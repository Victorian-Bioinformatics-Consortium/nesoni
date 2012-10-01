"""

Tools to do with variant calling.


TODO:

test-variant-call should use vcf-patch-in

"""


import random, os, re, sys, collections

import nesoni
from nesoni import config, legion, io, bio, workspace, working_directory, reference_directory, workflows, bowtie, sam, reporting
from nesoni.third_party import vcf

@config.help("""
Run FreeBayes on a set of samples.

For high quality variant calling in bacteria, \
we suggest setting --ploidy greater than 1 here, \
then reducing the ploidy to 1 using "nesoni vcf-filter".
""")
@config.Main_section('samples', 
    'Working directories or BAM files. '
    'The read group should be set '
    '(nesoni\'s wrappers of read aligners do this as of version 0.87).'
    )
@config.Int_flag('ploidy', 'Ploidy of genotype calls.')
class Freebayes(config.Action_with_prefix):
    ploidy = 2
    samples = [ ]

    def run(self):
        bams = [ ]
        reference = None
        reference2 = None
        
        for sample in self.samples:
            if sam.is_bam(sample):
                bams.append(sample)
            elif os.path.isdir(sample):
                working = working_directory.Working(sample,True)
                bams.append( working.get_filtered_sorted_bam() )
                if reference2 is None:
                    reference2 = working.get_reference().reference_fasta_filename()
            elif io.is_sequence_file(sample):
                assert reference is None, 'Only one reference FASTA file allowed.'
                reference = sample
        
        if reference is None:
            reference = reference2
        if reference is None:
            raise grace.Error('No reference FASTA file given.')
        
        command = [ 
            'freebayes',
            '-f', reference,
            #'--no-population-priors',
            '--ploidy',str(self.ploidy),
            '--pvar','0.0',
            ] + bams
        
        self.log.log('Running: '+' '.join(command)+'\n')
        
        io.execute(
            command,
            stdout = open(self.prefix + '.vcf','wb'),
            )

class _Filter(Exception): pass

@config.help("""\
Filter a VCF file, eg as produced by "nesoni freebayes:".

- reduce ploidy of genotype calls
  eg reducing ploidy from 4 to 2 would
     1/1/1/1 -> 1/1
     1/1/2/2 -> 1/2
     1/2/2/2 -> ./. (fail)

- genotype calls that don't pass filters (eg genotype quality) are failed

- variants in which genotyping failed in all samples are removed

Note:

- phased VCF files are not supported

""")
@config.Positional('vcf', 'VCF file, eg as produced by "nesoni freebayes:"')
@config.Float_flag('min_gq', 'Genotype quality cutoff, based on GQ as produced by Freebayes, phred scale.')
@config.Int_flag('ploidy', 
    'Reduce to this polidy. Should be a divisor of the ploidy of the input. '
    'If the genotype does not have an exact representation at the reduced ploidy, it is filtered.')
class Vcf_filter(config.Action_with_prefix):
    vcf = None
    min_gq = 20.0
    ploidy = 1
    
    def _blank_gt(self):
        return '/'.join(['.']*self.ploidy)
    
    def _modify_sample(self, variants, sample):
        try:
            if '.' in sample.data.GT: 
                raise _Filter()
            if self.min_gq is not None and sample.data.GQ < self.min_gq:
                raise _Filter()
        
            assert re.match('[0-9]+(/[0-9]+)*$', sample.data.GT), 'Unsupported genotype format: '%sample.data.GT
        
            gt = sorted(int(item) for item in sample.data.GT.split('/'))
            assert len(gt) % self.ploidy == 0, 'Can\'t reduce ploidy from %d to %d' % (len(gt),self.ploidy)
            stride = len(gt) // self.ploidy
        
            new_gt = [ ]
            for i in xrange(0,len(gt),stride):
                if len(set(gt[i:i+stride])) != 1:
                    raise _Filter()
                new_gt.append(gt[i])        

            new_gt = '/'.join(str(item) for item in new_gt)        
        except _Filter:
            new_gt = self._blank_gt()
        sample.data = sample.data._replace(GT=new_gt)
            
    
    def run(self):
        reader = vcf.Reader(open(self.vcf,'rU'))
        
        writer = vcf.Writer(open(self.prefix + '.vcf','wb'), reader)
        
        #print dir(reader)
        #print reader.formats
        #print 
        #print reader.infos
        #print 
        for record in reader:
            print record
            variants = [ record.REF ]
            if isinstance(record.ALT,list):
                variants.extend(str(item) for item in record.ALT)
            else:
                variants.append(str(record.ALT))
            print variants
            
            print record.INFO
            
            any = False
            
            for sample in record.samples:
                self._modify_sample(variants, sample)
                
                any = any or sample.data.GT != self._blank_gt()

                #print call.sample
                #for key in call.data._fields:
                #    print key, getattr(call.data,key), reader.formats[key].desc
                #    
                #counts = [ call.data.RO ]
                #if isinstance(call.data.QA,list):
                #    counts.extend(call.data.QA)
                #else:
                #    counts.append(call.data.QA)
                #print variants, counts
                #
                #
                #if self.min_gq is not None and call.data.GQ < self.min_gq:
                #    call.data = call.data._replace(GT='.')
                #    print call.data
                #else:
                #    any = True
            
            #if any:
            writer.write_record(record)
                
        writer.close()

@config.help("""\
Patch variants in a VCF file into the reference sequence, \
producing genome sequences for each sample in the VCF file.

Only haploid variants are supported.

Note that in parts of the sequence with insufficient coverage to call variants \
or in which the variant could not be called for some other reason, \
the reference sequence will incorrectly be retained.
""")
@config.Positional('reference', 'Reference directory')
@config.Positional('vcf', 'VCF file')
class Vcf_patch(config.Action_with_output_dir):
    reference = None
    vcf = None
    
    def run(self):
        workspace = self.get_workspace()
        
        reference = reference_directory.Reference(self.reference, must_exist=True)
        
        reader = vcf.Reader(open(self.vcf,'rU'))
        variants = collections.defaultdict(list)
        for record in reader:
            variants[record.CHROM].append(record)
        
        for chrom in variants:
            variants[chrom].sort(key=lambda item: item.POS)
        
        filenames = [ workspace/(item+'.fa') for item in reader.samples ]
        for filename in filenames:
            with open(filename,'wb'): pass
        
        for name, seq in io.read_sequences(reference.reference_fasta_filename()):
            for i, sample in enumerate(reader.samples):            
                revised = [ ]
                pos = 0
                for variant in variants[name]:
                    gt = variant.samples[i].data.GT
                    if gt is None: continue
                    assert gt.isdigit(), 'Unsupported genotype (can only use haploid genotypes): '+gt
                    gt_number = int(gt)
                    if gt_number == 0:
                        var_seq = variant.REF
                    else:
                        var_seq = str(variant.ALT[gt_number-1])
                        assert re.match('[ACGTN]*$', var_seq), 'Unsupported variant type: '+var_seq
                    new_pos = variant.POS-1
                    revised.append(seq[pos:new_pos])
                    pos = new_pos
                    revised.append(var_seq)
                    assert seq[pos:pos+len(variant.REF)].upper() == variant.REF, 'REF column in VCF does not match reference sequence'
                    pos += len(variant.REF)
                revised.append(seq[pos:])
                            
                with open(filenames[i],'ab') as f:
                    io.write_fasta(f, name, ''.join(revised))

            del variants[name]        
        assert not variants, 'Chromosome names in VCF not in reference: '+' '.join(variants)
        
        


def rand_seq(n):
    return ''.join( 
        random.choice('ACGT') 
        for i in xrange(n) 
    )

_analysis_presets = [
    ('shrimp', workflows.Analyse_sample(clip=None,reconsensus=None), 'Do analysis using SHRiMP'),
    ('bowtie', workflows.Analyse_sample(clip=None,align=bowtie.Bowtie(),reconsensus=None), 'Do analysis using Bowtie'),
]

@config.Positional('ref', 'Reference sequence\neg AAG')
@config.Main_section('variants', 
    'Variants, each with a number of reads, eg ACGx10 ATGx5. '
    'The first variant given is treated as the correct variant.')
@config.Configurable_section('analysis', presets=_analysis_presets)
class Test_variant_call(config.Action_with_output_dir):
    ref = None
    variants = [ ]
    analysis = _analysis_presets[0][1]()
    
    def run(self):
        workspace = self.get_workspace()
        
        read_length = 100
        left = rand_seq(read_length-1)
        while True:
            flank = rand_seq(1)
            if flank != self.ref[:1]: break
        left += flank
        
        right = rand_seq(read_length-1)
        while True:
            flank = rand_seq(1)
            if flank != self.ref[-1:]: break
        right = flank+right
        
        with open(workspace/'reference.fa','wb') as f:
            io.write_fasta(f,'chr1',left+self.ref+right)
        
        i = 0
        
        variants_used = [ ]
        
        with open(workspace/'reads.fq','wb') as f:
            for variant in self.variants:
                if 'x' in variant:
                    variant, count = variant.split('x')
                    count = int(count)
                else:
                    count = 10
                variants_used.append( (variant,count) )
                for j in xrange(count):
                    seq = left+variant+right
                    
                    pos = len(variant)+random.randrange(read_length-len(variant))
                    read = seq[pos:pos+read_length]
                    if random.randrange(2):
                        read = bio.reverse_complement(read)
                    i += 1
                    io.write_fastq(f,'read_%s_%d' % (variant,i),read,chr(64+30)*len(read))
        
        legion.remake_needed()
        
        self.analysis(
            workspace/'sample',
            workspace/'reference.fa',
            reads = [ workspace/'reads.fq' ],
        ).run()
        
        Freebayes(
            workspace/'freebayes',
            workspace/'sample',
        ).run()
        
        Vcf_filter(
            workspace/'filtered',
            workspace/'freebayes.vcf',
        ).run()
        
        #nesoni.Shrimp(
        #    workspace/'align',
        #    workspace/'reference.fa',
        #    reads=[ workspace/'reads.fq' ],
        #    ).run()
        #
        #nesoni.Consensus(
        #    workspace/'align',
        #    ).run()
        #
        #nesoni.Make_genome(
        #    workspace/'genome',
        #    workspace/'align/reference'
        #    ).run()
        #
        #output_dir = self.output_dir                
        
        #command = (
        #   'samtools mpileup -f %(output_dir)s/align/reference/reference.fa %(output_dir)s/align/alignments_filtered_sorted.bam'
        #   ' |java -jar VarScan.v2.3.2.jar mpileup2cns --variants 1 --output-vcf 1'
        #   % locals()
        #)
        #print command
        #assert 0 == os.system(command)
        
        with open(workspace/'variants.vcf','wb') as f:
            io.execute([
                'freebayes',
                workspace/('sample','alignments_filtered_sorted.bam'),
                '-f', workspace/('sample','reference','reference.fa'),
                '--no-population-priors',
                '--ploidy','4',
                '--pvar','0.0',
            ], stdout=f)
        
        #io.execute(['igvtools','index',workspace/'variants.vcf'])
        
        reader = vcf.Reader(open(workspace/'filtered.vcf','rU'))
        
        good_count = 0
        bad_count = 0
        
        for record in reader:
            variants = [ record.REF ]
            if isinstance(record.ALT,list):
                variants.extend(str(item) for item in record.ALT)
            else:
                variants.append(str(record.ALT))
            sample = record.samples[0]
            pos = record.POS-1
            if sample.data.GT and sample.data.GT.isdigit():
                call = variants[ int(sample.data.GT) ]
            else:
                call = None

            print pos, variants, call, variants_used[0][0]            
            while len(set( item[:1] for item in variants )) == 1:
                variants = [ item[1:] for item in variants ]
                pos += 1

            if sample.data.GT and sample.data.GT.isdigit():
                call = variants[ int(sample.data.GT) ]
            else:
                call = None
            print pos, variants, call, variants_used[0][0]            
            
            if call is not None:
                good = pos == read_length and variants[0] == self.ref and call == variants_used[0][0]
                if good:
                    good_count += 1
                else:
                    bad_count += 1

        self.log.log('\n')
        self.log.datum(workspace.name,'correct variant', good_count > 0)
        self.log.datum(workspace.name,'incorrect variants', bad_count)
        self.log.log('\n')


@config.help("""
Assess ability to call variants under a spread of different conditions.
""")
@config.Configurable_section('template','Setting for "nesoni test-variant-call".')
class Power_variant_call(config.Action_with_output_dir):
    template = Test_variant_call()

    def tryout(self, ref, variants):
        work = self.get_workspace()
        with workspace.tempspace(work.working_dir) as temp:
            job = self.template(temp.working_dir, ref=ref, variants=variants)            
            job.run()            
            
            result = dict( tuple(item.values()) for item in reporting.mine_logs([job.log_filename()]) )
            good = {'yes':True,'no':False}[result['correct variant']]
            bad = int(result['incorrect variants'])
            return good, bad

    def depth_test(self, ref, variant, others=[]):
        depths = range(1,30+1)
        results = [ self.tryout(ref, ['%sx%d'%(variant,i)] + others) for i in depths ]
        
        report = '%s -> %s' % (ref or '-', variant or '-')
        if others: report += '  with contamination  ' + ', '.join(others)
        report += '\n'
        
        report += 'good  [' + ''.join( '+' if item[0] else ' ' for item in results ) + ']\n'
        report += 'bad   [' + ''.join( 'x' if item[1] else ' ' for item in results ) + ']\n'

        depth_line = ''
        for i in depths:
            if i == 1 or i % 5 == 0:
                depth_line += ' '*(i-len(depth_line)) + str(i)
        report += 'depth '+depth_line+'\n'
        
        return report
        

    def run(self):
        work = self.get_workspace()
        
        #print [ self.tryout('A', ['CCCCCx%d'%i]) for i in xrange(1,11) ]
        
        report = ''
        
        for job in [
            ('A','C',[]),
            ('A','C',['Gx1']),
            ('A','',[]),
            ('','A',[]),
        ]:        
            report += self.depth_test(*job)+'\n'
        
        print report



