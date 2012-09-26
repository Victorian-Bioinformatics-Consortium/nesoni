"""

Tools to do with variant calling.

"""


import random, os

import nesoni
from nesoni import config, legion, io, bio, working_directory, workflows, bowtie
from nesoni.third_party import vcf

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
@config.Main_section('variants', 'Variants, each with a number of reads\neg ACGx10 ATGx5')
@config.Configurable_section('analysis', presets=_analysis_presets)
class Test_variant_call(config.Action_with_output_dir):
    ref = None
    variants = [ ]
    analysis = _analysis_presets[0][1]()
    
    def run(self):
        workspace = self.get_workspace()
        
        read_length = 100
        left = rand_seq(read_length)
        right = rand_seq(read_length)
        
        with open(workspace/'reference.fa','wb') as f:
            io.write_fasta(f,'chr1',left+self.ref+right)
        
        i = 0
        with open(workspace/'reads.fq','wb') as f:
            for variant in self.variants:
                if 'x' in variant:
                    variant, count = variant.split('x')
                else:
                    count = 10
                count = int(count)
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
        
        reader = vcf.Reader(open(workspace/'variants.vcf','rU'))
        print dir(reader)
        print reader.formats
        print
        print reader.infos
        print
        for record in reader:
            print record
            variants = [ record.REF ]
            if isinstance(record.ALT,list):
                variants.extend(str(item) for item in record.ALT)
            else:
                variants.append(str(record.ALT))
            print variants
            
            print record.INFO
            for call in record.samples:
                print call.sample
                for key in call.data._fields:
                    print key, getattr(call.data,key), reader.formats[key].desc
                    
                counts = [ call.data.RO ]
                if isinstance(call.data.QA,list):
                    counts.extend(call.data.QA)
                else:
                    counts.append(call.data.QA)
                print variants, counts
                    
            print
