
"""

Produce files for use with IGV



"""

from nesoni import config, io, grace, legion, trivia, annotation, working_directory, reference_directory

import itertools, math, os

def iter_add(a,b):
    for a,b in itertools.izip(a,b):
        yield a+b




@config.help("""\
Create an IGV .genome file from sequences and annotations, \
or from a nesoni reference directory created with "nesoni make-reference:".
""")
@config.String_flag('name', 'Descriptive name.')
@config.Main_section('filenames', 'Input filenames or reference directory.')
class Make_genome(config.Action_with_prefix):
    prefix = None
    name = None
    filenames = [ ]
    
    def run(self):
        base = os.path.split(self.prefix)[1]
        
        annotations = [ ]
        sequences = [ ]
        
        for filename in self.filenames:
            any = False
            if os.path.isdir(filename):
                reference = reference_directory.Reference(filename,must_exist=True)
                sequences.append(reference.reference_fasta_filename())
                if reference.annotations_filename():
                    annotations.append(reference.annotations_filename())
                any = True
            if io.is_sequence_file(filename):
                sequences.append(filename)
                any = True
            if annotation.is_annotation_file(filename):
                annotations.append(filename)
                any = True
            assert any, 'File is neither a recognized sequence or annotation file'

        cytoband_filename = os.path.join(self.prefix,base+'_cytoband.txt')
        property_filename = os.path.join(self.prefix,'property.txt')
        gff_filename = os.path.join(self.prefix,base+'.gff')
        output_filenames = [ cytoband_filename, property_filename, gff_filename ] 

        if not os.path.exists(self.prefix):
            os.mkdir(self.prefix)
            
        f = open(property_filename,'wb')
        print >> f, 'ordered=true'
        print >> f, 'id=%s' % base
        print >> f, 'name=%s' % (self.name or base)
        print >> f, 'cytobandFile=%s_cytoband.txt' % base
        print >> f, 'geneFile=%s.gff' % base
        print >> f, 'sequenceLocation=%s' % base
        f.close()
        
        trivia.As_gff(output=gff_filename,
               filenames=annotations,
               exclude=[ 'gene', 'source' ]
        ).run()
        
        f_cyt = open(cytoband_filename,'wb')
        for filename in sequences:
            for name, seq in io.read_sequences(filename):
                assert '/' not in name
                f = open(os.path.join(self.prefix, name + '.txt'), 'wb')
                f.write(seq)
                f.close()
                print >> f_cyt, '%s\t0\t%d' % (name, len(seq))
        f_cyt.close()
        
        genome_filename = self.prefix + '.genome'
        if os.path.exists(genome_filename):
            os.unlink(genome_filename)
        io.execute(
            ['zip', '-j', io.abspath(genome_filename)] +
            [ io.abspath(item) for item in output_filenames ]
        )
        for filename in output_filenames:
            if os.path.exists(filename):
                os.unlink(filename)
        



@config.help("""\
Produce IGV format plots.

Ambiguous mappings are including. Use the "mapping-ambiguity.tdf" plot \
to see where this might be an issue.

Note: "nesoni as-fasta:" and "nesoni as-gff:" can be used to create the files \
necessary to import a genome into IGV, if you only have a GENBANK file.
""")
@config.Bool_flag('strand_specific', 'Output strand-specific plots.') 
@config.Bool_flag('raw', 'Output un-normalized plots.')
@config.Bool_flag('norm', 'Output normalized plots.') 
@config.String_flag('genome', 'IGV ".genome" file. If specified, igvtools will be used to create ".tdf" files.')
@config.Bool_flag('delete_igv', 'Delete the ".igv" files after converting to ".tdf".')
@config.Main_section('working_dirs', 'Working directories containing the results of "filter:" or "consensus:".')
@config.String_flag('norm_file', 'File of normalizations produced by "norm-from-counts:". Defaults to equalizing the average unambiguous depth.')
#@config.Float_section('norm_mult', 
#    'Normalize by multiplying corresponding sample depths by these multipliers. Defaults to equalizing the average unambiguous depth. '
#    'You might alternately use EdgeR\'s Trimmed Mean Normalization, eg in R: library(nesoni); read.counts(\'counts.txt\', min.total=10)$samples$normalizing.multiplier')
class IGV_plots(config.Action_with_prefix):
    strand_specific = True
    raw = True
    norm = True
    genome = None
    delete_igv = True

    working_dirs = [ ]

    norm_file = None    
    
    def iter_over(self, access_func, zeros=True):
        for name in self.chromosome_names:
            this_depths = [ iter(access_func(item[name])) for item in self.depths ]
            for pos, item in enumerate(itertools.izip(*this_depths)):
                if zeros or any(item):
                    yield name, pos, item
    
    def iter_over_unstranded(self, access_func, zeros=True):
        return self.iter_over(lambda item: iter_add(access_func(item)[0],access_func(item)[1]), zeros)        

    def calculate_norm_mult(self):
        grace.status('Calculating normalization')
        
        totals = [ 0 ] * len(self.working_dirs)
        
        for i in xrange(len(self.workspaces)):
            for name in self.lengths:
                totals[i] += self.depths[i][name].depths[0].total() + self.depths[i][name].depths[1].total()
        
        #for name, pos, depths in self.iter_over_unstranded(lambda item: item.depths):
        #    if pos % 1000000 == 0:
        #        grace.status('Calculating normalization %s %s' % (name, grace.pretty_number(pos)))
        #    for i, depth in enumerate(depths):
        #        totals[i] += depth
        
        grace.status('')
        
        nonzero = [ item for item in totals if item ]
        geomean = math.exp(sum( math.log(item) for item in nonzero ) / len(nonzero))
        self.norm_mult = [
            1.0 if not item else geomean / item
            for item in totals
        ]

    def load_norm_mult(self):
        mults = { }
        for record in io.read_table(self.norm_file):
            mults[record['Sample']] = float(record['Normalizing.multiplier'])        
        self.norm_mult = [ mults[name] for name in self.sample_names ]
    
    def normalizer(self):
        #norm_factors = self.norm_factors
        #return lambda depths: [ a*b for a,b in zip(depths, norm_factors) ]
        
        expr = 'lambda depths: ('+''.join( 
            'depths[%d]*%f,' % item
            for item in enumerate(self.norm_mult)
        )+')'
        return eval(expr)
        
    def normalize_iter(self, iterator):
        normalize = self.normalizer()
        for name, pos, depths in iterator:
            yield name, pos, normalize(depths)

    def find_maximum_depth(self):
        grace.status('Finding maximum depth')        
        if self.strand_specific:
            iterator = legion.interleave([
                self.iter_over(lambda item: item.ambiguous_depths[0], zeros=False),
                self.iter_over(lambda item: item.ambiguous_depths[1], zeros=False),
            ])
        else:
            iterator = self.iter_over_unstranded(lambda item: item.ambiguous_depths, zeros=False)

        maximum = 1
        norm_maximum = 1.0
        normalize = self.normalizer()
        n = 0
        #for name, pos, depths in iterator:
        #    if n % 1000000 == 0:
        #        grace.status('Finding maximum depth %s %s' % (name, grace.pretty_number(pos)))
        #    n += 1
        #    
        #    maximum = max(maximum,max(depths))
        #    norm_maximum = max(norm_maximum,max(normalize(depths)))

        futures = [ ]        
        for name in self.chromosome_names:
            for i, depth in enumerate(self.depths):
                if self.strand_specific:
                    #this = futures.append(max(
                    #    max(depth[name].ambiguous_depths[0]),
                    #    max(depth[name].ambiguous_depths[1])
                    #)
                    
                    futures.append( (name,i,legion.future(max, depth[name].ambiguous_depths[0])) )
                    futures.append( (name,i,legion.future(max, depth[name].ambiguous_depths[1])) )
                else:
                    #this = max(iter_add(depth[name].ambiguous_depths[0],depth[name].ambiguous_depths[1]))
                    futures.append( (name,i,legion.future(lambda item: max(iter_add(item[0],item[1])), depth[name].ambiguous_depths)) )

        for name, i, future in futures:
            grace.status('Finding maximum depth %s %s' % (name, self.sample_names[i]))
            this = future()
            maximum = max(maximum, this)
            norm_maximum = max(maximum, self.norm_mult[i] * this)
        
        self.maximum = maximum
        self.norm_maximum = norm_maximum
        grace.status('')

    def make_plot(self, plots_name, plot_names, iterator, maximum, color='0,0,0', scale_type='log', windowing='maximum'):
        grace.status('Write '+plots_name)
        filename = self.prefix + plots_name + '.igv'
        f = open(filename, 'wb')
        
        height = max(10,int(100.0/math.sqrt(len(plot_names)))) #...
        
        print >> f, '#track viewLimits=0:%(maximum)f autoScale=off scaleType=%(scale_type)s windowingFunction=%(windowing)s maxHeightPixels=200:%(height)d:1 color=%(color)s' % locals()
        print >> f, '\t'.join(
            [ 'Chromosome', 'Start', 'End', 'Feature'] + plot_names
        )
        for name, pos, depths in iterator:
            print >> f, '\t'.join(
                [ name, str(pos), str(pos+1), 'F' ] +
                [ str(item) for item in depths ]
            )
        
        f.close()
        grace.status('')
        
        if self.genome:
            p = io.run([
                'igvtools', 'tile',
                filename, 
                self.prefix + plots_name + '.tdf',
                self.genome,
                '-f', 'max,mean'
            ], stdin=None, stdout=None)
            self.processes.append((p, filename))
            
    def wait_for_igv(self):
        while self.processes:
           p, filename = self.processes.pop()
           assert p.wait() == 0, 'igvtools tile failed'
           if self.delete_igv:
               os.unlink(filename)
                
    def iter_ambiguity(self):
        for (name, pos, depths), (name1, pos1, ambiguous_depths) in itertools.izip(
            self.iter_over_unstranded(lambda item: item.depths),
            self.iter_over_unstranded(lambda item: item.ambiguous_depths)):
            total_ambiguous = sum(ambiguous_depths)
            if not total_ambiguous:
                yield name, pos, (0.0,)
            else:
                yield name, pos, (1.0 - float(sum(depths)) / total_ambiguous,)

    def iter_total(self):
        for (name,pos,depths_fwd), (name1,pos1,depths_rev) in itertools.izip(
            self.iter_over(lambda item: item.ambiguous_depths[0]),
            self.iter_over(lambda item: item.ambiguous_depths[1])
        ):
            yield name, pos, (sum(depths_fwd), sum(depths_rev))

    def iter_5prime(self):
        for (name,pos,depths_fwd), (name1,pos1,depths_rev) in itertools.izip(
            self.iter_over(lambda item: item.ambiguous_depths[0].iter_starts()),
            self.iter_over(lambda item: item.ambiguous_depths[1].iter_ends())
        ):
            yield name, pos, (sum(depths_fwd), sum(depths_rev))

    def iter_3prime(self):
        for (name,pos,depths_fwd), (name1,pos1,depths_rev) in itertools.izip(
            self.iter_over(lambda item: item.ambiguous_depths[0].iter_ends()),
            self.iter_over(lambda item: item.ambiguous_depths[1].iter_starts())
        ):
            yield name, pos, (sum(depths_fwd), sum(depths_rev))

    def setup(self):
        grace.status('Load depths')
        self.sample_names = [ os.path.split(dirname)[1] for dirname in self.working_dirs ]
        self.workspaces = [ working_directory.Working(dirname, must_exist=True) for dirname in self.working_dirs ]
        
        self.depths = [ item.get_object('depths.pickle.gz') for item in self.workspaces ]
        #self.depths = list(legion.imap(lambda item: item.get_object('depths.pickle.gz'), self.workspaces, local=True))
        
        self.any_pairs = any(item.param['any_pairs'] for item in self.workspaces)
        grace.status('')
        
        lengths = self.workspaces[0].get_reference().get_lengths()
        self.chromosome_names = [ name for name, length in lengths ]
        self.lengths = dict(lengths)
        
        self.processes = [ ]

    def setdown(self):
        del self.sample_names
        del self.workspaces
        del self.depths
        del self.any_pairs
        del self.lengths
        self.wait_for_igv()
        del self.processes        

    def run(self):                
        assert self.working_dirs, 'No working directories given.'

        self.setup()
        
        if self.norm_file:
            self.load_norm_mult()
        else:
            self.calculate_norm_mult()
        
        self.find_maximum_depth()
        
        names_fwd = [ item + ' fwd' for item in self.sample_names ]
        names_rev = [ item + ' rev' for item in self.sample_names ]

        if self.raw:
            if self.strand_specific:
                self.make_plot('-reads-raw-fwd', names_fwd, self.iter_over(lambda item: item.ambiguous_depths[0]), self.maximum, '0,128,0')        
                self.make_plot('-reads-raw-rev', names_rev, self.iter_over(lambda item: item.ambiguous_depths[1]), self.maximum, '0,0,128')                
                if self.any_pairs:        
                    self.make_plot('-fragments-raw-fwd', names_fwd, self.iter_over(lambda item: item.ambiguous_pairspan_depths[0]), self.maximum, '0,128,0')        
                    self.make_plot('-fragments-raw-rev', names_rev, self.iter_over(lambda item: item.ambiguous_pairspan_depths[1]), self.maximum, '0,0,128')                
            else:
                self.make_plot('-reads-raw', self.sample_names, self.iter_over_unstranded(lambda item: item.ambiguous_depths), self.maximum, '0,128,128')
                if self.any_pairs:
                    self.make_plot('-fragments-raw', self.sample_names, self.iter_over_unstranded(lambda item: item.ambiguous_pairspan_depths), self.maximum, '0,128,128')

        if self.norm:
            if self.strand_specific:
                self.make_plot('-reads-normalized-fwd', names_fwd, self.normalize_iter(self.iter_over(lambda item: item.ambiguous_depths[0])), self.norm_maximum, '0,128,0')        
                self.make_plot('-reads-normalized-rev', names_rev, self.normalize_iter(self.iter_over(lambda item: item.ambiguous_depths[1])), self.norm_maximum, '0,0,128')
                if self.any_pairs:        
                    self.make_plot('-fragments-normalized-fwd', names_fwd, self.normalize_iter(self.iter_over(lambda item: item.ambiguous_pairspan_depths[0])), self.norm_maximum, '0,128,0')        
                    self.make_plot('-fragments-normalized-rev', names_rev, self.normalize_iter(self.iter_over(lambda item: item.ambiguous_pairspan_depths[1])), self.norm_maximum, '0,0,128')
            else:
                self.make_plot('-reads-normalized', self.sample_names, self.normalize_iter(self.iter_over_unstranded(lambda item: item.ambiguous_depths)), self.norm_maximum, '0,128,128')
                if self.any_pairs:
                    self.make_plot('-fragments-normalized', self.sample_names, self.normalize_iter(self.iter_over_unstranded(lambda item: item.ambiguous_pairspan_depths)), self.norm_maximum, '0,128,128')
        
        self.make_plot('-mapping-ambiguity', ['Ambiguity'], self.iter_ambiguity(), 1.0, '196,0,0', scale_type='linear', windowing='mean')

        maximum = self.maximum * len(self.working_dirs)
        self.make_plot('-total', ['Total fwd', 'Total rev'], self.iter_total(), maximum, '64,0,64')        
        self.make_plot('-5prime-ends', ['5p fwd', '5p rev'], self.iter_5prime(), maximum, '64,0,64')
        self.make_plot('-3prime-ends', ['3p fwd', '3p rev'], self.iter_3prime(), maximum, '64,0,64')
        
        self.setdown()
        

@config.help("""
Convert .igv files, for example as produced by "igv-plots: --delete-igv no" \
to .userplot files for use with Artemis.
""")
@config.Main_section('files', '.igv files.')
class As_userplots(config.Action_with_output_dir):
    files = [ ]
    
    def run(self):
        working = io.Workspace(self.output_dir, must_exist=False)
    
        for filename in self.files:
            reader = io.Table_reader(filename)
            
            name = os.path.splitext(os.path.split(filename)[1])[0]
            
            rname = None
            files = None
            for record in reader:
                if record['Chromosome'] != rname:
                    if files: 
                        for item in files: 
                            item.close()
                    rname = record['Chromosome']
                    grace.status('Convert '+name+' '+rname)
                    files = [
                        open(working / (
                            name + 
                            '-' + grace.filesystem_friendly_name(rname) + 
                            '-' + grace.filesystem_friendly_name(item) + '.userplot'
                        ), 'wb')
                        for item in reader.headings[4:]
                    ]
                    pos = 0
                assert int(record['Start']) == pos and int(record['End']) == pos + 1
                
                for val, f in zip(record.values()[4:], files):
                    print >> f, val
                
                pos += 1
            
            if files: 
                for item in files: 
                    item.close()
            grace.status('')

@config.help("""\
Open IGV with a specified .genome file. \
If the file is not in IGV's list of user defined genomes, it will be added.

.genome files can be created with "nesoni make-genome:".
""")
@config.Positional('genome','.genome file.')
@config.Main_section('args','command line arguments to IGV, eg files to load.', allow_flags=True)
class Run_igv(config.Action):
    genome = None
    args = [ ]

    def run(self):
        igv_dir = os.path.join(os.environ['HOME'],'igv')
        if not os.path.exists(igv_dir):
            os.mkdir(igv_dir)
        genomes_dir = os.path.join(igv_dir,'genomes')
        if not os.path.exists(genomes_dir):
            os.mkdir(genomes_dir)
        genomes_filename = os.path.join(
            os.environ['HOME'],
            'igv', 'genomes', 'user-defined-genomes.txt'
        )
        genomes = [ ]
        if os.path.exists(genomes_filename):
            with open(genomes_filename,'rU') as f:
                for line in f:
                    genomes.append(line.rstrip('\n').split('\t'))
    
        genome_filename = os.path.abspath(self.genome)
        
        name = None
        for item in genomes:
            if os.path.abspath(item[1]) == genome_filename:
                name = item[2]
                break
        
        if name is None:
            names = set(item[2] for item in genomes)
            i = 0
            while True:
                name = os.path.splitext(os.path.basename(genome_filename))[0]
                name += str(i).lstrip('0')
                if name and name not in names: break
                i += 1
            genomes.append([ name, genome_filename, name ])

            print 'Genome not in user defined genomes list.'
            print 'Adding it as:', name
        
            with open(genomes_filename,'wt') as f:
                for item in genomes:
                    print >> f, '\t'.join(item)
    
        #igv.sh does not begin with #!/bin/sh as at version 2.1.24
        #io.execute(['igv.sh','-g',name] + self.args)
        
        io.execute(['java','-jar',io.find_jar('igv.jar'),'-g',name] + self.args)
        
        



