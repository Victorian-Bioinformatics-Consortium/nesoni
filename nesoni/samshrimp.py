
from nesoni import io, bio, grace, sam, config, legion, working_directory

import os, sys, subprocess


@config.help("""\
Align reads to a reference using SHRiMP.

Paired end reads should either be given as two files in the "pairs" section, \
or as a single interleaved file in the "inerleaved" section.

To increase sensitivity, shrimp is by default invoked with the options:

  -n 2 -w 200% --half-paired --sam-unaligned

By default pairs are assumed to have:

  -p opp-in -I 0,500
  
You can override these in the shrimp-options section.

Useful shrimp options:
  
  -h threshold% - Reduce the hit threshold
                  (default is 55% of read length)

  -n 1          - For single reads: perform alignment on a single 
                  seed hit (default 2) 
                  (more sensitive, but slower)

  -p, -I        - Pairing parameters, see SHRiMP documentation
  
  --qv-offset   - Quality offset.
                  If not specified, nesoni will attempt to guess.
                  For Sanger format fastq and recent Illumina, use 33.
""")
@config.Main_section('references', 
    'Reference sequence filenames, '
    'or a directory created using "nesoni make-reference:" (recommended).')
@config.Bool_flag('cs', 'Are reads in colorspace?')
@config.Bool_flag('sam_unaligned', 'Pass --sam-unaligned to gmapper?')
@config.Bool_flag('half_paired', 'Pass --half-paired to gmapper?')
@config.Int_flag('cores', 'Maximum cores to use.', affects_output=False)
@config.Section('reads', 'Files containing unpaired reads.')
@config.Section('interleaved', 'Files containing interleaved read pairs.')
@config.Grouped_section('pairs', 'Pair of files containing read pairs.')
@config.Section('shrimp_options', 'Options to pass to SHRiMP.', allow_flags=True)
class Shrimp(config.Action_with_output_dir):
    cs = False
    sam_unaligned = True
    half_paired = True
    cores = 8
    references = []
    reads = []
    interleaved = []
    pairs = []
    shrimp_options = []
    
    _workspace_class = working_directory.Working

    def run(self):
        grace.require_shrimp_2()
        grace.require_samtools()
        assert self.references, 'No reference sequences given'
        assert self.reads or self.pairs or self.interleaved, 'No reads given'
        for pair in self.pairs:
            assert len(pair) == 2, 'Two files required in each pair: section'

        io.check_name_uniqueness(self.reads, self.pairs, self.interleaved)

        read_sets = [ ]
        for item in self.reads:
            read_sets.append( ([item], False) )
        for item in self.pairs:
            read_sets.append( (item, True) )
        for item in self.interleaved:
            read_sets.append( ([item], True) )

        #Create working directory
        
        workspace = self.get_workspace()
        workspace.setup_reference(self.references)        
        workspace.update_param(snp_cost=25)
        reference = workspace.get_reference()
        reference_filename = reference.reference_fasta_filename()

        cores = min(self.cores, legion.coordinator().get_cores())
                
        default_options = { 
            '-E' : None, 
            '-T' : None, 
            '-N' : str(cores), 
            '-n':'2', 
            '-w':'200%',
            '-p': 'opp-in', 
            '-I': '0,500', 
            '-X':None,
        }
        
        if self.sam_unaligned:
            default_options['--sam-unaligned'] = None
        
        if self.half_paired:
            default_options['--half-paired'] = None
        else:
            default_options['--no-half-paired'] = None

        cutoff = '55%' #Default changed in SHRiMP 2.0.2
        if '-h' in self.shrimp_options:
            cutoff = self.shrimp_options[ self.shrimp_options.index('-h')+1 ]
        
        #Run shrimp
        
        bam_filename = io.abspath(self.output_dir, 'alignments.bam')
        bam_prefix = io.abspath(self.output_dir, 'alignments')
        bam_sorted_prefix = io.abspath(self.output_dir, 'alignments_sorted')
        
        temp_filename = io.abspath(self.output_dir, 'temp.bam')
        
        log_filename = io.abspath(self.output_dir, 'shrimp_log.txt')
        log_file = open(log_filename, 'wb')
        
        sam_eater = sam.Bam_writer(temp_filename)
        
        sam_header_sent = [False]
        n_seen = [0]
        
        def eat(f):
            for line in f:
                if line.startswith('@'):
                    if sam_header_sent[0]: continue
                else:
                    n_seen[0] += 1
                    if n_seen[0] % 100000 == 0:
                        grace.status('%s alignments produced' % grace.pretty_number(n_seen[0]))
                sam_eater.write_raw(line)
            sam_header_sent[0] = True
        
        def remove_pair_options(options):
            for flag in ['-p','-I']:
                while flag in options:
                    pos = options.index(flag)
                    options = options[:pos] + options[pos+2:]
            for flag in ['--half-paired']:
                while flag in options:
                    pos = options.index(flag)
                    options = options[:pos] + options[pos+1:]
            return options
        
        for i, (filenames, is_paired) in enumerate(read_sets):
            options = self.shrimp_options[:]
               
            has_qualities = all(
                len( io.read_sequences(filename, qualities=True).next() ) == 3  #A little ugly
                for filename in filenames
            )
            if has_qualities:
                options.append( '--fastq' )
            
            if len(filenames) == 1:
                reads_parameters = [ filenames[0] ]
            else:
                reads_parameters = [ '-1', filenames[0], '-2', filenames[1] ]
            
            if '--qv-offset' not in self.shrimp_options:
                #guesses = [ ]
                #for filename in filenames:
                #    guesses.append(io.guess_quality_offset(filename))
                #assert len(set(guesses)) == 1, 'Conflicting quality offset guesses, please specify --qv-offset manually.'
                #default_options['--qv-offset'] = str(guesses[0])
                default_options['--qv-offset'] = str( io.guess_quality_offset(*filenames) )
                
            default_options['--read-group'] = '%s,%s' % (
                workspace.name.replace(',','_'),
                workspace.name.replace(',','_')
            )
            for flag in default_options:
                if flag not in options:
                    options.append(flag)
                    if default_options[flag] is not None:
                        options.append(default_options[flag])
            
            if not is_paired:
               options = remove_pair_options(options)
            
            grace.status('')
            
            full_param = reference.shrimp_command(self.cs, options + reads_parameters)
            
            print >> sys.stderr, 'Running', ' '.join(full_param)
            
            with io.pipe_from(full_param,
                    stderr=log_file,
                    cores=cores) as f:
                eat(f)
        
        log_file.close()
        
        sam_eater.close()
        
        grace.status('Sort')
        
        #io.execute([
        #    'samtools', 'sort', '-n', temp_filename, bam_prefix
        #])
        sam.sort_bam(temp_filename, bam_prefix, by_name=True, cores=self.cores)
        
        os.unlink(temp_filename)
        
        grace.status('')

