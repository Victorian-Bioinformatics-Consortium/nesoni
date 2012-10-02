
from nesoni import io, bio, grace, sam, config, working_directory

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
    'or a directory created using "nesoni make-reference: --bowtie yes" (recommended).')
@config.Bool_flag('cs', 'Are reads in colorspace?')
@config.Bool_flag('sam_unaligned', 'Pass --sam-unaligned to gmapper?')
@config.Bool_flag('half_paired', 'Pass --half-paired to gmapper?')
@config.Section('reads', 'Files containing unpaired reads.')
@config.Section('interleaved', 'Files containing interleaved read pairs.')
@config.Grouped_section('pairs', 'Pair of files containing read pairs.')
@config.Section('shrimp_options', 'Options to pass to SHRiMP.', allow_flags=True)
class Shrimp(config.Action_with_output_dir):
    cs = False
    sam_unaligned = True
    half_paired = True
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
        
        workspace = self.get_workspace() #working_directory.Working(self.output_dir, must_exist=False)
        workspace.setup_reference(self.references)        
        workspace.update_param(snp_cost=25)
        reference = workspace.get_reference()
        reference_filename = reference.reference_fasta_filename()


        
        default_options = { 
            '-E' : None, 
            '-T' : None, 
            '-N' : str(grace.how_many_cpus()), 
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
                
        #workspace = io.Workspace(self.output_dir)
        #
        #workspace.update_param( 
        #    shrimp_cutoff = cutoff
        #)
        #
        ##Make copy of reference sequences
        #
        #reference_filename = io.abspath(self.output_dir,'reference.fa')
        #reference_file = open(reference_filename,'wb')
        #
        #reference_genbank_filename = io.abspath(self.output_dir,'reference.gbk')
        #reference_genbank_file = open(reference_genbank_filename,'wb')
        #any_genbank = [ False ]
        #
        #def genbank_callback(name, record):
        #    """ Make a copy of any genbank files passed in. """
        #    from Bio import SeqIO
        #    
        #    SeqIO.write([record], reference_genbank_file, 'genbank')
        #    
        #    f = open(os.path.join(
        #        self.output_dir,
        #        grace.filesystem_friendly_name(name) + '.gbk'
        #    ), 'wb')
        #    SeqIO.write([record], f, 'genbank')
        #    f.close()
        #    
        #    any_genbank[0] = True
        #
        #for filename in self.references:
        #    for name, sequence in io.read_sequences(filename, genbank_callback=genbank_callback):
        #        #Don't retain any comment
        #        name = name.split()[0]
        #        io.write_fasta(reference_file, name, sequence.upper())
        #        
        #        f = open(os.path.join(
        #            self.output_dir,
        #            grace.filesystem_friendly_name(name) + '.fa'
        #        ), 'wb')
        #        io.write_fasta(f, name, sequence.upper())
        #        f.close()
        #        
        #
        #reference_file.close()
        #reference_genbank_file.close()
        #if not any_genbank[0]:
        #    os.unlink(reference_genbank_filename)
        #
        ## Create an index of the reference sequences
        #io.execute([
        #    'samtools', 'faidx', reference_filename
        #])
        
        #Run shrimp
        
        bam_filename = io.abspath(self.output_dir, 'alignments.bam')
        bam_prefix = io.abspath(self.output_dir, 'alignments')
        bam_sorted_prefix = io.abspath(self.output_dir, 'alignments_sorted')
        
        temp_filename = io.abspath(self.output_dir, 'temp.bam')
        
        log_filename = io.abspath(self.output_dir, 'shrimp_log.txt')
        log_file = open(log_filename, 'wb')
        
        sam_eater = sam.Bam_writer(temp_filename)
        
        #if self.cs:
        #    program = 'gmapper-cs'
        #else:
        #    program = 'gmapper-ls'
        
        sam_header_sent = [False]
        n_seen = [0]
        
        def eat(process):
            for line in process.stdout:
                if line.startswith('@'):
                    if sam_header_sent[0]: continue
                else:
                    n_seen[0] += 1
                    if n_seen[0] % 100000 == 0:
                        grace.status('%s alignments produced' % grace.pretty_number(n_seen[0]))
                sam_eater.write_raw(line)
                
            assert process.wait() == 0, 'shrimp failed'
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
        
        if '--qv-offset' not in self.shrimp_options:
            guesses = [ ]
            for filenames, is_paired in read_sets:
                for filename in filenames:
                    guesses.append(io.guess_quality_offset(filename))
            assert len(set(guesses)) == 1, 'Conflicting quality offset guesses, please specify --qv-offset manually.'
            default_options['--qv-offset'] = str(guesses[0])
                
        for i, (filenames, is_paired) in enumerate(read_sets):
            options = self.shrimp_options[:]
               
            has_qualities = all(
                len( io.read_sequences(filename, qualities=True).next() ) == 3  #A little ugly
                for filename in filenames
            )
            if has_qualities:
                options.append( '--fastq' )
            #    temp_read_filename = io.abspath(working_dir, 'temp.fa')
            #else:
            #    temp_read_filename = io.abspath(working_dir, 'temp.fq')
            
            #try:
            
            #if len(filenames) == 1: # gmapper can cope with gzipped     and filenames[0].endswith('.fa') or filenames[0].endswith('.fq'):
            #    actual_read_filename = filenames[0]
            #else:
            #    actual_read_filename = temp_read_filename
            #    grace.status('Copying reads')
            #    f = open(temp_read_filename, 'wb')
            #    if has_qualities:
            #        for reads in itertools.izip(*[ io.read_sequences(filename, qualities=True) for filename in filenames ]):
            #            for name, seq, qual in reads:
            #                io.write_fastq(f, name, seq, qual)
            #    else:
            #        for reads in itertools.izip(*[ io.read_sequences(filename) for filename in filenames ]):
            #            for name, seq in reads:
            #                io.write_fasta(f, name, seq)
            #    f.close()
            #    grace.status('')
            
            if len(filenames) == 1:
                reads_parameters = [ filenames[0] ]
            else:
                reads_parameters = [ '-1', filenames[0], '-2', filenames[1] ]
            
            default_options['--read-group'] = '%s-%d,%s' % (
                workspace.name.replace(',','_'),
                i+1,
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
            
            p = io.run(full_param,
                    stdout=subprocess.PIPE,
                    stderr=log_file)
            eat(p)
                        
            #finally:
            #    if os.path.exists(temp_read_filename):
            #        os.unlink(temp_read_filename)
        
        log_file.close()
        
        sam_eater.close()
        
        grace.status('Sort')
        
        io.execute([
            'samtools', 'sort', '-n', temp_filename, bam_prefix
        ])
        
        os.unlink(temp_filename)
        
        grace.status('')






#
#def samshrimp_main(args):
#    grace.require_shrimp_2()
#    grace.require_samtools()
#
#    original_args = args
#
#    use_cs, args = grace.get_option_value(args,'--cs', grace.as_bool, False)
#    use_sam_unaligned, args = grace.get_option_value(args,'--sam-unaligned', grace.as_bool, True)
#    use_half_paired, args = grace.get_option_value(args,'--half-paired', grace.as_bool, True)
#
#    if not args:
#        print SHRIMP_HELP
#        raise grace.Help_shown()
#
#    working_dir = args[0]
#    args = args[1:]
#    
#    reference_filenames = [ ]
#    read_sets = [ ]
#    
#    shrimp_options = [ ]
#    default_options = { '-E' : None, '-T' : None, '-N' : str(grace.how_many_cpus()), '-n':'2', '-w':'200%',
#                        '-p': 'opp-in', '-I': '0,500', '-X':None }
#
#    if use_sam_unaligned:
#        default_options['--sam-unaligned'] = None
#    
#    if use_half_paired:
#        default_options['--half-paired'] = None
#    else:
#        default_options['--no-half-paired'] = None
#        
#    def default_command(args):
#        grace.expect_no_further_options(args)
#        reference_filenames.extend(args)
#    def reads_command(args):
#        grace.expect_no_further_options(args)
#        read_sets.extend([ ([io.abspath(item)], False) for item in args ])
#    def pairs_command(args):
#        grace.expect_no_further_options(args)
#        assert len(args) == 2, 'Expected a pair of files for "pairs" section'
#        read_sets.append((args, True))
#    def interleaved_command(args):
#        grace.expect_no_further_options(args)
#        read_sets.extend([ ([io.abspath(item)], True) for item in args ])
#    def shrimp_options_command(args):
#        shrimp_options.extend(args)
#    
#    grace.execute(args, {
#        'reads' : reads_command,
#        'interleaved' : interleaved_command,
#        'pairs' : pairs_command,
#        'shrimp-options' : shrimp_options_command,
#    }, default_command)
#    
#    assert reference_filenames, 'No reference files given'
#    assert read_sets, 'No read files given'
#    
#    #Too clever: Might be remote file
#    #for filename in reference_filenames:
#    #    assert os.path.exists(filename.split('~~')[0]), '"%s" does not exist' % filename
#    #for item in read_sets:
#    #    for filename in item[0]:
#    #        assert os.path.exists(filename.split('~~')[0]), '"%s" does not exist' % filename
#
#    cutoff = '55%' #Default changed in SHRiMP 2.0.2
#    if '-h' in shrimp_options:
#        cutoff = shrimp_options[ shrimp_options.index('-h')+1 ]
#    
#    #Create working directory
#    
#    workspace = io.Workspace(working_dir)
#    
#    workspace.update_param( 
#        shrimp_options = original_args, 
#        shrimp_cutoff = cutoff
#    )
#    
#    #Make copy of reference sequences
#    
#    reference_filename = io.abspath(working_dir,'reference.fa')
#    reference_file = open(reference_filename,'wb')
#
#    reference_genbank_filename = io.abspath(working_dir,'reference.gbk')
#    reference_genbank_file = open(reference_genbank_filename,'wb')
#    any_genbank = [ False ]
#    
#    def genbank_callback(record):
#        """ Make a copy of any genbank files passed in. """
#        from Bio import SeqIO
#        
#        SeqIO.write([record], reference_genbank_file, 'genbank')
#        
#        f = open(os.path.join(
#            working_dir,
#            grace.filesystem_friendly_name(record.id) + '.gbk'
#        ), 'wb')
#        SeqIO.write([record], f, 'genbank')
#        f.close()
#        
#        any_genbank[0] = True
#    
#    for filename in reference_filenames:
#        for name, sequence in io.read_sequences(filename, genbank_callback=genbank_callback):
#            #Don't retain any comment
#            name = name.split()[0]
#            io.write_fasta(reference_file, name, sequence.upper())
#            
#            f = open(os.path.join(
#                working_dir,
#                grace.filesystem_friendly_name(name) + '.fa'
#            ), 'wb')
#            io.write_fasta(f, name, sequence.upper())
#            f.close()
#            
#    
#    reference_file.close()
#    reference_genbank_file.close()
#    if not any_genbank[0]:
#        os.unlink(reference_genbank_filename)
#    
#    # Create an index of the reference sequences
#    io.execute([
#        'samtools', 'faidx', reference_filename
#    ])
#    
#    #Run shrimp
#
#    bam_filename = io.abspath(working_dir, 'alignments.bam')
#    bam_prefix = io.abspath(working_dir, 'alignments')
#    bam_sorted_prefix = io.abspath(working_dir, 'alignments_sorted')
#    
#    temp_filename = io.abspath(working_dir, 'temp.bam')
#    
#    log_filename = io.abspath(working_dir, 'shrimp_log.txt')
#    log_file = open(log_filename, 'wb')
#
#    sam_eater = sam.Bam_writer(temp_filename)
#    
#    if use_cs:
#        program = 'gmapper-cs'
#    else:
#        program = 'gmapper-ls'
#    
#    sam_header_sent = [False]
#    n_seen = [0]
#    
#    def eat(process):
#        for line in process.stdout:
#            if line.startswith('@'):
#                if sam_header_sent[0]: continue
#            else:
#                n_seen[0] += 1
#                if n_seen[0] % 100000 == 0:
#                    grace.status('%s alignments produced' % grace.pretty_number(n_seen[0]))
#            sam_eater.write_raw(line)
#            
#        assert process.wait() == 0, 'shrimp failed'
#        sam_header_sent[0] = True
#
#    def remove_pair_options(options):
#        for flag in ['-p','-I']:
#            while flag in options:
#                pos = options.index(flag)
#                options = options[:pos] + options[pos+2:]
#        for flag in ['--half-paired']:
#            while flag in options:
#                pos = options.index(flag)
#                options = options[:pos] + options[pos+1:]
#        return options
#    
#    if '--qv-offset' not in shrimp_options:
#        guesses = [ ]
#        for filenames, is_paired in read_sets:
#            for filename in filenames:
#                guesses.append(io.guess_quality_offset(filename))
#        assert len(set(guesses)) == 1, 'Conflicting quality offset guesses, please specify --qv-offset manually.'
#        default_options['--qv-offset'] = str(guesses[0])
#            
#    for filenames, is_paired in read_sets:
#        options = shrimp_options[:]
#           
#        has_qualities = all(
#            len( io.read_sequences(filename, qualities=True).next() ) == 3  #A little ugly
#            for filename in filenames
#        )
#        if has_qualities:
#            options.append( '--fastq' )
#        #    temp_read_filename = io.abspath(working_dir, 'temp.fa')
#        #else:
#        #    temp_read_filename = io.abspath(working_dir, 'temp.fq')
#        
#        #try:
#        
#        #if len(filenames) == 1: # gmapper can cope with gzipped     and filenames[0].endswith('.fa') or filenames[0].endswith('.fq'):
#        #    actual_read_filename = filenames[0]
#        #else:
#        #    actual_read_filename = temp_read_filename
#        #    grace.status('Copying reads')
#        #    f = open(temp_read_filename, 'wb')
#        #    if has_qualities:
#        #        for reads in itertools.izip(*[ io.read_sequences(filename, qualities=True) for filename in filenames ]):
#        #            for name, seq, qual in reads:
#        #                io.write_fastq(f, name, seq, qual)
#        #    else:
#        #        for reads in itertools.izip(*[ io.read_sequences(filename) for filename in filenames ]):
#        #            for name, seq in reads:
#        #                io.write_fasta(f, name, seq)
#        #    f.close()
#        #    grace.status('')
#        
#        if len(filenames) == 1:
#            reads_parameters = [ filenames[0] ]
#        else:
#            reads_parameters = [ '-1', filenames[0], '-2', filenames[1] ]
#        
#        for flag in default_options:
#            if flag not in options:
#                options.append(flag)
#                if default_options[flag] is not None:
#                    options.append(default_options[flag])
#        
#        if not is_paired:
#           options = remove_pair_options(options)
#        
#        grace.status('')
#        
#        full_param = [ program ] + options + reads_parameters + [ reference_filename ]
#        
#        print >> sys.stderr, 'Running', ' '.join(full_param)
#        
#        p = io.run(full_param,
#                stdout=subprocess.PIPE,
#                stderr=log_file)
#        eat(p)
#                    
#        #finally:
#        #    if os.path.exists(temp_read_filename):
#        #        os.unlink(temp_read_filename)
#
#    log_file.close()
#
#    sam_eater.close()
#    
#    grace.status('Sort')
#
#    io.execute([
#        'samtools', 'sort', '-n', temp_filename, bam_prefix
#    ])
#    
#    os.unlink(temp_filename)
#    
#    grace.status('')
#    
#    ##TODO: make this optional
#    #sort_and_index(bam_filename, bam_sorted_prefix)
#
