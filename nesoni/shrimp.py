"""

Run shrimp, set up directory for nesoni consensus

"""

import sys, os, pprint, gzip, itertools

from nesoni import io, grace

USAGE = """\
Usage:

    nesoni shrimp: output_directory [options] \\
        reference.fa [...] \\
        [reads: single.fq | single.fa | interleaved.fq | interleaved.fa [...]] \\
        [pairs: left.fq right.fq | left.fa right.fa] \\
        [shrimp-options: ...options to pass directly to rmapper-xx... ]
  
Paired end reads should either be given as two files in the "pairs" section,
or as a single interleaved file in the "reads" section. 

Files may be gzipped. Reference files may be in FASTA or GENBANK format.

Options:
  
  --threshold N%%  - Hit score threshold (-h in rmapper-ls/cs)
                    Can either be a score or a percentage
                    (base match = 10 points, base mismatch = -15 points)  
                    Default 68%%

  --solid         - Reads are from a SOLiD sequencer, ie colorspace
                    (Default is FASTA or FASTQ format)
  
  --verbose       - Show output from SHRiMP
  --cpus N        - How many SHRiMPs to run in parallel 
                    (default: the number of CPUs in your system, %d)
  --batch-size N  - How many bases worth of reads to use in each SHRiMP 
                    invocation (default: 5000000)
  
  --stride N      - Only use 1 read in N

Useful shrimp-options (for more options, type "rmapper-ls"):

  -n 1            - Perform alignment on a single hash-table hit (default 4)
                    (more sensitive, but slow)
  -o 5            - Output a maximum of 5 hits per read (default 100)
                    (useful for repetitive genomes, 
                     but hit pairing may not work as well)
  -ungapped       - Only find ungapped alignments 
                    (faster)

"""

def load_config(dir_name):
    filename = os.path.join(dir_name, 'config.txt')
    f = open(filename, 'rb')
    config = eval(f.read())
    f.close()
    return config

def iter_reads(config, qualities=False):
    if 'stride' not in config:
        raise grace.Error('Please re-run nesoni shrimp, output format has changed')

    stride = config['stride']
    for reads_filename_set in config['reads']:
        if config['solid']:
            reader = [ io.read_solid(filename) for filename in reads_filename_set ]
        else:
            reader = [ io.read_sequences(filename, qualities) for filename in reads_filename_set ]
        reader = itertools.izip(*reader)
        
        for i, items in enumerate(reader):
            if i % stride == 0: 
               for item in items: 
                   yield item

def iter_read_hits(dir_name, hit_filename='shrimp_hits.txt.gz', qualities=False):
    """ yields ( (read_name, read_seq), [ hits ] ) """
    
    config = load_config(dir_name)
    hit_file = gzip.open(os.path.join(dir_name,hit_filename), 'rb')
    line = hit_file.readline()
    
    for read in iter_reads(config, qualities):
        prefix = '>' + read[0].split()[0] + '\t'
        hits = [ ]
        while line.startswith(prefix):
            hits.append(line)
            line = hit_file.readline()
        yield (read, hits)
    
    assert line == '', 'shrimp_hits.txt.gz contains unexpected hits!' + line    

    hit_file.close()

def calc_threshold(length, threshold):
    if threshold >= 0:
        return threshold
    else:
        return length * 10 * -threshold

def main(args):
    grace.require_shrimp_1()

    n_cpus = grace.how_many_cpus()
        
    solid, args = grace.get_flag(args, '--solid')
    verbose, args = grace.get_flag(args, '--verbose')

    threshold, args = grace.get_option_value(args, '--threshold', str, '68%')
    
    stride, args = grace.get_option_value(args, '--stride', int, 1)
    max_shrimps, args = grace.get_option_value(args, '--cpus', int, n_cpus)
    batch_size, args = grace.get_option_value(args, '--batch-size', int, 5000000)
        
    input_reference_filenames = [ ]
    reads_filenames = [ ]
    
    shrimp_options = [ '-h', threshold ]
    if threshold.endswith('%'):
        threshold = -float(threshold[:-1])/100.0
    else:
        threshold = int(threshold)
    
    output_dir = [ ]  #As list so can write to from function. Gah.
    
    def front_command(args):
        grace.expect_no_further_options(args)
        
        if len(args) < 1:
            return
        
        output_dir.append(args[0])        
        input_reference_filenames.extend(
            [ os.path.abspath(filename) for filename in args[1:] ])
    def reads_command(args):
        grace.expect_no_further_options(args)
        reads_filenames.extend([ [ os.path.abspath(filename) ] for filename in args])
    def pairs_command(args):
        grace.expect_no_further_options(args)
        assert len(args) == 2, 'Expected exactly two files in "pairs"'
        reads_filenames.append([ os.path.abspath(filename) for filename in args ])
    def shrimp_options_command(args):
        shrimp_options.extend(args)
    
    grace.execute(args, {
        'reads': reads_command,
        '--reads': reads_command,
        'pairs': pairs_command,
        'shrimp-options': shrimp_options_command,
        '--shrimp-options': shrimp_options_command,
    }, front_command)
    
    
    if not output_dir:
        print >> sys.stderr, USAGE % n_cpus
        return 1
    
    output_dir = output_dir[0]
    
    assert input_reference_filenames, 'No reference files given'
    assert reads_filenames, 'No read files given'
    
    for filename in itertools.chain(input_reference_filenames, *reads_filenames):
        assert os.path.exists(filename), '%s does not exist' % filename

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    if solid:
        shrimp = 'rmapper-cs'
    else:
        shrimp = 'rmapper-ls'
    
    
    reference_filename = os.path.join(output_dir,'reference.fa')
    reference_file = open(reference_filename,'wb')
    total_reference_sequences = 0
    total_reference_bases = 0
    for input_reference_filename in input_reference_filenames:
        for name, sequence in io.read_sequences(input_reference_filename):
            #Don't retain any comment
            name = name.split()[0]
            io.write_fasta(reference_file, name, sequence)
            
            total_reference_sequences += 1
            total_reference_bases += len(sequence)
            
    reference_file.close()
    
    print '%s base%s in %s reference sequence%s' % (
        grace.pretty_number(total_reference_bases), 's' if total_reference_bases != 1 else '',
        grace.pretty_number(total_reference_sequences), 's' if total_reference_sequences != 1 else '')
    
    assert total_reference_bases, 'Reference sequence file is empty' 
    
    config = {
        'references' : input_reference_filenames,
        'reads' : reads_filenames,
        'stride' : stride,
        'solid': solid,
        'threshold': threshold,
    }
    config_file = open(os.path.join(output_dir, 'config.txt'), 'wb')
    pprint.pprint(config, config_file)
    config_file.close()
    
    output_filename = os.path.join(output_dir, 'shrimp_hits.txt.gz')
    output_file = gzip.open(output_filename, 'wb')
    
    unmapped_filename = os.path.join(output_dir, 'unmapped.fa.gz')
    unmapped_file = gzip.open(unmapped_filename, 'wb')
    
    dirty_filenames = set()
    dirty_filenames.add(output_filename)
    dirty_filenames.add(unmapped_filename)
    
    #warn_low_threshold = True
    
    try: #Cleanup temporary files
        
        N = [0]
        def do_shrimp(read_set):
            my_number = N[0]
            N[0] += 1
            
            tempname = os.path.join(output_dir,'temp%d-%d.fa' % (os.getpid(),my_number))
            tempname_out = os.path.join(output_dir,'temp%d-%d.txt' % (os.getpid(),my_number))
            
            dirty_filenames.add(tempname)
            dirty_filenames.add(tempname_out)
            
            f = open(tempname,'wb')
            for read_name, read_seq in read_set:
                print >> f, '>' + read_name
                print >> f, read_seq
            f.close()
        
            command = shrimp + ' ' + ' '.join(shrimp_options) + ' ' + \
                      tempname + ' ' + reference_filename + ' >' + tempname_out
            if not verbose:
                command += ' 2>/dev/null'
            #f = os.popen(command, 'r')
            child_pid = os.spawnl(os.P_NOWAIT,'/bin/sh','/bin/sh','-c',command)
            #print 'SHRiMP %d running' % my_number
            
            def finalize():
                exit_status = os.waitpid(child_pid, 0)[1]
                assert exit_status == 0, 'Shrimp indicated an error'
                
                hits = { } # read_name -> [ hit line ]
                
                f = open(tempname_out,'rb')
                for line in f:
                    if line.startswith('>'):
                        read_name = line.split(None,1)[0][1:]
                        if read_name not in hits:
                            hits[read_name] = [ ]
                        hits[read_name].append(line)
                f.close()
                                
                for read_name, read_seq in read_set:
                    if read_name in hits:
                        for hit in hits[read_name]:
                            output_file.write(hit)
                    else:
                        print >> unmapped_file, '>' + read_name
                        print >> unmapped_file, read_seq

                output_file.flush()
                unmapped_file.flush()
        
                os.unlink(tempname)
                dirty_filenames.remove(tempname)
                os.unlink(tempname_out)
                dirty_filenames.remove(tempname_out)
                #print 'SHRiMP %d finished' % my_number
            return finalize
        
        
        shrimps = [ ]
        
        reader = iter_reads(config)
        read_count = 0
        
        while True:
            read_set = [ ]
            read_set_bases = 0

            #Read name should not include comment cruft
            # - SHRIMP passes this through
            # - might stuff up identification of pairs
            
            for read_name, read_seq in reader:
                read_name = read_name.split()[0]                
                read_set.append((read_name, read_seq))
                read_set_bases += len(read_seq)
                
                #if warn_low_threshold and len(read_seq)*7 < threshold: #Require 70% exact match
                #    sys.stderr.write('\n*** WARNING: Short reads, consider reducing --threshold ***\n\n')                    
                #    warn_low_threshold = False
            
                read_count += 1
                if read_set_bases >= batch_size: break
                
            if not read_set: break
        
            if len(shrimps) >= max_shrimps:
                shrimps.pop(0)()
            shrimps.append( do_shrimp(read_set) )
            
            grace.status('SHRiMPing %s' % grace.pretty_number(read_count))
        
        while shrimps:
            grace.status('Waiting for SHRiMPs to finish %d ' % len(shrimps) )
            shrimps.pop(0)()
        
        grace.status('')
        
        output_file.close()
        dirty_filenames.remove(output_filename)
        unmapped_file.close()
        dirty_filenames.remove(unmapped_filename)
        
        return 0

    finally:
        for filename in dirty_filenames:
            if os.path.exists(filename):
                os.unlink(filename)

