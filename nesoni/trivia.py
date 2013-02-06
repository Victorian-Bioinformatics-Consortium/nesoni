
import nesoni
from nesoni import grace, io, bio, config, annotation, selection

from nesoni.third_party import vcf

import sys, os, math, random

@config.help("""\
Test that nesoni is able to run correctly.
""")
class Test(config.Action):
    def hello(self):
        print 'Hello, world.'
        return 'Message from future.'

    def run(self):
        f = nesoni.future(self.hello)
        print f()

@config.help("""\
Output the contents of sequence files in FASTA format.
""")
@config.Main_section('filenames', 'Sequence files.', empty_is_ok=False)
@config.Section('only', 'Only output these sequences.')
class As_fasta(config.Action):
    filenames = [ ]
    only = [ ]

    def run(self):
        for filename in self.filenames:
            for name, seq in io.read_sequences(filename):
                if self.only and name.split()[0] not in self.only: continue
                io.write_fasta(sys.stdout, name, seq) 


@config.help("""\
Randomly sample from a sequence file.
""")
@config.Int_flag('n', 'How many sequences to output.')
@config.Main_section('filenames', 'Input sequence files.', empty_is_ok=False)
class Sample(config.Action):
    n = 10000
    filenames = [ ]
    
    def run(self):
        seqs = [ ]
        seen = 0
        for filename in self.filenames:
            for seq in io.read_sequences(filename, qualities=True):
                seen += 1
                if seen % 100000 == 0: grace.status('Scanned '+grace.pretty_number(seen))
                
                if len(seqs) < self.n:
                    seqs.append(seq)
                elif self.n <= random.random() * seen:
                    seqs[ random.randrange(self.n) ] = seq
        grace.status('')
        
        print >> sys.stderr, 'Sampled', grace.pretty_number(len(seqs)), 'of', grace.pretty_number(seen), 'sequences'
    
        if not seqs: return
        qualities = len(seqs[0])
        if qualities:
            for name, seq, qual in seqs:
                io.write_fastq(sys.stdout, name, seq, qual)
        else:
            for name, seq in seqs:
                io.write_fastq(sys.stdout, name, seq)


@config.help("""\
Show basic statistics of a sequence or annotation or VCF file.
""")
@config.Main_section('filenames', 'Sequence files.', empty_is_ok=False)
class Stats(config.Action_with_optional_output):
    filenames = [ ]
    
    def run(self):
        f = self.begin_output()
    
        for filename in self.filenames:
            info = io.get_file_info(filename)
            
            any = False
            
            name = os.path.splitext(os.path.split(filename)[1])[0]
            
            if info.matches('sequences'):
                total = 0
                total_length = 0
                for seq in io.read_sequences(filename, qualities=True):
                    total += 1
                    total_length += len(seq[1])
                print >> f, grace.datum(name, 'sequences', total)
                print >> f, grace.datum(name, 'total bases', total_length)
                if total:
                    print >> f, grace.datum(name, 'average length', float(total_length)/total)
                print >> f
                any = True
            
            if info.matches('annotations'):
                total = 0
                counts = { }
                for item in annotation.read_annotations(filename):
                    total += 1
                    counts[item.type] = counts.get(item.type,0)+1
                                
                print >> f, grace.datum(name, 'features', total)
                for key in sorted(counts):
                    print >> f, grace.datum(name, key + ' features', counts[key])
                print >> f
                any = True
            
            if info.matches('type-vcf'):
                reader_f = io.open_possibly_compressed_file(filename)
                reader = vcf.Reader(reader_f)
                n = 0
                for item in reader:
                    n += 1
                print >> f, grace.datum(name, 'variants', n)
                any = True
            
            if not any:
                raise grace.Error('Don\'t know what to do with ' + filename)

        self.end_output(f)

@config.help("""\
Output an annotation file in GFF format.
""")
@config.String_flag('output', 'Output file, defaults to STDOUT.')
@config.String_flag('select', 'What types of annotation to keep (selection expression).')
@config.Main_section('filenames', 'Annotation files.', empty_is_ok=False)
class As_gff(config.Action):
    filenames = [ ]
    select = 'all'
    
    output = None
    
    def ident(self):
        return config.Action.ident(self) + ('--' + self.output if self.output else '')
 
    def run(self):
        if self.output is not None:
           out_file = open(self.output,'wb')
        else:
           out_file = sys.stdout
    
        annotation.write_gff3_header(out_file)
        
        for filename in self.filenames:
            for item in annotation.read_annotations(filename):
                if not selection.matches(self.select, [item.type]): continue
                
                if 'ID' not in item.attr and 'locus_tag' in item.attr:
                    item.attr['ID'] = item.attr['locus_tag']
                    
                if 'color' not in item.attr:
                    if item.type == 'CDS':
                        item.attr['color'] = '#008800'
                    if item.type == 'rRNA':
                        item.attr['color'] = '#bb0000'
                    if item.type == 'tRNA':
                        item.attr['color'] = '#bb00bb'
                    if item.type == 'misc_feature':
                        item.attr['color'] = '#8888ff'

                print >> out_file, item.as_gff()
        
        if self.output is not None:
            out_file.close()


        

def plot(args):
    log_it, args = grace.get_option_value(args, '--log', grace.as_bool, False)

    grace.expect_no_further_options(args)
    
    import numpy, pylab
    
    pylab.rcParams['axes.formatter.limits'] = [ -20, 20 ]
    
    pylab.figure(figsize=(10,4))
    
    maximum = 0
    for filename in args:
        parts = filename.split('~~', 1)
        data = [ ]
        f = open(parts[0],'rb')
        for line in f:
            data.append(float(line.strip()))
        f.close()
        
        data = numpy.array(data)
        
        maximum = max(maximum,numpy.maximum.reduce(data))
        
        #if log_it:
        #    data = numpy.log(data + 1.0) / numpy.log(2.0)
        
        if log_it:
            pylab.semilogy( numpy.arange(1,len(data)+1), data, label=parts[-1] )
        else:
            pylab.plot( numpy.arange(1,len(data)+1), data, label=parts[-1] )
    
    if len(args) > 1:
        pylab.legend()
        
        if log_it:
            pylab.ylim( (1,maximum**1.2) )
        else:
            pylab.ylim( (0,maximum*1.2) )
    
    pylab.show()




def normalize_file(f_in, f_out, trim):
    headers = [ ]
    data = [ ]
    for line in f_in:
        if line.startswith('#'):
            headers.append(line)
        else:
            data.append( [ int(item) for item in line.strip().split() ] )

    is_multiplot = (len(headers) > 0)
    all_depths = [ ]
    for item in data:
        if is_multiplot:
            all_depths.extend(item[1:])
        else:
            all_depths.append(item[0])

    all_depths.sort()
    #median = all_depths[ len(all_depths) // 2]
    mean = float(sum(all_depths)) / len(all_depths)
    
    #import pylab
    #pylab.plot(all_depths)
    #pylab.show()
    
    #print >> sys.stderr
    #prev = mean
    #for cutoff in [ 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999 ]:
    #   trim_len = int(len(all_depths) * cutoff)
    #   trim_mean = float(sum(all_depths[:trim_len])) / trim_len
    #   print >> sys.stderr, cutoff, trim_mean, trim_mean / prev
    #   prev = trim_mean
    
    trim_len = int(len(all_depths) * (1.0-trim))
    trim_mean = float(sum(all_depths[:trim_len])) / trim_len
    
    if trim_mean > 0.0:
        print >> sys.stderr, ' Mean: %.2f  Trimmed mean: %.2f  Ratio: %.2f:1' % (mean, trim_mean, mean/trim_mean)
    else:
        print >> sys.stderr, ' Mean: %.2f  Trimmed mean: %.2f  Ratio: %.2f:1  Not normalized!' % (mean, trim_mean)
    
    norm = trim_mean
    if norm == 0: norm = 1.0 #Hmm

    for line in headers:
        f_out.write(line)
    
    for item in data:
        if is_multiplot:
            print >> f_out, item[0], ' '.join([ str(depth / norm) for depth in item[1:] ])
        else:
            print >> f_out, item[0] / norm


def read_userplot(filename):
    grace.status('Load '+filename)
    
    f = open(filename,'rb')
    lines = f.readlines()
    f.close()
    
    n = 0
    while n < len(lines) and lines[n].startswith('#'):
        n += 1
    headers = lines[:n]

    is_multiplot = (n > 0)

    data = [ tuple([ int(item) for item in lines[i].strip().split() ])
             for i in xrange(n,len(lines)) ]

    grace.status('')
    return headers, is_multiplot, data

def read_unstranded_userplot(filename):
    headers, is_multiplot, data = read_userplot(filename)
    if is_multiplot:
        return [ item[1]+item[2] for item in data ]
    else:
        return [ item[0] for item in data ]

def write_unstranded_userplot(filename, array):
    f = open(filename,'wb')
    for x in array:
        f.write( '%f\n' % float(x) )
    f.close()

def write_normalized_userplot((headers, is_multiplot, data), factor, filename):
    grace.status('Write '+filename)
    f = open(filename, 'wb')
    for line in headers:
        f.write(line)

    for item in data:
        if is_multiplot:
            print >> f, item[0], ' '.join([ str(depth * factor) for depth in item[1:] ])
        else:
            print >> f, item[0] * factor
    f.close()
    grace.status('')

def normalize_files(dirnames, prefix, min_depth):
    contents = [ read_userplot(os.path.join(item,prefix+'-depth.userplot')) for item in dirnames ]
    data = [ item[2] for item in contents ]

    if contents[0][1]:
        j_range = xrange(1,len(contents[2][0]))
    else:
        j_range = [0]
    
    totals = [ 0.0 ] * len(data)
    n = 0

    for i in xrange(len(contents[0][2])):
        for j in j_range:
            depths = [ item[i][j] for item in data ]
            good = True
            for item in depths: 
                if item < min_depth: 
                    good = False 
                    break
            if not good: continue
            
            for k, item in enumerate(depths):
                totals[k] += math.log(item)
            n += 1
    
    print prefix, 'sites used:', grace.pretty_number(n)
    
    if n == 0:
        print 'Can\'t normalize, skipping.'
        print
        return
    
    avg_total = sum( totals ) / len(data)
    norm_divisors = [ math.exp( (item - avg_total)/n ) for item in totals ]

    print 'Relative abundances:'
    for i in xrange(len(dirnames)):
        print '  %.3f %s' % (norm_divisors[i], dirnames[i])
    print
    
    #for i, item in enumerate(contents):
    #    write_normalized_userplot(item, 1.0/norm_divisors[i], os.path.join(dirnames[i],prefix+'-norm.userplot'))
    for i, dirname in enumerate(dirnames):
        for filename in os.listdir(dirname):
            if not filename.startswith(prefix): continue
            if not filename.endswith('.userplot'): continue
            if filename.endswith('-norm.userplot'): continue
            fullname = os.path.join(dirname, filename)
            full_outname = fullname[:-9] + '-norm.userplot'
            write_normalized_userplot(
                read_userplot(fullname),
                1.0/norm_divisors[i],
                full_outname)



NORMALIZE_HELP = """\

Usage:

    nesoni normalize: [options] working_dir1 working_dir2 [working_dir3 ...]

Creates xxx-norm.userplot files for each xxx.userplot file in the working directories
supplied.

Options:

    --min-depth NNN     - When calculating normalization factors, ignore
                          locations with less than this depth in at least one
                          of the working directories.
                          Default: 5

"""

def normalize(args):
    min_depth, args = grace.get_option_value(args, '--min-depth', int, 5) 
    grace.expect_no_further_options(args) 

    if len(args) < 2:
        print NORMALIZE_HELP
        raise grace.Help_shown()

    dirnames = args

    filenames = [ ]
    for dirname in dirnames:
        assert os.path.isdir(dirname), dirname + ' is not a directory'

        filenames.append(sorted(
            item for item in os.listdir(dirname)
            #if item.endswith('.userplot') and not item.endswith('-norm.userplot')
            if item.endswith('-depth.userplot')
            and not item.endswith('-ambiguous-depth.userplot')
            and not item.endswith('-pairspan-depth.userplot')
        ))
    
    for i in xrange(1,len(dirnames)):
        if filenames[i] != filenames[0]:
            raise grace.Error('Userplots in %s differ from those in %s' % (dirnames[i], dirnames[0]))
    filenames = filenames[0]

    for filename in filenames:
        normalize_files(dirnames, filename[:-15], min_depth)
        

# Triangular sliding window
def windower(thing, max_radius):
    import numpy

    thing_pad = numpy.concatenate((
        thing[-max_radius:], thing, thing[:max_radius]
        ))
    thing_sum = numpy.cumsum(numpy.cumsum(thing_pad))
    
    return (len(thing), thing_sum, max_radius) 

def use_windower(windower, window_radius):
    len_thing, thing_sum, max_radius = windower
    return (-2*thing_sum[max_radius:len_thing+max_radius]
            +thing_sum[max_radius+window_radius:][:len_thing] 
            +thing_sum[max_radius-window_radius:][:len_thing]) / float(window_radius**2)

def expected_depth(name, seq, depths, ambig_depths, radius=2):
    import numpy
    from numpy import linalg

    med = numpy.median(depths)
    sane = numpy.arange(len(depths))[ (depths > med*0.5) & (depths < med*2.0) & (depths*2.0 >= ambig_depths)]
    #print 'median', med, 'using', len(sane)
    
    if sum(sane) < 100:
        warn('Skipping depth correction on ' + name)
        return numpy.array( [numpy.average(depths)] * len(depths) ) 
    
    buckets = { }
    
    #radius = 2 # examine 5-mers
    n = radius*2+1
    
    sseqq = seq[len(seq)-radius:] + seq + seq[:radius]
    for i in sane:
        s = sseqq[i:i+n]
        if s not in buckets: buckets[s] = [ ]
        buckets[s].append( depths[i] )
    
    # Pool with reverse complement
    new_buckets = { }
    for kmer in buckets:
        rc = bio.reverse_complement(kmer)
        pool = buckets[kmer] + buckets.get(rc,[])
        new_buckets[kmer] = pool 
    buckets = new_buckets
    
    buckets_individual = buckets.copy()    
    for key in buckets:
        buckets[key] = numpy.average(buckets[key])
        
    
    avg_depth = numpy.average(depths[sane])
    listing = [ (key, numpy.log(value)-numpy.log(avg_depth) ) for key,value in buckets.items()
                 if key <= bio.reverse_complement(key) ]
    listing.sort(key=lambda x: abs(x[1]), reverse=True)
    print 'Top k-mer log2 fold change'
    for key,value in listing[:10]:
        print key, '% .2f' % (value / numpy.log(2.0)), '(%d)' % len(buckets_individual[key])
    print
    
    prediction = numpy.zeros(len(seq), 'float')
    for i in xrange(len(seq)):
        s = sseqq[i:i+n]
        prediction[i] = buckets.get(s,0.0)
    
    
    # selection of radii from 8 to 4096
    # TODO: make this configurable, or perhaps just larger
    radii = [ int(2**(0.5*i)) for i in xrange(3*2,12*2+1) ]
    
    prediction_windower = windower(numpy.log(prediction), radii[-1]) 
    
    a = numpy.arange(len(seq)) / float(len(seq))
    predictors = numpy.transpose(
    [   numpy.ones(len(seq), 'float'),
        numpy.cos(a * (2.0*numpy.pi)), 
        numpy.sin(a * (2.0*numpy.pi)),
    ] + [
        use_windower(prediction_windower, radius)
        for radius in radii
    ]
    )
    
    x = linalg.lstsq(predictors[sane], numpy.log(depths[sane]))[0]
    #print x
    prediction = numpy.sum(predictors * x[None,:], 1)
    
    change = numpy.median( numpy.abs( numpy.log(depths) - prediction ) )
    print 'Median log2 fold error:', change / numpy.log(2.0)
    print
    print
    
    return numpy.exp( prediction )

def debias(args):
    import numpy

    radius, args = grace.get_option_value(args, '--radius', int, 2) 

    dirs = args
    
    for dir_name in dirs:
        for name, seq in io.read_sequences(os.path.join(dir_name,'reference.fa')):
            for suffix, ambig_suffix in [
                ('-depth', '-ambiguous-depth'),
                ('-pairspan-depth', '-ambiguous-pairspan-depth'),
            ]:
                root = grace.filesystem_friendly_name(name)
                full_name = os.path.join(dir_name, root + suffix + '.userplot')
                full_ambig_name = os.path.join(dir_name, root + ambig_suffix + '.userplot')
                if not os.path.exists(full_name): continue
                if not os.path.exists(full_ambig_name): continue
                
                output_suffix = '-%d.userplot' % radius 

                print dir_name, root, output_suffix
                
                depths = numpy.array( read_unstranded_userplot(full_name) )
                ambig_depths = numpy.array( read_unstranded_userplot(full_ambig_name) )
                expect = expected_depth(root, seq, depths, ambig_depths, radius)
                
                write_unstranded_userplot(
                    os.path.join(dir_name, root + suffix + '-expected' + output_suffix),
                    expect) 
                
                corrected = depths / expect * numpy.median(expect)
                corrected[expect <= 5.0] = 0.0
                write_unstranded_userplot(
                    os.path.join(dir_name, root + suffix + '-corrected' + output_suffix),
                    corrected)                 
                
                ambig_corrected = ambig_depths / expect * numpy.median(expect)
                ambig_corrected[expect <= 0.0] = 0.0
                write_unstranded_userplot(
                    os.path.join(dir_name, root + ambig_suffix + '-corrected' + output_suffix),
                    ambig_corrected)                 
                
                



