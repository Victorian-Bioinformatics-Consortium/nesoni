
"""

Changes:

0.2  - Added --titleX options.

then incorporated into nesoni

"""

import sys, re, numpy, itertools, os

from nesoni import grace, io

USAGE = """\

Usage:

  nesoni fisher: [options] workingdir1 workingdir2 significance_cutoff

Options:

  --title1 "xyz"   - Column title for first evidence file in output
  --title2 "xyz"   - Column title for second evidence file in output

significance_cutoff can be in scientific notation, eg 1e-5
"""


def read_file(filename):
    f = open(filename,'rb')
    f.readline()
    for line in f:
        parts = line.rstrip('\n').split('\t')
        assert len(parts) >= 8, 'Old style evidence file, re-run samconsensus'
        yield int(parts[0]), parts[1], parts[2], parts[3], parts[6], parts[7]



LOG_FAC_CACHE = [ 0.0 ]
def log_fac(n):
    try:
        return LOG_FAC_CACHE[n]
    except IndexError:
        while len(LOG_FAC_CACHE) <= n:
            LOG_FAC_CACHE.append(LOG_FAC_CACHE[-1]+numpy.log(len(LOG_FAC_CACHE)))
        return LOG_FAC_CACHE[n]  

SUM_ERROR_MARGIN = 0.999999

class Cutoff_exceeded(Exception): pass


def fexact(matrix, significance_cutoff):
    matrix = numpy.asarray(matrix)
    n_row, n_col = matrix.shape
    row_sum = numpy.sum(matrix,1)
    col_sum = numpy.sum(matrix,0)
    n = numpy.sum(row_sum)

    cutoff = sum([ log_fac(item) for item in matrix.ravel() ]) * SUM_ERROR_MARGIN

    const_part = sum([ log_fac(item) for item in row_sum ]) + sum([ log_fac(item) for item in col_sum ]) - log_fac(n)
    
    significance = [ 0.0 ]
    
    row_remainders = row_sum.copy()
    def generate(row, col, col_remainder, total):
        cell_min = max(0, col_remainder - numpy.sum(row_remainders[row+1:]))
        cell_max = min(col_remainder, row_remainders[row])
        
        #next_row_remainders = row_remainders.copy()
        old_row_remainder = row_remainders[row]
        for i in xrange(cell_min,cell_max+1):
            next_total = total + log_fac(i)            
            row_remainders[row] = old_row_remainder - i
            if row+1 >= n_row:
                if col+1 >= n_col:
                    if next_total >= cutoff:
                        significance[0] += numpy.exp( const_part-next_total )
                        if significance[0] > significance_cutoff:
                            raise Cutoff_exceeded
                else:
                    generate(0, col+1, col_sum[col+1], next_total)
            else:
                generate(row+1, col, col_remainder-i, next_total)
        row_remainders[row] = old_row_remainder

    try:
        generate(0,0,col_sum[0],0.0)
    except Cutoff_exceeded:
        return None
    
    return significance[0]

SIG_CACHE = { }
def significance(ev1, ev2, cutoff):
    # Avoid duplicates in SIG_CACHE
    if ev2 < ev1:
        ev1,ev2 = ev2,ev1

    options = { }
    for item in ev1:
        options[item[0]] = len(options)
    for item in ev2:
        if item[0] not in options:
            options[item[0]] = len(options)
    
    n = len(options)
    
    matrix = numpy.zeros((n,2),'int')
    for item in ev1:
        matrix[ options[item[0]], 0 ] = item[1]
    for item in ev2:
        matrix[ options[item[0]], 1 ] = item[1]
    
    #s = fexact(matrix)
    #if s == 0:
    #    print matrix, s
    #return fexact(matrix)
    
    key = (tuple(matrix.ravel()),cutoff)
    if key not in SIG_CACHE:
        SIG_CACHE[key] = fexact(matrix,cutoff)    
    return SIG_CACHE[key]


def main(args):
    title1, args = grace.get_option_value(args, '--title1', str, None)
    title2, args = grace.get_option_value(args, '--title2', str, None)
    grace.expect_no_further_options(args)
        
    if len(args) != 3:
        print >> sys.stderr, USAGE
        return 1
    
    working_dir1 = args[0]
    working_dir2 = args[1]
    cutoff = float(args[2])
    
    sequence_names = [ name 
                       for name, sequence 
                       in io.read_sequences(os.path.join(working_dir1, 'reference.fa')) ]
    
    if title1 is None:
        title1 = working_dir1
    if title2 is None:
        title2 = working_dir2
        
    n = 1
    while significance([('A',n)],[('T',n)],1.0) > cutoff:
        n += 1
    
    print '%g\tsignificance cutoff' % cutoff
    print '%d\tdepth required to call substitution (greater if there are errors in the reads)' % n
        
    print 'Sequence\tPosition in reference\tChange type\tReference\t%s\t%s\tp-value (no correction for multiple testing)\t%s\t%s' % (title1, title2, title1, title2)

    for sequence_name in sequence_names:
        filename1 = os.path.join(working_dir1, grace.filesystem_friendly_name(sequence_name) + '-evidence.txt')
        filename2 = os.path.join(working_dir2, grace.filesystem_friendly_name(sequence_name) + '-evidence.txt')
    
        for (pos1, ins1, sub1, ref1, conins1, consub1), (pos2, ins2, sub2, ref2, conins2, consub2) in itertools.izip(read_file(filename1), read_file(filename2)):
            assert pos1 == pos2 and ref1 == ref2
        
            if pos1 % 1000 == 0:
                grace.status('Testing %s %d' % (sequence_name, pos1))
        
            dec_ins1 = io.decode_evidence(ins1)
            dec_ins2 = io.decode_evidence(ins2)
            if dec_ins1 and dec_ins2:
                sig = significance(io.decode_evidence(ins1), io.decode_evidence(ins2), cutoff)    
                if sig is not None and sig <= cutoff:
                    grace.status('')
                    print '%s\t%d\t%s\t\t%s\t%s\t%g\t%s\t%s' % (sequence_name, pos1, 'insertion-before', ins1, ins2, sig, conins1, conins2)
        
            dec_sub1 = io.decode_evidence(sub1)
            dec_sub2 = io.decode_evidence(sub2)
            if dec_sub1 and dec_sub2:
                sig = significance(dec_sub1, dec_sub2, cutoff)        
                if sig is not None and sig <= cutoff:
                    if dec_sub1[0][0] == '-' or dec_sub2[0][0] == '-':
                        what = 'deletion'
                    elif dec_sub1[0][0] != dec_sub2[0][0]:
                        what = 'substitution'
                    else:
                        what = 'different mix'
                    grace.status('')
                    print '%s\t%d\t%s\t%s\t%s\t%s\t%g\t%s\t%s' % (sequence_name, pos1, what, ref1, sub1, sub2, sig, consub1, consub2)
    
    grace.status('')
    return 0
