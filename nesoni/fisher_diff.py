
"""

Changes:

0.2  - Added --titleX options.

then incorporated into nesoni

"""

import sys, re, itertools, os, math

from nesoni import grace, io, config, working_directory


def read_file(filename):
    f = open(filename,'rb')
    f.readline()
    for line in f:
        parts = line.rstrip('\n').split('\t')
        assert len(parts) >= 8, 'Old style evidence file, re-run samconsensus'
        yield int(parts[0]), parts[1], parts[2], parts[3], parts[6], parts[7]



LOG_FAC_CACHE = [ 0.0 ]
def log_fac(n):
    n = int(n)
    try:
        return LOG_FAC_CACHE[n]
    except IndexError:
        while len(LOG_FAC_CACHE) <= n:
            LOG_FAC_CACHE.append(LOG_FAC_CACHE[-1]+math.log(len(LOG_FAC_CACHE)))
        return LOG_FAC_CACHE[n]  

SUM_ERROR_MARGIN = 0.999999

class Cutoff_exceeded(Exception): pass


def fexact(matrix, significance_cutoff):
    numpy = grace.get_numpy()

    matrix = numpy.array(matrix)
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
    numpy = grace.get_numpy()

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
    
    matrix = numpy.zeros((n,2),int)
    for item in ev1:
        matrix[ options[item[0]], 0 ] = item[1]
    for item in ev2:
        matrix[ options[item[0]], 1 ] = item[1]
    
    #s = fexact(matrix)
    #if s == 0:
    #    print matrix, s
    #return fexact(matrix)
    
    key = (tuple(int(i) for i in matrix.ravel()),cutoff)
    if key not in SIG_CACHE:
        SIG_CACHE[key] = fexact(matrix,cutoff)    
    return SIG_CACHE[key]

@config.help("""\
Compare results of two runs of nesoni consensus using Fisher's Exact Test for each site in the reference.
""")
@config.String_flag('title1', 'Title for first working directory in results table.')
@config.String_flag('title2', 'Title for second working directory in results table.')
@config.Float_flag('cutoff', 'Significance level cutoff.')
@config.Positional('working_dir1', 'First working directory for comparison.')
@config.Positional('working_dir2', 'Second working directory for comparison.')
class Fisher(config.Action_with_prefix):
    title1 = None
    title2 = None
    cutoff = 1e-5
    working_dir1 = None
    working_dir2 = None

    def run(self):
        title1 = self.title1
        title2 = self.title2
        
        working1 = working_directory.Working(self.working_dir1)
        working2 = working_directory.Working(self.working_dir2)
        
        cutoff = self.cutoff
        
        sequence_names = [ name 
                           for name, length 
                           in working1.get_reference().get_lengths() ]
        
        if title1 is None:
            title1 = working1.name
        if title2 is None:
            title2 = working2.name
            
        n = 1
        while significance([('A',n)],[('T',n)],1.0) > cutoff:
            n += 1

        f = open(self.prefix + '.txt','wb')        
        print >> f, '%g\tsignificance cutoff' % cutoff
        print >> f, '%d\tdepth required to call substitution (greater if there are errors in the reads)' % n
            
        print >> f, 'Sequence\tPosition in reference\tChange type\tReference\t%s\t%s\tp-value (no correction for multiple testing)\t%s\t%s' % (title1, title2, title1, title2)
    
        for sequence_name in sequence_names:
            filename1 = working1/(grace.filesystem_friendly_name(sequence_name) + '-evidence.txt')
            filename2 = working2/(grace.filesystem_friendly_name(sequence_name) + '-evidence.txt')
        
            for (pos1, ins1, sub1, ref1, conins1, consub1), (pos2, ins2, sub2, ref2, conins2, consub2) in itertools.izip(read_file(filename1), read_file(filename2)):
                assert pos1 == pos2 and ref1 == ref2
            
                if pos1 % 1000 == 0:
                    grace.status('Testing %s %d' % (sequence_name, pos1))
            
                dec_ins1 = io.decode_evidence(ins1)
                dec_ins2 = io.decode_evidence(ins2)
                if dec_ins1 and dec_ins2:
                    sig = significance(io.decode_evidence(ins1), io.decode_evidence(ins2), cutoff)    
                    if sig is not None and sig <= cutoff:
                        print >> f, '%s\t%d\t%s\t\t%s\t%s\t%g\t%s\t%s' % (sequence_name, pos1, 'insertion-before', ins1, ins2, sig, conins1, conins2)
                        f.flush()
            
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
                        print >> f, '%s\t%d\t%s\t%s\t%s\t%s\t%g\t%s\t%s' % (sequence_name, pos1, what, ref1, sub1, sub2, sig, consub1, consub2)
                        f.flush()
        
        f.close()
        
        grace.status('')
        return 0
