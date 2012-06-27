
"""

Core-genome tool.

"""

import sys, os
from os import path

from nesoni import grace, io, trivia, working_directory

HELP = """\

Usage:

    nesoni core: output_dir working_dir1 [working_dir2 ...] 

Identify a core genome present in all of a set of strains.
(Or alternatively sequence unique to the reference.)

Options:

    --mincov NNN       - Minimum coverage to declare sequence is present
                         Default: 1

    --maxdiff NNN      - Remove holes up to this size in the core genome
                         (morphological close operation)
                         Default: 16

    --minsize NNN      - Minimum size core and non-core chunks to output
                         Default: 200

    --what core/unique - What to compute
                           core   - sequence present in all samples and reference
                                    based on *-depth.userplot files
                           unique - sequence only present in reference
                                    based on *-ambiguous-depth.userplot files
                         Default: core

"""

def as_core_or_unique(string):
    string = string.lower()
    assert string in ('core', 'unique')
    return string

def main(args):
    mincov, args = grace.get_option_value(args, '--mincov', int, 1) 
    maxdiff, args = grace.get_option_value(args, '--maxdiff', int, 16) 
    minsize, args = grace.get_option_value(args, '--minsize', int, 200)
    what, args = grace.get_option_value(args, '--what', as_core_or_unique, 'core')    
    is_core = (what == 'core') 

    grace.expect_no_further_options(args)
    
    if len(args) < 2:
        print >> sys.stderr, HELP
        raise grace.Help_shown()
    
    output_dir, working_dirs = args[0], args[1:]
    
    #assert not path.exists(path.join(output_dir, 'reference.fa')), \
    assert not path.exists(nesoni.working_directory.Working(working_dirs[0]).get_reference().reference_fasta_filename()), \
            'Output directory not given'
    
    if not path.exists(output_dir):
        os.mkdir(output_dir)
    
    #for name, seq in io.read_sequences(path.join(working_dirs[0],'reference.fa')):
    for name, seq in io.read_sequences(nesoni.working_directory.Working(working_dirs[0]).get_reference().reference_fasta_filename())):
        print name
        friendly_name = grace.filesystem_friendly_name(name)
        
        good = [ True ] * len(seq)
        
        for working_dir in working_dirs:
            if is_core:
               suffix = '-depth.userplot'
            else:
               suffix = '-ambiguous-depth.userplot'
            data = trivia.read_unstranded_userplot(
                os.path.join(working_dir, friendly_name+suffix)
            )
            assert len(seq) == len(data)
            for i in xrange(len(seq)):
               if good[i]:
                   if is_core:
                       good[i] = data[i] >= mincov
                   else:
                       good[i] = data[i] < mincov

        #Close holes
        start = -maxdiff-1
        n_holes = 0
        for i in xrange(len(seq)):
            if good[i]:
                 if 0 < i-start <= maxdiff:
                     for j in xrange(start,i): good[j] = True
                     n_holes += 1
                 start = i+1
        print 'Closed', grace.pretty_number(n_holes), 'holes'
        
        
        f = open(path.join(output_dir, '%s-%s.fa' % (friendly_name,what)), 'wb')
        io.write_fasta(f, name,
            ''.join([ (seq[i] if good[i] else 'N')
                      for i in xrange(len(seq)) ])
        )
        f.close()

        f = open(path.join(output_dir, '%s-%s_masked.fa' % (friendly_name,what)), 'wb')
        io.write_fasta(f, name,
            ''.join([ (seq[i] if good[i] else seq[i].lower())
                      for i in xrange(len(seq)) ])
        )
        f.close()

        f_good = open(path.join(output_dir, '%s-%s_parts.fa' % (friendly_name,what)), 'wb')
        f_nongood = open(path.join(output_dir, '%s-non%s_parts.fa' % (friendly_name,what)), 'wb')
        start = 0
        n_good = [0]
        n_good_bases = [0]    
        def emit(i):
            if i-start < minsize: return
            if good[start]:
                n_good[0] += 1
                n_good_bases[0] += i-start
            io.write_fasta(
                f_good if good[start] else f_nongood,
                '%s:%d..%d' % (name, start+1,i),
                seq[start:i]
            )
        for i in xrange(1,len(seq)):
            if good[i] != good[start]:
                emit(i)
                start = i
        emit(len(seq))
        f_nongood.close()
        f_good.close()
        
        print grace.pretty_number(sum(good)), 'bases are '+what+', of', grace.pretty_number(len(seq)), 'in reference sequence'
        print grace.pretty_number(n_good[0]), 'parts at least', grace.pretty_number(minsize), 'bases long with', grace.pretty_number(n_good_bases[0]), 'total bases'

        print

