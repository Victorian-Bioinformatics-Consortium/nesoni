

import sys, os

from nesoni import io, grace


USAGE = """

Usage:

    nesoni shred: [options] sequence.fa >output.fa

Options:

    --size      Size of sequence to output.
                Default 200.
        
    --stride    Step size along input sequence.
                Default 50.

"""


def main(args):
    size, args = grace.get_option_value(args,'--size',int,200)
    stride, args = grace.get_option_value(args,'--stride',int,50)
    grace.expect_no_further_options(args)
    
    if not args:
        print USAGE
        return 1
    
    for filename in args:
        for name, seq in io.read_sequences(filename):
            name_parts = name.split(None, 1)
            name = name_parts[0]
            if len(name_parts) > 1:
               desc = ' ' + name_parts[1]
            else:
               desc = ''
            
            for i in xrange(-size+stride,len(seq),stride):
                start = max(0,min(len(seq),i))
                end = max(0,min(len(seq), i+size))
                io.write_fasta(
                    sys.stdout,
                    '%s:%d..%d' % (name,start+1,end) + desc,
                    seq[start:end]
                )
    
    return 0