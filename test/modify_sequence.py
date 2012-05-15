


# python2.6 modify_sequence.py data/velvet_test_reference.fa >data/velvet_test_reference_modified.fa

import sys, random

from nesoni import io

for name, seq in io.read_sequences(sys.argv[1]):
    j = 0    
    for i in xrange(0,len(seq)-100,100):
        original = seq[i]
        if j%3 == 0:
            while True:
                new = random.choice('ACGT')
                if new != original: break
            seq = seq[:i] + new + seq[i+1:]
        elif j%3 == 1:
            n = (j // 3) % 9 + 1
            seq = seq[:i] + seq[i+n:]
        else:
            n = (j // 3) % 9 + 1
            seq = seq[:i] + ''.join( random.choice('ACGT') for k in xrange(n) ) + seq[i:]
        j += 1
    io.write_fasta(sys.stdout, name, seq)
