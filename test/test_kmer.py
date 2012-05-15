
import sys, os

SCRIPT = '../nesoni_scripts/nesoni'

def run_nesoni(command):
    full_command = '%s %s %s' % (sys.executable, SCRIPT, command)
    print
    print 'Running:'
    print full_command
    print
    assert 0 == os.system(full_command)


run_nesoni('bag kmer_test new 30 reads data/reads_1.txt.gz data/reads_2.txt.gz')
run_nesoni('graph kmer_test/graph new with kmer_test build 1')