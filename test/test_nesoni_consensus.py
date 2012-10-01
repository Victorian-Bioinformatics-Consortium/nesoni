
import sys, os, glob

if not os.path.exists('output'):
    os.mkdir('output')

SCRIPT = '../nesoni_scripts/nesoni'

def run_nesoni(command, same_python=True):
    full_command = '%s %s' % (SCRIPT, command)
    if same_python:
        full_command = sys.executable + ' ' + full_command
    print
    print 'Running:'
    print full_command
    print
    assert 0 == os.system(full_command)

args = sys.argv[1:]

outer = [ 'shrimp', 'bowtie' ]

for aligner in outer:
    for prefix2, section in [('unpaired_','reads:'),('paired_','pairs:')]:
        for prefix3, monogamous_option in [('monogamous_',''), ('polygamous_',' --monogamous 0'),('random_',' --monogamous 0 --random 1')]:
            prefix1 = aligner + '_'
            name = 'output/'+prefix1+prefix2+prefix3+'consensus'
            
            print
            print '='*70
            print name
            print
            
            for filename in glob.glob(os.path.join(name, '*')):
                if not os.path.isdir(filename):
                    os.unlink(filename)            
            
            run_nesoni(r""" \
                %(aligner)s: %(name)s \
                data/NC_001422_modified.gbk \
                %(section)s \
                    data/reads_1.txt.gz \
                    data/reads_2.txt.gz
            """ % locals())
            
            run_nesoni('consensus: ' + name + ' ' + monogamous_option)
            
