
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

outer = [ ]
if 'shrimp1' in args: outer.append(('shrimp1_',False)) 
if 'noshrimp2' not in args: outer.append(('',True))

for prefix3, is_sam in outer:
    for prefix1, section, sep_option in [('unpaired_','reads:',''),('paired_','pairs:','' if is_sam else ' --max-pair-sep 300')]:
        for prefix2, monogamous_option in [('monogamous_',''), ('polygamous_',' --monogamous 0')]:
            name = 'output/'+prefix1+prefix2+prefix3+'consensus'
            
            print
            print '='*70
            print name
            print
            
            for filename in glob.glob(os.path.join(name, '*')):
                if not os.path.isdir(filename):
                    os.unlink(filename)            
            
            run_nesoni(r""" \
                %sshrimp: %s \
                data/NC_001422_modified.gbk \
                %s \
                    data/reads_1.txt.gz \
                    data/reads_2.txt.gz
            """ % ('sam' if is_sam else '', name, section))
#                data/NC_001422_modified.fna \

            
            #run_nesoni(r""" \
            #    shrimp %s \
            #    --threshold 50%% \
            #    data/velvet_test_reference_modified.fa \
            #    reads \
            #        data/velvet_test_reads.fa \
            #    shrimp-options -n 1
            #""" % name)
            
            run_nesoni(('sam' if is_sam else '') + 'consensus: ' + name + ' ' + sep_option + monogamous_option)
            
         #   run_nesoni(('sam' if is_sam else '') +'reconsensus: ' + name)
            
            if not is_sam:
                run_nesoni('tosam ' + name)
            
            #if prefix2 == 'monogamous_':
            #    run_nesoni('consequences --use-coverage data/NC_001422-malformed.gbk ' + name + ' >' + name + '/consequences_output.txt',
            #               False)

