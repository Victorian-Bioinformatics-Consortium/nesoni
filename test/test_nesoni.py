
import sys, os
sys.path.insert(0, os.path.join(os.path.split(__file__)[0], '../'))

import unittest

import nesoni
from nesoni import io

data = io.Workspace('data', must_exist=True)
output = io.Workspace('output', must_exist=False)


class Test_clip(unittest.TestCase):
    def test_single(self):
        nesoni.Clip(
            output / 'clip',
            reads = [ data/'reads_1.txt.gz' ],
        ).run()
    
    def test_paired(self):
        nesoni.Clip(
            output / 'clip',
            pairs = [ [ data/'reads_1.txt.gz', data/'reads_2.txt.gz' ] ],
        ).run()

    def test_custom_adapters(self):
        nesoni.Clip(
            output / 'clip',
            reads = [ data/'reads_1.txt.gz' ],
            adaptor_file = data/'test_adaptors.fa',
        ).run()

class Test_analyse_sample(unittest.TestCase):
    def setUp(self):
        nesoni.remake_needed()
        
        nesoni.Make_reference(
            output_dir=output/'reference', 
            filenames=[ data/'NC_001422.gbk' ],
            bowtie=True,
        ).run()

    def test_check_names(self):
        with self.assertRaises(AssertionError):
            nesoni.Analyse_sample(
                output_dir=output/'test-analyse', 
                reference=output/'reference', 
                clip=None,
                pairs=[[data/'reads_1.txt.gz',data/'reads_2.txt.gz'],
                       [data/'reads_1.txt.gz',data/'reads_2.txt.gz']]
            ).run()

        with self.assertRaises(AssertionError):
            nesoni.Analyse_sample(
                output_dir=output/'test-analyse', 
                reference=output/'reference', 
                clip=None,
                reads=[data/'reads_1.txt.gz',data/'reads_1.txt.gz']
            ).run()

        with self.assertRaises(AssertionError):
            nesoni.Analyse_sample(
                output_dir=output/'test-analyse', 
                reference=output/'reference', 
                clip=None,
                interleaved=[data/'reads_1.txt.gz']
            ).run()
    
    def test_analyse(self):        
        nesoni.Analyse_sample(
            output_dir=output/'test-analyse', 
            reference=output/'reference', 
            pairs=[[data/'reads_1.txt.gz',data/'reads_2.txt.gz']]
        ).run()

        nesoni.Nway(
            output=output/'test-nway.txt', 
            working_dirs=[output/'test-analyse']
        ).run()
        
        nesoni.Core(
            output_dir=output/'test-core',
            working_dirs=[output/'test-analyse'],
            what='core',
        ).run()

        nesoni.Core(
            output_dir=output/'test-unique',
            working_dirs=[output/'test-analyse'],
            what='unique',
        ).run()

    def test_bowtie(self):
        nesoni.Analyse_sample(
            output_dir=output/'test-analyse', 
            reference=output/'reference',
            pairs=[[data/'reads_1.txt.gz',data/'reads_2.txt.gz']],
            align=nesoni.Bowtie(),
        ).run()

if __name__ == '__main__':
    unittest.main()
