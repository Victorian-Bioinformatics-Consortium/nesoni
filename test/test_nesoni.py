
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

if __name__ == '__main__':
    unittest.main()