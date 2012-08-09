

import sys, os

from nesoni import io, grace, config

@config.help("""\
Break a sequence or sequences into small overlapping pieces.
""")
@config.Int_flag('size', 'Size of sequences to output.')
@config.Int_flag('stride', ' Step size along input sequence.')
@config.Main_section('filenames','Files containing sequences.', empty_is_ok=False)
class Shred(config.Action_with_optional_output):
    size = 200
    stride = 50
    filenames = [ ]
    
    def run(self):    
        f = self.begin_output()
    
        for filename in self.filenames:
            for name, seq in io.read_sequences(filename):
                name_parts = name.split(None, 1)
                name = name_parts[0]
                for i in xrange(-self.size+self.stride,len(seq),self.stride):
                    start = max(0,min(len(seq),i))
                    end = max(0,min(len(seq), i+self.size))
                    io.write_fasta(
                        f,
                        '%s:%d..%d' % (name,start+1,end),
                        seq[start:end]
                    )

        self.end_output(f)








