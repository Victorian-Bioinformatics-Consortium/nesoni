
import nesoni

from nesoni import config

class Analyse_sample(config.Action_with_output_dir):
    reads = [ ]
    interleaved = [ ]
    pairs = [ ]

    clip = None
    align = None
    filter = None
    consensus = None
    
    def run(self):
        reads = self.reads
        interleaved = self.interleaved
        pairs = self.pairs
                
        if self.clip:
            workspace = io.Workspace(self.output, must_exist=False)
            act = self.clip(
                workspace/'clipped',
                reads = reads,
                interleaved = interleaved,
                pairs = pairs,
            )
            act.make()
            reads = act.reads_out_filenames()
            interleaved = act.interleaved_out_filenames()
            pairs = [ ]
        
        self.align(
            self.output_dir,
            reads = reads,
            interleaved = interleaved,
            pairs = pairs, 
        ).make()
        
        if self.filter:
            self.filter(
                self.output_dir
            ).make()
        
        if self.consensus:
            self.consensus(
                self.output_dir
            ).make()

