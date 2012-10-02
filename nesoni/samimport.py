
from nesoni import io, bio, grace, sam, config, working_directory

import os, sys

@config.help("""\

Set up a nesoni working directory using a SAM or BAM file.

""")
@config.Float_flag('snp_cost', 'Cost per SNP in alignment scores of BAM.')
@config.Positional('input', 'SAM or BAM file.')
@config.Main_section('reference', 'Reference directory from "make-reference:", or reference sequences and annotations.')
class Import(config.Action_with_output_dir):
    snp_cost = 2.0
    input = None
    reference = [ ]
    
    def run(self):
        workspace = working_directory.Working(self.output_dir)        
        workspace.setup_reference(self.reference)
        workspace.update_param(snp_cost = self.snp_cost)
        
        #assert os.path.exists(self.reference), 'Reference file does not exist'
        #reference_filename = workspace._object_filename('reference.fa')
        #if os.path.exists(reference_filename):
        #   os.unlink(reference_filename)
        #os.symlink(os.path.relpath(self.reference, self.output_dir), reference_filename)
        
        bam_filename = io.abspath(self.output_dir, 'alignments.bam')
        bam_prefix = io.abspath(self.output_dir, 'alignments')

        if sam.is_bam(self.input):
            sort_input_filename = self.input
            temp_filename = None
        else:
            temp_filename = io.abspath(self.output_dir, 'temp.bam')
            sort_input_filename = temp_filename
            writer = io.Pipe_writer(temp_filename, ['samtools', 'view', '-S', '-b', '-'])
            f = open(self.input, 'rb')
            while True:
                data = f.read(1<<20)
                if not data: break
                writer.write(data)
            writer.close()
            f.close()
        
        grace.status('Sort')
        
        io.execute([
            'samtools', 'sort', '-n', sort_input_filename, bam_prefix
        ])
        
        if temp_filename is not None:
            os.unlink(temp_filename)
        
        grace.status('')




