
import subprocess

from nesoni import config, workspace, working_directory, reference_directory, io, grace


@config.help("""\

Align reads using Bowtie 2.

Paired end reads should either be given as two files in the "pairs"
section, or as a single interleaved file in the "interleaved" section.

""")
@config.Main_section('references', 
    'Reference sequence filenames, '
    'or a directory created using "nesoni make-reference: --bowtie yes" (recommended).')
@config.Section('reads', 'FASTQ files containing unpaired reads.')
@config.Section('interleaved', 'FASTQ files containing interleaved read pairs.')
@config.Grouped_section('pairs', 'Pair of FASTQ files containing read pairs.')
@config.Section('bowtie_options', 'Options to pass to bowtie.', allow_flags=True)
class Bowtie(config.Action_with_output_dir):
    cs = False
    sam_unaligned = True
    half_paired = True
    references = []
    reads = []
    interleaved = []
    pairs = []
    bowtie_options = []
    
    _workspace_class = working_directory.Working

    def run(self):
        assert self.reads or self.pairs or self.interleaved, 'No reads given'
    
        io.check_name_uniqueness(self.reads, self.pairs, self.interleaved)
        
        working = self.get_workspace()
        working.setup_reference(self.references, bowtie=True)        
        reference = working.get_reference()
        
        log_file = open(self.log_filename(),'wb')
              
        with workspace.tempspace(dir=working.working_dir) as temp:
            n = [ 0 ]
            def tempname():
                n[0] += 1
                return temp/('%d.fq'%n[0])
            def convert(filename):
                ok = io.get_compression_type(filename) in ('none','gzip','bzip2')
                if ok:
                    ok = io.get_file_type(filename) == 'fastq'
                if ok:
                    return filename            
                result_name = tempname()
                with open(result_name,'wb') as f:
                    for name, seq, qual in io.read_sequences(filename, qualities='required'):
                        io.write_fastq(f, name, seq, qual)
                return result_name
            
            ones = [ ]
            twos = [ ]
            singles = [ ]
            
            for pair in self.pairs:
                assert len(pair) == 2, 'Need two files in each "pair:" section.'
                ones.append(convert(pair[0]))
                twos.append(convert(pair[1]))
            
            for item in self.interleaved:
                left_name = tempname()
                right_name = tempname()
                ones.append(left_name)
                twos.append(right_name)
                with open(left_name,'wb') as left, \
                     open(right_name,'wb') as right:
                    reader = io.read_sequences(item, qualities='required')
                    while True:
                        try:
                            name, seq, qual = reader.next()
                        except StopIteration:
                            break
                        io.write_fastq(left, name,seq,qual)
                        
                        try:
                            name, seq, qual = reader.next()
                        except StopIteration:
                            raise grace.Error('Interleaved file contains odd number of sequences')
                        io.write_fastq(right, name,seq,qual)
            
            for item in self.reads:
                singles.append(convert(item))

            command = [ 
                'bowtie2', '-x', reference.get_bowtie_index_prefix(),
                '--rg-id', '1',
                '--rg', 'SM:'+working.name
            ]
            if ones:
                command.extend([ '-1', ','.join(ones), '-2', ','.join(twos) ])
            if singles:
                command.extend([ '-U', ','.join(singles) ])
            
            command.extend(self.bowtie_options)
            
            self.log.log('Running:\n' + ' '.join(command) + '\n')
            
            
            temp_bam_name = temp/'temp.bam'

            writer = io.run(['samtools', 'view', '-S', '-b', '-'],
                stdout=open(temp_bam_name,'wb'),
                stdin=subprocess.PIPE,
                stderr=log_file
            )            

            bowtie = io.run(command,
                stdout=writer.stdin,
                stderr=log_file
            )

            writer.stdin.close()
            
            assert bowtie.wait() == 0, 'bowtie2 failed, consult log: '+self.log_filename()
            assert writer.wait() == 0, 'samtools failed, consult log: '+self.log_filename()

            io.execute([
                'samtools', 'sort', '-n', temp_bam_name, working/'alignments'
            ])
            
        log_file.close()




