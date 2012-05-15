
import os

from nesoni import io, grace, config, annotation

class Reference(io.Workspace):
    def __init__(self, working_dir, must_exist):
        super(Reference,self).__init__(working_dir, must_exist)        
    
    def set_sequences(self, filenames):
        lengths = [ ]
        f = self.open('reference.fa', 'wb')
        for filename in filenames:
            for name, seq in io.read_sequences(filename):
                name = name.split()[0]
                lengths.append( (name, len(seq)) )
                io.write_fasta(f, name, seq)
        f.close()        
        self.set_object(lengths, 'reference-lengths.pickle.gz')
    
    def get_lengths(self):
        #Legacy working directory
        if not self.object_exists('reference-lengths.pickle.gz'):
            lengths = [ ]
            for name, seq in io.read_sequences(self.reference_fasta_filename()):
                name = name.split()[0]
                lengths.append( (name, len(seq)) )
            self.set_object(lengths, 'reference-lengths.pickle.gz')
                
        return self.get_object('reference-lengths.pickle.gz')
    
    def reference_fasta_filename(self):
        return self.object_filename('reference.fa')

    def build_shrimp_mmap(self, cs=False):
        suffix = '-cs' if cs else '-ls'
        
        grace.status('Building SHRiMP mmap')
        io.execute([
            'gmapper' + suffix,
            '--save', self.object_filename('reference' + suffix),
            self.reference_fasta_filename(),
        ])
        grace.status('')
    
    def shrimp_command(self, cs=False):
        suffix = '-cs' if cs else '-ls'
        result = [ 'gmapper' + suffix ]
        mmap_name = self.object_filename('reference' + suffix)
        if os.path.exists(mmap_name + '.genome'):
            result.extend([ '--load', mmap_name ])
        else:
            result.append( self.reference_fasta_filename() )
        return result
        
    def set_annotations(self, filenames):
        f = self.open('reference.gff','wb')
        print >> f, '##gff-version 3'
        for filename in filenames:
            for feature in annotation.read_annotations(filename):
                print >> f, feature.as_gff()
        f.close()
    
    def annotations_filename(self):
        name1 = self.object_filename('reference.gff')
        if os.path.exists(name1):
            return name1
        name2 = self.object_filename('reference.gbk')
        if os.path.exists(name2):
            return name2
        return None


@config.help("""\
Create a directory with a reference sequence, SHRiMP mmap files, \
and links to annotations.
""")
@config.Bool_flag('ls', 'Generate gmapper-ls mmap (faster SHRiMP startup for base-space reads).')
@config.Bool_flag('cs', 'Generate gmapper-cs mmap (faster SHRiMP startup for color-space reads).')
@config.Main_section('filenames', 'Sequence and annotation files.')
class Make_reference(config.Action_with_output_dir):
    ls = False
    cs = False
    filenames = [ ]

    def run(self):
        sequences = [ ]
        annotations = [ ]
        for filename in self.filenames:
            any = False
            if io.is_sequence_file(filename):
                sequences.append(filename)
                any = True
            if annotation.is_annotation_file(filename):
                annotations.append(filename)
                any = True            
            if not any:
                raise grace.Error(filename + ' is neither a sequence file nor an annotation file that nesoni can read.')
        
        reference = Reference(self.output_dir, must_exist=False)        
        reference.set_sequences(sequences)
        reference.set_annotations(annotations)
        if self.ls:
            reference.build_shrimp_mmap(False)
        if self.cs:
            reference.build_shrimp_mmap(True)
        
            
        
        
