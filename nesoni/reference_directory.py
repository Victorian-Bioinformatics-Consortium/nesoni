
import sys, os

from nesoni import io, grace, config, annotation

class Reference(io.Workspace):
    def __init__(self, working_dir, must_exist):
        super(Reference,self).__init__(working_dir, must_exist)
    
    def set_sequences(self, filenames):        
        reference_genbank_filename = self / 'reference.gbk'
        reference_filename = self / 'reference.fa'

        reference_genbank_file = open(reference_genbank_filename,'wb')
        any_genbank = [ False ]

        def genbank_callback(name, record):
            """ Make a copy of any genbank files passed in. """
            from Bio import SeqIO
            
            SeqIO.write([record], reference_genbank_file, 'genbank')
            
            f = open(self / (grace.filesystem_friendly_name(name) + '.gbk'), 'wb')
            SeqIO.write([record], f, 'genbank')
            f.close()
            
            any_genbank[0] = True
        
        lengths = [ ]
        f = open(reference_filename, 'wb')
        for filename in filenames:
            for name, seq in io.read_sequences(filename, genbank_callback=genbank_callback):
                name = name.split()[0]
                lengths.append( (name, len(seq)) )
                io.write_fasta(f, name, seq)
        f.close()        
        self.set_object(lengths, 'reference-lengths.pickle.gz')
        
        reference_genbank_file.close()
        if not any_genbank[0]:
            os.unlink(reference_genbank_filename)
            
        # Create an index of the reference sequences for samtools
        io.execute([
            'samtools', 'faidx', reference_filename
        ])

    
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

    def build_shrimp_mmap(self, cs=False, log_to=sys.stdout):
        suffix = '-cs' if cs else '-ls'
        
        grace.status('Building SHRiMP mmap')
        io.execute( [
                'gmapper' + suffix,
                '--save', self.object_filename('reference' + suffix),
                self.reference_fasta_filename(),            
                ],
            stdout=log_to
            )
        grace.status('')

    def build_bowtie_index(self, log_to=sys.stdout):
        io.execute([
                'bowtie2-build',
                self.reference_fasta_filename(),
                self/'bowtie',
                ],
            stdout = log_to,
            )
    
    def get_bowtie_index_prefix(self):
        assert os.path.exists(self/'bowtie.1.bt2'), 'bowtie2 index was not created. Please use "make-reference:" with "--bowtie yes".'
        return self/'bowtie'
    
    
    def shrimp_command(self, cs=False, parameters = [ ]):
        """ Parameters:
            First any parameters, then
            Single read file [ 'r.fa' ]
            Two paired reads file [ '-1', 'r1.fa', '-2', 'r2.fa' ]
        """
        suffix = '-cs' if cs else '-ls'
        result = [ 'gmapper' + suffix ] + parameters
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
@config.Bool_flag('bowtie', 'Generate bowtie2 index (necessary in order to use "nesoni bowtie:").')
@config.Main_section('filenames', 'Sequence and annotation files.')
class Make_reference(config.Action_with_output_dir):
    ls = False
    cs = False
    bowtie = False
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
        
        if not sequences:
            raise grace.Error('No reference sequence files given.')
        
        reference = Reference(self.output_dir, must_exist=False)        
        reference.set_sequences(sequences)
        reference.set_annotations(annotations)
        
        with open(self.log_filename(),'wb') as f:
            if self.ls:
                reference.build_shrimp_mmap(False, f)
            if self.cs:
                reference.build_shrimp_mmap(True, f)
            if self.bowtie:
                reference.build_bowtie_index(f)
            
        
        
