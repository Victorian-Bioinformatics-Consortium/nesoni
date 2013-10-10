
import sys, os

import nesoni
from nesoni import io, grace, config, annotation, legion

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
        seen = set()
        f = open(reference_filename, 'wb')
        for filename in filenames:
            for name, seq in io.read_sequences(filename, genbank_callback=genbank_callback):
                name = name.split()[0]
                assert name not in seen, 'Duplicate chromosome name: ' + name
                seen.add(name)
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
        
    def set_annotations(self, filenames):
        f = self.open('reference.gff','wb')
        annotation.write_gff3_header(f)
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

    def build_shrimp_mmap(self, cs=False):
        suffix = '-cs' if cs else '-ls'
        
        grace.status('Building SHRiMP mmap')
        io.execute( [
                'gmapper' + suffix,
                '--save', self.object_filename('reference' + suffix),
                self.reference_fasta_filename(),            
                ],
            )
        grace.status('')

    def build_bowtie_index(self):
        io.execute([
                'bowtie2-build',
                self.reference_fasta_filename(),
                self/'bowtie',
                ],
            )
    
    def get_bowtie_index_prefix(self):
        assert os.path.exists(self/'bowtie.1.bt2'), 'bowtie2 index was not created. Please use "make-reference:" with "--bowtie yes".'
        return self/'bowtie'

    def build_genome(self, select):
        nesoni.Make_genome(
            prefix = self/self.name,
            name = self.name,
            select = select,
            filenames = [ self.working_dir ],
            ).run()
    
    def get_genome_filename(self):
        filename = self / (self.name+'.genome')
        assert os.path.exists(filename), 'IGV .genome was not created. Please use "make-reference:" with "--genome yes".'
        return filename
    
    def get_genome_dir(self):
        return self / self.name
    
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

    def build_snpeff(self):
        jar = io.find_jar('snpEff.jar')
        
        with open(self/'snpeff.config','wb') as f:
            print >> f, 'data_dir = snpeff'
            print >> f, 'genomes : ' + self.name
            print >> f, self.name + '.genome : ' + self.name 
        
        snpwork = io.Workspace(self/'snpeff',must_exist=False)
        snpwork_genome = io.Workspace(snpwork/self.name,must_exist=False)
        snpwork_genomes = io.Workspace(snpwork/'genomes',must_exist=False)
        
        annotations = self.annotations_filename()
        assert annotations
        with open(snpwork_genome/'genes.gff','wb') as f:
            for record in annotation.read_annotations(annotations):
                if record.end <= record.start: continue
                if not record.attr:
                    record.attr['attributes'] = 'none'
                print >> f, record.as_gff()
        
        with open(snpwork_genomes/(self.name+'.fa'),'wb') as f:
            for name, seq in io.read_sequences(self.reference_fasta_filename()):
                io.write_fasta(f, name, seq)
                
        io.execute('java -jar JAR build NAME -gff3 -c CONFIG',
            JAR=jar, NAME=self.name, CONFIG=self/'snpeff.config')


@config.help("""\
Create a directory with a reference sequence, \
and optionally reference annotations, \
SHRiMP mmap files, Bowtie2 index files, \
and IGV .genome file.

If sequence filenames are not given, \
the directory must already exist, \
and files will be generated as requested \
using existing sequences and annotations.
""")
@config.Ifavailable_flag('ls', 'Generate gmapper-ls mmap (faster SHRiMP startup for base-space reads).')
@config.Ifavailable_flag('cs', 'Generate gmapper-cs mmap (faster SHRiMP startup for color-space reads).')
@config.Ifavailable_flag('bowtie', 'Generate bowtie2 index (necessary in order to use "nesoni bowtie:").')
@config.Bool_flag('genome', 'Create .genome file and directory for use with IGV.')
@config.String_flag('genome_select', 'What types of feature to use in IGV .genome (selection expression).')
@config.Ifavailable_flag('snpeff', 'Create snpEff files.')
@config.Main_section('filenames', 'Sequence and annotation files.')
class Make_reference(config.Action_with_output_dir):
    ls = 'ifavailable'
    cs = 'ifavailable'
    bowtie = 'ifavailable'
    genome = True
    genome_select = '-source'
    snpeff = 'ifavailable'
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
            assert not annotations, 'Annotations given without any reference sequences.'
            reference = Reference(self.output_dir, must_exist=True)        
        else:
            reference = Reference(self.output_dir, must_exist=False)        
            reference.set_sequences(sequences)
            reference.set_annotations(annotations)
        
        with legion.Stage() as stage:
            if config.apply_ifavailable_program(self.ls, 'gmapper-ls'):
                stage.process(reference.build_shrimp_mmap, False)
            if config.apply_ifavailable_program(self.cs, 'gmapper-cs'):
                stage.process(reference.build_shrimp_mmap, True)
            if config.apply_ifavailable_program(self.bowtie, 'bowtie2-build'):
                stage.process(reference.build_bowtie_index)
            if self.genome:
                stage.process(reference.build_genome, self.genome_select)
            if config.apply_ifavailable_jar(self.snpeff, 'snpEff.jar'):
                stage.process(reference.build_snpeff)
            
        
        
