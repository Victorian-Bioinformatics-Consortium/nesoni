
"""

Clip based on quality and presence of adaptor sequence

"""

import sys, os, itertools, collections

from nesoni import grace, io, bio, config

#python -c 'from nesoni import io;import pprint;pprint.pprint(list(io.read_sequences("...")),width=1000)'
ADAPTORS = [
 ('NexteraTM_DNA_Sample_Preparation_Kits_v2_transposomes_1', 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'),
 ('NexteraTM_DNA_Sample_Preparation_Kits_v2_transposomes_2', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'),
 ('TruSeq_Universal_Adapter', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('TruSeq_Adapter_Index_1', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_2', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_3', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_4', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_5', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_6', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_7', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_8', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_9', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_10', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_11', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_12', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_13', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_14', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_15', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_16', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_18_4', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_19', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_20', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_21', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_22', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_23', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_25', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG'),
 ('TruSeq_Adapter_Index_27', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG'),
 ('RNA_Adapter_RA5_part_15013205', 'GTTCAGAGTTCTACAGTCCGACGATC'),
 ('RNA_3p_Adapter_RA3_part_15013207', 'TGGAATTCTCGGGTGCCAAGG'),
 ('RNA_RT_Primer_RTP_part_15013981', 'GCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_RP1_part_15005505', 'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA'),
 ('RNA_PCR_Primer_Index_1_RPI1_5', 'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_2_RPI2', 'CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_3_RPI3', 'CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_4_RPI4', 'CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_5_RPI5', 'CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_6_RPI6', 'CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_7_RPI7', 'CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_8_RPI8', 'CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_9_RPI9', 'CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_10_RPI10', 'CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_11_RPI11', 'CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_12_RPI12', 'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_13_RPI13', 'CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_14_RPI14', 'CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_15_RPI15', 'CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_16_RPI16', 'CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_17_RPI17', 'CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_18_RPI18', 'CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_19_RPI19', 'CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_20_RPI20', 'CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_21_RPI21', 'CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_22_RPI22', 'CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_23_RPI23', 'CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_24_RPI24', 'CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_25_RPI25', 'CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_26_RPI26', 'CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_27_RPI27', 'CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_28_RPI28', 'CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_29_RPI29', 'CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_30_RPI30', 'CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_31_RPI31', 'CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_32_RPI32', 'CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_33_RPI33', 'CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_34_RPI34', 'CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_35_RPI35', 'CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_36_RPI36', 'CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_37_RPI37', 'CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_38_RPI38', 'CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_39_RPI39', 'CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_40_RPI40', 'CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_41_RPI41', 'CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_42_RPI42', 'CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_43_RPI43', 'CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_44_RPI44', 'CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_45_RPI45', 'CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_46_RPI46', 'CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_47_RPI47', 'CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('RNA_PCR_Primer_Index_48_RPI48', 'CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'),
 ('Oligonucleotide_sequences_for_Genomic_DNA_Adapters_1', 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'),
 ('Oligonucleotide_sequences_for_Genomic_DNA_Adapters_2', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('PCR_Primers_1', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('PCR_Primers_2', 'CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT'),
 ('Genomic_DNA_Sequencing_Primer', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('PE_Adapters_1', 'GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'),
 ('PE_Adapters_2', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('PE_PCR_Primer_1.0', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('PE_PCR_Primer_2.0', 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'),
 ('PE_Read_1_Sequencing_Primer', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('PE_Read_2_Sequencing_Primer', 'CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'),
 ('Multiplexing_Adapters_1', 'GATCGGAAGAGCACACGTCT'),
 ('Multiplexing_Adapters_2', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('Multiplexing_PCR_Primer_1.0', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('Multiplexing_PCR_Primer_2.0', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'),
 ('Multiplexing_Read_1_Sequencing_Primer', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'),
 ('Multiplexing_Index_Read_Sequencing_Primer', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'),
 ('Multiplexing_Read_2_Sequencing_Primer', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'),
 ('PCR_Primer_Index_1', 'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC'),
 ('PCR_Primer_Index_2', 'CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC'),
 ('PCR_Primer_Index_3', 'CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC'),
 ('PCR_Primer_Index_4', 'CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC'),
 ('PCR_Primer_Index_5', 'CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC'),
 ('PCR_Primer_Index_6', 'CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC'),
 ('PCR_Primer_Index_7', 'CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC'),
 ('PCR_Primer_Index_8', 'CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC'),
 ('PCR_Primer_Index_9', 'CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC'),
 ('PCR_Primer_Index_10', 'CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC'),
 ('PCR_Primer_Index_11', 'CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC'),
 ('PCR_Primer_Index_12', 'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC'),
 ('RT_Primer', 'CAAGCAGAAGACGGCATACGA'),
 ('5p_RNA_Adapter', 'GTTCAGAGTTCTACAGTCCGACGATC'),
 ('3p_RNA_Adapter', 'TCGUATGCCGTCTTCTGCTTGT'),
 ('v1.5_Small_RNA_3p_Adapter', 'ATCTCGTATGCCGTCTTCTGCTTG'),
 ('Small_RNA_PCR_Primer_1', 'CAAGCAGAAGACGGCATACGA'),
 ('Small_RNA_PCR_Primer_2', 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'),
 ('Small_RNA_Sequencing_Primer', 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC')
]



class Matcher(object):
    def __init__(self, seqs, names, max_error):
        self.max_error = max_error        
        self.names = names

        self.transitions = {'A':[], 'C':[], 'G':[], 'T':[] } # base -> id -> (new state id, emission)

        state_ids = { }

        def advance( state, base ):
            """ returns ( [new_state], [emit] ) """
        
            seq_no, pos, errors = state
            
            if seqs[seq_no][pos] != 'N' and seqs[seq_no][pos] != base:
                errors += 1
            
            if errors > max_error:
                return ( [], [] )
            
            pos += 1
            if pos == len(seqs[seq_no]):
                return ( [], [ (errors, seq_no) ] )
            else:
                return ( [ (seq_no, pos, errors) ], [] )


        todo = [ ]        
        def get_id(state):
            state = tuple(sorted(state))
            if state not in state_ids:
                state_ids[state] = len(state_ids)
                for base in 'ACGT':
                    self.transitions[base].append( None )
                todo.append(state)
            return state_ids[state]

        initial = [ ]
        for i in xrange(len(seqs)):
            for j in xrange(len(seqs[i])):
                initial.append( (i,j,0) )
        initial = tuple(sorted(initial))
        
        self.initial_id = get_id(initial)

        while todo:
            from_state = todo.pop(-1)
            from_id = get_id(from_state)
            
            for base in 'ACGT':
                 new_state = [ ]
                 emit = [ ]
                 for sub_state in from_state:
                     sub_new, sub_emit = advance(sub_state, base)
                     new_state.extend(sub_new)
                     emit.extend(sub_emit)
                 new_id = get_id(new_state)
                 self.transitions[ base ][ from_id ] = (new_id, sorted(emit))
    
    def match(self, seq):
        state = self.initial_id
        
        emit = None
        for i, base in enumerate(seq):        
            state, this_emit = self.transitions[ base ][ state ]
            if this_emit:
                emit = (i+1, this_emit)
        return emit


def deinterleave(iterator):
    for item1 in iterator:
        item2 = iterator.next()
        yield (item1, item2)


@config.help("""\
Clip adaptors and low quality bases from Illumina reads or read pairs.
""")
@config.Bool_flag('adaptor_clip',
    'Do adaptor clipping?'
)
@config.String_flag('adaptor_file',
    'FASTA file to read adaptors from. '
    'Defaults to a built-in list of Illumina adaptors.'
)
@config.Int_flag('match',
    'Minimum length adaptor match.'
)
@config.Int_flag('max_errors',
    'Maximum errors in adaptor match.\n'
    'Note: slow for N > 1'
)
@config.Bool_flag('clip_ambiguous', 
    'Clip ambiguous bases, eg N'
)
@config.Int_flag('quality',  
    'Quality cutoff'
)
@config.Int_flag('qoffset',
    'Quality character offset.\n'
    'sanger: 33, solexa: 59, illumina: 64\n'
    'Will guess if not given.'
)
@config.Int_flag('length',
    'Reads shorter than this will be discarded.'
)
@config.Bool_flag('homopolymers',  
    'Set to yes to *remove* reads containing all the same base.'
) 
@config.Int_flag('trim_start',
    'Trim NNN bases from start of each read irrespective of quality.'
)
@config.Int_flag('trim_end',                      
    'Trim NNN bases from end of each read irrespective of quality.'
)
@config.Bool_flag('revcom',
    'Reverse complement all reads.\n'
    'ie convert opp-out pairs to opp-in'
)
@config.Bool_flag('fasta',
    'Output in fasta rather than fastq format.'
)
@config.Bool_flag('gzip',
    'Gzip output.'
)
@config.Bool_flag('rejects',  
    'Output rejected reads in separate file.'
)
@config.Section('reads', 'Files containing unpaired reads.')
@config.Section('interleaved', 'Files containing interleaved read pairs.')
@config.Grouped_section('pairs', 'Pair of files containing read pairs.') 
class Clip(config.Action_with_prefix):
    quality = 10
    qoffset = None
    clip_ambiguous = True
    length = 24
    match = 10
    max_errors = 1
    adaptor_clip = True
    adaptor_file = None
    homopolymers = False
    revcom = False
    trim_start = 0
    trim_end = 0
    fasta = False
    gzip = True
    rejects = False
    
    reads = [ ]
    pairs = [ ]
    interleaved = [ ]


    def output_suffix(self):
        if self.fasta:
            suffix = '.fa'
        else:
            suffix = '.fq'
        if self.gzip:
            suffix += '.gz'
        return suffix

    def reads_output_filenames(self):
        return [ self.prefix + '_single' + self.output_suffix() ]
        
    def interleaved_output_filenames(self):
        if not self.pairs and not self.interleaved:
            return [ ]
        return [ self.prefix + '_paired' + self.output_suffix() ]

    def rejects_output_filenames(self):
        if not self.rejects:
            return [ ]
        return [ self.prefix + '_rejected' + self.output_suffix() ]


    def run(self):
        log = self.log
        
        #quality_cutoff, args = grace.get_option_value(args, '--quality', int, 10)
        #qoffset, args = grace.get_option_value(args, '--qoffset', int, None)
        #clip_ambiguous, args = grace.get_option_value(args, '--clip-ambiguous', grace.as_bool, True)
        #length_cutoff, args = grace.get_option_value(args, '--length', int, 24)
        #adaptor_cutoff, args = grace.get_option_value(args, '--match', int, 10)
        #max_error, args = grace.get_option_value(args, '--max-errors', int, 1)
        #adaptor_set, args = grace.get_option_value(args, '--adaptors', str, 'truseq-adapter,truseq-srna,genomic,multiplexing,pe,srna')
        #disallow_homopolymers, args = grace.get_option_value(args, '--homopolymers', grace.as_bool, False)
        #reverse_complement, args = grace.get_option_value(args, '--revcom', grace.as_bool, False)
        #trim_start, args = grace.get_option_value(args, '--trim-start', int, 0)
        #trim_end, args = grace.get_option_value(args, '--trim-end', int, 0)
        #output_fasta, args = grace.get_option_value(args, '--fasta', grace.as_bool, False)
        #use_gzip, args = grace.get_option_value(args, '--gzip', grace.as_bool, True)
        #output_rejects, args = grace.get_option_value(args, '--rejects', grace.as_bool, False)
        #grace.expect_no_further_options(args)
        
        prefix = self.prefix
        log_name = os.path.split(prefix)[1]
        
        quality_cutoff = self.quality
        qoffset = self.qoffset
        clip_ambiguous = self.clip_ambiguous
        length_cutoff = self.length
        adaptor_cutoff = self.match
        max_error = self.max_errors
        disallow_homopolymers = self.homopolymers
        reverse_complement = self.revcom
        trim_start = self.trim_start
        trim_end = self.trim_end
        output_fasta = self.fasta
        use_gzip = self.gzip
        output_rejects = self.rejects
    
        iterators = [ ]        
        filenames = [ ]
        any_paired = False
        
        for filename in self.reads:
            filenames.append(filename)
            iterators.append(itertools.izip(
                 io.read_sequences(filename, qualities=True)
            ))
        
        for pair_filenames in self.pairs:
            assert len(pair_filenames) == 2, 'Expected a pair of files for "pairs" section.'
            filenames.extend(pair_filenames)
            any_paired = True
            iterators.append(itertools.izip(
                io.read_sequences(pair_filenames[0], qualities=True),
                io.read_sequences(pair_filenames[1], qualities=True)
            ))
        
        for filename in self.interleaved:
            filenames.append(filename)
            any_paired = True
            iterators.append(deinterleave(
                io.read_sequences(filename, qualities=True)
            ))
        
        fragment_reads = (2 if any_paired else 1)
        read_in_fragment_names = [ 'read-1', 'read-2' ] if any_paired else [ 'read' ]
        
        assert iterators, 'Nothing to clip'
        
        io.check_name_uniqueness(self.reads, self.pairs, self.interleaved)
    
        if qoffset is None:
            guesses = [ io.guess_quality_offset(filename) for filename in filenames ]
            assert len(set(guesses)) == 1, 'Conflicting quality offset guesses, please specify manually.'
            qoffset = guesses[0]
            log.log('FASTQ offset seems to be %d\n' % qoffset)    
    
        quality_cutoff_char = chr(qoffset + quality_cutoff)
        
        #log.log('Minimum quality:        %d (%s)\n' % (quality_cutoff, quality_cutoff_char))
        #log.log('Clip ambiguous bases:   %s\n' % (grace.describe_bool(clip_ambiguous)))
        #log.log('Minimum adaptor match:  %d bases, %d errors\n' % (adaptor_cutoff, max_error))
        #log.log('Minimum length:         %d bases\n' % length_cutoff)
        
        adaptor_seqs = [ ]
        adaptor_names = [ ]
        if self.adaptor_clip:
            if self.adaptor_file:
                adaptor_iter = io.read_sequences(self.adaptor_file)
            else:
                adaptor_iter = ADAPTORS
            for name, seq in adaptor_iter:
                seq = seq.upper().replace('U','T')
                adaptor_seqs.append(seq)
                adaptor_names.append(name)
                adaptor_seqs.append(bio.reverse_complement(seq))
                adaptor_names.append(name)

        matcher = Matcher(adaptor_seqs, adaptor_names, max_error)
        
        start_clips = [ collections.defaultdict(list) for i in xrange(fragment_reads) ]
        end_clips = [ collections.defaultdict(list) for i in xrange(fragment_reads) ]
    
        if output_fasta:
            write_sequence = io.write_fasta_single_line
        else:
            write_sequence = io.write_fastq
    
        f_single = io.open_possibly_compressed_writer(self.reads_output_filenames()[0])
        if fragment_reads == 2:
            f_paired = io.open_possibly_compressed_writer(self.interleaved_output_filenames()[0])
        if output_rejects:
            f_reject = io.open_possibly_compressed_writer(self.rejects_output_filenames()[0])
        
        n_single = 0
        n_paired = 0
        
        n_in_single = 0
        n_in_paired = 0
        total_in_length = [ 0 ] * fragment_reads
        
        n_out = [ 0 ] * fragment_reads
        n_q_clipped = [ 0 ] * fragment_reads
        n_a_clipped = [ 0 ] * fragment_reads
        n_homopolymers = [ 0 ] * fragment_reads
        total_out_length = [ 0 ] * fragment_reads
        
        #log.attach(open(prefix + '_log.txt', 'wb'))
        
        for iterator in iterators:
          for fragment in iterator:
            if (n_in_single+n_in_paired) % 10000 == 0:
                grace.status('Clipping fragment %s' % grace.pretty_number(n_in_single+n_in_paired))
        
            if len(fragment) == 1:
                n_in_single += 1
            else:
                n_in_paired += 1
            
            graduates = [ ]
            rejects = [ ]
            for i, (name, seq, qual) in enumerate(fragment):
                seq = seq.upper()
                total_in_length[i] += len(seq)
                                    
                start = trim_start
                best_start = 0
                best_len = 0
                for j in xrange(len(seq)-trim_end):
                    if qual[j] < quality_cutoff_char or \
                       (clip_ambiguous and seq[j] not in 'ACGT'):
                        if best_len < j-start:
                            best_start = start
                            best_len = j-start
                        start = j + 1
                j = len(seq)-trim_end
                if best_len < j-start:
                    best_start = start
                    best_len = j-start
        
                clipped_seq = seq[best_start:best_start+best_len]
                clipped_qual = qual[best_start:best_start+best_len]
                if len(clipped_seq) < length_cutoff:
                    n_q_clipped[i] += 1
                    rejects.append( (name,seq,qual,'quality') ) 
                    continue
        
                match = matcher.match(clipped_seq)
                if match and match[0] >= adaptor_cutoff:
                    clipped_seq = clipped_seq[match[0]:]
                    clipped_qual = clipped_qual[match[0]:]
                    start_clips[i][match[0]].append( match[1][0] )
                    if len(clipped_seq) < length_cutoff:
                        n_a_clipped[i] += 1 
                        rejects.append( (name,seq,qual,'adaptor') ) 
                        continue
            
                match = matcher.match(bio.reverse_complement(clipped_seq))
                if match and match[0] >= adaptor_cutoff:
                    clipped_seq = clipped_seq[: len(clipped_seq)-match[0] ]    
                    clipped_qual = clipped_qual[: len(clipped_qual)-match[0] ]    
                    end_clips[i][match[0]].append( match[1][0] )
                    if len(clipped_seq) < length_cutoff:
                        n_a_clipped[i] += 1 
                        rejects.append( (name,seq,qual,'adaptor') ) 
                        continue
    
                if disallow_homopolymers and len(set(clipped_seq)) <= 1:
                    n_homopolymers[i] += 1
                    rejects.append( (name,seq,qual,'homopolymer') ) 
                    continue
        
                graduates.append( (name, clipped_seq, clipped_qual) )
                n_out[i] += 1
                total_out_length[i] += len(clipped_seq)
    
            if output_rejects:
                for name,seq,qual,reason in rejects:
                    write_sequence(f_reject, name + ' ' + reason, seq, qual)
             
            if graduates:
                if reverse_complement:
                    graduates = [
                        (name, bio.reverse_complement(seq), qual[::-1])
                        for name, seq, qual in graduates
                    ]
            
                if len(graduates) == 1:
                    this_f = f_single
                    n_single += 1
                else:
                    assert len(graduates) == 2
                    this_f = f_paired
                    n_paired += 1
                
                for name, seq, qual in graduates:
                    write_sequence(this_f, name, seq, qual)
        
        grace.status('')
        
        if output_rejects:
            f_reject.close()
        if fragment_reads == 2:
            f_paired.close()
        f_single.close()
        
        def summarize_clips(name, location, clips):
            total = 0
            for i in clips:
                total += len(clips[i])
            log.datum(log_name, name + ' adaptors clipped at ' + location, total) 
            
            if not clips:
                return
    
            for i in xrange(min(clips), max(clips)+1):
                item = clips[i]
                log.quietly_log('%3d bases: %10d ' % (i, len(item)))
                if item:
                    avg_errors = float(sum( item2[0] for item2 in item )) / len(item)
                    log.quietly_log(' avg errors: %5.2f  ' % avg_errors)
                    
                    counts = collections.defaultdict(int)
                    for item2 in item: counts[item2[1]] += 1
                    #print counts
                    for no in sorted(counts,key=lambda item2:counts[item2],reverse=True)[:2]:
                        log.quietly_log('%dx%s ' % (counts[no], matcher.names[no]))
                    if len(counts) > 2: log.quietly_log('...')
                    
                log.quietly_log('\n')
            log.quietly_log('\n')


        if n_in_paired:
            log.datum(log_name,'read-pairs', n_in_paired)                      
        if n_in_single:
            log.datum(log_name,'single reads', n_in_single)                      
        
        for i in xrange(fragment_reads):
            if start_clips:
                summarize_clips(read_in_fragment_names[i], 'start', start_clips[i])
        
            if end_clips:
                summarize_clips(read_in_fragment_names[i], 'end', end_clips[i])

                prefix = read_in_fragment_names[i]
                
            log.datum(log_name, prefix + ' too short after quality clip', n_q_clipped[i])
            log.datum(log_name, prefix + ' too short after adaptor clip', n_a_clipped[i])
            if disallow_homopolymers:
                log.datum(log_name, prefix + ' homopolymers', n_homopolymers[i])
            if fragment_reads > 1:
                log.datum(log_name, prefix + ' kept', n_out[i])
            log.datum(log_name, prefix + ' average input length',  float(total_in_length[i]) / (n_in_single+n_in_paired))                     
            if n_out[i]:
                log.datum(log_name, prefix + ' average output length', float(total_out_length[i]) / n_out[i])                     
        
        if fragment_reads == 2:
            log.datum(log_name,'pairs kept after clipping', n_paired)                      
        log.datum(log_name, 'reads kept after clipping', n_single)
        
        


