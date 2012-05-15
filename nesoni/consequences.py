
"""

Determine consequences of SNPs for amino acid coding

"""

from Bio import Seq, SeqIO
from Bio.Data import CodonTable
import gzip, string, bisect, sys, os, numpy
from numpy import linalg

from nesoni import io, bio, grace

def warn(message):
    sys.stderr.write('* Warning: %s\n' % message)


def band_limited_align(str1, str2, bandwidth):
    scores = { }   # (pos1, pos2, score)
    
    length = len(str1)
    assert length == len(str2) #Pad if necessary
    
    scores[(length,length)] = 0
    for pos in xrange(max(0, length-bandwidth),length):
        scores[(length,pos)] = 0
        scores[(pos,length)] = 0
    
    for pos1 in xrange(length-1,-1,-1):
        for pos2 in xrange(min(pos1+bandwidth, length-1),
                           max(pos1-bandwidth, 0)-1,
                           -1):
            score = scores[(pos1+1,pos2+1)]
            if str1[pos1] != str2[pos2]:
                score += 1
            
            below = (pos1,pos2+1)
            if below in scores:
                score = min( scores[below]+2, score )
            
            right = (pos1+1,pos2)
            if right in scores:
                score = min( scores[right]+2, score )
            
            scores[(pos1,pos2)] = score

    al1 = [ ]
    al2 = [ ]
    pos1 = 0
    pos2 = 0
    while pos1 < length and pos2 < length:
        diagonal = (pos1+1,pos2+1)
        below = (pos1+1,pos2)
        right = (pos1,pos2+1)
        
        if (below not in scores or scores[below] >= scores[diagonal]) and \
           (right not in scores or scores[right] >= scores[diagonal]):
           al1.append(str1[pos1])
           al2.append(str2[pos2])
           pos1 += 1
           pos2 += 1
        
        elif right not in scores or (below in scores and scores[right] >= scores[below]):
           al1.append( str1[pos1] )
           al2.append( '-' )
           pos1 += 1

        else:
           al1.append( '-' )
           al2.append( str2[pos2] )
           pos2 += 1
    
    return ''.join(al1), ''.join(al2)



COMPLEMENTER = string.maketrans('ABCDGHKMRSTVWY','TVGHCDMKYSABWR')
def reverse_complement(seq):
    return seq.translate(COMPLEMENTER)[::-1]

def sequence_slice(seq, start,end):
    if start > end:
        return reverse_complement(sequence_slice(seq,end,start))
    
    start = min(len(seq),max(0,start))
    end = min(len(seq),max(0,end))
    return seq[start:end]


class Half_alignment:
    def __init__(self, text):
        self.data = [0]
        
        for char in text:
            if char == '-':
                self.data.append(self.data[-1])
            else:
                self.data.append(self.data[-1] + 1) 
                
    def __getitem__(self,key):
        if key < 0:
            return key
        if key >= len(self.data):
            return self.data[-1]+key-(len(self.data)-1)
        return self.data[key]

    def find(self, value, left=True):
        if value <= 0:
            return value
        last_value = self.data[-1]
        if value >= last_value:
            return value-last_value+(len(self.data)-1)
        
        if left:
            return bisect.bisect_left(self.data, value)
        else:
            return bisect.bisect_right(self.data, value)-1

class Half_alignment_simple:
    def __getitem__(self, key): return key
    def find(self, value, left=True): return value
half_alignment_simple = Half_alignment_simple()


class Alignment:
    def __init__(self, start1, end1, forward1, text1, start2, end2, forward2, text2):
        self.start1 = start1
        self.end1 = end1
        self.forward1 = forward1
        if text1 is None:
            self.ali1 = half_alignment_simple
        else:
            self.ali1 = Half_alignment(text1)
        
        self.start2 = start2
        self.end2 = end2
        self.forward2 = forward2
        if text2 is None:
            self.ali2 = half_alignment_simple
        else:
            self.ali2 = Half_alignment(text2)

    def project(self, position, left=True):
        """ Positions lie between bases.
            start and end are positions. """
            
        if self.forward1:
            position = position - self.start1
        else:
            position = self.end1 - position
        
        position = self.ali2[ self.ali1.find(position, left) ]
        
        if self.forward2:
            position = position + self.start2
        else:
            position = self.end2 - position
        
        return position

    def back_project(self, position, left=True):
        """ Positions lie between bases.
            start and end are positions. """
            
        if self.forward2:
            position = position - self.start2
        else:
            position = self.end2 - position
        
        position = self.ali1[ self.ali2.find(position,left) ]
        
        if self.forward1:
            position = position + self.start1
        else:
            position = self.end1 - position
        
        return position


class Weird_alignment(Exception): pass

def half_alignment_from_feature(feature):
    if not feature.sub_features:
        xs = 'X' * (feature.location.nofuzzy_end-feature.location.nofuzzy_start)
        return xs
    
    result = ''
    position = feature.location.nofuzzy_start
    
    for sub_feature in feature.sub_features:
        assert sub_feature.strand == feature.strand 

        step = sub_feature.location.nofuzzy_start - position
        if step < 0:
            raise Weird_alignment()

        result += '-' * step
        
        result += half_alignment_from_feature(sub_feature)
        
        position = sub_feature.location.nofuzzy_end
    
    return result
        
def alignment_from_feature(seq, feature):
    hal2 = half_alignment_from_feature(feature)
    hal1 = 'X' * len(hal2)
    assert len(hal2) == feature.location.nofuzzy_end-feature.location.nofuzzy_start
    
    forward = feature.strand != -1
    if not forward: hal2 = hal2[::-1]
    return Alignment(feature.location.nofuzzy_start,feature.location.nofuzzy_end,forward,hal1,
                     0, hal2.count('X'),True, hal2)

class Half: pass

def strip_impure_insertions(half1, half2):
    """ Impure insertions are given in lowercase by nesoni consensus """
    new_half1 = [ ]
    new_half2 = [ ]
    for i in xrange(len(half1)):
        if not (half1[i] == '-' and half2[i].islower()):
            new_half1.append(half1[i])
            new_half2.append(half2[i])
    return ''.join(new_half1), ''.join(new_half2)

def load_alignments(filename):
    f = open(filename, 'rU')
    line = f.readline()
    
    result = [ ]
    while True:        
        while line and (line.startswith('#') or not line.strip()):
            line = f.readline()
        if not line: break
        
        assert line.startswith('a')
        
        line = f.readline()
        lines = [ ]
        while line.startswith('s'):
            lines.append(line)
            line = f.readline()
        
        assert len(lines) == 2, 'Only expecting two sequences in alignment'

        halves = [ ]
        for item in lines:
            half = Half()
            halves.append(half)
                    
            parts = item.rstrip().split()
            half.name = parts[1]
            half.ali = parts[6]
            assert int(parts[2]) == 0
            assert parts[4] == '+'
            assert int(parts[3]) == int(parts[5])
        
        halves[0].ali, halves[1].ali = strip_impure_insertions(halves[0].ali, halves[1].ali)
        
        for half in halves:
            half.seq = half.ali.replace('-','')
        
        result.append((
            halves[0].name,
            halves[0].seq, 
            halves[1].seq,
            Alignment(
                0,len(halves[0].seq),True, halves[0].ali,
                0,len(halves[1].seq),True, halves[1].ali
            )
        )) 
    
    f.close()
    
    return result    
        


def get_graph(path, name, suffix):
    filename = os.path.join(path, grace.filesystem_friendly_name(name) + '-' + suffix + '.userplot')
    
    result = [ ]
    for item in open(filename,'rb'):
        result.append( float(item.strip()) )
    return numpy.array(result)



# Triangular sliding window
def windower(thing, max_radius):
    thing_pad = numpy.concatenate((
        thing[-max_radius:], thing, thing[:max_radius]
        ))
    thing_sum = numpy.cumsum(numpy.cumsum(thing_pad))
    
    return (len(thing), thing_sum, max_radius) 

def use_windower(windower, window_radius):
    len_thing, thing_sum, max_radius = windower
    return (-2*thing_sum[max_radius:len_thing+max_radius]
            +thing_sum[max_radius+window_radius:][:len_thing] 
            +thing_sum[max_radius-window_radius:][:len_thing]) / float(window_radius**2)

def expected_depth(name, seq, depths, ambig_depths):
    med = numpy.median(depths)
    sane = numpy.arange(len(depths))[ (depths > med*0.5) & (depths < med*2.0) & (depths*2.0 >= ambig_depths)]
    #print 'median', med, 'using', len(sane)
    
    if sum(sane) < 100:
        warn('Skipping depth correction on ' + name)
        return numpy.array( [numpy.average(depths)] * len(depths) ) 
    
    buckets = { }
    
    radius = 2 # examine 5-mers
    n = radius*2+1
    
    sseqq = seq[len(seq)-radius:] + seq + seq[:radius]
    for i in sane:
        s = sseqq[i:i+n]
        if s not in buckets: buckets[s] = [ ]
        buckets[s].append( depths[i] )
    
    # Pool with reverse complement
    new_buckets = { }
    for kmer in buckets:
        rc = bio.reverse_complement(kmer)
        new_buckets[kmer] = buckets[kmer] + buckets.get(rc,[])
    buckets = new_buckets
    
    for key in buckets:
        buckets[key] = numpy.average(buckets[key])
    
    prediction = numpy.zeros(len(seq), 'float')
    for i in xrange(len(seq)):
        s = sseqq[i:i+n]
        prediction[i] = buckets.get(s,0.0)
    
    
    # selection of radii from 8 to 4096
    # TODO: make this configurable, or perhaps just larger
    radii = [ int(2**(0.5*i)) for i in xrange(3*2,12*2+1) ]
    
    prediction_windower = windower(prediction, radii[-1]) 
    
    a = numpy.arange(len(seq)) / float(len(seq))
    predictors = numpy.transpose(
    [   numpy.ones(len(seq), 'float'),
        numpy.cos(a * (2.0*numpy.pi)), 
        numpy.sin(a * (2.0*numpy.pi)),
    ] + [
        use_windower(prediction_windower, radius)
        for radius in radii
    ]
    )
    
    x = linalg.lstsq(predictors[sane], depths[sane])[0]
    #print x
    prediction = numpy.sum(predictors * x[None,:], 1)
    return prediction


def graphlet(data, expect):
    n = min(len(data),10)
    result = ''
    for i in xrange(n):
        value = numpy.average(data[i*len(data)//n:(i+1)*len(data)//n])
        value_expect = numpy.average(expect[i*len(data)//n:(i+1)*len(data)//n])
        
        if value_expect <= 0.0:
            result += '?'
        else:
            value_norm = value / value_expect    
            if value_norm <= 0.0:
                result += '_'
            else:
                value_norm = int(value_norm+0.5)
                if value_norm > 9:
                    result += '!'
                else:
                    result += str(value_norm) 
        
        #x = int(value*3.5)
        #
        #if value <= 0.0:
        #    result += unichr(0x2594)
        #elif value >= 2.0:
        #    result += unichr(0x2592)
        #else:
        #    result += unichr(0x2581+max(0,min(7,x)))
    return result #( unichr(0x2595) + result + unichr(0x258f) ).encode('utf-8') 

       
        

USAGE = """
Usage:

    nesoni consequences: [options] genbank_file.gbk[.gz] working_dir
    
Options:

    --use-coverage       - Use depth of coverage information.

    --coverage-cutoff    - If using --use-coverage, proportion of gene with 
                           strange depth required to be flagged as possible 
                           deletion or duplication. Strange depth means
                           depth > 1.5 or ambiguous depth < 0.5
                           Default: 0.1

    --transl_table NN    - Translation table to use, 
                           will be overridden by /transl_table qualifier if 
                           present. If not specified, /transl_table qualifier 
                           in each CDS is mandatory.
                           "--transl_table 1" is safe for NCBI GenBank files.
                           Default: 11 (bacterial/archaeal/plant)
    
    --tabular            - Spreadsheet friendly output.
    
    --noheader           - No title row.
    
    --verbose            - List amino acid changes, rather than just counting 
                           them.
    
    
    --band               - Width of banded needleman-wunsch for protein alignment,
                           default 20, which should be plenty.
        
"""

def main(args):
    default_transl_table, args = grace.get_option_value(args, '--transl_table', int, 11)
    use_coverage, args = grace.get_flag(args, '--use-coverage')
    coverage_cutoff, args = grace.get_option_value(args, '--coverage-cutoff', float, 0.1)
    tabular, args = grace.get_flag(args, '--tabular')
    noheader, args = grace.get_flag(args, '--noheader')
    verbose, args = grace.get_flag(args, '--verbose')
    bandwidth, args = grace.get_option_value(args, '--band', int, 20)
    grace.expect_no_further_options(args)

    if len(args) != 2:
        print USAGE
        return 1
    
    genbank_filename = args[0]
    alignment_filename = args[1]
    
    if os.path.isdir(alignment_filename):
        alignment_filename = os.path.join(alignment_filename, 'alignment.maf')
    
    working_dir = os.path.split(alignment_filename)[0]
    
    alignments = load_alignments(alignment_filename)
    
    summaries = [ ]
    details = [ ]
    
    if not noheader:
        fields = 'Sequence\tLocus tag\tOld length (aa)\tNew length (aa)\tAmino acid changes\t'
        if use_coverage: fields += 'Unambiguous coverage vs expected\t\tAmbiguous coverage vs expected\t\tAmbiguous percent with any hits\t'
        fields += 'Gene\tProduct'
        if tabular: fields += '\tChanges of note'
        print fields
    
    for record in SeqIO.parse(io.open_possibly_compressed_file(genbank_filename),'genbank'):
        sequence = record.seq.tostring()
    
        for name, seq1, seq2, alignment in alignments:
            if seq1 == sequence: break
        else:
            raise grace.Error('Genbank record %s sequence not identical to any reference sequence' % record.id)
             
        if use_coverage:       
            depth = get_graph(working_dir, name, 'depth')
            ambiguous_depth = get_graph(working_dir, name, 'ambiguous-depth')
            median_depth = numpy.median(depth)
            median_ambiguous_depth = numpy.median(ambiguous_depth)
            ambiguous_factor = float(median_ambiguous_depth) / median_depth
            depth_expect = expected_depth(name, sequence, depth, ambiguous_depth)
            
        
        for feature in record.features:
            if feature.type != 'CDS': continue
            
            if 'locus_tag' not in feature.qualifiers:
                locus_tag = '%d..%d' % (feature.location.nofuzzy_start+1,feature.location.nofuzzy_end)
            else:
                locus_tag = feature.qualifiers['locus_tag'][0]
            
            if 'transl_table' in feature.qualifiers:
                transl_table_no = int(feature.qualifiers['transl_table'][0])
            else:
                assert default_transl_table is not None, 'No /transl_table for CDS, and default transl_table not given'
                transl_table_no = default_transl_table
            
            transl_table = CodonTable.ambiguous_dna_by_id[transl_table_no]
            start_codons = transl_table.start_codons
            
            try:
                feature_alignment = alignment_from_feature(sequence, feature)
            except Weird_alignment:
                warn('%s has a location I could not handle, skipping, sorry' % locus_tag)
                continue
            
            dna = [ ]
            new_dna = [ ]
            shifts = [ ]
            for i in xrange(feature_alignment.end2):
                p1 = feature_alignment.back_project(i, left=False)
                p2 = feature_alignment.back_project(i+1, left=True)
                assert abs(p2-p1) < 2
                dna.append( sequence_slice(sequence,p1,p2) )
                
                p1a = alignment.project(p1, left=False)
                p2a = alignment.project(p2, left=False) #Hmm
                
                diff = (p2-p1)-(p2a-p1a)
                #if diff:
                #    if diff%3:
                #        frame_shift = True
                #    else:
                #        frame_preserving_shift = True
                new_dna.append( sequence_slice(seq2,p1a,p2a) )
                
                if diff:
                    shifts.append((i,dna[-1],new_dna[-1]))
                
            dna = ''.join(dna)
            new_dna = ''.join(new_dna)
            
            # This usually indicated a CDS truncated at the start?
            # in which case, will probably fail some way or other down the line.
            if 'codon_start' in feature.qualifiers:
                codon_start = int(feature.qualifiers['codon_start'][0]) - 1
            else:
                codon_start = 0
            dna = dna[codon_start:]
            new_dna = new_dna[codon_start:]
            
            if len(dna) % 3 != 0:
                warn(locus_tag + ' length not a multiple of 3')
            #assert len(new_dna) % 3 == 0
            
            protein = Seq.Seq(dna).translate(table=transl_table_no).tostring()            
            # http://en.wikipedia.org/wiki/Start_codon is always translated to M
            protein = 'M' + protein[1:]
            
            if dna[:3] not in start_codons:
                warn(locus_tag + ' has unknown start codon: ' + dna[:3])
                                    
            original_lacks_stop_codon = not protein.endswith('*')                 
            if original_lacks_stop_codon:
                warn(locus_tag + ' lacks end codon')
            original_stops_before_end = '*' in protein[:-1] 
            if original_stops_before_end:
                warn(locus_tag + ' contains stop codon before end')
                            
            if 'translation' in feature.qualifiers:
                expect = feature.qualifiers['translation'][0]
                if protein[:-1] != expect:
                    warn(locus_tag + ' translation given in feature does not match translation from DNA')                
        
            new_protein = Seq.Seq(new_dna).translate(table=transl_table_no).tostring()            
            new_protein = 'M' + new_protein[1:]
        
            # If end codon changed, find new end                
            # Don't bother if there are unknown amino acids or 
            # the original protein lacks a stop codon
            if 'X' not in new_protein and '*' not in new_protein and not original_lacks_stop_codon:
                #This is very inefficient
                i = feature_alignment.end2
                while True:
                    p1 = feature_alignment.back_project(i, left=False)
                    p2 = feature_alignment.back_project(i+1, left=True)
                    p1a = alignment.project(p1, left=False)
                    p2a = alignment.project(p2, left=False) #Hmm
                    if p1a < 0 or p2a < 0 or p1a > len(seq2) or p2a > len(seq2):
                        break
                        
                    new_dna += sequence_slice(seq2,p1a,p2a)                        
                    new_protein = Seq.Seq(new_dna).translate(table=transl_table_no).tostring()            
                    new_protein = 'M' + new_protein[1:]
                    if 'X' in new_protein or '*' in new_protein: break
                    
                    i += 1
            
            # Is the protein shorter?
            # Don't bother checking if the original protein has extra stop codons
            if '*' in new_protein and not original_stops_before_end:
                new_protein = new_protein[:new_protein.index('*')+1] 
        
            # If indels occurred, do an alignment
            # Don't bother otherwise
            if shifts:
                # Penalize gaps with cost 2 (vs 1 for mismatch)
                # If lengths don't match, pad with spaces (won't match longer seq),
                # aligner prefers mismatch to gaps
                
                #result = pairwise2.align.globalxs(protein      + ' '*max(0,len(new_protein)-len(protein)), 
                #                                  new_protein  + ' '*max(0,len(protein)-len(new_protein)), 
                #                                  -2.001,-2.000)[0]
                # 2.001 : very slightly prefer contiguous gaps. Also much faster!
        
                result = band_limited_align(protein      + ' '*max(0,len(new_protein)-len(protein)), 
                                            new_protein  + ' '*max(0,len(protein)-len(new_protein)), 
                                            bandwidth)
                
                
                protein_ali = result[0]
                new_protein_ali = result[1]
            else:
                protein_ali = protein
                new_protein_ali = new_protein
        
            diffs = [ ]
            j = 0
            k = 0
            for i in xrange(min(len(new_protein_ali),len(protein_ali))):
                if protein_ali[i] != ' ' and new_protein_ali[i] != ' ' and (
                      protein_ali[i] == '-' or 
                      new_protein_ali[i] == '-' or 
                      not bio.might_be_same_amino(protein_ali[i], new_protein_ali[i]) ):
                    diffs.append((i,j,k))
                if protein_ali[i] != '-': 
                    j += 1
                if new_protein_ali[i] != '-': 
                    k += 1
        
            diff_start = not bio.might_be_same_base(new_dna[0],dna[0]) or \
                         not bio.might_be_same_base(new_dna[1],dna[1]) or \
                         not bio.might_be_same_base(new_dna[2],dna[2]) 
        
            interesting_coverage = False
            if use_coverage:
                cds_depth = depth[feature_alignment.start1:feature_alignment.end1] #/ median_depth
                if not feature_alignment.forward1: cds_depth = cds_depth[::-1]
                cds_ambiguous_depth = ambiguous_depth[feature_alignment.start1:feature_alignment.end1] #/ median_ambiguous_depth
                if not feature_alignment.forward1: cds_ambiguous_depth = cds_ambiguous_depth[::-1]
                
                cds_depth_expect = depth_expect[feature_alignment.start1:feature_alignment.end1]
                if not feature_alignment.forward1: cds_depth_expect = cds_depth_expect[::-1]
                
                #cds_average_depth_ratio = numpy.average(depth[feature_alignment.start1:feature_alignment.end1]) / median_depth 
                #cds_average_ambiguous_depth_ratio = numpy.average(ambiguous_depth[feature_alignment.start1:feature_alignment.end1]) / median_ambiguous_depth                        
                #line += '%.1f\t' % cds_average_depth_ratio 
                #line += '%.1f\t' % cds_average_ambiguous_depth_ratio
                
                #line += '%.1f..%.1f\t' % (numpy.minimum.reduce(cds_depth)/median_depth, numpy.maximum.reduce(cds_depth)/median_depth) 
                #line += '%.1f+/-%.1f\t' % (numpy.average(cds_depth)/median_depth, numpy.var(cds_depth)**0.5/median_depth) 
                #line += '%.1f..%.1f\t' % (numpy.minimum.reduce(cds_ambiguous_depth)/median_ambiguous_depth, numpy.maximum.reduce(cds_ambiguous_depth)/median_ambiguous_depth)
                
                avg_expect = numpy.average(cds_depth_expect)
                if avg_expect > 0.0:
                    cds_avg_depth = numpy.average(cds_depth)/avg_expect
                    cds_avg_ambiguous_depth = numpy.average(cds_ambiguous_depth)/avg_expect/ambiguous_factor
                
                strange = (
                    (cds_depth >= cds_depth_expect*1.5) |
                    (cds_ambiguous_depth <= cds_depth_expect*(0.5*ambiguous_factor))
                )
                
                interesting_coverage = numpy.average(strange) >= coverage_cutoff
                     

            if interesting_coverage or diffs or diff_start or shifts or len(new_protein) != len(protein):
                line = name + '\t' + locus_tag + '\t' + \
                      '%d\t' % (len(protein)-1) + \
                      '%d\t' % (len(new_protein)-1) + \
                      '%d\t' % len(diffs)
                

                if use_coverage:
                    if avg_expect <= 0.0:
                        line += '\t\t\t'
                    else:
                        line += '%.1f\t' % (cds_avg_depth) + graphlet(cds_depth, cds_depth_expect)+'\t' 
                        line += '%.1f\t' % (cds_avg_ambiguous_depth) + graphlet(cds_ambiguous_depth, cds_depth_expect*ambiguous_factor)+'\t'
                        line += '%.1f%%\t' % (numpy.average(cds_ambiguous_depth > 0.0)*100.0)
                
                line += '%s\t' % feature.qualifiers.get('gene',[''])[0] + \
                        '%s' % feature.qualifiers.get('product',[''])[0]
                
                notes = [ ]
                
                if use_coverage and 'X' in new_protein:
                    xs = new_protein.count('X')
                    if xs == len(new_protein)-1: #First is M, so len-1
                        notes.append('\ No consensus')
                    else:
                        notes.append('\ No consensus for %d aa' % (new_protein.count('X')))
                                   
                if len(new_protein) < len(protein):
                    notes.append('\ Shorter by %d aa' % (len(protein)-len(new_protein)))
        
                if len(new_protein) > len(protein):
                    notes.append('\ Longer by %d aa' % (len(new_protein)-len(protein)))
                
                if diff_start:
                    notes.append('\ Start changed: %s -> %s' % (dna[:3], new_dna[:3]))
                    if new_dna[:3] not in start_codons:
                        notes.append('  No longer a start codon!')
                        
                if shifts:
                    notes.append('\ Indels:')
                
                    for pos, old, new in shifts:
                        notes.append('    base %5d / codon %5d   %s -> %s' % (pos+1,(pos//3)+1,old,new or '-'))
                    
                if diffs:
                    if verbose:
                        notes.append('\ Amino acid changes:')
                        for i, j, k in diffs:
                            notes.append('    codon %5d   %s->%s   (%s->%s)' % (
                                j+1, 
                                protein_ali[i], 
                                new_protein_ali[i], 
                                dna[j*3:j*3+3] if protein_ali[i] != '-' else '-', 
                                new_dna[k*3:k*3+3] if new_protein_ali[i] != '-' else '-'
                            ))
                
                #if len(new_protein) > len(protein):
                #    print 'New protein is longer:', new_protein[len(protein):]
                #if len(new_protein) < len(protein):
                #    print 'New protein is shorter:', protein[len(new_protein):]
                #print protein
                #print new_protein
                
                if tabular:
                    print line + '\t' + ' '.join([ ' '.join(note.strip().split()) for note in notes ])
                else:
                    print line
                    for note in notes:
                        print '\t' + note
    return 0
    


