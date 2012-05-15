"""

SHRiMP read consensus caller

"""


import sys, os, heapq, string, math, struct, cPickle, gzip

import nesoni
from nesoni import bio, io, grace, shrimp, statistics

#USAGE = """\
#
#Usage: 
#
#  nesoni consensus: working_dir [options]
#
#where working_dir was created by "nesoni shrimp".
#
#Consensus calling, short version: 
#   Adjust --cutoff to obtain reasonable depth cutoffs.
#   You should be able to leave --indel-prior and --prior-weight alone.
#
#Consensus calling, long version:   
#   A priori we believe each position in the reference to be some mixture of 
#   A, C, G, T, and deletion, but we are unsure what exact mixture it is.
#   This prior belief is expressed as a Dirichlet distribution. 
#   As we observe bases from various reads, our beliefs are updated.
#   A posteriori, if we believe with probability >= cutoff that there 
#   is a certain base (or a deletion) that forms an absolute majority in the 
#   mixture, it is called as the consensus.
#   
#   For insertions, a priori we believe there to be a mixture of insertion and 
#   no-insertion, and the same process occurs. 
#
#"""
#
#RECONSENSUS_USAGE = """\
#
#Usage:
#
#  nesoni reconsensus: working_dir [options]
#
#
#"""
#
#COMPLEMENTER = string.maketrans('ACGTacgt','TGCAtgca')
#
#def reverse_complement(seq):
#    """ Reverse complement of a string, works with ACGTN- """
#    
#    return seq.translate(COMPLEMENTER)[::-1]
#
#def bisect_left(a, key_func, x):
#    """Return the index where to insert item x in list a, assuming a is sorted.
#
#    The return value i is such that all e in a[:i] have e < x, and all e in
#    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
#    insert just before the leftmost x already there.
#
#    Optional args lo (default 0) and hi (default len(a)) bound the
#    slice of a to be searched.
#    """
#
#    lo = 0
#    hi = len(a)
#    while lo < hi:
#        mid = (lo+hi)//2
#        if key_func(a[mid]) < x: lo = mid+1
#        else: hi = mid
#    return lo
#    
#def edit_string_to_alignment(edit_string, corresp_seq):
#    """ Convert an edit string as output by SHRiMP and ELAND to
#        a conventional alignment. """
#
#    len_edit_string = len(edit_string)
#    read_ali = ''
#    contig_ali = ''
#    i = 0    
#    while i < len_edit_string:
#        if edit_string[i] == '(':
#            j = i + 1
#            while edit_string[j] != ')': j += 1
#            gap = edit_string[i+1:j].upper()
#            
#            #print 'insertion', gap
#            read_ali += gap
#            contig_ali += '-' * len(gap)
#            
#            i = j + 1
#            
#        elif edit_string[i].isdigit():
#            j = i
#            while j < len_edit_string and edit_string[j].isdigit(): j += 1
#            n_matches = int(edit_string[i:j])
#
#            #print 'match', n_matches
#            assert n_matches <= len(corresp_seq), "Bad edit string"   
#            read_ali += corresp_seq[:n_matches]
#            contig_ali += corresp_seq[:n_matches]
#            corresp_seq = corresp_seq[n_matches:]
#
#            i = j           
#             
#        elif edit_string[i] == '-':
#            #print 'deletion'
#            j = i
#            while j < len_edit_string and edit_string[j] == '-': j += 1
#            
#            n_deleted = j-i
#            
#            assert n_deleted <= len(corresp_seq), "Malformed edit string"
#            
#            contig_ali += corresp_seq[:n_deleted]
#            read_ali += '-' * n_deleted
#            corresp_seq = corresp_seq[n_deleted:]
#        
#            i = j
#            
#        elif edit_string[i] == 'x':
#            #Crossover in colorspace
#            i += 1
#            
#        else:
#            #print 'substitution', edit_string[i]
#            assert corresp_seq, "Strange edit string"
#            
#            contig_ali += corresp_seq[:1]
#            read_ali += edit_string[i].upper()
#            corresp_seq = corresp_seq[1:]
#            
#            i += 1
#
#    assert not corresp_seq, "Peculiar edit string"
#
#    return contig_ali, read_ali
#
#def roll_alignment(ali, ali_other):
#    """ Normalize one half of an alignment by shifting "-"s as far right as possible.
#    
#        (SHRiMP alignments can differ depending on whether the read is forward 
#         or reverse, when the location of an insertion or deletion is ambiguous.)
#    """
#    
#    if '-' not in ali: return ali
#
#    len_ali = len(ali)
#
#    i = 0
#    while i < len_ali:
#        if ali[i] != '-': 
#            i += 1
#            continue
#            
#        j = i+1
#        while j < len_ali and ali[j] == '-': j += 1
#        if j >= len_ali: 
#            break
#            
#        if ali_other[i] == ali_other[j] or ali[j] != ali_other[j]:
#            ali = ali[:i] + ali[j] + ali[i+1:j] + '-' + ali[j+1:]            
#            i = i + 1
#        else:
#            i = j
#            
#    return ali
#
#def _get0(x): return x[0]
#def _get1(x): return x[1]

EMPTY_EVIDENCE = ()

def evidence_add(evidence, item, n):
    pos = 0
    while pos < len(evidence) and evidence[pos] != item:
        pos += 2
    
    if pos < len(evidence):
        result = evidence[:pos+1] + (evidence[pos+1]+n,) + evidence[pos+2:]
    else:
        result = evidence + (item,n,)
    
    while pos and result[pos-1] < result[pos+1]:
        result = result[:pos-2] + result[pos:pos+2] + result[pos-2:pos] + result[pos+2:]
        pos -= 2
        
    return result

def evidence_merge(evidence1, evidence2):
    """ Could be more efficient """
    for pos in xrange(0,len(evidence2),2):
        evidence1 = evidence_add(evidence1, evidence2[pos], evidence2[pos+1])
    return evidence1

def evidence_total_count(evidence):
    return sum(evidence[1::2])

#
#def consensus(counts, min_depth,min_purity):
#    """ Call a consensus, if it meets minimum depth and purity 
#        (proportion of total) requirements. """
#
#    if not counts: return None
#
#    total = sum(counts[1::2])
#    if counts[1] < min_depth or \
#       counts[1] < min_purity * total:
#        return None
#    
#    return counts[0]
#


def bayesian_consensus(counts, p_cutoff,nonempty_prior,empty_prior,total_prior, proportion=0.5):
    """ Call a consensus, if it meets required probability of absolute majority. """

    if not counts: return None, 1.0

    total = sum(counts[1::2])
    
    if counts[0] == '-':
        x = empty_prior
    else:
        x = nonempty_prior 
    x += counts[1]
    
    p = statistics.probability_of_proportion_at_least(x, total+total_prior, proportion)
    if p < p_cutoff: return None, 1.0
    
    return counts[0], p


#AMBIGUITY_CODES = {
#    '-' : '-',
#    'A' : 'A',
#    'T' : 'T',
#    'G' : 'G',
#    'C' : 'C',
#    'GT' : 'K',
#    'AC' : 'M',
#    'AG' : 'R',
#    'CT' : 'Y',
#    'CG' : 'S',
#    'AT' : 'W',
#    'CGT' : 'B',
#    'ACG' : 'V',
#    'ACT' : 'H',
#    'AGT' : 'D',
#}
#
#AMBIGUITY_DECODE = { }
#for decode in AMBIGUITY_CODES:
#    AMBIGUITY_DECODE[ AMBIGUITY_CODES[decode] ] = decode
#
#def ambiguity_merge(*codes):
#    if None in codes: return None
#    bases = ''.join( AMBIGUITY_DECODE[code] for code in codes )
#    bases = ''.join(sorted(set(bases)))
#    return AMBIGUITY_CODES.get(bases, None)
#
#def ambiguity_code_consensus(evidence, min_depth,min_purity):
#    """ Call a consensus, if it meets minimum depth and purity 
#        (proportion of total) requirements. 
#        
#        Allow DNA ambiguity codes. 
#        
#        """
#
#    if not evidence: return None
#    
#    total = sum(evidence[1::2])
#    
#    total_count = 0
#    bases = [ ]
#    cutoff = None
#    
#    pos = 0
#    while pos < len(evidence):
#        if cutoff is not None and evidence[pos+1] < cutoff:
#            break
#    
#        total_count += evidence[pos+1]
#        bases.append( evidence[pos] )
#        
#        if total_count >= min_depth and \
#           total_count >= min_purity * total:
#            cutoff = evidence[pos+1]
#    
#        pos += 2
#    
#    if cutoff is None:
#        return None
#    
#    bases.sort()
#    bases = ''.join(bases)
#    return bio.AMBIGUITY_CODES.get(bases, None)

def bayesian_ambiguity_code_consensus(evidence, p_cutoff,nonempty_prior,empty_prior,total_prior, proportion=0.5):
    if not evidence: return None, 1.0
    
    total = sum(evidence[1::2])
    
    x = 0.0
    bases = [ ]
    cutoff = None
    
    pos = 0
    while pos < len(evidence):
        if cutoff is not None and evidence[pos+1] < cutoff:
            break
    
        if evidence[pos] == '-':
            x += empty_prior
        else:
            x += nonempty_prior
        x += evidence[pos+1]        
        
        bases.append( evidence[pos] )
        
        p = statistics.probability_of_proportion_at_least(x, total+total_prior, proportion)        
        if p >= p_cutoff:
            cutoff = evidence[pos+1]
    
        pos += 2
    
    if cutoff is None:
        return None, 1.0
    
    bases.sort()
    bases = ''.join(bases)
    return bio.AMBIGUITY_CODES.get(bases, None), p



    #total = sum(counts.values())
    #counts = counts.items()
    #counts.sort(key=_get1, reverse=True)
    #
    #total_count = 0
    #bases = [ ]
    #cutoff = None #<- if multiple items have the same count as that required to reach purity/depth
    #              #   be sure to include them all in the ambiguity code
    #for base, count in counts:
    #    if cutoff is not None and count < cutoff:
    #        break
    #
    #    total_count += count
    #    bases.append(base)
    #    
    #    if total_count >= min_depth and \
    #       total_count >= min_purity * total: 
    #        cutoff = count
    #
    #if cutoff is None:
    #    return None
    #
    #bases.sort()
    #bases = ''.join(bases)
    #
    #return AMBIGUITY_CODES.get(bases, None)

#def entropy(counts):
#    total = sum(counts.values())
#    result = 0.0
#    for value in counts.values():
#        if value == 0: continue
#        p = float(value) / total
#        result -= p * math.log(p)
#    return result / math.log(2.0)


#def iter_increasing_pairs(list1, list2, scorer=lambda a,b: a+b):
#    """ Yield pairs from list1 and list2, 
#        in increasing order of combined score.
#    
#        lists must be sorted and score combiner
#        must be monotonic increasing in both parameters. """
#
#    if not list1 or not list2: return
#    
#    heap = [ (scorer(list1[0],list2[0]),0,0) ]
#    while heap:
#        score, pos1, pos2 = heapq.heappop(heap)        
#        yield score, pos1, pos2
#        
#        if not pos2 and pos1+1 < len(list1):
#            heapq.heappush(heap,
#                (scorer(list1[pos1+1],list2[pos2]), pos1+1, pos2) )
#        if pos2+1 < len(list2):
#            heapq.heappush(heap,
#                (scorer(list1[pos1],list2[pos2+1]), pos1, pos2+1) )

#class Hit_pair_iter(object):
#    """ Yield pairs of hits from list1 and list2, 
#        in decreasing order of combined score.
#    
#        lists must be sorted. """
#        
#    def __init__(self, list1, list2):
#        if not list1 or not list2:
#            self.heap = [ ]
#        else:
#            self.heap = [ (-list1[0].score-list2[0].score,0,0) ]
#        self.list1 = list1
#        self.list2 = list2
#    
#    def __iter__(self):
#        return self
#    
#    def next(self):
#        if not self.heap:
#            raise StopIteration()
#    
#        neg_score, pos1, pos2 = heapq.heappop(self.heap)        
#        
#        if not pos2 and pos1+1 < len(self.list1):
#            heapq.heappush(self.heap,
#                (-self.list1[pos1+1].score-self.list2[pos2].score, pos1+1, pos2) )
#        if pos2+1 < len(self.list2):
#            heapq.heappush(self.heap,
#                (-self.list1[pos1].score-self.list2[pos2+1].score, pos1, pos2+1) )
#
#        return -neg_score, pos1, pos2


def pretty_evidence(evidence):
    return ' '.join([
        '%sx%d' % evidence[pos:pos+2]
        for pos in xrange(0,len(evidence),2)
    ])
    
    #counts = counts.items()
    #counts.sort(key=_get1, reverse=True)
    #return ' '.join([ '"%s"x%d' % item for item in counts ])


#class Picklable(object):
#    def __getstate__(self):
#        result = { }
#        for name in self.__class__.ATTRIBUTES:
#            result[name] = getattr(self, name)
#        return result
#
#    def __setstate__(self, state):
#        for name in self.__class__.ATTRIBUTES:
#            setattr(self, name, state[name])
#
#class Refseqset(Picklable):
#    """ A collection of reference seqeunces, to which reads are aligned.
#    
#        Properties:
#    
#            seqs - { name : Refseq }
#            hits - { read_name : [ Hit ] }
#
#    """
#    
#    #cdef public object seqs, seq_order, used_hit_file, unambiguous_seps, shrimp_threshold
#    #cdef public int infidelity
#    #cdef public bool is_paired_end, same_direction, is_monogamous, is_circular, only_pairs, whole_read_only
#    #cdef public int max_pair_sep, trim
#    #
#    #cdef public int n_single, n_single_with_hits, n_single_unambiguous, n_pair, n_pair_with_hits, n_pair_valid, n_pair_valid_unambiguous
#    
#    ATTRIBUTES = [ 'seqs', 'seq_order' ] # Only those relevant to reconsensusing
#
#    def __init__(self, is_paired_end, same_direction, max_pair_sep, trim, shrimp_threshold, infidelity, is_monogamous, is_circular, only_pairs, whole_read_only, used_hit_file):
#        self.seqs = { } 
#        self.seq_order = [ ]       
#        
#        self.is_paired_end = is_paired_end
#        self.same_direction = same_direction
#        self.max_pair_sep = max_pair_sep
#        self.trim = trim
#        self.shrimp_threshold = shrimp_threshold
#        self.infidelity = infidelity
#        self.is_monogamous = is_monogamous
#        self.is_circular = is_circular
#        self.only_pairs = only_pairs
#        self.whole_read_only = whole_read_only
#        self.used_hit_file = used_hit_file
#        
#        self.n_single = 0
#        self.n_single_with_hits = 0
#        self.n_single_unambiguous = 0
#        self.n_pair = 0
#        self.n_pair_with_hits = 0
#        self.n_pair_valid = 0
#        self.n_pair_valid_unambiguous = 0
#        
#        self.unambiguous_seps = [ ]
#        
#    def add_sequence(self, name, sequence):
#        assert name not in self.seqs, 'Duplicate reference sequence name'
#        
#        self.seqs[name] = Refseq(name, sequence)
#        self.seq_order.append(name)
#    
#    def run(self, directory):
#        reader = shrimp.iter_read_hits(directory)
#        
#        if not self.is_paired_end:
#            warn_about_paired_end = True
#        
#            for ((read_name, read_seq), lines) in reader:
#                if (self.n_single % 10000) == 0: 
#                    grace.status('Processing read %s' % grace.pretty_number(self.n_single))
#                self.n_single += 1
#                
#                if warn_about_paired_end and \
#                   (read_name[-3:] in ('_F3','_R3') or
#                    read_name[-2:] in ('/1','/2')):
#                    grace.status('')
#                
#                    sys.stderr.write('\n*** WARNING: read names look like paired end, but no --max-pair-sep given ***\n\n')
#                    warn_about_paired_end = False 
#            
#                hits = [ self.build_hit(line) for line in lines ]
#                
#                if self.whole_read_only:
#                    hits = filter(self.hit_uses_whole_read, hits)
#                
#                self.process_hit_set(read_name, read_seq, hits)
#        
#        else:
#            warn_about_uninterleaved = True
#        
#            while True:
#                try:
#                    ((read_name_1, read_seq_1), lines_1) = reader.next()
#                    ((read_name_2, read_seq_2), lines_2) = reader.next()
#                except StopIteration:
#                    break
#                
#                if warn_about_uninterleaved and read_name_1[-1:] == read_name_2[-1:]:
#                    sys.stderr.write('\n*** WARNING: alleged read pairs have same suffix ***\n')
#                    sys.stderr.write('Pairs should be passed to "nesoni shrimp" in a "pairs" section\n')
#                    sys.stderr.write('or interleaved before being given in a "reads" section.\n\n')
#                    warn_about_uninterleaved = False
#                
#                if (self.n_pair % 10000) == 0: grace.status('Processing pair %s' % grace.pretty_number(self.n_pair))
#                self.n_pair += 1
#                
#                hits_1 = [ self.build_hit(line) for line in lines_1 ]
#                hits_2 = [ self.build_hit(line) for line in lines_2 ]
#                
#                if self.whole_read_only:
#                    hits_1 = filter(self.hit_uses_whole_read, hits_1)
#                    hits_2 = filter(self.hit_uses_whole_read, hits_2)
#                
#                self.process_paired_hit_set(read_name_1, read_seq_1, hits_1, read_name_2, read_seq_2, hits_2)
#    
#    def hit_uses_whole_read(self, hit):
#        # But allow partial match at start/end of references
#        ref_length = len(self.seqs[hit.ref_name].reference)
#        return (hit.read_start == 0 or (hit.forward and hit.ref_start == 0) or (not hit.forward and hit.ref_end == ref_length)) and \
#               (hit.read_end == hit.read_length or (hit.forward and hit.ref_end == ref_length) or (not hit.forward and hit.ref_start == 0))
#    
#    def stats_text(self):
#        import numpy
#    
#        text = ''
#        if not self.n_single and not self.n_pair:
#            text += '\nNo reads.\n'
#            return text
#        
#        if self.n_pair:
#            text += '%20s read pairs\n' % grace.pretty_number(self.n_pair)
#        
#        if self.n_single:
#            text += '%20s ' % grace.pretty_number(self.n_single)
#            if self.n_pair: text += 'further '
#            text += 'single reads\n'
#        
#        if self.is_paired_end:
#            text += '\n%20s pairs where both reads hit something\n' % grace.pretty_number(self.n_pair_with_hits)
#            text += '%20s pairs validly oriented and spaced\n' % grace.pretty_number(self.n_pair_valid)
#            if self.is_monogamous:
#                text += '%20s unambiguously\n' % grace.pretty_number(self.n_pair_valid_unambiguous)
#            else:
#                text += '(%19s unambiguously, however all sufficiently good hits were used)\n' % grace.pretty_number(self.n_pair_valid_unambiguous)
#
#        if self.unambiguous_seps:
#            text += grace.pretty_number(int(numpy.median(self.unambiguous_seps)),20) + \
#                    ' median insert size (limit %d)\n' % self.max_pair_sep
#
#        
#        if not self.only_pairs or self.n_single or self.n_single_with_hits:
#            text += '\n%20s reads hit something on their own\n' % grace.pretty_number(self.n_single_with_hits)
#            if self.is_monogamous:
#                text += '%20s unambiguously\n' % grace.pretty_number(self.n_single_unambiguous)
#            else:
#                text += '(%19s unambiguously, however all sufficiently good hits were used)\n' % grace.pretty_number(self.n_single_unambiguous)
#        
#        return text
#        
#    
#    def build_hit(self, line):
#        result = [ ]
#
#        #(0 read_name, 1 contig_name, 2 strand, 
#        # 3 contig_start, 4 contig_end, 
#        # 5 read_start, 6 read_end, 7 read_length, 
#        # 8 score, 9 edit_string) = line.rstrip().split('\t')[:10]
#        
#        parts = line.rstrip('\n').split('\t')
#        
#        hit = Hit()
#        
#        hit.line = line
#        
#        hit.read_name = parts[0]
#        hit.ref_name = parts[1]
#        hit.ref_start = int(parts[3])-1
#        hit.ref_end = int(parts[4])
#        hit.read_start = int(parts[5])-1
#        hit.read_end = int(parts[6])
#        hit.read_length = int(parts[7])
#        hit.score = int(parts[8])
#        hit.forward = (parts[2] == '+')
#        hit.edit_string = parts[9]            
#        
#        return hit
#        
#        #if self.paired_end and not self.same_direction and hit.read_name.endswith(self.suffix_2):
#        #    hit.source_strand = 1 if hit.forward else 0
#        #else:
#        #    hit.source_strand = 0 if hit.forward else 1
#    
#            
#
#    def process_hit_set(self, read_name, read_seq, hits):
#        """ Update counts in Refseq objects based on hits read in. """
#        
#        if not hits: return
#        
#        actual_threshold = shrimp.calc_threshold(len(read_seq), self.shrimp_threshold)
#        
#        for hit in hits:
#            if hit.forward:
#                hit.source_strand = 0
#            else:
#                hit.source_strand = 1
#        
#        self.n_single_with_hits += 1
#        
#        hits.sort(key=_hit_score, reverse=True)
#        
#        n = 1
#        while n < len(hits) and hits[n].score >= hits[0].score-self.infidelity:
#            n += 1
#            
#        if n == 1:
#            self.n_single_unambiguous += 1
#            
#        if self.is_monogamous:
#            if n == 1 and \
#               hits[0].score >= actual_threshold+self.infidelity: #Don't allow infidelity to be hidden by threshold
#                self.seqs[ hits[0].ref_name ].process_hit(hits[0], self.trim, self.used_hit_file)                
#            
#            for hit in hits[:n]:
#                self.seqs[ hit.ref_name ].process_ambiguous_hit(hit, self.trim)
#        else:
#            for hit in hits[:n]:
#                self.seqs[ hit.ref_name ].process_hit(hit, self.trim, self.used_hit_file)
#
#        #stats_text = (
#        #   grace.pretty_number(n_total,20) + ' reads hit something\n' +
#        #   grace.pretty_number(n_unambiguous,20) + ' unambiguously\n'
#        #)
#
#        #grace.status('')        
#        #return stats_text
#    
#    def process_paired_hit_set(self, read_name_1, read_seq_1, hits_1, read_name_2, read_seq_2, hits_2):
#        """ Update counts in Refseq objects based on hits read in. """
#        import numpy
#        
#        actual_threshold = shrimp.calc_threshold(len(read_seq_1), self.shrimp_threshold) + \
#                           shrimp.calc_threshold(len(read_seq_2), self.shrimp_threshold) 
#        
#        #weird_pair_report = [ ]
#        
#        if not hits_1 and not hits_2:
#            return 
#        
#        if not hits_2:
#            if not self.only_pairs:
#                self.process_hit_set(read_name_1, read_seq_1, hits_1)
#            return
#        
#        if not hits_1:
#            if not self.only_pairs:
#                self.process_hit_set(read_name_2, read_seq_2, hits_2)
#            return
#
#        self.n_pair_with_hits += 1
#
#
#        for hit in hits_1:
#            if hit.forward:
#                hit.source_strand = 0
#            else:
#                hit.source_strand = 1
#        if self.same_direction:
#            for hit in hits_2:
#                if hit.forward:
#                    hit.source_strand = 0
#                else:
#                    hit.source_strand = 1
#        else:
#            for hit in hits_2:
#                if hit.forward:
#                    hit.source_strand = 1
#                else:
#                    hit.source_strand = 0
#            
#        #hits_1.sort(key=_hit_score, reverse=True)
#        #hits_2.sort(key=_hit_score, reverse=True)
#                        
#        top_hits = valid_pairs(hits_1, hits_2, self.same_direction, self.max_pair_sep, self.is_circular, self.seqs)
#        j = 0
#        while j < len(top_hits) and top_hits[j][0] >= top_hits[0][0]-self.infidelity:
#            j += 1
#        top_hits = top_hits[:j]
#            
#        if not top_hits:
#            #if (len(hits_1) == 1 or hits_1[0].score*self.fidelity > hits_1[1].score) and \
#            #   (len(hits_2) == 1 or hits_2[0].score*self.fidelity > hits_2[1].score):
#            #    hit_1 = hits_1[0]
#            #    hit_2 = hits_2[0]
#            #    if hit_1.forward:
#            #        hit_1_end = hit_1.ref_end
#            #    else:
#            #        hit_1_end = hit_1.ref_start
#            #    if hit_2.forward:
#            #        hit_2_end = hit_2.ref_end
#            #    else:
#            #        hit_2_end = hit_2.ref_start
#            #    #out by one error above?
#            #    
#            #    weird_pair_report.append((hit_1.ref_name, hit_1.forward, hit_1_end, hit_2.ref_name, hit_2.forward, hit_2_end))
#            #
#            return
#        
#        self.n_pair_valid += 1
#        
#        if len(top_hits) == 1:
#            self.n_pair_valid_unambiguous += 1
#            
#        if self.is_monogamous:
#            if len(top_hits) == 1 and \
#               top_hits[0][0] >= actual_threshold+self.infidelity: #Don't allow infidelity to be hidden by threshold
#                hits_to_use = top_hits
#            else:
#                hits_to_use = [ ]
#                 
#            for score, hit_1, hit_2 in top_hits:
#                self.seqs[ hit_1.ref_name ].process_ambiguous_hit(hit_1, self.trim)
#                self.seqs[ hit_2.ref_name ].process_ambiguous_hit(hit_2, self.trim)
#            
#        else:
#            hits_to_use = top_hits
#        
#        for score, hit_1, hit_2 in hits_to_use:
#            hit_1 = top_hits[0][1]
#            hit_2 = top_hits[0][2]
#            if hit_1.forward:
#                #sep = hit_2.ref_start - hit_1.ref_end
#                sep = hit_2.ref_end - hit_1.ref_start
#                span_start = hit_1.ref_start
#                span_end = hit_2.ref_end
#            else:
#                #sep = hit_1.ref_start - hit_2.ref_end
#                sep = hit_1.ref_end - hit_2.ref_start
#                span_start = hit_2.ref_start
#                span_end = hit_1.ref_end
#                
#            if len(top_hits) == 1:
#                self.unambiguous_seps.append(sep)
#            
#            if span_start <= span_end: #Not wrapped
#                span = self.seqs[ hit_1.ref_name ].depth_pairspan[hit_1.source_strand][span_start:span_end]
#                numpy.add(span,1,span)
#            #TODO: handle wrapping at end for circular genomes
#
#            self.seqs[ hit_1.ref_name ].process_hit( hit_1, self.trim, self.used_hit_file )
#        for score, hit_1, hit_2 in hits_to_use:
#            self.seqs[ hit_2.ref_name ].process_hit( hit_2, self.trim, self.used_hit_file )
#        
##
##        grace.status('')
##        
##        if not only_pairs:
##            single_text = ''
##            self.process_hits(trim, fidelity, is_monogamous, used_hit_file, orphans + unpaired)
##        else:
##            single_text = ' (ignored)'
##        
##        stats_text = (        
##           grace.pretty_number(n_total,20) + ' read pairs where both reads hit something\n' +
##           grace.pretty_number(n_valid_total,20) + ' read pairs validly oriented and spaced\n' +
##           grace.pretty_number(n_valid_unambiguous,20) + ' unambiguously\n'
##        )
##        
##        if unambiguous_seps:
##            stats_text += grace.pretty_number(int(numpy.median(unambiguous_seps)),20) + \
##                          ' median insert size (limit %d)\n' % max_pair_sep
##
##        stats_text += '\n' + grace.pretty_number(len(orphans),20) + ' reads with no hits to their pair%s\n' % single_text
##        
##        if unpaired:
##            stats_text += '\n' + grace.pretty_number(len(unpaired),20) + ' without a read-pair suffix, treated as single reads%s\n' % single_text
##        
##        return stats_text, weird_pair_report
#
#class Hit(object):
#    """ Structure to store hits """
#    pass
#
#    #__slots__ = (
#    #    'ref_name',
#    #    'read_name',
#    #    'score',
#    #    'forward',
#    #    'ref_start',
#    #    'ref_end',
#    #    'ref_ali',
#    #    'read_start',
#    #    'read_end',
#    #    'read_ali',
#    #    'read_length',
#    #)
#    
#    #cdef public object line
#    #
#    #cdef public object ref_name
#    #cdef public object read_name
#    #cdef public int score
#    #cdef public int forward
#    #cdef public int ref_start
#    #cdef public int ref_end
#    ##cdef public object ref_ali
#    #cdef public int read_start
#    #cdef public int read_end
#    ##cdef public object read_ali
#    #cdef public int read_length
#    #
#    #cdef public object edit_string
#    #
#    #cdef public int source_strand    # strand of *pair*, 0=forward, 1=reverse, for depth plots
#
#def _hit_score(hit): return hit.score
#def _hit_start_location(hit): return (hit.ref_name, hit.ref_start)
#def _hit_end_location(hit): return (hit.ref_name, hit.ref_end)
#    
#def valid_pairs(hits1, hits2, same_direction, max_pair_sep, circular_reference, seqs):
#    hits2_forward = [ ]
#    hits2_reverse = [ ]
#        
#    for hit in hits2:
#        if hit.forward:
#            hits2_forward.append(hit)
#        else:
#            hits2_reverse.append(hit)
#
#    if same_direction:
#        hits2_for_hit1_forward = hits2_forward
#        hits2_for_hit1_reverse = hits2_reverse
#    else:
#        hits2_for_hit1_forward = hits2_reverse
#        hits2_for_hit1_reverse = hits2_forward
#    
#    hits2_for_hit1_forward.sort(key=_hit_end_location)
#    hits2_for_hit1_reverse.sort(key=_hit_start_location)
#    
#    result = [ ] # (score, hit1, hit2)
#    
#    for hit1 in hits1:
#        if hit1.forward:
#            options = hits2_for_hit1_forward[
#                bisect_left(hits2_for_hit1_forward, _hit_end_location, (hit1.ref_name, hit1.ref_start)) :
#                bisect_left(hits2_for_hit1_forward, _hit_end_location, (hit1.ref_name, hit1.ref_start + max_pair_sep+1))
#            ]
#            if circular_reference:
#                ref_length = len(seqs[hit1.ref_name].reference)
#                options.extend(hits2_for_hit1_forward[
#                    bisect_left(hits2_for_hit1_forward, _hit_end_location, (hit1.ref_name, hit1.ref_start-ref_length)) :
#                    bisect_left(hits2_for_hit1_forward, _hit_end_location, (hit1.ref_name, hit1.ref_start-ref_length+max_pair_sep+1)) :
#                ])
#        else:
#            options = hits2_for_hit1_reverse[
#                bisect_left(hits2_for_hit1_reverse, _hit_start_location, (hit1.ref_name, hit1.ref_end -max_pair_sep)) :
#                bisect_left(hits2_for_hit1_reverse, _hit_start_location, (hit1.ref_name, hit1.ref_end +1))
#            ]
#            if circular_reference:
#                ref_length = len(seqs[hit1.ref_name].reference)
#                options.extend(hits2_for_hit1_reverse[
#                    bisect_left(hits2_for_hit1_reverse, _hit_start_location, (hit1.ref_name, hit1.ref_end+ref_length -max_pair_sep)) :
#                    bisect_left(hits2_for_hit1_reverse, _hit_start_location, (hit1.ref_name, hit1.ref_end+ref_length +1))
#                ])
#    
#        for hit2 in options:
#            result.append( (hit1.score+hit2.score, hit1, hit2) )
#    
#    result.sort(key=_get0, reverse=True)
#    return result
#
#    
#class Refseq(Picklable):
#    """ Reference sequence and statistics on alignments thereto.
#    
#        Properties:
#
#            name
#            reference - reference sequence
#        
#            depth - unambiguous hit depth
#            depth_ambiguous - depth counting all top-scoring hits
#        
#            base_counts - { letter : count } (includes deletions as '-')
#            insertions - { position : [ sequence ] }
#                
#    """
#    
#    #cdef public object name, reference, depth, depth_ambiguous, depth_pairspan, base_counts, insertions
#    
#    ATTRIBUTES = [ 'name', 'reference', 'base_counts', 'insertions', 'depth' ] # Only those relevant to reconsensusing
#
#    def __init__(self, name, reference):
#        import numpy
#    
#        self.name = name
#        self.reference = reference
#
#        #Keep stranded depth statistics 0 == forward, 1 == reverse
#        self.depth = [ numpy.zeros(len(self.reference), 'int') for strand in (0,1) ]
#        self.depth_ambiguous = [ numpy.zeros(len(self.reference), 'int') for strand in (0,1) ]
#        
#        self.depth_pairspan = [ numpy.zeros(len(self.reference), 'int') for strand in (0,1) ] 
#        
#        #self.incomplete_ends_count = numpy.zeros(len(self.reference), 'int')
#        #self.bias = numpy.zeros(len(self.reference), 'float')
#
#        #self.base_counts = { }        
#        #for base in 'ACGT-':
#        #    self.base_counts[base] = numpy.zeros(len(self.reference), 'int')
#            
#        #self.insertions = { } # index *after* insertion -> { inserted sequence : count }
#        
#        #self.base_counts = [ {} for i in xrange(len(reference)) ]
#        #self.insertions = [ {} for i in xrange(len(reference)) ]
#        self.base_counts = [ EMPTY_EVIDENCE ] * len(reference)
#        self.insertions = [ EMPTY_EVIDENCE ] * len(reference)
#
#    def process_hit(self, hit, trim, used_hit_file):
#        if used_hit_file is not None:
#            used_hit_file.write(hit.line)
#
#        # Convert edit string to alignment
#        
#        corresp_seq = self.reference[hit.ref_start:hit.ref_end]
#        
#        if not hit.forward:
#            corresp_seq = reverse_complement(corresp_seq)
#        
#        hit_ref_ali, hit_read_ali = edit_string_to_alignment(hit.edit_string, corresp_seq)	    
#        
#        if not hit.forward:
#            hit_ref_ali = reverse_complement(hit_ref_ali)
#            hit_read_ali = reverse_complement(hit_read_ali)
#        
#        #Normalization -- move "-"s as far right as possible
#        hit_read_ali = roll_alignment(hit_read_ali, hit_ref_ali)
#        hit_ref_ali = roll_alignment(hit_ref_ali, hit_read_ali)
#        
#        i = 0
#        position = hit.ref_start
#        
#        #hit_ref_ali = hit.ref_ali
#        #hit_read_ali = hit.read_ali
#        len_hit_ref_ali = len(hit_ref_ali)
#        
#        self_insertions = self.insertions
#        self_base_counts = self.base_counts
#        
#        self_depth = self.depth[hit.source_strand]
#        
#        if position != 0: #Don't trim at start of reference! 
#            while i < trim:
#                if hit_ref_ali[i] != '-':
#                    position += 1
#                i += 1
#
#        if hit.ref_end == len(self.reference): #Don't trim at end of reference!
#            trim = 0
#        
#        while i < len_hit_ref_ali-trim:
#            if hit_ref_ali[i] == '-':
#                j = i + 1
#                while j < len_hit_ref_ali and hit_ref_ali[j] == '-': j += 1
#                what = hit_read_ali[i:j]
#                
#                #if position not in self.insertions: self.insertions[position] = { }                
#                
#                #counter = self_insertions[position]
#                #counter[what] = counter.get(what,0) + 1
#                self_insertions[position] = evidence_add(self_insertions[position], what, 1)
#                
#                i = j	        
#            else:
#                base = hit_read_ali[i]
#                if base != 'N' and base != 'X': #X is a shrimp colorspace confuzzle
#                    #self.base_counts[ hit_read_ali[i] ][position] += 1
#                    
#                    #counter = self_base_counts[position]
#                    #counter[base] = counter.get(base,0)+1
#                    
#                    self_base_counts[position] = evidence_add(self_base_counts[position], base, 1)
#                    
#                self_depth[position] += 1
#                
#                position += 1        
#                i += 1
#    
#    def process_ambiguous_hit(self, hit, trim):
#        #TODO: don't ignore trim!
#        self.depth_ambiguous[hit.source_strand][hit.ref_start:hit.ref_end] = \
#            self.depth_ambiguous[hit.source_strand][hit.ref_start:hit.ref_end] + 1
#
#    def consensus(self, p_cutoff,indel_prior,prior_weight,use_ambiguity_codes):
#        """ 
#        Returns:        
#          consensus sequence        
#          snp/indel report        
#          For each position in the reference, whether there was a consensus        
#        """
#        total_prior = prior_weight
#        
#        empty_prior    = indel_prior * total_prior
#        nonempty_prior = (total_prior-empty_prior) / 4.0 # Each base equally likely
#        
#        insertion_present_prior = indel_prior * total_prior
#        insertion_absent_prior = total_prior - insertion_present_prior
#    
#        result = [ ]
#        result_masked_only = [ ]
#        
#        alignment_result    = [ ]
#        alignment_reference = [ ]
#        
#        #insertion_evidence = [ ]
#        #substitution_evidence = [ ]
#        evidence = [ ] # [(insertion_evidence, substitution_evidence, insertion_call, substitution_call)]
#        
#        has_consensus = [ ]
#        
#        grace.status('Consensus %s' % self.name)
#        
#        report = [ ]
#        
#        for i in xrange(len(self.reference)):
#            if i % 100000 == 0:
#                grace.status('Consensus %s %s' % (self.name, grace.pretty_number(i)))
#            #if i in self.insertions:
#            # Need to count absense of insertions in consensus
#            #if i in self.insertions:
#            #insertions = self.insertions.get(i,{}).copy()
#            insertions = self.insertions[i]
#            #total = sum(insertions.values())
#            total = evidence_total_count(insertions)
#            #else:
#            #    insertions = { }
#            #    total = 0
#            depth = self.depth[0][i] + self.depth[1][i]
#            if i: depth = min(depth, self.depth[0][i-1] + self.depth[1][i-1])
#            if depth > total: 
#                #insertions['-'] = depth-total
#                insertions = evidence_add(insertions, '-', depth-total)
#                
#            insertion_call_with_purity = bayesian_consensus(insertions, p_cutoff,insertion_present_prior,insertion_absent_prior,total_prior)[0]
#            insertion_call_without_purity = bayesian_consensus(insertions, 0.0,insertion_present_prior,insertion_absent_prior,total_prior)[0]
#            
#            #Indicate lack of purity with lower case in alignment and consensus files
#            if insertion_call_without_purity is not None and insertion_call_with_purity is None:
#                insertion_call_without_purity = insertion_call_without_purity.lower()
#            
#            if insertion_call_without_purity is not None and insertion_call_without_purity != '-': 
#                result.append(insertion_call_without_purity)
#                result_masked_only.append(insertion_call_without_purity)
#                
#                alignment_result.append(insertion_call_without_purity)
#                alignment_reference.append('-' * len(insertion_call_without_purity))
#
#            if insertion_call_with_purity is not None and insertion_call_with_purity != '-':
#                report.append(('insertion-before', i, '-', insertion_call_with_purity, insertions))
#            
#            #insertion_evidence.append(insertions)
#            
#            
#            
#            #counts = { }
#            #for base in self.base_counts:
#            #    if self.base_counts[base][i]: 
#            #        counts[base] = int(self.base_counts[base][i]) #Typecast to strip weird slow printing numpy type
#            counts = self.base_counts[i]
#            
#            if use_ambiguity_codes:
#                substitution_call = bayesian_ambiguity_code_consensus(counts, p_cutoff,nonempty_prior,empty_prior,total_prior)[0]
#            else:
#                substitution_call = bayesian_consensus(counts, p_cutoff,nonempty_prior,empty_prior,total_prior)[0]
#            
#            if substitution_call is None:
#                result.append('N')
#                result_masked_only.append( self.reference[i].lower() )
#                has_consensus.append(False)
#            elif substitution_call == '-':
#                report.append(('deletion',i, self.reference[i], '-', counts))
#                has_consensus.append(True)
#            else:
#                result.append(substitution_call)
#                result_masked_only.append(substitution_call)
#                has_consensus.append(substitution_call in 'ACGT') #Exclude ambiguity codes
#                
#                #if c != self.reference[i]:
#                if self.reference[i] not in AMBIGUITY_DECODE[substitution_call]:
#                    report.append(('substitution', i, self.reference[i], substitution_call, counts))
#                    
#            #substitution_evidence.append(counts)
#            
#            evidence.append((
#                insertions,
#                counts,
#                insertion_call_with_purity or 'N',
#                substitution_call or 'N'
#            ))
#            
#            alignment_result.append(substitution_call or 'N')
#            alignment_reference.append(self.reference[i])
#            
#        
#        grace.status('')
#        
#        return (
#            ''.join(result),
#            ''.join(result_masked_only), 
#            report, 
#            has_consensus, 
#            evidence,
#            ''.join(alignment_reference), 
#            ''.join(alignment_result)
#        )


def table_to_text(table, left_pad, alignment):
    if not table: return
    
    widths = [ 0 ] * len(table[0])
    for line in table:
        for i, item in enumerate(line):
            widths[i] = max(widths[i], len(item))

    output = [ ]

    for line in table:
        out_line = left_pad 
        for i, item in enumerate(line):
            pad = ' '*(widths[i]-len(item))
            if alignment[i] == 'L':
                item = item + pad
            else:
                assert alignment[i] == 'R'
                item = pad + item
            out_line += item + '  '
        output.append( out_line.rstrip(' ') + '\n' )
    
    return ''.join(output)


def write_userplot(filename, array):
    f = open(filename,'wb')
    for x in array:
        f.write( '%d\n' % int(x) )
    f.close()

def write_userplots(prefix, arrays, output_strand_specific_depths):
    if not output_strand_specific_depths:
        write_userplot(prefix + '.userplot', arrays[0]+arrays[1])
        return
        
    #if output_strand_specific_depths:
    #    write_userplot(prefix + '-forward.userplot', arrays[0])
    #    write_userplot(prefix + '-reverse.userplot', arrays[1])
    
    f = open(prefix + '.userplot', 'wb')
    f.write('# BASE forward reverse\n')
    f.write('# colour 0:196:0 196:0:0\n')
    for i in xrange(len(arrays[0])):
        f.write('%d %d %d\n' % (i+1,int(arrays[0][i]),int(arrays[1][i])))
    f.close()



def consensus_calling_advice(p_cutoff, indel_prior, prior_weight, 
                             title='With these options, the required depth of coverage will be:\n',
                             suffix='\nErrors can occur due to miscalled bases, incorrect alignment,\nor the sample might be truly mixed.\n',
                             proportion=0.5
                             ):
    result = [ ]

    total_prior = prior_weight
    
    empty_prior    = indel_prior * total_prior
    nonempty_prior = (total_prior-empty_prior) / 4.0 # Each base equally likely
    
    insertion_present_prior = indel_prior * total_prior
    insertion_absent_prior = total_prior - insertion_present_prior

    n = 1
    while bayesian_consensus( ('A',n), p_cutoff, nonempty_prior, empty_prior, total_prior, proportion )[0] is None:
        n += 1    
    result.append([ str(n), 'to call a base with no errors' ])

    n = 1
    while bayesian_consensus( ('A',n-1,'C',1), p_cutoff, nonempty_prior, empty_prior, total_prior, proportion )[0] is None:
        n += 1    
    result.append([ str(n), 'to call a base with 1 error' ])

    n = 2
    while bayesian_consensus( ('A',n-2,'C',2), p_cutoff, nonempty_prior, empty_prior, total_prior, proportion )[0] is None:
        n += 1    
    result.append([ str(n), 'to call a base with 2 errors' ])

    n = 3
    while bayesian_consensus( ('A',n-3,'C',3), p_cutoff, nonempty_prior, empty_prior, total_prior, proportion )[0] is None:
        n += 1    
    result.append([ str(n), 'to call a base with 3 errors' ])

    n = 1
    while bayesian_consensus( ('-',n), p_cutoff, nonempty_prior, empty_prior, total_prior, proportion )[0] is None:
        n += 1    
    result.append([ str(n), 'to call a deletion with no errors' ])

    n = 1
    while bayesian_consensus( ('A',n), p_cutoff, insertion_present_prior, insertion_absent_prior, total_prior, proportion )[0] is None:
        n += 1    
    result.append([ str(n), 'to call an insertion present with no errors' ])

    n = 1
    while bayesian_consensus( ('-',n), p_cutoff, insertion_present_prior, insertion_absent_prior, total_prior, proportion )[0] is None:
        n += 1    
    result.append([ str(n), 'to call no insertion with no errors' ])
    
    return title + table_to_text(result, '    ', 'RL') + suffix
       
#
#
#def main(args):
#    p_cutoff, args = grace.get_option_value(args,'--cutoff', float, 0.99)
#    indel_prior, args = grace.get_option_value(args,'--indel-prior', float, 0.2)
#    prior_weight, args = grace.get_option_value(args,'--prior-weight', float, 1.0)
#
#    whole_read_only, args = grace.get_option_value(args, '--whole-read-only', grace.as_bool, False)
#    trim, args = grace.get_option_value(args,'--trim', int, 5)
#    infidelity, args = grace.get_option_value(args,'--infidelity', int, 50)
#    is_monogamous, args = grace.get_option_value(args,'--monogamous', grace.as_bool, True)
#    use_ambiguity_codes, args = grace.get_option_value(args,'--ambiguity-codes', grace.as_bool, True)
#    
#    max_pair_sep, args = grace.get_option_value(args, '--max-pair-sep', int, None)
#    only_pairs, args = grace.get_option_value(args, '--only-pairs', grace.as_bool, False)
#    same_direction, args = grace.get_option_value(args, '--same-dir', grace.as_bool, False)
#
#    circular_reference, args = grace.get_option_value(args, '--circular', grace.as_bool, True)
#    
#    output_strand_specific_depths, args = grace.get_option_value(args, '--strand-specific', grace.as_bool, False)
#    save_hits, args = grace.get_option_value(args, '--save-hits', grace.as_bool, False)
#    
#    grace.expect_no_further_options(args)
#
#    print 'Options:\n'
#    
#    if max_pair_sep is None:
#         max_pair_sep_text = 'not given'
#    else:
#         max_pair_sep_text = '%d' % max_pair_sep
#    table = [
#        ['--cutoff', '%.6f' % p_cutoff, '(probability of absolute majority required to'],
#        ['','',' call consensus)'],
#        
#        ['--indel-prior', '%.6f' % indel_prior, '(prior expected proportion of'],
#        ['','', ' insertion or deletion at each base)'],
#        ['--prior-weight', '%.1f' % prior_weight, '(our prior belief is as though we have'],
#        ['','',' this many existing observations)'],
#        
#        ['','',''],
#        ['--infidelity', '%d' % infidelity, '(any runner-up hits scoring this close to'],
#        ['','',' the best hit\'s score are considered possibly valid hits.'],
#        ['','', ' Note: each extra base mismatch loses 25 points)'],
#        
#        ['--monogamous', grace.describe_bool(is_monogamous), '(discard reads with more than'],
#        ['','',                ' one possibly valid hit)'],
#        ['--whole-read-only', grace.describe_bool(whole_read_only), '(only accept alignments that cover'],
#        ['','',' the whole read)'],
#        ['--trim', '%d' % trim, '(amount to trim from start/end of alignments)'],
#        ['','',''],
#        ['--ambiguity-codes', grace.describe_bool(use_ambiguity_codes), '(use IUPAC ambiguity codes)'],
#        ['','',''],
#        ['--max-pair-sep', max_pair_sep_text, '(maximum size of a read pair *including'],
#        ['', '', ' the length of the reads themselves*,'],
#        ['', '', ' will treat reads as unpaired if not given)'],
#        ['--only-pairs', grace.describe_bool(only_pairs), '(ignore reads without a hit to their pair)'],
#        ['--same-dir',  grace.describe_bool(same_direction) , '(read pairs have the same orientation)'],
#        ['','',''],
#        ['--circular', grace.describe_bool(circular_reference) , '(assume reference sequences are circular'],
#        ['', '', ' when checking for valid pairing)'],
#        ['','',''],
#        ['--strand-specific', grace.describe_bool(output_strand_specific_depths), '(output strand-specific depths)'],
#        ['--save-hits', grace.describe_bool(save_hits) , '(save file of hits ultimately used)'],
#    ]
#    table_text = table_to_text(table, '  ', 'LLL')
#    print table_text
#
#    consensus_advice = consensus_calling_advice(p_cutoff, indel_prior, prior_weight)
#    print consensus_advice    
#    
#    if len(args) != 1:
#        sys.stderr.write( USAGE )
#        return 1
#    
#    output_dir = args[0]
#
#    reference_filename = os.path.join(output_dir,'reference.fa')
#        
#    shrimp_filename = os.path.join(output_dir,'shrimp_hits.txt.gz')
#    
#    for filename in [reference_filename, shrimp_filename]:
#        assert os.path.exists(filename), filename + ' does not exist'
#    
#    if not os.path.isdir(output_dir):
#        os.mkdir(output_dir)
#    
#    if save_hits:
#        used_hit_file = gzip.open(os.path.join(output_dir, 'used_shrimp_hits.txt.gz'), 'wb')
#    else:
#        used_hit_file = None
#
#    shrimp_config = shrimp.load_config(output_dir)
#    assert 'threshold' in shrimp_config, 'You need to re-run nesoni shrimp with the latest version'
#    shrimp_threshold = shrimp_config['threshold']
#    
#    if infidelity and is_monogamous:
#        if max_pair_sep is not None:
#           print 'Hit-pairs scoring below the threshold supplied to shrimp plus %d' % infidelity
#        else:
#           print 'Hits scoring below the threshold supplied to shrimp plus %d' % infidelity
#        print 'will be discarded because we can\'t be certain there was not another'
#        print 'sufficiently good hit just below the threshhold given to SHRiMP.'
#        print
#
#    seqset = Refseqset(max_pair_sep is not None, same_direction, max_pair_sep or 0, trim, shrimp_threshold, infidelity, is_monogamous, circular_reference, only_pairs, whole_read_only, used_hit_file)
#    
#    for name, seq in io.read_fasta(reference_filename):
#        seqset.add_sequence(name, seq.upper())
#
#    seqset.run(output_dir)
#
#    stats_text = seqset.stats_text()
#
#    log_file = open(os.path.join(output_dir, 'consensus_log.txt'), 'wb')
#    log_file.write( table_text )
#    log_file.write( '\n' )
#    log_file.write( consensus_advice )
#    log_file.write( '\n' )
#    log_file.write( stats_text )
#    log_file.close()
#
#    print stats_text
#
#    if save_hits:
#        used_hit_file.close()
#        
#    grace.status('Save information for reconsensus')
#    pickle_file = gzip.open(os.path.join(output_dir, 'evidence.pickle.gz'), 'wb')
#    cPickle.dump(seqset, pickle_file, 2)
#    pickle_file.close()    
#    grace.status('')
#
#    for name in seqset.seq_order:
#        grace.status('Write depth for ' + name)
#        
#        write_userplots(
#            os.path.join(output_dir, grace.filesystem_friendly_name(name) + '-depth'),
#            seqset.seqs[name].depth,
#            output_strand_specific_depths
#        )
#        
#        if is_monogamous:
#            grace.status('Write ambiguous depth for ' + name)
#            
#            write_userplots(
#                os.path.join(output_dir, grace.filesystem_friendly_name(name) + '-ambiguous-depth'),
#                seqset.seqs[name].depth_ambiguous,
#                output_strand_specific_depths
#            )
#        
#        if max_pair_sep is not None:
#            grace.status('Write pair-span depth for ' + name)
#            
#            write_userplots(
#                os.path.join(output_dir, grace.filesystem_friendly_name(name) + '-pairspan-depth'),
#                seqset.seqs[name].depth_pairspan,
#                output_strand_specific_depths
#            )
#        
#        grace.status('')
#    
#    do_consensus(seqset, output_dir, p_cutoff, indel_prior, prior_weight, use_ambiguity_codes)
#
#def reconsensus(args):
#    p_cutoff, args = grace.get_option_value(args,'--cutoff', float, 0.99)
#    indel_prior, args = grace.get_option_value(args,'--indel-prior', float, 0.2)
#    prior_weight, args = grace.get_option_value(args,'--prior-weight', float, 1.0)
#
#    use_ambiguity_codes, args = grace.get_option_value(args,'--ambiguity-codes', grace.as_bool, True)
#
#    grace.expect_no_further_options(args)
#
#    print 'Options:\n'
#    
#    table = [
#        ['--cutoff', '%.6f' % p_cutoff],        
#        ['--indel-prior', '%.6f' % indel_prior],
#        ['--prior-weight', '%.1f' % prior_weight],
#        ['--ambiguity-codes', grace.describe_bool(use_ambiguity_codes)],
#    ]
#    table_text = table_to_text(table, '  ', 'LLL')
#    print table_text
#    consensus_advice = consensus_calling_advice(p_cutoff, indel_prior, prior_weight)
#    print consensus_advice        
#    
#    if len(args) != 1:
#        sys.stderr.write( RECONSENSUS_USAGE )
#        return 1
#    
#    output_dir = args[0]
#
#    grace.status('Load evidence')
#    pickle_file = gzip.open(os.path.join(output_dir, 'evidence.pickle.gz'), 'rb')
#    seqset = cPickle.loads(pickle_file.read())  #Note: inefficient memory use, but loading a gzip file is very slow for some reason 
#    pickle_file.close()
#    grace.status('')
#
#    log_file = open(os.path.join(output_dir, 'consensus_log.txt'), 'ab')
#    log_file.write( '\nReconsensus:\n' )
#    log_file.write( table_text )
#    log_file.write( '\n' )
#    log_file.write( consensus_advice )
#    log_file.write( '\n' )
#    log_file.close()
#     
#    do_consensus(seqset, output_dir, p_cutoff, indel_prior, prior_weight, use_ambiguity_codes)
#
#
#def do_consensus(seqset, output_dir, p_cutoff, indel_prior, prior_weight, use_ambiguity_codes):        
#    consensus_file = open(os.path.join(output_dir, 'consensus.fa'), 'wb')
#    consensus_masked_file = open(os.path.join(output_dir, 'consensus_masked.fa'), 'wb')
#    reference_having_consensus_file = open(os.path.join(output_dir, 'reference_having_consensus.fa'), 'wb')
#    report_file = open(os.path.join(output_dir, 'report.txt'), 'wb')
#    report_gff_file = open(os.path.join(output_dir, 'report.gff'), 'wb')
#    alignment_file = open(os.path.join(output_dir, 'alignment.maf'), 'wb')
#    
#    report_file.write( 'Sequence\tPosition in reference\tChange type\tOld\tNew\tEvidence\n' )
#    
#    alignment_file.write( '##maf version=1\n' )
#    alignment_file.write( '#nesoni %s\n' % nesoni.VERSION )
#
#    report_gff_file.write( '##gff-version 3\n' ) 
#    
#    for name in seqset.seq_order:
#        consensus, consensus_masked_only, report, has_consensus, evidence, \
#        alignment_reference, alignment_result = \
#            seqset.seqs[name].consensus(p_cutoff, indel_prior, prior_weight, use_ambiguity_codes)
#        
#        grace.status('Write results for ' + name)
#            
#        io.write_fasta(consensus_file, name, consensus)
#        
#        io.write_fasta(consensus_masked_file, name, consensus_masked_only)
#        
#        ref = seqset.seqs[name].reference
#        seq = ''.join([ item and ref[i] or 'n' for i, item in enumerate(has_consensus) ])
#        io.write_fasta(reference_having_consensus_file, name, seq)
#        
#        alignment_file.write( '\n' )
#        alignment_file.write( 'a'+'\n')
#        alignment_file.write( 's %s           0 %d + %d %s\n' % (name, len(ref), len(ref), alignment_reference))
#        alignment_file.write( 's %s-consensus 0 %d + %d %s\n' % (name, len(ref), len(ref), alignment_result))
#                
#        grace.status('Write report for ' + name)
#            
#        for change_type, position, old, new, counts in report:
#            report_file.write( '%s\t%d\t%s\t%s\t%s\t%s\n' % (
#                name,
#                position+1,
#                change_type,
#                old,
#                new,
#                pretty_evidence(counts)
#            ))
#            
#            if change_type == 'deletion':
#                start = position+1
#                end = position+1
#                product = 'Base deleted: %s' % (old)
#            elif change_type == 'insertion-before':
#                #Bracket insertion
#                start = position
#                end = position+1
#                product = 'Insertion: .' + new + '.'
#            else:
#                start = position+1
#                end = position+1
#                product = 'Substitution: %s became %s' % (old,new)
#            
#            product += ' ('+pretty_evidence(counts)+')'
#            
#            report_gff_file.write( '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n' % (
#                name,
#                'nesoni_consensus_'+nesoni.VERSION,
#                'variation',
#                start,
#                end,
#                '.', #score
#                '+', #strand
#                '.', #frame
#                'product='+product
#            ))
#            
#        grace.status('Write evidence for ' + name)
#            
#        f = open(os.path.join(output_dir, grace.filesystem_friendly_name(name) + '-evidence.txt'),'wb')
#        f.write( 'Position\tInsertion-before evidence\tSubstitution evidence\tReference\tInsertion call\tSubstitution call\n' )
#        for i in xrange(len(evidence)):
#            f.write( '%d\t%s\t%s\t%s\t%s\t%s\n' % (
#                i+1,
#                pretty_evidence(evidence[i][0]),
#                pretty_evidence(evidence[i][1]),
#                seqset.seqs[name].reference[i],
#                evidence[i][2],
#                evidence[i][3],
#            ))
#        
#        grace.status('')
#
#    consensus_file.close()
#    consensus_masked_file.close()
#    reference_having_consensus_file.close()
#    report_file.close()
#    report_gff_file.close()
#    alignment_file.close()
#    
#    return 0
#    