"""

Utilities for dealing with sequence strings

"""

import string


AMBIGUOUS_AMINO_DECODE = {
    'X' : 'ACDEFGHIKLMNPQRSTVWY*UO',
    'B' : 'DN',
    'Z' : 'EQ',
    'J' : 'LI',
}
for letter in 'ACDEFGHIKLMNPQRSTVWY*UO': AMBIGUOUS_AMINO_DECODE[letter] = letter
del letter

def might_be_same_amino(ambiguous_amino_1, ambiguous_amino_2):
    if ambiguous_amino_1 == ambiguous_amino_2: return True

    for amino1 in AMBIGUOUS_AMINO_DECODE[ambiguous_amino_1.upper()]:
        if amino1 in AMBIGUOUS_AMINO_DECODE[ambiguous_amino_2.upper()]:
            return True
    return False


COMPLEMENTER = string.maketrans('ABCDGHKMRSTVWYabcdghkmrstvwy', 'TVGHCDMKYSABWRtvghcdmkysabwr')
def reverse_complement(seq):
    """ Reverse complement,
        works with IUPAC ambiguity codes. """

    return seq.translate(COMPLEMENTER)[::-1]


AMBIGUITY_CODES = {
    '-' : '-',
    'A' : 'A',
    'T' : 'T',
    'G' : 'G',
    'C' : 'C',
    'GT' : 'K',
    'AC' : 'M',
    'AG' : 'R',
    'CT' : 'Y',
    'CG' : 'S',
    'AT' : 'W',
    'CGT' : 'B',
    'ACG' : 'V',
    'ACT' : 'H',
    'AGT' : 'D',
}

AMBIGUITY_DECODE = { }
for decode in AMBIGUITY_CODES:
    AMBIGUITY_DECODE[ AMBIGUITY_CODES[decode] ] = decode
del decode

def might_be_same_base(ambiguous_base_1, ambiguous_base_2):
    if ambiguous_base_1 == ambiguous_base_2: return True

    ambiguous_base_1 = ambiguous_base_1.upper()
    ambiguous_base_2 = ambiguous_base_2.upper()
    if ambiguous_base_1 == 'N' or ambiguous_base_2 == 'N': 
        return True
    for base1 in AMBIGUITY_DECODE[ambiguous_base_1]:
        if base1 in AMBIGUITY_DECODE[ambiguous_base_2]:
            return True
    return False

def might_be_same_bases(string1, string2):
    if len(string1) != len(string2): 
        return False
    for i in xrange(len(string1)):
        if not might_be_same_base(string1[i],string2[i]):
            return False
    return True


def represent_evidence(counts):
    counts = counts.items()
    counts.sort(key=lambda x:x[1], reverse=True)
    return ' '.join([ '"%s"x%d' % item for item in counts ])

def consensus(counts, min_depth,min_purity):
    """ Call a consensus, if it meets minimum depth and purity 
        (proportion of total) requirements. """

    if not counts: return None

    total = sum(counts.values())
    counts = counts.items()
    counts.sort(key=lambda x:x[1])
    
    if counts[-1][1] < min_depth or \
       counts[-1][1] < min_purity * total: 
        return None
    
    return counts[-1][0]

def ambiguity_code_consensus(counts, min_depth,min_purity):
    """ Call a consensus, if it meets minimum depth and purity 
        (proportion of total) requirements. 
        
        Allow DNA ambiguity codes. 
        
        """

    if not counts: return None

    total = sum(counts.values())
    counts = counts.items()
    counts.sort(key=lambda x:x[1], reverse=True)
    
    total_count = 0
    bases = [ ]
    cutoff = None #<- if multiple items have the same count as that required to reach purity/depth
                  #   be sure to include them all in the ambiguity code
    for base, count in counts:
        if cutoff is not None and count < cutoff:
            break
    
        total_count += count
        bases.append(base)
        
        if total_count >= min_depth and \
           total_count >= min_purity * total: 
            cutoff = count
    
    if cutoff is None:
        return None

    bases.sort()
    bases = ''.join(bases)
    
    return AMBIGUITY_CODES.get(bases, None)

