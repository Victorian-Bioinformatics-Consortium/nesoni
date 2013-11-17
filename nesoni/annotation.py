"""

Some code taken from https://github.com/chapmanb/bcbb

"""

from nesoni import grace, io

import sys, re, urllib, os, os.path

strand_to_gff = { 1:'+', -1:'-', 0:'.', None:'?' }
strand_from_gff = { '+':1, '-':-1, '.':0, '?':None }

class Annotation(object):
    """ 
    Structure to store a GFF line or BioPython Feature
    
    Note: Coordinates are zero-based! Because this is the correct way to do coordinates, damnit.
    
    Attributes:
        seqid   - chromosome name
        source  - generating program
        type    - feature type
        start   - start, zero based
        end     - end, zero based
        strand  - '+' (gff) = 1 (here), '-' = -1, '.' = 0, '?' = None
                  ie 0 is unstranded, None is unknown
                  Note: (strand or 1) will give strand with default to 1 if unknown/unstranded 
        score   - None or float
        phase   - None or 0,1,2
        attr    - { 'key':'value',... }
    """
    
    def __init__(self, seqid='', source='nesoni', type='feature', start=0, end=0, strand=0, score=None, phase=None, attr={}):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.strand = strand
        self.score = score
        self.phase = phase        
        self.attr = dict(attr)
    
    def __repr__(self):
        return '%s%s[%d,%d) %s %s' % (
             self.seqid,
             strand_to_gff[self.strand],
             self.start,
             self.end,
             self.type,
             ' '.join( key+'='+val for key,val in self.attr.items() )
        ) 

    def get_id(self):
        for key in ('ID','id','locus_tag'):
            if key in self.attr:
                return self.attr[key]

        return '%s/%d..%d' % (self.seqid,self.start+1,self.end)

    def overlaps(self, other, allowance=0, check_strand=True):
        return (
            self.seqid == other.seqid and
            (not check_strand or self.strand == other.strand) and 
            other.start < self.end+allowance and 
            self.start < other.end+allowance
        )

    def contains(self, other):
        return (
            self.seqid == other.seqid and
            self.strand == other.strand and 
            self.start <= other.start and 
            other.end <= self.end
        )

    def as_gff(self):
        return '\t'.join([
            self.seqid,
            self.source,
            self.type,
            str(self.start+1),
            str(self.end),
            '.' if self.score is None else str(self.score),
            strand_to_gff[self.strand],
            '.' if self.phase is None else str(self.phase),
            encode_keyvals(self.attr),
        ])


def link_up_annotations(annotations):
    index = { }
    for item in annnotations:
        id = item.get_id()
        assert id not in index, 'Annotations contain a duplicated ID.'
        index[id] = item
        item.children = [ ]
        
    for item in annotations:
        if 'Parent' not in item.attr:
            item.parents = [ ]
        else:
            item.parents = [ index[parent] for parent in item.attr['Parent'].split(',') ]
            for parent in item.parents:
                index[parent].children.append(item)
                
    for item in annotations:
        if item.strand == -1:
            item.children.sort(key=lambda i: i.end, reverse=True)
        else:
            item.children.sort(key=lambda i: i.start)



gff3_kw_pat = re.compile("\w+=")

def split_keyvals(keyval_str):
    """Split key-value pairs in a GFF2, GTF and GFF3 compatible way.

    GFF3 has key value pairs like:
    count=9;gene=amx-2;sequence=SAGE:aacggagccg
    GFF2 and GTF have:
    Sequence "Y74C9A" ; Note "Clone Y74C9A; Genbank AC024206"
    name "fgenesh1_pg.C_chr_1000003"; transcriptId 869
    """
    quals = { }
    if keyval_str is None:
        return quals
    
    # ensembl GTF has a stray semi-colon at the end
    if keyval_str[-1] == ';':
        keyval_str = keyval_str[:-1]
    
    # GFF2/GTF has a semi-colon with at least one space after it.
    # It can have spaces on both sides; wormbase does this.
    # GFF3 works with no spaces.
    # Split at the first one we can recognize as working
    parts = keyval_str.split(" ; ")
    if len(parts) == 1:
        parts = keyval_str.split("; ")
        if len(parts) == 1:
            parts = keyval_str.split(";")
    
    # check if we have GFF3 style key-vals (with =)
    is_gff2 = True
    if gff3_kw_pat.match(parts[0]):
        is_gff2 = False
        key_vals = [p.split('=') for p in parts]
    
    # otherwise, we are separated by a space with a key as the first item
    else:
        pieces = []
        for p in parts:
            # fix misplaced semi-colons in keys in some GFF2 files
            if p and p[0] == ';':
                p = p[1:]
            pieces.append(p.strip().split(" "))
        key_vals = [(p[0], " ".join(p[1:])) for p in pieces]
    
    for item in key_vals:
        # standard in-spec items are key=value
        if len(item) == 2:
            key, val = item

        # out-of-spec files can have just key values. We set an empty value
        # which will be changed to true later to standardize.
        else:
            assert len(item) == 1, item
            key = item[0]
            val = ''

        # remove quotes in GFF2 files
        if (len(val) > 0 and val[0] == '"' and val[-1] == '"'):
            val = val[1:-1]

        assert key not in quals
        quals[key] = val

    for key, val in quals.items():
        quals[key] = urllib.unquote(val)

    return quals


def quote(s):
    result = []
    for char in s:
        if char in '=;' or char < ' ':
            result.append('%%%02X' % ord(char))
        else:
            result.append(char)
    return ''.join(result)

def encode_keyvals(quals):
    return ';'.join( key+'='+quote(val) for key,val in quals.items() )



def read_gff(filename):
    f = io.open_possibly_compressed_file(filename)
    for line in f:
        line = line.rstrip()
        
        if line == '##FASTA':
            break
        
        if not line or line.startswith('#'): 
            continue
        
        parts = line.split('\t')
        assert len(parts) >= 8, parts

        result = Annotation()
        
        result.seqid = parts[0]
        result.source = parts[1]
        result.type = parts[2]
        result.start = int(parts[3])-1
        result.end = int(parts[4])
        result.score = None if parts[5] == '.' else float(parts[5])
        result.strand = strand_from_gff[parts[6]]
        result.phase = None if parts[7] == '.' else int(parts[7])
        result.attr = { } if len(parts) < 9 else split_keyvals(parts[8])
        
        yield result        

    f.close()


def read_genbank(filename):
    from Bio import Seq, SeqIO
    f = io.open_possibly_compressed_file(filename)
    
    id_counter = 0
    
    for record in SeqIO.parse(f,'genbank'):
        name = record.id
        if name == '' or name == 'unknown':
            name = record.name

        for root_feature in record.features:
            todo = [ root_feature ]
            while todo:
                feature = todo.pop()            
                result = Annotation()
                result.seqid = name
                result.source = 'genbank-file'
                result.type = feature.type
                result.start = feature.location.nofuzzy_start
                result.end = feature.location.nofuzzy_end
                result.score = None
                result.strand = feature.strand
                result.phase = 0 #FIXME
                result.attr = { }
                for key in feature.qualifiers:
                    result.attr[key] = ', '.join(feature.qualifiers[key])
                yield result
                
                if 'ID' not in result.attr:
                    id_counter += 1
                    result.attr['ID'] = '%d' % id_counter                 
                for sub_feature in feature.sub_features:
                    feature.qualifiers['Parent'] = [ result.attr['ID'] ]
                
                todo.extend(feature.sub_features[::-1])
                
    f.close()


def read_annotations(filename):
    f = io.open_possibly_compressed_file(filename)
    peek = f.read(1024)
    f.close()
    
    if peek.startswith('LOCUS'):
        return read_genbank(filename)
    elif peek.startswith('##gff') or peek.split('\n')[0].count('\t') in (7,8):
        return read_gff(filename)
    else:
        raise grace.Error('Not an annotation file.')


def is_annotation_file(filename):
    if not os.path.isfile(filename): 
        return False
    try:
        read_annotations(filename)
        return True
    except grace.Error:
        return False


def write_gff3_header(f):    
    print >> f, '##gff-version 3'




