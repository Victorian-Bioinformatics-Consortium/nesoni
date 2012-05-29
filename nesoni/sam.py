
"""

SAM-based reboot

"""

import sys, os, subprocess, itertools, array, datetime, socket


from nesoni import grace, bio, io, consensus, legion

FLAG_PAIRED = 1
FLAG_PROPER = 2
FLAG_UNMAPPED = 4
FLAG_MATE_UNMAPPED = 8
FLAG_REVERSE = 16
FLAG_MATE_REVERSE = 32
FLAG_FIRST = 64
FLAG_SECOND = 128
FLAG_NONPRIMARY = 256
FLAG_FAIL = 512 #Read failed platform/vendor quality checks
FLAG_DUP = 1024 #PCR or optical duplicate

class Alignment(object):
    #__slots__ = [
    #    'qname',
    #    'flag',
    #    'rname',
    #    'pos',
    #    'mapq',
    #    'cigar',
    #    'mrnm',
    #    'mpos',
    #    'isize',
    #    'seq',
    #    'qual',
    #    'extra',
    #    'length',
    #]
    
    #cdef public str qname
    #cdef public int flag
    #cdef public str rname
    #cdef public int pos
    #cdef public int mapq
    #cdef public str cigar
    #cdef public str mrnm
    #cdef public int mpos
    #cdef public int isize
    #cdef public str seq
    #cdef public str qual
    #cdef public list extra
    #cdef public int length

    def __init__(self, line):
        parts = line.rstrip('\n').split('\t')
        self.extra = parts[11:]
        (self.qname, 
         flag, 
         self.rname, 
         pos, 
         mapq, 
         self.cigar, 
         self.mrnm, 
         mpos, 
         isize, 
         self.seq, 
         self.qual) = parts[:11]
        
        self.flag = int(flag)
        self.pos = int(pos) # 1-based
        self.mapq = int(mapq)
        self.mpos = int(mpos)
        self.isize = int(isize)
                
        self.length = 0
        n = 0
        for value in array.array('B', self.cigar):
            if 48 <= value <= 57:
                n = n*10+(value-48)
            else:
                #if char in 'MDNP=X':
                if value == 77 or value == 68 or value == 78 or value == 80 or value == 61 or value == 88:
                    self.length += n
                n = 0
    
    def original_name(self):
        #Assuming it was Illumina
        if self.flag&FLAG_PAIRED:
            if self.flag&FLAG_FIRST:
                return self.qname+'/1'
            elif self.flag&FLAG_SECOND:
                return self.qname+'/2'
            else:
                raise grace.Error('Confused by SAM file')                
        else:
            return self.qname

    def get_qual(self):
        if self.qual == '*':
            return '^' * len(self.seq)
        else:
            return self.qual

    def __repr__(self):
        return '%s @ %d..%d %s %s' % (self.original_name(),self.pos,self.pos+self.length-1, '-' if self.flag&FLAG_REVERSE else '+',self.rname,) 

    def set_flag(self, mask, value):
        if value:
            self.flag |= mask
        else:
            self.flag &= ~mask

    def get_mrnm(self):
        if self.mrnm == '=':
            return self.rname
        else:
            return self.mrnm
    
    def get_AS(self):
        """ Get alignment score.
            Prefer AS if present. Meaning of mapq is complicated. """
        for item in self.extra:
            if item.startswith('AS:i:'):
                return int(item[5:])
        
        return self.mapq        
        #raise grace.Error('SAM line lacks AS')
    
    #def get_length(self):
    #    length = 0
    #    i = 0
    #    n = 0
    #    for c in self.cigar:
    #        if c.isdigit():
    #            n = n*10 + ord(c)-48
    #        else:
    #            if c in 'MDNP': #TODO: checkme
    #                length += n
    #            n = 0
    #    return length


def is_bam(filename):
    import gzip
    try:
        f = gzip.open(filename,'rb')
        magic = f.read(4)
        f.close()
    except IOError:
        return False
    return magic == 'BAM\x01'

def open_bam(filename):
    process = io.run([
        'samtools', 'view', '-h',
        io.abspath(filename),
    ])    
    return process.stdout



def bam_headers(filename):
    process = io.run([
        'samtools',
        'view',
        '-H',
        io.abspath(filename),
    ])
    
    headers = process.stdout.read()                
    assert process.wait() == 0, '"samtools view -H ..." failed'
    
    return headers


class Bam_reader(object):
    def __init__(self, filename):
        assert os.path.exists(filename), filename + ' does not exist'
        
        if is_bam(filename):
            self.process = io.run([
                'samtools',
                'view',
                io.abspath(filename),
            ])
            
            ## Godawful hack
            #self.process.stdout = io.process_buffer(self.process.stdout)
            self.file = self.process.stdout
        else:
            self.process = None
            self.file = io.open_possibly_compressed_file(filename)
    
    def __iter__(self):
        return self
        
    def next(self):
        line = self.file.readline()
        if not line:
            if self.process is not None:
                assert self.process.wait() == 0, '"samtools view ..." failed'
            raise StopIteration()
        
        return Alignment(line)


def bam_iter_fragments(filename, status_text='Processing'):
    reader = Bam_reader(filename)

    n = 0
    n_ambiguous = 0
    for read_name, alignment_iter in itertools.groupby(reader, lambda read: read.qname):
        if n % 100000 == 0:
            grace.status(status_text + ' fragment %s' % grace.pretty_number(n))
        n += 1
        
        unpaired = [ ]
        first = [ ]
        second = [ ]
        unmapped = [ ]
        for al in alignment_iter:
            if al.flag&FLAG_UNMAPPED:
                unmapped.append(al)            
            elif not al.flag&FLAG_PAIRED or al.flag&FLAG_MATE_UNMAPPED:
                unpaired.append((al,))
            elif al.flag&FLAG_FIRST:
                first.append(al)
            elif al.flag&FLAG_SECOND:
                second.append(al)
            else:
                assert False, 'Read in pair that is neither first nor second'
        
        pairs = [ ]
        
        unused = set(first + second)
        
        second_index = { }
        for al in second:
            key = (al.rname, al.pos)
            if key not in second_index: second_index[key] = [ ]
            second_index[key].append(al)
        
        for al1 in first:
            key = (al1.get_mrnm(), al1.mpos)
            for al2 in second_index.get(key, ()):
                if al2.get_mrnm() != al1.rname or \
                   al2.mpos != al1.pos:
                    continue
                
                if al1 not in unused or al2 not in unused:
                    #TODO: this is fucked
                    n_ambiguous += 1
                    continue
                
                pairs.append( (al1, al2) )
                unused.remove(al1)
                unused.remove(al2)
        
        assert not unused, 'Alignment pairing not even pretending to make sense'
        
        yield read_name, pairs + unpaired, unmapped
    
    grace.status('')
    if n_ambiguous:
        print >> sys.stderr 
        print >> sys.stderr, 'The alignment pairing was unclear %s times, and alignments were paired arbitrarily.' % grace.pretty_number(n_ambiguous)
        print >> sys.stderr, 'Blame the SAM format and/or SHRiMP.'
        print >> sys.stderr 

        
class Bam_writer(object):
    def __init__(self, filename, headers=''):
        self.writer = io.Pipe_writer(filename, ['samtools', 'view', '-S', '-b', '-'])
        self.writer.write(headers) 

    def write(self, al):
        line = '\t'.join([
                al.qname,
                str(al.flag),
                al.rname,
                str(al.pos),
                str(al.mapq),
                al.cigar,
                al.mrnm,
                str(al.mpos),
                str(al.isize),
                al.seq,
                al.qual
            ] + al.extra)
        self.writer.write(line + '\n')
    
    def write_raw(self, text):
        self.writer.write(text)
    
    def close(self):
        self.writer.close()




def sort_and_index(in_filename, out_prefix):
    io.execute([
        'samtools', 'sort', in_filename, out_prefix
    ])
    
    io.execute([
        'samtools', 'index', out_prefix + '.bam'
    ])



