"""

Read and write various types of file / directory

"""

from nesoni import grace, legion

import os, cPickle, sys, re, weakref, collections, csv, subprocess, gzip, bz2, itertools

STREAM_PROCESS = weakref.WeakKeyDictionary()
def close(f):
    """ Violently close pipes """
    if f in STREAM_PROCESS:
        STREAM_PROCESS[f].kill()
    f.close()



class Workspace(object):
    """ Directory containing pickled objects, etc """

    def __init__(self, working_dir, must_exist=False):
        self.working_dir = os.path.normpath(working_dir)
        if not os.path.exists(self.working_dir):
            assert not must_exist, working_dir + ' does not exist'
            os.mkdir(self.working_dir)
        else:
            assert os.path.isdir(self.working_dir), self.working_dir + ' exists and is not a directory'
        
        self.name = os.path.split(os.path.abspath(working_dir))[1]

        self._load_param()        

    def _load_param(self):
        if self.object_exists('parameters'):
            self.param = self.get_object('parameters', plain_text=True)
        else:
            self.param = { }

    def update_param(self, remove=[], **updates):
        self._load_param()
        for item in remove:
            if item in self.param:
                del self.param[item]
        self.param.update(updates)
        self.set_object(self.param, 'parameters', plain_text=True)

    def open(self, path, mode):
        return open(self.object_filename(path), mode)
    
    def object_exists(self, path):
        return os.path.exists(self.object_filename(path))

    def get_object(self, path, plain_text=False):
        f = open_possibly_compressed_file(self._object_filename(path))
        if plain_text:
            result = eval(f.read())
        else:
            result = cPickle.load(f)
        f.close()
        return result

    def set_object(self, obj, path, plain_text=False):
        temp_filename = self._object_filename('tempfile')
        if plain_text:
            f = open(temp_filename, 'wb')
            f.write(repr(obj))
            f.close()
        else:
            f = Pipe_writer(temp_filename, ['gzip'])
            cPickle.dump(obj, f, 2)
            f.close()
        
        os.rename(temp_filename, self._object_filename(path))

    def path_as_relative_path(self, path):        
        if os.path.isabs(path):
            return path

        assert os.path.sep == '/' #Someone else can work this out on windows and mac   
        me = os.path.abspath(self.working_dir).strip('/').split('/')
        it = os.path.abspath(path).strip('/').split('/')
        n_same = 0
        while n_same < len(me) and n_same < len(it) and me[n_same] == it[n_same]:
            n_same += 1        
        return os.path.normpath( '/'.join( ['..']*(len(me)-n_same) + it[n_same:] ) )

    def relative_path_as_path(self, path):
        if os.path.isabs(path):
            return path
        
        return os.path.normpath(os.path.join(self.working_dir,path))    
    object_filename = relative_path_as_path
    _object_filename = object_filename
    
    # A single step down the path of darkness, what's the worst that can happen?
    __div__ = relative_path_as_path


def run(args, stdin=None, stdout=subprocess.PIPE, stderr=None):
    return legion.subprocess_Popen(
        args,
        bufsize=1<<24,
        stdin=stdin,        
        stdout=stdout,
        stderr=stderr,
        close_fds=True,
    )

def execute(args, stdin=None, stdout=None):
    p = run(args, stdin=stdin, stdout=stdout)
    assert p.wait() == 0, 'Failed to execute "%s"' % ' '.join(args)


#def peek_and_pipe(f, n_peek):
#    try:
#        f.seek(0)
#        peek = f.read(n_peek)
#        f.seek(0)
#        f.flush() #Necessary in CPython but noy PyPy to actually seek file
#        return peek, f
#    except IOError: #Non-seekable
#        pass
#
#    script = """
#
#import os, sys, select
#n = %d
#f_in = os.fdopen(os.dup(sys.stdin.fileno()),'rb',0)
#f_out = os.fdopen(os.dup(sys.stdout.fileno()),'wb',0)
#try:
#    peek = f_in.read(n)
#    f_out.write(peek + chr(0) * (n-len(peek)))
#    f_out.write(peek)
#    chunks = [ ]
#    done = False
#    while not done or chunks:
#        rlist, wlist, elist = select.select([ f_in ] if len(chunks) < 16 else [ ], [ f_out ] if chunks else [ ], [ ])
#        if wlist:
#            f_out.write(chunks.pop(0))
#        elif rlist:
#            chunk = f_in.read(1<<16)
#            if not chunk:
#                done = True                
#            else:
#                chunks.append( chunk )
#except IOError, error:
#    if error.errno != 32: #Broken pipe
#        raise error
#except KeyboardInterrupt:
#    sys.exit(1)
#
#""" % n_peek
#
#    p = legion.subprocess_Popen(
#         [sys.executable,'-u','-c',script],
#         bufsize=0, #Random crashes if set to, eg, 1<<20
#         stdin=f,       
#         stdout=subprocess.PIPE,
#         close_fds=True,
#         )
#    f.close()
#    
#    peek = p.stdout.read(n_peek).rstrip(chr(0))
#    d = os.dup(p.stdout.fileno())
#    p.stdout.close()
#    f_out = os.fdopen(d,'rb', 1<<20)
#    STREAM_PROCESS[f_out] = p
#    
#    return peek, f_out
#
#def process_buffer(f):
#    return peek_and_pipe(f, 0)[1]  


def is_remote_filename(filename):
    return bool( re.match('\w+:', filename) )

def abspath(*components):
    """ Absolute path of a filename.
        Note: this also ensures the filename does not look like a flag.
    """    
    filename = os.path.join(*components)

    if is_remote_filename(filename):
        return filename
    else:
        return os.path.abspath(filename)

#def open_possibly_remote_file(filename):
#    """ Use lftp to read remote files.
#    """
#
#    if not is_remote_filename(filename):
#        # Doesn't look like a URL
#        return open(filename, 'rb')
#    
#    p = legion.subprocess_Popen(
#         ['lftp', '-c', 'cat', filename],
#         stdout=subprocess.PIPE,
#         close_fds=True,
#         )
#    return p.stdout
    
def open_possibly_compressed_file(filename):
    """ Read a file. It might be compressed. """
    #if isinstance(filename, file):
    #    return filename #Whatever
        
    #peek, f = peek_and_pipe(open_possibly_remote_file(filename), 4)
    
    f = open(filename,'rb')
    peek = f.read(4)
    f.close()
    #f = open(filename,'rb')

    if peek.startswith('\x1f\x8b'):
        #command = 'gunzip'
        return gzip.open(filename, 'rb')
    elif peek.startswith('BZh'):
        #command = 'bunzip2'
        return bz2.BZFile(filename, 'rb')
    else:
        #command = None
        return open(filename, 'rb')
    
    #if command:
    #    p = legion.subprocess_Popen(
    #        [command],
    #        stdin=f,       
    #        stdout=subprocess.PIPE,
    #        close_fds=True,
    #        )
    #    f.close()        
    #    f = p.stdout        
    #
    #return f

def copy_file(source, dest):
    f_in = open_possibly_remote_file(source)
    f_out = open(dest, 'wb')
    while True:
        text = f_in.read(1<<20)
        if not text: break
        f_out.write(text)
    f_out.close()
    f_in.close()


class Pipe_writer(object):
    """ Write to a file via another process. 
    
        Buffering to avoid slowness in pypy.
        
        Drop-in replacement for a file, assuming
        you only use write and close.
    """

    def __init__(self, filename, command):
        self.command = command
        f_out = open(filename,'wb')
        self.process = legion.subprocess_Popen(
            command,
            stdin = subprocess.PIPE,
            stdout = f_out,
    #        bufsize = 1<<24,
            close_fds = True
        )
        f_out.close()
        
        #self.buffer = [ ]
        #self.buf_size = 0

    #def write(self, text):
    #    self.buffer.append( text )
    #    self.buf_size += len(text)
    #    if self.buf_size >= 1<<20: self.flush_buffer()
    
    #def flush_buffer(self):
    #    self.process.stdin.write( ''.join(self.buffer) )
    #    self.buffer = [ ]
    #    self.buf_size = 0
    
    def fileno(self):
        return self.process.stdin.fileno()
    
    def write(self, text):
        self.process.stdin.write( text )
    
    def flush_buffer(self):
        pass
    
    def close(self):
        self.flush_buffer()
        self.process.stdin.close()
        assert self.process.wait() == 0, ' '.join(self.command) + ' failed'


def open_gzip_writer(filename):
    return Pipe_writer(filename, ['gzip'])

def open_bzip2_writer(filename):    
    return Pipe_writer(filename, ['bzip2'])

def open_possibly_compressed_writer(filename):
    if filename[-3:].lower() == '.gz':
        return open_gzip_writer(filename)
    elif filename[-4:].lower() == '.bz2':
        return open_bzip2_writer(filename)
    else:
        return open(filename, 'wb')


def read_solid(filename):
    reads_file = open_possibly_compressed_file(filename)
    
    while True:
        line1 = reads_file.readline()
        while line1.startswith('#'):
            line1 = reads_file.readline()
        if not line1: break
        assert line1.startswith('>'), 'Not a SOLiD CSFASTA file?'
        line2 = reads_file.readline()

        read_name = line1.rstrip('\n')[1:]
        read_seq = line2.rstrip('\n')
        yield read_name, read_seq



def read_illumina(filename):
    reads_file = open_possibly_compressed_file(filename)
    
    while True:
        line1 = reads_file.readline()
        if not line1: break
        line2 = reads_file.readline()
        line3 = reads_file.readline()
        line4 = reads_file.readline()
            
        assert line1.startswith('@'), 'Not an Illumina FASTQ file?'
        assert line3.startswith('+'), 'Not an Illumina FASTQ file?'
            
        read_name = line1.rstrip('\n')[1:]
        assert read_name, 'FASTQ file contains record with no name'
        
        read_seq = line2.rstrip('\n')
        yield read_name, read_seq

def read_illumina_with_quality(filename):
    reads_file = open_possibly_compressed_file(filename)
    
    while True:
        line1 = reads_file.readline()
        if not line1: break
        line2 = reads_file.readline()
        line3 = reads_file.readline()
        line4 = reads_file.readline()
            
        assert line1.startswith('@'), 'Not an Illumina FASTQ file?'
        assert line3.startswith('+'), 'Not an Illumina FASTQ file?'
            
        read_name = line1.rstrip('\n')[1:]
        read_seq = line2.rstrip('\n')
        read_qual = line4.rstrip('\n')
        yield read_name, read_seq, read_qual

def read_fasta(filename):
    reads_file = open_possibly_compressed_file(filename)
    
    line = reads_file.readline()
    while line:
        line = line.rstrip()
        assert line.startswith('>'), 'Not a FASTA file?'
        read_name = line[1:] #.split()[0]
        assert read_name, 'FASTA file contains record with no name'
        
        line = reads_file.readline()
        parts = [ ]
        while line and not line.startswith('>'):
            parts.append(line.rstrip())
            line = reads_file.readline()
        
        yield read_name, ''.join(parts)


def read_gff3_sequence(filename):
    f = open_possibly_compressed_file(filename)
    
    for line in f:
        if line.rstrip() == '##FASTA':
            break
    else:
        raise grace.Error('Tried reading file as a GFF3 but it contains no ##FASTA section')
    
    return read_fasta(f)


def read_genbank_sequence(filename, genbank_callback=None):
    from Bio import SeqIO

    f = open_possibly_compressed_file(filename)
    
    for record in SeqIO.parse(f,'genbank'):
        name = record.id
        if name == '' or name == 'unknown':
            name = record.name
        assert name, 'GENBANK file contains record with no accession or name' 
        
        #Hideous hack: samshrimp makes a copy of each record in a genbank file
        if genbank_callback: genbank_callback(name, record)
        
        yield name, record.seq.tostring()    
    f.close()

def read_empty(filename):
    f = open_possibly_compressed_file(filename)
    f.close()
    return
    yield


def _filter_name(iterator, argument):
    for read in iterator:
        yield (read[0].split()[0],) + read[1:]

def _filter_bar(iterator, argument):
    nth = int(argument)

    for read in iterator:
        parts = read[0].split('|')
        assert len(parts) > nth, 'Not enough parts in: '+name
        yield (parts[nth].strip(),) + read[1:]

def _filter_rename(iterator, argument):
    for read in iterator:
        yield (argument,) + read[1:]

def _filter_select(iterator, argument):
    good = set(argument.split(','))
    for read in iterator:
        if read[0] in good:
            yield read

def _filter_lengthatleast(iterator, argument):
    cutoff = int(argument)
    for read in iterator:
        if len(read[1]) >= cutoff:
            yield read

def _filter_pfilter(iterator, argument):
    for read in iterator:
        variables = {'name':name,'seq':seq}
        if len(read)>2: variables['quality'] = read[2]
        result = eval(argument, variables)
        if not result: continue
        if isinstance(result, str): read = (result,)+read[1:]
        yield read

def _filter_qclip(iterator, argument):
    cutoff = chr(64 + int(argument))
    for name, seq, qual in iterator:
        best_start = 0
        best_end = 0
        start = 0
        for i in xrange(len(qual)):
            if qual[i] < cutoff:
                start = i+1 
            else:
                if i+1-start > best_end-best_start:
                    best_start = start
                    best_end = i+1
        yield name, seq[best_start:best_end], qual[best_start:best_end] 

def _filter_first(iterator, argument):
    n = int(argument)
    for item in iterator:
        if n <= 0: break
        n -= 1
        yield item

FILTERS = {
    'name': _filter_name,
    'bar:': _filter_bar,
    'rename:': _filter_rename,
    'select:': _filter_select,
    'lengthatleast:': _filter_lengthatleast,
    'pfilter:': _filter_pfilter,
    'qclip:': _filter_qclip,
    'first:': _filter_first,
}

def filter_no_qualities(iterator):
    for read in iterator:
        yield read[:2]    

def read_sequences(filename, qualities=False, genbank_callback=None):
    """ Read fasta or illumina sequences, possibly compressed 
    
        Post reading filters can be applied.
    """
    
    parts = filename.split('~~')

    f = open_possibly_compressed_file(parts[0])
    peek = f.read(8)
    f.close()
    
    have_qualities = False
    
    if not peek:
        result = read_empty(parts[0])
    elif peek.startswith('>'):
        result = read_fasta(parts[0])
    elif peek.startswith('LOCUS'):
        result = read_genbank_sequence(parts[0], genbank_callback)
    elif peek.startswith('@'):
        have_qualities = True
        result = read_illumina_with_quality(parts[0])    
    elif peek.startswith('##gff'):
        result = read_gff3_sequence(parts[0])
    elif peek.startswith('.sff'):
        f.close()
        grace.require_sff2fastq()
        have_qualities = True
        process = run(['sff2fastq', parts[0]])
        result = read_illumina_with_quality(process.stdout)
    else:
        raise grace.Error('Unrecognized file format for '+filename)
    
    for part in parts[1:]:
        for prefix in FILTERS:
            if part.lower().startswith(prefix):
                result = FILTERS[prefix](result, part[len(prefix):])
                break
        else:
            raise grace.Error('Unrecognized filter: '+part)

    if have_qualities and not qualities:
        result = filter_no_qualities(result)
    
    return result   


def is_sequence_file(filename):
    if not os.path.isfile(filename): 
        return False
    try:
        read_sequences(filename)
    except grace.Error:
        return False
    return True


def guess_quality_offset(filename):
    grace.status('Guessing quality offset')
    try:
        min_value = chr(255)
        #max_value = chr(0)
        for i, item in enumerate(read_sequences(filename, qualities=True)):
            if len(item) == 2: return 33 #Not fastq
            
            min_value = min(min_value, min(item[2]))
            #max_value = max(max_value, max(item[2]))
            
            if i >= 100000: break
        
        low = ord(min_value)
        #high = ord(max_value)
        #print 'Quality chars in range %d-%d in %s' % (low,high,filename)
        if low < 59: return 33 #Sanger and Illumina 1.8+
        return 64 #Illumina pre 1.8
    finally:
        grace.status('')


def check_name_uniqueness(filenames):
    """ Check first few read names are unique """
    names = set()
    for filename in filenames:
        for i, (name, seq) in io.read_sequences(filename):
            assert name not in names, 'Duplicate sequence name: '+name
            names.add(name)
            if i >= 1000: break


def write_fasta(f, name, sequence, qual=None):
    print >> f, '>' + name
    for i in xrange(0,len(sequence),70):
        print >> f, sequence[i:i+70] 

def write_fasta_single_line(f, name, sequence, qual=None):
    print >> f, '>' + name
    print >> f, sequence 

def write_fastq(f, name, sequence, qual):
    print >> f, '@' + name
    print >> f, sequence
    #for i in xrange(0,len(sequence),70):
    #    print >> f, sequence[i:i+70]
    print >> f, '+' 
    print >> f, qual
    #for i in xrange(0,len(qual),70):
    #    print >> f, qual[i:i+70]



def decode_evidence(desc):
    result = [ ]
    for item in desc.split():
        match = re.match('(.*)x([0-9]+)$', item)
        result.append( (match.group(1).replace('"',''), int(match.group(2))) )
        # Old evidence format example: "A"x42
        # New evidence format example: Ax42
        # Handle both
    return result

def read_evidence_file(filename):
    f = open(filename,'rb')
    f.readline()
    for line in f:
        parts = line.rstrip('\n').split('\t')
        yield int(parts[0]), parts[1], parts[2], parts[3]



class Table_reader(object):
    def __init__(self, filename):
        self.f = open(filename, 'rb')
        line = self.f.readline()
        while line and line.startswith('#') or not line.strip():
            line = self.f.readline()    
        
        assert line, 'Table has not even a heading'
        
        if '\t' in line:
            self.parse = lambda line: line.rstrip('\n').split('\t')
        elif ',' in line:
            self.parse = lambda line: csv.reader([line]).next()
        else:
            assert False, 'Strange table'
        
        self.headings = self.parse(line)
    
    def __iter__(self):
        return self
        
    def next(self):
        while True:
            line = self.f.readline()
            if not line: raise StopIteration()        
            if line.startswith('#') or not line.strip(): continue
            
            values = self.parse(line)
            if not values: continue
            
            assert len(values) == len(self.headings)
            return collections.OrderedDict(zip(self.headings,values))

read_table = Table_reader


def write_csv(filename, iterable):
    """ Write a sequence of OrderedDicts of strings as a CSV """
    f = open(filename, 'wb')
    writer = csv.writer(f)
    keys = None
    
    for record in iterable:
        if keys is None:
           keys = record.keys()
           writer.writerow(keys)
        assert record.keys() == keys
        writer.writerow(record.values())
    f.close()







