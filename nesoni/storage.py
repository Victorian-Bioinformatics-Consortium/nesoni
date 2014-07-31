
import mmap, struct, cPickle, threading, platform
import numpy as np

LOCAL = threading.local()
LOCAL.storage = None

pypy = platform.python_implementation() == 'PyPy'
if pypy:
    import cffi
    ffi = cffi.FFI()
    ffi.cdef('typedef size_t off_t;')
    ffi.cdef('void *mmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset);')
    ffi.cdef('int munmap(void *addr, size_t length);')
    C = ffi.dlopen(None)

    def get_mmap(fileno, size):
        addr = C.mmap(ffi.NULL, size, mmap.PROT_READ|mmap.PROT_WRITE, mmap.MAP_SHARED, fileno, 0)
        assert addr != -1, 'mmap failed'
        def close():
            C.munmap(addr, size)
        
        return ffi.buffer(addr, size), close

else:
    def get_mmap(fileno, size):
        result = mmap.mmap(fileno, size, access=mmap.ACCESS_WRITE)
        return result, result.close


class Disk_array(np.ndarray):
    def __new__(cls, offset, shape, dtype):
        result = np.ndarray.__new__(
            Disk_array,
            shape=shape,
            dtype=dtype,
            offset=offset,
            buffer=LOCAL.storage._mmap)
        result._disk_array_storage = LOCAL.storage
        result._disk_array_new = (offset, shape, dtype)
        return result
    
    def __reduce__(self):
        return (Disk_array, self._disk_array_new)

    

class _Store(object): 
    """ Object to pickle """
    pass
    
class _Value(object): 
    """ Default value type """
    pass


class Storage(object):
    """
    By means of pickle-fu,
    store data on disk, 
    including arbitrary python objects
    and memory mapped numpy arrays.
    
    The idea is to work with a fairly small amount of python objects
    which are fully unpickled on load,
    and some very large memory mapped arrays,
    the memory for which is only loaded on demand.
    
    Access data via .data

    Create memory-mapped arrays with .empty, .zeros, .ones,
    store them in .data. Call .save()
    
    Garbage collection is presently not implemented.
    It would be relatively easy to implement as part of ._load()
    (ie a GC pass would require saving and loading).
    
    As the file grows it is mmapped multiple times.
    Old mmaps are retained as old arrays point into them.
    The file grows by (default) 1.5x each time, 
    so not too many mappings are required.
    This seems to work.

    Call .close() to unmap mmaps. 
    Things will explode 
    if you attempt to access memory mapped arrays after this.        
    """
    
    
    def __init__(self, filename, clear=False, preallocate=0, expand=1.5):
        self._closers = [ ] # Functions to close all mmaps
        self._close_size = None # Truncate any unused space on close
        
        self._filename = filename
        self._file = open(filename, 'ab+')
        self._file.seek(0, 2)
        if clear or self._file.tell() < 24:
            self._file.truncate(0)
            self._file.write('nesoni00' + chr(0)*16)
        self._file.seek(0, 2)
        self._size = self._file.tell()
        self._expand = expand
        
        self._mmap = None
        self._require(24 + preallocate)
        
        self._load()
    
    def close(self):
        del self._obj
        while self._closers:
            self._closers.pop()()
        if self._close_size is not None:
            self._file.truncate(self._close_size)
        self._file.close()
    
    @property
    def data(self):
        return self._obj.value
    
    @data.setter
    def data(self, value):
        self._obj.value = value

            
    def save(self):
        dump = self._in_context(cPickle.dumps, self._obj, 2)
        offset = self._obj.alloc
        size = len(dump)
        self._require(offset+size)
        self._mmap[8:24] = struct.pack('<QQ', offset,size)
        self._mmap[offset:offset+size] = dump
        self._close_size = offset+size
        
    def _require(self, size):
        if self._size < size:
            size = max(int(self._size*self._expand), size)
            self._file.truncate(size)
            self._size = size
            self._mmap = None

        if self._mmap is None:
            self._mmap, closer = get_mmap(self._file.fileno(), self._size)
            self._closers.append(closer)
    
    def _load(self):
        self._require(24)
        offset, size = struct.unpack('<QQ', self._mmap[8:24])
        if size == 0:
            self._obj = _Store()
            self._obj.alloc = 24
            self._obj.value = _Value()
        else:
            self._file.seek(offset)
            self._obj = self._in_context(cPickle.loads,self._mmap[offset:offset+size])
    
    def _in_context(self, func, *args, **kwargs):
        old = LOCAL.storage
        try:
            LOCAL.storage = self
            return func(*args, **kwargs)
        finally:
            LOCAL.storage = old

    def _allocate(self, size):
        result = self._obj.alloc
        self._obj.alloc += size
        self._require(self._obj.alloc)
        return result

    def empty(self, shape, dtype):
        if isinstance(shape, int):
            shape = (shape,)
        
        size = np.dtype(dtype).itemsize
        for item in shape:
            size *= item
                
        offset = self._allocate(size)
        return self._in_context(Disk_array, offset, shape, dtype)

    def zeros(self, shape, dtype):
        result = self.empty(shape,dtype)
        result[()] = 0
        return result

    def ones(self, shape, dtype):
        result = self.empty(shape,dtype)
        result[()] = 1
        return result

    def array(self, data):
        # Slight inefficiency
        data = np.asarray(data)
        result = self.empty(data.shape, data.dtype)
        result[()] = data
        return result




