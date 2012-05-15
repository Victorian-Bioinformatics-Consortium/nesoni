"""

    Key value store, with fast writing, semigroup update semantics

"""

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"
    ctypedef void** const_void_ptr_ptr "const void**"

cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object obj, const_void_ptr_ptr buffer, Py_ssize_t *buffer_len)
    object PyString_FromStringAndSize(char *v, int len)
    
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)
    void  fprintf(FILE* f, char* s, char* s)
    
    ctypedef unsigned long size_t
    FILE *fopen(char *path, char *mode)
    int fclose(FILE *fp)
    size_t fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
    int fputc(int c, FILE *stream)
    long ftell(FILE *stream)
    
    int memcmp(void *s1, void *s2, size_t n)

ctypedef FILE *FILE_ptr
ctypedef char *char_ptr
    
import sys, os, mmap, heapq, itertools

from treemaker import worker_pool

cdef inline int compare_strings(char *str1, long len1, char *str2, long len2):
    cdef long min_length
    if len1 < len2:
        min_length = len1
    else:
        min_length = len2
        
    cdef int result = memcmp(
        str1,
        str2,
        min_length
    )
        
    if not result:
        if len1 < len2:
            result = -1
        elif len1 > len2:
            result = 1
        
    return result


cdef inline void write_int64(FILE_ptr f, long value):
    fwrite( &value, 8, 1, f )

cdef inline str string_from_int64(long value):
    return PyString_FromStringAndSize(<char*>&value, 8)

cdef inline long int64_from_string(string): #TODO: restore to str string when Cython is fixed
    return (<long*><char*>string)[0]

cdef inline void write_int(FILE_ptr f, long value):
    """ Write an unsigned integer, concisely """

    while value >= 128:
        fputc( (value&127)|128, f )
        value >>= 7
    fputc( value, f )

cdef inline long read_int(char **data):
    """ Read a concise unsigned integer """

    cdef long result = 0
    cdef unsigned char x
    cdef int shift = 0
    
    while True:
        x = data[0][0]
        data[0] += 1        
        result += (<long>(x&127)) << shift
        if (x&128) == 0: break
        shift += 7
    return result


cpdef remove_tree(path):
    if not os.path.exists(path): return
    
    for filename in os.listdir(path):
        assert filename.isdigit()
        os.unlink(os.path.join(path, filename))
    
    os.rmdir(path)


cpdef make_tree(path, sequence):
    cdef long i, j, ones, left_loc, right_loc
    
    #TODO: restore when Cython is fixed
    #cdef str key, value
    
    cdef int n_levels = 0
    cdef FILE_ptr files[64]    
    cdef long last[64]

    remove_tree(path)
    os.mkdir(path)
    
    i = 0
    for key, value in sequence:
        # The run of least significant ones determines the depth in the balanced tree
        ones = 0
        j = i
        while i & (1<<ones):
            ones += 1
        
        while n_levels <= ones:
            filename = os.path.join(path, str(n_levels))
            files[n_levels] = fopen(filename, 'wb')
            last[n_levels] = 0
            n_levels += 1            

        # last[ones] = levels[ones].tell()
        
        if ones:
            left_loc = last[ones-1]
            for j in xrange(ones):            
                last[j] = ftell(files[j])
            right_loc = last[ones-1]
            
            write_int(files[ones], left_loc)
            
            # Right should not be far from left,
            # write as a relative position to safe space
            write_int(files[ones], right_loc - left_loc)
            
        write_int(files[ones], len(key))        
        fwrite( <char*>key, len(key), 1, files[ones] )

        write_int(files[ones], len(value))
        fwrite( <char*>value, len(value), 1, files[ones] )
        
        i += 1

    if n_levels: 
        # The last little bits of the sequence might not be a part of the binary tree yet
        # so need to store their location
        # See also: Cursor.right()
    
        write_int64(files[0], 0)
        for i in xrange(0,n_levels-1):
            write_int64(files[i+1], last[i])
            
    for i in xrange(n_levels): 
        fclose(files[i])


cdef class Cursor

cdef class Tree:
    cdef public object path, levels
    cdef Py_ssize_t sizes[64]
    cdef char_ptr maps[64]

    def __init__(self, path):
        self.path = path
        self.levels = [ ]
        while True:
             filename = os.path.join(path, str(len(self.levels)))
             if not os.path.exists(filename): break
             f = open(filename,'rb')
             self.levels.append(mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ))
             f.close()
      
        for i in xrange(len(self.levels)):
            #self.sizes[i] = len(self.levels[i])
            PyObject_AsReadBuffer(self.levels[i], <const_void_ptr_ptr>&(self.maps[i]),
                <Py_ssize_t*>&(self.sizes[i]))

    def __dealloc__(self):
        for level in self.levels:
            level.close()

    cpdef size(self):
        return sum([ len(item) for item in self.levels ])

    cpdef Cursor cursor(self):
        if not self.levels: return None
        return Cursor(self, len(self.levels)-1, 0, None)
        

cdef class Cursor:
    cdef public Tree tree
    cdef public Cursor successor
    cdef public long level

    cdef public long left_loc, right_loc, key_length
    cdef char *data_ptr

    def __cinit__(Cursor self, Tree tree, long level, long position, Cursor successor):
        self.tree = tree
        self.level = level
        #self.position = position
        self.successor = successor
        
        self.data_ptr = tree.maps[level] + position
        if level:
            self.left_loc = read_int(&self.data_ptr)
            self.right_loc = self.left_loc + read_int(&self.data_ptr)
        self.key_length = read_int(&self.data_ptr)

    def __cmp__(Cursor self, Cursor other):
        return compare_strings(self.data_ptr,self.key_length,other.data_ptr,other.key_length)
    cpdef int compare(Cursor self, Cursor other):
        return compare_strings(self.data_ptr,self.key_length,other.data_ptr,other.key_length)

    # This is weirdly broken.
    #cpdef int compare_to_string(Cursor self, char *string, long length):        
    #    return compare_strings(self.data_ptr,self.key_length,string,length)

    cpdef int compare_to_str(Cursor self, string): #TODO: restore to str string when Cython is fixed        
        cdef char *str_ptr = string
        return compare_strings(self.data_ptr,self.key_length,str_ptr, len(string))
        
    cpdef str key(Cursor self):
        return PyString_FromStringAndSize(self.data_ptr, self.key_length)
        
    cpdef str value(Cursor self):
        cdef char *value_ptr = self.data_ptr + self.key_length
        cdef long value_length = read_int(&value_ptr)
    
        return PyString_FromStringAndSize(value_ptr, value_length)

    cpdef Cursor left(Cursor self):
        if self.level == 0: return None
        
        return Cursor(self.tree, self.level-1, self.left_loc, self)

    cpdef Cursor right(Cursor self):
        if self.level == 0: return None

        cdef long new_level = self.level-1
        cdef long new_position = self.right_loc
        
        while new_position == self.tree.sizes[new_level] - 8:
            if new_level == 0: return None
        
            new_position = (<long*>(self.tree.maps[new_level]+new_position))[0]            
            new_level -= 1
        
        return Cursor(self.tree, new_level, new_position, 
                      self.successor)
    
    cpdef Cursor next(Cursor self):
        cdef Cursor right = self.right()
        if right is None: return self.successor
        return left_most(right)


cpdef Cursor left_most(Cursor cursor):
    if not cursor: return None
    while cursor.level > 0:
        cursor = cursor.left()
    return cursor


cdef class Cursor_merge(object):
    cdef Cursor cursor1, cursor2
    cdef object reducer

    def __cinit__(self, cursor1, cursor2, reducer):
        self.cursor1 = cursor1
        self.cursor2 = cursor2
        self.reducer = reducer

    def __iter__(self):
        return self

    def __next__(self):        
        if self.cursor1 is None:
            if self.cursor2 is None:
                raise StopIteration()
                
            result = (self.cursor2.key(), self.cursor2.value())
            self.cursor2 = self.cursor2.next()
            return result

        if self.cursor2 is None:
            result = (self.cursor1.key(), self.cursor1.value())
            self.cursor1 = self.cursor1.next()
            return result
        
        cdef int comparison = self.cursor1.compare(self.cursor2)
        
        if comparison < 0:
            result = (self.cursor1.key(), self.cursor1.value())
            self.cursor1 = self.cursor1.next()
            return result

        elif comparison > 0:
            result = (self.cursor2.key(), self.cursor2.value())
            self.cursor2 = self.cursor2.next()
            return result

        else:
            result = (self.cursor1.key(), self.reducer(self.cursor1.value(),self.cursor2.value()))
            self.cursor1 = self.cursor1.next()
            self.cursor2 = self.cursor2.next()
            return result


cdef class Cursor_merge_all:
    cdef list cursor_heap # heap of (cursor, priority)
    cdef object reducer

    def __cinit__(self, cursors, reducer):
        # Cursors should be passed as a list from newest to oldest

        self.reducer = reducer    
        self.cursor_heap = [ ]
        
        for i, cursor in enumerate(cursors):
            heapq.heappush(self.cursor_heap, (cursor, i))
    
    def __iter__(self):
        return self

    def __next__(self):
        if not self.cursor_heap:
            raise StopIteration()
        
        heap_items = [ heapq.heappop(self.cursor_heap) ]
        while self.cursor_heap and heap_items[0][0].compare(self.cursor_heap[0][0]) == 0:
            heap_items.append(heapq.heappop(self.cursor_heap))
        
        for cursor, oldness in heap_items:
            next_cursor = cursor.next()
            if next_cursor is not None:
                heapq.heappush(self.cursor_heap, (next_cursor, oldness))
                
        result = heap_items[0][0].value()
        for cursor, oldness in heap_items[1:]:
            result = self.reducer( cursor.value(), result )
        return (heap_items[0][0].key(), result)


cdef class Skip_discardables:
    cdef object cursor_merger, is_discardable

    def __cinit__(self, cursor_merger, is_discardable):
        self.cursor_merger = cursor_merger
        self.is_discardable = is_discardable
    
    def __iter__(self):
        return self
    
    def __next__(self):
        result = self.cursor_merger.next() #May throw Stop_iteration, that's fine
        while self.is_discardable(result[1]):
            result = self.cursor_merger.next()
        return result




cpdef _default_is_final(item): return True
cpdef _default_is_discardable(item): return False
cpdef _default_reducer(item_older, item_newer): return item_newer

cdef class Store:
    """ Key-value store with per-key reductions
    
        Default behaviour is to over-write existing keys,
        and to return None if a key is not present
    
        Order is important!
        reducer takes arguments like this: reducer(older item, newer item)
    
    """

    cdef public object path, trees, cursors, next_tree_id, read_only
    cdef public object memory, memory_size
    cdef public object lock, busy_paths,
    # merger_processes
    #cdef public object manager, worker_pool, unbusy_queue
    cdef public object pool

    cdef public object default_value, is_discardable, is_final, reducer

    def __init__(self, 
                 path, 
                 default_value=None, 
                 is_final=_default_is_final, 
                 is_discardable=_default_is_discardable, 
                 reducer=_default_reducer, 
                 memory_size=100000, #Formerly: 500000
                 read_only=False
                 ):
        self.read_only = read_only
        
        self.path = path        
        if not os.path.exists(path):
            os.mkdir(path)
        
        self.default_value = default_value
        self.is_final = is_final
        self.is_discardable = is_discardable
        self.reducer = reducer

        self.memory_size = memory_size
        self.memory = { }
        
        #self.manager = multiprocessing.Manager()
        #self.lock = self.manager.Lock()
        #self.merger_processes = [ ]
        self.busy_paths = set()
        
        #self.unbusy_queue = self.manager.Queue()
        #self.unbusy_inpipe, self.unbusy_outpipe = multiprocessing.Pipe(False)
        #self.worker_pool = multiprocessing.Pool()
        
        self.lock = worker_pool.Lock()
        self.pool = worker_pool.Worker_pool(do_mergedown, (self.lock,self.reducer))
        
        self.trees = [ ]
        
        self._clean_house()
        self._remap()

    cpdef close(self):
        self.flush()
        self.finish_merging()
    
    cpdef destroy(self):
        """ Delete entire store directory """
        self.clear()
        self.close()
        os.rmdir(self.path)
    
    cpdef _clean_house(self):
        """ Called by init. """
        if self.read_only: return
        
        for filename in os.listdir(self.path):
            full_name = os.path.join(self.path, filename)
            
            if filename == 'newtree' or \
               filename.endswith('-build'):
                remove_tree(full_name)
            
            elif filename.endswith('-ready'):
               tree1, tree2, ready_string = filename.split('-')
               tree1 = os.path.join(self.path, tree1)      
               tree2 = os.path.join(self.path, tree2)
               remove_tree(tree1)
               remove_tree(tree2)
               os.rename(full_name, tree1)      
                
    cpdef _remap(self):
        self.lock.acquire()
    
        list_numbers = [ int(filename) 
                         for filename in os.listdir(self.path)
                         if filename.isdigit() ]
        list_numbers.sort()
        
        if not list_numbers:
            self.next_tree_id = 0
        else:
            self.next_tree_id = list_numbers[-1] + 1
    
        self.trees = [ ]
        for i in list_numbers:
            filename = os.path.join(self.path, str(i))            
            if not os.path.exists(filename): break
            self.trees.append( Tree(filename) )
                
        self.lock.release()
        
        self.cursors = [ ]
        for tree in self.trees[::-1]: #Note reverse order!
            cursor = tree.cursor()
            if cursor is not None: self.cursors.append(cursor)
        
    
    cpdef clear(self):
        self.finish_merging()
    
        for tree in self.trees:
            remove_tree(tree.path)
        self.memory = { }
                
        self._remap()

    
    cpdef finish_merging(self, wait=True):
        #still_waiting = [ ]
        #any = False
        #for item in self.merger_processes:
        #    if not wait and item[0].is_alive():
        #        still_waiting.append(item)
        #        continue
        #    
        #    #if item[0].is_alive(): print 'wait!'
        #    
        #    item[0].join()
        #    self.busy_paths.remove(item[1])
        #    self.busy_paths.remove(item[2])
        #    #print 'done', item[1], item[2]
        #    any = True
        #
        #self.merger_processes = still_waiting
        #
        
        #any = False
        #while self.busy_paths and (wait or not self.unbusy_queue.empty()):
        #    self.busy_paths.remove( self.unbusy_queue.get() )
        #    any = True

        any = False
        while self.busy_paths:
            try:
                result = self.pool.get(wait)
            except worker_pool.Empty:
                break
            
            for path in result:
                self.busy_paths.remove(path)
            any = True
        
        if any:
            self._remap()
            
        return any
        
    
    cpdef start_merging(self, all=False):
        self.finish_merging(False)
        
        pos = len(self.trees)-2
        while pos >= 0:
            if self.trees[pos].path not in self.busy_paths and \
               self.trees[pos+1].path not in self.busy_paths and \
               (all or self.trees[pos].size() <= self.trees[pos+1].size() * 2):
                path1 = self.trees[pos].path
                path2 = self.trees[pos+1].path
                self.busy_paths.add(path1)
                self.busy_paths.add(path2)
                #print 'merge', path1, path2
                
                if pos == 0:
                    is_discardable_if_needed = self.is_discardable
                else:
                    is_discardable_if_needed = None
                
                #process = multiprocessing.Process(
                #    target=do_mergedown,
                #    args=(self.path, path1, path2,
                #          self.reducer, is_discardable_if_needed,
                #          self.lock, self.unbusy_queue))
                #process.start()
                #self.merger_processes.append((process,path1,path2))
                
                #self.worker_pool.apply_async(
                #    do_mergedown,
                #    (self.lock, self.reducer, self.path, path1, path2,
                #     self.reducer, is_discardable_if_needed))

                self.pool.put((
                     self.path, path1, path2,
                     is_discardable_if_needed))
                
                pos -= 2
            else:
                pos -= 1


    def optimize(self):
        self.flush()
    
        #Not optimal? Need a way to wait for any one process to finish.
        while len(self.trees) > 1:
            self.start_merging(True)
            self.finish_merging(True)        

    
    cpdef flush(self):
        if not self.memory: return
        
        tempname = os.path.join(self.path, 'newtree')
        filename = os.path.join(self.path, str(self.next_tree_id))

        # This seems to be about the fastest way to do this, surprisingly
        items = [ (item, self.memory[item]) for item in sorted(self.memory.iterkeys()) ]
        self.memory = {}
        
        try:
            make_tree(tempname, items)
        except:
            remove_tree(tempname)
            raise
        
        os.rename(tempname, filename)
        self._remap()        
        self.start_merging()
        
    def _put(self, key, value):
        if key in self.memory:
            self.memory[key] = self.reducer(self.memory[key], value)
        else:
            self.memory[key] = value
            if len(self.memory) >= self.memory_size: 
                self.flush()       

    def __getitem__(self, key):
        return self._get(key)
    
    def _get(self, key):
        cdef Cursor cursor

        if self.busy_paths:
            if self.finish_merging(False):
                self.start_merging()
        
        if key in self.memory:
            result = self.memory[key]
            if self.is_final(result): return result
        else:
            result = None
            
        #cdef char *key_ptr = key
        #cdef long key_len = len(key)
        cdef int comparison
        
        for cursor in self.cursors:
            while cursor:
                comparison = cursor.compare_to_str(key)
                if comparison == 0: 
                    this_value = cursor.value()
                    if result is None:
                        result = this_value
                    else:
                        result = self.reducer(this_value, result)
                    if self.is_final(result):
                        return result
                    break
                    
                if comparison > 0:
                    cursor = cursor.left()
                else:
                    cursor = cursor.right()
        
        if result is None:
            result = self.default_value
        
        return result
        
    def __contains__(self, key):
        return self._get(key) != self.default_value

    def iter_from(self, key):
        self.flush() # self.memory is a dict, not ordered, need to flush to disk to sort it
    
        #cdef char *key_ptr = key
        #cdef long key_len = len(key)
        
        cursors = [ ]
        
        for cursor in self.cursors:
            while True:
                comparison = cursor.compare_to_str(key)
                if comparison == 0:
                    break
                
                if comparison > 0:
                    new_cursor = cursor.left()
                else:
                    new_cursor = cursor.right()
                    
                if new_cursor is None: break
                
                cursor = new_cursor
            
            if cursor.compare_to_str(key) < 0:
                cursor = cursor.next()
            if cursor is not None:
                cursors.append(cursor)
        
        return Skip_discardables(
            Cursor_merge_all(cursors, self.reducer), 
            self.is_discardable)
            
            

def do_mergedown(lock, reducer, path, path1, path2, is_discardable):
    name1 = os.path.split(path1)[1]
    name2 = os.path.split(path2)[1]
    
    new_filename_base = os.path.join(path, name1+'-'+name2)
    new_filename_build = new_filename_base + '-build'
    new_filename_ready = new_filename_base + '-ready'
    
    lock.acquire()
    tree1 = Tree(path1)
    tree2 = Tree(path2)
    lock.release()
    
    merger = Cursor_merge(
                left_most(tree1.cursor()),
                left_most(tree2.cursor()),
                reducer)
        
    if is_discardable is not None:
        merger = Skip_discardables(merger, is_discardable)
    
    try:
        make_tree(new_filename_build, merger)
    except:
        remove_tree(new_filename_build)
        raise
        
    lock.acquire()
    os.rename(new_filename_build, new_filename_ready)
    remove_tree(path1)
    remove_tree(path2)
    os.rename(new_filename_ready, path1)
    lock.release()
    
    #unbusy_queue.put(path1)
    #unbusy_queue.put(path2)
    return path1, path2    


_counting_store_default = string_from_int64(0)

cpdef _counting_store_is_final(item): 
    return False

cpdef _counting_store_is_discardable(item): 
    return item == _counting_store_default

cpdef _counting_store_reduce(item_older, item_newer):
    return string_from_int64(int64_from_string(item_older) + int64_from_string(item_newer))

cpdef _counting_store_iter_filter(pair):
    return (pair[0], int64_from_string(pair[1]))

cdef class Counting_store(Store):
    """ Store for counting strings with """
    
    def __init__(self, path, **kwargs):
        Store.__init__(
            self,
            path,
            default_value = _counting_store_default,
            is_final = _counting_store_is_final,
            is_discardable = _counting_store_is_discardable,
            reducer = _counting_store_reduce,
            **kwargs
        )

    def __getitem__(self, key):
        return int64_from_string(self._get(key))

    def add(self, key, amount=1):
        self._put(key, string_from_int64(amount))

    def iter_from(self, key):
        return itertools.imap(
            _counting_store_iter_filter,
            Store.iter_from(self,key))



_count_and_shadow_store_default = string_from_int64(0) * 2

cpdef _count_and_shadow_store_is_final(item): 
    return False

# Discard shadow if count == 0
cpdef _count_and_shadow_store_is_discardable(item): 
    return item[:8] == _counting_store_default

cpdef _count_and_shadow_store_reduce(item_older, item_newer):
    older_count = int64_from_string(item_older[0:8])
    older_shadow = int64_from_string(item_older[8:16])
    newer_count = int64_from_string(item_newer[0:8])
    newer_shadow = int64_from_string(item_newer[8:16])
    return string_from_int64(older_count+newer_count) + string_from_int64(max(older_shadow,newer_shadow))

cpdef _count_and_shadow_decode(item):
    count = int64_from_string(item[0:8])
    shadow = int64_from_string(item[8:16])
    return (count, shadow)

cpdef _count_and_shadow_store_iter_filter(pair):
    return (pair[0], _count_and_shadow_decode(pair[1]))

cdef class Count_and_shadow_store(Store):
    """ Store for counting strings with """
    
    def __init__(self, path, **kwargs):
        Store.__init__(
            self,
            path,
            default_value = _count_and_shadow_store_default,
            is_final = _count_and_shadow_store_is_final,
            is_discardable = _count_and_shadow_store_is_discardable,
            reducer = _count_and_shadow_store_reduce,
            **kwargs
        )

    def __getitem__(self, key):
        return _count_and_shadow_decode(self._get(key))

    def add(self, key, amount, shadow_amount):
        self._put(key, string_from_int64(amount) + string_from_int64(shadow_amount))

    def iter_from(self, key):
        return itertools.imap(
            _count_and_shadow_store_iter_filter,
            Store.iter_from(self,key))




cpdef _general_store_is_discardable(item): 
    return item == ''
    
cpdef _general_store_is_final(item): 
    return item == '' or item[0] == 'S'

cpdef _general_store_reduce(item_older, item_newer):
    #Delete
    if item_newer == '': 
        return item_newer
    
    #Set
    if item_newer[0] == 'S':
        return item_newer

    #Append
    if item_newer[0] == 'A':
        #to deleted
        if item_older == '':
            return 'S' + item_newer[1:]
        
        #to set or append
        return item_older + item_newer[1:]

cpdef _general_store_iter_filter(pair):
    return (pair[0], pair[1][1:])        

cpdef _general_store_iter_keys_filter(pair):
    return pair[0]        

cpdef _identity(x): return x

cdef class General_store(Store):
    """ Store that allows setting, deleting and appending 
    
        You can specify functions to encode and decode values.
        If you will be using append_to, the decoder
        should cope with encoded values being concatenated together.
    """
    
    cdef public object value_encoder, value_decoder

    def __init__(self, path,
                 value_encoder=_identity,
                 value_decoder=_identity, 
                 **kwargs):
        
        self.value_encoder = value_encoder
        self.value_decoder = value_decoder
        
        Store.__init__(
            self,
            path,
            default_value = '',
            is_discardable = _general_store_is_discardable,
            is_final = _general_store_is_final,
            reducer = _general_store_reduce,
            **kwargs
        )

    def append_to(self, key, value):
        self._put(key, 'A' + self.value_encoder(value))

    def __setitem__(self, key, value):
        self._put(key, 'S' + self.value_encoder(value))

    def __delitem__(self, key):
        self._put(key, '')

    def __getitem__(self, key):
        result = self._get(key)        
        if result == '': 
            raise KeyError(key)

        return self.value_decoder( result[1:] )

    def get(self, key, default):
        result = self._get(key)        
        if result == '': 
            return default

        return self.value_decoder( result[1:] )

    def _iter_filter(self, pair):
        return (pair[0], self.value_decoder( pair[1][1:] ))
                
    def iter_from(self, key):
        return itertools.imap(
            self._iter_filter,
            Store.iter_from(self,key))

    def iter_keys_from(self, key):
        return itertools.imap(
            _general_store_iter_keys_filter,
            Store.iter_from(self,key))



cpdef _int_list_store_encode(values):
    return ''.join([ string_from_int64(item) for item in values ])
    
cpdef _int_list_store_decode(str string):
    return [ int64_from_string(string[i:i+8]) for i in xrange(0,len(string),8) ]

cdef class Int_list_store(General_store):
    def __init__(self, path, **kwargs):
        General_store.__init__(self, path,
            value_encoder=_int_list_store_encode,
            value_decoder=_int_list_store_decode,
            **kwargs)



