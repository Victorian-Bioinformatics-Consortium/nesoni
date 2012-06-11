"""

    Key value store, with fast writing, semigroup update semantics

"""
    
import sys, os, mmap, heapq, itertools, struct


def write_int64(f, value):
    f.write(struct.pack('<q', value))

def string_from_int64(value):
    return struct.pack('<q', value)

def int64_from_string(string): #TODO: restore to str string when Cython is fixed
    return struct.unpack('<q', string)[0]

def write_int(f, value):
    """ Write an unsigned integer, concisely """

    while value >= 128:
        f.write( chr( (value&127)|128 ) )
        value >>= 7
    f.write(chr(value))

def read_int(data, offset):
    """ Read a concise unsigned integer """

    result = 0
    shift = 0
    while True:
        x = ord(data[offset])
        offset += 1
        result += (x&127) << shift
        if (x&128) == 0: break
        shift += 7
    return result, offset


def remove_tree(path):
    if not os.path.exists(path): return
    
    for filename in os.listdir(path):
        assert filename.isdigit()
        os.unlink(os.path.join(path, filename))
    
    os.rmdir(path)


def make_tree(path, sequence):
    n_levels = 0
    files = [ ]
    last = [ ]

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
            files.append( open(filename, 'wb') )
            last.append( 0 )
            n_levels += 1            

        if ones:
            left_loc = last[ones-1]
            for j in xrange(ones):            
                last[j] = files[j].tell()
            right_loc = last[ones-1]
            
            write_int(files[ones], left_loc)
            
            # Right should not be far from left,
            # write as a relative position to save space
            write_int(files[ones], right_loc - left_loc)
            
        write_int(files[ones], len(key))
        files[ones].write(key)        

        write_int(files[ones], len(value))
        files[ones].write(value)
        
        i += 1

    if n_levels: 
        # The last little bits of the sequence might not be a part of the binary tree yet
        # so need to store their location
        # See also: Cursor.right()
    
        write_int64(files[0], 0)
        for i in xrange(0,n_levels-1):
            write_int64(files[i+1], last[i])
            
    for i in xrange(n_levels): 
        files[i].close()


class Tree:
    def __init__(self, path):
        self.path = path
        self.levels = [ ]
        while True:
             filename = os.path.join(path, str(len(self.levels)))
             if not os.path.exists(filename): break
             f = open(filename,'rb')
             self.levels.append(mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ))
             f.close()

        self.sizes = [ len(item) for item in self.levels ]      

    def __dealloc__(self):
        for level in self.levels:
            level.close()

    def size(self):
        return sum(self.sizes)

    def cursor(self):
        if not self.levels: return None
        return Cursor(self, len(self.levels)-1, 0, None)
        

class Cursor:
    def __init__(self, tree, level, position, successor):
        self.tree = tree
        self.level = level
        #self.position = position
        self.successor = successor
        
        #self.data_ptr = tree.maps[level] + position
        self.data = tree.levels[level]
        self.offset = position
        
        if level:
            self.left_loc, self.offset = read_int(self.data, self.offset)
            self.right_loc, self.offset = read_int(self.data, self.offset)
            self.right_loc += self.left_loc
        self.key_length, self.offset = read_int(self.data, self.offset)

    def __cmp__(self, other):
        return cmp(self.key(), other.key())
        
    compare = __cmp__

    def compare_to_str(self, string):        
        return cmp(self.key(), string)
        
    def key(self):
        return self.data[self.offset:self.offset+self.key_length]
        
    def value(self):
        value_length, value_offset = read_int(self.data,self.offset+self.key_length)
        return self.data[value_offset:value_offset+value_length]

    def left(self):
        if self.level == 0: return None
        
        return Cursor(self.tree, self.level-1, self.left_loc, self)

    def right(self):
        if self.level == 0: return None

        new_level = self.level-1
        new_position = self.right_loc
        
        while new_position == self.tree.sizes[new_level] - 8:
            if new_level == 0: return None
        
            new_position = int64_from_string(self.tree.levels[new_level][new_position:new_position+8])
            new_level -= 1
        
        return Cursor(self.tree, new_level, new_position, 
                      self.successor)
    
    def next(self):
        right = self.right()
        if right is None: return self.successor
        return left_most(right)


def left_most(cursor):
    if not cursor: return None
    while cursor.level > 0:
        cursor = cursor.left()
    return cursor


class Cursor_merge(object):
    def __init__(self, cursor1, cursor2, reducer):
        self.cursor1 = cursor1
        self.cursor2 = cursor2
        self.reducer = reducer

    def __iter__(self):
        return self

    def next(self):        
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
        
        comparison = self.cursor1.compare(self.cursor2)
        
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


class Cursor_merge_all:
    # cursor_heap = heap of (cursor, priority)
    
    def __init__(self, cursors, reducer):
        # Cursors should be passed as a list from newest to oldest

        self.reducer = reducer    
        self.cursor_heap = [ ]
        
        for i, cursor in enumerate(cursors):
            heapq.heappush(self.cursor_heap, (cursor, i))
    
    def __iter__(self):
        return self

    def next(self):
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


class Skip_discardables:
    def __init__(self, cursor_merger, is_discardable):
        self.cursor_merger = cursor_merger
        self.is_discardable = is_discardable
    
    def __iter__(self):
        return self
    
    def next(self):
        result = self.cursor_merger.next() #May throw Stop_iteration, that's fine
        while self.is_discardable(result[1]):
            result = self.cursor_merger.next()
        return result




def _default_is_final(item): return True
def _default_is_discardable(item): return False
def _default_reducer(item_older, item_newer): return item_newer

class Store:
    """ Key-value store with per-key reductions
    
        Default behaviour is to over-write existing keys,
        and to return None if a key is not present
    
        Order is important!
        reducer takes arguments like this: reducer(older item, newer item)
    
    """

    #cdef public object path, trees, cursors, next_tree_id, read_only
    #cdef public object memory, memory_size
    #cdef public object lock, busy_paths,
    ## merger_processes
    ##cdef public object manager, worker_pool, unbusy_queue
    #cdef public object pool
    #
    #cdef public object default_value, is_discardable, is_final, reducer

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
        
        self.trees = [ ]
        
        self._clean_house()
        self._remap()

    def close(self):
        self.flush()
        self.finish_merging()
    
    def destroy(self):
        """ Delete entire store directory """
        self.clear()
        self.close()
        os.rmdir(self.path)
    
    def _clean_house(self):
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
                
    def _remap(self):
        #self.lock.acquire()
    
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
                
        #self.lock.release()
        
        self.cursors = [ ]
        for tree in self.trees[::-1]: #Note reverse order!
            cursor = tree.cursor()
            if cursor is not None: self.cursors.append(cursor)
        
    
    def clear(self):
        self.finish_merging()
    
        for tree in self.trees:
            remove_tree(tree.path)
        self.memory = { }
                
        self._remap()

    
    def finish_merging(self, wait=True):
        return False
        
    
    def start_merging(self, all=False):
        self.finish_merging(False)
        
        pos = len(self.trees)-2
        while pos >= 0:
            if (all or self.trees[pos].size() <= self.trees[pos+1].size() * 2):
                path1 = self.trees[pos].path
                path2 = self.trees[pos+1].path
                
                if pos == 0:
                    is_discardable_if_needed = self.is_discardable
                else:
                    is_discardable_if_needed = None
                
                do_mergedown(self.reducer, self.path, path1, path2, is_discardable_if_needed)
                
                pos -= 2
            else:
                pos -= 1
        
        self._remap()


    def optimize(self):
        self.flush()
    
        #Not optimal? Need a way to wait for any one process to finish.
        while len(self.trees) > 1:
            self.start_merging(True)
            self.finish_merging(True)        

    
    def flush(self):
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
        #if self.busy_paths:
        #    if self.finish_merging(False):
        #        self.start_merging()
        
        if key in self.memory:
            result = self.memory[key]
            if self.is_final(result): return result
        else:
            result = None
            
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
            
            

def do_mergedown(reducer, path, path1, path2, is_discardable):
    name1 = os.path.split(path1)[1]
    name2 = os.path.split(path2)[1]
    
    new_filename_base = os.path.join(path, name1+'-'+name2)
    new_filename_build = new_filename_base + '-build'
    new_filename_ready = new_filename_base + '-ready'
    
    #lock.acquire()
    tree1 = Tree(path1)
    tree2 = Tree(path2)
    #lock.release()
    
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
        
    #lock.acquire()
    os.rename(new_filename_build, new_filename_ready)
    remove_tree(path1)
    remove_tree(path2)
    os.rename(new_filename_ready, path1)
    #lock.release()


_counting_store_default = string_from_int64(0)

def _counting_store_is_final(item): 
    return False

def _counting_store_is_discardable(item): 
    return item == _counting_store_default

def _counting_store_reduce(item_older, item_newer):
    return string_from_int64(int64_from_string(item_older) + int64_from_string(item_newer))

def _counting_store_iter_filter(pair):
    return (pair[0], int64_from_string(pair[1]))

class Counting_store(Store):
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

def _count_and_shadow_store_is_final(item): 
    return False

# Discard shadow if count == 0
def _count_and_shadow_store_is_discardable(item): 
    return item[:8] == _counting_store_default

def _count_and_shadow_store_reduce(item_older, item_newer):
    older_count = int64_from_string(item_older[0:8])
    older_shadow = int64_from_string(item_older[8:16])
    newer_count = int64_from_string(item_newer[0:8])
    newer_shadow = int64_from_string(item_newer[8:16])
    return string_from_int64(older_count+newer_count) + string_from_int64(max(older_shadow,newer_shadow))

def _count_and_shadow_decode(item):
    count = int64_from_string(item[0:8])
    shadow = int64_from_string(item[8:16])
    return (count, shadow)

def _count_and_shadow_store_iter_filter(pair):
    return (pair[0], _count_and_shadow_decode(pair[1]))

class Count_and_shadow_store(Store):
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




def _general_store_is_discardable(item): 
    return item == ''
    
def _general_store_is_final(item): 
    return item == '' or item[0] == 'S'

def _general_store_reduce(item_older, item_newer):
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

def _general_store_iter_filter(pair):
    return (pair[0], pair[1][1:])        

def _general_store_iter_keys_filter(pair):
    return pair[0]        

def _identity(x): return x

class General_store(Store):
    """ Store that allows setting, deleting and appending 
    
        You can specify functions to encode and decode values.
        If you will be using append_to, the decoder
        should cope with encoded values being concatenated together.
    """

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



def _int_list_store_encode(values):
    return ''.join([ string_from_int64(item) for item in values ])
    
def _int_list_store_decode(string):
    return [ int64_from_string(string[i:i+8]) for i in xrange(0,len(string),8) ]

class Int_list_store(General_store):
    def __init__(self, path, **kwargs):
        General_store.__init__(self, path,
            value_encoder=_int_list_store_encode,
            value_decoder=_int_list_store_decode,
            **kwargs)



