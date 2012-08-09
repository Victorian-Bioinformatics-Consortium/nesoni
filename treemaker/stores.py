"""

    Key value store, with fast writing, semigroup update semantics

"""
    
import sys, os, mmap, heapq, itertools, struct


def write_int64(f, value):
    f.write(struct.pack('<q', value))

def string_from_int64(value):
    return struct.pack('<q', value)

def int64_from_string(string):
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


def write_string(f, value):
    write_int(f, len(value))
    f.write(value)

def read_string(data, offset):
    length, offset = read_int(data, offset)
    return data[offset:offset+length], offset+length


def remove_tree(path):
    if not os.path.exists(path): return
    
    for filename in os.listdir(path):
        assert filename.isdigit()
        os.unlink(os.path.join(path, filename))
    
    os.rmdir(path)


def make_tree(path, sequence, key_writer, value_writer):
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
            
        #write_int(files[ones], len(key))
        #files[ones].write(key)        

        #write_int(files[ones], len(value))
        #files[ones].write(value)
        
        key_writer(files[ones], key)
        value_writer(files[ones], value)
        
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


class Tree(object):
    def __init__(self, path, key_reader, value_reader):
        self.path = path
        self.key_reader = key_reader
        self.value_reader = value_reader
        
        self.levels = [ ]
        while True:
             filename = os.path.join(path, str(len(self.levels)))
             if not os.path.exists(filename): break
             f = open(filename,'rb')
             self.levels.append(mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ))
             f.close()

        self.sizes = [ len(item) for item in self.levels ]      

    def __del__(self):
        for level in self.levels:
            level.close()

    def size(self):
        return sum(self.sizes)

    def cursor(self):
        if not self.levels: return None
        return Cursor(self, len(self.levels)-1, 0, None)
        

class Cursor(object):
    def __init__(self, tree, level, position, successor):
        self.tree = tree
        self.level = level
        self.successor = successor
        
        self.data = tree.levels[level]
        self.offset = position        
        
        if level:
            self.left_loc, self.offset = read_int(self.data, self.offset)
            self.right_loc, self.offset = read_int(self.data, self.offset)
            self.right_loc += self.left_loc

        self.key, self.offset = tree.key_reader(self.data, self.offset)
        
        self._left = None
        self._right = None

    def __cmp__(self, other):
        return cmp(self.key, other.key)
        
    compare = __cmp__

    def compare_to_key(self, key):        
        return cmp(self.key, key)
        
    def get_value(self):
        return self.tree.value_reader(self.data, self.offset)[0]

    def left(self):
        if self.level == 0: return None
        if self._left: return self._left
        
        result = Cursor(self.tree, self.level-1, self.left_loc, self)
        if self.level >= len(self.tree.levels)-1-20:
            self._left = result
        return result

    def right(self):
        if self.level == 0: return None
        if self._right: return self._right

        new_level = self.level-1
        new_position = self.right_loc
        
        while new_position == self.tree.sizes[new_level] - 8:
            if new_level == 0: return None
        
            new_position = int64_from_string(self.tree.levels[new_level][new_position:new_position+8])
            new_level -= 1
        
        result = Cursor(self.tree, new_level, new_position, 
                      self.successor)
        if self.level >= len(self.tree.levels)-1-20:
            self._right = result
        return result
    
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
                
            result = (self.cursor2.key, self.cursor2.get_value())
            self.cursor2 = self.cursor2.next()
            return result

        if self.cursor2 is None:
            result = (self.cursor1.key, self.cursor1.get_value())
            self.cursor1 = self.cursor1.next()
            return result
        
        comparison = self.cursor1.compare(self.cursor2)
        
        if comparison < 0:
            result = (self.cursor1.key, self.cursor1.get_value())
            self.cursor1 = self.cursor1.next()
            return result

        elif comparison > 0:
            result = (self.cursor2.key, self.cursor2.get_value())
            self.cursor2 = self.cursor2.next()
            return result

        else:
            result = (self.cursor1.key, self.reducer(self.cursor1.get_value(),self.cursor2.get_value()))
            self.cursor1 = self.cursor1.next()
            self.cursor2 = self.cursor2.next()
            return result


class Cursor_merge_all(object):
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
                
        result = heap_items[0][0].get_value()
        for cursor, oldness in heap_items[1:]:
            result = self.reducer( cursor.get_value(), result )
        return (heap_items[0][0].key, result)


class Skip_discardables(object):
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

class Store(object):
    """ Key-value store with per-key reductions
    
        Default behaviour is to over-write existing keys,
        and to return None if a key is not present
    
        Order is important!
        reducer takes arguments like this: reducer(older item, newer item)
    
    """
    
    _default_value = None
    _key_writer = staticmethod( write_string )
    _key_reader = staticmethod( read_string )
    _value_writer = staticmethod( write_string )
    _value_reader = staticmethod( read_string )
    @staticmethod
    def _reducer(older, newer):
        return newer
    @staticmethod
    def _is_final(item):
        return True
    @staticmethod
    def _is_discardable(item):
        return False

    def __init__(self, 
                 path, 
                 memory_size=1<<20, #Formerly: 100000, #Formerly: 500000
                 read_only=False
                 ):
        self.read_only = read_only
        
        self.path = path        
        if not os.path.exists(path):
            os.mkdir(path)
        
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
            self.trees.append( Tree(filename, self._key_reader, self._value_reader) )
                
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
                    is_discardable_if_needed = self._is_discardable
                else:
                    is_discardable_if_needed = None
                
                do_mergedown(self._reducer, self.path, path1, path2, is_discardable_if_needed, self._key_reader, self._key_writer, self._value_reader, self._value_writer)
                
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
            make_tree(tempname, items, self._key_writer, self._value_writer)
        except:
            remove_tree(tempname)
            raise
        
        os.rename(tempname, filename)
        self._remap()        
        self.start_merging()
        
    def _put(self, key, value):
        if key in self.memory:
            self.memory[key] = self._reducer(self.memory[key], value)
        else:
            self.memory[key] = value
            if len(self.memory) >= self.memory_size: 
                self.flush()       

    def __getitem__(self, key):
        return self._get(key)
    
    def _get(self, key):
        have_value = key in self.memory
        if have_value:
            result = self.memory[key]
            if self._is_final(result): return result
            
        for cursor in self.cursors:
            while cursor:
                comparison = cursor.compare_to_key(key)
                if comparison == 0: 
                    this_value = cursor.get_value()
                    if not have_value:
                        result = this_value
                        have_value = True
                    else:
                        result = self._reducer(this_value, result)
                    if self._is_final(result):
                        return result
                    break
                    
                if comparison > 0:
                    cursor = cursor.left()
                else:
                    cursor = cursor.right()
        
        if not have_value:
            result = self._default_value
        
        return result
        
    def __contains__(self, key):
        return self._get(key) != self.default_value

    def iter_from(self, key):
        self.flush() # self.memory is a dict, not ordered, need to flush to disk to sort it
    
        cursors = [ ]
        for cursor in self.cursors:
            while True:
                comparison = cursor.compare_to_key(key)
                if comparison == 0:
                    break
                
                if comparison > 0:
                    new_cursor = cursor.left()
                else:
                    new_cursor = cursor.right()
                    
                if new_cursor is None: break
                
                cursor = new_cursor
            
            if cursor.compare_to_key(key) < 0:
                cursor = cursor.next()
            if cursor is not None:
                cursors.append(cursor)
        
        return Skip_discardables(
            Cursor_merge_all(cursors, self._reducer), 
            self._is_discardable)
            
            

def do_mergedown(reducer, path, path1, path2, is_discardable, key_reader, key_writer, value_reader, value_writer):
    name1 = os.path.split(path1)[1]
    name2 = os.path.split(path2)[1]
    
    new_filename_base = os.path.join(path, name1+'-'+name2)
    new_filename_build = new_filename_base + '-build'
    new_filename_ready = new_filename_base + '-ready'
    
    tree1 = Tree(path1, key_reader, value_reader)
    tree2 = Tree(path2, key_reader, value_reader)
    
    merger = Cursor_merge(
                left_most(tree1.cursor()),
                left_most(tree2.cursor()),
                reducer)
        
    if is_discardable is not None:
        merger = Skip_discardables(merger, is_discardable)
    
    try:
        make_tree(new_filename_build, merger, key_writer, value_writer)
    except:
        remove_tree(new_filename_build)
        raise
        
    os.rename(new_filename_build, new_filename_ready)
    remove_tree(path1)
    remove_tree(path2)
    os.rename(new_filename_ready, path1)


class Counting_store(Store):
    """ Store for counting strings with """

    _default_value = 0
    _value_writer = staticmethod( write_int )
    _value_reader = staticmethod( read_int )
    
    @staticmethod
    def _is_final(value):
        return False
    
    @staticmethod
    def _reducer(older, newer):
        return older + newer
    
    @staticmethod
    def _is_discardable(value):
        return value == 0

    def add(self, key, amount=1):
        self._put(key, amount)


class Count_and_shadow_store(Store):
    """ Store for counting strings with """
    
    _default_value = (0,0)
    
    @staticmethod
    def _value_writer(f, (count,shadow)):
        write_int(f, count)
        write_int(f, shadow)
    
    @staticmethod
    def _value_reader(data, offset):
        count, offset = read_int(data, offset)
        shadow, offset = read_int(data, offset)
        return (count,shadow), offset
    
    @staticmethod
    def _reducer((older_count,older_shadow),(newer_count,newer_shadow)):
        return (older_count+newer_count,max(older_shadow,newer_shadow))
    
    @staticmethod
    def _is_final(value):
        return False
    
    @staticmethod
    def _is_discardable(value):
        return value == (0,0) 
    
    def add(self, key, amount, shadow_amount):
        self._put(key, (amount, shadow_amount))




class General_store(Store):
    """ Store that allows setting, deleting and appending 
    
        You can override _value_encode and _value_decode to encode and decode values.
        If you will be using append_to, the decoder
        should cope with encoded values being concatenated together.
    """

    @staticmethod    
    def _is_discardable(item): 
        return item == ''
    
    @staticmethod
    def _is_final(item): 
        return item == '' or item[0] == 'S'
    
    @staticmethod
    def _reduce(item_older, item_newer):
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

    _value_encoder = staticmethod( lambda x:x )
    _value_decoder = staticmethod( lambda x:x )

    def append_to(self, key, value):
        self._put(key, 'A' + self._value_encoder(value))

    def __setitem__(self, key, value):
        self._put(key, 'S' + self._value_encoder(value))

    def __delitem__(self, key):
        self._put(key, '')

    def __getitem__(self, key):
        result = self._get(key)        
        if result == '': 
            raise KeyError(key)

        return self._value_decoder( result[1:] )

    def get(self, key, default):
        result = self._get(key)        
        if result == '': 
            return default

        return self._value_decoder( result[1:] )
                
    def iter_from(self, key):
        return itertools.imap(
            lambda (key,value): (key,self._value_decoder(value[1:])),
            Store.iter_from(self,key))

    def iter_keys_from(self, key):
        return itertools.imap(
            lambda (key,value): key,
            Store.iter_from(self,key))


class Int_list_store(General_store):
    @staticmethod
    def _value_encoder(values):
        return ''.join([ string_from_int64(item) for item in values ])
    
    @staticmethod    
    def _value_decoder(string):
        return [ int64_from_string(string[i:i+8]) for i in xrange(0,len(string),8) ]


