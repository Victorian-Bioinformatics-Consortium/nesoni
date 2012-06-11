
import sys, os
sys.path.insert(0, os.path.join(os.path.split(__file__)[0], '../'))

#import pyximport
#pyximport.install()

import time, sys, os

try: 
    import bsddb
except:
    bsddb = None

import treemaker

def run_test(title, n, setter, flusher, getter):
    """ Insert pseudo-random values,
        then check that they can be retrieved """
    
    print '%30s' % title,
    sys.stdout.flush()

    store = treemaker.Store('test_store')
    store.clear()
    
    t1 = time.time()
    for i in xrange(n):
        key = str( (i * 302513) % 1000000007 )
        value = key[-1] * 8
        #store[ key ] = key[-1] * 8
        setter(key, value)

    # This is not necessary
    # it just doesn't seem quite fair not to do it
    flusher()

    t2 = time.time()
    
    print '   set: %8.1f sec' % (t2-t1),
    sys.stdout.flush()

    for i in xrange(n):
        key = str( (i * 302513) % 1000000007 )
        value = key[-1] * 8
        assert getter(key) == value
    t3 = time.time()

    print '   get: %8.1f sec' % (t3-t2)
    sys.stdout.flush()
    
    #p.join()
    
    return 0

def do_nothing(): pass

def main():
    n = 250000

    while n <= 1024000000:
        print
        print n, 'items'
        print
        store = treemaker.Store('test_store')
        store.clear()
        run_test('treemaker.Store', n, store._put, store.flush, store.__getitem__)
        store.close()
        store = None

        store = treemaker.General_store('test_store')
        store.clear()
        run_test('treemaker.General_store', n, store.__setitem__, store.flush, store.__getitem__)
        store.close()
        store = None
        
        if n <= 32000000:
            db = { }
            run_test('in memory dict', n, db.__setitem__, do_nothing, db.__getitem__)
            db = None
        
        if n <= 1000000 and bsddb is not None:
            if os.path.exists('test.db'): os.unlink('test.db')
            db = bsddb.btopen('test.db')
            run_test('bsddb btree', n, db.__setitem__, do_nothing, db.__getitem__)
            db.close()
            db = None

            if os.path.exists('test.db'): os.unlink('test.db')
            db = bsddb.hashopen('test.db')
            run_test('bsddb  hash', n, db.__setitem__, do_nothing, db.__getitem__)
            db.close()
            db = None

        n *= 2

if __name__ == '__main__':
    main()


