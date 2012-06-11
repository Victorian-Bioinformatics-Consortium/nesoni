
import sys, os
sys.path.insert(0, os.path.join(os.path.split(__file__)[0], '../'))

import treemaker

import tokenize, random

def test_counting_store():
    print 'Testing Counting_store'

    counter = treemaker.Counting_store('test_count', memory_size = 10)
    counter.clear()

    check = { }

    print ' - Adding some words'
    for token in tokenize.generate_tokens(open(__file__,'rU').readline):
        if token[1].strip():
            x = token[1]
            counter.add( x )
            check[x] = check.get(x,0)+1

    print ' - Checking it has the right counts'
    for key in check:
        assert check[key] == counter[key]
    
    print ' - Checking it iterates properly'
    for key, value in counter.iter_from(''):
        assert value == counter[key]
        assert value == check[key]
        #print key, value, check[key]

    print ' - Checking for absence of "hammertime"'
    assert counter['hammertime'] == 0
    assert 'hammertime' not in counter
    assert 'for' in counter

    counter.close()

def test_general_store():
    print 'Testing General_store'

    store = treemaker.General_store('test_general_store', memory_size = 10)
    
    check = { }

    print ' - Manipulating the store'

    for i in xrange(10000):
        key = str(random.randrange(100))
        action = random.randrange(3)
        if action == 0:
            del store[key]
            if key in check: del check[key]
        elif action == 1:
            what = str(random.random())
            store.append_to(key, what)
            check[key] = check.get(key,'') + what
        else:
            what = str(random.random())
            store[key] = what
            check[key] = what

    print ' - Checking it has all the right keys'

    for key in check:
        assert store[key] == check[key]
        assert key in store
    
    print ' - Checking it iterates properly'
    
    for key, value in store.iter_from(''):
        assert value == store[key]
        assert value == check[key]
        assert key in store
        del check[key]
    assert not check
    
    print ' - Checking for absence of "hammertime"'
    
    assert 'hammertime' not in store
    try:
        store['hammertime']
        assert False
    except KeyError:
        pass

if __name__ == '__main__':
    test_counting_store()
    test_general_store()



