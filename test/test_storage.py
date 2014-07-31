
import os, random

import numpy as np

from nesoni import storage

#s = storage.Storage('test.store', clear=True, preallocate=1<<20)
#
#val = s.data
#val.foo = s.empty(5, 'int')
#val.foo[1] = 2
##val.bar = s.array([1.0,2.0,3.0])
#val.x = 'hi'
#s.save()
#
#s.close()
#
#print '==='
#
#s = storage.Storage('test.store')
#val = s.data
##print repr( s.get().bar )
#print s.data.x
#print repr( s.data.foo )
#
#s.close()
#

#print s.get().x
#print repr( s.get().foo )

s = storage.Storage('output/test.store', clear=True)

a = [ ]
b = [ ]
s.data = a

for i in xrange(10000):
    j = random.choice([0] + ([1] if a else []))
    if j == 0:
        size = random.randrange(1,1000)
        def op(alist, allocator):
            alist.append( allocator.zeros(size, 'int32') )
    elif j == 1:
        k = random.randrange(len(a))
        l = random.randrange(len(a[k]))
        m = random.randrange(l,len(a[k]))
        def op(alist, allocator):
            alist[k][l:m] += 1
    
    #print a
    #print b  
    op(a,s)
    op(b,np)
    for aa,bb in zip(a,b):
        assert np.all(aa == bb)
        
    if i % 1000 == 0:
        print 'reload'
        s.save()
        s.close()
        s = storage.Storage('output/test.store')
        a = s.data



s.save()
s.close()


