
import time

import nesoni
from nesoni import legion, grace

class Timer(object):
    def __init__(self, name):
        self.name = name
    
    def __enter__(self):
        self.t1 = time.time()
    
    def __exit__(self, *etc):
        t2 = time.time()
        print self.name, t2-self.t1

def hello():
    print 'Hello world.'


def func(name, deps):
    for item in deps:
        item()
    print name 

def do_nothing(*args):
    pass

def do_status(i):
    grace.status(str(i))
    grace.status('')

def main():
    #print dir(hello)
    import sys
    print sys.modules['__main__'].__file__
    
    a = nesoni.future(func,'a',[])
    b = nesoni.future(func,'b',[])
    c = nesoni.future(func,'c',[a,b])
    d = nesoni.future(func,'d',[a,b])
    
    c()
    d()

    #print legion.coordinator().get_cores()
    #legion.coordinator().job(hello)

    with Timer('100 status updates'):
        for i in xrange(100):
            grace.status(str(i))
        grace.status('')

    with Timer('100 parallel threads'):
        @nesoni.thread_for(xrange(100))
        def loop(i):
            pass

    with Timer('100 parallel processes'):
        nesoni.parallel_for(xrange(100))(do_nothing)

    with Timer('100 parallel processes updating statuses'):
        nesoni.parallel_for(xrange(100))(do_status)
    
    #with Timer('Nested processes'):
    #    @nesoni.parallel_for(xrange(3))
    #    def _(item):
    #        print item
    #        @nesoni.parallel_for(xrange(3))
    #        def _(item2):
    #            print item, item2

if __name__ == '__main__':
    nesoni.run_script(main)
