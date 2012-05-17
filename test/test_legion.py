
import time

from nesoni import *
from nesoni import legion, grace

class Timer(object):
    def __init__(self, name):
        self.name = name
    
    def __enter__(self):
        self.t1 = time.time()
    
    def __exit__(self, *etc):
        t2 = time.time()
        print self.name, t2-self.t1

def main():
    with Timer('100 status updates'):
        for i in xrange(100):
            grace.status(str(i))
        grace.status('')

    with Timer('100 parallel processes'):
        @parallel_for(xrange(100))
        def loop(i):
            pass

    with Timer('100 parallel processes updating statuses'):
        @parallel_for(xrange(100))
        def loop(i):
            grace.status(str(i))
            grace.status('')

    with Timer('Nested processes'):
        @parallel_for(xrange(3))
        def _(item):
            print item
            @parallel_for(xrange(3))
            def _(item2):
                print item, item2
    
if __name__ == '__main__':
    run_script(main)
