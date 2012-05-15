"""

    multiprocessing uses threading in Queue and Pool, raising GIL issues.

"""

try:
    import multiprocessing
except ImportError:
    import processing as multiprocessing
    
    #Hmm
    multiprocessing.Process.is_alive = multiprocessing.Process.isAlive

import time

Lock = multiprocessing.Lock

from Queue import Empty

TIMEOUT = 10.0

def work(func, inherited_args, queue_in, queue_out, semaphore, task):
    """ Process at least one task,
    
        for further tasks, you'll need to take a token. """
    try:
        queue_out.put( func(*(inherited_args+task)) )
        
        while True:
            #Here's a token. Grab it and give me a job.
            semaphore.release()
            try:
                task = queue_in.get(True, TIMEOUT)
            except Empty:
                #Meh, I give up, better get my token back
                if semaphore.acquire(False):
                    break
                else:
                    #I can't have it back? There must be a job.
                    task = queue_in.get()
            
            queue_out.put( func(*(inherited_args+task)) )
    except KeyboardInterrupt:
        pass #ssh

class Worker_pool:
    def __init__(self, func, inherited_args):
        self.func = func
        self.inherited_args = inherited_args
        self.semaphore = multiprocessing.Semaphore(0)
        self.queue_in = multiprocessing.Queue()
        self.queue_out = multiprocessing.Queue()
        self.out_list = [ ]
        self.n_ready = 0

    def get(self, block=True, timeout=None):
        return self.queue_out.get(block,timeout)

    def put(self, task):
        #Grab a token
        if not self.semaphore.acquire(False):
            #No token available right now?
            multiprocessing.Process(
                target=work,
                args=(self.func, self.inherited_args,
                      self.queue_in, self.queue_out, self.semaphore,
                      task),
            ).start()
            return

        self.queue_in.put(task)
        
        
def test_func():
    print 'Hello'
    time.sleep(1.0)
    print 'world'

if __name__ == '__main__':        
    w = Worker_pool(test_func, ())
    for i in xrange(100):
        w.put(())
        time.sleep(0.1)
    w.get()



