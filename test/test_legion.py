

from nesoni import *
from nesoni import legion

print legion.__name__

legion.coordinator().set_cores(1)

def main():
    @parallel_for(xrange(3))
    def _(item):
        print item
        @parallel_for(xrange(3))
        def _(item2):
            print item, item2
    
if __name__ == '__main__':
    run_script(main)
