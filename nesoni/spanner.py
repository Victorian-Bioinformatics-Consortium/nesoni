
import heapq


def iter_to_spanner(an_iterable):
    an_iter = iter(an_iterable)

    pos = 0
    start = 0
    try:
        current = an_iter.next()
        pos += 1
    except StopIteration:
        return
    
    for item in an_iter:
        if current != item:
            yield (pos-start,current)
            current = item
            start = pos
        pos += 1
    if pos != start:
        yield (pos-start,current)

def spanner_to_iter(a_spanner):
    for length, value in a_spanner:
        for i in xrange(length):
            yield value

def zip_spanners(*spanners):
    n = len(spanners)
    
    start = 0
    current = [ ]
    next = [ ]
    
    for i in xrange(n):
        try:
            length, value = spanners[i].next()
            current.append(value)
            heapq.heappush(next, (length,i))
        except StopIteration:
            current.append(None)
    
    while next:
        pos, i = heapq.heappop(next)
        if pos != start:
            yield (pos-start, tuple(current))
            start = pos
        try:
            length, value = spanners[i].next()
            current[i] = value
            heapq.heappush(next,(pos+length, i))
        except StopIteration:
            current[i] = None

def map_spanner(a_func, a_spanner):
    for length, value in a_spanner:
        yield length, a_func(value)


 
