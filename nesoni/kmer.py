"""

A kmer collection contains

- kmer counts
- read / read pair indexing

tree of counts

file of <class>\tread\t[read]\n


Allows

- depth stats, histogram

- SplitTree output - present/absent/not sure call
- depth cutoff
- tail trimming

- contig output
  - pair seperation estimation per class

- contig output with hints from eg PCR

- graph layout
  - pickle or such
- graph interaction


"""

import sys, os, string, re, cPickle, time, struct, math

from nesoni import io, bio, grace

#try:
#    import pyximport
#    pyximport.install()
#except:
#    raise grace.Error('Couldn\'t install pyximport (from Cython), these tools still need Cython, sorry')
import treemaker

BAD_BASES = re.compile('[^ACGT]')

##Too slow
#def norm(kmer):
#    #rkmer = reverse_complement(kmer)
#    #kmer = min(kmer, rkmer)
#    
#    kp1 = len(kmer) + 1
#    kmer2 = kmer+'-'+kmer
#    best = kmer[:kp1]
#    for i in xrange(1,kp1):
#        if kmer2[i] <= best:    
#            best = min(best, kmer2[i:i+kp1])
#    return best

def norm(kmer):
    rkmer = bio.reverse_complement(kmer)
    return min(kmer, rkmer)

def norm_all(thing):
    if isinstance(thing, str):
        return norm(thing)
    else:
        return [ norm_all(item) for item in thing ]

def reverse_kmer_list(thing):
    return [ bio.reverse_complement(kmer) for kmer in thing[::-1] ]


def kmers_from_sequence(k, sequence):
    result = [ ]
    for good_sequence in BAD_BASES.split(sequence):
        for i in xrange(len(good_sequence)-k+1):
            result.append(good_sequence[i:i+k])
    return result

def sequence_from_kmer_list(kmers):
    for i in xrange(1,len(kmers)):
        assert kmers[i-1][1:] == kmers[i][:-1]
    
    return ''.join( kmer[0] for kmer in kmers ) + kmers[-1][1:]


def successors(kmer):
    km1mer = kmer[1:]
    return [ norm(km1mer+letter) for letter in 'ACGT' ]
def predecessors(kmer):
    km1mer = kmer[:-1]
    return [ norm(letter+km1mer) for letter in 'ACGT' ]
def friends(kmer):
    return successors(kmer) + predecessors(kmer)

def norm_and_sign(kmer):
    rkmer = bio.reverse_complement(kmer)
    if rkmer < kmer:
        return rkmer, -1
    elif rkmer == kmer:
        return kmer, 0
    else:
        return kmer, 1

def friends_with_sign_travel(kmer):
    result = [ ]
    head = kmer[:-1]
    tail = kmer[1:]
    for base in 'ACGT':
        successor, successor_sign = norm_and_sign(tail+base)
        result.append( (successor, successor_sign, 1) )
        predecessor, predecessor_sign = norm_and_sign(base+head)
        result.append( (predecessor, predecessor_sign, -1) )
    return result


def is_hairtip(kmers, kmer):
    n_succ = 0
    for kmer2 in successors(kmer):
        if kmer2 in kmers: n_succ += 1
    n_pred = 0
    for kmer2 in predecessors(kmer):
        if kmer2 in kmers: n_pred += 1
    return n_succ + n_pred < 2
    
def trace_forward(kmers, root, n, equivalent_getter=lambda x:[x]):
    result = [ (0, root) ]
    result_set = { }
    
    switchbacks = [ ]
    
    i = 0
    while i < len(result):
        distance, kmer = result[i]
        if distance > n: break
        i += 1
        
        for kmer2 in equivalent_getter(kmer):
            if kmer2 in result_set: continue

            if kmer2 == bio.reverse_complement(kmer2):
                #print kmer, 'is a palindrome!'
                switchbacks.append(kmer2)
        
            result_set[kmer2] = distance
            
            head = kmer2[1:]
            for ext in 'ACGT':
                kmer3 = head+ext
                if norm(kmer3) in kmers:
                    result.append( (distance+1, kmer3) )
    
    return result_set, switchbacks

def single_strand_forward(kmer, kmers, check_predecessors=True):
    result = [ ]
    result_set = set()
    while True:
        if kmer in result_set: break #Loop!
    
        successors = [ ]
        predecessors = [ ]
        for char in 'ACGT':
            successor = kmer[1:] + char
            if norm(successor) in kmers:
                successors.append(successor)
            
            predecessor = char + kmer[:-1]
            if norm(predecessor) in kmers:
                predecessors.append(predecessor)
        if len(successors) != 1: break
        if check_predecessors and len(predecessors) != 1: break
        
        #***Don't include branch point***
        result.append(kmer)
        result_set.add(kmer)
        kmer = successors[0]
    return result[1:]

def single_strand_from_seed(kmer, kmers):
    f = single_strand_forward(kmer, kmers)
    b = single_strand_forward(bio.reverse_complement(kmer), kmers)
    
    return reverse_kmer_list(b) + [kmer] + f


def greedy_single_strand_forward(kmer, kmers, check_predecessors=True):
    result = [ ]
    result_set = set()
    previous = None
    while True:
        if kmer in result_set: break #Loop!
        result.append(kmer)
        result_set.add(kmer)
    
        successors = [ ]
        predecessors = [ ]
        for char in 'ACGT':
            successor = kmer[1:] + char
            nsuccessor = norm(successor)
            if nsuccessor in kmers:
                successors.append( (kmers[nsuccessor],successor) )
                        
        
        successors.sort(reverse=True)
        if not successors: break
        if len(successors) > 1 and successors[1][0] >= successors[0][0]: break
        
        if check_predecessors:
            predecessors = [ ]
            for char in 'ACGT':
                predecessor = char + kmer[:-1]
                npredecessor = norm(predecessor)
                if npredecessor in kmers:
                    predecessors.append( (kmers[npredecessor],predecessor) )
            
            predecessors.sort(reverse=True)
            if len(predecessors) > 1:
                if previous is not None and predecessors[0][1] != previous: break
                if predecessors[1][0] >= predecessors[0][0]: break
        
        previous = kmer
        kmer = successors[0][1]
        
    return result[1:]

def greedy_single_strand_from_seed(kmer, kmers, check_predecessors=True):
    f = greedy_single_strand_forward(kmer, kmers, check_predecessors)
    b = greedy_single_strand_forward(bio.reverse_complement(kmer), kmers, check_predecessors)
    
    return reverse_kmer_list(b) + [kmer] + f

def greedy_strands(kmers):
    kmers_copy = kmers.copy()
    
    strands = [ ]
    
    for kmer in sorted(kmers, key=lambda kmer: kmers[kmer], reverse=True):
        if kmer not in kmers_copy: continue
        
        strand = greedy_single_strand_from_seed(kmer, kmers_copy, False)
        for kmer2 in strand:
            nkmer2 = norm(kmer2)
            if nkmer2 in kmers_copy:
                del kmers_copy[norm(kmer2)]
            
        if len(strand) > 1:
            strands.append(strand)
    
    strands.sort(key=lambda x: len(x), reverse=True)
    return strands


def find_path(start, end, kmers):
    steps = [ (start, None) ]
    used = set([start])
    i = 0
    while i < len(steps):
        kmer = steps[i][0]
        if kmer == end:
            j = i
            collection = [ ]
            while j != None:
                collection.append(steps[j][0])
                j = steps[j][1]
            return collection
        
        for kmer2 in friends(kmer):
            if kmer2 in kmers and kmer2 not in used:
                used.add(kmer2)
                steps.append( (kmer2, i) )
    
        i += 1
    
    return [ ]



class Kmer_bag(io.Workspace):
    def __init__(self, working_dir, read_only=False):
        io.Workspace.__init__(self, working_dir)

        self.k = self.param.get('k', None)

        sequence_filename = self._object_filename('sequences')
        if not os.path.exists(sequence_filename):
            open(sequence_filename,'wb').close()
        self.sequence_file = open(sequence_filename,'rb+')

        k_to_s_filename = self._object_filename('kmer_to_sequences')
        self.kmer_to_sequences = treemaker.Int_list_store(k_to_s_filename, read_only=read_only)
        
        count_filename = self._object_filename('kmer_counts')
        self.counts = treemaker.Counting_store(count_filename, read_only=read_only)


    def set_k(self, k):
        assert self.k is None, 'k already set, can not be changed'
        
        self.update_param(k=k)
        self.k = k

        
    def clear(self):
        self.counts.clear()
        self.kmer_to_sequences.clear()
        self.sequence_file.seek(0)
        self.sequence_file.truncate()
        
        self.update_param(remove=['k'])
        self.k = None


    def close(self):
        self.sequence_file.close()
        self.counts.close()
        self.kmer_to_sequences.close()

    def load_reads_count_only(self, filenames):
        for filename in filenames:
            grace.status('Loading '+filename)
            
            nth = 0
            for name, sequence in io.read_sequences(filename):
                nth += 1
                if nth % 100000 == 0:
                    grace.status('Loading '+filename+' read '+grace.pretty_number(nth))
            
                for good_sequence in BAD_BASES.split(sequence):
                    for i in xrange(len(good_sequence)-self.k+1):
                        kmer = good_sequence[i:i+self.k]
                        nkmer = norm(kmer)
                        self.counts.add(nkmer, 1)
        grace.status('')

    def get_sequences(self, location):
        self.sequence_file.seek(location)
        line = self.sequence_file.readline().rstrip('\n')
        items = line.split('\t')
        return items

    def insert_sequences(self, sequences):
        """ Insert a read or read pair """
        line = '\t'.join(sequences) + '\n'        

        self.sequence_file.seek(0,2)
        pos = self.sequence_file.tell()        
        self.sequence_file.write(line)
        
        for sequence in sequences:
            for good_sequence in BAD_BASES.split(sequence):
                for i in xrange(len(good_sequence)-self.k+1):
                    kmer = good_sequence[i:i+self.k]
                    nkmer = norm(kmer)
                    self.counts.add(nkmer, 1)
                    self.kmer_to_sequences.append_to(nkmer, [pos])
        

    def load_reads(self, filenames):
        #n = 0
        for filename in filenames:
            grace.status('Loading '+filename)
            for name, sequence in io.read_sequences(filename):
                self.insert_sequences([sequence.upper()])
                
                #n += 1
                #if n%1000==0: print n
                #if n >= 100000: break
        grace.status('')


    def load_pairs(self, filenames):
        temp_store = treemaker.General_store(self._object_filename('temp'))
        try:
            temp_store.clear()

            n = 0
            for filename in filenames:
                grace.status('Reading '+filename)
                for name, sequence in io.read_sequences(filename):
                    temp_store[name] = sequence.upper()
                    n += 1
                    
                    #if n >= 1000000: break
                    
            grace.status('Sorting')
            temp_store.optimize()

            i = 0
            for name, sequence in temp_store.iter_from(''):
                if i % 1000 == 0:
                    grace.status('Loading read %s/%s' % (grace.pretty_number(i),grace.pretty_number(n)))
                i += 1
                    
                if name.endswith('/2'):
                    other_name = name[:-2] + '/1'
                    if other_name not in temp_store:
                        self.insert_sequences([sequence])
                    continue
                
                assert name.endswith('/1'), 'Unexpected pair naming'
                
                other_name = name[:-2] + '/2'
                other_sequence = temp_store.get(other_name, None)
                
                if other_sequence is None:
                    self.insert_sequences([sequence])
                    continue
                
                self.insert_sequences([sequence, 
                                       bio.reverse_complement(other_sequence)]) 
            
        finally:
            grace.status('Cleaning up')
            temp_store.destroy()
        grace.status('')


    def kmers_reaching_cutoff(self, cutoff):
        kmers = { }
        for kmer, count in self.counts.iter_from(''):
            if count >= cutoff:
                kmers[kmer] = count
        return kmers
    

    def stats(self, plot=False):
        print
        print 'k =', self.k
    
        try: import numpy
        except ImportError: import numpypy as numpy
        
        if plot:
            import pylab
    
        histogram = [ ]
        for i, (kmer, count) in enumerate(self.counts.iter_from('')):
            while len(histogram) <= count: histogram.append(0)
            histogram[count] += 1
            
            if (i&65535) == 0:
                grace.status('Scanning ' + kmer + ' max multiplicity: %d' % (len(histogram)-1))
        grace.status('')
        
        unique_total = numpy.sum(histogram)
        
        weighted_histogram = [ histogram[i] * i for i in xrange(len(histogram)) ]
        total = sum(weighted_histogram)
        
        remainder = total*0.5
        N50 = 0
        while remainder > weighted_histogram[N50]:
            remainder -= weighted_histogram[N50]
            N50 += 1
        
        total_multiplicity = sum( weighted_histogram[i]*i for i in xrange(len(histogram)) )
         
        print
        print 'Unique k-mers:        %20s' % grace.pretty_number(unique_total) 
        print 'Total k-mers:         %20s' % grace.pretty_number(total)
        print
        print 'Average multiplicity: %22.1f' % (float(total_multiplicity)/total)
        print 'Median multiplicity:  %20d' % N50
        print 'Maximum multiplicity: %20d' % (len(histogram)-1)
        print
        
        print 'At or above this multiplicity, there are this many unique k-mers'
        cutoff = 1
        while cutoff < len(histogram):
            print '%25d %25d' % (cutoff, sum(histogram[cutoff:]))
            
            f = str(cutoff)[0]
            if f in '15': cutoff *= 2
            else: cutoff = cutoff*5//2
        
        if plot:
            #pylab.loglog(range(1,len(histogram)), weighted_histogram[1:], '-+')
            #pylab.text(N50, weighted_histogram[N50], 'Median multiplicity: %d'%N50)
            #pylab.axvline(N50)
            #pylab.xlabel('Multiplicity')
            #pylab.ylabel('Number of kmers (non-unique)')
            #pylab.show()
            
            pylab.loglog(range(len(histogram)-1,0,-1), numpy.cumsum(histogram[:0:-1]), '-+')
            pylab.xlabel('Multiplicity cutoff')
            pylab.ylabel('Number of unique k-mers')
            pylab.show()

def depth_stats(bags, filenames, format='text'):
    import numpy
    
    assert format in ['text', 'jal'], 'Unknown format'
    
    if format == 'jal':
        colors = [
           ('A', 0x64,0xf7,0x3f),
           ('T', 0x3c,0x88,0xee),
           ('C', 0xff,0xb3,0x40),
           ('G', 0xeb,0x41,0x3c),
           ('N', 0xff,0xff,0xff),
        ]
        for base,r,g,b in colors:
            print base+'depth0\tff0000'
            for i in xrange(31):
                c = 127-min(127.0,i*128 / 24.0)
                print base+'depth%d\t%02x%02x%02x' % (1 << i, int(r*0.5+c), int(g*0.5+c), int(b*0.5+c))
        #print 'depth0\tff0000'
        #for i in xrange(31):
        #    print 'depth%d\t%02x%02x%02x' % (1 << i, max(0,255-i*32), max(0,255-i*16), max(0,255-i*8))
        #print
        
        print 'STARTGROUP\tDepths'
    elif format == 'text':
        print
        print 'k-mer depths:'
        print

    k = bags[0].k

    for filename in filenames:
        for name, seq in io.read_sequences(filename):
            counts = [ ]
            for i in xrange(len(seq)-k+1):
                kmer = seq[i:i+k]
                nkmer = norm(kmer)
                counts.append(sum( bag.counts[nkmer] for bag in bags ))
            counts = numpy.array(counts)
            if format == 'text':
                print name, '%.1f%% above zero, mean %.1f, median %d, min %d, max %d' % (
                    100.0*numpy.average(counts > 0),
                    numpy.average(counts), 
                    numpy.median(counts), 
                    numpy.minimum.reduce(counts), 
                    numpy.maximum.reduce(counts)
                )
                plot = [ ]
                peak = numpy.maximum.reduce(counts)
                for item in counts:
                    plot.append( '_0123456789'[min(10, (item * 10 + (peak-1)) // peak)] )
                print ''.join(plot)
                print
            
            elif format == 'jal':
                dilation = [
                    numpy.average(
                        counts[ max(0,i-k+1) : min(len(counts),i+1) ]
                    )
                    for i in xrange(len(seq))
                ]
                
                for i in xrange(len(seq)):
                    if dilation[i] == 0:
                       n = 0
                    else:
                       n = 1 << min(30, int(numpy.log(dilation[i])/numpy.log(2.0)))
                
                    print 'Average kmer depth %.0f\t%s\t-1\t%d\t%d\t%sdepth%d' % (
                        dilation[i], name.split()[0], i+1, i+1, seq[i], n 
                    )
                
    if format == 'jal':
        print 'ENDGROUP\tDepths'
                    
                    
                


def get_root(kmer, parent):
    root = kmer
    reversed = False
    while root in parent: 
        root, this_reversed = parent[root]
        reversed = reversed ^ this_reversed
    #while kmer in parent:
    #    #kmer, parent[kmer] = parent[kmer], root
    #    k2 = parent[kmer]
    #    parent[kmer] = root
    #    kmer = k2
        
    return root, reversed

def mutate(kmer, n=1, start=0):
    if n < 1: return

    for pos in xrange(start, len(kmer)):
        for base in 'ACGT':
            if base == kmer[pos]: continue
            new_kmer = kmer[:pos]+base+kmer[pos+1:]
            
            yield new_kmer
            if n > 1:
                for new_kmer2 in mutate(new_kmer, n-1, pos+1):
                    yield new_kmer2 


class Kmer_graph(io.Workspace):
    def __init__(self, working_dir):
        io.Workspace.__init__(self, working_dir)

        if 'bags' not in self.param:
            self.update_param(bags=[])
        
        self.bags = { } # path -> instance
        for path in self.param.get('bags',()):
            assert os.path.exists(self.relative_path_as_path(path)), path + ' does not exist'
            self.bags[path] = Kmer_bag(self.relative_path_as_path(path), read_only=True)

    def close(self):
        for bag in self.bags.values():
            bag.close()
    
    def clear(self):
        self.bags = { }        
        self.update_param(bags=[])
    
    def add_bag(self, path):
        path = self.path_as_relative_path(path)

        assert path not in self.bags
        assert os.path.exists(self.relative_path_as_path(path))
        
        self.update_param(bags = self.param['bags'] + [path])
        self.bags[path] = Kmer_bag(self.relative_path_as_path(path))
    
    def make_graph(self, cutoff, seed_filenames=[], minimum_object_size=None, merge_snps=False, merge_indels=False):
        # Graph layout depends on GTK, etc, so not imported at top
        from nesoni import graph_layout
    
        bags = self.bags.values()
        k = bags[0].k #TODO: check k consistency
        
        if len(bags) == 1:
           kmers = bags[0].kmers_reaching_cutoff(cutoff)
        else:
           kmers = { }
           possible_cutoff = (cutoff+len(bags)-1) // len(bags)
           print possible_cutoff
           for bag in bags:
               print 'scanning', bag.working_dir
               for kmer, count in bag.counts.iter_from(''):
                   if count >= possible_cutoff:
                       kmers[kmer] = 0
           
           print grace.pretty_number(len(kmers)), 'possibilities'

           possible = kmers.keys()
           possible.sort()
           for bag in bags:
               for kmer in possible:
                   kmers[kmer] += bag.counts[kmer]
        
           for kmer in kmers.keys():
               if kmers[kmer] < cutoff:
                   del kmers[kmer] 

        print grace.pretty_number(len(kmers)), 'at or above cutoff'
        


        if seed_filenames:
            reachable_kmers = { }
            
            seeds = set()
            for filename in seed_filenames:
                for seq_name, seq in io.read_sequences(filename):
                    seeds.update( norm(kmer) for kmer in kmers_from_sequence(k, seq) )
            
            todo = list(seeds)
            while todo:
                kmer = todo.pop(-1)
                if kmer not in kmers: continue
                if kmer in reachable_kmers: continue
                reachable_kmers[kmer] = kmers[kmer]
                todo.extend(friends(kmer))
            
            kmers = reachable_kmers
            print grace.pretty_number(len(kmers)), 'linked to seed'
        
        if minimum_object_size is not None:
            good_kmers = { }
            seen = set()
            
            n = 0
            for kmer in kmers.keys():
                if kmer in seen: continue
                
                n += 1
                if (n&65535) == 0:
                    grace.status('object sizing %d' % len(seen))
                
                connected_object = set()
                connected_object_size = 0
                todo = [ kmer ]
                while todo:
                    kmer2, kmer2_reversed = parent[ todo.pop(-1) ]
                    if kmer2 not in kmers: continue
                    if kmer2 in connected_object: continue
                    
                    connected_object_size += 1
                    for kmer3 in equivalent_kmers(kmer2):
                        connected_object.add(kmer3)
                        seen.add(kmer3)
                        todo.extend(friends(kmer3))
                
                if connected_object_size >= minimum_object_size:
                    for kmer2 in connected_object:
                        good_kmers[kmer2] = kmers[kmer2]
            
            grace.status('')
            
            kmers = good_kmers
            print grace.pretty_number(len(kmers)), 'in object above cutoff size', minimum_object_size    
            
            #todo = list(kmers)
            #while todo:
            #    kmer = todo.pop(-1)
            #    if kmer not in kmers: continue                
            #    if not is_hairtip(kmers, kmer): continue
            #    
            #    del kmers[kmer]
            #    todo.extend(friends(kmer))
            #
            #print len(kmers), 'when dehaired'

        assert kmers, 'No k-mers'
        
        #==== Refactor me! =====
        
        parent = { } # kmer -> (kmer, reversed)

        if merge_snps:
            grace.status('Merging kmers by SNP')
            
            n = 0
            for kmer in kmers:
                n += 1
                if n % 1000 == 0:
                    grace.status('Merging kmer %s by SNP' % grace.pretty_number(n))
                  
                #for pos in xrange(len(kmer)):
                #    for base in 'ACGT':
                #        if base == kmer[pos]: continue
                #        new_kmer = kmer[:pos]+base+kmer[pos+1:]
                for new_kmer in mutate(kmer, 1):                    
                    new_kmer_r = bio.reverse_complement(new_kmer)
                    new_reversed = new_kmer_r < new_kmer
                    if new_reversed:
                        new_kmer = new_kmer_r
                    
                    if new_kmer not in kmers: continue
                    
                    this_root, this_root_reversed = get_root(kmer, parent)
                    new_root, new_root_reversed = get_root(new_kmer, parent)
                    if this_root != new_root:
                        parent[this_root] = (new_root, new_reversed ^ new_root_reversed ^ this_root_reversed)

        if merge_indels:
            grace.status('Merging kmers by indel')

            n = 0
            for kmer in kmers:
                n += 1
                if n % 1000 == 0:
                    grace.status('Merging kmer %s by indel' % grace.pretty_number(n))
                    
                for pos in xrange(1,k-1):
                    if pos < k//2:
                        a = kmer[1:pos]
                        b = kmer[pos:]
                    else:
                        a = kmer[:pos]
                        b = kmer[pos:-1]
                    for base in 'ACGT':
                        new_kmer = a+base+b                        
                        new_kmer_r = bio.reverse_complement(new_kmer)
                        new_reversed = new_kmer_r < new_kmer
                        norm_new_kmer = new_kmer_r if new_reversed else new_kmer
                        if norm_new_kmer not in kmers: continue
                        
                        if new_kmer[1:] == kmer[:-1] or new_kmer[:-1] == kmer[1:]: continue
                    
                        this_root, this_root_reveresed = get_root(kmer, parent)
                        new_root, new_root_reversed = get_root(norm_new_kmer, parent)
                        if this_root != new_root:
                            parent[this_root] = (new_root, new_reversed ^ new_root_reversed ^ this_root_reversed) 
        
        grace.status('Finding kmer-sets')
        children = { }  # kmer -> [ (kmer, reversed) ]
        for kmer in kmers:
            root, root_reversed = get_root(kmer, parent)
            if root not in children: children[root] = [ ]
            children[root].append((kmer, root_reversed))
            
        for kmer in children:
            children[kmer].sort(key=lambda x: kmers[x[0]], reverse=True) #Sort highest to lowest depth
            for child, child_reversed in children[kmer]:
                parent[child] = (kmer, child_reversed)
        
        # We now have that:
        # - children maps from roots to a sorted list of kmers, in reverse order of depth
        # - parent points directly to the root, even root points to root        
        
        grace.status('')
        
        
        def equivalent_kmers(kmer):
            return [ item[0] for item in children[get_root(kmer, parent)] ]
        
        #if 0:
        #    indominable_kmers = { }
        #
        #    print 'puff'
        #    puff = kmers.copy()
        #    for kmer in kmers:
        #        bad = False
        #        for pos in xrange(len(kmer)):
        #            for base in 'ACGT':
        #                if base == kmer[pos]: continue
        #                new_kmer = norm(kmer[:pos]+base+kmer[pos+1:])
        #                puff[norm(new_kmer)] = max(puff.get(new_kmer,0),kmers[kmer])
        #
        #    print 'check'
        #    for kmer in kmers:
        #        bad = False
        #        for pos in xrange(len(kmer)):
        #            for base in 'ACGT':
        #                if base == kmer[pos]: continue
        #                new_kmer = kmer[:pos]+base+kmer[pos+1:]
        #                if puff.get(norm(new_kmer),0) > kmers[kmer]:
        #                    bad = True
        #                    break
        #                
        #                #for pos2 in xrange(pos+1,self.k):
        #                #    for base2 in 'ACGT':
        #                #        if new_kmer[pos2] == base2: continue
        #                #        new_kmer2 = new_kmer[:pos]+base2+new_kmer[pos+1:]
        #                #        if kmers.get(norm(new_kmer2),0) > kmers[kmer]:
        #                #            bad = True
        #                #            break
        #                if bad: break
        #                
        #                
        #            if bad: break
        #        if not bad:
        #            indominable_kmers[kmer] = kmers[kmer]
        #    #kmers = indominable_kmers
        #    print len(indominable_kmers), 'undominated by 2NP kmer'
        #    
        #    knowaguy_kmers = { }
        #    todo = list(indominable_kmers.keys())
        #    while todo:
        #        kmer = todo.pop(-1)
        #        if kmer in knowaguy_kmers: continue
        #        knowaguy_kmers[kmer] = kmers[kmer]
        #        
        #        for kmer2 in (kmer, bio.reverse_complement(kmer)):
        #            successors_present = [ ]
        #            for kmer3 in successors(kmer2):
        #                if kmer3 in kmers:
        #                    successors_present.append(kmer3)
        #            if len(successors_present) == 1:
        #                todo.append(successors_present[0])
        #    
        #    kmers = knowaguy_kmers
        #    print len(kmers), 'undominated, or simply linked to undominated'
        
            
            
        #names = kmers.keys()
        #names.sort()
        #graph = graph_layout.Graph(names)
        #for kmer1 in kmers:
        #    for kmer2 in friends(kmer1):
        #        if kmer2 > kmer1 and kmer2 in kmers:
        #            graph.link(kmer1, kmer2)
        names = children.keys()
        names.sort()
        
        graph = graph_layout.Graph(
            names,
            [ sum( kmers[kmer2] for kmer2, kmer2_reversed in children[kmer] ) for kmer in names ]
        )
        
        for kmer1 in children:
            friend_set = set()
            for kmer1a, kmer1a_reversed in children[kmer1]:
                for kmer2, kmer2_sign, kmer2_travel in friends_with_sign_travel(kmer1a):
                    if kmer2 in kmers:
                        kmer2_root, kmer2_root_reversed = parent[kmer2]
                        if kmer2_root > kmer1:
                            if kmer1a_reversed:
                                kmer2_root_travel = -kmer2_travel
                            else:
                                kmer2_root_travel = kmer2_travel
                            if kmer2_root_reversed ^ kmer1a_reversed:
                                kmer2_root_sign = -kmer2_sign
                            else:
                                kmer2_root_sign = kmer2_sign 
                            friend_set.add( (kmer2_root, kmer2_root_sign, kmer2_root_travel) )
            for kmer2, kmer2_sign, kmer2_travel in friend_set:
                graph.link(kmer1, kmer2, kmer2_sign, kmer2_travel)
                        


        #Link SNPs!
        #for kmer in kmers:
        #    for pos in xrange(len(kmer)):
        #        for base in 'ACGT':
        #            if base == kmer[pos]: continue
        #            new_kmer = norm(kmer[:pos]+base+kmer[pos+1:])
        #            if new_kmer < kmer: continue
        #            if new_kmer not in kmers: continue
        #            
        #            graph.link(kmer, new_kmer)
                    
        self.set_object(graph, 'graph')
        self.set_object(kmers, 'graph_kmer_counts')
        
        self.set_object(parent, 'graph_parents')
        self.set_object(children, 'graph_children')
        
        
        #TODO: fixme
        bodgy_kmers = { }
        for kmer in children:
            kmer2 = children[kmer][0][0]
            #for k,r in children[kmer]:
            #    print k, kmers[k]
            #print
            bodgy_kmers[kmer2] = kmers[kmer2]
        self.create_strands(bodgy_kmers)

    def create_strands(self, kmers):
        import numpy
    
        filename = self._object_filename('strands.fa')
        f = open(filename, 'wb')
        
        strands = [ ]
        
        strands = greedy_strands(kmers)        
        
        #kmers_copy = kmers.copy()
        #
        #for kmer in sorted(kmers, key=lambda kmer: kmers[kmer], reverse=True):
        #    if kmer not in kmers_copy: continue
        #    
        #    strand = greedy_single_strand_from_seed(kmer, kmers_copy, False)
        #    for kmer2 in strand:
        #        nkmer2 = norm(kmer2)
        #        if nkmer2 in kmers_copy:
        #            del kmers_copy[norm(kmer2)]
        #        
        #    if len(strand) > 1:
        #        strands.append(strand)
        #        grace.status('%d strands' % len(strands))
        #
        #grace.status('')        
        #
        #strands.sort(key=lambda x: len(x), reverse=True)
        #
        for i, strand in enumerate(strands):
            depths = [ kmers[norm(kmer2)] for kmer2 in strand ]
            avg_depth = numpy.average(depths)
            min_depth = numpy.minimum.reduce(depths)
            max_depth = numpy.maximum.reduce(depths)
            
            seq = sequence_from_kmer_list(strand)
            io.write_fasta(f, 
                           'strand%d %dbp depth avg %.1f min %d max %d' % 
                               (i+1, len(seq), avg_depth, min_depth, max_depth), 
                           seq)        
        
        f.close()
        
    def graph_layout(self):
        graph = self.get_object('graph')
        for i in xrange(10000):
            grace.status('Iterating %d %f' % (i,graph.update_amount))
            graph.iterate()
            self.set_object(graph,'graph')
        grace.status('')
    
    def graph_show(self, overlay_filenames, bag_to_compare=None):
        # Graph layout depends on GTK, etc, so not imported at top
        from nesoni import graph_layout
    
        graph = self.get_object('graph')
        
        #weights = None #self.get_object('graph_kmer_counts')
        
        self.kmer_counts = self.get_object('graph_kmer_counts')
        self.children = self.get_object('graph_children')
        self.parent = self.get_object('graph_parents')
        
        #weights = { }
        #for kmer in self.children:
        #    weights[kmer] = sum( self.kmer_counts[kmer2] for kmer2, kmer2_reversed in self.children[kmer] )

        viewer = graph_layout.Graph_viewer(graph,
            actions=['read/pair connections', 'graph trace'],
            callback=self._graph_viewer_callback)

        bags = self.bags.values()
        #print bags[0].working_dir
        #print bags[1].working_dir
        #
        #weights = { }
        #for name in graph.names:
        #    ac = bags[0].counts[name]
        #    bc = bags[1].counts[name]
        #    
        #    n = ac+bc
        #    
        #    #r = float(ac)/(ac+bc)
        #    #viewer.annotate(name, 50.0, r,1-r,0)
        #    if ac>=10 and not bc:
        #        viewer.annotate(name, 50.0, 1,0,0)
        #    if bc>=10 and not ac:
        #        viewer.annotate(name, 50.0, 0,1,0)
        #viewer.refresh_annotation()
        
        
        n = 0
        for filename in overlay_filenames:
            for name, seq in io.read_sequences(filename):
                kmers = kmers_from_sequence(bags[0].k, seq.upper())
                nkmers = norm_all(kmers)
                good_nkmers = [ self.parent[nkmer][0] for nkmer in nkmers if nkmer in self.parent ]
                viewer.arrow(good_nkmers, name)
                
                angle = 3.88322208 * n #2*pi/golden ratio
                n += 1
                r = math.sin(angle)*0.5+0.5
                g = math.sin(angle+math.pi/3.0)*0.5+0.5
                b = math.sin(angle-math.pi/3.0)*0.5+0.5
                
                #n_snp = 0
                for nkmer in good_nkmers:
                    viewer.annotate(self.parent[nkmer][0], 4.0, r,g,b)
                    #else:
                    #    for i in xrange(len(nkmer)):
                    #        for base in 'ACGT':
                    #            if nkmer[i] != base:
                    #                nkmer2 = norm(nkmer[:i] + base + nkmer[i+1:])
                    #                if graph.has(nkmer2):
                    #                    viewer.annotate(nkmer2, 2.0, 1.0,0.0,0.0)
                    #                    n_snp += 1
                
                print name
                print '  ', grace.pretty_number(len(nkmers)), 'kmers'
                print '  ', grace.pretty_number(len(good_nkmers)), 'present'
                #print '  ', grace.pretty_number(n_good+n_snp), 'allowing 1 change'                        
                                        
                
                #if len(kmers) < 100: continue
                
                #for i in [ len(kmers)//2 ]: #xrange(0, len(kmers), 100):
                #    nkmer = norm(kmers[i])
                #    if nkmer in graph.name_to_ident:
                #        viewer.label(nkmer, name+' %d'%i)



        if bag_to_compare is not None:
            bag = Kmer_bag(bag_to_compare, read_only=True)
            
            kmers_copy = self.kmer_counts.copy()
            for kmer, count in bag.counts.iter_from(''):
                if kmer in kmers_copy:
                    del kmers_copy[kmer]
            for kmer in kmers_copy:
                viewer.annotate(kmer, 8.0, 1.0,0.0,0.0)
            
            bag.close()


        viewer.refresh_annotation()        
    
        viewer.run()

    def _equivalents(self, kmer):
        rkmer = bio.reverse_complement(kmer)
        reverse = rkmer < kmer
        if reverse:
            nkmer = rkmer
        else:
            nkmer = kmer
        
        root, root_reverse = self.parent[nkmer]
        for child, child_reverse in self.children[root]:
            if child_reverse ^ root_reverse ^ reverse:
                yield bio.reverse_complement(child)
            else:
                yield child

    def _root_norm(self, kmer):
        return self.parent[norm(kmer)][0]
    
    def _graph_viewer_callback(self, viewer, action, kmer1, kmer2):
        graph = viewer.graph

        print
        print kmer1 + ' -- ' + kmer2
        collection = find_path(kmer1, kmer2, self.kmer_counts)
        strands = greedy_strands(dict( (kmer,self.kmer_counts[kmer]) for kmer in collection ))
        for strand in strands:
            norm_strand = norm_all(strand)
            direction_evidence = 0
            if kmer1 in norm_strand:
                pos = norm_strand.index(kmer1)
                if pos*2 < len(strand): direction_evidence += 1
                if pos*2 > len(strand): direction_evidence -= 1
            if kmer2 in norm_strand:
                pos = norm_strand.index(kmer2)
                if pos*2 < len(strand): direction_evidence -= 1
                if pos*2 > len(strand): direction_evidence += 1
            seq = sequence_from_kmer_list(strand)
            if direction_evidence < 0:
                seq = bio.reverse_complement(seq)
            print '    ', seq, '>>>' if direction_evidence else ''
    
        if action == 'graph trace':
            n = 100
            f_set, f_switchbacks = trace_forward(self.parent, kmer1, n, self._equivalents)
            b_set, b_switchbacks = trace_forward(self.parent, bio.reverse_complement(kmer1), n, self._equivalents)
            
            for item in f_switchbacks + b_switchbacks:
                viewer.annotate(self._root_norm(item), 100.0, 1.0,1.0,0.0)            
            for item in f_set:
                viewer.annotate(self._root_norm(item), float(n-f_set[item])/n*25.0, 1.0,0.0,0.0)
            for item in b_set:
                viewer.annotate(self._root_norm(item), float(n-b_set[item])/n*25.0, 0.0,1.0,0.0)

        if action == 'read/pair connections':
            for bag in self.bags.values():
                locations = bag.kmer_to_sequences.get(kmer1,[])
                
                print bag.working_dir, len(locations)
                if not locations: continue
                
                rkmer = bio.reverse_complement(kmer1)
                weight = 50.0 / len(locations)
                for location in locations:
                    sequences = bag.get_sequences(location)
                    for sequence in sequences:
                        origin = kmer1 in sequence or rkmer in sequence
                        for i in xrange(len(sequence)-bag.k+1):
                            kmer2 = norm(sequence[i:i+bag.k])
                            if kmer2 in graph.name_to_ident:
                                if origin:
                                    viewer.annotate(kmer2, weight, 0.9,0.0,0.9)
                                else:
                                    viewer.annotate(kmer2, weight, 0.9,0.9,0.0)
        
    
                

def histogram_from_values(values):
    histogram = [ ]
    for n in values:
        while len(histogram) <= n: histogram.append(0)
        histogram[n] += 1
    return histogram
                    

def compare(bag1, bag2):
     print 'Scan 1'
     kmers1 = bag1.kmers_reaching_cutoff(10)
     print 'Scan 2'
     kmers2 = bag2.kmers_reaching_cutoff(10)
     
     set_kmers1 = set(kmers1)
     set_kmers2 = set(kmers2)
     
     print
     
     print len(set_kmers1), bag1.working_dir
     print len(set_kmers2), bag2.working_dir
     
     print
     
     print len(set_kmers1 | set_kmers2), 'total kmers'
     print len(set_kmers1 & set_kmers2), 'present in both'
     print len(set_kmers1 - set_kmers2), 'unique to', bag1.working_dir
     print len(set_kmers2 - set_kmers1), 'unique to', bag2.working_dir

     print
     all_kmers1 = set(bag1.kmers_reaching_cutoff(1))
     all_kmers2 = set(bag2.kmers_reaching_cutoff(1))
     print len(set_kmers1 - all_kmers2), 'good and unique to', bag1.working_dir
     print len(set_kmers2 - all_kmers1), 'good and unique to', bag2.working_dir
     
     
     kmers = set_kmers1 & set_kmers2
     
     import pylab
     pylab.subplot(1,2,1)
     
     print 'h1'
     h1 = histogram_from_values(kmers1.values())
     pylab.loglog(range(1,len(h1)), h1[1:], label=bag1.working_dir)
     
     print 'h2'
     h2 = histogram_from_values(kmers2.values())
     pylab.loglog(range(1,len(h2)), h2[1:], label=bag2.working_dir)
     pylab.legend()
     
     pylab.subplot(1,2,2)
     pylab.loglog(
         [ kmers1.get(kmer,0) for kmer in kmers ],
         [ kmers2.get(kmer,0) for kmer in kmers ],
         ','
     )
     pylab.xlabel(bag1.working_dir)
     pylab.ylabel(bag2.working_dir)
     pylab.show()

def coverage_plot(sequence_filename, bag_filenames):
    sequence = io.read_sequences(sequence_filename).next()[1]
    
    import pylab
    for bag_filename in bag_filenames:
        bag = Kmer_bag(bag_filename)
        
        kmers = kmers_from_sequence(bag.k, sequence)
        coverage = [ bag.counts[norm(item)] for item in kmers ]
        pylab.plot(coverage, label=bag_filename)
        
        bag.close()
    
    pylab.legend()
    pylab.show()

def pair_plot(sequence_filename, bag_filenames):
    import pylab, numpy

    sequence = io.read_sequences(sequence_filename).next()[1]
    
    for bag_filename in bag_filenames:
        bag = Kmer_bag(bag_filename)
        
        kmers = [ norm(kmer) for kmer in kmers_from_sequence(bag.k, sequence) ]
        kmer_pos = { }
        for i, kmer in enumerate(kmers):
            kmer_pos[kmer] = i
        
        counts = numpy.zeros((len(kmers),len(kmers)), 'float64')
        for i in xrange(len(kmers)):
            print i, 'of', len(kmers)
            loci = bag.kmer_to_sequences.get(kmers[i],[]) [::100]
            if not loci: continue
            weight = 1.0 #/len(loci)
            for location in loci:
                for read in bag.get_sequences(location)[1]:
                    for kmer in kmers_from_sequence(bag.k, read):
                        kmer = norm(kmer)
                        if kmer in kmer_pos:
                            j = kmer_pos[kmer]
                            #if abs(i-j) > 200:
                            counts[i,j] += weight
        
        pylab.imshow(counts, interpolation='nearest')
        pylab.show()


def chimera_plot(sequence_filename, bag_filenames):
    sequence = io.read_sequences(sequence_filename).next()[1]
    
    import pylab
    for bag_filename in bag_filenames:
        bag = Kmer_bag(bag_filename)

        kmers = [ kmer for kmer in kmers_from_sequence(bag.k, sequence) ]
        
        counts = numpy.zeros((len(kmers),len(kmers)), 'float64')
        
        split = bag.k//2 + 1
        #for split in xrange(5,bag.k-5):
        for i in xrange(len(kmers)):
            for j in xrange(len(kmers)):
                if abs(i-j) < 10: continue
                chimera = kmers[i][:split] + kmers[j][split:]
                #chimera = kmers[i][:split] + bio.reverse_complement(kmers[j])[split:]
                
                n = bag.counts[norm(chimera)]
                if n:
                    print i,j,n
                    
                    counts[i,j] += n
        
        pylab.imshow(counts, interpolation='nearest')
        pylab.show()



def clean(k, min_cutoff, max_cutoff, shadow_cutoff, filenames):
    """
        For each read:
            Index nkmers->read#,position in read (encode as read#<<32 + position)
       
        for each nkmer:
            If depth good,
                save (read#,pos)->''
       
        For each read,
            get good
            clip, output
            
        
        
        Note clipping chooses the longest good subrange.
        This may discard some kmers, pushing others below the cutoff in the output... don't panic about this!
        
    """
    
    def scan_reads():
        i = 0
        for filename in filenames:
            for name, seq in io.read_sequences(filename):
                yield i, name, seq
                i += 1
    
    def pack(integer):
        return struct.pack('>Q', integer)
        #return '%060d' % integer
    def unpack(string):
        return struct.unpack('>Q', string)[0]
        #return int(string)
    
    temp_prefix = 'temp_%s_' % os.getpid()
    
    has_shadow_cutoff = (shadow_cutoff is not None)
    has_max_cutoff = (max_cutoff is not None)
    
    kmer_to_readpos = treemaker.Int_list_store(temp_prefix + 'kmer')
    kmer_count_and_shadow = treemaker.Count_and_shadow_store(temp_prefix + 'count_and_shadow')
    good_kmers_in_reads = treemaker.General_store(temp_prefix + 'good')
    try:
        grace.status('Indexing kmers')
        for i, name, seq in scan_reads():
            if i % 10000 == 0:
                grace.status('Indexing kmers, read %s' % grace.pretty_number(i))
            for good_seq in BAD_BASES.split(seq):
                for j in xrange(len(good_seq)-k+1):
                    nkmer = norm(good_seq[j:j+k])
                    kmer_to_readpos.append_to(nkmer, [(i<<32) + j] )
                
                    kmer_count_and_shadow.add(nkmer, 1,0)
        
        if has_shadow_cutoff:
            grace.status('Casting shadows')
            i = 0
            for nkmer, (count, shadow) in kmer_count_and_shadow.iter_from(''):
                if i % 1000 == 0:
                    grace.status('Casting shadows, kmer %s' % nkmer)
            
                for j in xrange(k):
                    for char in 'ACGT':
                        if char != nkmer[j]:
                            shadowed_nkmer = norm( nkmer[:j] + char + nkmer[j+1:] )
                            kmer_count_and_shadow.add(shadowed_nkmer, 0, count)
                
                i += 1
        
        grace.status('Selecting good kmers')
        
        #for nkmer, values in kmer_to_readpos.iter_from(''):
        #    n = len(values)
        #    if n < min_cutoff or (max_cutoff is not None and n > max_cutoff): continue <-- note bug
        #    for value in values:
        #        good_kmers_in_reads[pack(value)] = ''

        for nkmer, (count, shadow) in kmer_count_and_shadow.iter_from(''):
            if count < min_cutoff or \
               (has_max_cutoff and count > max_cutoff) or \
               (has_shadow_cutoff and count < shadow_cutoff * shadow):
                continue
                
            for value in kmer_to_readpos[nkmer]:
                good_kmers_in_reads[pack(value)] = ''

    
        #Don't need these any more
        kmer_to_readpos.clear()
        kmer_count_and_shadow.clear()
    
        grace.status('Writing clipped reads')
        for i, name, seq in scan_reads():
            good_positions = [ ]
            for packed_value, empty in good_kmers_in_reads.iter_from(pack(i<<32)):
                value = unpack(packed_value)
                if (value>>32) != i: break
                good_positions.append(value & 0xffffffff)
            
            best_start = 0
            best_length = 0
            current_start = 0
            current_length = 0
            for position in good_positions:
                if position != current_start+current_length:
                    current_start = position
                    current_length = 0
                
                current_length += 1                
                if current_length > best_length:
                    best_start = current_start
                    best_length = current_length
            
            if best_length:
                end = best_start+best_length+k-1
                clipped = seq[best_start:end]
            else:
                clipped = ''
            io.write_fasta(sys.stdout, name, clipped)    
    finally:
        grace.status('Cleaning up')
        kmer_to_readpos.destroy()
        kmer_count_and_shadow.destroy()
        good_kmers_in_reads.destroy()
        grace.status('')

CLEAN_USAGE = """\

Usage:

    nesoni clean [options] [<reads.fa> ...] >cleaned_reads.fa

Options:

    --k   NNN     k-mer size, default 30
    
    --min NNN     Minimum number of times a k-mer must be seen, 
                  default: 1
    
    --max NNN     Maximum number of times a k-mer can be seen, 
                  default: no upper limit
                
    --shadow N.NN k-mer must be seen this proportion of the 
                  times k-mers one SNP away were seen. 
                  Expensive but powerful!
                  eg --shadow 0.5
                  default: not enabled
                 


Reads are clipped or removed so as to remove all k-mers seen too few
times or too many times, or seen much less frequently than SNP-neighbouring 
k-mers.

That is, good k-mers are identified in each read, and the longest contiguous
run of such k-mers is chosen.

"""

def clean_main(args):
    k, args = grace.get_option_value(args, '--k', int, 30)
    min_cutoff, args = grace.get_option_value(args, '--min', int, 1)
    max_cutoff, args = grace.get_option_value(args, '--max', int, None)
    shadow_cutoff, args = grace.get_option_value(args, '--shadow', float, None)
    filenames = args

    if len(filenames) < 1:
        print CLEAN_USAGE
        return 1    
    
    clean(k, min_cutoff, max_cutoff, shadow_cutoff, filenames)    

        

BAG_USAGE = """\

Usage:

    nesoni bag <name> new <k>
    
where <k> is the k-mer size will set up a new k-mer bag in a directory called <name>.


    nesoni bag <name> reads <reads.fa> [<reads.fa> ...]

adds a set of reads to a bag.


    nesoni bag <name> pairs <reads.fa> [<reads.fa> ...]
    
adds a set of read-pairs to a bag. Illumina naming convention is assumed.


    nesoni bag <name> kmers <reads.fa> [<reads.fa> ...]

will count kmers in a set of reads, but not index from kmers to reads.
This is faster and uses less space, but you won't be able to examine
read or pair mapping in a graph layout.  


    nesoni bag <name> stats [--plot yes]

displays some k-mer frequency statistics. If you specify --plot yes,
it will plot some stuff.


    nesoni bag <name> depth [--format text/jal] <sequences.fa>

displays depth statistics for the sequences given.


You can combine multiple commands.

Typical usage:

    nesoni bag foo new 25 pairs foo_reads_1.fa foo_reads_2.fa 

"""

def bag_main(args):
    if not args:
        print BAG_USAGE
        return 1

    kmer_bags = [ ]
    try:
        def add_bags(args):
            grace.expect_no_further_options(args)
            for working_dir in args:
                kmer_bags.append( Kmer_bag(working_dir) )
            for bag in kmer_bags:
                assert bag.k == kmer_bags[0].k, 'Inconsistent kmer size in bag collection'

        def new(args):
            need_exactly_one()
            k = int( args.pop(0) )
            kmer_bags[0].clear()
            kmer_bags[0].set_k(k)
        
        #Commands that follow assume k has been set
        def need_k():
            if kmer_bags[0].k is None:
                raise grace.Error('First set up the working directory with the "new" command')
        
        def need_exactly_one():
            assert len(kmer_bags) == 1, 'This command needs exactly one bag'
        
        def kmers(args):
            need_exactly_one()
            need_k()
            filenames = args
            kmer_bags[0].load_reads_count_only(filenames)
        
        def reads(args):
            need_exactly_one()
            need_k()
            filenames = args
            kmer_bags[0].load_reads(filenames)
        
        def pairs(args):
            need_exactly_one()
            need_k()
            filenames = args
            kmer_bags[0].load_pairs(filenames)
        
        def stats(args):
            need_exactly_one()
            need_k()
            show_plot, args = grace.get_option_value(args,'--plot',grace.as_bool,False)
            kmer_bags[0].stats(show_plot)
        
        def depth(args):
            need_k()
            format, args = grace.get_option_value(args,'--format',str,'text')
            depth_stats(kmer_bags, args, format)
                
        grace.execute(args, [
            new, kmers, reads, pairs, stats, depth
        ], default_command = add_bags)
    finally:
        for kmer_bag in kmer_bags:
            kmer_bag.close()    



GRAPH_USAGE = """\

Usage:

    nesoni graph <name> new
    
will set up a new graph directory called <name>.


    nesoni graph <name> with <bag_name> [<bag_name>...]
    
will add a kmer-bag to the graph.


    nesoni graph <name> build [options] <cutoff> [<seed_sequences.fa> ...]
    
    Options:
        --min-size NNN        - minimum object size
                                default: no minimum
        --merge-snps yes/no   - merge kmers differing by a SNP
                                default: no
        --merge-indels yes/no - merge kmers differing by an indel
                                default: no

construct and lay out a deBruijn graph of kmers occurring at 
least <cutoff> times in the bags associated with this graph.

If <minimum connected size> is given, only connected objects containing
this many k-mers or more will be used.

If <seed_sequences.fa> is given, only connected objects containing at least
one k-mer from one of the seed sequences will be used.


    nesoni graph <name> show [<overlay_sequences.fa> ...]

will display the deBruijn graph and allow you to interact with it.
You may do this while the graph is still "build"ing. If you specify
a file of overlay sequences, these will be highlighted in the graph
and annotated with the sequence names.


You can combine multiple commands (unless the number of parameters
is variable).

Typical usage:

    nesoni graph bar new with foo baz build 10 

    nesoni graph bar show

"""

def graph_main(args):
    if not args:
        print GRAPH_USAGE
        return 1
        
    args = args[:]
    working_dir = args.pop(0)
    
    kmer_graph = Kmer_graph(working_dir)
    try:
        def do_new(args):
            kmer_graph.clear()
        
        def do_with(args):
            for path in args:
                kmer_graph.add_bag(path)
        
        def do_build(args):
            minimum_object_size, args = grace.get_option_value(args,'--min-size',int,None)
            merge_snps, args = grace.get_option_value(args,'--merge-snps',grace.as_bool,False)
            merge_indels, args = grace.get_option_value(args,'--merge-indels',grace.as_bool,False)
            
            cutoff = int(args.pop(0))
            
            #minimum_object_size = None
            #if len(args) and args[0].isdigit():
            #    minimum_object_size = int(args.pop(0))
            
            seed_filenames = args
            kmer_graph.make_graph(cutoff, seed_filenames, minimum_object_size, merge_snps, merge_indels)
            kmer_graph.graph_layout()

        def do_layout(args):
            kmer_graph.graph_layout()
        
        def do_show(args):
            bag_to_compare, args = grace.get_option_value(args, '--compare-bag', str, None) 
            overlay_filenames = args
            kmer_graph.graph_show(overlay_filenames, bag_to_compare)
                
        grace.execute(args, {
            'new' : do_new,
            'with' : do_with,
            'build' : do_build,
            'layout' : do_layout,
            'show' : do_show,
        })
        
        #while args:
        #    command = args.pop(0)
        #    
        #    if command == 'new':
        #        kmer_graph.clear()
        #        continue
        #
        #    if command == 'with':
        #        path = args.pop(0)
        #        kmer_graph.add_bag(path)
        #        continue
        #    
        #    if command == 'build':
        #        cutoff = int(args.pop(0))
        #        
        #        minimum_object_size = None
        #        if len(args) and args[0].isdigit():
        #            minimum_object_size = int(args.pop(0))
        #        
        #        seed_filenames = args
        #        args = []
        #        kmer_graph.make_graph(cutoff, seed_filenames, minimum_object_size)
        #        kmer_graph.graph_layout()
        #        continue
        #    
        #    if command == 'layout':
        #        kmer_graph.graph_layout()
        #        continue
        #    
        #    if command == 'show':
        #        overlay_filenames = args
        #        args = []
        #        kmer_graph.graph_show(overlay_filenames)
        #        continue
        #        
        #    print GRAPH_USAGE
        #    return 1
        
    finally:
        kmer_graph.close()


#def main(args):
#    args = args[:]
#    what = args.pop(0)
#    
#    if what == 'bag':
#        bag_main(args)
#    elif what == 'graph':
#        graph_main(args)
#    elif what == 'coverage':
#        coverage_plot(args[0], args[1:])
#    elif what == 'pair':
#        pair_plot(args[0], args[1:])
#    elif what == 'chimera':
#        chimera_plot(args[0], args[1:])
#    else:
#        assert False


