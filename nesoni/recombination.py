
import sys, os, itertools, copy, numpy

from nesoni import io, bio, shrimp, grace

def pad_slice(string, start, end):
    if start < 0:
        return ' '*(-start) + string[:end]
    return string[start:end]

class Alignment:
   """
      Properties
      
      target_name
      target_seq
      target_start
      target_end
      
      query_name
      query_seq
      query_forward
      query_start
      query_end
   
   """
   
   def reversed(self): # Reverse query sequence
       result = copy.copy(self)
       result.query_seq = bio.reverse_complement(result.query_seq)
       result.query_start, result.query_end = \
           len(self.query_seq)-result.query_end, len(self.query_seq)-result.query_start
       result.query_forward = not result.query_forward
       return result
   
   #def flipped(self): # Exchange query and target
   #    result = Alignment()
   #    result.target_name = self.query_name
   #    result.target_seq = self.query_seq
   #    result.target_start = self.query_start
   #    result.target_end = self.query_end
   #    result.query_name = self.target_name
   #    result.query_seq = self.target_seq
   #    result.q
   #    return result

def alignment_from_shrimp(line, targets, read_name, read_seq):
    parts = line.rstrip().split()
    result = Alignment()
    result.query_name = parts[0][1:]
    result.target_name = parts[1]
    result.target_seq = targets[result.target_name]
    result.query_forward = (parts[2] == '+')
    result.target_start = int(parts[3])-1
    result.target_end = int(parts[4])
    result.query_start = int(parts[5])-1
    result.query_end = int(parts[6])
    result.query_seq = read_seq
    assert read_name == result.query_name
    assert len(read_seq) == int(parts[7])
    return result



USAGE = """\

Usage:

    nesoni recombination working-dir sequence-name >output.txt

Uses partial read alignments to find points of divergence from 
the reference sequence.

"""

def recombination(args):
    grace.expect_no_further_options(args)
    if len(args) != 2:
        print >> sys.stderr, USAGE
        raise grace.Help_shown()

    working_dir, seq_name = args

    references = dict(io.read_sequences(os.path.join(working_dir, 'reference.fa')))
    
    depth = { }
    prefixes = { }
    suffixes = { }
    for name in references:
        depth[name] = numpy.zeros(len(references[name]), 'int64')
        prefixes[name] = [ [] for base in references[name] ]
        suffixes[name] = [ [] for base in references[name] ]
    def register_divergence(hit):
        if not hit.query_forward:
            hit = hit.reversed()
        
        margin = 20
        
        if hit.target_end - hit.target_start < 20: 
            return False
        
        depth[hit.target_name][hit.target_start : hit.target_end] += 1

        any = False
        
        if hit.query_end <= len(hit.query_seq)-margin: # and hit.target_end < len(hit.target_seq):
            suffixes[hit.target_name][hit.target_end-1].append( hit.query_seq[hit.query_end:] )
            any = True
        
        if hit.query_start >= margin: # and hit.target_start > 0:
            prefixes[hit.target_name][hit.target_start].append( hit.query_seq[:hit.query_start] )
            any = True
        
        return any

    n = 0
    for (read_name, read_seq), hits in shrimp.iter_read_hits(working_dir):
        # Skip reads containing Ns
        if 'N' in read_seq: continue
    
        for line in hits:
            register_divergence(alignment_from_shrimp(line, references, read_name, read_seq))
        
        n += 1
        #if n > 100000:
        #    break
            
        if n%10000 == 0:
            grace.status('Processing read %s' % grace.pretty_number(n))

    grace.status('')

    
    def show_items(items):
        original_length = len(items)
        cut = 0
        while len(items) > 80:
            cut += 1
            items = [ item for item in items if item[0] >= cut ]
        for item in items:
            print item[1]
        if len(items) < original_length:
            print '(and %d more occurring %d times or less)' % (original_length-len(items), cut-1) 
    
    def score(items):
        if not items: return 1.0
        return float(sum( item[0] * item[0] for item in items )) / (sum( item[0] for item in items )**2)
    
    def summarize_prefixes(seqs, pad):
        seqs = sorted(seqs, key=lambda seq: seq[::-1])
        
        cut = 100
        while True:
            items = [ ] 
            for (seq, iterator) in itertools.groupby(seqs, key = lambda x: x[-cut:]):
                ss = list(iterator)
                anylong = any( item != seq for item in ss )            
                n = len(ss)
                items.append( (n, ('%'+str(pad)+'s')%(('...' if anylong else '') + seq) + ' x %d' % n) )
            
            if score(items) >= 1.0/20: break
            cut -= 1
        
        show_items(items)
        
    def summarize_suffixes(seqs, pad):
        seqs = sorted(seqs)
        
        cut = 100
        while True:
            items = [ ] 
            for (seq, iterator) in itertools.groupby(seqs, key = lambda x: x[:cut]):
                ss = list(iterator)            
                anylong = any( item != seq for item in ss )            
                n = len(ss)
                items.append( (n, ('%'+str(pad)+'s')%('%d x '%n) + seq + ('...' if anylong else '')) )
            
            if score(items) >= 1.0/20: break
            cut -= 1
        show_items(items)
    
    print 'Position        Depth        Changed prefixes             Changed suffixes'
    print '                          Count    % of depth       Count    % of depth'
    for i in xrange(len(references[seq_name])):
        print '%8d   %10d %9d  %11s       %9d  %11s' % (
            i+1,
            depth[seq_name][i],
            len(prefixes[seq_name][i]),
            '%.3f%%' % (len(prefixes[seq_name][i])*100.0/depth[seq_name][i]) if prefixes[seq_name][i] else '',  
            len(suffixes[seq_name][i]),
            '%.3f%%' % (len(suffixes[seq_name][i])*100.0/depth[seq_name][i]) if suffixes[seq_name][i] else '')
        #summarize_suffixes(suffixes[name][i], references[name][i+1:], references[name], suffix_depth[name][i])

    print
    print 'Details'
    print
    for i in xrange(len(references[seq_name])):
        print '%-80s*' % ('Base %d' % (i+1))
        
        print pad_slice(references[seq_name], i-80,i+1+80)        
        summarize_prefixes(prefixes[seq_name][i], 80)
        summarize_suffixes(suffixes[seq_name][i], 81)
        print


