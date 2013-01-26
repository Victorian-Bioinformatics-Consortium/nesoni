"""

Use all top scoring alignments

Trim start and end by n bases
Call regions above a threshold
Discard if less than min size


"""

from __future__ import division

import os, collections, math

import nesoni
from nesoni import config, io, annotation, sam, grace


class Explorer(object):
    """ Generate regions of various sizes and positions.
    
        Initially, explore a set of regions following a
        roughly wavelet-like decompositions: many small regions,
        fewer large regions.
    
        If a region is found to be interesting, investigate
        similar regions (nearby and of similar size).
        """
    
    def __init__(self, size, min_step):
        self.size = size
        self.queue = [ (0,self.size) ]
        self.queue_pos = 0
        
        initial_step = 1
        while initial_step*2 < size: initial_step *= 2
        self.step = { (0,self.size) : initial_step }
        
        i = 0
        while i < len(self.queue):
            if self.step[self.queue[i]] >= max(min_step,(self.queue[i][1]-self.queue[i][0])//8):
                self.elaborate(self.queue[i][0],self.queue[i][1], True)
            else:
                i += 1
            
    def elaborate(self, start, end, once=False):
        while True:
            step = self.step[(start,end)]
            if not step: break
            all = True
            for new_start, new_end in [
                     (start+step,end),
                     (start-step,end),
                     (start,end+step),
                     (start,end-step),
                     (start+step,end+step),
                     (start+step,end-step),
                     (start-step,end+step),
                     (start-step,end-step),
                     ]:
                if new_start >= new_end or new_end > self.size or new_start < 0 or \
                   (new_start,new_end) in self.step: 
                    all = False
                    continue
                self.queue.append((new_start,new_end))
                self.step[(new_start,new_end)] = step
            self.step[(start,end)] = step//2
            if once or all: break
            #if once: break
    
    def __iter__(self):
        return self
    
    def next(self):
        if self.queue_pos >= len(self.queue):
            raise StopIteration()
        result = self.queue[self.queue_pos]
        self.queue_pos += 1
        return result

#e = Explorer(65536, 8)
#print e.queue
#import pylab
#pylab.plot([x for x,y in e.queue],[y for x,y in e.queue],'.')
#pylab.show()
#foo

def find_peaks(depth, min_depth, power, depth_power=1.0, min_initial_size=8):
    """ Find peaks from a depth of coverage profile (list of integers). 
        """
    integral = [ 0.0 ]
    for i in xrange(len(depth)):
        integral.append(integral[i]+depth[i]**(power * depth_power))

    dominance = [ 0.0 ] * len(depth)            
    candidates = [ ]
    
    exp = Explorer(len(depth), min_initial_size)
    for start,end in exp:
        width = end-start
        mean = (integral[end]-integral[start]) / width
        #mean = sum(pdepth[start:end]) / width
        if mean <= 0.0: continue
        
        score = math.log(mean) + math.log(width)*power
        
        all = True
        n = 0
        for i in xrange(start,end):
            if score > dominance[i]:
                dominance[i] = score
                n += 1
            else:
                all = False
        if n*16 > end-start:
            exp.elaborate(start,end)
        if all:
            if max(depth[start:end]) >= min_depth:
                candidates.append((start,end,score))
    
    peaks = [ ]
    for start,end,score in candidates:
        if score >= max(dominance[start:end]):
            peaks.append((start,end))
    peaks.sort()
    return peaks



@config.help(
    'Call peaks that are higher than their surrounding coverage.\n\n'
    
    'Peaks are often surrounded by a low level of noise.'
    ' The peak caller seeks to ignore this noise,'
    ' in a manner somewhat similar to human auditory masking.'
    ' This is controlled by the --power parameter,'
    ' and is independant of scale and depth of coverage.'
    '\n\n'
    
    'Algorithm details:'
    ' Potential peaks are scored as the mean of the depths aveaged using a power mean'
    ' multiplied by the peak width.'
    ' A peak is reported if there is no higher scoring potential peak that overlaps it.'
    ' (A heuristic method is used to choose potential peaks to examine, the search is moderately thorough but'
    ' not guaranteed to be exhaustive.)\n\n'
    
    'Note: Not intended for calling RNA transcripts.'
    )
@config.Float_flag(
    'power',
    '0 < power < 1\n'
    'Larger values will suppress calling of peaks near taller peaks.'
    )
@config.Float_flag(
    'depth_power',
    'Smaller values will allow more variability within a peak.'
    )
@config.Int_flag(
    'min_depth',
    'Peaks must have at least this average depth of coverage.'
    )
@config.String_flag(
    'filter',
    'Filtering mode:\n'
    'poly     - Use all top scoring alignments\n'
    'mono     - Use top scoring alignment if unique\n'
    'existing - Use alignments from "filter:" or "consensus:"\n'
    )
@config.Bool_flag(
    'strand_specific',
    'Are the reads strand specific?'
    )
@config.String_flag(
    'type',
    'Type of feature to produce in GFF output.'
    )
@config.Main_section(
    'filenames',
    'Working directories or BAM files (sorted by read name).'
    )
class Peaks(config.Action_with_prefix):
    min_depth = 25
    power = 0.1
    depth_power = 0.5
    trim = 0
    filter = 'poly'
    strand_specific = True
    type = 'peak'
    filenames = [ ]
    
    def run(self):
        if self.filter == 'poly':
            use_bam_filename = 'alignments.bam'
            use_only_top = True
            use_only_monogamous = False
            expect_multiple_alignments = True
        elif self.filter == 'mono': 
            use_bam_filename = 'alignments.bam'
            use_only_top = True
            use_only_monogamous = True
            expect_multiple_alignments = True
        else:
            assert self.filter == 'existing', 'Unrecognized filtering mode'
            use_bam_filename = 'alignments_filtered.bam'
            use_only_top = False
            use_only_monogamous = False
            expect_multiple_alignments = False
                    
        spans = collections.defaultdict(list)
        
        for i, filename in enumerate(self.filenames):
            if os.path.isdir(filename):
                filename = os.path.join(filename, use_bam_filename)
            
            n = 0
            for read_name, fragment_alignments, unmapped in \
                    sam.bam_iter_fragments(
                        filename, 
                        'Scanning sample %d of %d' % (i+1,len(self.filenames))):
                if not fragment_alignments:
                    continue
                    
                if use_only_top:
                    fragment_scores = [ sum( al.get_AS() for al in item ) for item in fragment_alignments ]            
                    best_score = max(fragment_scores)
                    fragment_alignments = [ 
                        item 
                        for item, score in zip(fragment_alignments, fragment_scores)
                        if score >= best_score ]            
                
                for alignments in fragment_alignments:
                    if self.strand_specific:
                        strand = -1 if alignments[0].flag&sam.FLAG_REVERSE else 1
                    else:
                        strand = 0
                
                    start = min(item.pos-1 for item in alignments)
                    end = max(item.pos+item.length-1 for item in alignments)
                    if end-start <= self.trim*2: continue
                    
                    rname = alignments[0].rname                    
                    spans[(rname, strand)].append((start+self.trim,end-self.trim))
                
                n += 1
                #if n > 100000: break

        grace.status('Calling peaks')

        f = open(self.prefix+'.gff', 'wb')
        annotation.write_gff3_header(f)

        for (rname, strand), span_list in spans.items():
            depth = [ 0 ] * (1+max( item[1] for item in span_list ))
            for start, end in span_list:
                depth[start] += 1
                depth[end] -= 1
            for i in xrange(1,len(depth)):
                depth[i] += depth[i-1]

            #import pylab
            #pylab.plot(depth)
            
            n = 0
            for start, end in find_peaks(depth, self.min_depth, self.power, self.depth_power):
                #pylab.axvspan(start-0.5,end-0.5,alpha=0.25)
                
                n += 1
                
                if strand == -1:
                    id = '%s-%d..%d' % (rname,start,end+1)
                elif strand == 0:
                    id = '%s.%d..%d' % (rname,start,end+1)
                else:
                    id = '%s+%d..%d' % (rname,start,end+1)
                
                ann = annotation.Annotation()
                ann.source = 'nesoni'
                ann.type = self.type
                ann.seqid = rname
                ann.start = start - self.trim
                ann.end = end + self.trim
                ann.strand = strand
                ann.score = None
                ann.phase = None
                ann.attr = { 'id' : id }
                print >> f, ann.as_gff()
            f.flush()
            
            #pylab.show()

        f.close()
        
        grace.status('')

            
            
if __name__ == '__main__':
    nesoni.run_tool(Peaks)            

            
    

