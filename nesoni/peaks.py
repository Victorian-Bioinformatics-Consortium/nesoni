"""

Use all top scoring alignments

Trim start and end by n bases
Call regions above a threshold
Discard if less than min size


"""

from __future__ import division

import os, collections, math

import nesoni
from nesoni import config, legion, io, annotation, sam, grace


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
            if self.step[self.queue[i]] >= max(min_step,(self.queue[i][1]-self.queue[i][0])//4):
                self.elaborate(self.queue[i][0],self.queue[i][1], True)
            else:
                i += 1
        self.queue.sort(key=lambda item: item[1]-item[0])
            
    def elaborate(self, start, end, once=False):
        while True:
            step = self.step[(start,end)]
            if not step: break
            all = True
            for new_start, new_end in [
                     (start-step,end+step),
                     (start,end+step),
                     (start-step,end),
                     (start+step,end+step),
                     (start-step,end-step),
                     (start+step,end),
                     (start,end-step),
                     (start+step,end-step),
                     ]:
                if new_start >= new_end or new_end > self.size or new_start < 0 or \
                   (new_start,new_end) in self.step: 
                    all = False
                    continue
                self.queue.append((new_start,new_end))
                self.step[(new_start,new_end)] = step
            self.step[(start,end)] = step//2
            #if once or all: break
            if once: break
    
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


def find_peaks(depth, power, moderation=5.0, width_power=1.0, min_initial_size=8):
    """ Find peaks from a depth of coverage profile (list of integers). 
        """
    #integral = [ 0.0 ]
    #for i in xrange(len(depth)):
    #    integral.append(integral[i]+(depth[i]+moderation)**-power)
    
    pow_depth = [ max(item,moderation)**-power for item in depth ]

    dominance = [ 0.0 ] * len(depth)            
    candidates = [ ]
    
    exp = Explorer(len(depth), min_initial_size)
    for start,end in exp:
        width = end-start
        #total = integral[end]-integral[start]
        total = sum( pow_depth[i] for i in xrange(start,end) )
        score = math.log(width)*width_power - math.log(total)
        
        all_good = True
        #n = 0
        for i in xrange(start,end):
            if score > dominance[i]:
                dominance[i] = score
                #n += 1
            else:
                all_good = False
        #if n*16 > end-start:
        if all_good:
            exp.elaborate(start,end)
            candidates.append((start,end,score))
    
    peaks = [ ]
    for start,end,score in candidates:
        if ( all( dominance[i] <= score for i in xrange(start,end) ) and
             (power < 0 or any( depth[i] > moderation for i in xrange(start,end) )) ):
            peaks.append((start,end))

    peaks.sort()
    return peaks



@config.String_flag(
    'filter',
    'Filtering mode:\n'
    'poly     - Use all top scoring alignments\n'
    'mono     - Use top scoring alignment if unique\n'
    'existing - Use alignments from "filter:" or "consensus:"\n'
    )
@config.Bool_flag(
    'deduplicate',
    'Count fragment alignments with the exact same start and end positions as a single alignment.'
    )
@config.Bool_flag(
    'strand_specific',
    'Are the reads strand specific?'
    )
@config.String_flag(
    'what',
    'fragment - Peaks based on entire fragment.\n'
    '3prime - Peaks based on last base of fragment.\n'
    '5prime - Peaks based on first base of fragment.\n'
    )
@config.Int_flag(
    'lap',
    'Add this many bases to the end of each fragment when calculating depth, '
    'then subtract this many bases from called peaks, back to original size. '
    'Positive values will tend to join up small gaps. '
    'Negative values may allow calling of slightly overlapping transcripts, '
    'and enhance calling of transcripts that are close together. '
    )
@config.Float_flag(
    'crosstalk',
    'Before calling peaks, subtract the depth of the reverse strand depth scaled by this amount.'
    )
@config.String_flag(
    'type',
    'Type of feature to produce in GFF output.'
    )
@config.Main_section(
    'filenames',
    'Working directories or BAM files (sorted by read name).'
    )
class Span_finder(config.Action_with_prefix):
    filter = 'poly'
    deduplicate = False
    strand_specific = True
    what = 'fragment'
    lap = 0
    crosstalk = 0.0
    type = 'peak'
    filenames = [ ]
    
    def run(self):
        assert self.what in ('fragment','5prime','3prime'), 'Unknown option for --what.'
        #assert self.moderation > 0.0, '--moderation must be greater than zero.'
        #assert self.power > 0.0, '--power must be greater than zero.'
        #assert self.width_power >= 1.0, '--width-power must be greater than or equal to one.'
    
        #if self.filter == 'poly':
        #    use_bam_filename = 'alignments.bam'
        #    use_only_top = True
        #    use_only_monogamous = False
        #    expect_multiple_alignments = True
        #elif self.filter == 'mono': 
        #    use_bam_filename = 'alignments.bam'
        #    use_only_top = True
        #    use_only_monogamous = True
        #    expect_multiple_alignments = True
        #else:
        #    assert self.filter == 'existing', 'Unrecognized filtering mode'
        #    use_bam_filename = 'alignments_filtered.bam'
        #    use_only_top = False
        #    use_only_monogamous = False
        #    expect_multiple_alignments = False
                    
        spans = collections.defaultdict(list)
        
        for item in legion.parallel_imap(self._load_bam, self.filenames):
            for key,value in item.items():
                spans[key].extend(value)
        
        #for i, filename in enumerate(self.filenames):
        #    if os.path.isdir(filename):
        #        filename = os.path.join(filename, use_bam_filename)
        #    
        #    n = 0
        #    for read_name, fragment_alignments, unmapped in \
        #            sam.bam_iter_fragments(
        #                filename, 
        #                'Scanning sample %d of %d' % (i+1,len(self.filenames))):
        #        if not fragment_alignments:
        #            continue
        #            
        #        if use_only_top:
        #            fragment_scores = [ sum( al.get_AS() for al in item ) for item in fragment_alignments ]            
        #            best_score = max(fragment_scores)
        #            fragment_alignments = [ 
        #                item 
        #                for item, score in zip(fragment_alignments, fragment_scores)
        #                if score >= best_score ]            
        #        
        #        for alignments in fragment_alignments:
        #            if self.strand_specific:
        #                strand = -1 if alignments[0].flag&sam.FLAG_REVERSE else 1
        #            else:
        #                strand = 0
        #        
        #            start = min(item.pos-1 for item in alignments)
        #            end = max(item.pos+item.length-1 for item in alignments)
        #            if end-start <= self.trim*2: continue
        #            
        #            rname = alignments[0].rname                    
        #            spans[(rname, strand)].append((start+self.trim,end-self.trim))
        #        
        #        n += 1
        #        #if n > 100000: break
        #
        #if self.deduplicate:
        #    for key in spans:
        #        spans[key] = list(set(spans[key]))

        grace.status('Calling peaks')

        f = open(self.prefix+'.gff', 'wb')
        annotation.write_gff3_header(f)
        
        n = 0

        for (rname, strand), span_list in spans.items():
            depth = [ 0.0 ] * (1+max( item[1] for item in span_list ))
            for start, end in span_list:
                depth[start] += 1.0
                depth[end] -= 1.0
            
            if self.crosstalk and strand and (rname,-strand) in spans:
                for start, end in spans[(rname,-strand)]:
                    if start < len(depth): depth[start] -= self.crosstalk
                    if end < len(depth): depth[end] += self.crosstalk
            
            for i in xrange(1,len(depth)):
                depth[i] += depth[i-1]

            if self.crosstalk:
                for i in xrange(len(depth)):
                    depth[i] = max(0.0,depth[i])

            #import pylab
            #pylab.plot(depth)
            
            for start, end in self._find_spans(depth):
                #pylab.axvspan(start-0.5,end-0.5,alpha=0.25)
                
                if end-self.lap-start <= 0: continue
                
                n += 1
                
                id = 'peak%d' % n
                
                #if strand == -1:
                #    id = '%s-%d..%d' % (rname,start,end+1)
                #elif strand == 0:
                #    id = '%s.%d..%d' % (rname,start+1,end)
                #else:
                #    id = '%s+%d..%d' % (rname,start+1,end)
                
                ann = annotation.Annotation()
                ann.source = 'nesoni'
                ann.type = self.type
                ann.seqid = rname
                ann.start = start
                ann.end = end - self.lap
                ann.strand = strand
                ann.score = None
                ann.phase = None
                ann.attr = { 
                    'id' : id,
                    'color' : '#00ff00' if strand > 0 else '#0000ff' if strand < 0 else '#008080',
                    }
                print >> f, ann.as_gff()
            f.flush()
            
            #pylab.show()

        f.close()
        
        self.log.datum('-','called peaks',n)
        
        grace.status('')


    def _load_bam(self, filename):
        spans = { }

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
        
        if os.path.isdir(filename):
            filename = os.path.join(filename, use_bam_filename)
        
        for read_name, fragment_alignments, unmapped in \
                sam.bam_iter_fragments(
                    filename, 
                    'Scanning'):
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
                
                if self.what == '5prime':
                   if strand >= 0:
                       end = start+1
                   else:
                       start = end-1
                elif self.what == '3prime':
                   if strand >= 0:
                       start = end-1
                   else:
                       end = start+1
                
                if end+self.lap-start <= 0: continue
                
                rname = alignments[0].rname
                if (rname,strand) not in spans: 
                    spans[(rname,strand)] = [ ]          
                spans[(rname, strand)].append((start,end+self.lap))
                
        if self.deduplicate:
            for key in spans:
                spans[key] = list(set(spans[key]))
        
        return spans


@config.help(
    'Call peaks that are higher than their surrounding coverage.\n\n'    
    'Peaks are often surrounded by a low level of noise.'
    ' The peak caller seeks to ignore this noise,'
    ' in a manner somewhat similar to human auditory masking.'
    ' This is controlled by the --power parameter,'
    ' and is independant of scale and depth of coverage.'
    '\n\n'    
    'Algorithm details: '
    'The peak caller is intended to be invariant under rescaling of the size and depth of peaks. '
    'The basic idea is to score potential peaks using a the area of a rectangle with '
    'width equal to the width of the span, and height equal to the harmonic mean of some function of the depths within the span. '
    ' A peak is reported if there is no higher scoring potential peak that overlaps it.'
    ' (A heuristic method is used to choose potential peaks to examine, the search is moderately thorough but'
    ' not guaranteed to be exhaustive.) '
    'The actual scoring function is: \n\n'
    '  width^width_power / sum( (depth+moderation)^-power )'    
    )
@config.Float_flag(
    'moderation',
    'moderation > 0\n'
    'The water line. Depth is clipped to be at least this much. '
    'Peaks will not be called with depth less than this.'
    )
@config.Float_flag(
    'power',
    'power > 0\n'
    'Smaller values will encourage more and shorter peaks. '
    'Larger values will encourage less and longer peaks. '
    )
@config.Float_flag(
    'width_power',
    'width-power >= 1\n'
    'This probably isn\'t an important parameter to tweak. '
    'Larger values will encourage calling wider peaks. '
    'Smaller values will encourage calling shorter peaks.'
    )
class Islands(Span_finder):
    #min_depth = 25
    moderation = 5.0
    power = 10.0
    width_power = 2.0

    def _find_spans(self, depth):            
        return find_peaks(depth, power=self.power, moderation=self.moderation, width_power=self.width_power)




@config.help(
    'Call transcripts.\n'
    'Potential transcripts are chosen as runs having at least a given minimum depth, '
    'then filtered on having at least a given median depth.'    
    )
@config.Int_flag(
    'min_depth',
    'All of the transcript must have at least this depth.'
    )
@config.Int_flag(
    'median_depth',
    'Potential transcripts identified using --min-depth must also have at least this median depth.'
    )    
class Transcripts(Span_finder):
    min_depth = 1
    median_depth = 100

    def _find_spans(self, depth):
        result = [ ]
        def consider(start,end):
            if start < end:
                items = depth[start:end]
                items.sort()
                if items[len(items)//2] >= self.median_depth:
                    result.append((start,end))
            
        start = 0
        for i, item in enumerate(depth):
            if item < self.min_depth:
                consider(start,i)
                start = i+1                
        consider(start,len(depth))
        
        return result

@config.help(
    'Call peaks that are higher than everything within "radius", '
    'and higher than some minimum depth.\n'
    '\n'
    'Note: The --lap parameter can be used to apply some smoothing, '
    'if you are looking for modes that have a bit of width as well as pure height.\n'
    )
@config.Int_flag(
    'min_depth',
    'Minimum depth.'
    )
@config.Int_flag(
    'radius',
    'Smaller mode suppression radius.'
    )
class Modes(Span_finder):
    min_depth = 5
    radius = 20

    def _find_spans(self, depth):
        result = [ ]
        
        i = 0
        while i < len(depth):
            j = i+1
            while j < len(depth) and depth[j] == depth[i]:
                j += 1
            
            lap = max(0,self.lap-(j-i)+1)
            lap_back = lap//2
            lap_forward = lap-lap_back

            if depth[i] >= self.min_depth:
                for k in xrange(max(0,i-self.radius),min(len(depth),j+self.radius)):
                    if depth[k] > depth[i] or (depth[k] == depth[i] and k < i): #Resolve ties arbitrarily
                        break
                else:
                    result.append((i-lap_back,j+lap_forward))                
            i = j
            
        return result


if __name__ == '__main__':
    nesoni.run_tool(Peaks)            

            
    

