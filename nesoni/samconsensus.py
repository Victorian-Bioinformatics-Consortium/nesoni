
import nesoni
from nesoni import io, bio, grace, config, sam, samshrimp, consensus, working_directory

import os, sys, subprocess, array, datetime, socket, bisect, random, itertools

FILTER_HELP = """\

Usage:

    nesoni samfilter: working_dir [options]

"""

class Depth(object):
    def __init__(self, size):
        self.size = size
        self.starts = { }
        self.ends = { }
    
    def increment(self, start, end):
        self.starts[start] = self.starts.get(start,0)+1
        self.ends[end] = self.ends.get(end,0)+1
    
    def __len__(self):
        return self.size
    
    def __iter__(self):
        delta = [ 0 ] * (self.size+1)
        for key, value in self.starts.iteritems():
            delta[key] += value
        for key, value in self.ends.iteritems():
            delta[key] -= value
        depth = 0
        for i in xrange(self.size):
            depth += delta[i]
            yield depth
        
        #value = 0
        #starts = self.starts
        #ends = self.ends
        #for i in xrange(self.size):
        #    value += starts.get(i,0)
        #    value -= ends.get(i,0)
        #    yield value

    def iter_starts(self):
        """ How many start at this base position """
        result = [ 0 ] * self.size
        for key, value in self.starts.iteritems():
            result[key] += value
        for i in xrange(self.size):
            yield result[i]
            
        #starts = self.starts
        #for i in xrange(self.size):
        #    yield starts.get(i,0)
    
    def iter_ends(self):
        """ How many end at this base position """
        result = [ 0 ] * (self.size+1)
        for key, value in self.ends.iteritems():
            result[key] += value
        for i in xrange(self.size):
            yield result[i]
            
        #ends = self.ends
        #for i in xrange(self.size):
        #    yield ends.get(i+1,0)

    def total(self):
        total = 0
        for pos, value in self.starts.iteritems():
            total += value*(self.size-pos)
        for pos, value in self.ends.iteritems():
            total -= value*(self.size-pos)
        return total


class Ref_depth(object): pass

def write_userplot(filename, array):
    f = open(filename,'wb')
    for x in array:
        f.write( '%d\n' % int(x) )
    f.close()

def write_userplots(prefix, arrays, output_strand_specific_depths):
    f = open(prefix + '.userplot', 'wb')

    if not output_strand_specific_depths:
        for a,b in itertools.izip(arrays[0],arrays[1]):
            f.write( '%d\n' % (a+b) )
    else:    
        f.write('# BASE forward reverse\n')
        f.write('# colour 0:196:0 196:0:0\n')
        for i, (a,b) in enumerate(itertools.izip(arrays[0],arrays[1])):
            f.write('%d %d %d\n' % (i+1,a,b))
    f.close()




@config.help("""\
Filter BAM file in directory created by "shrimp:" or "import:", and calculate depth of coverage.
""")
@config.Float_flag('infidelity', 
    'Any runner-up alignments scoring this close to the best alignment\'s score are considered possibly valid alignments. '
    'Unit is number of SNPs, allowance for indels depends on aligner\'s scoring scheme relative to SNPs.')
@config.Bool_flag('monogamous', 'Discard fragments with more than one possibly valid alignment.')
@config.Bool_flag('userplots', 'Output Artemis userplots.')
@config.Bool_flag('strand_specific', 'Output strand-specific Artemis userplots.')
@config.Bool_flag('random', 
    'Instead of using all ambiguous alignments, pick one at random from the equal top alignments. '
    'Affects depth of coverage plots, and base calling if --monogamous is no.')
class Filter(config.Action_with_working_dir):
    infidelity = 2.0
    monogamous = True
    userplots = True
    strand_specific = False
    random = False
    
    _workspace_class = working_directory.Working
    
    def run(self, log=None):
        grace.require_samtools()
        filter(self.working_dir, self.infidelity, self.monogamous, self.userplots, self.strand_specific, self.random, self.log)    



#def filter_parse_args(args, log):
#    #Parse arguments
#
#    log.log('\nFilter and depth plot options:\n\n') 
#    
#    infidelity, args = grace.get_option_value(args,'--infidelity', int, 50,
#        log, '%d',
#        'Any runner-up hits scoring this close to the best hit\'s score are considered possibly valid hits. '
#        'Note: each extra base mismatch loses 25 points')
#    is_monogamous, args = grace.get_option_value(args,'--monogamous', grace.as_bool, True,
#        log, grace.describe_bool,
#        'Discard reads with more than one possibly valid hit.')    
#    output_userplots, args = grace.get_option_value(args, '--userplots', grace.as_bool, True,
#        log, grace.describe_bool,
#        'Output Artermis userplots.')
#    output_strand_specific_depths, args = grace.get_option_value(args, '--strand-specific', grace.as_bool, False,
#        log, grace.describe_bool,
#        'Output strand-specific depths.')
#    use_random, args = grace.get_option_value(args, '--random', grace.as_bool, False,
#        log, grace.describe_bool,
#        'Instead of using all ambiguous hits, pick one at random from the equal top alignments. Affects ambiguous depth plots, and base calling if --monogamous is no.')
#
#    # Secret option: Output to an alternate directory
#    # This is just a way to produce alternate filterings, it won't work downstream
#    output_dir, args = grace.get_option_value(args, '--output', str, None)
#
#    return (infidelity,is_monogamous,output_userplots,output_strand_specific_depths,use_random,output_dir), args
#
#def filter_main(args):
#    grace.require_samtools()
#
#    log = grace.Log()
#
#    filter_args, args = filter_parse_args(args, log)
#
#    grace.expect_no_further_options(args)                    
#    if len(args) != 1:
#        log.log(FILTER_HELP)
#        raise grace.Help_shown()
#
#    log.attach(open(os.path.join(args[0], 'consensus_log.txt'), 'wb'))
#    
#    return filter(args[0], *filter_args, **{'log':log})

#def zeros(size):
#    return array.array('l', chr(0)*(8*size))
#
#def increment_range(array, start, end):
#    #for i in xrange(start,end):
#    #    array[i] += 1
#    
#    #For speed we will store the derivative of the depth, and integrate it at the end
#    array[start] += 1
#    if end < len(array):
#        array[end] -= 1    
#
#def integrate(a):
#    total = 0
#    for i in xrange(len(a)):
#        total += a[i]
#        a[i] = total

class Hit_already_seen(object): pass

def filter(working_dir, infidelity_snps, is_monogamous, output_userplots, output_strand_specific_depths, use_random, log):    
    workspace = working_directory.Working(working_dir, must_exist=True)
    reference = workspace.get_reference()
    
    infidelity = infidelity_snps * workspace.param.get('snp-cost',25)
    
    #if 'shrimp_cutoff' in workspace.param:
    #    calc_cutoff = samshrimp.cutoff_interpreter(workspace.param['shrimp_cutoff'])
    #else:
    #    calc_cutoff = lambda length: -1

    depths = { }
    for name, length in reference.get_lengths():
        depths[name] = Ref_depth()
        depths[name].ambiguous_depths = [ Depth(length) for direction in (0,1) ]
        depths[name].ambiguous_pairspan_depths = [ Depth(length) for direction in (0,1) ]
        depths[name].depths = [ Depth(length) for direction in (0,1) ]
        depths[name].pairspan_depths = [ Depth(length) for direction in (0,1) ]

    input_filename = workspace.object_filename('alignments.bam')
    
    writer = sam.Bam_writer(
        workspace.object_filename('alignments_filtered.bam'), 
        sam.bam_headers(input_filename)
    )
    
    f_unmapped_single = workspace.open('unmapped_single.fq', 'wb')
    f_unmapped_paired = workspace.open('unmapped_paired.fq', 'wb')

    n = 0
    n_kept_single = 0
    n_kept_paired = 0
    
    fragment_count = 0
    fragment_size_total = 0.0
    fragment_size_total2 = 0.0
    fragment_insert_total = 0.0
    
    n_polygamous = 0
    #n_possibly_polygamous = 0
    any_pairs = False
    
    n_single_unmapped = 0
    n_paired_unmapped = 0

    for read_name, fragment_alignments, unmapped in sam.bam_iter_fragments(input_filename,'Filtering'):
        if unmapped:
            # Shrimp can output multiple copies of an unmapped read, if its pair maps
            true_unmapped = [ ]
            seen = set()
            for item in fragment_alignments:
                for item2 in item:
                    seen.add( item2.original_name() )
            for item in unmapped:
                original_name = item.original_name()
                if original_name in seen: continue
                seen.add(original_name)
                true_unmapped.append(item)
            unmapped = true_unmapped
        
        if unmapped:
            for alignment in unmapped:
                writer.write(alignment)
        
            if len(unmapped) == 1:
                io.write_fastq(f_unmapped_single, unmapped[0].original_name(), unmapped[0].seq, unmapped[0].get_qual())
                n_single_unmapped += 1
            else:
                assert len(unmapped) == 2
                unmapped.sort(key=lambda al: al.flag & sam.FLAG_SECOND)
                io.write_fastq(f_unmapped_paired, unmapped[0].original_name(), unmapped[0].seq, unmapped[0].get_qual())
                io.write_fastq(f_unmapped_paired, unmapped[1].original_name(), unmapped[1].seq, unmapped[1].get_qual())
                n_paired_unmapped += 1
        
        if not fragment_alignments: continue

        n += 1
        
        fragment_alignments = [ (sum( al.get_AS() for al in item ), item) for item in fragment_alignments ]
        
        fragment_alignments.sort(key=lambda x:x[0], reverse=True)
        
        # Apply infidelity
        i = 0
        while i < len(fragment_alignments) and \
              fragment_alignments[i][0] >= fragment_alignments[0][0]-infidelity:
            i += 1 
        fragment_alignments = fragment_alignments[:i]        
        fragment_is_monogamous = len(fragment_alignments) == 1
        if use_random and not fragment_is_monogamous:
            #Pick randomly from equal-top alignments
            i = 0
            while i < len(fragment_alignments) and \
                  fragment_alignments[i][0] >= fragment_alignments[0][0]:
                i += 1
            fragment_alignments = [ random.choice(fragment_alignments[:i]) ]

        for fragment_alignment in fragment_alignments:
            fragment_strand = (1 if fragment_alignment[1][0].flag & sam.FLAG_REVERSE else 0)
        
            for alignment in fragment_alignment[1]:
                depths[alignment.rname].ambiguous_depths[fragment_strand].increment(
                    alignment.pos-1,alignment.pos-1+alignment.length)
                      #was: [1 if alignment.flag & FLAG_REVERSE else 0]
            if len(fragment_alignment[1]) == 2:
                any_pairs = True
                start = min(fragment_alignment[1][0].pos-1,
                            fragment_alignment[1][1].pos-1)
                end = max(fragment_alignment[1][0].pos-1+fragment_alignment[1][0].length,
                          fragment_alignment[1][1].pos-1+fragment_alignment[1][1].length)                
                depths[alignment.rname].ambiguous_pairspan_depths[fragment_strand].increment(start,end)
        
        if not fragment_is_monogamous:
            n_polygamous += 1
            #if is_monogamous: continue
        
        fragment_bases = sum( len(item.seq) for item in fragment_alignments[0][1] )
        #if fragment_alignments[0][0] < calc_cutoff(fragment_bases) + infidelity:
        #    fragment_is_monogamous = False
        #    n_possibly_polygamous += 1
        #    #if is_monogamous: continue
        
        #Emit
        #TODO: mark non-primary alignments
        pair = False
        single = False
        for fragment_alignment in fragment_alignments:
            fragment_strand = (1 if fragment_alignment[1][0].flag & sam.FLAG_REVERSE else 0)
            
            if not is_monogamous or fragment_is_monogamous: 
                for alignment in fragment_alignment[1]:
                        writer.write(alignment)
                if len(fragment_alignment[1]) == 1:
                    single = True
                else:
                    pair = True
                    any_pairs = True
                    
            if fragment_is_monogamous:
                for alignment in fragment_alignment[1]:
                        depths[alignment.rname].depths[fragment_strand].increment(
                            alignment.pos-1,alignment.pos-1+alignment.length)
            
            if fragment_is_monogamous and len(fragment_alignment[1]) == 2:
                start = min(fragment_alignment[1][0].pos-1,
                            fragment_alignment[1][1].pos-1)
                end = max(fragment_alignment[1][0].pos-1+fragment_alignment[1][0].length,
                          fragment_alignment[1][1].pos-1+fragment_alignment[1][1].length)                
                depths[alignment.rname].pairspan_depths[fragment_strand].increment(start,end)
                
                fragment_count += 1
                fragment_size_total += end-start
                fragment_size_total2 += (end-start)**2
                fragment_insert_total += end-start - fragment_alignment[1][0].length - fragment_alignment[1][1].length
                
        if pair:                
            n_kept_paired += 1
        elif single:
            n_kept_single += 1

    writer.close()
    
    f_unmapped_single.close()
    f_unmapped_paired.close()

    #for name, seq in io.read_fasta(reference_filename):
    #    for i in (0,1):
    #        if is_monogamous:
    #            integrate(depths[name].ambiguous_depths[i])
    #            integrate(depths[name].ambiguous_pairspan_depths[i])
    #        integrate(depths[name].depths[i])
    #        integrate(depths[name].pairspan_depths[i])

    log.log('\n')    
    if n_paired_unmapped:
       #log.log('%15s' % grace.pretty_number(n_paired_unmapped) + ' unmapped pairs\n')
       log.datum(workspace.name, 'unmapped pairs', n_paired_unmapped)
    if n_single_unmapped:
       #log.log('%15s' % grace.pretty_number(n_single_unmapped) + ' unmapped reads\n')
       log.datum(workspace.name, 'unmapped reads', n_single_unmapped)
    #log.log(
    #   '%15s' % grace.pretty_number(n) + ' reads' + ('/pairs' if any_pairs else '') + ' with alignments\n' +
    #   '%15s' % grace.pretty_number(n_polygamous) + ' hit multiple locations' + (' (discarded)' if is_monogamous else '') + '\n'
    #)
    log.datum(workspace.name, 'reads' + ('/pairs' if any_pairs else '') + ' with alignments', n)
    log.datum(workspace.name, 'hit multiple locations' + (' (discarded)' if is_monogamous else ''), n_polygamous)

    if n_kept_paired:
       #log.log('%15s' % grace.pretty_number(n_kept_paired) + ' pairs kept\n')
       log.datum(workspace.name, 'pairs kept', n_kept_paired)
    if n_kept_single:
       #log.log('%15s' % grace.pretty_number(n_kept_single) + ' reads kept\n')
       log.datum(workspace.name, 'reads kept', n_kept_single)

    total_length = 0
    total_depth = 0
    total_ambiguous_depth = 0
    for name, length in reference.get_lengths():
        total_length += length
        total_depth += depths[name].depths[0].total() + depths[name].depths[1].total()
        total_ambiguous_depth += depths[name].ambiguous_depths[0].total() + depths[name].ambiguous_depths[1].total()

    log.datum(workspace.name, 'average depth of coverage, ambiguous', float(total_ambiguous_depth)/total_length)
    log.datum(workspace.name, 'average depth of coverage, unambiguous', float(total_depth)/total_length)
    
    log.log('\n')
    
    if fragment_count:
        #log.log('%15.1f' % (float(fragment_insert_total)/fragment_count) + ' mean insert size (unsequenced middle of fragment)\n')
        mean = float(fragment_size_total)/fragment_count
        var = float(fragment_size_total2)/fragment_count - mean*mean
        #log.log('%15.1f' % mean + ' mean fragment size\n')
        #log.log('%15.1f' % (var**0.5) + ' s.d. fragment size\n')
        log.datum(workspace.name, 'mean fragment size', mean)
        log.datum(workspace.name, 's.d. fragment size', (var**0.5))
        log.log('\n')    
    
    sam.sort_and_index_bam(
        workspace.object_filename('alignments_filtered.bam'),
        workspace.object_filename('alignments_filtered_sorted')
    )
    
    # Write depths
    grace.status('Write depth pickle')
    workspace.set_object(depths, 'depths.pickle.gz')
    workspace.update_param(any_pairs = any_pairs)
    grace.status('')
    
    if output_userplots:
        for name in depths:
            grace.status('Write depth for ' + name)
            
            write_userplots(
                os.path.join(working_dir, grace.filesystem_friendly_name(name) + '-depth'),
                depths[name].depths,
                output_strand_specific_depths
            )
            
            grace.status('Write ambiguous depth for ' + name)
            
            write_userplots(
                os.path.join(working_dir, grace.filesystem_friendly_name(name) + '-ambiguous-depth'),
                depths[name].ambiguous_depths,
                output_strand_specific_depths
            )
            
            if any_pairs:        
                grace.status('Write pair-span depth for ' + name)
                
                write_userplots(
                    os.path.join(working_dir, grace.filesystem_friendly_name(name) + '-pairspan-depth'),
                    depths[name].pairspan_depths,
                    output_strand_specific_depths
                )
                
                grace.status('Write ambiguous pair-span depth for ' + name)
                
                write_userplots(
                    os.path.join(working_dir, grace.filesystem_friendly_name(name) + '-ambiguous-pairspan-depth'),
                    depths[name].ambiguous_pairspan_depths,
                    output_strand_specific_depths
                )
            
            grace.status('')


# ===================================================================================
# Handling of BioPython sequence features        

def iter_feature_positions(feature):
    if feature.sub_features:
        if feature.strand >= 0:
            sub_features = feature.sub_features
        else:
            sub_features = feature.sub_features[::-1]
            
        for sub_feature in feature.sub_features:
            for item in iter_feature_positions(sub_feature):
                yield item
    else:
        loc = feature.location
        if feature.strand >= 0:
            yield loc.nofuzzy_start, 1, False
            for i in xrange(loc.nofuzzy_start+1, loc.nofuzzy_end):
                yield i, 1, True
                yield i, 1, False
        else:
            yield loc.nofuzzy_end-1, -1, False
            for i in xrange(loc.nofuzzy_end-2, loc.nofuzzy_start-1, -1):
                yield i+1, -1, True
                yield i, -1, False

TRUMPS = {
    'mRNA' : ('CDS'),
    'gene' : ('mRNA', 'CDS'),
}

# =====================================================
# Reference sequence object


class Ref_seq(object):
    def __init__(self, seq):
        self.seq = seq
        #self.depth = [ 0 ] * len(seq)
        self.stranded_base_counts = [ [ consensus.EMPTY_EVIDENCE ] * len(seq) for strand in (0,1) ]
        self.stranded_insertions_before  = [ [ consensus.EMPTY_EVIDENCE ] * (len(seq)+1) for strand in (0,1) ]
    
        self.base_annotation = [ () ]*len(seq)
        self.insertion_annotation = [ () ]*len(seq) 
    
    def process_alignment(self, alignment, trim, whole_read_only):    
         #Note: orientation of read, not of fragment
         if alignment.flag & sam.FLAG_REVERSE:
             strand = 1
         else:
             strand = 0

         pos = alignment.pos - 1

         if whole_read_only:
              #No soft clipping
              if pos > 0 and alignment.cigar.lstrip('0123456789').startswith('S'):
                  return
              if pos+alignment.length < len(self.seq) and alignment.cigar.endswith('S'):
                  return

         seq = alignment.seq
         seq_pos = 0
         n = 0         
         
         trim_start = pos
         if trim_start > 0:
             trim_start += trim
                          
         trim_end = min(pos + alignment.length, len(self.seq))
         if trim_end < len(self.seq):
             trim_end -= trim
         
         base_counts = self.stranded_base_counts[strand]
         insertions_before = self.stranded_insertions_before[strand]
                  
         evidence_add = consensus.evidence_add
         
         for value in array.array('B', alignment.cigar):
            if 48 <= value <= 57:
                n = n*10+(value-48)
            else:
                if value == ord('M') or value == ord('=') or value == ord('X'):
                    #Match or substitution
                    for i in xrange(n):
                        if trim_start <= pos < trim_end:
                            base_counts[pos] = evidence_add(base_counts[pos], seq[seq_pos], 1)
                        pos += 1
                        seq_pos += 1
                elif value == ord('D'):
                    for i in xrange(n):
                        if trim_start <= pos < trim_end:
                            base_counts[pos] = evidence_add(base_counts[pos], '-', 1)
                        pos += 1                    
                elif value == ord('I'):
                    if trim_start <= pos <= trim_end:
                        insertions_before[pos] = evidence_add(
                            insertions_before[pos],
                            seq[seq_pos:seq_pos+n],
                            1
                        )
                    seq_pos += n
                elif value == ord('S'):
                    #Soft clipping
                    seq_pos += n
                elif value == ord('H'):
                    pass
                else:
                    raise grace.Error('Unhandled cigar character: %s' % char)
                
                n = 0    
    
    def call_consensus(self, p_cutoff, stranded_p_cutoff, indel_prior, prior_weight, use_ambiguity_codes, proportion):
        total_prior = prior_weight
        
        empty_prior    = indel_prior * total_prior
        nonempty_prior = (total_prior-empty_prior) / 4.0 # Each base equally likely
        
        insertion_present_prior = indel_prior * total_prior
        insertion_absent_prior = total_prior - insertion_present_prior
        
        self.base_counts = [ ]
        self.insertions_before = [ ]

        self.base_calls = [ ] # Uppercase ambiguity code or N
        self.insertion_calls = [ ]  # Uppercase / - or N
        #self.insertion_calls_without_purity = [ ]  #Uppercase (pure) or lowercase / - (impure) or N
        
        self.consensus = [ ]
        self.consensus_masked = [ ]
        
        self.alignment_reference = [ ]
        self.alignment_result = [ ]
        
        self.report = [ ]
        
        self.expected_miscalls = 0.0
        self.expected_miscalled_changes = 0.0
        
        for i in xrange(len(self.seq)):
            for j in (0,1):
                # How many alignments didn't have insertions?
                total = consensus.evidence_total_count(self.stranded_insertions_before[j][i])
                depth = consensus.evidence_total_count(self.stranded_base_counts[j][i])
                if i: depth = min(depth, consensus.evidence_total_count(self.stranded_base_counts[j][i-1]))
                if depth > total: 
                    self.stranded_insertions_before[j][i] = consensus.evidence_add(self.stranded_insertions_before[j][i], '-', depth-total)
        
            self.base_counts.append( consensus.evidence_merge( self.stranded_base_counts[0][i],self.stranded_base_counts[1][i] ) )
            self.insertions_before.append( consensus.evidence_merge( self.stranded_insertions_before[0][i],self.stranded_insertions_before[1][i] ) )
        
            insertion_call_with_purity, p = consensus.bayesian_consensus(self.insertions_before[i], p_cutoff,insertion_present_prior,insertion_absent_prior,total_prior, proportion)
            if insertion_call_with_purity is not None:
                self.expected_miscalls += 1.0-p
                if insertion_call_with_purity != '-':
                    self.expected_miscalled_changes += 1.0-p
            
            if stranded_p_cutoff and insertion_call_with_purity is not None:
                call1, p1 = consensus.bayesian_consensus(self.stranded_insertions_before[0][i], stranded_p_cutoff,insertion_present_prior,insertion_absent_prior,total_prior, proportion)
                call2, p2 = consensus.bayesian_consensus(self.stranded_insertions_before[1][i], stranded_p_cutoff,insertion_present_prior,insertion_absent_prior,total_prior, proportion)
                if insertion_call_with_purity != call1 or insertion_call_with_purity != call2:
                    insertion_call_with_purity = None
            
            #2/10/2012 - this isn't useful
            #
            # insertion_call_without_purity, p = consensus.bayesian_consensus(self.insertions_before[i], 0.0,insertion_present_prior,insertion_absent_prior,total_prior, proportion)
            #
            # #Indicate lack of purity with lower case in alignment and consensus files
            # if insertion_call_without_purity is not None and insertion_call_with_purity is None:
            #    insertion_call_without_purity = insertion_call_without_purity.lower()
            #
            # if insertion_call_without_purity is not None and insertion_call_without_purity != '-': 
            #     self.consensus.append(insertion_call_without_purity)
            #     self.consensus_masked.append(insertion_call_without_purity)
            #     
            #     self.alignment_result.append(insertion_call_without_purity)
            #     self.alignment_reference.append('-' * len(insertion_call_without_purity))

            if insertion_call_with_purity is not None and insertion_call_with_purity != '-':
                self.consensus.append(insertion_call_with_purity)
                self.consensus_masked.append(insertion_call_with_purity)
                
                self.alignment_result.append(insertion_call_with_purity)
                self.alignment_reference.append('-' * len(insertion_call_with_purity))

                self.report.append(('insertion-before', i, '-', insertion_call_with_purity, self.insertions_before[i],self.stranded_insertions_before[0][i],self.stranded_insertions_before[1][i]))

            self.insertion_calls.append(insertion_call_with_purity or 'N')

            #insertion_call = bio.consensus(self.insertions_before[i], min_depth,min_purity) or 'N'
            #
            #if insertion_call != 'N':
            #    insertion_call_without_purity = insertion_call
            #else:
            #    insertion_call_without_purity = (bio.consensus(self.insertions_before[i], min_depth, 0.0) or '-').lower()
            #    # NOTE: defaults to saying no insertion if we can't even get consensus here
            #self.insertion_calls.append(insertion_call)
            #if insertion_call_without_purity != '-':
            #     consensus.append(insertion_call_without_purity)
            #     consensus_masked.append(insertion_call_without_purity)
            #     alignment_reference.append('-' * len(insertion_call_without_purity))
            #     alignment_result.append(insertion_call_without_purity)

            if use_ambiguity_codes:
                substitution_call, p = consensus.bayesian_ambiguity_code_consensus(self.base_counts[i], p_cutoff,nonempty_prior,empty_prior,total_prior, proportion)
                if stranded_p_cutoff and substitution_call is not None:
                    call1, p1 = consensus.bayesian_ambiguity_code_consensus(self.stranded_base_counts[0][i], stranded_p_cutoff,nonempty_prior,empty_prior,total_prior, proportion)
                    call2, p2 = consensus.bayesian_ambiguity_code_consensus(self.stranded_base_counts[1][i], stranded_p_cutoff,nonempty_prior,empty_prior,total_prior, proportion)
                    substitution_call = consensus.ambiguity_merge(substitution_call, call1, call2)
            
            else:
                substitution_call, p = consensus.bayesian_consensus(self.base_counts[i], p_cutoff,nonempty_prior,empty_prior,total_prior, proportion)
                if stranded_p_cutoff and substitution_call is not None:
                    call1, p1 = consensus.bayesian_consensus(self.stranded_base_counts[0][i], stranded_p_cutoff,nonempty_prior,empty_prior,total_prior, proportion)
                    call2, p2 = consensus.bayesian_consensus(self.stranded_base_counts[1][i], stranded_p_cutoff,nonempty_prior,empty_prior,total_prior, proportion)
                    if substitution_call != call1 or substitution_call != call2:
                        substitution_call = None
            
            if substitution_call is not None:
                self.expected_miscalls += 1.0-p
                if not bio.might_be_same_base(self.seq[i], substitution_call):
                    self.expected_miscalled_changes += 1.0-p
            
            if substitution_call is None:
                self.consensus.append('N')
                self.consensus_masked.append( self.seq[i].lower() )
                #has_consensus.append(False)
            elif substitution_call == '-':
                self.report.append(('deletion',i, self.seq[i], '-', self.base_counts[i],self.stranded_base_counts[0][i],self.stranded_base_counts[1][i]))
                #has_consensus.append(True)
            else:
                self.consensus.append(substitution_call)
                self.consensus_masked.append(substitution_call)
                #has_consensus.append(substitution_call in 'ACGT') #Exclude ambiguity codes
                
                if not bio.might_be_same_base(self.seq[i], substitution_call):
                    self.report.append(('substitution', i, self.seq[i], substitution_call, self.base_counts[i],self.stranded_base_counts[0][i],self.stranded_base_counts[1][i]))
            
            self.base_calls.append(substitution_call or 'N')        
            self.alignment_result.append(substitution_call or 'N')
            self.alignment_reference.append(self.seq[i])

            #if use_ambiguity_codes:
            #    base_call = bio.ambiguity_code_consensus(self.base_counts[i], min_depth,min_purity) or 'N'
            #else:
            #    base_call = bio.consensus(self.base_counts[i], min_depth,min_purity) or 'N'
            #            
            #self.base_calls.append(base_call)            
            #                        
            #consensus.append(base_call)
            #if base_call == 'N':
            #    consensus_masked.append(self.seq[i])
            #else:
            #    consensus_masked.append(base_call)            
            #alignment_reference.append(self.seq[i])
            #alignment_result.append(base_call)
        
        self.consensus = ''.join(self.consensus)
        self.consensus_masked = ''.join(self.consensus_masked)
        self.alignment_reference = ''.join(self.alignment_reference)
        self.alignment_result = ''.join(self.alignment_result)


    def annotate_consensus(self, record, default_transl_table):
        assert record.seq.tostring().upper() == self.seq.upper(), 'Annotation sequence does not match. Huh?'
        annotations = { } 
        for feature in record.features:
            self.annotate_consensus_feature(annotations, feature, default_transl_table)
        for (pos,is_insertion,feature_type),text in annotations.items():
            for trump in TRUMPS.get(feature_type,()):
                if (pos,is_insertion,trump) in annotations: break
            else:
                array = self.insertion_annotation if is_insertion else self.base_annotation
                array[pos] = array[pos] + (text,)

    def annotate_consensus_feature(self, annotations, feature, default_transl_table):
        from Bio import Seq
        
        if feature.type == 'source': return

        codon_start = 0

        if feature.type == 'CDS':
            if 'transl_table' in feature.qualifiers:
                transl_table_no = int(feature.qualifiers['transl_table'][0])
            else:
                assert default_transl_table is not None, 'No /transl_table for CDS, and default transl_table not given'
                transl_table_no = default_transl_table
            
            if 'codon_start' in feature.qualifiers:
                codon_start = int(feature.qualifiers['codon_start'][0]) - 1
            
            def translate(seq, is_first):
                result = Seq.Seq(seq).translate(table=transl_table_no).tostring()
                if is_first and result: result = 'M'+result[1:]
                return result
            codon_size = 3
        else:
            def translate(seq, is_first):
                return seq    
            codon_size = 1

        if 'locus_tag' not in feature.qualifiers:
            locus_tag = '%d..%d' % (feature.location.nofuzzy_start+1,feature.location.nofuzzy_end)
        else:
            locus_tag = feature.qualifiers['locus_tag'][0]

        #if feature.type != 'CDS': locus_tag = feature.type+':'+locus_tag

        shift = 0
        phase = 0
        
        block_changes = []
        block_start = 0
        dna_block = []
        new_dna_block = []
        
        i = codon_start        
        position_iter = iter_feature_positions(feature)        
        for j in xrange(codon_start): position_iter.next()
        finished = False
        
        while True:
            try:
                p, strand, insertion = position_iter.next()
            except StopIteration:
                finished = True
        
            if (phase == 0 and shift == 0) or finished:
                #Common codon boundary
                dna_block = ''.join(dna_block)
                new_dna_block = ''.join(new_dna_block)
                
                if not bio.might_be_same_bases(dna_block, new_dna_block):
                    protein = translate(dna_block, block_start==codon_start)
                    new_protein = translate(new_dna_block, block_start==codon_start)
                    
                    cod_start = (block_start-codon_start)//codon_size
                    cod_end = (i-codon_start)//3
                    if codon_size == 1 or cod_start+1 == cod_end:
                        cod_desc = '' #'codon %d' % (cod_start+1)
                    else:
                        cod_desc = ' of codons %d..%d' % (cod_start+1,cod_end)
                    
                    extra = ''
                    if 'product' in feature.qualifiers:
                        extra += ' ' + (' '.join(feature.qualifiers['product']).replace(',',';'))
                    
                    if new_protein == protein:
                        what = 'synonymous'
                    elif shift != 0:
                        what = 'frame-shift'
                    else:
                        what = '%s=>%s' % (protein or '-', new_protein or '-')                 
                    
                    for pos, is_insertion, ref_pos in block_changes:
                        where = '%sbase %d' % (
                            'before ' if is_insertion else '',
                            pos+1
                        )
                        if codon_size != 1: 
                            where += '%s codon %d' % (
                                ' before' if is_insertion and (pos-codon_start)%3 == 0 else '', 
                                (pos-codon_start)//codon_size+1
                            )
                        
                        item = feature.type + ' ' + what + ' ' + locus_tag + ' ' + where + cod_desc + extra
                                                
                        annotations[(ref_pos, is_insertion, feature.type)] = item                     
                
                block_start = i
                block_changes = [ ]
                dna_block = [ ]
                new_dna_block = [ ]
            
            if finished: break
            
            if insertion:
                dna = ''
                new_dna = self.insertion_calls[p]
                if new_dna in ('-','N'): 
                    new_dna = ''
                else:
                    block_changes.append( (i,True,p) )
            else:
                dna = self.seq[p]
                new_dna = self.base_calls[p]
                if new_dna == '-': 
                    new_dna = ''
                    block_changes.append( (i,False,p) )
                elif not bio.might_be_same_base(dna, new_dna):
                    block_changes.append( (i,False,p) )
            
            if strand < 0:
                dna = bio.reverse_complement(dna)
                new_dna = bio.reverse_complement(new_dna)
            
            dna_block.append(dna)
            new_dna_block.append(new_dna)
            shift = (shift + len(new_dna)-len(dna))%codon_size
            phase = (phase+len(dna))%codon_size
            i += len(dna)
            
            #
            #if strand >= 0:
            #    dna_block.append( self.seq[p] )
            #
            #    new = self.base_calls[p]
            #    if new == '-':
            #        shift += 1
            #        block_changes.append( (i, 'deletion', p) )                    
            #    else:
            #        new_dna_block.append(new)
            #        if not bio.might_be_same_bases(new_dna_block[-1], dna_block[-1]):
            #            block_changes.append( (i, 'substitution', p) )
            #        
            #    if p+1 < len(self.seq) and self.insertion_calls[p+1] not in ('N','-'):
            #        new_dna_block.append( self.insertion_calls[p+1] )
            #        shift -= len(self.insertion_calls[p+1])
            #        block_changes.append( (i+1, 'insertion', p+1) )
            #        
            #else:
            #    dna_block.append( bio.reverse_complement(self.seq[p]) )
            #
            #    new = self.base_calls[p]
            #    if new == '-':
            #        shift += 1
            #        block_changes.append( (i, 'deletion', p) )
            #    else:
            #        new_dna_block.append(bio.reverse_complement(new))
            #        if not bio.might_be_same_bases(new_dna_block[-1], dna_block[-1]):
            #            block_changes.append( (i, 'substitution', p) )
            #
            #    if p >= 0 and self.insertion_calls[p-1] not in ('N','-'):
            #        new_dna_block.append( bio.reverse_complement(self.insertion_calls[p-1]) )
            #        shift -= len(self.insertion_calls[p-1])
            #        block_changes.append( (i+1, 'insertion', p-1) )
            #
            #shift = shift%codon_size
            #phase = (phase+1)%codon_size
            #i += 1 




CONSENSUS_BLURB = """\
Consensus calling, short version: 

Adjust --cutoff to obtain reasonable depth cutoffs. \
You should be able to leave --indel-prior and --prior-weight alone.

If you might have clonal reads, also set --strand-cutoff. This lets \
you specify an aditional strand specific cutoff.

Consensus calling, long version:

A priori we believe each position in the reference to be some mixture of \
A, C, G, T, and deletion, but we are unsure what exact mixture it is. \
This prior belief is expressed as a Dirichlet distribution. \
As we observe bases from various reads, our beliefs are updated. \
A posteriori, if we believe with probability >= cutoff that there \
is a certain base (or a deletion) that forms a majority in the \
mixture, it is called as the consensus.

(By default the required majority is 50%, but this can also be adjusted. \
For example if you want to avoid calling a consensus for mixed populations, \
you might set --majority higher.)

For insertions, a priori we believe there to be a mixture of insertion and \
no-insertion, and the same process occurs. 
"""

@config.help("""\

Make consensus calls on a working directory previously filtered by
"filter:" or "consensus:".


""" + CONSENSUS_BLURB)
@config.Bool_flag('whole_read_only', 'Only use alignments to whole read.')
@config.Int_flag('trim', 'Trim alignment start and end by this many bases. '
                         '(Reduces SNP/indel misclassifications.)')
@config.Float_flag('cutoff', 'Probability of majority required to call consensus.')
@config.Float_flag('strand_cutoff', 'Additional cutoff applied to each strand on its own.')
@config.Float_flag('indel_prior', 'Prior expected proportion of insertion or deletion at each base.')
@config.Float_flag('prior_weight', 'Our prior belief is as though we have this many existing observations.')
@config.Float_flag('majority', 'Required majority, higher is stricter, must be less than 1.')
@config.Bool_flag('ambiguity_codes', 'Use IUPAC ambiguity codes.')
@config.Int_flag('transl_table', 'Default codon translation table.') 
class Reconsensus(config.Action_with_working_dir):
    whole_read_only = False
    trim = 5
    cutoff = 0.99
    strand_cutoff = 0.0
    indel_prior = 0.2
    prior_weight = 1.0
    majority = 0.5
    ambiguity_codes = False
    transl_table = 11
    
    _workspace_class = working_directory.Working

    def describe(self, *args, **kwargs):
        desc = super(Reconsensus, self).describe(*args, **kwargs)
        
        if not kwargs.get('brief',False):
            desc += '\n' + consensus.consensus_calling_advice(self.cutoff, self.indel_prior, self.prior_weight, proportion=self.majority)    
            if self.strand_cutoff:
                desc += '\n' + consensus.consensus_calling_advice(self.strand_cutoff, self.indel_prior, self.prior_weight, 'The per-strand coverage required is:\n', proportion=self.majority)
        return desc    
    
    def run(self):
        invocation = config.strip_color(self.describe())
        
        consensus_run(
            invocation=invocation, filter_needed=False, log=self.log, whole_read_only=self.whole_read_only, 
            trim=self.trim, p_cutoff=self.cutoff, stranded_p_cutoff=self.strand_cutoff, 
            indel_prior=self.indel_prior, prior_weight=self.prior_weight, proportion=self.majority, 
            use_ambiguity_codes=self.ambiguity_codes, default_transl_table=self.transl_table, 
            working_dir=self.working_dir
        )


@config.help("""\

Make consensus calls for SNPs and indels based on the directory created by "shrimp:" or "import:".

""" + CONSENSUS_BLURB)
class Consensus(Filter, Reconsensus):
    def run(self):
        Filter.run(self)
        Reconsensus.run(self)



#CONSENSUS_HELP = """\
#
#Consensus calling, short version: 
#   Adjust --cutoff to obtain reasonable depth cutoffs.
#   You should be able to leave --indel-prior and --prior-weight alone.
#   
#   If you might have clonal reads, also set --strand-cutoff. This lets
#   you specify an aditional strand specific cutoff.
#
#Consensus calling, long version:   
#   A priori we believe each position in the reference to be some mixture of 
#   A, C, G, T, and deletion, but we are unsure what exact mixture it is.
#   This prior belief is expressed as a Dirichlet distribution. 
#   As we observe bases from various reads, our beliefs are updated.
#   A posteriori, if we believe with probability >= cutoff that there 
#   is a certain base (or a deletion) that forms a majority in the 
#   mixture, it is called as the consensus. 
#   
#   (By default the required majority is 50%, but this can also be adjusted.
#   For example if you want to avoid calling a consensus for mixed populations,
#   you might set --majority higher.)
#   
#   For insertions, a priori we believe there to be a mixture of insertion and 
#   no-insertion, and the same process occurs. 
#
#Usage:
#
#    nesoni samconsensus: working_dir [options]
#
#where working_dir was created by "nesoni samshrimp".
#
#"""
#
#def consensus_main(args, filter_needed=True):        
#    grace.require_samtools()
#
#    log = grace.Log()
#    
#    #original_args = args    
#    if filter_needed:
#        invocation = 'samconsensus ' + ' '.join(args)
#    else:
#        invocation = 'samreconsensus ' + ' '.join(args)
#    
#    if filter_needed:
#        filter_args, args = filter_parse_args(args, log)
#    
#    log.log('\nConsensus calling options:\n\n')
#
#    whole_read_only, args = grace.get_option_value(args, '--whole-read-only', grace.as_bool, False,
#        log, grace.describe_bool,
#        'Only use alignments to whole read.')
#
#    trim, args = grace.get_option_value(args,'--trim', int, 5,
#        log, '%d',
#        'Trim alignment start and end by this many bases. (Reduces SNP/indel misclassification.)')
#
#    p_cutoff, args = grace.get_option_value(args,'--cutoff', float, 0.99,
#        log, '%.3f',
#        'Probability of majority required to call consensus.')
#    assert 0.0 <= p_cutoff < 1.0, '--cutoff specified is not in range [0,1)'
#
#    stranded_p_cutoff, args = grace.get_option_value(args,'--strand-cutoff', float, 0.0,
#        log, '%.3f',
#        'Additional cutoff applied to each strand on its own.')
#    assert 0.0 <= stranded_p_cutoff < 1.0, '--strand-cutoff specified is not in range [0,1)'
#
#    indel_prior, args = grace.get_option_value(args,'--indel-prior', float, 0.2,
#        log, '%.2f',
#        'Prior expected proportion of insertion or deletion at each base.')
#
#    prior_weight, args = grace.get_option_value(args,'--prior-weight', float, 1.0,
#        log, '%.2f',
#        'Our prior belief is as though we have this many existing observations.')
#    assert prior_weight > 0.0, '--prior-weight is not greater than zero'
#
#    proportion, args = grace.get_option_value(args,'--majority', float, 0.5,
#        log, '%.3f',
#        'Required majority, higher is stricter, must be less than 1.')    
#    assert 0.0 < proportion < 1.0, '--majority specified is not in range (0,1)' 
#    
#    use_ambiguity_codes, args = grace.get_option_value(args,'--ambiguity-codes', grace.as_bool, True,
#        log, grace.describe_bool,
#        'Use IUPAC ambiguity codes.')
#
#    default_transl_table, args = grace.get_option_value(args, '--transl_table', int, 11,
#        log, '%d',
#        'Default codon translation table.')
#
#    log.log('\n' + consensus.consensus_calling_advice(p_cutoff, indel_prior, prior_weight, proportion=proportion))
#    
#    if stranded_p_cutoff:
#        log.log('\n' + consensus.consensus_calling_advice(stranded_p_cutoff, indel_prior, prior_weight, 'The per-strand coverage required is:\n', proportion=proportion))    
#    
#    grace.expect_no_further_options(args)                    
#    if len(args) != 1:
#        log.log(CONSENSUS_HELP)
#        raise grace.Help_shown()
#
#    working_dir = args[0]
#    
#    log.attach(open(os.path.join(working_dir, 'consensus_log.txt'), 'wb'))
#
#    if filter_needed:
#        filter(working_dir, *filter_args, **{'log':log})
#
#    return consensus_run( 
#      invocation, filter_needed, log, whole_read_only, trim, p_cutoff, stranded_p_cutoff, 
#      indel_prior, prior_weight, proportion, use_ambiguity_codes, default_transl_table, working_dir
#    )


def consensus_run(    
    invocation, filter_needed, log, whole_read_only, trim, p_cutoff, stranded_p_cutoff, 
    indel_prior, prior_weight, proportion, use_ambiguity_codes, default_transl_table, working_dir
    ):

    working = working_directory.Working(working_dir, must_exist=True)
    reference = working.get_reference()
    
    references = { }
    for name, seq in io.read_fasta(reference.reference_fasta_filename()):
        references[name] = Ref_seq( seq.upper() )
    
    bam_filename = working.object_filename('alignments_filtered_sorted.bam')
    
    grace.status('Processing')
    reader = sam.Bam_reader(bam_filename)
    nth = 0
    
    for alignment in reader:
         if alignment.flag & sam.FLAG_UNMAPPED: continue
         
         nth += 1
         if nth % 100000 == 0:
             grace.status('Processing alignment %s' % grace.pretty_number(nth))       
    
         references[alignment.rname].process_alignment(alignment, trim, whole_read_only)
    
    grace.status('Writing results')
    
    consensus_filename = os.path.join(working_dir, 'consensus.fa')
    consensus_file = open(consensus_filename, 'wb')
    
    consensus_masked_filename = os.path.join(working_dir, 'consensus_masked.fa')
    consensus_masked_file = open(consensus_masked_filename, 'wb')
    
    alignment_filename = os.path.join(working_dir, 'alignment.maf')
    alignment_file = open(alignment_filename, 'wb')

    gff_header = (
        '##gff-version 3\n' +
        '##feature-ontology http://song.cvs.sourceforge.net/viewvc/*checkout*/song/ontology/so.obo?revision=1.263\n' + 
        '# date %s\n' % datetime.date.today().isoformat() +
        '# cmdline nesoni\n' +
        config.wrap(invocation, 100, '# ')+'\n' +
        '# nesoni version %s\n' % nesoni.VERSION +
        '# run %s@%s\n' % (os.environ.get('USER','?'), socket.gethostname())
    )

    report_file = open(os.path.join(working_dir, 'report.txt'), 'wb')
    report_gff_file = open(os.path.join(working_dir, 'report.gff'), 'wb')
    report_gff_file.write(gff_header)
    
    report_file.write('Sequence\tPosition in reference\tChange type\tOld\tNew\tEvidence\tConsequences')
    if stranded_p_cutoff: report_file.write('\tForward\tReverse')
    report_file.write('\n')
    
    n_substitutions = 0
    n_insertions = 0
    n_deletions = 0
    
    expected_miscalls = 0.0
    expected_miscalled_changes = 0.0
    
    for rname in sorted(references.keys()):
        ref = references[rname]
        ref.call_consensus(p_cutoff, stranded_p_cutoff, indel_prior, prior_weight, use_ambiguity_codes, proportion)
        
        expected_miscalls += ref.expected_miscalls
        expected_miscalled_changes += ref.expected_miscalled_changes
    
        #annotation_filename = io.abspath(working_dir, grace.filesystem_friendly_name(rname) + '.gbk')
        annotation_filename = reference / (grace.filesystem_friendly_name(rname) + '.gbk')
        if os.path.exists(annotation_filename):
            from Bio import SeqIO
            annotation_file = open(annotation_filename, 'rb')
            annotation = SeqIO.parse(annotation_file, 'genbank').next()
            annotation_file.close()
            ref.annotate_consensus(annotation, default_transl_table)
        
        io.write_fasta(consensus_file, rname, ref.consensus)
        io.write_fasta(consensus_masked_file, rname, ref.consensus_masked)

        alignment_file.write( '\n' )
        alignment_file.write( 'a\n')
        alignment_file.write( 's %s           0 %d + %d %s\n' % (rname, len(ref.seq), len(ref.seq), ref.alignment_reference))
        result_size = len(ref.alignment_result.replace('-',''))
        alignment_file.write( 's %s-consensus 0 %d + %d %s\n' % (rname, result_size, result_size, ref.alignment_result))
        
        friendly_name = grace.filesystem_friendly_name(rname)
        evidence_filename = os.path.join(working_dir, friendly_name + '-evidence.txt')
        individual_report_gff_filename = os.path.join(working_dir, friendly_name + '.gff')
        
        f = open(evidence_filename,'wb')
        f.write( 'Position\tInsertion-before evidence\tSubstitution evidence\t'
                 'Reference\tInsertion call\tSubstitution call\t'
                 'Insertion consequences\tSubstitution consequences\n' )
        for i in xrange(len(ref.seq)):
            f.write( '%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                i+1,
                consensus.pretty_evidence(ref.insertions_before[i]),
                consensus.pretty_evidence(ref.base_counts[i]),
                ref.seq[i],
                ref.insertion_calls[i],
                ref.base_calls[i],
                ', '.join(ref.insertion_annotation[i]),
                ', '.join(ref.base_annotation[i]),
            ))
        f.close()
        
        individual_report_gff_file = open(individual_report_gff_filename, 'wb')
        individual_report_gff_file.write(gff_header) 

        for change_type, position, old, new, counts, counts_forward, counts_reverse in ref.report:
            if change_type == 'insertion-before':
                consequence = ref.insertion_annotation[position]
            else:
                consequence = ref.base_annotation[position]
            consequence = ', '.join(consequence)
            
            report_file.write( '%s\t%d\t%s\t%s\t%s\t%s\t%s' % (
                rname,
                position+1,
                change_type,
                old,
                new,
                consensus.pretty_evidence(counts),
                consequence
            ))
            
            if stranded_p_cutoff:
                report_file.write( '\t%s\t%s' % (
                    consensus.pretty_evidence(counts_forward),
                    consensus.pretty_evidence(counts_reverse),
                ))
            
            report_file.write('\n')
            
            if change_type == 'deletion':
                n_deletions += 1
                start = position+1
                end = position+1
                so_term = 'deletion'
                product = 'Base deleted: %s' % (old)                
            elif change_type == 'insertion-before':
                n_insertions += 1
                #Bracket insertion
                start = position
                end = position+1
                so_term = 'insertion'
                product = 'Insertion: .' + new + '.'
            else:
                n_substitutions += 1
                start = position+1
                end = position+1
                so_term = 'point_mutation'
                product = 'Substitution: %s became %s' % (old,new)

            product += ' ('+consensus.pretty_evidence(counts)+')'

            if consequence:            
                product += ' '+consequence
                        
            gff_line = '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n' % (
                rname,
                'nesoni_samconsensus_' + nesoni.VERSION,
                so_term,
                start,
                end,
                '.', #score
                '+', #strand
                '.', #frame
                'product='+product
            )
            
            report_gff_file.write(gff_line)
            individual_report_gff_file.write(gff_line)

        individual_report_gff_file.close()        

        #depth_filename = os.path.join(working_dir, friendly_name + '-depth.txt')
        #f = open(depth_filename, 'wb')
        #for d in ref.depth:
        #    print >> f, d
        #f.close()

    consensus_file.close()
    consensus_masked_file.close()
    alignment_file.close()
    report_file.close()
    report_gff_file.close()

    grace.status('')

    #log.log('%15s' % grace.pretty_number(n_substitutions) + ' SNPs called\n')
    #log.log('%15s' % grace.pretty_number(n_deletions) + ' deletions called\n')
    #log.log('%15s' % grace.pretty_number(n_insertions) + ' insertions called\n')
    #log.log('%15.3f' % expected_miscalled_changes + ' expected number of false SNP/indel calls\n') 
    log.datum(working.name,'SNPs called',n_substitutions)
    log.datum(working.name,'deletions called',n_deletions)
    log.datum(working.name,'insertions called',n_insertions)
    log.datum(working.name,'expected number of false SNP/indel calls',expected_miscalled_changes) 

    #If you must know, uncomment
    #log.log('%15.2f' % (expected_miscalls-expected_miscalled_changes) + ' expected false no-change calls\n') 

    log.log('\n')

