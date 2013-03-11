
from nesoni import io, bio, grace, config, sam, span_index, annotation, working_directory

import os, sys, collections, random

def subsample(n, p):
    #Note that this is rather inefficient
    result = 0
    for i in xrange(n):
        if random.random() <= p:
            result += 1
    return result

def add_lists(list1, list2):
    assert len(list1) == len(list2)
    return [ a+b for a,b in zip(list1,list2) ]

def add_defdicts(dict1, *other_dicts):
    result = dict1.copy()
    for other_dict in other_dicts:
        for item in other_dict:
            result[item] += other_dict[item]
    return result

class Feature(object):
   """ Fields
   
       name
       count       -    strand -> column_number -> count
       
       common      -   (strand,strand) -> feature -> count    -    hit overlaps multiple features
       ambiguous   -   (strand,strand) -> feature -> count    -    fragment has multiple hits
   """

   def __init__(self, n_samples):
       self.count = {-1: [0]*n_samples, 1:[0]*n_samples }
       #self.count = {-1: [ frag_histogram() for i in xrange(n_samples) ], 
       #               1: [ frag_histogram() for i in xrange(n_samples) ] }
       self.ambiguous = { 
           (-1,-1):collections.defaultdict(int), 
           (-1,1):collections.defaultdict(int), 
           (1,-1):collections.defaultdict(int), 
           (1,1):collections.defaultdict(int), 
       }
       self.common = { 
           (-1,-1):collections.defaultdict(int), 
           (-1,1):collections.defaultdict(int), 
           (1,-1):collections.defaultdict(int), 
           (1,1):collections.defaultdict(int), 
       }
   
   def __repr__(self): 
       return self.name
   
   #def add(self, other):
   #    for strand in (-1,1):
   #        self.count[strand] = [ a+b for a,b in zip(self.count[strand],other.count[strand]) ]
   #        for key in other.ambiguous[strand]:
   #            self.ambiguous[strand][key] += other.ambiguous[strand][key]
   #        for key in other.common[strand]:
   #            self.common[strand][key] += other.common[strand][key]



class Span_entry(object):
    """
        start
        end
        strand
        feature    
    """
    
    def __init__(self, start,end,strand,feature):
        self.start = start
        self.end = end
        self.strand = strand
        self.feature = feature
    
    def __repr__(self):
        return '%d:%d %d %s' % (self.start,self.end,self.strand,self.feature)



@config.help("""\
Count alignments to annotated features.

Working directories given should be the result of "shrimp:". It is not necessary to use samfilter or \
samconsensus. Alternatively a BAM file can be used, this should sorted in read-name order.

If a GENBANK file was given as reference to "shrimp:", this will be used as the reference, \
otherwise you need to specify an annotation file here.

By default, if a fragment aligns to multiple locations, only top scoring hits will be counted. \
If a fragment has several equal top alignments, it will be counted multiple times. \
Output contains columns detailing where this has happened. See the "--filter" flag for other options.
""")
@config.String_flag('filter', """\
Filtering mode:
poly     - Use all top scoring alignments
mono     - Use top scoring alignment if unique
existing - Use alignments from "filter:" or "consensus:"\
""")
@config.String_flag('types', 'Comma separated list of feature types to use.')
@config.String_flag('locii', 'Only use features with these locus_tags. Note: --types filter is still applied.')
@config.String_flag('strand', """\
How to use the relative strand of fragment and feature:
pool    - data is not single stranded, use both strands
forward - only count fragments on the same strand
reverse - only count fragments on the opposite strand
both    - output forward and reverse counts separately

This option MUST be specified explicitly.\
""")
@config.String_flag('qualifiers', 'Comma separated list of feature qualifiers to take from annotation file.') 
@config.Int_flag('min_score', 'Minimum fragment score. Note: 10 points equals 1 matched base.')
@config.Int_flag('min_size', 'Minimum fragment size.')
@config.Int_flag('max_size', 'Maximum fragment size.')
@config.Bool_flag('equalize', 'Subsample counts down to the size of the smallest sample.')
@config.Main_section('filenames', 'Working directories and/or BAM files and/or annotation files.')
class Count(config.Action_with_prefix):
    filter = 'poly'
    types = 'CDS'
    locii = None
    strand = None  #It is not good to have a default, this needs to be specified explicitly
    qualifiers = 'gene,product'
    min_score = None
    min_size = None
    max_size = None
    equalize = False
    filenames = [ ]

    def run(self):
        grace.require_samtools()
        count_run(
            min_score=self.min_score, min_size=self.min_size, max_size=self.max_size, 
            filter_mode=self.filter, equalize=self.equalize, types=self.types, locii=self.locii, 
            qualifiers=self.qualifiers, use_strand=self.strand, merge_filename=None, limit=None, 
            output_prefix=self.prefix, filenames=self.filenames, log=self.log
        )

#
#COUNT_HELP = """\
#
#Usage:
#
#    nesoni samcount: [options] output_prefix [working_dir ...] [bam_file.bam ...] [annotations.gbk/gff ...] 
#
#Count hits to CDS features.
#
#working_dirs should be the result of nesoni samshrimp. It is not necessary to use samfilter or 
#samconsensus. Or a BAM file can be used, this should sorted in read-name order.
#
#If a GENBANK file was given as reference to nesoni samshrimp, this will be used as the reference,
#otherwise you need to specify it here.
#
#If a fragment hits multiple locations, only top scoring hits will be counted.
#If a fragment has several equal top hits, it will be counted multiple times.
#Output contains columns detailing where this has happened.
#
#Options:
#
#    --min-score NNNN       - Minimum total fragment score
#                             Note 10 points == 1 matched base
#                             Default: 0
#
#    --min-size NNNN        - Minimum fragment size.
#                             Default: 0
#
#    --max-size NNNN        - Maximum fragment size.
#                             Default: 0
#
#    --filter poly/mono/existing
#                           - Filtering mode
#                             poly - Use all top scoring alignments
#                             mono - Use top scoring alignment if unique
#                             existing - Use alignments from "samfilter" or "samconsensus"
#                             Default: poly
#
#    --equalize yes/no      - Subsample counts down to the size of 
#                             the smallest sample.
#                             Default: no
#
#    --types type1,type2,...
#                      Comma separated list of feature types to use.
#                      Not case sensitive.
#                      Default: CDS
#
#    --locii locus_tag,locus_tag,...
#                      Only use features with these locus_tags.
#                      Note: --types filter is still applied
#                      Default: all features
#
#    --strand pool/forward/reverse/both
#                      How to use relative strand of read and feature:
#                        pool    - data is not single stranded, use both strands
#                        forward - only output count of forward strand
#                        reverse - only output count of reverse strand
#                        both    - output forward and reverse separately                      
#                      Default: both
#
#    --qualifiers qualifier,qualifier,...
#                      Comma separated list of feature qualifiers to take 
#                      from Genbank file.
#                      Default: gene,product
#
#"""
#    #--merge filename
#    #                  Merge features (for example splice variants).                      
#    #                  File should contain lines consisting of:
#    #                  <merged-name> <feature-name> [<feature-name>...]
#    #
#
#
#def count_main(args):
#    grace.require_samtools()
#
#    min_score, args = grace.get_option_value(args,'--min-score', int, 0)
#    min_size, args = grace.get_option_value(args,'--min-size', int, 0)
#    max_size, args = grace.get_option_value(args,'--max-size', int, None)
#    #monogamous, args = grace.get_option_value(args,'--monogamous', grace.as_bool, False)
#
#    filter_mode, args = grace.get_option_value(args,'--filter', str, 'poly')
#
#    equalize, args = grace.get_option_value(args,'--equalize', grace.as_bool, False)
#
#    types, args = grace.get_option_value(args, '--types', str, 'CDS')
#    locii, args = grace.get_option_value(args, '--locii', str, '')
#
#    qualifiers, args = grace.get_option_value(args, '--qualifiers', str, 'gene,product')
#    
#    use_strand, args = grace.get_option_value(args, '--strand', str, 'both')
#    
#    merge_filename, args = grace.get_option_value(args, '--merge', str, None)
#
#    limit, args = grace.get_option_value(args, '--limit', int, None)
#
#    grace.expect_no_further_options(args)
#
#    output_prefix = args[0]    
#    filenames = args[1:]
#
#    if len(args) < 2:
#        print >> sys.stderr, COUNT_HELP
#        raise grace.Help_shown()
#
#    count_run(
#        min_score, min_size, max_size, filter_mode, equalize, types, locii,
#        qualifiers, use_strand, merge_filename, limit, output_prefix, filenames
#    )

def tab_encode(listing):
    return '\t'.join( item.replace('\t','        ').replace('\n',' ') for item in listing )

def count_encode(counts):
    return ' '.join(
        '%dx%s' % (item[1],str(item[0]).replace(' ','-'))
        for item in sorted(counts.items(), key=lambda item:item[1], reverse=True)
        ) 

def count_parse(count_str):
    result = { }
    if count_str:
        for item in count_str.split(' '):
            x = item.index('x')
            result[item[x+1:]] = int(item[:x])
    return result


def count_run(
    min_score, min_size, max_size, filter_mode, equalize, types, locii,
    qualifiers, use_strand, merge_filename, limit, output_prefix, filenames, log):
    
    if filter_mode == 'poly':
        use_bam_filename = 'alignments.bam'
        use_only_top = True
        use_only_monogamous = False
        expect_multiple_alignments = True
    elif filter_mode == 'mono': 
        use_bam_filename = 'alignments.bam'
        use_only_top = True
        use_only_monogamous = True
        expect_multiple_alignments = True
    else:
        assert filter_mode == 'existing', 'Unrecognized filtering mode'
        use_bam_filename = 'alignments_filtered.bam'
        use_only_top = False
        use_only_monogamous = False
        expect_multiple_alignments = False

    types = types.lower().split(',')

    qualifiers = qualifiers.split(',')

    if locii:
        locii = locii.lower().split(',')
    else:
        locii = None

    assert use_strand is not None, 'You must now explicitly specify --strand'
    assert use_strand in ('pool','forward','reverse','both'), "Can't understand --strand specification."
    
    from Bio import Seq, SeqIO

    annotation_filenames = [ ]
    bam_filenames = [ ]
    for arg in filenames:
        if annotation.is_annotation_file(arg):
            annotation_filenames.append(arg)
        else:
            bam_filenames.append(arg)    

    n_samples = len(bam_filenames)
    titles = bam_filenames[:]
    tags = [ ]
    for i in xrange(len(bam_filenames)):
        if os.path.isdir(bam_filenames[i]):
            working = working_directory.Working(bam_filenames[i])            
            titles[i] = working.name
            tags.append(working.get_tags())
            if not annotation_filenames:
                reference_filename = working.get_reference().annotations_filename()
                if reference_filename is not None:
                    annotation_filenames.append(reference_filename)
            
            bam_filenames[i] = os.path.join(bam_filenames[i], use_bam_filename)
    
    assert bam_filenames, 'No reference alignments given' 

    merge = { }
    merge_qualifiers = { }
    if merge_filename is not None:
        #First line gives qualifiers
        #remaining lines give <qualifier> <qualifier...> <gene> <transcript> <transcript...>
    
        f = open(merge_filename,'rU')
        qualifiers = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if not parts: continue
            for name in parts[len(qualifiers)+1:]:
                assert name not in merge, 'Duplicate feature name in merge file'
                merge[name] = parts[len(qualifiers)]
                merge_qualifiers[name] = parts[:len(qualifiers)]
        f.close()

    
    genes = { }  # reference name -> gene index
    
    feature_names = { }   # feature_name -> number of occurrences
    
    features = [ ]
    
    n_features = 0

    chromosome_length = { }
    for filename in bam_filenames:
        headers = sam.bam_headers(filename)
        for line in headers.split('\n'):
            if not line: continue
            parts = line.split('\t')
            if parts[0] != '@SQ': continue
            
            name = None
            length = None
            for part in parts[1:]:
                if part.startswith('SN:'): name = part[3:]
                if part.startswith('LN:'): length = int(part[3:])
            assert name is not None and length is not None
            
            if name in chromosome_length:
                assert chromosome_length[name] == length
            else:
                chromosome_length[name] = length
    
    for name in chromosome_length:
        genes[name] = span_index.Span_index()

    
    if annotation_filenames:
        assert not merge, 'Merging not supported with annotation files'
    
        for filename in annotation_filenames:
            for feature in annotation.read_annotations(filename):
                if feature.type.lower() not in types: continue
                
                if (locii is not None and
                    ('locus_tag' not in feature.attr or
                     feature.attr['locus_tag'].lower() not in locii)):
                    continue
                
                f = Feature(n_samples)
                f.name = feature.get_id()                
                if feature.type.lower() != 'cds' and len(types) > 1:
                    f.name = feature.type + ':' + f.name

                feature_names[f.name] = feature_names.get(f.name,0)+1
                if feature_names[f.name] > 1:
                   f.name += '/%d' % feature_names[f.name]
                   
                f.qualifiers = [ feature.attr.get(item,'') for item in qualifiers ]
                
                f.length = feature.end - feature.start

                assert feature.seqid in genes, 'Annotation for sequence that is not in BAM files'
                genes[feature.seqid].insert(Span_entry(feature.start, feature.end, feature.strand or 1, f))
                features.append(f)

    else:
        # Sequences as features        
        log.log('No annotation files given or found, using sequences as features\n')
        
        name_feature = { } # (merged)name -> feature
        
        for name in chromosome_length:
            merged_name = merge.get(name, name)
            
            if merged_name not in name_feature: 
                f = Feature(n_samples)
                f.name = merged_name
                f.length = length
                f.qualifiers = merge_qualifiers.get(name, ('',)*len(qualifiers))
                n_features += 1
                name_feature[merged_name] = f
                features.append(f)
            else:
                f = name_feature[merged_name]
                f.length = max(f.length, length) #...
            
            genes[name].insert(Span_entry(0, chromosome_length[name], 1, f))
                
                

    
    log.log('%d features\n\n' % len(features))

    for name in genes: 
        genes[name].prepare()
        
    n_fragments = [ 0 ] * n_samples
    n_fragments_aligned = [ 0 ] * n_samples
    n_low_score = [ 0 ] * n_samples
    n_something = [ 0 ] * n_samples
    n_multiple = [ 0 ] * n_samples
    n_span = [ 0 ] * n_samples

    for i in xrange(n_samples):
        for read_name, fragment_alignments, unmapped in sam.bam_iter_fragments(bam_filenames[i], 'Counting sample %d of %d' % (i+1,n_samples)):
            n_fragments[i] += 1
            
            if not fragment_alignments: 
                continue
            
            n_fragments_aligned[i] += 1
            
            feature_hits = [ ] # [ [ (feature, strand) ] ]
            
            # Use only top scoring alignments
            fragment_scores = [ sum( al.get_AS() for al in item ) for item in fragment_alignments ]
            
            best_score = max(fragment_scores)
            
            if min_score is not None and best_score < min_score:
                n_low_score[i] += 1
                continue
            
            if use_only_top:
                cutoff = max(best_score, min_score)
            else:
                cutoff = min_score            
            fragment_alignments = [ item 
                                    for item, score in zip(fragment_alignments, fragment_scores)
                                    if score >= cutoff ]            
            
            for alignments in fragment_alignments:
                strand = -1 if alignments[0].flag&sam.FLAG_REVERSE else 1
            
                start = min(item.pos-1 for item in alignments)
                end = max(item.pos+item.length-1 for item in alignments)
                length = end-start
                if min_size is not None and length < min_size: continue
                if max_size is not None and length > max_size: continue
                
                rname = alignments[0].rname
                strand = -1 if alignments[0].flag&sam.FLAG_REVERSE else 1
                assert alignments[0].rname in genes, 'Alignment refers to sequence not present in GENBANK file'
            
                this_feature_hits = [ ]    
                for item in genes[rname].get(start, end):
                    rel_strand = strand * item.strand
                    key = (item.feature, rel_strand)
                    if key in this_feature_hits: continue
                    this_feature_hits.append( key )
                    if not use_only_monogamous or len(fragment_alignments) == 1:
                        item.feature.count[rel_strand][i] += 1
                
                if this_feature_hits: 
                    feature_hits.append( this_feature_hits )
                
                if len(this_feature_hits) > 1:
                    for a in this_feature_hits:
                        for b in this_feature_hits:
                            if a[0] is b[0]: continue
                            a[0].common[(a[1],b[1])][b[0]] += 1
                                
                
            if len(feature_hits) > 0: 
                n_something[i] += 1
            #else:
            #    print fragment_alignments
            #    print genes[fragment_alignments[0][0].rname].indexes
            #    print
            
            
            if len(feature_hits) > 1: 
                n_multiple[i] += 1
                for j in xrange(len(feature_hits)):
                    for k in xrange(len(feature_hits)):
                        if j == k: continue
                        for a in feature_hits[j]:
                            for b in feature_hits[k]:
                                if a[0] is b[0]: continue
                                a[0].ambiguous[(a[1],b[1])][b[0]] += 1
            
            if any(len(item) > 1 for item in feature_hits): n_span[i] += 1
            
            if limit is not None and n_fragments[i] >= limit: break

        grace.status('')

        #log.log('%s\n' % titles[i])
        #log.log('%20s fragments\n' % grace.pretty_number(n_fragments[i]))
        #log.log('%20s fragments aligned to the reference\n' % grace.pretty_number(n_fragments_aligned[i]))
        #if n_low_score[i]:
        #    log.log('%20s had too low an alignment score, discarded\n' % grace.pretty_number(n_low_score[i]))
        #log.log('%20s aligned to an annotated gene\n' % grace.pretty_number(n_something[i]))
        #if expect_multiple_alignments or n_multiple[i]:
        #    log.log('%20s aligned to multiple genes\n' % grace.pretty_number(n_multiple[i]))
        #log.log('%20s had an alignment that spanned multiple genes\n' % grace.pretty_number(n_span[i]))
        #log.log('\n')

        log.datum(titles[i], 'fragments', n_fragments[i])
        log.datum(titles[i], 'fragments aligned to the reference', n_fragments_aligned[i])
        if n_low_score[i]:
            log.datum(titles[i], 'had too low an alignment score, discarded', n_low_score[i])
        log.datum(titles[i], 'aligned to an annotated gene', n_something[i])
        if expect_multiple_alignments or n_multiple[i]:
            log.datum(titles[i], 'aligned to multiple genes', n_multiple[i])
        log.datum(titles[i],'had an alignment that spanned multiple genes', n_span[i])
        log.log('\n')

    strandedness = [ ]
    for feature in features:
        n_forward = sum(feature.count[1])
        n_reverse = sum(feature.count[-1])
        if n_forward+n_reverse < 5: continue
        strandedness.append( (n_forward-n_reverse)*100.0 / (n_forward+n_reverse) )
    strandedness = sum(strandedness) / max(1,len(strandedness))
    log.log('Strand specificity score: %.0f\n'
            '  (~ -100 reverse strand, ~ 0 non-specific, ~ 100 forward strand\n'
            '   Average over all features with at least 5 hits.)\n'
            % strandedness)


    if use_strand == 'pool':
        getters = [ lambda f: (feature.name, 
                               add_lists(feature.count[1],feature.count[-1]),
                               add_defdicts(feature.common[(1,1)], 
                                            feature.common[(1,-1)], 
                                            feature.common[(-1,1)], 
                                            feature.common[(-1,-1)]),
                               add_defdicts(feature.ambiguous[(1,1)], 
                                            feature.ambiguous[(1,-1)], 
                                            feature.ambiguous[(-1,1)], 
                                            feature.ambiguous[(-1,-1)])) ]
    elif use_strand == 'forward':
        getters = [ lambda f: (feature.name, feature.count[1], feature.common[(1,1)], feature.ambiguous[(1,1)]) ]
    elif use_strand == 'reverse':
        getters = [ lambda f: (feature.name, feature.count[-1], feature.common[(-1,-1)], feature.ambiguous[(-1,-1)]) ]
    elif use_strand == 'both':
        getters = [ lambda f: (feature.name, feature.count[1], feature.common[(1,1)], feature.ambiguous[(1,1)]),
                    lambda f: (feature.name + 'r', feature.count[-1], feature.common[(-1,-1)], feature.ambiguous[(-1,-1)]) ]

    total_hits = [0] * n_samples
    for feature in features:
        for getter in getters:
            total_hits = add_lists(total_hits, getter(feature)[1])

    if equalize:
        min_hits = min(total_hits)
        p = [ float(min_hits)/item for item in total_hits ]
        total_hits = [ min_hits ] * n_samples

    
    
    comments = [ '#Counts' ] + [
        '#sampleTags='+','.join(item)
        for item in tags
        ]
    
    names = [ ]
    
    count_type = io.named_list_type(titles)    
    counts = [ ]
    
    #rpkm_type = io.named_list_type(titles)
    #rpkms = [ ]
    
    annotation_type = io.named_list_type([ 'Length' ] + qualifiers)
    annotations = [ ]
    
    alignment_type = io.named_list_type(
        [ 'On same fragment' ] + 
            [ 'Ambiguous alignment' ] 
                if expect_multiple_alignments 
                else [ ]
        )
    alignments = [ ]
    
    for feature in features:
        for getter in getters:
            feature_name, count, common, ambiguous = getter(feature)
            
            if equalize:
                count = [
                    subsample(count[i], p[i])
                    for i in xrange(n_samples)
                    ]
            
            #rpkm = [ count[i] * 1.0e9 / feature.length / total_hits[i] for i in xrange(n_samples) ]
            
            #common_str = ' '.join(
            #    '%dx%s' % (item[1],item[0])
            #    for item in sorted(common.items(), key=lambda item:item[1], reverse=True)
            #    ) 
            #ambiguous_str = ' '.join(
            #    '%dx%s' % (item[1],item[0])
            #    for item in sorted(ambiguous.items(), key=lambda item:item[1], reverse=True)
            #    )
            common_str = count_encode(common)
            ambiguous_str = count_encode(ambiguous)
            
            names.append(feature_name)
            counts.append(count_type(count))
            #rpkms.append(rpkm_type(rpkm))
            annotations.append(annotation_type([ str(feature.length) ] + list(feature.qualifiers)))
            alignments.append(alignment_type([ common_str ] + [ ambiguous_str ] if expect_multiple_alignments else [ ]))

    groups = [
        ('Count', io.named_list_type(names,count_type)(counts)),
        #('RPKM', io.named_list_type(names,rpkm_type)(rpkms)),
        ('Annotation', io.named_list_type(names,annotation_type)(annotations)),
        ('Alignment', io.named_list_type(names,alignment_type)(alignments)),
        ]
    
    io.write_grouped_csv(output_prefix + '.csv', groups, rowname_name='Feature', comments=comments)

#
#
#
#
#    f = open(output_prefix + '.txt', 'wb')
#    #log.attach(open(output_prefix + '_log.txt', 'wb'))
#
#    print >> f, tab_encode(
#        [ 'Feature' ] +
#        titles +
#        [ 'RPKM ' + item for item in titles ] +
#        [ 'Length' ] +
#        qualifiers +
#        [ 'On same fragment' ] +
#        ([ 'Ambiguous alignment' ] if expect_multiple_alignments else [ ])
#    )
#    
#    for feature in features:
#        for getter in getters:
#            feature_name, count, common, ambiguous = getter(feature)
#            
#            if equalize:
#                count = [
#                    subsample(count[i], p[i])
#                    for i in xrange(n_samples)
#                ]
#            
#            rpkm = [ count[i] * 1.0e9 / feature.length / total_hits[i] for i in xrange(n_samples) ]
#            
#            common_str = ' '.join(
#                '%dx%s' % (item[1],item[0])
#                for item in sorted(common.items(), key=lambda item:item[1], reverse=True)
#            ) 
#            ambiguous_str = ' '.join(
#                '%dx%s' % (item[1],item[0])
#                for item in sorted(ambiguous.items(), key=lambda item:item[1], reverse=True)
#            ) 
#            
#            print >> f, tab_encode(
#                [ feature_name ] +
#                [ str(item) for item in count ] +
#                [ '%.2f' % item for item in rpkm ] +
#                [ str(feature.length) ] +
#                list(feature.qualifiers) +
#                [ common_str ] +
#                ([ ambiguous_str ] if expect_multiple_alignments else [ ]) 
#            )
#            
#
#    f.close()




###class frag_histogram(collections.defaultdict):
###    def __init__(self):
###        collections.defaultdict.__init__(self, lambda: 0)
###
###    def total(self):
###        return sum(self.values())
###
###    def add(self, other):
###        for item in other:
###            self[item] += other[item]
###
###    def __add__(self, other):
###        result = frag_histogram()
###        result.add(self)
###        result.add(other)
###        return result
###
###def count_main(args):
###    grace.require_samtools()
###    
###    log = grace.Log()
###
###    min_score, args = grace.get_option_value(args,'--min-score', int, 0)
###    monogamous, args = grace.get_option_value(args,'--monogamous', grace.as_bool, False)
###    fraggalize, args = grace.get_option_value(args,'--fraggalize', grace.as_bool, False)
###
###    types, args = grace.get_option_value(args, '--types', str, 'CDS')
###    types = types.lower().split(',')
###    locii, args = grace.get_option_value(args, '--locii', str, '')
###    if locii:
###        locii = locii.lower().split(',')
###    else:
###        locii = None
###    
###    qualifiers, args = grace.get_option_value(args, '--qualifiers', str, 'gene,product')
###    qualifiers = qualifiers.split(',')
###    
###    use_strand, args = grace.get_option_value(args, '--strand', str, 'both')
###    
###    merge_filename, args = grace.get_option_value(args, '--merge', str, None)
###
###    limit, args = grace.get_option_value(args, '--limit', int, None)
###
###    assert use_strand in ('pool','forward','reverse','both'), "Can't understand --strand specification."
###
###    if len(args) < 2:
###        print >> sys.stderr, COUNT_HELP
###        raise grace.Help_shown()
###    
###    from Bio import Seq, SeqIO
###
###    output_prefix = args[0]    
###
###    genbank_filenames = [ ]
###    bam_filenames = [ ]
###    for arg in args[1:]:
###        if os.path.isdir(arg):
###            peek = ''
###        else:
###            f = io.open_possibly_compressed_file(arg)
###            peek = f.read(5)
###            f.close()
###        if peek == 'LOCUS':
###            genbank_filenames.append(arg)
###        else:
###            bam_filenames.append(arg)    
###
###    n_samples = len(bam_filenames)
###    titles = bam_filenames[:]    
###    for i in xrange(len(bam_filenames)):
###        if os.path.isdir(bam_filenames[i]):
###            titles[i] = os.path.basename(bam_filenames[i])
###            if not genbank_filenames:
###                reference_filename = os.path.join(bam_filenames[i], 'reference.gbk')
###                if os.path.exists(reference_filename):
###                    genbank_filenames.append(reference_filename)
###            
###            bam_filenames[i] = os.path.join(bam_filenames[i], 'alignments.bam')
###    
###    #assert genbank_filenames, 'No reference GENBANK files given or found in reference alignment directories.'
###    
###    assert bam_filenames, 'No reference alignments given' 
###
###    merge = { }
###    if merge_filename is not None:
###        for line in open(merge_filename,'rU'):
###            parts = line.strip().split()
###            if not parts: continue
###            for name in parts[1:]:
###                assert name not in merge, 'Duplicate feature name in merge file'
###                merge[name] = parts[0]
###
###    
###    genes = { }  # reference name -> gene index
###    chromosome_names = [ ]
###    
###    feature_names = { }   # feature_name -> number of occurrences
###    
###    features = [ ]
###    
###    n_features = 0
###    
###    if genbank_filenames:
###        assert not merge, 'Merging not supported with GENBANK files'
###    
###        for genbank_filename in genbank_filenames:
###          for record in SeqIO.parse(io.open_possibly_compressed_file(genbank_filename),'genbank'):
###            name = record.id
###            seq = record.seq.tostring()
###            assert name not in genes, 'Duplicate record name in genbank file'
###            chromosome_names.append(name)
###            genes[name] = span_index.Span_index()
###        
###            for feature in record.features:
###                if feature.type.lower() not in types: continue
###                
###                if locii is not None:
###                    if 'locus_tag' not in feature.qualifiers or \
###                       feature.qualifiers['locus_tag'][0].lower() not in locii:
###                        continue 
###                
###                f = Feature(n_samples)
###                if 'locus_tag' not in feature.qualifiers:
###                    f.name = '%d..%d' % (feature.location.nofuzzy_start+1,feature.location.nofuzzy_end)
###                else:
###                    f.name = feature.qualifiers['locus_tag'][0]
###                    
###                if feature.type.lower() != 'cds':
###                    f.name = feature.type + ':' + f.name
###                
###                feature_names[f.name] = feature_names.get(f.name,0)+1
###                if feature_names[f.name] > 1:
###                   f.name += '/%d' % feature_names[f.name]
###                
###                #f.gene = feature.qualifiers.get('gene',[''])[0].strip()
###                #f.product = feature.qualifiers.get('product',[''])[0].strip()
###                f.feature = feature
###                
###                start = feature.location.nofuzzy_start
###                end = feature.location.nofuzzy_end
###                strand = feature.strand
###                if strand is None:
###                    strand = 1
###                
###                f.length = end-start
###                
###                #f.sequence = record.seq[start:end]
###                #if strand < 0:
###                #   f.sequence = bio.reverse_complement(f)
###                    
###                genes[name].insert(start, end, strand, f)
###                features.append(f)
###
###    else:
###        # Sequences as features
###        qualifiers = [ ]
###        
###        log.log('No GENBANK files given or found, using sequences as features\n')
###        
###        name_feature = { } # (merged)name -> feature
###        
###        for filename in bam_filenames:
###            headers = bam_headers(filename)
###            for line in headers.split('\n'):
###                if not line: continue
###                parts = line.split('\t')
###                if parts[0] != '@SQ': continue
###                
###                name = None
###                length = None
###                for part in parts[1:]:
###                    if part.startswith('SN:'): name = part[3:]
###                    if part.startswith('LN:'): length = int(part[3:])
###                assert name is not None and length is not None
###                
###                if name in genes: continue
###                genes[name] = span_index.Span_index()
###                chromosome_names.append(name)
###                
###                merged_name = merge.get(name, name)
###                
###                if merged_name not in name_feature: 
###                    f = Feature(n_samples)
###                    f.name = merged_name
###                    f.length = length
###                    n_features += 1
###                    name_feature[merged_name] = f
###                    features.append(f)
###                else:
###                    f = name_feature[merged_name]
###                    f.length = max(f.length, length) #...
###                
###                #f.start = 0
###                #f.end = length
###                #f.strand = 1                
###                genes[name].insert(0, length, 1, f)
###                
###                
###
###    
###    log.log('%d features\n\n' % len(features))
###
###    for name in genes: 
###        genes[name].prepare()
###        
###    n_fragments = [ 0 ] * n_samples
###    n_fragments_aligned = [ 0 ] * n_samples
###    n_low_score = [ 0 ] * n_samples
###    n_something = [ 0 ] * n_samples
###    n_multiple = [ 0 ] * n_samples
###    n_span = [ 0 ] * n_samples
###
###    for i in xrange(n_samples):
###        for read_name, fragment_alignments, unmapped in sam.bam_iter_fragments(bam_filenames[i], 'Counting sample %d of %d' % (i+1,n_samples)):
###            n_fragments[i] += 1
###            
###            if not fragment_alignments: 
###                continue
###            
###            n_fragments_aligned[i] += 1
###            
###            feature_hits = [ ] # [ [ (feature, strand) ] ]
###            
###            # Use only top scoring alignments
###            fragment_scores = [ sum( al.get_AS() for al in item ) for item in fragment_alignments ]
###            best_score = max(fragment_scores)
###            
###            if best_score < min_score:
###                n_low_score[i] += 1
###                continue
###            
###            fragment_alignments = [ item 
###                                    for item, score in zip(fragment_alignments, fragment_scores)
###                                    if score >= best_score ]            
###            
###            for alignments in fragment_alignments:
###                strand = -1 if alignments[0].flag&sam.FLAG_REVERSE else 1
###            
###                start = min(item.pos-1 for item in alignments)
###                end = max(item.pos+item.length-1 for item in alignments)
###                length = end-start
###                rname = alignments[0].rname
###                strand = -1 if alignments[0].flag&sam.FLAG_REVERSE else 1
###                assert alignments[0].rname in genes, 'Alignment refers to sequence not present in GENBANK file'
###            
###                this_feature_hits = [ ]    
###                for item in genes[rname].get(start, end):
###                    rel_strand = strand * item.strand
###                    key = (item.feature, rel_strand)
###                    if key in this_feature_hits: continue
###                    this_feature_hits.append( key )
###                    if not monogamous or len(fragment_alignments) == 1:
###                        item.feature.count[rel_strand][i][length] += 1
###                
###                if this_feature_hits: 
###                    feature_hits.append( this_feature_hits )
###                
###                if len(this_feature_hits) > 1:
###                    for a in this_feature_hits:
###                        for b in this_feature_hits:
###                            if a[0] is b[0]: continue
###                            a[0].common[a[1]][b[0]] += 1
###                                
###                
###            if len(feature_hits) > 0: 
###                n_something[i] += 1
###            #else:
###            #    print fragment_alignments
###            #    print genes[fragment_alignments[0][0].rname].indexes
###            #    print
###            
###            
###            if len(feature_hits) > 1: 
###                n_multiple[i] += 1
###                for j in xrange(len(feature_hits)):
###                    for k in xrange(len(feature_hits)):
###                        if j == k: continue
###                        for a in feature_hits[j]:
###                            for b in feature_hits[k]:
###                                if a[0] is b[0]: continue
###                                a[0].ambiguous[a[1]][b[0]] += 1
###            
###            if any(len(item) > 1 for item in feature_hits): n_span[i] += 1
###            
###            if limit is not None and n_fragments[i] >= limit: break
###
###        grace.status('')
###
###        log.log('%s\n' % titles[i])
###        log.log('%20s fragments\n' % grace.pretty_number(n_fragments[i]))
###        log.log('%20s fragments aligned to the reference\n' % grace.pretty_number(n_fragments_aligned[i]))
###        if n_low_score[i]:
###            log.log('%20s had too low an alignment score, discarded\n' % grace.pretty_number(n_low_score[i]))
###        log.log('%20s aligned to an annotated gene\n' % grace.pretty_number(n_something[i]))
###        log.log('%20s aligned to multiple genes\n' % grace.pretty_number(n_multiple[i]))
###        log.log('%20s had an alignment that spanned multiple genes\n' % grace.pretty_number(n_span[i]))
###        log.log('\n')
###
###    strandedness = [ ]
###    for feature in features:
###        n_forward = sum(feature.count[1], start=frag_histogram()).total()
###        n_reverse = sum(feature.count[-1], start=frag_histogram()).total()
###        if n_forward+n_reverse < 5: continue
###        strandedness.append( (n_forward-n_reverse)*100.0 / (n_forward+n_reverse) )
###    strandedness = sum(strandedness) / len(strandedness)
###    log.log('Strand specificity: %.0f%%\n'
###            '  (~ -100%% reverse strand, ~ 0%% non-specific, ~ 100%% forward strand\n'
###            '   Average over all features with at least 5 hits.)\n'
###            % strandedness)
###
###
###    if use_strand == 'pool':
###        getters = [ lambda f: (feature.name, 
###                               add_lists(feature.count[1],feature.count[-1]),
###                               add_defdicts(feature.common[1], feature.common[-1]),
###                               add_defdicts(feature.ambiguous[1], feature.ambiguous[-1])) ]
###    elif use_strand == 'forward':
###        getters = [ lambda f: (feature.name, feature.count[1], feature.common[1], feature.ambiguous[1]) ]
###    elif use_strand == 'reverse':
###        getters = [ lambda f: (feature.name, feature.count[-1], feature.common[-1], feature.ambiguous[-1]) ]
###    elif use_strand == 'both':
###        getters = [ lambda f: (feature.name, feature.count[1], feature.common[1], feature.ambiguous[1]),
###                    lambda f: (feature.name + 'r', feature.count[-1], feature.common[-1], feature.ambiguous[-1]) ]
###
###    total_hist = [ frag_histogram() for i in xrange(n_samples) ]
###    for feature in features:
###        for getter in getters:
###            total_hist = add_lists(total_hist, getter(feature)[1])
###
###    total_hits = [ item.total() for item in total_hist ]
###
###    if fraggalize:
###        min_hist = frag_histogram()
###        for item in total_hist:
###            for i in item:
###                if i not in min_hist: 
###                    min_hist[i] = item[i]
###                else:
###                    min_hist[i] = min(item[i], min_hist[i])
###        
###        print min_hist
###        
###        scaling = [ frag_histogram() for i in xrange(n_samples) ]
###        for i in xrange(n_samples):
###            for j in total_hist[i]:
###                scaling[i][j] = float(min_hist[j])/total_hist[i][j]
###
###        total_hits = [ min_hist.total() ] * n_samples
###
###
###    f = open(output_prefix + '.txt', 'wb')
###    log.attach(open(output_prefix + '_log.txt', 'wb'))
###
###    print >> f, 'Feature\t%s\t%s\tLength%s\tOn same fragment\tAmbiguous alignment' % (
###        '\t'.join( titles ),
###        '\t'.join( 'RPKM ' + item for item in titles ),
###        ''.join( '\t'+item for item in qualifiers ),
###    )
###    
###    for feature in features:
###        for getter in getters:
###            feature_name, count, common, ambiguous = getter(feature)
###            
###            if fraggalize:
###               scaled_count = [ frag_histogram() for i in xrange(n_samples) ]
###               for i in xrange(n_samples):
###                   for j in count[i]:
###                       scaled_count[i][j] = count[i][j] * scaling[i][j]
###               count = scaled_count 
###            
###            count_str = '\t'.join( str(item.total()) for item in count )
###            
###            rpkm = [ count[i].total() * 1.0e9 / feature.length / total_hits[i] for i in xrange(n_samples) ]
###            rpkm_str = '\t'.join( '%.2f' % item for item in rpkm )
###            
###            common_str = ' '.join(
###                '%dx%s' % (item[1],item[0])
###                for item in sorted(common.items(), key=lambda item:item[1], reverse=True)
###            ) 
###            ambiguous_str = ' '.join(
###                '%dx%s' % (item[1],item[0])
###                for item in sorted(ambiguous.items(), key=lambda item:item[1], reverse=True)
###            ) 
###            
###            this_qualifiers = ''.join([
###                '\t' + ' '.join(feature.feature.qualifiers.get(item,[''])).strip()
###                for item in qualifiers
###            ])
###            
###            print >> f, '%s\t%s\t%s\t%d%s\t%s\t%s' % (
###                feature_name, count_str, rpkm_str,
###                feature.length, this_qualifiers,
###                common_str, ambiguous_str
###            )
###
###    f.close()
###

def merge_require_equal(a,b):
    assert a == b, 'Values aren\'t equal while merging matrices'
    return a

def merge_counts(a,b):
    counts = count_parse(a)
    for key, value in count_parse(b).items():
        counts[key] = counts.get(key,0) + value
    return count_encode(counts)
    

def matrix_merge(matrices, reducer=merge_require_equal):
    """ Merge several matrices. 
        Matrices should have the same row names.
        Where multiple matrices have the same column, the values are merged using the reducer function.
        """
    columns = collections.OrderedDict()
    row_names = matrices[0].keys()
    for i, matrix in enumerate(matrices):
        assert row_names == matrix.keys(), 'Row names don\'t match in matrix merge'
        for column in matrix.value_type().keys():
            if column not in columns:
                columns[column] = [ ]
            columns[column].append( i )
    
    return io.named_matrix_type(row_names, columns.keys())([
        [ reduce(reducer, [ matrices[i][name][column] for i in columns[column] ])
            for column in columns
            ]
        for name in row_names
        ])
    

@config.help(
    'Merge several count files as produced by "count:".'
    )
@config.Main_section('filenames', 'CSV files to merge.')
class Merge_counts(config.Action_with_prefix):
    filenames = [ ]
    
    def run(self):
        assert self.filenames, 'No files given to merge.'
        
        tables = [ ]
        for filename in self.filenames:
            tables.append(io.read_grouped_table(
                filename,
                [('Count',str), ('Annotation',str), ('Alignment',str)],
                'Count',
                ))
                
        result = io.Grouped_table()
        result.comments = [ '#Counts' ]
        for table in tables:
            for comment in table.comments:
                if comment != '#Counts':
                    result.comments.append(comment)
        
        result['Count'] = matrix_merge([ table['Count'] for table in tables ])
        result['Annotation'] = matrix_merge([ table['Annotation'] for table in tables ])
        result['Alignment'] = matrix_merge([ table['Alignment'] for table in tables ], merge_counts)        
        result.write_csv(self.prefix + '.csv')














