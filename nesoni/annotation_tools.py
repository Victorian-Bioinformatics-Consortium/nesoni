
import collections

from nesoni import config, annotation, span_index

def join_descriptions(seq):
    result = [ ]
    for item in seq:
        #parts = item.split('isoform')        
        #if len(parts) == 2 and len(parts[1].strip()) <= 1:
        #   item = parts[0].rstrip()
        if item not in result: 
            result.append(item)
    return ', '.join(result)


@config.help(
    'Modify annotated features.'
    '\n\n'
    'Note: Features without a strand will not be shifted.'
    )
@config.Int_flag('shift_start', 'Bases to shift feature start.')
@config.Int_flag('shift_end', 'Bases to shift feature end.')
@config.String_flag('type', 'Output feature type.\nDefault: retain existing type.')
@config.String_flag('select', 'What types of annotation to use (selection expression).')
@config.Main_section('filenames', 'Annotation files.',empty_is_ok=False)
class Modify_features(config.Action_with_prefix):
    type = None
    shift_start = 0
    shift_end = 0
    select = 'all'
    filenames = [ ]
 
    def run(self):
        out_file = open(self.prefix+'.gff','wb')    
        annotation.write_gff3_header(out_file)
        
        for filename in self.filenames:
            for item in annotation.read_annotations(self.filename):
                if not selection.matches(self.select, [item.type]): continue
                
                if item.strand == 1:
                    item.start += self.shift_start
                    item.end += self.shift_end
                elif item.strand == -1:
                    item.end -= self.shift_start
                    item.start -= self.shift_end
                
                if self.type:
                    item.type = self.type
            
                print >> out_file, item.as_gff()

        out_file.close()


@config.help(
    'Merge overlapping features.'
    )
@config.Int_flag(
    'overlap', 
    'Allowed overlap in bases.'
    ' Negative to merge features that aren\'t quite touching.'
    )
@config.String_flag('type', 'Output feature type.\nDefault: retain existing type.')
@config.String_flag('select', 'What types of annotation to use (selection expression).')
@config.Main_section('filenames', 'Annotation files.',empty_is_ok=False)
class Collapse_features(config.Action_with_prefix):
    overlap = 0
    type = None
    select = 'all'
    filenames = [ ]

    def run(self):
        annotations = [ ]
        for filename in self.filenames:
            for item in annotation.read_annotations(self.filename):
                if not selection.matches(self.select, [item.type]): continue
                if self.type:
                    item.type = self.type
                annotations.append(item)
        
        annotations.sort(key=lambda item: (item.seqid, item.strand, item.start))
        
        group = [ ]
        groups = [ ]
        def emit():
            if not group: return
            groups.append(group[:])
            del group[:]        
        seqid = None
        strand = None
        end = 0
        for item in annotations:
            if item.seqid != seqid or item.strand != strand or item.start >= end:
                emit()
                seqid = item.seqid
                strand = item.strand
                end = item.end-self.overlap
            group.append(item)
            end = max(item.end-self.overlap, end)
        emit()

        out_file = open(self.prefix+'.gff','wb')
        annotation.write_gff3_header(out_file)

        for group in groups:
            item = annotation.Annotation()
            item.source = group[0].source
            item.type = join_descriptions( item2.type for item2 in group )
            item.seqid = group[0].seqid
            item.strand = group[0].strand
            item.start = min( item2.start for item2 in group )
            item.end = max( item2.end for item2 in group )
            item.score = None
            item.phase = None
            item.attr = { }
            
            for item2 in group:
                for key in item2.attr:
                    if key in item.attr: continue
                    item.attr[key] = join_descriptions( item3.attr[key] for item3 in group if key in item3.attr )
            
            print >> out_file, item.as_gff()
            
        out_file.close()


class _Related_feature(collections.namedtuple(
        '_Related_feature', 
        'feature start end relations')):
    def modify_with_relations(self, parent_selection):
        buckets = collections.defaultdict(list)
        
        my_strand = self.feature.strand or 0
        for item in self.relations:
            their_strand = item.feature.strand or 0
            overlaps = self.feature.overlaps(item,check_strand=False)
            if my_strand * their_strand == -1:
                if overlaps:
                    relation = 'opposite'
                elif item.feature.start*my_strand < self.feature.start*my_strand:
                    relation = 'upstrand-opposite'
                else:
                    relation = 'downstrand-opposite'
            elif overlaps:
                relation = 'in'
            else:
                strand = my_strand or their_strand
                if not strand:
                    relation = 'near'                
                else:
                    if item.feature.start*strand < self.feature.start*strand:
                        relation = 'upstrand'
                    else:
                        relation = 'downstrand'
            
            buckets[relation].append(item)
                
        

@config.help(
    'Annotate features with their relation to a second set of features.'
    ' For example, relate transcribed regions to annotated genes.'
    '\n\n'
    'Although we refer to the feature sets as parent and child,'
    ' the relation is treated as many-to-many.'
    '\n\n'
    'Possible relations for unstranded features are:\n'
    'in near'
    '\n\n'
    'Possible relations for stranded features are:\n'
    'in opposite upstrand downstrand upstrand-opposite downstrand-opposite\n\n'
    'Features that overlap are also called "in" (or "opposite").'
    )
@config.Int_flag('upstrand', 'Number of bases upstrand of features in features-1 to look.')
@config.Int_flag('downstrand', 'Number of bases downstrand of features in features-1 to look.')
@config.String_flag('use', 'What relations to set "Parent" attribute for (selection expression).')
@config.String_flag('select_parent', 'What types of annotation to use from parent features file (selection expression).')
@config.String_flag('select_child', 'What types of annotation to use from child features file (selection expression).')
@config.Positional('parent', 'File containing "parent" features.')
@config.Positional('child', 'File containing "child" features.')
class Relate_features(config.Action_with_prefix):
    upstrand = 0
    downstrand = 0
    use = 'all'
    select_parent = 'all'
    select_child = 'all'
    parent = None
    child = None

    # Also output un-related features.
    
    def run(self):
        features_parent = [ 
            _Related_feature(item,item.start,item.end,[]) 
            for item in read_annotations(self.parent) 
            if selection.matches(self.select_parent, [item.type]) 
            ]
        features_child = [ 
            _Related_feature(item,item.start,item.end,[]) 
            for item in read_annotations(self.child) 
            if selection.matches(self.select_child, [item.type])
            ]
        
        index = Span_index()
        for item in features_child:
            index.insert(item)
        index.prepare()
        
        for item_1 in features_parent:
            if item_1.strand == 1:
                start = item_1.start - self.upstrand
                end = item_1.end + self.downstrand
            elif item_1.strand == -1:
                start = item_1.start - self.downstrand
                end = item_1.end + self.upstrand
            else:
                start = item_1.start - max(self.upstrand,self.downstrand)
                end = item_1.end + max(self.upstrand,self.downstrand)
            for item_2 in index.get(start,end):
                item_1.relations.append(item_2)
                item_2.relations.append(item_1)

        for item in features_parent:
            item.modify_with_relations(self.use)
        for item in features_child:
            item.modify_with_relations('-all')
        
        f_1 = open(self.prefix + '-parent.gff','wb')
        f_2 = open(self.prefix + '-child.gff','wb')




