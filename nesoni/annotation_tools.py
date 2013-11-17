
import collections, math

from nesoni import config, annotation, span_index, selection

def decode_shift(expression):
    absolute = 0.0
    proportion = 0.0
    i = 0
    while i < len(expression):
        j = i
        while j < len(expression) and expression[j] in '+-':
            j += 1
        sign_chunk = expression[i:j]        

        i = j 
        while j < len(expression) and expression[j] not in '+-':
            j += 1
        chunk = expression[i:j]
        i = j
        
        assert chunk, 'Couldn\'t parse shift'
        if sign_chunk in ('+',''):
            sign = 1
        elif sign_chunk == '-':
            sign = -1
        else:
            assert False, 'Couldn\'t parse shift'
            
        if chunk.endswith('%'):
            proportion += sign*float(chunk[:-1])/100.0
        else:
            absolute += sign*float(chunk)
    return absolute, proportion


def join_descriptions(seq, joiner=', '):
    result = [ ]
    for item in seq:
        #parts = item.split('isoform')        
        #if len(parts) == 2 and len(parts[1].strip()) <= 1:
        #   item = parts[0].rstrip()
        if item not in result: 
            result.append(item)
    return joiner.join(result)


STRAND_CHANGE = {
   'no'      : {-1:-1,0:0,1:1,None:None},
   'flip'    : {-1:1,0:0,1:-1,None:None},
   'clear'   : {-1:0,0:0,1:0,None:0},
   'forward' : {-1:1,0:1,1:1,None:1},
   'reverse' : {-1:-1,0:-1,1:-1,None:-1},
}


@config.help(
    'Modify annotated features.'
    '\n\n'
    'Note: Features without a strand will not be shifted.'
    )
@config.String_flag('shift_start', 
    'Bases to shift feature start. '
    'Can be absolute or a percentage, or a combination, eg '
    '5 50% 50%+20 10-100%')
@config.String_flag('shift_end', 'Bases to shift feature end. Format as for --shift-start.')
@config.String_flag(
    'change_strand',
    'no / flip / clear / forward / reverse\n'
    'What to do with the strand of features. '
    '(Applied after shifting the feature.)'
    )
@config.String_flag('type', 'Output feature type.\nDefault: retain existing type.')
@config.String_flag('rename', 'Rename attributes. Comma separated list of newname=oldname.')
@config.String_flag('select', 'What types of annotation to use (selection expression).')
@config.Main_section('filenames', 'Annotation files.',empty_is_ok=False)
class Modify_features(config.Action_with_prefix):
    type = None
    shift_start = '0'
    shift_end = '0'
    rename = ''
    change_strand = 'no'
    select = 'all'
    filenames = [ ]
 
    def run(self):
        assert self.change_strand in STRAND_CHANGE, 'Unknown way to change strand.'
        strand_changer = STRAND_CHANGE[self.change_strand]
        
        shift_start_absolute, shift_start_proportion = decode_shift(self.shift_start)
        shift_end_absolute, shift_end_proportion = decode_shift(self.shift_end)
        
        renames = [ ]
        if self.rename:
            for item in self.rename.split(','):
                new, old = item.split('=')
                if new != old:
                    renames.append((new,old))
    
        out_file = open(self.prefix+'.gff','wb')    
        annotation.write_gff3_header(out_file)
        
        for filename in self.filenames:
            for item in annotation.read_annotations(filename):
                if not selection.matches(self.select, [item.type]): continue
                
                if self.type:
                    item.type = self.type
                
                length = item.end-item.start
                shift_start = int(math.floor(0.5+shift_start_absolute+shift_start_proportion*length))
                shift_end = int(math.floor(0.5+shift_end_absolute+shift_end_proportion*length))
                
                if item.strand == 1:
                    item.start += shift_start
                    item.end += shift_end
                elif item.strand == -1:
                    item.end -= shift_start
                    item.start -= shift_end
                
                item.strand = strand_changer[item.strand]
                
                old_attr = item.attr.copy()
                for new,old in renames:
                    if old in item.attr:
                       del item.attr[old]
                for new,old in renames:
                    if old in old_attr:
                       item.attr[new] = old_attr[old]
            
                print >> out_file, item.as_gff()

        out_file.close()


@config.help(
    'Merge overlapping features. Only features of the same type are merged. "Parent" relationships are preserved.'
    )
@config.Int_flag(
    'overlap', 
    'Allowed overlap in bases.'
    ' Negative to merge features that aren\'t quite touching.'
    )
@config.String_flag('type', 'Output feature type.\nDefault: retain existing type.')
@config.String_flag('select', 'What types of annotation to use (selection expression).')
@config.String_flag('joiner', 'Separator for joined fields.')
@config.Main_section('filenames', 'Annotation files.',empty_is_ok=False)
class Collapse_features(config.Action_with_prefix):
    overlap = 0
    type = None
    select = 'all'
    joiner = '/'
    filenames = [ ]

    def run(self):
        annotations = [ ]
        for filename in self.filenames:
            for item in annotation.read_annotations(filename):
                if not selection.matches(self.select, [item.type]): continue
                if self.type:
                    item.type = self.type
                annotations.append(item)
        
        annotations.sort(key=lambda item: (item.type, item.seqid, item.strand, item.start))
        
        group = [ ]
        groups = [ ]
        def emit():
            if not group: return
            groups.append(group[:])
            del group[:]        
        type = None
        seqid = None
        strand = None
        end = 0
        for item in annotations:
            if item.type != type or item.seqid != seqid or item.strand != strand or item.start >= end:
                emit()
                type = item.type
                seqid = item.seqid
                strand = item.strand
                end = item.end-self.overlap
            group.append(item)
            end = max(item.end-self.overlap, end)
        emit()


        items = [ ]
        
        id_map = { }

        for group in groups:
            item = annotation.Annotation()
            item.source = group[0].source
            item.type = group[0].type
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
                    item.attr[key] = join_descriptions([ item3.attr[key] for item3 in group if key in item3.attr ], self.joiner )

            item.parents = [ ]
            for item2 in group:
                if 'ID' in item2.attr:
                    assert item2.attr['ID'] not in id_map, 'Duplicate ID: '+item2.attr['ID']
                    id_map[item2.attr['ID']] = item.attr['ID']
                if 'Parent' in item2.attr:
                    item.parents.append(item2.attr['Parent'])
            
            items.append(item)
        
        for item in items:
            if item.parents:
                item.attr['Parent'] = join_descriptions([ id_map.get(parent,parent) for parent in item.parents ], ',')
        
        with open(self.prefix+'.gff','wb') as out_file:
            annotation.write_gff3_header(out_file)
            for item in items:
                print >> out_file, item.as_gff()


class _Related_feature(collections.namedtuple(
        '_Related_feature', 
        'feature start end relations')):
    def __hash__(self):
        return id(self)
    
    def add_to_attr(self, name, value):
        if name in self.feature.attr:
            self.feature.attr[name] += ','+value
        else:
            self.feature.attr[name] = value
    
    def modify_with_relations(self, use, to_child, to_parent):
        buckets = collections.defaultdict(list)
        
        my_strand = self.feature.strand or 0
        for item in self.relations:
            their_strand = item.feature.strand or 0
            overlaps = self.feature.overlaps(item.feature,check_strand=False)
            if my_strand * their_strand == -1:
                if overlaps:
                    relation = 'opposite'
                elif item.feature.start*my_strand < self.feature.start*my_strand:
                    relation = 'upstrand_opposite'
                else:
                    relation = 'downstrand_opposite'
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
        
        for name,relatives in buckets.items():        
            if selection.matches(use, [name]):
                for relative in relatives:
                    self.add_to_attr('has_'+name, relative.feature.get_id())
                    relative.add_to_attr('is_'+name, self.feature.get_id())
                    relative.add_to_attr('Parent', self.feature.get_id())
                    
                    for key in self.feature.attr:
                        if selection.matches(to_child,[key]):
                            relative.add_to_attr(key, self.feature.attr[key])
                    for key in relative.feature.attr:
                        if selection.matches(to_parent,[key]):
                            self.add_to_attr(key, relative.feature.attr[key])
                
        

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
    'in opposite upstrand downstrand upstrand_opposite downstrand_opposite\n\n'
    'Features that overlap are also called "in" (or "opposite").'
    )
@config.Int_flag('upstrand', 'Number of bases upstrand of parent features to look.')
@config.Int_flag('downstrand', 'Number of bases downstrand of parent features to look.')
@config.String_flag('use', 'What relationship types to use (selection expression).')
@config.String_flag('to_child', 'What attributes to copy to children.')
@config.String_flag('to_parent', 'What attributes to copy to parent.')
@config.String_flag('select_parent', 'What types of annotation to use from parent features file (selection expression).')
@config.String_flag('select_child', 'What types of annotation to use from child features file (selection expression).')
@config.Positional('parent', 'File containing "parent" features.')
@config.Positional('child', 'File containing "child" features.')
class Relate_features(config.Action_with_prefix):
    upstrand = 0
    downstrand = 0
    use = 'all'
    to_child = '-all'
    to_parent = '-all'
    select_parent = 'all'
    select_child = 'all'
    parent = None
    child = None

    # TODO: Also output un-related features.
    
    def run(self):
        features_parent = [ 
            _Related_feature(item,item.start,item.end,[]) 
            for item in annotation.read_annotations(self.parent) 
            if selection.matches(self.select_parent, [item.type]) 
            ]
        features_child = [ 
            _Related_feature(item,item.start,item.end,[]) 
            for item in annotation.read_annotations(self.child) 
            if selection.matches(self.select_child, [item.type])
            ]
        
        index = { }
        for item in features_child:
            if item.feature.seqid not in index:
                index[item.feature.seqid] = span_index.Span_index()
            index[item.feature.seqid].insert(item)
        for value in index.values():
            value.prepare()
        
        for item_1 in features_parent:
            if item_1.feature.strand == 1:
                start = item_1.start - self.upstrand
                end = item_1.end + self.downstrand
            elif item_1.feature.strand == -1:
                start = item_1.start - self.downstrand
                end = item_1.end + self.upstrand
            else:
                start = item_1.start - max(self.upstrand,self.downstrand)
                end = item_1.end + max(self.upstrand,self.downstrand)
            if item_1.feature.seqid in index:
                for item_2 in index[item_1.feature.seqid].get(start,end):
                    item_1.relations.append(item_2)
                    item_2.relations.append(item_1)

        for item in features_parent:
            item.modify_with_relations(self.use, self.to_child, self.to_parent)
        
        with open(self.prefix + '-parent.gff','wb') as f:
            annotation.write_gff3_header(f)
            for item in features_parent:
                print >> f, item.feature.as_gff()
        
        with open(self.prefix + '-child.gff','wb') as f:
            annotation.write_gff3_header(f)
            for item in features_child:
                print >> f, item.feature.as_gff()




