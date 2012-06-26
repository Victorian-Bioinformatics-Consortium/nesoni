
import sys, os, collections, itertools

from nesoni import io, grace, working_directory, config

USAGE = """\

Usage:

    nesoni nway: [options] workingdir1 workingdir2 [...] \\
                 [splitting: ... from: ...]

Options:

    --gbk file.gbk       annotate with features from GENBANK file 

    --as table/compact/nexus/counts
                         output format, see below
                         default: table

    --evidence yes/no    include evidence in table output
                         default: yes
                         
    --consequences yes/no 
                         include protein-level consequences in table output
                         default: yes

    --reference yes/no   include reference
                         default: yes
    
    --indels yes/no      use insertions and deletions
                         default: yes

    --require-all yes/no require all inputs to have a consensus
                         default: no

    --require-bisect yes/no   
                         only consider positions that split the
                         sample set in two (implies --require-all)
                         default: no
    
    --full yes/no        output uninteresting positions
                         default: no
    
    splitting: workingdirA1 [workingdirA2 ...]
    from: workingdirB1 [workingdirB2 ...]
                         consider only sites that partition As from Bs
                         use "reference" to specify the reference
    
Output formats:

    table      - output as table with consensus calls and evidence
    
    compact    - compact output of consensus only
    
    nexus      - NEXUS format, suitable for SplitsTree, etc
                 Note: Will output missing sites unless
                 you use --require-all. Your phylogenetic
                 software may or may not do something sane
                 with this.
    
    counts     - output partition counts
                 (implies --require-all) 

Lack of consensus or ambiguity codes in any data set will cause differences 
at that position to be ignored.

Three or more different bases/insertions at a given position will cause it to be
ignored.

"""


AMBIGUOUS = set('KMRYSWBVHDN')

Call = collections.namedtuple('Call', 'consensus evidence consequences')

Calls = collections.namedtuple('Calls', 'ref_name ref_pos is_insertion calls features')

def do_fasta_output(names, interesting):
    for i, name in enumerate(names):
        sequence = [ ]
        #for refname in substitution_calls:
        #    sequence.extend(substitution_calls[refname][i])
        #for i in xrange(len(sequence)):
        #    if sequence[i] not in 'ACGT':
        #        sequence[i] = 'N'

        for refname, position, change_type, values, has_ambiguous, evidence in interesting:
            if not has_ambiguous and change_type == 'substitution':
                sequence.append(values[i]) 
            
        io.write_fasta(sys.stdout, name + ' variable sites with consensus in all', ''.join(sequence))


def reference_reader(sequence):
    for i in xrange(len(sequence)):
        yield Call('-', '', '')
        yield Call(sequence[i], '', '') 

def evidence_reader(working_dir, name):
    filename = os.path.join(working_dir, grace.filesystem_friendly_name(name) + '-evidence.txt')
    f = open(filename,'rb')
    
    header = f.readline()
    if header.count('\t') != 7:
        raise grace.Error('Old style evidence file. Please re-run nesoni consensus.')
    
    for line in f:
        fields = line.rstrip('\n').split('\t')
        yield Call(fields[4], fields[1], fields[6])
        yield Call(fields[5], fields[2], fields[7])

    f.close()

def reader(working_dirs, references, use_reference, annotations={}):
    for name, sequence in references:
        features = annotations.get(sequence, [])
    
        if use_reference:
            readers = [ reference_reader(sequence) ]
        else:
            readers = [ ]
        
        readers.extend( evidence_reader(working_dir, name) for working_dir in working_dirs )
        
        active_features = [ ]
        feature_pos = 0        
        
        for i in xrange(len(sequence)):
            if i % 10000 == 0:
                grace.status('%s %s' % (name, grace.pretty_number(i)))
            
            active_features = [ item for item in active_features if item.location.nofuzzy_end > i ]
            while feature_pos < len(features) and \
                  features[feature_pos].location.nofuzzy_start <= i:
                active_features.append(features[feature_pos])
                feature_pos += 1
        
            for is_insertion in (True, False):
                 yield Calls(name, i, is_insertion, [ item.next() for item in readers ], active_features)
        
        for reader in readers:
            for item in reader:
                raise grace.Error('Unexpected extra data in evidence file')

    grace.status('')

def is_interesting(calls):
    unambiguous = set( item.consensus for item in calls.calls if item.consensus not in AMBIGUOUS )
    return len(unambiguous) >= 2

def is_binary_partition(calls):
    return len(set( item.consensus for item in calls.calls )) == 2

def is_sane_partition(calls): #Can be converted into 0123456789
    return 2 <= len(set( item.consensus for item in calls.calls )) <= 10

def is_split(split_a, split_b):
    all = split_a + split_b
    def func(calls):
        for i in all:
            if calls.calls[i].consensus in AMBIGUOUS: return False
        set1 = set( calls.calls[i].consensus for i in split_a )
        set2 = set( calls.calls[i].consensus for i in split_b )
        return len(set1 & set2) == 0
    return func


def fully_unambiguous(calls):
    for call in calls.calls:
        if call.consensus in AMBIGUOUS: return False
    return True
 
def has_no_indels(calls):
    if calls.is_insertion: return False
    for call in calls.calls:
        if call.consensus == '-': return False
    return True

def not_boring_insertion(calls):
    if not calls.is_insertion: return True
    for call in calls.calls:
        if call.consensus not in ('-','N'): return True
    return False

def partition_string(calls):
    renaming = { }
    for char in 'ACGT-': 
        renaming[char] = char
    n = 0        
    for call in calls.calls:
        if call.consensus not in renaming and call.consensus not in AMBIGUOUS:
            renaming[call.consensus] = chr(48 + n)
            n += 1
    
    return ''.join( renaming.get(item.consensus,'N') for item in calls.calls )

def change_type(calls):
    unambiguous = set( item.consensus for item in calls.calls if item.consensus not in AMBIGUOUS )
    if len(unambiguous) <= 1: return ''
    if calls.is_insertion:
        return 'insertion-before'
    if '-' in unambiguous:
        return 'deletion'
    return 'substitution'

def transpose_strings(strings, pad=' ', min_n=0, pad_width=40):
    joiner = ' '*min(2,pad_width // len(strings))
    n = max(min_n, max( len(item) for item in strings ))
    return [ joiner.join( item[i] if i < len(item) else ' ' for item in strings )
             for i in xrange(n) ]


def norm_name(working_dir):
    return os.path.basename(os.path.normpath(working_dir))

def describe_features(features):
    result = [ ]
    for feature in features:
        info = [feature.type] + feature.qualifiers.get('locus_tag',[]) + \
               feature.qualifiers.get('gene',[]) + \
               feature.qualifiers.get('product',[])
        result.append(' '.join(info))

    return ', '.join(result)


@config.help('Compare results of several runs of nesoni consensus, amongst themselves and optionally with the reference.',
"""\
"splitting: A1 A2... from: B1 B2..." can be used to consider only sites that partition As from Bs. \
Use "reference" to specify the reference.

Output formats:

    table      - output as table with consensus calls and evidence
    
    compact    - compact output of consensus only
    
    nexus      - NEXUS format, suitable for SplitsTree, etc
                 Note: Will output missing sites unless
                 you use --require-all. Your phylogenetic
                 software may or may not do something sane
                 with this.
    
    counts     - output partition counts
                 (implies --require-all) 
""")
@config.String_flag('as_','Output format, see below.\n'
                          'Options: table compact nexus counts\n')
@config.String_flag('gbk','annotate with features from GENBANK file')
@config.Bool_flag('evidence','include evidence in table output')
@config.Bool_flag('consequences','include protein-level consequences in table output')
@config.Bool_flag('reference','include reference in comparison')
@config.Bool_flag('indels','use insertions and deletions')
@config.Bool_flag('require_all','require all inputs to have a consensus')
@config.Bool_flag('require_bisect','only consider positions that split the sample set in two (implies --require-all)')
@config.Bool_flag('full','output uninteresting positions')
@config.Main_section('working_dirs','Working directories', empty_is_ok=False)
@config.Section('splitting')
@config.Section('from_')
class Nway(config.Action_with_optional_output):
    as_ = 'table'
    gbk = None
    indels = True
    reference = True
    evidence = True
    consequences = True
    require_all = False
    require_bisect = False
    full = False

    working_dirs = [ ]
    splitting = [ ]
    from_ = [ ]
    
    def run(self):
        f = self.begin_output()
        nway_main(gbk_filename=self.gbk, use_indels=self.indels, use_reference=self.reference, 
                  give_evidence=self.evidence, give_consequences=self.consequences,
                  require_all=self.require_all, require_bisect=self.require_bisect, 
                  full_output=self.full, format=self.as_, 
                  working_dirs=self.working_dirs, 
                  split_a=self.splitting, split_b=self.from_, f=f)
        self.end_output(f)

#def main(args):
#    gbk_filename, args = grace.get_option_value(args,'--gbk',str,None)
#    use_indels, args = grace.get_option_value(args,'--indels',grace.as_bool,True)
#    use_reference, args = grace.get_option_value(args,'--reference',grace.as_bool,True)
#    give_evidence, args = grace.get_option_value(args,'--evidence',grace.as_bool,True)
#    give_consequences, args = grace.get_option_value(args,'--consequences',grace.as_bool,True)
#    require_all, args = grace.get_option_value(args,'--require-all',grace.as_bool,False)
#    require_bisect, args = grace.get_option_value(args,'--require-bisect',grace.as_bool,False)
#    full_output, args = grace.get_option_value(args,'--full',grace.as_bool,False)
#    format, args = grace.get_option_value(args,'--as',str,'table')
#    
#    ## Secret option!
#    #limit, args = grace.get_option_value(args,'--limit',int,None)
#    
#    grace.expect_no_further_options(args)
#
#    if len(args) < 1:
#        sys.stderr.write(USAGE)
#        return 1
#
#    working_dirs = [ ]
#    split_a = [ ]
#    split_b = [ ]
#    def default(args):
#        working_dirs.extend(args)
#    def splitting(args):
#        split_a.extend(args)
#    def splitting_from(args):
#        split_b.extend(args)
#        
#    grace.execute(args, {
#        'splitting' : splitting,
#        'from' : splitting_from 
#    }, default
#    )

def nway_main(gbk_filename, use_indels, use_reference, give_evidence, give_consequences,
              require_all, require_bisect, full_output, format, working_dirs, split_a, split_b, f=sys.stdout):
    assert working_dirs, 'Need at least one working directory.'
    workspaces = [ working_directory.Working(dirname, must_exist=True) for dirname in working_dirs ]
    reference = workspaces[0].get_reference()
    #if not annotation_filename:
    #    annotation_filename = reference.annotations_filename() #May still be None
    
    if use_reference:
        names = ['reference']
        evidence_start = 1
    else:
        names = [ ]
        evidence_start = 0
        
    names.extend( norm_name(item) for item in  working_dirs )
    
    references = io.read_sequences(reference.reference_fasta_filename())
    
    annotations = { }
    if gbk_filename:
        from Bio import SeqIO
        for record in SeqIO.parse(io.open_possibly_compressed_file(gbk_filename),'genbank'):
            sequence = record.seq.tostring()
            features = [ item for item in record.features if item.type != 'source' ]
            features.sort(key=lambda item: item.location.nofuzzy_start)
            annotations[sequence] = features
    
    iterator = reader(working_dirs, references, use_reference, annotations)
    
    if not use_indels:
        iterator = itertools.ifilter(has_no_indels, iterator)

    if require_all or require_bisect or format == 'counts':
        iterator = itertools.ifilter(fully_unambiguous, iterator)
    
    if require_bisect:
        iterator = itertools.ifilter(is_binary_partition, iterator)

    if not require_bisect:
        if full_output:
            iterator = itertools.ifilter(not_boring_insertion, iterator)
        else:
            iterator = itertools.ifilter(is_interesting, iterator)

    if split_a or split_b:
        assert len(names) == len(set(names)), 'Two samples with the same name'
        try:
            split_a = [ names.index(norm_name(item)) for item in split_a ]
            split_b = [ names.index(norm_name(item)) for item in split_b ]
        except ValueError:
            raise grace.Error('Sample to be split is not amongst samples given')
        iterator = itertools.ifilter(is_split(split_a, split_b), iterator)

    #if limit:
    #    iterator = itertools.islice(iterator, limit)
    
    if format == 'table':
        line = 'Reference\tPosition\tChange type'
        line +=  '\t' + '\t'.join(names)
        if give_evidence:
            line += '\t' + '\t'.join(names[evidence_start:])
        if give_consequences:
            line += '\t' + '\t'.join(names[evidence_start:])
        if annotations:
            line += '\tAnnotations'
        print >> f, line
        for calls in iterator:
            line = '%s\t%d\t%s\t%s' % (
                calls.ref_name, 
                calls.ref_pos+1, 
                change_type(calls), 
                '\t'.join(item.consensus for item in calls.calls))
            if give_evidence:
                line += '\t' + '\t'.join(item.evidence for item in calls.calls[evidence_start:])
            if give_consequences:
                line += '\t' + '\t'.join(item.consequences for item in calls.calls[evidence_start:])
            if annotations:
                line += '\t' + describe_features(calls.features)
            print >> f, line

    elif format == 'compact':
        for line in transpose_strings(names):
            print >> f, line
        print >> f
        
        for calls in iterator:
            if calls.is_insertion:
                footer = '%12d.5 %s' % (calls.ref_pos, calls.ref_name)
            else: 
                footer = '%12d   %s' % (calls.ref_pos+1, calls.ref_name)
            
            t = transpose_strings([ item.consensus for item in calls.calls ], '-', 1)
            top = t[0] + ' ' + footer
            if give_consequences:
                consequences = [ ]
                for call in calls.calls:
                    if call.consequences:
                        for item in call.consequences.split(', '):
                            item = ' '.join(item.split()[:3])
                            if item not in consequences: consequences.append(item)
                        
                if consequences:
                    top += '  ' + ' / '.join(sorted(consequences))
            top += '  ' + describe_features(calls.features)
            print >> f, top
            for line in t[1:]:
                print >> f, line            
    
    elif format == 'nexus':
        buckets = [ [ ] for name in names ]
        for calls in iterator:
            for i, char in enumerate(partition_string(calls)):
                buckets[i].append(char)
        
        print >> f, '#NEXUS'
        print >> f, 'begin taxa;'
        print >> f, 'dimensions ntax=%d;' % len(names)
        print >> f, 'taxlabels'
        for name in names:
            print >> f, name
        print >> f, ';'
        print >> f, 'end;'

        print >> f, 'begin characters;'
        print >> f, 'dimensions nchar=%d;' % len(buckets[0])
        print >> f, 'format datatype=STANDARD symbols="ACGT-0123456789" missing=N;'
        print >> f, 'matrix'
        for name, bucket in itertools.izip(names, buckets):
            print >> f, name, ''.join(bucket)
        print >> f, ';'
        print >> f, 'end;'
    
    elif format == 'counts':
        for line in transpose_strings(names):
            print >> f, line
        print >> f

        counts = { }
        for calls in iterator:
            count_str = partition_string(calls)
            if count_str not in counts:
                counts[count_str] = 1
            else:
                counts[count_str] += 1
        
        for count_str in sorted(counts, key=lambda x: (counts[x], x), reverse=True):
            print >> f, '%s   %d' % (transpose_strings(count_str)[0], counts[count_str])
    
    else:
        raise grace.Error('Unknown output format: ' + format)
        
    #for calls in iterator:
    #    print calls.ref_name, calls.ref_pos+1, ''.join( item.consensus for item in calls.calls )




#def find_interesting(what, calls, evidences):
#    result = [ ]
#    for refname in calls:
#        for i, row in enumerate(zip(*calls[refname])):
#            unambiguous = [ item for item in row if item not in AMBIGUOUS ]
#            if len(set(unambiguous)) <= 1: continue
#            
#            evidence = [ evidences[refname][j][i] for j in xrange(len(row)) ]
#            
#            result.append( (refname, i+1, what, row, len(row) != len(unambiguous), evidence) )
#
#    return result
#
#
#def old_main(args):
#    use_indels, args = grace.get_option_value(args,'--indels',int,1)
#    use_reference, args = grace.get_option_value(args,'--reference',int,1)
#    make_list, args = grace.get_option_value(args,'--list',int,0)
#    fasta_output, args = grace.get_option_value(args,'--fasta',int,0)
#    grace.expect_no_further_options(args)
#    
#    if len(args) < 1:
#        sys.stderr.write(USAGE)
#        return 1
#        
#    if fasta_output and use_indels:
#        print >> sys.stderr, 'Indels will not be included in FASTA output'
#        use_indels = 0
#    
#    working_dirs = args
#    
#    #reference_data = { } # (ref_name, position, change_type) -> string
#    #strain_data = { } # working_dir -> (ref_name, position, change_type) -> string
#    
#    names = ['reference'] + working_dirs
#    
#    substitution_calls = { } # ref_name -> [ [ call ] ]
#    insertion_calls = { } # ref_name -> [ [ call ] ]
#    substitution_evidence = { }
#    insertion_evidence = { }
#    
#    for name, sequence in io.read_sequences(os.path.join(working_dirs[0], 'reference.fa')):
#        substitution_calls[name] = [ list(sequence.upper()) ]
#        insertion_calls[name] = [ [ '-' ] * len(sequence) ]
#        substitution_evidence[name] = [ [ '' ] * len(sequence) ]    
#        insertion_evidence[name] = [ [ '' ] * len(sequence) ]    
#    
#    for working_dir in working_dirs:
#        for name in substitution_calls:
#            filename = os.path.join(working_dir, grace.filesystem_friendly_name(name) + '-evidence.txt')
#            f = open(filename,'rb')
#            
#            this_substitution_calls = [ ]
#            this_insertion_calls = [ ]
#            this_substitution_evidence = [ ]
#            this_insertion_evidence = [ ]
#            
#            header = f.readline()
#            if header.count('\t') != 5:
#                print >> sys.stderr, 'Old style evidence file. Please re-run nesoni consensus.'
#                return 1
#            
#            for line in f:
#                fields = line.rstrip('\n').split('\t')
#                this_substitution_calls.append(fields[5])
#                this_insertion_calls.append(fields[4])
#                this_substitution_evidence.append(fields[2])
#                this_insertion_evidence.append(fields[1])
#            
#            substitution_calls[name].append(this_substitution_calls)
#            insertion_calls[name].append(this_insertion_calls)
#            substitution_evidence[name].append(this_substitution_evidence)
#            insertion_evidence[name].append(this_insertion_evidence)
#    
#    if not use_reference:
#        names.pop(0)
#        for name in substitution_calls:
#            substitution_calls[name].pop(0)
#            insertion_calls[name].pop(0)
#            substitution_evidence[name].pop(0)
#            insertion_evidence[name].pop(0)
#
#    interesting = find_interesting('substitution', substitution_calls, substitution_evidence)
#    if use_indels:
#        interesting.extend( find_interesting('insertion-before', insertion_calls, insertion_evidence) )
#
#    if not use_indels:
#        interesting = [ item for item in interesting if '-' not in item[3] ]
#    
#    interesting.sort()
#
#
#    if fasta_output:
#        do_fasta_output(names, interesting)
#        return 0 
#
#    
#    #strain_reference_having_consensus = { } # working_dir -> ref_name -> string
#    #
#    #for working_dir in working_dirs:
#    #    assert working_dir not in strain_data, 'Working directory given twice'
#    #    strain_data[working_dir] = { }
#    #    
#    #    report_file = open(os.path.join(working_dir, 'report.txt'), 'rU')
#    #    report_file.readline()
#    #    for line in report_file:
#    #        ref_name, position, change_type, old, new, evidence = \
#    #            line.rstrip('\n').split('\t')
#    #        
#    #        if change_type == 'deletion':
#    #            change_type = 'substitution'
#    #        
#    #        if not use_indels and \
#    #           (change_type == 'insertion-before' or new == '-'):
#    #            continue
#    #        
#    #        key = (ref_name, int(position), change_type)
#    #        if key in reference_data:
#    #            assert reference_data[key] == old
#    #        else:
#    #            reference_data[key] = old
#    #        
#    #        strain_data[working_dir][key] = new
#    #    report_file.close()
#    #    
#    #    strain_reference_having_consensus[working_dir] = { }
#    #    ref_have_con_filename = os.path.join(working_dir, 'reference_having_consensus.fa')
#    #    for name, sequence in io.read_fasta(ref_have_con_filename):
#    #        strain_reference_having_consensus[working_dir][name] = sequence
#    #
#    #keys = sorted(reference_data)
#    #
#    ##Fill in any blanks
#    #for working_dir in working_dirs:
#    #    for key in keys:
#    #        if key in strain_data[working_dir]: continue
#    #    
#    #        # - Positions in report files start from 1 not 0
#    #        # - Insertions must be bracketed
#    #        lacks_consensus = (
#    #            strain_reference_having_consensus[working_dir][key[0]][key[1]-1] == 'N' or
#    #            (key[2] == 'insertion-before' and key[1] > 1 and
#    #             strain_reference_having_consensus[working_dir][key[0]][key[1]-2] == 'N')
#    #        )
#    #        
#    #        #If there's no consensus, record it as ambiguous
#    #        if lacks_consensus:
#    #            strain_data[working_dir][key] = 'N'                
#    #        else:
#    #            strain_data[working_dir][key] = reference_data[key]
#
# 
#    #all_data_names = ([ 'reference' ] if use_reference else []) + working_dirs
#    #all_data = ([ reference_data ] if use_reference else []) + \
#    #           [ strain_data[working_dir] for working_dir in working_dirs ] 
#    
#
#    #all_data_names = ([ 'reference' ] if use_reference else []) + working_dirs
#    
#    
#
#    
#    
#    ones = ( 1 << len(names) )-1
#    
#    total_differences = 0
#    
#    if make_list:
#        print '\t'.join(['Partition','Sequence','Position in reference','Change type'] + names + names) 
#    
#    for i in xrange(1,(1<<len(names))-1,2):
#        set1 = [ ]
#        set2 = [ ]
#        for j in xrange(len(names)):
#            if i & (1<<j):
#                set1.append(j)
#            else:
#                set2.append(j)
#
#        if make_list:
#            print
#            print ', '.join( names[i] for i in set1 ) + '   vs   ' + \
#                  ', '.join( names[i] for i in set2 )
#            print
#                
#        n = 0
#        for refname, position, change_type, values, has_ambiguous, evidence in interesting: 
#            #Skip if *any* ambiguity
#            if has_ambiguous:
#                continue
#            
#            if any( values[i] != values[set1[0]] for i in set1[1:] ) or \
#               any( values[i] != values[set2[0]] for i in set2[1:] ):
#                continue
#            
#            if make_list:
#                if change_type == 'substitution' and '-' in values: change_type = 'deletion'
#                print '\t%s\t%d\t%s\t' % (refname,position,change_type) + '\t'.join(values) + '\t' + '\t'.join(evidence) 
#            
#            n += 1
#
#        total_differences += n
#
#        if not make_list:
#            print ', '.join( names[i] for i in set1 ) + '   vs   ' + \
#                  ', '.join( names[i] for i in set2 ) + \
#                  ': %d differences' %n            
#
#    if not make_list:
#        print
#        print 'Total: %d' % total_differences
#
#
#    if make_list:
#        print
#        print 'Ignored'
#        print
#    
#    n_multiway = 0
#    n_ambiguous = 0    
#    for refname, position, change_type, values, has_ambiguous, evidence in interesting: 
#        confusing = False
#        if has_ambiguous:
#            n_ambiguous += 1
#            confusing = True
#        elif len(set(values)) > 2:
#            n_multiway += 1
#            confusing = True
#        
#        if make_list and confusing:
#            print '\t%s\t%d\t%s\t' % (refname,position,change_type) + '\t'.join(values) + '\t' + '\t'.join(evidence) 
#
#    if not make_list:
#        print
#        print 'Ambiguities ignored: %d' % n_ambiguous
#        print 'Multi-way changes ignored: %d' % n_multiway
#    
#    assert total_differences + n_ambiguous + n_multiway == len(interesting)
#    
#    return 0
#
