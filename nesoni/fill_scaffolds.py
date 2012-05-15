
from nesoni import io, bio, grace

import sys, os, collections

USAGE = """

Fill gaps in a 454 scaffold with paths through the contig
graph having maximum minimum coverage along their length.

This might be useful if you want a *very* *rough* idea
of what might be lurking in the gaps. 

Usage:

    nesoni fill-scaffolds: [options] output-dir 454-assembly-dir \\
        [scaffold: [--circular yes/no] +/-n +/-n ...]...

Options:

    --max-filler nnn   - Maximum length contig to use when filling
                         gaps in a scaffold.
                         Default: 4000                      

You must have specified minimum contig length 0 in Newbler, so that
it creates the file 454ContigGraph.txt.

The "scaffold" section lets you specify a custom scaffold, eg:

    scaffold: --circular yes +5 -10 +3           

specifies a circular molecule consisting of contig00005 forward,
contig00010 reversed, then contig00003 forward.

"""

def fill_scaffolds(args):
    max_filler_length, args = grace.get_option_value(args, '--max-filler', int, 4000)
    
    if len(args) < 2:
        print USAGE
        return 1
    
    (output_dir, graph_dir), args = args[:2], args[2:]

    scaffolds = [ ]
    
    def scaffold(args):
        circular, args = grace.get_option_value(args, '--circular', grace.as_bool, False)
        
        scaffold = [ ]
        for item in args:
            scaffold.append( ('contig', int(item)) )
            scaffold.append( ('gap', None) )
        
        if not circular: scaffold = scaffold[:-1]
        
        name = 'custom_scaffold_%d' % (len(scaffolds)+1)
        scaffolds.append( (name, scaffold) )
            
    grace.execute(args, [scaffold])
    
    custom_scaffolds = (len(scaffolds) != 0)    
    
    sequences = dict( 
        (a.split()[0], b.upper()) 
          for a,b in 
            io.read_sequences(os.path.join(
              graph_dir, '454AllContigs.fna')))
    
    sequence_names = sorted(sequences)
    sequence_ids = dict(zip(sequence_names, xrange(1,len(sequence_names)+1)))
    
    contexts = { }
    context_names = { }
    context_depths = { }
    for i in xrange(1,len(sequence_names)+1):
        seq = sequences[sequence_names[i-1]]
        contexts[ i ] = seq
        context_names[ i ] = sequence_names[i-1]+'-fwd'
        contexts[ -i ] = bio.reverse_complement(seq)
        context_names[ -i ] = sequence_names[i-1]+'-rev'
    
    links = collections.defaultdict(list)
    
    for line in open(
      os.path.join(graph_dir, '454ContigGraph.txt'),
      'rU'):
        parts = line.rstrip('\n').split('\t')
        
        if parts[0].isdigit():
            seq = sequence_ids[parts[1]]
            context_depths[ seq] = float(parts[3])
            context_depths[-seq] = float(parts[3])
        
        if parts[0] == 'C':    
            name1 = 'contig%05d' % int(parts[1])
            dir1 = {"3'" : 1, "5'" : -1 }[parts[2]]
            name2 = 'contig%05d' % int(parts[3])
            dir2 = {"5'" : 1, "3'" : -1 }[parts[4]]
            depth = int(parts[5])
            #print name1, dir1, name2, dir2, depth
            
            links[ sequence_ids[name1] * dir1 ].append( (depth, sequence_ids[name2] * dir2) )
            links[ sequence_ids[name2] * -dir2 ].append( (depth, sequence_ids[name1] * -dir1) )
    
        if parts[0] == 'S' and not custom_scaffolds:  
            name = 'scaffold%05d' % int(parts[2])  
            components = parts[3].split(';')
            scaffold = [ ]
            for component in components:
                a,b = component.split(':')
                if a == 'gap':
                    scaffold.append( ('gap',int(b)) )
                else:
                    strand = { '+': +1, '-': -1 }[ b ]
                    scaffold.append( ('contig', sequence_ids['contig%05d'%int(a)] * strand) )
            scaffolds.append( (name, scaffold) )
    
    
    
    #paths = { }
    #
    #todo = [ ]
    #for i in contexts:
    #    for depth_left, neg_left in links[-i]:
    #        left = -neg_left
    #        for depth_right, right in links[i]:
    #            todo.append( ( max(-depth_left,-depth_right,-context_depths[i]), left, right, (i,)) )
    #
    #heapq.heapify(todo)
    #while todo:
    #    score, source, dest, path = heapq.heappop(todo)
    #    if (source,dest) in paths: continue
    #    
    #    paths[(source,dest)] = path
    #    
    #    if len(contexts[dest]) > max_filler_length: continue
    #    
    #    for depth, next in links[dest]:
    #        heapq.heappush(todo,
    #            ( max(score,-depth,-context_depths[dest]), source, next, path+(dest,))
    #        )
    
    
    path_source_dest = collections.defaultdict(dict) # source -> dest -> next
    path_dest_source = collections.defaultdict(dict) # dest -> source -> next
    
    
    # Use links, in order to depth of coverage, to construct paths between contigs
    # Thus: paths have maximum minimum depth
    #       subsections of paths also have this property
    
    todo = [ ]
    for i in contexts:    
        for depth_link, right in links[i]:
            todo.append( ( depth_link, i, right) )
    todo.sort(reverse=True)
    for score, left, right in todo:
        if right in path_source_dest[left]: continue
        
        sources = [(left,right)]
        if len(contexts[left]) <= max_filler_length:
            sources += path_dest_source[left].items()
        destinations = [right]
        if len(contexts[right]) <= max_filler_length:
            destinations += path_source_dest[right].keys()
        
        for source, next in sources:
            for dest in destinations:
                if dest in path_source_dest[source]: continue
                path_source_dest[source][dest] = next
                path_dest_source[dest][source] = next
    
    
    workspace = io.Workspace(output_dir)
    scaffold_f = workspace.open('scaffolds.fa','wb')
    
    #comments = [ ]
    features = [ ]
    
    used = set()
    previous_total = 0
    
    for i, (name, scaffold) in enumerate(scaffolds):
        result = '' # Inefficient. Meh.
        n_filled = 0
        n_failed = 0
        for j, item in enumerate(scaffold):
            if item[0] == 'contig':
                result += contexts[item[1]]
                used.add(abs(item[1]))
            else:
                left = scaffold[j-1]
                right = scaffold[ (j+1) % len(scaffold) ] #If gap at end, assume circular
                assert left[0] == 'contig'
                assert right[0] == 'contig'
                
                gap_start = len(result)
    
                can_fill = right[1] in path_source_dest[left[1]]
                if can_fill:
                    n = 0
                    k = path_source_dest[left[1]][right[1]]
                    while k != right[1]:
                        n += len(contexts[k])
                        result += contexts[k].lower()
                        used.add(abs(k))
                        
                        k = path_source_dest[k][right[1]]
                    
                    n_filled += 1
                        
                    if item[1] is not None and max(n,item[1]) > min(n,item[1])*4:
                        print >> sys.stderr, 'Warning: gap size changed from %d to %d in scaffold %d' % (item[1],n,i+1)
                else:
                    n_failed += 1
                    
                    #print >> sys.stderr, 'Warning: No path to fill a gap in scaffold %d' % (i+1)
                    result += 'n' * (9 if item[1] is None else item[1])
    
                gap_end = len(result)
                
                #features.append( '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s' % (
                #    'all-scaffolds',
                #    'fill-scaffolds',
                #    'gap',
                #    previous_total + gap_start+1,
                #    previous_total + max(gap_end, gap_start+1), #Allow for zeroed out gaps. Hmm.
                #    '.', #score
                #    '+', #strand
                #    '.', #frame
                #    '' #properties
                #))
                features.append( '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s' % (
                    name,
                    'fill-scaffolds',
                    'gap',
                    gap_start+1,
                    max(gap_end, gap_start+1), #Allow for zeroed out gaps. Hmm.
                    '.', #score
                    '+', #strand
                    '.', #frame
                    '' #properties
                ))
                    
    
        io.write_fasta(scaffold_f, name, result)
        previous_total += len(result)
        #comments.append('##sequence-region    %s %d %d' % (name, 1, len(result)))
        print >> sys.stderr, 'Scaffold%05d: %d gaps filled, %d could not be filled' % (i+1, n_filled, n_failed)
    
    scaffold_f.close()
    
    gff_f = workspace.open('scaffolds.gff', 'wb')
    #print >>gff_f, '##gff-version    3'
    #for comment in comments:
    #    print >>gff_f, comment
    for feature in features:
        print >>gff_f, feature
    gff_f.close()
    
    
    leftovers_f = workspace.open('leftovers.fa', 'wb')
    for name in sequence_names:
        if sequence_ids[name] not in used:
            io.write_fasta(leftovers_f, name, sequences[name])
    leftovers_f.close()
    
    ends = { }
    for i, (name, scaffold) in enumerate(scaffolds):
        if scaffold[-1][0] == 'gap': continue
        ends[ '%s start' % name ] = scaffold[-1][1]
        ends[ '%s end  ' % name ] = -scaffold[0][1] 
    
    for end1 in sorted(ends):
        options = [ end2 for end2 in ends if -ends[end2] in path_source_dest[ends[end1]] ]
        if len(options) == 1:
            print >> sys.stderr, 'Note: from', end1, 'only', options[0], 'is reachable'
    
