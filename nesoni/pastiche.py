
"""

Paste together contigs to cover a reference

"""


import sys, os, collections, subprocess

from nesoni import io, bio, grace

USAGE = """

Use MUMMER to plaster a set of contigs over reference sequences.

Usage:

    nesoni pastiche: [options] <output-dir> <reference.fa> [...] contigs: <contigs.fa> [...]

Options:

    --mask yes/no  - What to do if nothing matches part of the reference.
                     If yes, retain reference in lowercase.
                     If no, output "n"s.
                     Default: no 

    --min-leftover nnn
                   - Discard leftover sequence shorter than this.
                     Default: 20

"""


def run(args, stdin=None, stdout=subprocess.PIPE, stderr=None, bufsize=1<<20):
    return subprocess.Popen(
        args,
        bufsize=bufsize,
        stdin=stdin,        
        stdout=stdout,
        stderr=stderr,
        close_fds=True,
    )

def execute(args, stdin=None, stdout=None):
    p = run(args, stdin=stdin, stdout=stdout, stderr=subprocess.PIPE)
    for line in p.stderr:
        #sys.stdout.write(line)
        pass
    assert p.wait() == 0, 'Failed to execute "%s"' % ' '.join(args)

    
def pastiche(args):
    if len(args) < 4:
        print USAGE
        return 1

    mask_only, args = grace.get_option_value(args, '--mask', grace.as_bool, False)
    min_leftover, args = grace.get_option_value(args, '--min-leftover', int, 20)
        
    output_dir, args = args[0], args[1:]
    
    #, ref_filename, contig_filenames = args[0], args[1], args[2:]
    
    ref_filenames = [ ]
    contig_filenames = [ ]
    grace.execute(args, {
        'contigs' : lambda args: contig_filenames.extend(args)
    }, lambda args: ref_filenames.extend(args))
    
    assert ref_filenames, 'No reference sequences given'
    assert contig_filenames, 'No contig sequences given'
    
    contigs = dict([ 
                 (name.split()[0], seq) 
                 for filename in contig_filenames 
                 for name, seq in io.read_sequences(filename) 
              ])
    dir_contigs = { }
    for name in contigs:
        dir_contigs[name + '+'] = contigs[name]
        dir_contigs[name + '-'] = bio.reverse_complement(contigs[name])
    
    dir_contigs_used = { }
    for name in dir_contigs:
        dir_contigs_used[name] = [ False ] * len(dir_contigs[name])


    workspace = io.Workspace(output_dir)
    temp_prefix = workspace._object_filename('temp-pastiche')
    
    out_f = workspace.open('pastiche.fa', 'wb')
    
    for ref_filename in ref_filenames:
      for ref_name, ref_seq in io.read_sequences(ref_filename):
        ref_name = ref_name.split()[0]
        
        grace.status(ref_name)
        
        f = open(temp_prefix + '.fa','wb')
        io.write_fasta(f, 'ref', ref_seq)
        f.close()
    
        scores = [ -1 ] * (len(ref_seq)*2)
        strings = [ 'N', '' ] * (len(ref_seq))
        contexts = [ None for i in xrange(len(ref_seq)*2) ]
        
        #MAXSCORE = len(ref_seq)+1
        #for i in xrange(len(ref_seq)):
        #    if ref_seq[i].upper() != 'N':
        #        strings[i*2] = ref_seq[i]
        #        scores[i*2] = MAXSCORE
        #for i in xrange(len(ref_seq)-1):
        #    if ref_seq[i].upper() != 'N' and ref_seq[i+1].upper() != 'N':
        #        scores[i*2+1] = MAXSCORE

        if mask_only:        
            for i in xrange(len(ref_seq)):
                strings[i*2] = ref_seq[i].lower()
        
        
        def put(position, dir_contig_name, start, end, score):
            if scores[position] < score:
                scores[position] = score
                strings[position] = dir_contigs[dir_contig_name][start:end]
                contexts[position] = (dir_contig_name, start, end, score)

        for contig_filename in contig_filenames:
            execute(['nucmer',
                     '--prefix', temp_prefix,
                     #'--maxmatch', #Very slow
                     '--nosimplify',
                     '--minmatch', '9',
                     '--mincluster', '50',
                     #'--maxgap', '1000',
                     #'--breaklen', '1000', # Increasing this reduces Ns, but is slow
                     #'--diagfactor', '1.0',
                     temp_prefix+'.fa',
                     contig_filename])
            
            for contig_name, contig_seq in io.read_sequences(contig_filename):
                contig_name = contig_name.split()[0]
                grace.status(ref_name + ' vs ' + contig_name)
                p = run(['show-aligns', temp_prefix+'.delta', 'ref', contig_name],
                        stderr=subprocess.PIPE)
                
                alignments = [ ]
                
                while True:
                    line = p.stdout.readline()
                    if not line: break
                    if not line.startswith('-- BEGIN'):
                        continue
                    
                    parts = line.split()
                    
                    ref_start = int(parts[5])
                    ref_end = int(parts[7])
                    query_start = int(parts[10])
                    query_end = int(parts[12])
                    
                    #assert ref_start < ref_end
                    #ref_start -= 1 #Zero based coordinates
                    
                    al_ref = [ ]
                    al_query = [ ]
                    
                    while True:
                        block = [ ]
                        end = False
                        while True:
                            line = p.stdout.readline()
                            if line.startswith('--   END'): 
                                end = True
                                break
                            if line == '\n':
                                if block: 
                                    break
                                else:
                                    continue
                            block.append(line)
                        
                        if end: break
                        
                        al_ref.append(block[0].split()[1])
                        al_query.append(block[1].split()[1])
                        
                    al_ref = ''.join(al_ref)
                    al_query = ''.join(al_query)            
                    
                    if ref_start > ref_end:
                       al_ref = bio.reverse_complement(al_ref)
                       al_query = bio.reverse_complement(al_query)
                       ref_start, ref_end = ref_end, ref_start
                       query_start, query_end = query_end, query_start
                    
                    if query_start > query_end:
                       dir_contig_name = contig_name + '-'
                       query_start = len(contig_seq)+1-query_start
                       query_end = len(contig_seq)+1-query_end
                    else:
                       dir_contig_name = contig_name + '+'
                       
                    ref_start -= 1 #Zero based coordinates
                    query_start -= 1
                    
                    #print al_ref
                    #print al_query
                    
                    #Pretty dumb scoring scheme
                    al_score = 0
                    for i in xrange(len(al_ref)):
                        if al_ref[i] == al_query[i]:
                            al_score += 1
                        #else:
                        #    al_score -= 1
                    
                    #Pastiche alignment over reference
                    ref_pos = ref_start
                    query_pos = query_start
                    al_pos = 0
                    while al_pos < len(al_ref):
                        assert al_ref[al_pos] != '.'                
                        if al_query[al_pos] == '.':
                            put(ref_pos*2, dir_contig_name, query_pos, query_pos, al_score)
                        else:
                            assert al_query[al_pos].lower() == dir_contigs[dir_contig_name][query_pos].lower()
                            put(ref_pos*2, dir_contig_name, query_pos, query_pos+1, al_score)
                            query_pos += 1
                        al_pos += 1
                        
                        al_pos_end = al_pos
                        query_pos_end = query_pos
                        while al_pos_end < len(al_ref) and al_ref[al_pos_end] == '.':
                            al_pos_end += 1
                            query_pos_end += 1
                        #put(ref_pos*2+1, al_query[al_pos:al_pos_end], al_score)
                        assert al_query[al_pos:al_pos_end].lower() == dir_contigs[dir_contig_name][query_pos:query_pos_end].lower() 
                        put(ref_pos*2+1, dir_contig_name, query_pos,query_pos_end, al_score)
                        al_pos = al_pos_end
                        query_pos = query_pos_end
                        ref_pos += 1
                    
                    
                p.wait()
            
        grace.status(ref_name)
        
        result = ''.join(strings)    
        io.write_fasta(out_f, ref_name, result)
        
        
        for context in contexts:
            if context is None: continue
            name,start,end,score = context
            for i in xrange(start,end):
                dir_contigs_used[name][i] = True
        
        
        #Interpolation
        #result = [ ]
        #i = 0
        #while i < len(ref_seq):
        #    if strings[i*2].upper() != 'N':
        #        result.append(strings[i*2])
        #        result.append(strings[i*2+1])
        #        i += 1
        #        continue
        #    
        #    j = i
        #    while strings[j*2].upper() == 'N':
        #        j += 1
        #    
        #    grace.status('')
        #    print >> sys.stderr, 'interpolating', i+1,'..',j
        #    
        #    window = 20 #!!!!!!!!!!!
        #    left_contexts = collections.defaultdict(lambda:0)
        #    for i1 in xrange(max(0,i-window),i):
        #        for context_name, context_start, context_end, context_score in contexts[i1*2]:
        #            key = (context_name, context_end + i - i1)
        #            left_contexts[key] = max(left_contexts[key],context_score)
        #        
        #    right_contexts = collections.defaultdict(lambda:0)
        #    for j1 in xrange(j,min(j+window,len(ref_seq))):
        #        for context_name, context_start, context_end, context_score in contexts[j1*2]:
        #            key = (context_name, context_start + j - j1)
        #            right_contexts[key] = max(left_contexts[key],context_score)
        #    
        #    #print >> sys.stderr, left_contexts
        #    #print >> sys.stderr, right_contexts
        #    
        #    options = [ ]
        #    
        #    for (left_name, left_pos), left_score in left_contexts.items():
        #        for (right_name, right_pos), right_score in right_contexts.items():
        #            if left_name != right_name: continue
        #            if right_pos < left_pos: continue
        #            
        #            if right_pos-left_pos > (j-i) * 4.0 + 10: continue   #!!!!!!!!!!!!!!!!!!!!!!1
        #            if right_pos-left_pos < (j-i) * 0.25 - 10: continue
        #            
        #            score = float(min(right_pos-left_pos,j-i))/max(right_pos-left_pos,j-i)                  
        #            score *= left_score + right_score
        #            #print >> sys.stderr, left_name, right_pos-left_pos, j-i, score
        #            options.append( (score, left_name, left_pos, right_pos) )
        #    
        #    if options:
        #        best = max(options, key=lambda option: option[0])
        #        print >> sys.stderr, '->', best
        #        result.append( dir_contigs[best[1]][best[2]:best[3]].lower() )
        #    else:
        #        print >> sys.stderr, '-> no good interpolation'
        #        result.append( ref_seq[i:j] )
        #    
        #    i = j
        #
        #result = ''.join(result)    
        #io.write_fasta(sys.stdout, ref_name, result)
        
        
        #print >> sys.stderr, len(result), result.count('N')
        #for pos, size in N_runs:
        #    out_size = len(''.join( strings[pos*2:pos*2+2] ))
        #    print >> sys.stderr, pos, size, '->', out_size        
    
    out_f.close()
    
    grace.status('')
    
    #for name, seq in io.read_sequences(ref_filename):
    #    result = pastiche(seq, contigs_filename)
    #    io.write_fasta(sys.stdout, name, result)
    
    
    leftover_f = workspace.open('leftovers.fa','wb')

    for name in sorted(contigs):
        used = [ (a or b) for a,b in zip(dir_contigs_used[name+'+'],dir_contigs_used[name+'-'][::-1]) ]

        i = 0
        while i < len(used):
            j = i
            while j < len(used) and not used[j]: 
                j += 1
            if j-i > min_leftover:
                if i == 0 and j == len(used):
                    out_name = name
                else:
                    out_name = name + ':%d..%d' % (i+1,j)
                io.write_fasta(leftover_f, out_name, contigs[name][i:j])
            
            i = j+1        

    leftover_f.close()

    for suffix in ['.fa', '.delta']:
        os.unlink(temp_prefix + suffix)



