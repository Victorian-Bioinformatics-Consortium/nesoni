
import sys, os, itertools, subprocess

import shrimp, consensus, bio, io, grace

def safe_filename(*filename_parts):
    """ Should not look like a flag. """
    return os.path.abspath(os.path.join(*filename_parts))

def run(args, stdin=None, stdout=subprocess.PIPE):
    return subprocess.Popen(
        args,
        stdin=stdin,
        stdout=stdout,
        close_fds=True,
    )

def execute(args, stdin=None, stdout=None):
    p = run(args, stdin=stdin, stdout=stdout)
    assert p.wait() == 0, 'Failed to execute "%s"' % ' '.join(args)

USAGE = """

Usage:

    nesoni tosam working_dir [options]

Options:

    --used-only     Only output hits used by nesoni consensus.
                    You must have specified "--save-hits yes" in nesoni consensus.
                    Default: no

Convert shrimp hits to BAM format.

"""

def main(args):
    used_only, args = grace.get_option_value(args,'--used-only', grace.as_bool, False)
        
    grace.expect_no_further_options(args)

    if len(args) != 1:
        sys.stderr.write( USAGE )
        return 1
    
    working_dir = args[0]
    
    print
    print 'Note: This is still under development'
    print '      Pairing information is not included'
    print '      Only the part of the read that was aligned is included'
    print
    
    if used_only:
        hit_filename = 'used_shrimp_hits.txt.gz'
        output_prefix = 'used_hits'
    else:
        hit_filename = 'shrimp_hits.txt.gz'
        output_prefix = 'hits'

    reference_filename = os.path.join(working_dir,'reference.fa')
    references = dict( io.read_fasta(reference_filename) )
    for name in references:
        references[name] = references[name].upper()
    
    bam_filename = safe_filename(working_dir, output_prefix+'_unsorted.bam')
    bam_sorted_prefix = safe_filename(working_dir, output_prefix)

    f = open(bam_filename, 'wb')
    sam_eater = run(['samtools', 'view', '-S', '-b', '-'],
                    stdin=subprocess.PIPE,
                    stdout=f.fileno())
    f.close()
    
    
    for name in references:
        print >> sam_eater.stdin, '@SQ\tSN:%s\tLN:%d' % (name, len(references[name]))
        
    for i, (read, hits) in enumerate(shrimp.iter_read_hits(working_dir, hit_filename, qualities=True)):
        if (i % 10000) == 0: 
            grace.status('Processing read %s' % grace.pretty_number(i))
    
        for line in hits:
            parts = line.rstrip('\n').split('\t')
            read_name = parts[0]
            ref_name = parts[1]
            ref_start = int(parts[3])-1
            ref_end = int(parts[4])
            read_start = int(parts[5])-1
            read_end = int(parts[6])
            read_length = int(parts[7])
            score = int(parts[8])
            forward = (parts[2] == '+')
            edit_string = parts[9]

            corresp_seq = references[ref_name][ref_start:ref_end]
            
            if not forward:
                corresp_seq = bio.reverse_complement(corresp_seq)
            
            hit_ref_ali, hit_read_ali = consensus.edit_string_to_alignment(edit_string, corresp_seq)	    
            
            if not forward:
                hit_ref_ali = bio.reverse_complement(hit_ref_ali)
                hit_read_ali = bio.reverse_complement(hit_read_ali)
            
            #Normalization -- move "-"s as far right as possible
            hit_read_ali = consensus.roll_alignment(hit_read_ali, hit_ref_ali)
            hit_ref_ali = consensus.roll_alignment(hit_ref_ali, hit_read_ali)
            
            if len(read) > 2:
                qual = read[2][read_start:read_end]
            else:
                qual = '*'
            
            if forward:
                forwardized_seq = read[1][read_start:read_end]
                forwardized_qual = qual
            else:
                forwardized_seq = bio.reverse_complement(read[1][read_start:read_end])
                forwardized_qual = qual[::-1]
            
            #if forwardized_seq != hit_read_ali.replace('-',''):
            #    print line
            #    print read
            #    print corresp_seq
            assert forwardized_seq == hit_read_ali.replace('-',''), forwardized_seq + ' ' + hit_read_ali
            
            cigar_items = [ ]
            for i in xrange(len(hit_ref_ali)):
                if hit_ref_ali[i] == '-':
                    cigar_items.append('I')
                elif hit_read_ali[i] == '-':
                    cigar_items.append('D')
                else:
                    cigar_items.append('M')
            
            cigar = ''.join(
                str(len(list(subiterator))) + key
                for key, subiterator in itertools.groupby(cigar_items)
            )
            
            flags = 0
            if not forward: flags += 0x0010
            
            sam_line = '\t'.join([
                read_name,
                str(flags),
                ref_name,
                str(ref_start+1),
                '255',
                cigar,
                '=',
                '0',
                '0',
                forwardized_seq,
                forwardized_qual,
            ])
            
            print >> sam_eater.stdin, sam_line

    sam_eater.stdin.close()
    assert sam_eater.wait() == 0, 'samtools failed'
    
    grace.status('')
    
    execute([
        'samtools', 'sort', bam_filename, bam_sorted_prefix
    ])
    
    os.unlink(bam_filename)
    
    execute([
        'samtools', 'index', bam_sorted_prefix + '.bam'
    ])
    
    print
    print working_dir + '/' + output_prefix + '.bam and index created'
    print


