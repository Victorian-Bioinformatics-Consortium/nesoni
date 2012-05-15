
from nesoni import grace, io

import sys, os, hashlib, pipes
from os.path import join

def quote_param(items):
    return ' '.join( pipes.quote(item.replace('\n', ' ')) for item in items )

def make_quote(string):
    string = string.replace('$','$$')
    return string

def require_dir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    assert os.path.isdir(dirname)


class Batch:
    def __init__(self, dirname, submit='%'):
        self.dirname = dirname
        self.submit = submit
        require_dir(dirname)
        require_dir(join(dirname, 'state'))
        
        self.lines = [ ]
        self.all = [ ]
        
    def require_dir(self, name):
        require_dir(join(self.dirname, name))
    
    def target(self, path, dep, *commands):
        """ path is a directory or prefix or directory/prefix
            command is a command to execute to produce it
        """

        state_prefix = join(
            'state',
            grace.filesystem_friendly_name(path) + '_'
        )
        
        state_name = state_prefix + hashlib.sha1(
            '\n'.join(item.strip() for item in commands)
        ).hexdigest()
        self.all.append(state_name)
        
        self.lines.extend([
           '',
           '%s : %s' % (state_name, ' '.join(dep)),
           '\t@rm -f %s%s' % (state_prefix, '?'*40),
        ] + [
           '\t%s' % make_quote(self.submit.replace('%',command)) for command in commands 
        ] + [
           '\t@touch %s' % state_name,
        ])
        
        return state_name
    
    def virtual_target(self, name, dep, *commands):
        """ Target not included in "all", and with a fixed name. """    
        self.lines.extend([
            '',
            '%s : %s' % (name, ' '.join(dep)),            
        ] + [
           '\t%s' % make_quote(self.submit.replace('%',command)) for command in commands 
        ])
            
    def close(self):
        f = open(join(self.dirname, 'Makefile'), 'wb')
        print >> f, 'all :', ' '.join(self.all)
        print >> f, '\t@echo All done.'
        for line in self.lines:
            print >> f, line
        f.close()

        
class Options: pass


BATCH_HELP = r"""
Usage:

    nesoni batch: <batch_dir> [options] \
        reference: <reference.fa/gbk/etc> [...] \
        [sample: <sample_name> \ 
            [reads: <reads.fa> [...] ] \
            [pairs: <reads1.fa> <reads2.fa>] \
            [interleaved: <interleaved.fa> [...] ] ] \
        [import: <working_dir> [...]] \
        [do-clip: <nesoni clip options>] \
        [do-shrimp: <nesoni samshrimp options>] \
        [do-consensus: <nesoni samconsensus options> OR no ] \
        [do-count: <nesoni samcount options>] \
        [do-test-counts: <nesoni test-counts parameters>] [...] \
        [do-report: <nesoni report options>]

Construct a Makefile that will run other nesoni commands. 

For a set of samples, the Makefile will run:
1. nesoni clip:
2. nesoni samshrimp:
3. nesoni samconsensus:

If desired:
4. nesoni count:
5. nesoni test-counts:

Finally:
6. nesoni report:

Options:

    --damp-run yes/no     - Process only the first 10,000 reads/read-pairs,
                            to test for any problems.
                            Default: no 

    --run NNN             - Run now with specified level of parallelism (make -j NNN)
                            Default: just create directories and Makefile

    --nesoni COMMAND      - Command line to invoke nesoni.
                            Default: %s
    
    --pypy-nesoni COMMAND - Command line to invoke nesoni modules 
                            that work with pypy.
                            Default: same as --nesoni 

    --submit PATTERN      - How to run a command.
                            Should wait until command completes.
                            %% is replaced by the actual command.
                            Example: "jobrunningscript %%"
                            Default: %%

    --input-prefix PREFIX - Prefix for input filenames.
                            eg sftp://user@machine/home/user/foo/bar/baz/
                            Default: use absolute path of filename

"""

def batch_main(args):
    options = Options()
    
    options.references = [ ]
    
    options.clip_options = [ ]
    options.shrimp_options = [ ]
    options.do_consensus = True
    options.consensus_options = [ ]
    options.samples = [ ]
    
    options.do_count = False
    options.count_options = [ ]
    options.tests = [ ]
    
    options.report_options = [ ]

    default_nesoni = sys.executable + ' ' + sys.argv[0]
    options.nesoni, args = grace.get_option_value(args, '--nesoni', str, default_nesoni)    
    options.pypy_nesoni, args = grace.get_option_value(args, '--pypy-nesoni', str, options.nesoni)

    options.prefix, args = grace.get_option_value(args,'--input-prefix', str, None)
    options.submit, args = grace.get_option_value(args,'--submit', str, '%')
    assert '%' in options.submit, 'Bad submit pattern'
    
    options.damp, args = grace.get_option_value(args, '--damp-run', grace.as_bool, False)
    
    options.run, args = grace.get_option_value(args, '--run', int, None)
    
    def absolutize(filename):
        if options.prefix is not None:
            return options.prefix + filename
        else:
            return io.abspath(filename)
    def path_param(filenames, damp=False):
        if damp:
           filenames = [ item+'~~first:10000' for item in filenames ]
    
        return ' '.join(absolutize(filename) for filename in filenames)
       
    def default(args):
        grace.expect_no_further_options(args)
        if len(args) != 1:
            print >> sys.stderr, BATCH_HELP % default_nesoni
            raise grace.Help_shown()
        options.dirname = args[0]

    def reference(args):
        grace.expect_no_further_options(args)
        options.references.extend(args)

    def do_clip(args):
        options.clip_options.extend(args)

    def do_shrimp(args):
        options.shrimp_options.extend(args)
    
    def do_consensus(args):
        if args == ['no']:
            options.do_consensus = False
        else:
            options.consensus_options.extend(args)

    def sample(args):
        sample = Options()
        sample.imported = False
        sample.reads = [ ]
        sample.pairs = [ ]
        sample.interleaved = [ ]
        options.samples.append(sample)        
        def default(args):
            assert len(args) == 1, 'Expected a sample name in "sample:"'
            sample.name = args[0]
        def reads(args):
            grace.expect_no_further_options(args)
            sample.reads.extend(args)
        def pairs(args):
            grace.expect_no_further_options(args)
            assert len(args) == 2, 'Expected exactly two files in "pairs:"'
            sample.pairs.append(args)
        def interleaved(args):
            grace.expect_no_further_options(args)
            sample.interleaved.extend(args)
        grace.execute(args, [reads,pairs,interleaved], default)
        assert sample.reads or sample.pairs or sample.interleaved, 'No reads for sample'

    def import_(args):
        grace.expect_no_further_options(args)
        for item in args:
            sample = Options()
            options.samples.append(sample)
            sample.imported = True
            sample.clip_dest = None
            sample.align_dest = absolutize(item)

    def do_count(args):
        options.do_count = True
        options.count_options.extend(args)

    def do_test_counts(args):
        assert len(args) > 1, 'Incorrect parameters for test-counts'
        test = Options()
        options.tests.append(test)
        test.args = args
    
    def do_report(args):
        options.report_options.extend(args)

    grace.execute(args, [
        reference,
        do_clip,
        do_shrimp,
        do_consensus,
        sample,
        import_,
        do_count,
        do_test_counts,
        do_report,
    ], default)

    if options.damp:
        options.dirname += '-damp'
    
    if options.tests:
        options.do_count = True
    
    batch = Batch(options.dirname, options.submit) 

    for sample in options.samples:
        if sample.imported: continue
    
        # CLIP ===========================================
        batch.require_dir('clip')
        
        sample.clip_dest = join('clip', sample.name)
        command = (
            options.pypy_nesoni + 
            ' clip: ' +
            sample.clip_dest
        )
        if options.clip_options:
            command += ' ' + quote_param(options.clip_options)         
        if sample.reads:
            command += ' reads: ' + path_param(sample.reads, options.damp)
        for pair in sample.pairs:
            command += ' pairs: ' + path_param(pair, options.damp)
        if sample.interleaved:
            command += ' interleaved: ' + path_param(sample.interleaved, options.damp)
        
        sample.has_pairs = bool(sample.pairs) or bool(sample.interleaved)
        sample.clip_state = batch.target(
            sample.clip_dest,
            [],
            command
        )

        # ALIGN ==========================================
        batch.require_dir('align')
        
        sample.align_dest = join('align', sample.name) 
        command = options.pypy_nesoni + ' samshrimp: ' + sample.align_dest
        command += ' ' + path_param(options.references)
        command += ' reads: ' + sample.clip_dest + '_single.fq.gz'
        if sample.has_pairs:
            command += ' interleaved: ' + sample.clip_dest + '_paired.fq.gz'
        command += ' ' + quote_param(options.shrimp_options) 
        
        sample.align_state = batch.target(
            sample.align_dest,
            [ sample.clip_state ],
            command
        )
        
        # CONSENSUS =======================================
        if options.do_consensus:
            command = (
                options.pypy_nesoni + 
                ' samconsensus: ' + 
                sample.align_dest +
                ' ' + quote_param(options.consensus_options)
            )
            
            sample.consensus_state = batch.target(
                join(sample.align_dest, 'consensus'),
                [ sample.align_state ],
                command 
            )

        
    batch.virtual_target(
        'clip',
        [ sample.clip_state for sample in options.samples if not sample.imported ]
    )
    batch.virtual_target(
        'align',
        [ sample.align_state for sample in options.samples if not sample.imported ]
    )
    
    # COUNT ==========================================
    if options.do_count:
        command = options.pypy_nesoni + ' samcount: counts ' + quote_param(options.count_options)
        command += ' ' + ' '.join( sample.align_dest for sample in options.samples )
        
        options.counts_state = batch.target(
            'count',
            [ (sample.consensus_state if options.do_consensus else sample.align_state) for sample in options.samples if not sample.imported ],
                # count: --filter existing can depend on consensus
            command
        )
        batch.virtual_target('count', [ options.counts_state ])
        
        command = options.pypy_nesoni + ' plot-counts: scatter-plots counts.txt'
        options.plot_state = batch.target(
            'plot',
            [ options.counts_state ],
            command
        )
        

    # TEST ============================================
    for test in options.tests:
        batch.require_dir('test')
        test.dest = join('test', test.args[0])
        param = test.args[1:]
        
        command = options.pypy_nesoni + ' test-counts: ' + test.dest + ' counts.txt'
        command += ' ' + quote_param(param)
        
        if options.damp: 
            command += ' --min-count 1'
        
        test.state = batch.target(
            test.dest,
            [ options.counts_state ],
            command
        )
    
    if options.tests:
        command1 = 'rm -f differential-expression-tests.zip'
        command2 = (
            'zip -j differential-expression-tests.zip ' +
            ' '.join( test.dest + '*' for test in options.tests )
        )
        
        options.edger_zip_state = batch.target(
            'differential-expression-tests',
            [ test.state for test in options.tests ],
            command1,
            command2, 
        )
        
    
    # REPORT ===========================================

    command = options.pypy_nesoni + ' report: report ' + quote_param(options.report_options)
    command += ' reference: ' + path_param(options.references)
    command += ' clips: ' + ' '.join( sample.clip_dest+'_log.txt' for sample in options.samples if sample.clip_dest is not None )
    if options.do_consensus:
        command += ' aligns: ' + ' '.join( sample.align_dest for sample in options.samples )
    if options.do_count:
        command += ' count-log: counts_log.txt'

    if options.do_count:
        command += ' file: counts.txt \'Table of raw counts, RPKMs, and statistics on alignments spanning multiple genes.\''
        command += ' file: scatter-plots-count.png \'Pairwise scatter plots of number of reads aligning to each gene.\''
        command += ' file: scatter-plots-RPKM.png \'Pairwise scatter plots of RPKM values.\''
    
    if options.tests:
        command += ' file: differential-expression-tests.zip \'Differential gene expression analysis\''     
    
    options.report_state = batch.target(
        'report',
        batch.all[:], #Meh
        command
    )
    batch.virtual_target('report', [ options.report_state ])
    batch.virtual_target(
        'view', 
        [ options.report_state ],
        'firefox -no-remote report/index.html'
    )
        
    batch.close()
    
    if options.run is None:
        print
        print 'Now type:'
        print
        print 'make -C %s' % pipes.quote(options.dirname)
        print
    else:
        command = 'make -C %s -j %d' % (pipes.quote(options.dirname), options.run)
        print
        print command
        print
        assert 0 == os.system(command)





