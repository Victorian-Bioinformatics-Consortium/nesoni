
import os, glob

import nesoni

from nesoni import config, clip, samshrimp, bowtie, samconsensus, working_directory, reference_directory, \
                   workspace, io, reporting

class Context(object): pass

@config.help('Analyse reads from a single sample.')
@config.Positional('reference', 'Reference directory created by "make-reference:".')
@config.Section('tags', 'Tags for this sample. (See "nesoni tag:".)')
@config.Section('reads', 'Files containing unpaired reads.')
@config.Section('interleaved', 'Files containing interleaved read pairs.')
@config.Grouped_section('pairs', 'Pair of files containing read pairs.')
@config.Configurable_section('clip', 'Options for clip', presets=[
    ('clip', lambda obj: clip.Clip(), 'Clip reads'),
    ('none', lambda obj: None, 'Don\'t clip reads'),
    ])
@config.Configurable_section('align', 'Options for aligner', presets=[
    ('shrimp', lambda obj: nesoni.Shrimp(), 'Align using SHRiMP 2'),
    ('bowtie', lambda obj: nesoni.Bowtie(), 'Align using Bowtie 2'),
    ])
@config.Configurable_section('filter', 'Options for filter', presets=[
    ('filter', lambda obj: nesoni.Filter(), ''),
    ])
@config.Configurable_section('reconsensus', 'Options for reconsensus', presets=[
    ('reconsensus', lambda obj: nesoni.Reconsensus(), 'Do consensus call'),
    ('none',        lambda obj: None,                       'Do not do consensus call'),
    ])
@config.Configurable_section('count', 'Options for count, if doing RNA-seq.', presets=[
    ('none',    lambda obj: None,                           'Not RNA-seq, do not count expression levels.'),
    ('pool',    lambda obj: nesoni.Count(strand='pool'),    'Pool forward and reverse strands.'),
    ('forward', lambda obj: nesoni.Count(strand='forward'), 'Count on forward strand.'),
    ('reverse', lambda obj: nesoni.Count(strand='reverse'), 'Count on reverse strand.'),
    ('both',    lambda obj: nesoni.Count(strand='both'),    'Count both strands.'),
    ])
class Analyse_sample(config.Action_with_output_dir):
    reference = None
    tags = [ ]
    reads = [ ]
    interleaved = [ ]
    pairs = [ ]

    #clip = clip.Clip()
    #align = samshrimp.Shrimp()
    #filter = samconsensus.Filter()
    #reconsensus = samconsensus.Reconsensus()
    #count = samcount.Count()
    
    _workspace_class = working_directory.Working
    
    def get_context(self):
        context = Context()
        context.action = self
        context.space = self.get_workspace()

        reads = self.reads
        interleaved = self.interleaved
        pairs = self.pairs
        
        if not self.clip:
            context.clip = None
        else:
            context.clip = self.clip(
                context.space/context.space.name,
                reads = reads,
                interleaved = interleaved,
                pairs = pairs,
            )
            reads = context.clip.reads_output_filenames()
            pairs = context.clip.pairs_output_filenames()
            interleaved = context.clip.interleaved_output_filenames()
        
        context.align = self.align(
            self.output_dir,
            self.reference,
            reads = reads,
            interleaved = interleaved,
            pairs = pairs, 
            )
        
        if not self.filter:
            context.filter = None
        else:
            context.filter = self.filter(self.output_dir)
        
        if not self.reconsensus:
            context.reconsensus = None
        else:
            context.reconsensus = self.reconsensus(self.output_dir)
        
        if not self.count:
            context.count = None
        else:
            context.count = self.count(os.path.join(self.output_dir,'counts'), filenames = self.count.filenames + [ self.output_dir ])
        
        return context
        
    
    def run(self):
        context = self.get_context()
        
        if context.clip: 
            context.clip.make()
        
        context.align.make()
        
        if context.filter: 
            context.filter.make()
        
        with nesoni.Stage() as stage:
            if context.reconsensus: 
                context.reconsensus.process_make(stage)        
        
            nesoni.Tag(self.output_dir,tags=self.tags).make()

            if context.count:
                context.count.make()


@config.Positional('reference',
    'Reference directory created with "make-reference:".'
    )
@config.Main_section('samples', 
    'Working directories.'
    )
@config.Configurable_section('freebayes', 'Options for "freebayes:".', presets=[
    ('freebayes', lambda obj: nesoni.Freebayes(), ''),
    ])
@config.Configurable_section('vcf_filter', 'Options for "vcf-filter:".', presets=[
    ('vcf-filter', lambda obj: nesoni.Vcf_filter(), ''),
    ])
@config.Configurable_section('snpeff', 'Options for "snpeff:".', presets=[
    ('snpeff', lambda obj: nesoni.Snpeff(), ''),
    ('none', lambda obj: None, 'Do not use snpeff'),
    ])
@config.Configurable_section('analysis', '"analyse-sample:" option to use with "power-variant-call:". Does not affect actual results.', presets=[
    ('shrimp', lambda obj: Analyse_sample(align='shrimp'), 'Align using SHRiMP 2'),
    ('bowtie', lambda obj: Analyse_sample(align='bowtie'), 'Align using Bowtie 2'),
    ('none', lambda obj: None, 'Don\'t test power'),
    ])
class Analyse_variants(config.Action_with_output_dir):
    reference = None
    samples = [ ]
    
    #freebayes = 
    #vcf_filter = 
    #snpeff = 
    #align = 

    def run(self):
        assert self.reference is not None, 'No reference directory given.'
        space = self.get_workspace()
        
        if self.analysis:
            nesoni.Power_variant_call(
                space/'power',
                template__analysis   = self.analysis,
                template__freebayes  = self.freebayes,
                template__vcf_filter = self.vcf_filter,
                legacy = False,
                ).make()
        
        self.freebayes(
            space / 'variants-raw',
            samples=self.samples,
            ).make()
        
        self.vcf_filter(
            space / 'variants-filtered',
            space / 'variants-raw.vcf',
            ).make()        
        filename = space/'variants-filtered.vcf'
        
        if self.snpeff:
            self.snpeff(
                space / 'variants-filtered-annotated',
                self.reference,
                space / 'variants-filtered.vcf'
                ).make()
            filename = space / 'variants-filtered-annotated.vcf'
        
        io.symbolic_link(source=filename, link_name=space / 'variants.vcf')
        if os.path.exists(filename+'.idx'):
            io.symbolic_link(source=filename+'.idx', link_name=space / 'variants.vcf.idx')
        
        nesoni.Vcf_patch(
            space / 'patched',
            self.reference,
            space / 'variants.vcf'
            ).make()
        
        nesoni.Vcf_nway(
            space / 'net',
            space / 'variants.vcf',
            require='all',
            as_='splitstree',
            ).make()

        reporter = reporting.Reporter(space / 'report', 'Variants analysis')
                
        reporter.report_logs(None, [ space / 'variants-filtered_log.txt' ], 
            renaming = {'input':'Found by freebayes', 'kept':'Kept after quality filtering'})
        
        reporter.p(reporter.get(filename))
        if os.path.exists(filename+'.idx'):
            reporter.p(reporter.get(filename + '.idx') + ' (needed to view VCF file in IGV)')
        
        reporter.p(reporter.get(space / 'net.svg', title='Phylogenetic net'))
        
        if self.analysis:
            reporter.p(reporter.get(space / 'power_log.txt', title='Power report') +
                       '<br/>(Test of the ability of the pipeline to call various variants at various depths of coverage and in the presence of errors, using synthetic reads.)'
                       )
        
        reporter.close()
        


#@config.Positional('reference',
#    'Reference directory created with "make-reference:".'
#    )
@config.Main_section('samples', 
    'Working directories.'
    )
#@config.Configurable_section('count',
#    'Options for "count:".',
#    presets=[('default',lambda obj: nesoni.Count(),'')],
#    )
@config.Configurable_section('norm_from_counts',
    'Options for "norm-from-counts:" depth of coverage normalization.',
    presets=[('default',lambda obj: nesoni.Norm_from_counts(),'')],
    )
@config.Grouped_configurable_section('heatmap',
    'One or more heatmaps to produce.',
    template_getter=lambda obj: nesoni.Heatmap(),
    )
@config.Grouped_configurable_section('test',
    'One or more "test-counts:" to perform. '
    'Note: on the command line, '
    'give the counts file as ".", this is slightly broken and will be fixed in a future revision.',
    template_getter=lambda obj: nesoni.Test_counts(),
    )
class Analyse_expression(config.Action_with_output_dir):
    #reference = None
    samples = [ ]
    
    #count = 
    #norm_from_counts = 
    heatmap = [ ]
    test = [ ]
    
    def run(self):
        #assert self.reference is not None, 'No reference directory given.'
        space = self.get_workspace()
            
        #self.count(
        #    space / 'counts',
        #    filenames=self.samples,
        #    ).make()
        
        nesoni.Merge_counts(
            space / 'counts',
            filenames = [ os.path.join(item,'counts.csv') for item in self.samples ]
            ).make()
        
        self.norm_from_counts(
            space / 'norm',
            space / 'counts.csv'
            ).make()
        
        similarity = nesoni.Similarity(
            space / 'similarity',
            space / 'counts.csv',
            norm_file = space / 'norm.csv',
            )
        
        heatmaps = [
            heatmap(
                space / 'heatmap-'+heatmap.prefix,
                space / 'counts.csv',
                norm_file = space / 'norm.csv',
                )
            for heatmap in self.heatmap ]
        
        tests = [
            test(
                space / 'test-'+test.prefix,
                space / 'counts.csv',
                norm_file = space / 'norm.csv',
                )
            for test in self.test ]

        with nesoni.Stage() as stage:
            similarity.process_make(stage)
            for heatmap in heatmaps: 
                heatmap.process_make(stage)            
            for test in tests: 
                test.process_make(stage)       
        
        reporter = reporting.Reporter(space / 'report', 'Expression analysis')
        
        similarity.report(reporter)

        #reporter.heading('Sample similarity')
        #
        #reporter.p(
        #    'The following plots attempt to summarize the similarity/differences in expression patterns between samples, '
        #    'based on the glog2-transformed normalized read counts. '
        #    'Samples from the same experimental group should cluster together.'
        #    )
        #
        #reporter.p(
        #    reporter.get(space / 'similarity-plotMDS.png',
        #        title = 'limma\'s "plotMDS" Multi-Dimensional Scaling plot of sample similarity',
        #        image = True
        #        )
        #    )
        #
        #reporter.p(
        #    reporter.get(space / 'similarity.svg',
        #        title = 'Split Network visualization of sample similarity.',
        #        image = True
        #        ) +
        #    '<br>(Visualization of euclidean distances as a split network. '
        #    'Note: This is <i>not</i> a phylogenetic network.)'
        #    )
                
        if heatmaps:
            reporter.heading('Heatmaps')
            for heatmap in heatmaps:
                reporter.report_heatmap(heatmap)
        
        if tests:
            reporter.heading('Differential expression analysis')
            for test in tests:
                #reporter.report_test(test)
                test.report(reporter)
        
        reporter.heading('Raw data')
        
        reporter.p(reporter.get(space / 'counts.csv'))
        reporter.p(reporter.get(space / 'norm.csv'))
        
        reporter.close()
        


@config.help("""\
Analyse multiple samples. \
"analyse-sample:" is performed on a set of different samples, \
then some or all of \
"analyse-variants:", "analyse-expression:", "igv-plots:".

If "analyse-expression:" is used, \
"igv-plots:" will use the same normalization as used in \
"analyse-expression:".

For your sanity \
I recommend invoking this from a Python script \
rather than the command line. \
However, if you insist, it can be done. \
Example:

  nesoni analyse-samples: my_analysis my_reference_dir \\
      template: align: bowtie \\
      sample: sam1 tags: wildtype reads: sam1reads.fq \\
      sample: sam2 tags: mutant   reads: sam2reads.fq
      
""")
@config.Positional('reference', 'Reference directory created by "make-reference:".')
@config.Configurable_section('template', 
    'Common options for each "sample:". '
    'Put this section before the actual "sample:" sections.',
    presets = [ ('default', lambda obj: Analyse_sample(), 'default "analyse-sample:" options') ],
    )
#@config.Grouped_configurable_section('sample', 
#    'Sample for analysis. Give one "sample:" section for each sample. '
#    'Give a name for the sample and reads or read pairs as in "analyse-sample:". '
#    'Also, any options in "analyse-sample:" may be overridden. See example below.',
#    template_getter=lambda obj: obj.template
#    )
@config.Configurable_section_list('samples',
    'Samples for analysis. Give one "sample:" section for each sample. '
    'Give a name for the sample and reads or read pairs as in "analyse-sample:". '
    'Also, any options in "analyse-sample:" may be overridden. See example below.',
    templates = [ ],
    sections = [ ('sample', lambda obj: obj.template, 'A sample, parameters as per "analyse-sample:".') ],
    )
@config.Configurable_section('variants',
    'Options for "analyse-variants:".',
    presets = [ 
        ('default', lambda obj: Analyse_variants(analysis=None), 'default "analyse-variants:" options'),
        ('none', lambda obj: None, 'don\'t analyse variants'),
        ],
    )
@config.Configurable_section('expression',
    'Options for "analyse-expression:".',
    presets = [ 
        ('none', lambda obj: None, 'don\'t analyse expressions'),
        ('default', lambda obj: Analyse_expression(), 'default "analyse-expression:" options'),
        ],
    )
@config.Configurable_section('igv_plots',
    presets = [
        ('default', lambda obj: nesoni.IGV_plots(), 'default "igv-plots:" options'),
        ('none', lambda obj: None, 'don\'t produce IGV plots'),
        ],
    )
@config.String_flag('report_title')
@config.Bool_flag('include_genome', 'Include .genome with IGV plots in report.')
@config.Bool_flag('include_bams', 'Include BAM files in report.')
class Analyse_samples(config.Action_with_output_dir):
    reference = None
    #template = Analyse_sample()
    samples = [ ]
    
    #variants = 
    #expression =
    #igv-plots =
    
    report_title = 'High-throughput sequencing analysis'
    include_genome = True
    include_bams = True

    def get_context(self):
        assert self.reference is not None, 'reference not given.'

        context = Context()
        context.action = self
        context.space = self.get_workspace()
        context.name = context.space.name
        context.sample_space = workspace.Workspace(context.space/'samples', False)
        context.reference = reference_directory.Reference(self.reference, True)

        for sample in self.samples:
            assert sample.reference is None, 'reference should not be specified within sample.'
            assert sample.output_dir, 'sample given without name.'
            assert os.path.sep not in sample.output_dir, 'sample name contains '+os.path.sep
        context.samples = [
            sample(
                output_dir = context.sample_space / sample.output_dir,
                reference = self.reference,                
                )
            for sample in self.samples
            ]
        context.sample_dirs = [ item.output_dir for item in context.samples ]
        
        if not self.variants:
            context.variants = None
        else:
            context.variants = self.variants(
                context.space/'variants',
                self.reference,
                samples=context.sample_dirs,
                analysis=self.variants.analysis or self.template, #Use aligner from template unless aligner explicitly given
                )

        if not self.expression:
            context.expression = None
        else:
            context.expression = self.expression(
                context.space/'expression',
                samples=context.sample_dirs,
                )
        
        return context

    def run(self):
        context = self.get_context()

        with nesoni.Stage() as stage:
            for sample in context.samples:
                sample.process_make(stage)
            
        with nesoni.Stage() as stage:
            if context.variants:
                context.variants.process_make(stage)
    
            if context.expression:
                context.expression.process_make(stage)

        if self.igv_plots:
            plot_space = workspace.Workspace(context.space/'plot',False)
            self.igv_plots(
                prefix = plot_space / ('plot'),
                genome = context.reference.get_genome_filename(),
                norm_file = context.space/('expression','norm.csv') if context.expression else None,
                working_dirs = context.sample_dirs,
                ).make()

        # =================================================================================
        # =================================================================================
        # =================================================================================

        reporter = reporting.Reporter(context.space / 'report', self.report_title, context.name)
        
        reporter.report_logs('alignment-statistics',
            [ sample.get_context().clip.log_filename() 
                for sample in context.samples
                if sample.clip
                ] +
            [ sample.get_context().filter.log_filename() 
               if not sample.count
               else sample.get_context().count.log_filename()
               for sample in context.samples 
               if sample.filter or sample.count
               ],
            filter=lambda sample,field: field != 'fragments',
            )
        
        if self.expression:
            io.symbolic_link(source=context.space/('expression','report'),link_name=context.space/('report','expression'))
            reporter.heading('<a href="expression/index.html">&gt; Expression analysis</a>')
        
        if self.variants:
            io.symbolic_link(source=context.space/('variants','report'),link_name=context.space/('report','variants'))
            reporter.heading('<a href="variants/index.html">&gt; Variants analysis</a>')
                
        if self.igv_plots:
            reporter.heading('IGV plots')            
            reporter.p('These files show the depth of coverage. They can be viewed with the IGV genome browser.')

            genome_files = [ ]
            if self.include_genome:
                genome_filename = context.reference.get_genome_filename()
                genome_dir = context.reference.get_genome_dir()
                genome_files.append(genome_filename)
                if genome_dir:
                    base = os.path.split(genome_dir)[1]
                    for filename in os.listdir(genome_dir):
                        genome_files.append((
                            os.path.join(genome_dir, filename),
                            os.path.join(base, filename)
                            ))
            
            reporter.p(reporter.tar('igv-plots',
                genome_files +
                glob.glob(plot_space/'*.tdf')
                ))

        if self.include_bams:
            reporter.heading('BAM files')
            
            reporter.p('These BAM files contain the alignments of reads to the reference sequences.'
                       ' They can also be viewed using IGV.')
            
            bam_files = [ ]
            for sample in self.samples:
                name = sample.output_dir
                bam_files.append( (context.space/('samples',name,'alignments_filtered_sorted.bam'),name+'.bam') )
                bam_files.append( (context.space/('samples',name,'alignments_filtered_sorted.bam.bai'),name+'.bam.bai') )
            reporter.p(reporter.tar('bam-files', bam_files))
        
        reporter.write('<p/><hr/>\n')
        reporter.p('nesoni version '+nesoni.VERSION)
        reporter.close()




