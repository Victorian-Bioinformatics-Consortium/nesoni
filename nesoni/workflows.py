
import os

import nesoni

from nesoni import config, clip, samshrimp, bowtie, samconsensus, working_directory, reference_directory, \
                   workspace, io, reporting

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
    ('shrimp', lambda obj: samshrimp.Shrimp(), 'Align using SHRiMP 2'),
    ('bowtie', lambda obj: bowtie.Bowtie(), 'Align using Bowtie 2'),
    ])
@config.Configurable_section('filter', 'Options for filter', presets=[
    ('filter', lambda obj: samconsensus.Filter(), ''),
    ])
@config.Configurable_section('reconsensus', 'Options for reconsensus', presets=[
    ('reconsensus', lambda obj: samconsensus.Reconsensus(), 'Do consensus call'),
    ('none',        lambda obj: None,                       'Do not do consensus call'),
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
    
    _workspace_class = working_directory.Working
    
    def run(self):
        reads = self.reads
        interleaved = self.interleaved
        pairs = self.pairs
        
        workspace = self.get_workspace()
        
        if self.clip:
            act = self.clip(
                workspace/workspace.name,
                reads = reads,
                interleaved = interleaved,
                pairs = pairs,
            )
            act.make()
            reads = act.reads_output_filenames()
            pairs = act.pairs_output_filenames()
            interleaved = act.interleaved_output_filenames()
        
        self.align(
            self.output_dir,
            self.reference,
            reads = reads,
            interleaved = interleaved,
            pairs = pairs, 
        ).make()
        
        if self.filter:
            self.filter(
                self.output_dir
            ).make()
        
        if self.reconsensus:
            self.reconsensus(
                self.output_dir
            ).make()
        
        nesoni.Tag(
            self.output_dir, 
            tags=self.tags
            ).make()


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
class Analyse_variants(config.Action_with_output_dir):
    reference = None
    samples = [ ]
    
    #freebayes = 
    #vcf_filter = 
    #snpeff = 

    def run(self):
        assert self.reference is not None, 'No reference directory given.'
        space = self.get_workspace()
        
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


#@config.Positional('reference',
#    'Reference directory created with "make-reference:".'
#    )
@config.Main_section('samples', 
    'Working directories.'
    )
@config.Configurable_section('count',
    'Options for "count:".',
    presets=[('default',lambda obj: nesoni.Count(),'')],
    )
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
            
        self.count(
            space / 'counts',
            filenames=self.samples,
            ).make()
        
        self.norm_from_counts(
            space / 'norm',
            space / 'counts.csv'
            ).make()
        
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
                norm_file = space / 'norm.csv'
                )
            for test in self.test ]
        
        with nesoni.Stage() as stage:
            for heatmap in heatmaps: 
                heatmap.process_make(stage)            
            for test in tests: 
                test.process_make(stage)        
        
        
        reporter = reporting.Reporter(space / 'report', 'Expression analysis')
        
        if heatmaps:
            reporter.heading('Heatmaps')
            for heatmap in heatmaps:
                reporter.report_heatmap(heatmap)
        
        if tests:
            reporter.heading('Differential expression analysis')
            for test in tests:
                reporter.report_test(test)
        
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
@config.Grouped_configurable_section('sample', 
    'Sample for analysis. Give one "sample:" section for each sample. '
    'Give a name for the sample and reads or read pairs as in "analyse-sample:". '
    'Also, any options in "analyse-sample:" may be overridden. See example below.',
    template_getter=lambda obj: obj.template
    )
@config.Configurable_section('variants',
    'Options for "analyse-variants:".',
    presets = [ 
        ('default', lambda obj: Analyse_variants(), 'default "analyse-variants:" options'),
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
class Analyse_samples(config.Action_with_output_dir):
    reference = None
    #template = Analyse_sample()
    sample = [ ]
    
    #variants = 
    #expression =
    #igv-plots =

    def get_sample_space(self):
       work = self.get_workspace()
       return workspace.Workspace(work/'samples',must_exist=False)

    def get_samples(self):
        """ Fill in given samples with reference. """
        assert self.reference is not None, 'reference not given.'
        for sample in self.sample:
            assert sample.reference is None, 'reference should not be specified within sample.'
            assert sample.output_dir, 'sample given without name.'
            assert os.path.sep not in sample.output_dir, 'sample name contains '+os.path.sep
        sample_space = self.get_sample_space()
        return [
            sample(
                output_dir = sample_space / sample.output_dir,
                reference = self.reference,                
            )
            for sample in self.sample
        ]

    def run(self):
        assert self.reference is not None, 'reference not given.'
        reference = reference_directory.Reference(self.reference, True)
        samples = self.get_samples()
        sample_dirs = [ sample.output_dir for sample in samples ]
        work = self.get_workspace()

        with nesoni.Stage() as stage:
            for sample in samples:
                sample.process_make(stage)
            
        with nesoni.Stage() as stage:
            if self.variants:
                self.variants(
                    work/'variants',
                    self.reference,
                    samples=sample_dirs,
                    ).process_make(stage)
    
            if self.expression:
                self.expression(
                    work/'expression',
                    samples=sample_dirs,
                    ).process_make(stage)

        if self.igv_plots:
            plot_space = workspace.Workspace(work/'plot',False)
            self.igv_plots(
                prefix = plot_space / ('plot'),
                genome = reference.get_genome_filename(),
                norm_file = work/('expression','norm.csv') if self.expression else None,
                working_dirs = sample_dirs,
                ).make()



