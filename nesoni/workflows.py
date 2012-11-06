
import os

import nesoni

from nesoni import config, clip, samshrimp, bowtie, samconsensus, working_directory, workspace

@config.help('Analyse reads from a single sample.')
@config.Positional('reference', 'Reference directory created by "make-reference:".')
@config.Section('tags', 'Tags for this sample. (See "nesoni tag:".)')
@config.Section('reads', 'Files containing unpaired reads.')
@config.Section('interleaved', 'Files containing interleaved read pairs.')
@config.Grouped_section('pairs', 'Pair of files containing read pairs.')
@config.Configurable_section('clip', 'Options for clip', presets=[
    ('clip', clip.Clip(), 'Clip reads'),
    ('none', None, 'Don\'t clip reads'),
])
@config.Configurable_section('align', 'Options for aligner', presets=[
    ('shrimp', samshrimp.Shrimp(), 'Align using SHRiMP 2'),
    ('bowtie', bowtie.Bowtie(), 'Align using Bowtie 2'),
])
@config.Configurable_section('filter', 'Options for filter')
@config.Configurable_section('reconsensus', 'Options for reconsensus', presets=[
    ('reconsensus', samconsensus.Reconsensus(), 'Do consensus call'),
    ('none',        None,                       'Do not do consensus call'),
])
class Analyse_sample(config.Action_with_output_dir):
    reference = None
    tags = [ ]
    reads = [ ]
    interleaved = [ ]
    pairs = [ ]

    clip = clip.Clip()
    align = samshrimp.Shrimp()
    filter = samconsensus.Filter()
    reconsensus = samconsensus.Reconsensus()
    
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
class Analyse_variants(config.Action_with_output_dir):
    reference = None
    samples = [ ]

    def run(self):
        assert self.reference is not None, 'No reference directory given.'
        space = self.get_workspace()
        
        nesoni.Freebayes(
            space / 'variants-raw',
            samples=self.samples,
            ).make()
        
        nesoni.Vcf_filter(
            space / 'variants-filtered',
            space / 'variants-raw.vcf',
            ).make()
        
        nesoni.Snpeff(
            space / 'variants-filtered-annotated',
            self.reference,
            space / 'variants-filtered.vcf'
            ).make()
        
        nesoni.Vcf_patch(
            space / 'patched',
            self.reference,
            space / 'variants-filtered.vcf'
            ).make()


@config.help("""\
Analyse multiple samples. \
"analyse-sample:" is performed on a set of different samples, \
then "analyse-variants:".

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
    presets = [ ('default', Analyse_sample(), 'default "analyse-sample:" options') ],
    )
@config.Grouped_configurable_section('sample', 
    'Sample for analysis. Give one "sample:" section for each sample. '
    'Give a name for the sample and reads or read pairs as in "analyse-sample:". '
    'Also, any options in "analyse-sample:" may be overridden. See example below.',
    template_getter=lambda obj: obj.template
    )
class Analyse_samples(config.Action_with_output_dir):
    reference = None
    template = Analyse_sample()
    sample = [ ]

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
        samples = self.get_samples()
        work = self.get_workspace()

        with nesoni.Stage() as stage:
            for sample in samples:
                sample.process_make(stage)

        Analyse_variants(
            work/'variants',
            self.reference,
            samples=[ sample.output_dir for sample in samples ]
            ).make()




