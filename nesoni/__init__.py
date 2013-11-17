VERSION='0.113'
#^ Note: this first line is read by the setup.py script to get the version

import sys

import grace, config

from reference_directory import Make_reference
from clip import Clip
from working_directory import Tag
from samimport import Import
from samshrimp import Shrimp
from bowtie import Bowtie
from samconsensus import Filter, Reconsensus, Consensus
from nway_diff import Nway
from fisher_diff import Fisher
from core import Core
from samcount import Count, Merge_counts
from runr import Test_counts, Plot_counts, Test_power, Heatmap, Similarity, Compare_tests, Norm_from_counts, Glog, NMF
from normalize import Norm_from_samples
from trivia import Test, As_fasta, As_gff, Sample, Stats
from shred import Shred
from igv import Make_genome, IGV_plots, As_userplots, Run_igv
from variant import Freebayes, Vcf_filter, Snpeff, Vcf_nway, Vcf_patch, Test_variant_call, Power_variant_call
from peaks import Islands,Transcripts,Modes
from annotation_tools import Modify_features, Collapse_features, Relate_features
from workflows import Analyse_sample, Analyse_variants, Analyse_expression, Analyse_samples
from changes import Changes

from legion import *

BOLD = '\x1b[1m'
END = '\x1b[m'

USAGE = """
%(BOLD)snesoni%(END)s %(VERSION)s - high-throughput sequencing data analysis toolset

Usage:

    %(BOLD)snesoni <tool>: %(END)s...

Give <tool>: without further arguments for help on using that tool.


%(BOLD)sAlignment to reference -- core tools:%(END)s

    make-reference:
                  - Set up a directory containing a reference sequence,
                    annotations, and files for SHRiMP and/or Bowtie.

    shrimp:       - Run SHRiMP 2 on a read set to set up a working
                    directory.
    
    bowtie:       - Run Bowtie 2 on a read set to set up a working
                    directory.
    
    consensus     - Filter read hits, and try to call a consensus for each 
                    position in reference.
    
    (import:      - Pipe SAM alignments to set up a working directory)
    (filter:      - Filter read hits, but do not call consensus)
    (reconsensus: - Re-call consensus, using previously filtered hits)


%(BOLD)sAlignment to reference -- VCF based tools: (under development)%(END)s

These provide an alternative to consensus calling using "nesoni consensus:"
- better handling of complicated Multi-Nucleotide Polymorphisms
- can't distinguish between absence of a variant and insufficient data
  (but can distinguish absence of a variant from insufficient data 
   in a single sample if variant present in other samples)

    freebayes:    - Run FreeBayes to produce a VCF file.
    
    vcf-filter:   - Filter a VCF file, eg as produced by "nesoni freebayes:".
    
    snpeff:       - Run snpEff to annotate variants with their effects.
    
    vcf-nway:     - Summarize a VCF file in a variety of possible ways.
    
    vcf-patch:    - Patch in variants to produce genome of samples.
                    (similar to consensus_masked.fa produced by "nesoni consensus:")
    
    test-variant-call:
                  - Generate synthetic reads, see what variant is called.
    
    power-variant-call:
                  - Apply "neosni test-variant-call:" to a variety of
                    different variants over a range of depths.


%(BOLD)sAlignment to reference -- analysis tools:%(END)s

    igv-plots:    - Generate plots for IGV.

    nway:         - Compare results of two or more runs of nesoni consensus,
                    amongst themselves and optionally with the reference.
                 
                    Can produce output suitable for phylogenetic analysis
                    in SplitsTree4.
        
    fisher:       - Compare results of two runs of nesoni consensus using
                    Fisher's Exact Test for each site in the reference.

    core:         - Infer core genome present in a set of strains.

    (consequences: 
                  - Determine effects at the amino acid level of SNPs and INDELs
                    called by nesoni consensus. Most of the features of this tool
                    are now a part of "samconsensus:".)


%(BOLD)sAlignment to reference -- differential expression:%(END)s

    count:        - Count number of alignments to genes, using output from
                    "shrimp:".

    test-counts:  - Use edgeR or limma from BioConductor to detect differentially
                    expressed genes, using output from samcount.
    
    test-power:   - Test the statistical power of "nesoni test-counts:" with
                    simulated data.

    plot-counts:  - Plot counts against each other.
    
    norm-from-counts:
                  - Calculate normalizing multipliers from counts using TMM.
    
    norm-from-samples:
                  - Calculate normalizing multipliers from working directories.
    
    heatmap:      - Draw a heat map of counts.
    
    nmf:          - Perform a Non-negative Matrix Factorization of counts.
                    NMF is a type of fuzzy clustering.
    
    compare-tests:
                  - Compare the output from two runs of "test-counts:"
                    eg to compare the results of different "--mode"s
    
    similarity:
                  - Compare samples in a counts file in various ways.
    
    glog:
                  - Obtain glog2 RPM values from a counts file.

An R+ module is included with nesoni which will help load the output from
samcount, for analysis with BioConductor packages.


%(BOLD)sPeak calling and annotation manipulation tools:%(END)s

    islands:
    transcripts:
    modes:
                  - Various peak and transcript calling algorithms.
    
    modify-features:
                  - Shift start or end position of features,
                    filter by type, change type.
    
    collapse-features:
                  - Merge overlapping features.
    
    relate-features:
                  - Find features from one set that are near or overlapping
                    features from another set.
    
    as-gff:       - Output an annotation in GFF format,
                    optionally filtering by annotation type.

    
%(BOLD)sk-mer tools:%(END)s (experimental)

    bag:          - Create an index of kmers in a read set for analysis with
                    nesoni graph.
    
    graph:        - Use a bag or bags to lay out a deBruijn graph.
                    Interact with the graph in various ways.


%(BOLD)sUtility tools:%(END)s

    clip:         - Remove Illumina adaptor sequences and low quality bases
                    from reads.

    shred:        - Break a sequence into small overlapping pieces.
                    In case you want to run an existing sequence through the 
                    above tools. Yes, this isn't ideal.

    as-fasta:     - Output a sequence file in FASTA format.
    
    as-userplots: - Convert a .igv file to a set of .userplot files
                    for viewing in Artemis.
    
    make-genome:  - Make an IGV .genome file.
    
    run-igv:      - Run IGV with a specified .genome file.
    
    sample:       - Randomly sample from a sequence file.
    
    stats:        - Show some statistics about a sequence or annotation file.

    fill-scaffolds: - Guess what might be in the gaps in a 454 scaffold.
    
    pastiche:     - Use MUMMER to plaster a set of contigs over reference 
                    sequences.

    changes:      - Prints out change log file.
                    

%(BOLD)sPipeline tools:%(END)s
    
    analyse-sample:
                  - Clip, align, and call consensus on a set of reads.
                  
    analyse-variants:
                  - Produce a VCF file listing SNPs and other variants in
                    a set of samples.
    
    analyse-expression:
                  - Count alignments of fragments to genes,
                    then perform various types of statistics and 
                    visualization on this.

    analyse-samples:
                  - Run "analyse-sample:" on a set of different samples,
                    then run "analyse-variants:" and/or "analyse-expression".

If a pipeline tool is run again, it restarts only from the point affected 
by the changed parameters. The following global flags control pipeline tool 
behaviour:
%(MAKE)s

                    
%(BOLD)sInput files:%(END)s
- sequence files can be in FASTA, FASTQ, or GENBANK format.
- annotation files can be in GENBANK or GFF format 
  (GFF is not yet supported by all tools).
- nesoni is able to read files compressed with gzip or bzip2.


%(BOLD)sSelections and sorts:%(END)s
Working directories can be given a set of tags using "tag:". They also 
implicitly have a tag for the name of the directory, and a tag "all".

A %(BOLD)sselection expression%(END)s is a logical expression used to select a subset 
of working directories. It may consist of (grouped by precedence):

  tag        - Working directories with tag

  [exp]      - exp

  -exp       - not exp

  exp1:exp2  - exp1 and exp2
  exp1/exp2  - exp1 or exp2
  exp1^expr2 - exp1 xor exp2

Example:

  [strain1:time1]/[-strain1:-time1]   
             - Samples either from strain1 at time1, 
               or not from strain1 and not from time1.
               Equivalently: strain1^-time1

A %(BOLD)ssort expression%(END)s is a comma separated list of selection expressions,
used to sort a list of working directories.

Example:

  strain1,strain2,time1,time2,time3,replicate1,replicate2
             - Sort, grouping by strain, then by time, then by replicate

""" % { 'BOLD' : '\x1b[1m', 'END' : '\x1b[m', 'VERSION' : VERSION, 'MAKE' : Make().describe('', show_help=True, escape_newlines=False) }

def get_actions():
    for item in globals().values():
        if isinstance(item, type) and issubclass(item, config.Action):            
            yield item

def get_all():
    return (
        legion.__all__ +
        [ key for key,value in globals().items()
          if isinstance(value, type) and issubclass(value, config.Action)
        ]
    )

__all__ = get_all()

__doc__ = """

Nesoni tools may be accessed via Python in addition to access 
via command line. The Python interface closely follows the command
line interface.

Simple usage would be:

  import nesoni
  
  if __name__ == '__main__':
      nesoni.Tool_name(
          arg1, arg2, ..., 
          flag_name=value,
          section_name=[values],...
          ).run()


More advanced usage would make use of the facilities in nesoni.legion:
  
  import nesoni
  
  def main():
      with nesoni.Stage() as stage:
          nesoni.Tool_1(...).process_make(stage)
          nesoni.Tool_2(...).process_make(stage)
          stage.barrier()
          nesoni.Tool_3(...).make()
            
  if __name__ == '__main__':
      nesoni.run_script(main)

"make" only starts re-running tools when it reaches a tool whos parameters have changed.
(Command line options --make-do and --make-done can be used to override this.)

Allowing tools to run in parallel eliminates the implicit dependency between them.
barrier() waits for processes to finish, and declares that what follows depends 
on the output of these processes.

Here Tool_1 and Tool_2 run in parallel, and are forced to complete by barrier(),
then Tool_3 runs. Tool_3 is taken to depend on the output of Tool_1 and Tool_2
since it must execute after them, and will be re-run if Tool_1 or Tool_2 need to
be re-run.

A fully and correctly parallelized script will also have the correct dependency 
structure for recomputation when parameters are changed.


The "run_script" function gives your script various command line options to control 
making.


The following tools are available:

""" + '\n'.join(sorted(
    item.__name__
    for item in get_actions()
))


def get_commands():
    commands = { }
    
    for item in get_actions():
        name = item.__name__.lower().replace('_','-')
        def func(args, action=item,name=name):
            config.shell_run(action(), args, 'nesoni %s:' % name)
        commands[name] = func
    
    def add(func):
        commands[ func.__name__.replace('_','-').strip('-') ] = func
            
    @add
    def consequences(args):
        grace.load('consequences').main(args)
        
    @add
    def bag(args):
        grace.load('kmer').bag_main(args)
    
    @add
    def graph(args):
        grace.load('kmer').graph_main(args)
    
    @add
    def clean(args):
        grace.load('kmer').clean_main(args)
    
    commands['samimport'] = commands['import']
    commands['samshrimp'] = commands['shrimp']
    commands['samfilter'] = commands['filter']
    commands['samconsensus'] = commands['consensus']
    commands['samreconsensus'] = commands['reconsensus']
    commands['samcount'] = commands['count']
        
    @add
    def fill_scaffolds(args):
        grace.load('fill_scaffolds').fill_scaffolds(args)

    @add
    def pastiche(args):
        grace.load('pastiche').pastiche(args)

    @add
    def plot(args):
        grace.load('trivia').plot(args)

    @add
    def normalize(args):
        grace.load('trivia').normalize(args)
    
    @add
    def debias(args):
        grace.load('trivia').debias(args)
        
    @add
    def recombination(args):
        grace.load('recombination').recombination(args)

    return commands



def main():
    try:
        args = sys.argv[1:]
        args = configure_making(args)
    
        if not args:
            config.write_colored_text(sys.stdout, USAGE)
            grace.check_installation()
            sys.exit(1)
    
        commands = get_commands()
        
        command, args = args[0], args[1:]
        
        mangled_command = command.lower().rstrip(':')
        if mangled_command not in commands:
            raise grace.Error("Don't know how to "+command)
        
        commands[mangled_command](args)

    except grace.Help_shown:
        sys.exit(1)
    
    except Exception:
        config.report_exception()
        sys.exit(1)     
    


