
import grace, config

import sys

VERSION='0.73'

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
                    annotations, and files for SHRiMP.

    shrimp:       - Run SHRiMP 2 on a read set to set up a working
                    directory.
    
    consensus     - Filter read hits, and try to call a consensus for each 
                    position in reference.
    
    (import:      - Pipe SAM alignments to set up a working directory)
    (filter:      - Filter read hits, but do not call consensus)
    (reconsensus: - Re-call consensus, using previously filtered hits)

%(BOLD)sAlignment to reference -- analysis tools:%(END)s

    igv-plots:    - Generate plots for IGV.

    nway:         - Compare results of two or more runs of nesoni consensus,
                    amongst themselves and optionally with the reference.
                 
                    Can produce output suitable for phylogenetic analysis
                    in SplitsTree4.
        
    fisher:       - Compare results of two runs of nesoni consensus using
                    Fisher's Exact Test for each site in the reference.

    normalize:    - Create normalized Artemis depth plots.
                    See also "igv-plots:".

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
                  - Calculate normalizing multipliers from counts.
    
    heatmap:      - Draw a heat map of counts.
    
    nmf:          - Perform a Non-negative Matrix Factorization of counts.
                    NMF is a type of fuzzy clustering.
    
    compare-tests:
                  - Compare the output from two runs of "test-counts:"
                    eg to compare the results of different "--mode"s

An R+ module is included with nesoni which will help load the output from
samcount, for analysis with BioConductor packages.

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
    
    as-gff:       - Output an annotation in GFF format,
                    optionally filtering by annotation type.
    
    as-userplots: - Convert a .igv file to a set of .userplot files
                    for viewing in Artemis.
    
    make-genome:  - Make an IGV .genome file.
    
    sample:       - Randomly sample from a sequence file.
    
    stats:        - Show some statistics about a sequence or annotation file.

    fill-scaffolds: - Guess what might be in the gaps in a 454 scaffold.
    
    pastiche:     - Use MUMMER to plaster a set of contigs over reference 
                    sequences.
                    
%(BOLD)sInput files:%(END)s
- sequence files can be in FASTA, FASTQ, or GENBANK format.
- annotation files can be in GENBANK or GFF format (GFF is not yet supported by all tools).
- nesoni is able to read files compressed with gzip or bzip2.
- remote files can be specified as a URL (will be streamed using lftp).

""" % { 'BOLD' : '\x1b[1m', 'END' : '\x1b[m', 'VERSION' : VERSION }

from reference_directory import Make_reference
from clip import Clip
from samimport import Import
from samshrimp import Shrimp
from samconsensus import Filter, Reconsensus, Consensus
from samcount import Count
from runr import Test_counts, Heatmap, Compare_tests, Norm_from_counts, NMF
from trivia import As_fasta, As_gff, Sample, Stats
from igv import Make_genome, IGV_plots, As_userplots 
from workflows import Analyse_sample

from legion import *

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

  from nesoni import *
  
  Tool_name(
      arg1, arg2, ..., 
      flag_name=value,
      section_name=[values],...
  ).run()


More advanced usage would make use of the facilities in nesoni.legion:
  
  from nesoni import *
  
  def main():
      process_make( Tool_1(...) )
      process_make( Tool_2(...) )
      barrier()
      make( Tool_3(...) )
            
  if __name__ == '__main__':
      run_script(main)

"make" only starts re-running tools when it reaches a tool whos parameters have changed. 
A directory ".state" in the current directory is created to store the parameters from 
previous runs. You can also delete the appropriate file from .state to force re-running.

Allowing tools to run in parallel eliminates the implicit dependency between them.
barrier() waits for processes to finish, and declares that what follows depends 
on the output of these processes.

Here Tool_1 and Tool_2 run in parallel, and are forced to complete by barrier(),
then Tool_3 runs. Tool_3 is taken to depend on the output of Tool_1 and Tool_2
since it must execute after them, and will be re-run if Tool_1 or Tool_2 need to
be re-run.

A fully and correctly parallelized script will also have the correct dependency 
structure for recomputation when parameters are changed.

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
    def nway(args):
        grace.load('nway_diff').main(args)
    
    @add
    def fisher(args):
        grace.load('fisher_diff').main(args)
        
    @add
    def core(args):
        grace.load('core').main(args)
    
    @add
    def bag(args):
        grace.load('kmer').bag_main(args)
    
    @add
    def graph(args):
        grace.load('kmer').graph_main(args)
    
    @add
    def clean(args):
        grace.load('kmer').clean_main(args)
    
    @add
    def shred(args):
        grace.load('shred').main(args)

    #@add
    #def samshrimp(args):
    #    grace.load('samshrimp').samshrimp_main(args)
    #    #config.shell_run( grace.load('samshrimp').Samshrimp(), args, 'nesoni samshrimp:')
    
    commands['samimport'] = commands['import']
    commands['samshrimp'] = commands['shrimp']
    commands['samfilter'] = commands['filter']
    commands['samconsensus'] = commands['consensus']
    commands['samreconsensus'] = commands['reconsensus']
    commands['samcount'] = commands['count']
        
    #@add
    #def samfilter(args):
    #    grace.load('samconsensus').filter_main(args)
    #
    #@add
    #def samconsensus(args):
    #    grace.load('samconsensus').consensus_main(args, True)
    #
    #@add
    #def samreconsensus(args):
    #    grace.load('samconsensus').consensus_main(args, False)
    #
    #@add
    #def samcount(args):
    #    grace.load('samcount').count_main(args)
        
    @add
    def fill_scaffolds(args):
        grace.load('fill_scaffolds').fill_scaffolds(args)

    @add
    def pastiche(args):
        grace.load('pastiche').pastiche(args)

    #@add
    #def clip(args):
    #    grace.load('clip').clip(args)
    
    @add
    def plot_counts(args):
        grace.load('runr').plot_counts_main(args)

    #@add
    #def test_counts(args):
    #    grace.load('runr').test_counts_main(args)

    @add
    def test_power(args):
        grace.load('runr').test_power_main(args)
        
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
        
    @add
    def report(args):
        grace.load('report').report_main(args)
    
    #@add
    #def batch(args):
    #    grace.load('batch').batch_main(args)

    return commands



def main(args):
    if not args:
        config.write_colored_text(sys.stdout, USAGE)
        return 1

    commands = get_commands()
    
    command, args = args[0], args[1:]
    
    mangled_command = command.lower().rstrip(':')
    if mangled_command not in commands:
        raise grace.Error("Don't know how to "+command)
    
    commands[mangled_command](args)
    return 0
    



