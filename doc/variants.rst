
Variant calling for a set of samples
====================================

This is how to produce:

* A VCF file containing a list of variants in a set of samples.
* Various filtered VCF files and summaries.
* Phylogenetic trees or nets.
* A modified version of the reference sequence with the variants seen in a sample patched in.

Variant calling is actually performed automatically 
by the ``analyse-samples:`` tool we saw in :doc:`rnaseq`,
the idea being that there's no reason not to, and it might pick up a bodgy sample.
However here we'll see how to run this pipeline manually, 
and the steps and assumptions in it.

Preliminaries
-------------

This recipe assumes you have a set of working directories,
either created using the ``analyse-samples:`` command we saw in
:doc:`rnaseq`,
or individually by ``analyse-sample:`` as in :doc:`analyse-sample`.

By default the ``analyse-samples:`` command actuall already does variant calling!

Variant calling in a single command
-----------------------------------

::

  nesoni analyse-variants: myvariants sample1 sample2 sample3

Variant calling as individual steps
-----------------------------------

The first step is to use the FreeBayes program developed by the Marth Lab 
to identify variants.
By default Nesoni invokes FreeBayes with options that casts a wide net, 
because the next step will be to filter what it finds using our own statistical method.
This includes by default specifying an high ploidy, 4, to encourage FreeBayes to
find more alelles.
::

   nesoni freebayes: rawvariants sample1 sample2 sample3

In default operation, the next step reduces the ploidy down to one,
substitutes in a much more paranoid, and we believe realistic, genotype quality value,
then filters on the basis of this genotype quality.
::

   nesoni vcf-filter: filteredvariants rawvariants.vcf

This new genotype quality score is based on the following model: 
Alleles are drawn from an unknown mixture.
Our belief about this mixture is represented as a Dirichlet distribution,
the parameters of which are updated from some initial belief
by observing the allele seen in each relevent read in turn.
The genotype quality for a particular allele is then the probability that 
it makes up at least 50% (parameter ``--dirichlet-majority``) of the mixture.

We can then annotate the variants using the snpEFF program.
This requires specifying a reference directory (see :doc:`analyse-sample`).
::

   nesoni snpeff: annotatedvariants myref filteredvariants.vcf

We can also patch the observed allele into each sample
to produced "fixed" sequences.
::

   nesoni vcf-patch: patched filteredvariants.vcf

Variant calling in a single command, with more options
------------------------------------------------------

The options available for ``freebayes:``, ``vcf-filter:``, and ``snpeff:``
can all be given in subsections of the ``analyse-variants:`` command.

Slicing and dicing variants
---------------------------

The tool ``vcf-nway:`` provides a way to select subsets of variants,
for example based on snpEFF annotations,
or subsets of samples.
It can also output them in a variety of formats,
and even invoke SplitsTree4 to construct a phylogenetic net.
