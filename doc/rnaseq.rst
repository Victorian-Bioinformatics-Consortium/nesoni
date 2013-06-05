
Analyse gene expression patterns of a set of RNA-seq samples
============================================================

This is how to produce:

* Everything from :doc:`analyse-sample` for each of a set of samples.
* Depth of coverage plots for the IGV browser.
* A web page summarizing the results.

This page assumes you already have a reference directory, see :ref:`make-reference`.


Basic RNA-seq analysis in one command
-------------------------------------

Example:
::

  nesoni analyse-samples: myanalysis myref \
      template: count: forward \
      sample: foo1 tags: foo rep1 reads: foo1-reads.fa \
      sample: foo2 tags: foo rep2 reads: foo2-reads.fa \
      sample: bar1 tags: bar rep1 reads: bar1-reads.fa \
      sample: bar2 tags: bar rep2 reads: bar2-reads.fa \
      expression: default

Basic RNA-seq analysis as individual steps
------------------------------------------

The first step is to generate a working directory for each of the samples,
as described in :doc:`analyse-sample`.

It will be useful at this stage to give each sample a set of tags, using the ``tag:`` tool.
For example, we might want to tag sample foo1 with the tags foo and rep[licate]1:
::

  nesoni tag: foo1 foo rep1

We now need to count the number of reads (or fragments if we have read pairs) aligning
to each gene in each sample.
This is achived with the ``count:`` command.
::

  nesoni count: mycounts foo1 foo2 bar1 bar2

This produces a file ``mycounts.csv``.

*Optional:* It's also possible at this stage to create a normalization file.
::

  nesoni norm-from-counts: mynorm mycounts.csv
  
The resultant normalization file ``mynorm.csv`` can be passed to 
downstream tools using the ``--norm`` flag.
If this step is omitted,
the normalization will be recalculated each time it is needed.

Basic RNA-seq analysis in one command, with more options
--------------------------------------------------------

The above steps, and production of a default set of heatmaps, are encapsulated in 
a pipeline command called ``analyse-expression:``.

This ``analyse-expression:`` pipeline is itself an optional part of the
``analyse-samples:`` pipeline, 
all of its options and sections can be given after ``expression: default``
and an ``analyse-samples:`` command.
It's entirely possible that this is too clever by half.

One difference between this method and the method just described is that counting
is done for each sample individually, and these are then merely merged by 
``analyse-expression:`` pipeline (using the ``merge-counts:`` tool).
Therefore if using ``analyse-sample:`` with ``analyse-expression:`` be sure to include
a ``count:`` section, for example:
::

  nesoni analyse-sample: foo1 tags: foo rep1 reads: foo1.fq count: forward

``count: forward`` here means your reads are strand-specific and on the same strand as the gene.
Other options can be seen in the help for ``analyse-sample:`` and ``count:``.
For example, 
by default CDS features are counted, 
but this can be changed with the ``--types`` flag.

Heatmaps
--------

Heatmaps with hierarchical clustering provide a basic way to
assess the quality of your data, and the patterns of differential expression it contains.

For example, 
this will generate a heatmap of all genes
for which there is at least a 8-fold difference between two samples
(log2(8) = 3).
::
  
  nesoni heatmap: myheatmap mycounts.csv --min-span 3

The colors in the heatmap are based on log2 transformed normalized counts,
with moderation of values near zero 
(the log transforms of which are necessarily noisy).
The moderation amount is controlled by the parameter ``--glog-moderation``,
we believe it has a reasonable default value and shouldn't need adjusting.
  
Differential expression testing
-------------------------------

Nesoni provides a convenient interface to the limma and edgeR Bioconductor packages 
in a tool called ``test-counts:``.

A basic knowledge of linear model based hypothesis testing is assumed.
`Please, please, know what you are doing with this tool`.
Read the documentation displayed by ``nesoni test-counts:``,
read the documentation for limma and edgeR 
and work through some examples in R+ as well as using this tool.

Here we will test which genes in group bar differ from those in group foo.
::

  nesoni test-counts: mytest mycounts.csv bar

Here ``bar`` is a selection expression on the tags we gave to samples earlier,
resulting in a term for the linear model that is 1 for samples tagged as "bar" and 0 otherwise.
The syntax for selection expression is described in Nesoni's main help blurb, 
which is displayed by typing ``nesoni``.
We are comparing the linear model containing this term and a constant term to the linear
model containing just a constant term.
Don't use this tool unless and until this makes sense.

This uses the default expression testing method ``--mode voom`` which uses voom and limma.

A .csv file or results is produced, with an accompanying heatmap.
The heatmap can be used to visually assess how well the analysis method is treating 
within-group and between-group variation when calling significant differential expression.


Non-negative Matrix Factorization
---------------------------------

A further, experimental, perspective on your data is provided by the
Non-negative Matrix Factorization (NMF) tool.
::

  nesoni nmf: mynmf-5 mycounts.csv --rank 5

NMF is a fuzzy clustering technique.
Nesoni adds a further twist or two on top of standard NMF,
the supreme awesomeness of which will be described in a future paper.
In the meantime it may suggest things to test with other, less magical, tools.
Try it with a few different ranks (number of classes).


