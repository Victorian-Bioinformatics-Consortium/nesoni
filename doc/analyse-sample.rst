
Analyse a single sample
=======================

This is how to produce:

* BAM file containing read alignments.
* Depth of coverage plots for the Artemis genome browser.
* A list of SNPs and indels.
* A modified version of the reference sequence with the SNPs and indels patched into it.
* A sample "working directory" used by other Nesoni tools.


.. _make-reference:

First make a reference directory
--------------------------------

Our first step is to make a "reference directory" for nesoni,
which collects all the information about the reference sequence that nesoni needs,
and any indexes needed by various tools.

The input is some sequences and, if available, annotations of those sequences.
This could be just FASTA sequences, 
a GENBANK file containing sequence and annotations, 
or FASTA sequences and a GFF file.

::

  nesoni make-reference: myref ref-sequences.fa

This will generate SHRiMP indicies if SHRiMP is installed
(and a Bowtie2 index if Bowtie2 is installed),
a .genome file for the IGV browser,
and files needed to use snpEFF if snpEFF is installed.

Analysis in one command
-----------------------

::

  nesoni analyse-sample: mysample myref reads: reads.fastq

This clips reads, aligns, filters alignments, and calls SNPs and indels, all in one command!

If you have read-pairs, you can use a ``pairs: reads1.fastq reads2.fastq`` section
instead of the ``reads:`` section.

Depth of coverage plots viewable in Artemis are also produced.
We'll see how to create depth of coverage plots for IGV 
when we look at analysing multiple samples.

Analysis as individual steps
----------------------------

Let's break it down into individual steps.
Then you can specify options for each of those steps.

First clip reads based on quality
and (by default) Illumina adaptor sequences.
::

  nesoni clip: myclipped reads: reads.fastq

Next align reads to the reference.
This creates what Nesoni refers to as a `working directory`.
::

  nesoni shrimp: mysample myref reads: myclipped-single.fq.gz

``bowtie:`` can be used instead of ``shrimp:``, 
or a BAM file from another read aligner can be imported with ``import:``.

Next, filter the read alignments and call SNPs and indels:
::

  nesoni filter: mysample
  nesoni reconsensus: mysample

``consensus:`` can be used as a shorthand for ``filter:`` followed by ``reconsensus:``.

The default filtering mode is to only retain reads that map to a single location.
You can retain all possibly valid alignments with ``--monogamous no``, 
or choose randomly amongst possibly valid alignments with
``--monogamous no --random yes``.

Read pairs
----------

Nesoni tools that take reads as inputs generally have three possible sections
for giving reads: 
``reads:`` for single reads, 
``interleaved:`` for paired reads in a single file,
or one or more ``pairs:`` sections giving read pairs as a pair of files.

Analysis in one command, with more options
------------------------------------------

All of these options are available in the one-command version as well.
The options for each command are given in subsections
``clip:``, ``align:``, ``filter:``, and ``reconsensus:``.
``align: shrimp`` aligns using SHRiMP (the default) and 
``align: bowtie`` aligns using Bowtie2.

If you run the pipeline again with some options changed,
it will run from the first stage where an option was changed.
Also note Nesoni is not smart enough to noticed changed input files,
just changed options.
Use ``--make-do all`` to force the entire pipeline to be re-run.

From python
-----------

We can create an object representing the action of analysing the sample,
and then tell it to "make".

.. code-block:: python
  
  import nesoni
  
  action = nesoni.Analyse_sample(
      'mysample',
      reference = 'myref',
      reads = [ 'myreads.fastq' ],
      )

  def main():  
      action.make()

  if __name__ == '__main__': 
      nesoni.run_script(main)

Wrapping ``main`` in ``nesoni.run_script`` gives access to the make options,
such as ``--make-do all`` to force the entire pipeline to re-run.


