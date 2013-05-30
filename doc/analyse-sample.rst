
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

  nesoni make-reference: myref ref-sequences.fa \
      --ls yes --genome yes

This says to generate a SHRiMP base-space index,
so read alignment will start up faster,
and a .genome file for the IGV browser.

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

::

  nesoni clip:
  nesoni shrimp:
  nesoni filter:
  nesoni reconsensus:


Analysis in one command, with more options
------------------------------------------

All of these options are available in the one-command version as well.

From python
-----------

We create an object representing the action of analysing the sample,
and then tell it to "make".
::
  
  import nesoni
  
  action = nesoni.Analyse_sample(
    )

  def main():  
      action.make()

  if __name__ == '__main__': 
      nesoni.run_script(main)
  
