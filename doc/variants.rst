
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

Variant calling in a single command
-----------------------------------

Variant calling as individual steps
-----------------------------------

The first step is to use freebayes to identify variants.
By default this casts a wide net, 
because the next step will be to filter what it finds using our own statistical method.


Variant calling in a single command, with more options
------------------------------------------------------

Slicing and dicing variants
---------------------------