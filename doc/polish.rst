
Polish a de-novo assembly
=========================

Suppose you have assembled a genome, 
but suspect it has some minor errors.
For example, 
it might have been assembled from 454 or IonTorrent reads
and have homopolymer length errors.
The assembly can then be polished using, say, an Illumina read set
(higher quality but shorter).

*Method 1:* Follow the instructions in :doc:`analyse-sample`. 
A file ``consensus_masked.fa`` is produced
with SNPs and indels patched into the reference sequences.
In this file, 
parts of the genome with insufficient coverage
to be sure of the base calling are given in lower case.

You may wish to reduce the parameter ``reconsensus: --cutoff``,
in order to favour the Illumina read set even if the evidence from it is fairly weak.

*Method 2:* Follow the instructions in :doc:`variants`.
A file with the variants detected patched into the reference sequences is produced
as part of the output.

Again, you may wish to reduce ``vcf-filter: --min-gq``
to favour using whatever evidence there is in the Illumina reads.