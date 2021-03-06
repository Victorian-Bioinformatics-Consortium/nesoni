Nesoni 
======

With the closure of the Victorian Bioinformatics Consortium, I anticipate little further development on Nesoni.

Nesoni still serves as the basis for "Tail Tools", which remains under development. The parallel processing code is also used by "Demakein".

- Paul Harrison, January 2015


What Nesoni does
----------------

Nesoni is a high-throughput sequencing data analysis toolset, which the Victorian Bioinformatics Consortium developed to cope with the flood of Illumina, 454, and SOLiD data being produced. 

The VBC's work was largely with bacterial genomes, and the design tradeoffs in Nesoni reflect this.

Nesoni focusses on analysing the alignment of reads to a reference genome. Use of the SHRiMP and Bowtie2 read aligners is automated by nesoni. We use SHRiMP as it is able to detect small insertions and deletions in addition to SNPs. Output from other aligners may be imported in SAM format.

Nesoni can call a consensus of read alignments, taking care to indicate ambiguity. This can then be used in various ways: to determine the protein level changes resulting from SNPs and indels, to find differences between multiple strains, or to produce n-way comparison data suitable for phylogenetic analysis in SplitsTree4.


Requirements
============

Python 2.7 or higher. Use of PyPy where possible is highly recommended 
for performance.

Python libraries
* Recommended:
  * BioPython [1]
* Optional (used by non-core nesoni tools):
  * numpy
  * matplotlib

External programs:
* Required:
  * SHRiMP or Bowtie2
  * samtools (0.1.19 or higher)
* Required for VCF based variant calling:
  * Picard [2]
  * Freebayes (9.9.2 or higher)
* Optional for VCF based variant calling:
  * SplitsTree4

R libraries required by R-based tools (mostly for RNA-seq):
* Required:
  * varistran from https://github.com/MonashBioinformaticsPlatform/varistran
  * limma, edgeR (3.2.4 or higher) from BioConductor
  * seriation from CRAN
* Optional:
  * goseq from BioConductor


[1] BioPython is used for reading GenBank files.

[2] There does not seem to be a standard way to install .jar files. 
    Nesoni will search for .jar files in directories listed in 
    environment variables $PATH and $JARPATH.


Installation
============

The easy way to install or upgrade:

    pip install -I nesoni

Then type "nesoni" and follow the command to install the R module.

See below for more ways to install nesoni.


Advanced Installation 
---------------------

From source, download and untar the source tarball, then:

    python setup.py install

Optional:

    R CMD INSTALL nesoni/nesoni-r


For PyPy it seems to be currently easiest to install nesoni in 
a virtualenv:

    virtualenv -p pypy my-pypy-env
    my-pypy-env/bin/pip install -I git+https://bitbucket.org/pypy/numpy.git
    my-pypy-env/bin/pip install -I biopython 
    my-pypy-env/bin/pip install -I nesoni

You can also set up a CPython virtualenv like this:

    virtualenv my-python-env
    my-python-env/bin/pip install -I numpy 
    my-python-env/bin/pip install -I matplotlib 
    my-python-env/bin/pip install -I biopython 
    my-python-env/bin/pip install -I nesoni


Using nesoni
============

Type

    nesoni

for usage instructions.

nesoni can also be used without installing, from the directory in
which you unpacked it:

    python -m nesoni




