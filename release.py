#!/usr/bin/env python


import os, datetime, sys

import nesoni
from nesoni import config

def cleanup(text):
    lines = text.split('\n')    
    lines[0] = '<b style="font-size: 125%">'+lines[0].strip()+'</b>'
    del lines[1]
    return '\n'.join(lines)

VERSION = nesoni.VERSION

README = open('README','rb').read().split('~')
REQUIREMENTS = cleanup(README[1])
INSTALL = cleanup(README[2])

HELP = config.strip_color(nesoni.USAGE).replace('<','&lt;').replace('>','&gt;')
HELP = HELP[:HELP.index('If a pipeline tool is run again')]

PAGE = r"""<!--#include virtual="top.html" -->

<style>
pre { line-height: 100%%; font-size: 100%%; }
.box { background: #eeffee; border: 1px solid #008800; margin: 1em; padding: 1em; float: left; width: 25em; }
h1, h2, h3 { clear: both; }
</style>

<div style="float: left;">

<h2>Nesoni</h2>

<h3>Download</h3>

<div style="float: right; margin: 1em; margin-right: 20%%; border: 2px solid black; padding: 1em; font-size: 150%%;">
    <a href="nesoni-cookbook/">Read the Nesoni Cookbook</a>
</div>

<p>%(date)s:

<ul>
<li> <a href="%(release_tarball_name)s">%(release_tarball_name)s</a> </li>

<li> <a href="https://pypi.python.org/pypi/nesoni/">PyPI page</a>

<li> <a href="https://github.com/Victorian-Bioinformatics-Consortium/nesoni">Development version on github</a>

</ul>

<p>Nesoni is free software, released under the GPL (version 2).

<h3>Description</h3>

<p>
Nesoni is a high-throughput sequencing data analysis toolset,
which the VBC has developed to cope with the flood of Illumina, 
454, and SOLiD data now being produced. 

<p>
Our work is largely with bacterial genomes, and the design tradeoffs 
in nesoni reflect this.

<div class="box">
  <h3>Alignment to reference</h3>
  
  <p>
  Nesoni focusses on analysing the alignment of reads to a reference 
  genome. We use the SHRiMP read aligner, as it is able to detect 
  small insertions and deletions in addition to SNPs.
  
  <p>
  Nesoni can call a consensus of read alignments, taking care to
  indicate ambiguity. This can then be used in various ways: to determine
  the protein level changes resulting from SNPs and indels, to find 
  differences between multiple strains, or to produce n-way comparison data 
  suitable for phylogenetic analysis in SplitsTree4. 
  
  </div>

<div class="box">
  <h3>Explore your data in IGV</h3>
  
  <p>
  Besides producing BAM files,
  Nesoni can produce multi-sample depth of coverage plots for
  <a href="http://www.broadinstitute.org/igv/">the IGV browser</a>.
  Display and quickly navigate depth of coverage plots from tens of samples
  simultaneously, from whole chromosome overview to the individual base level.
  
  <p>
  Tracks indicating the presence of read multi-mapping, 
  locations of 5' and 3' read ends, and total depth are also produced.
  
  </div>

<div class="box">
  <h3>Multiple perspectives on expression data</h3>
  
  <p>
  Nesoni provides multiple tools for exploring expression data.
  Don't just blindly reach into your data and pick out a handful of genes,
    actually understand the story your data is telling you.
  
  <ul>
  <li>Differential expression testing with limma and edgeR.
      Accompanying heatmaps make it easy to understand the tradeoffs of within group and between group variation
      that drive these tools.
      </li>
  <li>Heatmaps with hierarchical clustering.
      </li>
  <li>Non-negative Matrix Factorization, a fuzzy clustering technique.
      </li>
  </ul>
  
  </div>

<div class="box">
  <h3>Workflow engine</h3>
  
  <p>
  Nesoni includes tools that make it easy to write parallel processing pipelines
  in Python.
  
  <p>
  Pipelines are expressed as Python functions. 
  The translation of a serial program with for-loops and function calls
  into a parallel program requires only simple localized modifications to the code. 
  
  <p>
  Pipelines expressed in this way are <i>composable</i>, just like ordinary functions.
  
  <p>
  Much like <tt>make</tt>, the resultant program will only re-run tools as necessary 
  if parameters are changed.
  The dependancy structure is implicit from the parallel program,
  if a tool needs to be re-run, 
  only things that <i>must</i> execute after that tool also need re-running.
  
  <p>
  <a href="https://docs.google.com/document/pub?id=1vt__lbYwnoMGU0m1XTqw4j-H0URez201i9JLLtN_aVA">Further documentation</a>
  </p>
  
  </div>

<!--
<div class="box">
  <h3>k-mer and De Bruijn graph tools</h3>
  
  <p>
  Nesoni also includes some highly experimental tools for working with sets of k-mers and De Bruijn graphs. You can:
  
  <ul>
  <li>Produce a 2D layout of a De Bruijn graph, and interact with it: zoom in to examine details, overlay read-pair data, overlay sequences on top of it (for example, to examine the behaviour of a de-novo assembler such as Velvet).
  <li>Clip a set of reads to remove low-frequency or high-frequency k-mers, or k-mers where there is a more frequent k-mer differing by one SNP.
  </ul>
  
  </div>
-->

<div style="clear: both"/>
<br/><br/>

<div style="float: right; margin: 2em; width: 400px;">
<a href="nesoni_ba2009_poster.pdf"><img src="nesoni_poster.png"></a>
<p><a href="nesoni_ba2009_poster.pdf">This poster</a>, presented at BA2009, describes some of Nesoni's basic capabilities.
</div>

<pre>
%(REQUIREMENTS)s
%(INSTALL)s
</pre>

<h3>Documentation</h3>

<ul>
<li><a href="nesoni-cookbook/">Cookbook</a>
</ul>

<p>
Nesoni provides the following specific usage information when run with no parameters:

<pre>
%(HELP)s
</pre>

<h3>Contact</h3>
<ul>
<li><a href='mailto:paul.harrison@monash.edu'>Paul Harrison</a>  
<li><a href='mailto:torsten.seemann@monash.edu'>Torsten Seemann</a>
</ul>

<h3>Older versions</h3>

<ul>
%(OLD)s
</ul>

<p>

</div>

<!--#include virtual="bot.html" -->
"""





FITNOISE_PAGE = r"""<!--#include virtual="top.html" -->

<style>
pre { padding: 1em; margin: 1em; border: 1px solid #888; width: 60em; }
</style>

<h2>Fitnoise</h2>

<div style="float:right; margin-right:2em;">
<a href="documents/Fitnoise-poster-ABiC-October-2014.pdf">
<div style="border: 1px solid black; padding: 1em">
<img width="150" 
     src="documents/Fitnoise-poster-ABiC-October-2014.png">
</div>
Poster describing Fitnoise 
</a>
<br/>presented at <a href="http://bioinformatics.net.au/abic2014/">ABiC 2014</a>
</div>

Fitnoise is R+ software for Empirical Bayesian linear modelling.
It has capabilities very similar to 
<a href="http://bioinf.wehi.edu.au/limma/">Limma</a>
but allows parametric noise models with an arbitrary number of parameters.
Parameters are fitted using <a href="https://en.wikipedia.org/wiki/Restricted_maximum_likelihood">REML</a>.

<p>%(date)s, from Nesoni version %(VERSION)s:

<ul>
<li><a href="fitnoise.R">fitnoise.R</a>
</ul>

<p>
fitnoise.R is distributed as part of 
<a href="software.nesoni.shtml">Nesoni</a>, 
but also works as a standalone R+ module.

<p>
A unique feature of Fitnoise is its ability to
detect differential poly(A) tail length in PAT-Seq data. 
See also <a href="software.tail-tools.shtml">Tail Tools</a>.

<h3>Basic usage</h3>

<p>
Say you have a dgelist, mydgelist, and a design matrix mydesign, and a set of coefficients, testcoefs, to test:

<p>
H1 is all columns in the design matrix. H0 is those columns not named in testcoefs.

<p>Run R+ with:

<pre>
OPENBLAS_NUM_THREADS=1 R
</pre>

<p>
...

<pre>
source("fitnoise.R")

# voom from limma to calculate weights
myelist <- voom(mydgelist, mydesign)

myfit <- fit.elist(myelist, mydesign, model=model.t.standard, cores=8)

#examine noise fit
myfit

#produce a toptable-like result data-frame
result <- test.fit(myfit, coefs=testcoefs)
</pre>

<p>
<b>When examining myfit, if the "noise combined p-value" is low this indicates a poor
fit by the noise model.</b> A more accurate noise model may be required.

<h3>Noise models</h3>

<p>
The novel feature of Fitnoise is the availability of different noise models. 
Above we used model.t.standard, which performs moderated t tests or F tests similarly to limma.
Available are:
<ul>

<li>model.t.standard, attempts to closely emulate limma.
The weights matrix is used if present in the EList.

<li>model.t.independent, which simply performs independent t tests or F tests.

<li>model.normal.standard, which performs z tests or chi-square tests.
The weights matrix is used if present in the EList.

<li>model.t.patseq and model.normal.patseq, for PAT-Seq poly(A) tail length data.

<li><b>Experiemental:</b> model.t.per.sample.var and model.normal.per.sample.var
fit per-sample variances, similar to arrayWeights in limma.
The weights matrix is used if present in the EList.
Would suggest only using this with an absolute minimum of three samples per group, 
or with negative controls (see below).

</ul>

<p>
Writing new noise models is straightforward.

<p>
Fitnoise handles missing data gracefully.

<h3>Use without replicates</h3>

Say mydesign represents a null hypothesis H0.
We seek genes that are outliers given a noise model.
t-distribution based models allow for outliers,
so we instead use a normal-distribution based model.
This is effectively asserting that genes have fixed biological variability.
This is a poor assertion to make.

<pre>
myfit <- fit.elist(myelist, mydesign, model=model.normal.standard, cores=16)

#produce a toptable-like result data-frame based on surprising deviation 
#from the noise model in individual genes
result <- test.fit(myfit)
</pre>

<p>
Note how we never specified an alternative hypothesis H1.

<p>
Equivalently, you can call fit.elist with H1 as the design and an extra parameter "noise.design"
giving H0, then call test.fit normally. This will give identical p-values,
and the output will include coefficient estimates.

<p>
This produces a ranking of genes by interest, taking into account uncertainty
in their expression levels. <b>p-values / FDR-values are for ranking purposes only.</b>
This method can not distinguish highly variable genes from truly differentially expressed genes.
Now go and tell your collaborator to produce replicates next time.


<h3>Experimental: negative control features</h3>

<p>
A negative control feature is one in which the noise distribution
is typical, and which is believed to not have been changed by the experiment.
Constitutively expressed housekeeping genes would make good negative controls.
Spiked in RNA would <i>not</i> make good negative controls as they
will not have typical biological noise.
Empirically chosen negative controls 
might underestimate or misestimate the noise distribution, use with caution
(eg features chosen from non-significant features from a differential test without negative controls).

<p>
A logical vector of negative control features can be passed to fit.elist
as parameter "controls". 
You can also pass your null hypothesis design as "control.design"
(this defaults to a single column of ones).

<p>
Negative controls, if known, allow noise to be fitted more accurately.
For example they may be useful for model.t.per.sample.var if you have few samples.
If the noise model is covariant between samples
with covariance spanning experimental groups,
negative controls are absolutely necessary
(none of the current noise models use this.)

<br/><br/><br/><br/><br/>

<!--#include virtual="bot.html" -->
"""








old = [ ]
for item in os.listdir('/home/websites/vicbioinformatics.com'):
    if item.startswith('nesoni-') and item.endswith('.tar.gz'):
       old.append(item)
       
old.sort(key=lambda item:
    [ int(item2) if item2.isdigit() else item2 for item2 in item.replace('-','.').split('.') ]
)

OLD = '\n'.join([
    '<li><a href="%s">%s</a>' % (item,item) for item in old[::-1]
])


os.environ['PATH'] = '/bio/sw/python/bin:' + os.environ['PATH']

# Force MANIFEST rebuild
os.system('rm MANIFEST')

do_install = 'noinstall' not in sys.argv[1:]
do_upload = 'noupload' not in sys.argv[1:]
do_rebuild = 'rebuild' in sys.argv[1:]

release_tarball_name = 'nesoni-%s.tar.gz' % nesoni.VERSION
if do_upload:
    assert 'force' in sys.argv[1:] or not os.path.exists('dist/'+release_tarball_name), release_tarball_name + ' already exists'

def sh(cmd): assert 0 == os.system(cmd)

date = datetime.date.today().strftime('%e %B %Y')

try:
    #assert 0 == os.system('cd test && pypy test_nesoni.py')
    
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_scripts --install-dir /bio/sw/python/bin/')
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_lib')
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/python2.6 setup.py install_lib')
    
    #assert 0 == os.system('sudo pypy setup.py install --home /bio/sw/python')
    #assert 0 == os.system('sudo PYTHONPATH=/bio/sw/python/lib/python '
    #                      'python setup.py install --home /bio/sw/python')
    
    if do_install:
        if do_rebuild:
            os.system('sudo -H rm -r /bio/sw/python/env-pypy')
            sh('sudo -H virtualenv -p /bio/sw/python/download/pypy/bin/pypy /bio/sw/python/env-pypy')
            sh('sudo -H ln -s pypy /bio/sw/python/env-pypy/bin/pypy-bio')
            sh('sudo -H /bio/sw/python/env-pypy/bin/pip install --upgrade distribute')
            sh('sudo -H /bio/sw/python/env-pypy/bin/pip install --upgrade biopython')
        sh('sudo -H /bio/sw/python/env-pypy/bin/python setup.py install')
        sh('sudo -H /bio/sw/python/env-pypy/bin/pip install --upgrade tail-tools')

        if do_rebuild:        
            os.system('sudo -H rm -r /bio/sw/python/env-python')
            sh('sudo -H virtualenv -p python /bio/sw/python/env-python')
            sh('sudo -H ln -s python /bio/sw/python/env-python/bin/python-bio')
            sh('sudo -H /bio/sw/python/env-python/bin/pip install --upgrade distribute')
            sh('sudo -H /bio/sw/python/env-python/bin/pip install --upgrade numpy')
            sh('sudo -H /bio/sw/python/env-python/bin/pip install --upgrade matplotlib')
            sh('sudo -H /bio/sw/python/env-python/bin/pip install --upgrade biopython')
        sh('sudo -H /bio/sw/python/env-python/bin/python setup.py install')        
        sh('sudo -H /bio/sw/python/env-python/bin/pip install --upgrade tail-tools')
        
        sh('sudo -H R CMD INSTALL --library=/bio/sw/R nesoni/nesoni-r')
        
        os.system('sudo -H rm -r nesoni.egg-info build')
        
        print
        print
        
        sh('python setup.py sdist')
        
        sh('cp dist/%s /home/websites/vicbioinformatics.com' % release_tarball_name)
        f = open('/home/websites/vicbioinformatics.com/software.nesoni.shtml','wb')
        f.write(PAGE % locals())
        f.close()

    sh('cd doc && make html')
    sh('cp -r doc/_build/html/* /home/websites/vicbioinformatics.com/nesoni-cookbook/')    

    sh('cp nesoni/nesoni-r/R/fitnoise.R /home/websites/vicbioinformatics.com/')
    f = open('/home/websites/vicbioinformatics.com/software.fitnoise.shtml','wb')
    f.write(FITNOISE_PAGE % locals())
    f.close()
    
    #assert os.path.exists(release_tarball_name), release_tarball_name + ' failed to build'
    
    #os.system('rsync -rltvz /opt/python/ msgln1.its.monash.edu.au:opt-python/')

except:
    os.system('rm dist/%s' % release_tarball_name)
    raise    

if do_upload:
    os.system('python setup.py sdist upload')        


