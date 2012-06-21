#!/opt/python/bin/python2.6


import os, datetime

import nesoni
from nesoni import config

README = open('README','rb').read().split('~')
REQUIREMENTS = README[1]
INSTALL = README[2]

HELP = config.strip_color(nesoni.USAGE).replace('<','&lt;').replace('>','&gt;')

PAGE = r"""<!--#include virtual="top.html" -->

<style>
pre { line-height: 100%%; font-size: 130%%; }
</style>

<h2>Nesoni</h2>

<p>
Nesoni is a high-throughput sequencing data analysis toolset,
which the VBC has developed to cope with the flood of Illumina, 
454, and SOLiD data now being produced. 

<p>
Our work is largely with bacterial genomes, and the design tradeoffs 
in nesoni reflect this.

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

<p>
Alternatively, the raw counts of bases at each position in the reference
seen in two different sequenced strains can compared using Fisher's Exact Test.

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

<h3>k-mer and De Bruijn graph tools</h3>

<p>
Nesoni also includes some highly experimental tools for working with sets of k-mers and De Bruijn graphs. You can:

<ul>
<li>Produce a 2D layout of a De Bruijn graph, and interact with it: zoom in to examine details, overlay read-pair data, overlay sequences on top of it (for example, to examine the behaviour of a de-novo assembler such as Velvet).
<li>Clip a set of reads to remove low-frequency or high-frequency k-mers, or k-mers where there is a more frequent k-mer differing by one SNP.
</ul>

<h3>Documentation</h3>

<table cellpadding="0" cellspacing="0"><tr>
<td valign="top" style="padding-right: 1em"><a href="nesoni_ba2009_poster.pdf"><img src="nesoni_poster.png"></a></td>
<td valign="top">
<p>
<a href="nesoni_ba2009_poster.pdf">This poster</a>, presented at BA2009, gives an overview of Nesoni's capabilities.
</td></tr></table>

<p>
Nesoni provides the following specific usage information when run with no parameters:

<pre>
%(HELP)s
</pre>


<h3>Download</h3>

<p>%(date)s:

<ul>
<li> <a href="%(release_tarball_name)s">%(release_tarball_name)s</a> </li>
</ul>

<p>Nesoni is free software, released under the GPL (version 2).

<pre>
%(REQUIREMENTS)s
%(INSTALL)s
</pre>

<h3>Contact</h3>
<ul>
<li><a href='mailto:paul.harrison@monash.edu'>Paul Harrison</a>  
<li><a href='mailto:torsten.seemann@monash.edu'>Torsten Seemann</a>
</ul>


<!--#include virtual="bot.html" -->
"""

os.environ['PATH'] = '/bio/sw/python/bin:' + os.environ['PATH']

# RAGE
os.system('rm MANIFEST')

release_tarball_name = 'nesoni-%s.tar.gz' % nesoni.VERSION
assert not os.path.exists('dist/'+release_tarball_name), release_tarball_name + ' already exists'

try:
    assert 0 == os.system('cd test && pypy test_nesoni.py')
    
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_scripts --install-dir /bio/sw/python/bin/')
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_lib')
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/python2.6 setup.py install_lib')
    assert 0 == os.system('sudo pypy setup.py install --home /bio/sw/python')
    assert 0 == os.system('sudo R CMD INSTALL --library=/bio/sw/R nesoni-r')
    
    print
    print
    
    date = datetime.date.today().strftime('%e %B %Y')
    
    assert 0 == os.system('python setup.py sdist')
    
    f = open('/home/torsten/public_html/vicbioinformatics.com/software.nesoni.shtml','wb')
    f.write(PAGE % locals())
    f.close()
    
    assert 0 == os.system('cp dist/%s /home/torsten/public_html/vicbioinformatics.com' % release_tarball_name)
    assert 0 == os.system('cd /home/torsten/public_html/vicbioinformatics.com ; make install')
    
    #assert os.path.exists(release_tarball_name), release_tarball_name + ' failed to build'
    
    #os.system('rsync -rltvz /opt/python/ msgln1.its.monash.edu.au:opt-python/')

except:
    os.system('rm dist/%s' % release_tarball_name)
    raise    

