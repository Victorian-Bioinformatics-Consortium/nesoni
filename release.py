#!/usr/bin/env python


import os, datetime, sys

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

<h3>Download</h3>

<p>%(date)s:

<ul>
<li> <a href="%(release_tarball_name)s">%(release_tarball_name)s</a> </li>
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

<p style="font-size: 125%%">
<a href="https://docs.google.com/document/pub?id=1vt__lbYwnoMGU0m1XTqw4j-H0URez201i9JLLtN_aVA">Further documentation</a>
</p>

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

<pre>
%(REQUIREMENTS)s
%(INSTALL)s
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

release_tarball_name = 'nesoni-%s.tar.gz' % nesoni.VERSION
assert 'force' in sys.argv[1:] or not os.path.exists('dist/'+release_tarball_name), release_tarball_name + ' already exists'

def sh(cmd): assert 0 == os.system(cmd)

try:
    #assert 0 == os.system('cd test && pypy test_nesoni.py')
    
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_scripts --install-dir /bio/sw/python/bin/')
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_lib')
    #assert 0 == os.system('sudo -E /bio/sw/python/bin/python2.6 setup.py install_lib')
    
    #assert 0 == os.system('sudo pypy setup.py install --home /bio/sw/python')
    #assert 0 == os.system('sudo PYTHONPATH=/bio/sw/python/lib/python '
    #                      'python setup.py install --home /bio/sw/python')
    
    os.system('sudo rm -r /bio/sw/python/env-pypy')
    sh('sudo virtualenv -p pypy /bio/sw/python/env-pypy')
    sh('sudo ln -s pypy /bio/sw/python/env-pypy/bin/pypy-bio')
    sh('sudo /bio/sw/python/env-pypy/bin/pip install biopython')
    sh('sudo /bio/sw/python/env-pypy/bin/python setup.py install')

    os.system('sudo rm -r /bio/sw/python/env-python')
    sh('sudo virtualenv -p python /bio/sw/python/env-python')
    sh('sudo ln -s python /bio/sw/python/env-python/bin/python-bio')
    sh('sudo /bio/sw/python/env-python/bin/pip install numpy biopython')
    sh('sudo /bio/sw/python/env-python/bin/python setup.py install')        
    
    sh('sudo R CMD INSTALL --library=/bio/sw/R nesoni/nesoni-r')
    
    os.system('sudo rm -r nesoni.egg-info')
    
    print
    print
    
    date = datetime.date.today().strftime('%e %B %Y')
    
    assert 0 == os.system('python setup.py sdist')
    
    assert 0 == os.system('cp dist/%s /home/websites/vicbioinformatics.com' % release_tarball_name)
    
    f = open('/home/websites/vicbioinformatics.com/software.nesoni.shtml','wb')
    f.write(PAGE % locals())
    f.close()
    
    #assert os.path.exists(release_tarball_name), release_tarball_name + ' failed to build'
    
    #os.system('rsync -rltvz /opt/python/ msgln1.its.monash.edu.au:opt-python/')

except:
    os.system('rm dist/%s' % release_tarball_name)
    raise    

os.system('python setup.py sdist upload')        


