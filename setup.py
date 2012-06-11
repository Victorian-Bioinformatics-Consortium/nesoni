#!/opt/python/bin/python2.6

from distutils.core import setup

import nesoni

import sys, os
if not os.path.isfile(sys.executable):
    print
    print 'Python thinks this is its executable:', sys.executable
    print
    print 'I disagree.'
    print
    print 'If using sudo, try sudo -E'
    print
    sys.exit(1)

setup(name='nesoni',
      version=nesoni.VERSION,
      
      packages = [
          'nesoni', 
          'treemaker'
      ],
      
      package_data = {
          'treemaker' : ['*.pyx'],
          'nesoni' : ['*.pyx'],
      },
      
      scripts=['nesoni_scripts/nesoni'],

      classifiers = [
          'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      ],
      url = 'http://bioinformatics.net.au/software.nesoni.shtml',
      author = 'Paul Harrison',
      author_email = 'paul.harrison@monash.edu',
)
