#!/usr/bin/env python

from distutils.core import setup

import nesoni

setup(name='nesoni',
      version=nesoni.VERSION,
      
      packages = [
          'nesoni', 
          'treemaker'
      ],
      
      package_data = {
          'treemaker' : ['*.pyx'],
          'nesoni' : ['*.pyx','nesoni-r/DESCRIPTION','nesoni-r/R/*'],                    
      },
      
      scripts=['nesoni_scripts/nesoni'],

      classifiers = [
          'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      ],
      url = 'http://bioinformatics.net.au/software.nesoni.shtml',
      author = 'Paul Harrison',
      author_email = 'paul.harrison@monash.edu',
)
