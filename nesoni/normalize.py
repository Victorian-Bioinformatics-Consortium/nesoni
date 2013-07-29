
from . import grace, config, working_directory, io

import os, itertools, collections
    
@config.help(
    'Create a normalization file based on depth of coverage.',
    'Normalization is based on total depth of coverage of those locations in the genome where read mapping is entirely unambiguous.'
    )
@config.Main_section('working_dirs', 'Working directories containing the results of "filter:" or "consensus:".')
class Norm_from_samples(config.Action_with_prefix):
    working_dirs = [ ]
    
    def run(self):
        sample_names = [ os.path.split(dirname)[1] for dirname in self.working_dirs ]
        workspaces = [ working_directory.Working(dirname, must_exist=True) for dirname in self.working_dirs ]        
        depths = [ item.get_object('depths.pickle.gz') for item in workspaces ]

        lengths = workspaces[0].get_reference().get_lengths()
        chromosome_names = [ name for name, length in lengths ]
        lengths = dict(lengths)
        
        totals = [ 0 ] * len(sample_names)
        sites = 0
        good_sites = 0
        for name in chromosome_names:
            grace.status('Normalization '+name)
            iter_unambiguous0 = itertools.izip(*[ item[name].depths[0] for item in depths ])   
            iter_unambiguous1 = itertools.izip(*[ item[name].depths[1] for item in depths ])   
            iter_ambiguous0 = itertools.izip(*[ item[name].ambiguous_depths[0] for item in depths ])   
            iter_ambiguous1 = itertools.izip(*[ item[name].ambiguous_depths[1] for item in depths ])
            for unambiguous0, unambiguous1, ambiguous0, ambiguous1 in itertools.izip(
                iter_unambiguous0, iter_unambiguous1, iter_ambiguous0, iter_ambiguous1):
                sites += 1
                if unambiguous0 != ambiguous0 or unambiguous1 != ambiguous1:
                    continue
                good_sites += 1
                for i in xrange(len(sample_names)):
                    totals[i] += ambiguous0[i]+ambiguous1[i]
        grace.status('')
        
        self.log.log('Locations used: %s of %s\n' % (grace.pretty_number(good_sites),grace.pretty_number(sites)))
        
        for name, value in zip(sample_names, totals):
            self.log.datum(name, 'mean depth', float(value) / sites)
        
        mult = 1
        for item in totals: mult *= item
        geomean = mult ** (1.0/len(sample_names))
        
        self.log.log('Geometric mean: %.3f\n' % (geomean/sites))

        def writer():
            for i, name in enumerate(sample_names):
                row = collections.OrderedDict()
                row['Name'] = name
                row['Normalizing.multiplier'] = geomean / totals[i]
                yield row
        io.write_csv(self.prefix+'.csv', writer(), comments=['Normalization'])
        



