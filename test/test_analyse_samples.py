
import sys, os
sys.path.insert(0, os.path.join(os.path.split(__file__)[0], '../'))

import os
from Bio import Entrez

import nesoni

from nesoni import config, io

# ref.: http://wilke.openwetware.org/Parsing_Genbank_files_with_Biopython.html

# replace with your real email (optional):
#Entrez.email = 'whatever@mail.com'

@config.Positional('accession')
@config.Positional('email',affects_output=False)
class Fetch_genbank(config.Action_with_output_dir):
    accession = None
    email = None
    
    def state_filename(self):
        return os.path.join(self.output_dir, self.accession + '.state')
    
    def run(self):
        work = self.get_workspace()
        acc = self.accession
        
        Entrez.email = self.email

        #handle=Entrez.efetch(db='nucleotide',id=acc,rettype='gb')
        handle=Entrez.efetch(db='nuccore',id=acc,rettype='gbwithparts')
        with open(work/(acc+'.gbk'),'wb') as f:
            f.write(handle.read())
        handle.close()

@config.Positional('run_accession')
class Fetch_sra(config.Action_with_output_dir):
    run_accession = None
    
    def state_filename(self):
        return os.path.join(self.output_dir, self.run_accession + '.state')
    
    def run(self):
        work = self.get_workspace()
        acc = self.run_accession
        io.execute(
            'wget -c URL',
            #URL='http://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra'
            URL='http://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra'
                % (acc[:3],acc[:6],acc,acc),
            cwd=work.working_dir,
            )
        
        io.execute(
            'fastq-dump --split-files --bzip2 FILENAME',
            FILENAME='./'+acc+'.sra',
            cwd=work.working_dir,
            )



@config.Positional('email', 'Email address to give when accessing NCBI', affects_output=False)
class Test_analyse_samples(config.Action_with_output_dir):
    email = None
    genbanks = 'NC_017331.1'.split()
    accessions = 'ERR048386 ERR048396 ERR048414 ERR086164 ERR086220 ERR086222'.split()

    def run(self):
        work = self.get_workspace()
        with nesoni.Stage() as stage:
            for accession in self.genbanks:
                Fetch_genbank(
                    work/'genbank',
                    accession,
                    self.email,
                    ).process_make(stage)
        
            for accession in self.accessions:
                Fetch_sra(
                    work/'sra', 
                    accession
                    ).process_make(stage)
                    
        nesoni.Make_reference(
            output_dir=work/'TW20', 
            filenames=[ work/('genbank',accession+'.gbk') for accession in self.genbanks ],
            genome=True, bowtie=True, ls=True, snpeff=True,
            ).make()
            
        analyser = nesoni.Analyse_samples
        analyser(
            work/'analysis',
            work/'TW20',
            samples=[
                nesoni.Analyse_sample(
                    accession,
                    pairs=[[ work/('sra',accession+'_1.fastq.bz2'),work/('sra',accession+'_2.fastq.bz2') ]]
                    )
                for accession in self.accessions
                ]
            ).make()


if __name__ == '__main__':
    nesoni.run_tool(Test_analyse_samples(output_dir='output'))        
    
