from itertools import izip_longest
import logging
import multiprocessing as mp
import os
import subprocess as sp
import sys
import tempfile as tf

import ROOT

from sample import Sample
import settings


class Process(object):

    def __init__(self, name = '', process_cut = '', samples = [], sample_cuts = [], **kwargs):
        
        self.logger = logging.getLogger('Process')
        self.logger.info('Initialized for {}'.format(name))

        self.name = name
        self.process_cut = process_cut
        self.samples = samples
        self.sample_cuts = sample_cuts

    def make(self):
    
        # Output Directory
        try:
            os.makedirs(settings.PROCESS_DIR)
        except OSError:
            if not os.path.isdir(settings.PROCESS_DIR):
                raise

        # Prepare Samples
        sample_files = self._check_samples()

        inputfiles = []

        n_tasks = 0
        tasks = mp.Queue()
        results = mp.Queue()

        for sample_file, sample_cut in izip_longest(sample_files, self.sample_cuts, fillvalue = ''):
            if (sample_cut is '') and (self.process_cut is ''):
                inputfiles.append(sample_file)
            else:
                tasks.put((sample_file, sample_cut))
                n_tasks += 1

        # Parallel Cut 
        if (n_tasks > 0):
            
            tmpdir = tf.mkdtemp(prefix = self.name, dir = settings.PROCESS_DIR)
            self.tmpdir = tmpdir + '/'

            _processes = [
                mp.Process(target = self._cut_sample, args = (tasks, results))
                for i in range(min(n_tasks, mp.cpu_count()))
            ]

            for p in _processes:
                tasks.put(None)
                p.start()

            for p in _processes:
                p.join()

            results.put(None)

            for r in iter(results.get, None):
                inputfiles.append(r)

        # hadd Files
        outputfile = settings.PROCESS_DIR + self.name + '.root'
        
        if (len(inputfiles) == 1) and (settings.SAMPLE_DIR in inputfiles[0]):
            # Make symbolic link instead to save memory.
            sp.check_call(['ln', '-s', inputfiles[0], outputfile])
        else:
            sp.check_call(['hadd', '-f', outputfile] + inputfiles)

        sp.check_call(['rm', '-r', self.tmpdir])
         
    def _check_samples(self):
    
        sample_files = []
     
        for sample in self.samples:

            fname = settings.SAMPLE_DIR + sample + '.root'

            if os.path.isfile(fname):
                sample_files.append(fname)
            else:
                self.logger.info('Getting missing sample {}'.format(sample))
                Sample(sample, **settings.SAMPLES[sample]).make()
                sample_files.append(fname)

        return sample_files

    def _cut_sample(self, tasks = None, results = None):
    
        for sample_file, sample_cut in iter(tasks.get, None):

            outname = self.tmpdir + os.path.basename(sample_file)

            infile = ROOT.TFile(sample_file, 'read')
            outfile = ROOT.TFile(outname, 'recreate')

            intree = infile.Get('tree')
            outtree = intree.CopyTree(self.process_cut & sample_cut)

            outtree.Write()
            
            for key in infile.GetListOfKeys():
                if key.GetName() == 'tree':
                    continue
                obj = key.ReadObj()
                obj.Write()

            outfile.Close()
            infile.Close()

            results.put(outname)

#------
# Main 
#------

if __name__ == '__main__':

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    for name in sys.argv[1:]:

        logging.basicConfig(level = logging.INFO,
                            format = '%(name)s(%(levelname)s) - %(message)s')

        Process(name, **settings.PROCESSES[name]).make()

