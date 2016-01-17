import glob
import logging
import multiprocessing as mp
import os
import subprocess as sp
import sys
import tempfile as tf

import ROOT

from process import Process, PROCESS_DIR
from settings import WORK_DIR, PROCESSES, REGIONS


# Output Directory
REGION_DIR = WORK_DIR + 'regions/'

class Region(object):

    def __init__(self, name = '', cuts = [], **kwargs):
    
        self.logger = logging.getLogger('Region')
        self.logger.info('Initialized for {}'.format(name))

        self.name = name
        self.cuts = cuts

    def make(self):

        # Output Directory
        try:
            os.makedirs(REGION_DIR)
        except OSError:
            if not os.path.isdir(REGION_DIR):
                raise

        # Prepare Processes
        self._check_processes()

        # Parallel Cut
        tasks = mp.Queue()
        results = mp.Queue()

        tmpdir = tf.mkdtemp(prefix = self.name, dir = REGION_DIR)
        self.tmpdir = tmpdir + '/'

        # _processes refers to newly spawned Python processes
        _processes = [
            mp.Process(target = self._cut_process, args = (tasks, results))
            for cpu in xrange(mp.cpu_count())
        ]

        for process in PROCESSES:
            tasks.put(process)

        for p in _processes:
            tasks.put(None)
            p.start()

        for p in _processes:
            p.join()

        results.put(None)

        for r in iter(results.get, None):
            self.logger.info('Selected {!s} out of {!s} entries in {}.'.format(*r))

        # hadd Files
        inputfiles = glob.glob(self.tmpdir + '*.root')
        outputfile = REGION_DIR + self.name + '.root'

        sp.check_call(['hadd', '-f', outputfile] + inputfiles)
        sp.check_call(['rm', '-r', self.tmpdir])

    def _check_processes(self):
        
        for process in PROCESSES:
        
            if os.path.isfile(PROCESS_DIR + process + '.root'):
                continue
            
            self.logger.info('Getting missing process {}'.format(process))
            Process(process, **PROCESSES[process]).make()

    
    def _cut_process(self, tasks = None, results = None):

        for process in iter(tasks.get, None):

            infile = ROOT.TFile(PROCESS_DIR + process + '.root', 'read')
            outfile = ROOT.TFile(self.tmpdir + process + '.root', 'recreate')

            intree = infile.Get('tree')
            intree.SetName(process)
            
            for i, cut in enumerate(self.cuts):
               intree.Draw('>>{0!s}_elist_{1!s}'.format(process,i), cut)
               eventlist = ROOT.gDirectory.Get('{0!s}_elist_{1!s}'.format(process,i))
               intree.SetEventList(eventlist)

            outtree = intree.CopyTree('')

            result = (outtree.GetEntriesFast(), intree.GetEntriesFast(), process)

            outtree.Write()
            
            outfile.Close()
            infile.Close()

            results.put(result)

#------
# Main
#------

if __name__ == '__main__':

    # Set ROOT to run in batch mode.
    ROOT.gROOT.SetBatch(1)

    for name in sys.argv[1:]:

        logging.basicConfig(level = logging.INFO,
                            format = '%(name)s(%(levelname)s) - %(message)s')

        Region(name, **REGIONS[name]).make()

