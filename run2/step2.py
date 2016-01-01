import logging
import multiprocessing as mp
import os
import subprocess as sp
import tempfile as tf

import numpy as np
import ROOT

import settings


EOS = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'

class Step2(object):
    
    def __init__(self, name = '', path = '', cuts = [], xsec = None):

        self.logger = logging.getLogger('Step2')
        self.logger.info('Initialized for {}'.format(name))

        self.name = name
        self.path = path
        self.cuts = '&&'.join('({})'.format(x) for x in SKIM + cuts)
        if xsec is None:
            self.logger.info('No cross section provided, assuming data sample')
        else:
            self.xsec = xsec 

    def run(self):

        self.path = self._find_dir()
        self.logger.info('EOS Path: {}'.format(self.path))

        files = self._find_files()
        self.logger.info('Number of ".root" files: {!s}'.format(len(files)))

        # Output Directory
        try:
            os.makedirs(settings.STEP2_DIR)
        except OSError:
            if not os.path.isdir(settings.STEP2_DIR):
                raise

        tmpdir = tf.mkdtemp(prefix = self.name, dir = settings.STEP2_DIR)
        self.tmpdir = tmpdir + '/'

        # Parallel Copy
        tasks = mp.Queue()
        results = mp.Queue()
 
        processes = [
            mp.Process(target = self._copy_file, args = (tasks, results))
            for cpu in range(mp.cpu_count())
        ]

        for f in files:
            tasks.put(f)

        for p in processes:
            tasks.put(None)
            p.start()
        
        for p in processes:
            p.join()

        results.put(None)

        for r in iter(results.get, None):
            self.logger.info('{} skimmed from {!s} to {!s} entries.'.format(*r))

        # hadd Files
        inputfiles = glob.glob(self.tmpdir + '*.root')
        outputfile = settings.STEP2DIR + self.name + '.root'

        hadd_log = tf.TemporaryFile(dir = settings.STEP2_DIR)
        sp.check_call(['hadd', '-f', outputfile] + inputfiles, stdout = hadd_log, stderr = hadd_log)
        sp.check_call(['rm', '-r', self.tmpdir])

    def _find_dir(self):
        
        eos_find_d = sp.check_output([EOS, 'find', '-d', self.path])
 
        dirtree = (line.split() for line in eos_find_d.splitlines())
 
        for level in dirtree:
            path, n_dir, n_files = (int(x.split('=')[1]) if '=' in x else x for x in level)
            if n_dir > 1:
                raise ValueError, 'Unable to find a unique path: "{}" contains {!s} subdirectories'.format(path, n_dir)
        else:
            return path
                
    def _find_files(self):

        eos_ls = sp.check_output([EOS, 'ls', self.path])
       
        files = [line for line in eos_ls.splitlines() if line.endswith('.root')]
        
        if files:
            return files
        else:
            raise ValueError, 'Unable to find any ".root" files'

    def _copy_file(self, tasks = None, results = None):

        xrd_prefix = 'root://eoscms.cern.ch/'
 
        for fname in iter(tasks.get, None):

            infile = ROOT.TFile.Open(xrd_prefix + self.path + fname, 'read')
            outfile = ROOT.TFile(self.tmpdir + fname, 'recreate')

            intree = infile.Get('tree')
            outtree = intree.CopyTree(self.cuts)
            
            result = (fname, intree.GetEntriesFast(), outtree.GetEntriesFast())
            
            outtree.Write()

            for key in infile.GetListOfKeys():
                if key.GetName() == 'tree':
                    continue
                obj = key.ReadObj()
                obj.Write()
           
            outfile.Close()
            infile.Close()

            results.put(result)

#------
# Main
#------

if __name__ == '__main__':

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    for name in sys.argv[1:]:

        logging.basicConfig(level = logging.INFO,
                            format = '%(name)s(%(levelname)s) - %(message)s')

        step2 = Step2(name, **settings.SAMPLES[name])
        step2.run()
 


