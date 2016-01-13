import glob
import logging
import multiprocessing as mp
import os
import subprocess as sp
import sys
import tempfile as tf

import numpy as np
import ROOT

import settings


EOS = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'

class Sample(object):
    
    def __init__(self, name = '', path = '', xsec = None, **kwargs):

        self.logger = logging.getLogger('Sample')
        self.logger.info('Initialized for {}'.format(name))

        self.name = name
        self.path = path
        if xsec is None:
            self.logger.info('No cross section provided. Assuming data rather than MC sample.')
        else:
            self.xsec = xsec 

    def make(self):

        self.path = self._find_dir()
        self.logger.info('EOS Path: {}'.format(self.path))

        files = self._find_files()
        self.logger.info('Number of ".root" files: {!s}'.format(len(files)))

        # Output Directory
        try:
            os.makedirs(settings.SAMPLE_DIR)
        except OSError:
            if not os.path.isdir(settings.SAMPLE_DIR):
                raise

        # Temporary Work Directory
        tmpdir = tf.mkdtemp(prefix = self.name, dir = settings.SAMPLE_DIR)
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
        outputfile = settings.SAMPLE_DIR + self.name + '.root'

        hadd_log = tf.TemporaryFile(dir = settings.SAMPLE_DIR)
        sp.check_call(['hadd', '-f', outputfile] + inputfiles, stdout = hadd_log, stderr = hadd_log)
        sp.check_call(['rm', '-r', self.tmpdir])

        # Add Sample Luminosity Branch
        if hasattr(self, 'xsec'):
            sample_lumi = self._write_sample_lumi(outputfile)
            self.logger.info('Sample Luminosity: {} pb-1'.format(sample_lumi))

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
            outtree = intree.CopyTree(settings.SKIM)
            
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

    @staticmethod
    def _write_sample_lumi(fname):
        
        infile = ROOT.TFile(fname, 'update')
        tree = infile.Get('tree')

        sample_lumi_address = np.array([-999], np.float32)
        sample_lumi_branch = tree.Branch('sample_lumi', sample_lumi_address, 'sample_lumi/F')
        
        n_pos = infile.Get('CountPosWeight').GetBinContent(1)
        n_neg = infile.Get('CountNegWeight').GetBinContent(1)
        sample_lumi_address[0] = (n_pos - n_neg) / self.xsec

        for i in range(0, tree.GetEntriesFast()):
            tree.GetEntry(i)
            sample_lumi_branch.Fill()

        tree.Write()
        infile.Close()

        return sample_lumi_address[0]

#------
# Main
#------

if __name__ == '__main__':

    # Set ROOT to batch mode.
    ROOT.gROOT.SetBatch(1)

    for name in sys.argv[1:]:

        logging.basicConfig(level = logging.INFO,
                            format = '%(name)s(%(levelname)s) - %(message)s')

        Sample(name, **settings.SAMPLES[name]).make()

