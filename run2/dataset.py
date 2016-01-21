import glob
import logging
import multiprocessing as mp
import os
import subprocess as sp
import tempfile as tf

import numpy as np
import ROOT

from process import PROCESS_DIR
from settings import WORK_DIR, PROCESSES, CLASSIFICATION


# Output Directory
DATASET_DIR = WORK_DIR + 'classification/datasets/'

class Dataset(object):

    def __init__(self, name = ''):

        self.logger = logging.getLogger('Dataset')

        self.name = name
        self.preselection = CLASSIFICATION.get('preselection', [''])

    def make(self):

        self.logger.info('Generating dataset "{}"'.format(self.name))

        # Output Directory
        try:
            os.makedirs(DATASET_DIR)
        except OSError:
            if not os.path.isdir(DATASET_DIR):
                raise

        # Temporary Work Directory
        tmpdir = tf.mkdtemp(dir = DATASET_DIR)
        self.tmpdir = tmpdir + '/'

        # Parallel Split
        n_tasks = 0
        tasks = mp.Queue()
        results = mp.Queue()

        for process, properties in PROCESSES.iteritems():
            types = set(properties['types'].lower().split(':'))
            if 'data' in types:
                continue
            if 'sig' in types:
                tasks.put((process, 'sig'))
            elif 'bkg' in types:
                tasks.put((process, 'bkg'))
            n_tasks += 1

        _processes = [
            mp.Process(target = self._split_process, args = (tasks, results))
            for i in xrange(min(n_tasks, mp.cpu_count()))
        ]

        for p in _processes:
            tasks.put(None)
            p.start()

        for p in _processes:
            p.join()

        results.put(None)

        for r in iter(results.get, None):
            self.logger.info('{!s} entries in {!s} passed preselection and were split into {!s}({!s}) for training(testing).'.format(*r))

        # hadd Files
        inputfiles = glob.glob(self.tmpdir + '*.root')
        outputfile = DATASET_DIR + self.name + '.root'

        sp.check_call(['hadd', '-f', outputfile] + inputfiles)
        sp.check_call(['rm', '-r', self.tmpdir])

    def _split_process(self, tasks = None, results = None):
    
        for process, types in iter(tasks.get, None):

            infile = ROOT.TFile(PROCESS_DIR + process + '.root', 'read')
            outfile = ROOT.TFile(self.tmpdir + process + '.root', 'recreate')

            intree = infile.Get('tree')
            intree.SetName(process)

            # Apply preselection cuts and shuffle the entries which pass.
            for i, cut in enumerate(self.preselection):
                n_entries = intree.Draw('>>{0!s}_skim_{1!s}'.format(process, i), cut)
                eventlist = ROOT.gDirectory.Get('{0!s}_skim_{1!s}'.format(process, i))
                intree.SetEventList(eventlist)

            entries = np.zeros(n_entries, dtype = np.int64)
            for i in xrange(n_entries):
                entries[i] = eventlist.GetEntry(i)

            #np.random.seed(0)
            np.random.shuffle(entries)
 
            # Split the entries into a training and test set.
            n_test = n_entries / 2
            n_train = n_test if (n_entries % 2 == 0) else n_test + 1

            train_elist = ROOT.TEventList(process + '_train', '', n_train)
            for entry in np.sort(entries[:n_train]):
                train_elist.Enter(entry)

            test_elist = ROOT.TEventList(process + '_test', '', n_test)
            for entry in np.sort(entries[n_train:]):
                test_elist.Enter(entry)

            # Create the training and test set trees.
            intree.SetEventList(train_elist)
            train_tree = intree.CopyTree('')
            train_tree.SetName('{}_{}_train'.format(process, types))
            train_tree.Write()

            intree.SetEventList(test_elist)
            test_tree = intree.CopyTree('')
            test_tree.SetName('{}_{}_test'.format(process, types))
            test_tree.Write()

            # Remove any autosaved TTrees.
            outfile.Delete(process + ';*')

            result = (n_entries, process, train_tree.GetEntriesFast(), test_tree.GetEntriesFast())

            outfile.Close()
            infile.Close()

            results.put(result)
