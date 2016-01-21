import logging
import multiprocessing as mp
import time

from scipy import stats
import numpy as np

from BDT import BDT


HYPERPARAM_SPACE = {
    'NTrees': stats.randint(200, 1001),
    'MaxDepth': stats.randint(3, 11),
    'MinNodeSize': 0.05,
    'nCuts': stats.randint(10, 51),
    'BoostType': 'AdaBoost',
    'AdaBoostBeta': 0.5,
    'SeparationType': ['CrossEntropy', 'GiniIndex', 'GiniIndexWithLaplace', 'MisClassificationError'],
    #'Shrinkage': stats.expon(loc = 0.001, scale = 0.1),
}

def sample_hyperparams():
    hyperparams = {}
    for name, value in HYPERPARAM_SPACE.iteritems():
        if hasattr(value, 'rvs'):
            hyperparams[name] = value.rvs()
        elif isinstance(value, list):
            hyperparams[name] = value[np.random.randint(len(value))]
        else:
            hyperparams[name] = value
    return hyperparams

class Worker(mp.Process):

    """
    github.com/qbuat/tauperf/blob/master/tauperf/parallel.py
    Thank you, Quentin! Otherwise I'd be stuck in TMVA hell.
    """

    def __init__(self, func, *args, **kwargs):
        super(Worker, self).__init__()
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.result = mp.Queue()

    def run(self):
        self.func(*self.args, **self.kwargs)

#------
# Main
#------

if __name__ == '__main__':

    # Logging setup.
    logging.basicConfig(level = logging.INFO,
                        format = '%(name)s(%(levelname)s) - %(message)s')

    logger = logging.getLogger('RandomSearch')

    # Random Search Configuration
    n_trials = 20

    # BDT Configuration
    options = {
        'dataset': '20Jan2016',
    }

    # Initialize Random Search Trials
    bdts = []
    for i in xrange(n_trials):
        options['job_id'] = 'RandomSearch_Trial{!s}'.format(i)
        options['hyperparams'] = sample_hyperparams()
        bdts.append(BDT(**options))
    trials = [Worker(bdt.train) for bdt in bdts]

    # Parallel Training
    n_processes = min(n_trials, mp.cpu_count())
    _processes = []
    p = None

    try:
        while True:
            n_active = len(mp.active_children())
            while n_active < n_processes and len(trials) > 0:
                p = trials.pop(0)
                p.start()
                _processes.append(p)
                n_active = len(mp.active_children())
            if len(trials) == 0 and n_active == 0:
                break
            time.sleep(0.1)
    except KeyboardInterrupt, SystemExit:
        if p is not None:
            p.terminate()
        for p in _processes:
            p.terminate()
        raise

