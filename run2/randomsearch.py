import logging
import multiprocessing as mp

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

def worker(options):
    bdt = BDT(**options)
    return bdt.train()

#------
# Main
#------

if __name__ == '__main__':

    # Logging setup.
    logging.basicConfig(level = logging.INFO,
                        format = '%(name)s(%(levelname)s) - %(message)s')

    logger = logging.getLogger('RandomSearch')

    # Random Search Options
    n_trials = 5

    options = {
        'job_name': 'test',
        'dataset': '20Jan2016',
    }

    # Parallel Training
    pool = mp.Pool(min(n_trials, mp.cpu_count()))
    for i in xrange(n_trials): 
        options['name'] = 'BDT_Trial{!s}'.format(i)
        options['hyperparams'] = sample_hyperparams()
        pool.apply_async(worker, (options,))
    pool.close()
    pool.join()

