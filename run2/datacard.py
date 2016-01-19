from operator import itemgetter
import os

from settings import WORK_DIR, PROCESSES


# Output Directory
DATACARD_DIR = WORK_DIR + 'datacards/'

try:
    os.makedirs(DATACARD_DIR)
except OSError:
    if not os.path.isdir(DATACARD_DIR):
        raise

with open(DATACARD_DIR + 'vhbb_Znn_13TeV.txt', 'w') as datacard:

    # The number of observables, or the number of bins in a binned shape fit.
    imax = 1

    # The number of background sources, to be determined.
    jmax = 0

    # Load the processes.
    processes = []

    for process, properties in PROCESSES.iteritems():
        types = set(properties['types'].lower().split(':'))
        if 'pid' in properties:
            processes.append((process, properties['pid']))
        if 'bkg' in types:
            jmax += 1

    processes = sorted(processes, key = itemgetter(1))
    
    
    # Write to datacard.
    datacard.write('imax {!s:<4} number of bins\n'.format(imax))
    datacard.write('jmax {!s:<4} number of backgrounds\n'.format(jmax))
    datacard.write('kmax *    number of nuisance parameters\n')
    datacard.write('---------------------------------------\n')
    datacard.write('shapes * * FIX ME\n')
    datacard.write('---------------------------------------\n')
    datacard.write('bin FIX ME\n')
    datacard.write('observation FIX ME\n')
    datacard.write('---------------------------------------\n')
    datacard.write('bin FIX ME\n')
    datacard.write('process          ' + ''.join('{!s:<10}'.format(p[0]) for p in processes) + '\n')
    datacard.write('process          ' + ''.join('{!s:<10}'.format(p[1]) for p in processes) + '\n')
    
