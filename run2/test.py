import multiprocessing as mp
import subprocess as sp
import time


with open('/afs/cern.ch/project/eos/installation/cms/etc/setup.sh', 'r') as f:
    for line in f:
        if 'alias eos=' in line:
            eos = line.split('"')[1]

Step1_src = '/store/group/cmst3/user/degrutto/ZnnHbbV13/'
Step1_dir = 'ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/'

# Generate a list of names for the constituent ntuples.
eos_ls = sp.Popen([eos, 'ls', Step1_src + Step1_dir], stdout = sp.PIPE, stderr = sp.PIPE)
out, err = eos_ls.communicate()
ntuples = [x for x in out.split('\n')[:-1]]

parallel_1 = time.time()
pool = mp.Pool(processes=6)
for x in ntuples:
    pool.apply_async(sp.check_call(['cmsStage', Step1_src + Step1_dir + x, '.']))
pool.close()
pool.join()
parallel_2 = time.time()

for x in ntuples:
    sp.check_call(['rm', x])

serial_1 = time.time()
for x in ntuples:
    sp.check_call(['cmsStage', Step1_src + Step1_dir + x, '.'])
serial_2 = time.time()

parallel_time = parallel_2 - parallel_1
serial_time = serial_2 - serial_1

print 'Parallel Time: %s' % parallel_time
print 'Serial Time: %s' % serial_time

