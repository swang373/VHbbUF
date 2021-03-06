import time
import sys,os
import getopt
import subprocess

usage = """
maketiert2list.py [options]
--input-dir, -i                    origin directory
--wildcard, -w                     selects only directories with w to transfer
--source, -s                       pisa,fnal,cern
--nolsmode, -n                     don't check if files already exist, default off
--checkmode, -c                    check if all files are transferred, default off
"""

pisapre='srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/'
cernpre='srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms/store/user/'
fnalpre='srm://cmseos.fnal.gov:8443/srm/v2/server?SFN=/eos/uscms/store/user/'

sourcepre=''

source=''

inputDir=''

wildcard=''

checkmode=False
nolsmode=False

(opts, args) = getopt.getopt(sys.argv[1:], 'd:i:o:w:s:hcn', ['source=', 'input-dir=',  'help', 'wildcard=','checkmode'])

for opt,argm in opts:
    #print opt,argm
    if (opt == "--help" or opt == "-h"):
        print 'Usage: %s' % (usage)
        sys.exit(0)
    elif (opt == "-s" or opt == "--source"):
        source = argm

    elif (opt == "-i" or opt == "--input-dir"):
        inputDir = argm
    elif (opt == "-w" or opt == "--wildcard"):
        wildcard = argm
    elif (opt == "-c" or opt == "--checkmode"):
        checkmode=True
    elif (opt == "-n" or opt == "--nolsmode"):
        nolsmode=True
    else:
        print 'Wrong options: %s' % (opt)
        sys.exit(3)

if source.find('pisa') != -1:
    sourcepre=pisapre
elif source.find('cern') != -1:
    sourcepre=cernpre
elif source.find('fnal') != -1:
    sourcepre=fnalpre
else:
    print "source",source,"is not setup"
    sys.exit(0)



def LCG_LS_WRAP(args):
    cmd=['lcg-ls','-b','-D','srmv2']
    cmd.extend(args)
    checkOutput = subprocess.check_output(cmd)
    return checkOutput


def LCG_CP_WRAP(input, output):
    cmd=['lcg-cp','-b','-D','srmv2',input,output]
    #print cmd
    subprocess.Popen(cmd)


def DoesFileExist(inpath,outpath):
    try:
        #print inpath
        outls=LCG_LS_WRAP(["-l",outpath])
        #check that file is the right size
        inls=LCG_LS_WRAP(["-l",inpath])
        #print "insize, outsize,",inls.split()[4],",",outls.split()[4]
        if inls.split()[4] == outls.split()[4]:
            #same size
            return True
        else:
            return False
    except Exception as e:
        #print "ERROR",e
        return False


def MakeFileList(rootpath):
    #print rootpath
    output=LCG_LS_WRAP([rootpath])  
    files=[]
    list=output.split("\n")
    for item in list:
        if item is not "":
            item=item.split("user/")[-1]
            # removing things here could be problematic
            # if there are many subdirs without datasetnames
            if wildcard is not "":
                if item.find(wildcard) == -1:
                    print "skipping",item
                    continue
            if item.find(".root") != -1:
                files.append(item)
            else:
                try:
                    files.extend(MakeFileList(sourcepre+item))
                except:
                    print "The following item fails "+item
    return files


def SetMaxThreads():
    cmd=['cat','/proc/cpuinfo']
    output=subprocess.check_output(cmd)
    corelines=[line for line in output.split('\n') if line.find("cpu cores") != -1]
    cores=0
    for coreline in corelines:
        cores=cores+int(coreline.split()[-1])

    return cores


MaxThreads=SetMaxThreads()

#inputDir='/cms/store/user/arizzi/VHBBHeppyTest/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/'

print "sourcepre",sourcepre
print "inputDir",inputDir

files=MakeFileList(sourcepre+inputDir)


ifile=0
isleep=0
newsleep=1
untransfiles=[]

if len(files) is 0:
    print "There are no files to transfer."
    sys.exit(0)

print "MaxThreads files ",MaxThreads,len(files)

file = open('myfile.dat', 'w+')
for f in files:
   file.write(f)
   file.write("\n")


#while ifile!=len(files):
#    out=[line for line in subprocess.check_output(["ps"]).split("\n") if line.find("lcg-cp") !=-1 ]
#   
#    #if ifile>13:
#    #    break

#    if len(out)< MaxThreads:
#        if files[ifile].find("user/") is not -1:
#            files[ifile]=files[ifile].split("user/")[-1]
#        outpath=files[ifile].replace(inputDir,"")
#        if dest is "local":
#            outputDir=outpath.replace(outpath.split("/")[-1], "")
#            if not os.path.exists(outputDir):
#                os.makedirs(outputDir)
#        outpath=destpre+outputDir+outpath
#        inpath=sourcepre+files[ifile]#

#        transfer=True
#        if nolsmode is False:
#            transfer= not DoesFileExist(inpath, outpath)
        
#        if transfer is True:
#            #print "File does not exist", outpath
#            if checkmode is True:
#               untransfiles.append(outpath) 
#            else:
#                #LCG_CP_WRAP(inpath, outpath)
#                newsleep=1
#        
#        ifile=ifile+1
#    else:
#        isleep=isleep+1
#        if newsleep == 1:
#            print "sleeping",isleep
#            print "ifile",ifile,"of",len(files)
#            newsleep=0
#        time.sleep(2)

#if checkmode is True:
#    print "Missing:"
#    for file in untransfiles:
#        print file
#    print len(untransfiles),"of",len(files),"to go"
      
#else:
#    out=[line for line in subprocess.check_output(["ps"]).split("\n") if line.find("lcg-cp") !=-1 ]
#    print "Number of on-going lcg-cps:",len(out)
