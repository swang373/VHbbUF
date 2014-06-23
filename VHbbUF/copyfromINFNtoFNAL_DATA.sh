#!/bin/bash

#INFN="srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/arizzi/Step1V33_Step2_V3"
#FNAL="srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/jiafu/ZbbHinv/Step1V33_Step2_V3"
#INFN="srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/arizzi/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V4a_MC/"
#FNAL="srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/jiafu/ZnunuHbb/Step1V42_Step2Tag_EDMV42_Step2_V4a_MC/"
#INFN="srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/arizzi/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V4a_Data/"
#FNAL="srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/jiafu/ZnunuHbb/Step1V42_Step2Tag_EDMV42_Step2_V4a_Data/"
#INFN="srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/arizzi/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V4a_MC_BESTCSV/"
#FNAL="srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/jiafu/ZnunuHbb/Step1V42_Step2Tag_EDMV42_Step2_V4a_MC_BESTCSV/"
INFN="srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/arizzi/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_Data/"
FNAL="srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/lpchbb/apana/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_Data/"

CUT="${INFN#srm*\=}"

#voms-proxy-init -voms cms
lcg-ls -b -D srmv2 $INFN > COPYME_DATA.txt
sed -e "s@$CUT\(.*root\)@$INFN\1 $FNAL\1@" -i COPYME_DATA.txt
#sed -e "s@$CUT\(.*root\)@'$INFN\1' '$FNAL\1'@" -i COPYME_DATA.txt

echo "To start transferring, do:"
echo "srmcp -2 -streams_num=20 -copyjobfile=COPYME_DATA.txt"
