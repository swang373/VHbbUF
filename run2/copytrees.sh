#!/bin/bash

dir="/store/cmst3/user/degrutto/ZnnHbbV13skimmedMet150/MET__Run2015C-PromptReco-v1/"

for i in $(eos ls $dir) 
do
    if [ "$i" == "log" ]
    then
        continue
    fi

    cmsStage $dir$i .

done
