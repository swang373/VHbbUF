#!/bin/bash

groups=(Data_MET ZnnH125 ggZH125 WlnH125 
        WJets ZJets TT ST QCD VV)

(for group in "${groups[@]}"
do
    echo "$group"
done) | parallel -j2 python step3.py >> step3.log 2>&1

