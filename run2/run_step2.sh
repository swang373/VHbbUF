#!/bin/bash

labels=(MET_RunC MET_RunD MET_RunD_Prompt 
        ZnnH125 ggZH125 WlnH125 
        WJetsHT100 WJetsHT200 WJetsHT400 WJetsHT600 WJetsIncl
        ZJetsHT100 ZJetsHT200 ZJetsHT400 ZJetsHT600 
        TTPow T_tW Tbar_tW T_s_comb_lep T_t_lep Tbar_t_lep
        QCDHT100 QCDHT200 QCDHT300 QCDHT500 QCDHT700 QCDHT1000 QCDHT1500 QCDHT2000
        WW WZ ZZ)

(for label in "${labels[@]}"
do
    echo "$label"
done )| parallel -j2 python step2.py >> step2.log 2>&1

