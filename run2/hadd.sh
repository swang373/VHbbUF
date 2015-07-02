#!/bin/bash

hadd -f Step3_QCD.root   Step3_QCDHT100.root   Step3_QCDHT250.root   Step3_QCDHT500.root   Step3_QCDHT1000.root 
hadd -f Step3_s_Top.root Step3_T_s.root        Step3_Tbar_s.root     Step3_T_t.root        Step3_Tbar_t.root   Step3_T_tW.root Step3_Tbar_tW.root
hadd -f Step3_WJets.root Step3_WJetsHT100.root Step3_WJetsHT200.root Step3_WJetsHT400.root Step3_WJetsHT600.root
hadd -f Step3_ZJets.root Step3_ZJetsHT100.root Step3_ZJetsHT200.root Step3_ZJetsHT400.root Step3_ZJetsHT600.root
