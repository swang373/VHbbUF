#!/bin/sh

DIR=Step4_20130404
mkdir -p $DIR/stitch
cd $DIR/stitch

# ZH
ln -s ../Step4_ZnnH125.root Step4_ZH125.root
# WH
ln -s ../Step4_WlnH125.root Step4_WH125.root
# WjLF + WjHF
hadd -f Step4_Wj.root ../Step4_WJetsPtW70.root ../Step4_WJetsPtW100.root ../Step4_WJetsPtW180.root 
ln -s ../Step4_WJetsHW.root Step4_WjHW.root
# ZjLF + ZjHF
hadd -f Step4_Zj.root ../Step4_ZJetsHT50.root ../Step4_ZJetsHT100.root ../Step4_ZJetsHT200.root ../Step4_ZJetsHT400.root ../Step4_ZJetsPtZ100.root
ln -s ../Step4_ZJetsHW.root Step4_ZjHW.root
# TT
ln -s ../Step4_TTPowheg.root .
hadd -f Step4_TT.root ../Step4_TTFullLeptMG.root ../Step4_TTSemiLeptMG.root ../Step4_TTHadronicMG.root
# s_Top
hadd -f Step4_s_Top.root ../Step4_T_s.root ../Step4_T_t.root ../Step4_T_tW.root ../Step4_Tbar_s.root ../Step4_Tbar_t.root ../Step4_Tbar_tW.root
# VV
ln -s ../Step4_ZZ.root
hadd -f Step4_VV.root ../Step4_WW.root ../Step4_WZ.root ../Step4_ZZ.root
# QCD
hadd -f Step4_QCD.root ../Step4_QCDPt80.root ../Step4_QCDPt120.root ../Step4_QCDPt170.root ../Step4_QCDPt300.root ../Step4_QCDPt470.root ../Step4_QCDPt600.root ../Step4_QCDPt1000.root ../Step4_QCDPt1400.root ../Step4_QCDPt1800.root
# data_obs
hadd -f Step4_data_obs.root ../Step4_Data_R.root ../Step4_Data_P.root

cd -
