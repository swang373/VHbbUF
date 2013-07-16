#ifndef __CINT__
#include <TROOT.h>
#include <Rtypes.h>
#include <TString.h>
#endif



#include "TSystem.h"
#include <stdlib.h>
#include <stdio.h>
#include "TFile.h"
#ifdef __CINT__
  /// usage getNWZ file.root run

void getNWZ(){


  const char *file = gSystem->Getenv("filenamev1");
  const char *run = gSystem->Getenv("run");
  TFile * root_file = new TFile(file, "read");
  std::string dirW_ewkMuDqm = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuDQM/MT_LASTCUT";
  std::string dirZ_ewkMuDqm = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuDQM/DIMUONMASS_AFTERZCUTS";
  std::string dirW = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuLumiMonitorDQM/TMASS"; 
  std::string dirZ1HLT = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuLumiMonitorDQM/Z_1HLT_MASS"; 
  std::string dirZ2HLT = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuLumiMonitorDQM/Z_2HLT_MASS"; 
  std::string dirZNotIso = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuLumiMonitorDQM/Z_NOTISO_MASS"; 
  std::string dirZMSta = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuLumiMonitorDQM/Z_GLBSTA_MASS";
  std::string dirZMTrk = "/DQMData/Run " + string(run) + "/Physics/Run summary/EwkMuLumiMonitorDQM/Z_GLBTRK_MASS";  
  
  if (!root_file->IsZombie()) {
  TH1D * histoW_ewkMuDqm = root_file->Get(dirW_ewkMuDqm.c_str());
  std::cout << run << " NofW[50-200]fromEWKMuDQMv1 " << histoW_ewkMuDqm->Integral(25,100) << endl;
  TH1D * histoZ_ewkMuDqm = root_file->Get(dirZ_ewkMuDqm.c_str());
  std::cout << run << " NofZ[60-120]fromEWKMuDQMv1 " << histoZ_ewkMuDqm->Integral(30,60) << endl;
  TH1D * histoW = root_file->Get(dirW.c_str());
  std::cout << run << " NofW[50-200]fromEWKMuLumiMonitorDQMv1 " << histoW->Integral(50,200) << endl;
  TH1D * histoZ1HLT = root_file->Get(dirZ1HLT.c_str());
  std::cout << run << " NofZ1HLT[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZ1HLT->Integral(60,120) << endl;
  TH1D * histoZ2HLT = root_file->Get(dirZ2HLT.c_str());
  std::cout << run << " NofZ2HLT[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZ2HLT->Integral(60,120) << endl;
  TH1D * histoZNotIso = root_file->Get(dirZNotIso.c_str());
  std::cout << run << " NofZNotIso[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZNotIso->Integral(60,120) << endl;
  TH1D * histoZMSta = root_file->Get(dirZMSta.c_str());
  std::cout << run << " NofZMuSta[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZMSta->Integral(60,120) << endl;
  TH1D * histoZMTrk = root_file->Get(dirZMTrk.c_str());
  std::cout << run << " NofZMuTrk[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZMTrk->Integral(60,120) << endl;
  }

  const char *file2 = gSystem->Getenv("filenamev2");
   TFile * root_file2 = new TFile(file2, "read");
if  (!root_file2->IsZombie()) {
 TH1D * histoW_ewkMuDqm = root_file2->Get(dirW_ewkMuDqm.c_str());
  std::cout << run << " NofW[50-200]fromEWKMuDQMv1 " << histoW_ewkMuDqm->Integral(25,100) << endl;
  TH1D * histoZ_ewkMuDqm = root_file2->Get(dirZ_ewkMuDqm.c_str());
  std::cout << run << " NofZ[60-120]fromEWKMuDQMv1 " << histoZ_ewkMuDqm->Integral(30,60) << endl;
  TH1D * histoW2 = root_file2->Get(dirW.c_str());
  std::cout << run <<" NofW[50-200]fromEWKMuLumiMonitorDQMv2 " << histoW2->Integral(50,200) << endl;
TH1D * histoZ1HLT = root_file2->Get(dirZ1HLT.c_str());
  std::cout << run << " NofZ1HLT[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZ1HLT->Integral(60,120) << endl;
  TH1D * histoZ2HLT = root_file2->Get(dirZ2HLT.c_str());
  std::cout << run << " NofZ2HLT[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZ2HLT->Integral(60,120) << endl;
  TH1D * histoZNotIso = root_file2->Get(dirZNotIso.c_str());
  std::cout << run << " NofZNotIso[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZNotIso->Integral(60,120) << endl;
  TH1D * histoZMSta = root_file2->Get(dirZMSta.c_str());
  std::cout << run << " NofZMuSta[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZMSta->Integral(60,120) << endl;
  TH1D * histoZMTrk = root_file2->Get(dirZMTrk.c_str());
  std::cout << run << " NofZMuTrk[60-120]fromEWKMuLumiMonitorDQMv1 " << histoZMTrk->Integral(60,120) << endl;
 }
  

}

 


