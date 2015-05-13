#include "plotHistos_Wen.h"

#include "TROOT.h"
#include "TStyle.h"


// To run, do "root -l plotHistos_Wen.C+"

void plotHistos_Wen() {
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    TH1::SetDefaultSumw2(1);
    //gROOT->SetBatch(1);

    TString plotdir = "plots/";

    if (gSystem->AccessPathName(plotdir))
        gSystem->mkdir(plotdir);


    ////////////////////////////////////////////////////////////////////////////
    // Task 1 (a)                                                             //
    // - please comment out other tasks                                       //
    ////////////////////////////////////////////////////////////////////////////

    // Read from ntuples
    Events * ev = new Events();
    TCut cutmc_all   = "Vtype==3";   // change this to your channel
    TCut cutdata_all = "Vtype==3";   // change this to your channel
    TString process  = "VH";         // can try other processes: VH, WJ, ZJ, TT
    ev->read(cutmc_all, cutdata_all, process);

    TString var      = "H.mass";                // the variable to plot
    TCut cut         = "V.pt>120";              // the selection cut
    TString title    = ";m(jj) [GeV]";          // the title of the histogram
    TString plotname = process + "_Hmass";      // the name of the image file
    int nbinsx       = 15;                      // number of bins
    double xlow      = 30.0;                    // the low edge of x-axis
    double xup       = 255.0;                   // the upper edge of x-axis
    TString options  = "!plotLog:!plotNorm";    // use "plotLog" to plot on log-y scale,
                                                // use "!plotLog" to plot on linear scale;
                                                // use "plotNorm" to plot normalized plots,
                                                // use "!plotNorm" otherwise.

    // Use "ev->WH" for WH, or "ev->WjLF", "ev->WjHF", "ev->ZjLF", "ev->ZjHF", "ev->TT" for other processes
    MakePlot(ev->WH, var, cut, title, nbinsx, xlow, xup, plotname, plotdir, options);

    // You can put in the parameters directly as in the following commented out line:
    //MakePlot(ev->WH, "H.pt", cut, "; p_{T}(jj) [GeV]", 16, 0, 240., process + "_Hpt", plotdir, options);


    ////////////////////////////////////////////////////////////////////////////
    // Task 1 (b)                                                             //
    // - please comment out other tasks                                       //
    ////////////////////////////////////////////////////////////////////////////
/*
    // Read from ntuples
    Events * ev = new Events();
    TCut cutmc_all   = "Vtype==3";   // change this to your channel
    TCut cutdata_all = "Vtype==3";   // change this to your channel
    ev->read(cutmc_all, cutdata_all, "VH:WJ");  // read both VH and WJ processes

    TString var      = "H.mass";
    //TCut cut         = "V.pt>120";
    TCut cut         = "V.pt>120 && hJet_pt[0]>30 && hJet_pt[1]>30";  // for Wln, Znn, change to tighter cut
    TString title    = ";m(jj) [GeV]";
    TString plotname = "WH_vs_WJ_Hmass";
    int nbinsx       = 15;
    double xlow      = 30.0;
    double xup       = 255.0;
    TString options  = "!plotLog:plotNorm";

    // Using "ev->WH" for WH and "ev->WjHF" for W+HF
    MakePlot2(ev->WH, ev->WjHF, var, cut, title, nbinsx, xlow, xup, plotname, plotdir, options);
*/

    ////////////////////////////////////////////////////////////////////////////
    // Task 2                                                                 //
    // - please comment out other tasks                                       //
    ////////////////////////////////////////////////////////////////////////////
/*
    // Zmm______________________________________________________________________
    //TString channel  = "Zmm";

    // These are loose cuts for all plots in this particular channel
    //TCut cutmc_all   = "Vtype==0 && V.pt>120 && hJet_pt[0]>20 && hJet_pt[1]>20 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && vLepton_pt[0]>20 && vLepton_pt[1]>20 && abs(vLepton_eta[0])<2.4 && abs(vLepton_eta[1])<2.4 && METtype1corr.et<60 && 75<V.mass && V.mass<105 && H.dR<1.6 && hbhe==1";
    //cutmc_all       += "min(hJet_csv_nominal[0], hJet_csv_nominal[1])>0.4";  // tighter cut
    //TCut cutdata_all = cutmc_all;
    //cutmc_all       *= "weightTrig2012";  // apply trigger weight for MC

    // Scale factors in order of: WjLF, WjHF, ZjLF, ZjHF, TT
    // NOTE: WjLF, WjHF are not needed for Zll
    //double scalefactors[5] = {1.00, 1.00, 1.00, 1.00, 1.00};


    // Zee______________________________________________________________________
    //TString channel  = "Zee";

    // These are loose cuts for all plots in this particular channel
    //TCut cutmc_all   = "Vtype==1 && V.pt>120 && hJet_pt[0]>20 && hJet_pt[1]>20 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && vLepton_pt[0]>20 && vLepton_pt[1]>20 && abs(vLepton_eta[0])<2.5 && abs(vLepton_eta[1])<2.5 && METtype1corr.et<60 && 75<V.mass && V.mass<105 && H.dR<1.6 && hbhe==1";
    //cutmc_all       += "min(hJet_csv_nominal[0], hJet_csv_nominal[1])>0.4";  // tighter cut
    //TCut cutdata_all = cutmc_all;
    //cutmc_all       *= "weightTrig2012";  // apply trigger weight for MC

    // Scale factors in order of: WjLF, WjHF, ZjLF, ZjHF, TT
    // NOTE: WjLF, WjHF are not needed for Zll
    //double scalefactors[5] = {1.00, 1.00, 1.00, 1.00, 1.00};


    // Wmn______________________________________________________________________
    //TString channel  = "Wmn";

    // These are loose cuts for all plots in this particular channel
    //TCut cutmc_all   = "Vtype==2 && H.pt>80 && hJet_pt[0]>30 && hJet_pt[1]>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && vLepton_pt[0]>30 && abs(vLepton_eta[0])<2.4 && METtype1corr.et >45 && nalep==0 && Sum$(aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)==0 && hbhe==1";
    //cutmc_all       += "min(hJet_csv_nominal[0], hJet_csv_nominal[1])>0.4 && abs(HVdPhi)>2.0";  // tighter cut
    //TCut cutdata_all = cutmc_all;
    //cutmc_all       *= "weightTrig2012";  // apply trigger weight for MC

    // Scale factors in order of: WjLF, WjHF, ZjLF, ZjHF, TT
    // NOTE: ZjLF, ZjHF are not needed for Wln
    //double scalefactors[5] = {1.00, 1.00, 1.00, 1.00, 1.00};


    // Wen______________________________________________________________________
    TString channel  = "Wen";

    // These are loose cuts for all plots in this particular channel
    TCut cutmc_all   = "Vtype==3 && H.pt>80 && hJet_pt[0]>30 && hJet_pt[1]>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && vLepton_pt[0]>30 && abs(vLepton_eta[0])<2.5 && METtype1corr.et >45 && nalep==0 && Sum$(aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)==0 && hbhe==1";
    cutmc_all       += "min(hJet_csv_nominal[0], hJet_csv_nominal[1])>0.4 && abs(HVdPhi)>2.0";  // tighter cut
    TCut cutdata_all = cutmc_all;
    cutmc_all       *= "weightTrig2012";  // apply trigger weight for MC

    // Scale factors in order of: WjLF, WjHF, ZjLF, ZjHF, TT
    // NOTE: ZjLF, ZjHF are not needed for Wln
    double scalefactors[5] = {1.00, 1.00, 1.00, 1.00, 1.00};


    // Znn______________________________________________________________________
    //TString channel  = "Znn";

    // These are loose cuts for all plots in this particular channel
    //TCut cutmc_all   = "Vtype==4 && H.pt>130 && hJet_pt[0]>80 && hJet_pt[1]>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && nalep==0 && Sum$(aJet_pt>25 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)==0 && min(Min$(abs(deltaPhi(METtype1corr.phi,hJet_phi))),Min$(abs(deltaPhiMETjets(METtype1corr.phi,aJet_phi,aJet_pt,aJet_eta)))+999*(Sum$(aJet_pt>25 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)==0) )>0.5 && hbhe==1";
    //cutmc_all       += "min(hJet_csv_nominal[0], hJet_csv_nominal[1])>0.4 && abs(HVdPhi)>2.0";  // tighter cut
    //TCut cutdata_all = cutmc_all;
    //cutmc_all       *= "(triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1) && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent"; // apply trigger bits and MET cleaning for MC (they are already applied on data)
    //cutmc_all       *= "triggercorrMET(METtype1corr.et)";

    // Scale factors in order of: WjLF, WjHF, ZjLF, ZjHF, TT
    //double scalefactors[5] = {1.00, 1.00, 1.00, 1.00, 1.00};


    // All channels_____________________________________________________________
    // Read from ntuples
    Events * ev = new Events();
    ev->read(cutmc_all, cutdata_all);

    // Set the scale factors
    ev->set_sf(scalefactors);

    // Optimize these five variables (default: recommendations for Zll)
    //double vpt    = 150.;
    //double hpt    = 0.;
    //double maxcsv = 0.679;
    //double mincsv = 0.5;
    //double dPhi   = 0.;

    // Optimize these five variables (default: recommendations for Wln)
    double vpt    = 150.;
    double hpt    = 100.;
    double maxcsv = 0.898;
    double mincsv = 0.5;
    double dPhi   = 2.95;

    // Optimize these five variables (default: recommendations for Znn)
    //double vpt    = 170.;  // for Znn, pT(V) = MET
    //double hpt    = 170.;
    //double maxcsv = 0.898;
    //double mincsv = 0.5;
    //double dPhi   = 2.95;

    // If doing cut and count analysis, cut on H.mass by changing the values of minhmass and maxhmass
    //double minhmass = 0.;
    //double maxhmass = 9999.;
    double minhmass = 110.;
    double maxhmass = 140.;

    // These are tight cuts for this particular plot
    TCut cutmc = Form("V.pt>%.2f && H.pt>%.2f && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>%.3f && min(hJet_csv_nominal[0], hJet_csv_nominal[1])>%.3f && abs(HVdPhi)>%.2f && %.2f<H.mass && H.mass<%.2f", vpt, hpt, maxcsv, mincsv, dPhi, minhmass, maxhmass);
    TCut cutdata = cutmc;

    TString var      = "H.mass";
    TString title    = ";m(jj) [GeV]";
    TString plotname = channel + "_Hmass";
    int nbinsx       = 15;
    double xlow      = 30.0;
    double xup       = 255.0;
    TString options  = "printStat:plotSig:!plotData:!plotLog";

    MakePlots(ev, var, cutmc, cutdata, title, nbinsx, xlow, xup, plotname, plotdir, options);

    // Or, just put in them directly as in the following commented out line:
    //MakePlots(ev, "H.mass", cutmc, cutdata, "m(jj) [GeV]", 15, 30.0, 255.0, channel+"_Hmass", plotdir, "printStat:plotSig:!plotData:!plotLog");
*/

    ////////////////////////////////////////////////////////////////////////////
    // Task 3                                                                 //
    // - please comment out other tasks, but keep Task 2                      //
    ////////////////////////////////////////////////////////////////////////////
/*
    TString dcname    = Form("vhbb_%s_8TeV.txt", channel.Data());   // the datacard name
    TString wsname    = plotdir + plotname +".root";                // the workspace name
    bool    useshapes = false;
    TString options1  = "!unblind:SplusB";

    // For cut-and-count analysis, apply H.mass cut before calling MakeDatacard(...)
    MakeDatacard(channel, dcname, wsname, useshapes, options1);


    // For shape analysis, remove H.mass cut before calling MakeDatacard(...)
    //cutmc = Form("V.pt>%.2f && H.pt>%.2f && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>%.3f && min(hJet_csv_nominal[0], hJet_csv_nominal[1])>%.3f && abs(HVdPhi)>%.2f", vpt, hpt, maxcsv, mincsv, dPhi);
    //cutdata = cutmc;
    //plotname = channel + "_Hmass_shapes";
    //MakePlots(ev, var, cutmc, cutdata, title, nbinsx, xlow, xup, plotname, plotdir, options);

    //dcname    = Form("vhbb_shapes_%s_8TeV.txt", channel.Data());    // the datacard name
    //wsname    = plotdir + plotname +".root";                        // the workspace name
    //useshapes = true;
    //options1  = "unblind:SplusB";
    //MakeDatacard(channel, dcname, wsname, useshapes, options1);
*/

}
