#include "plotHistos_Znn_13TeV.h"
#include "TROOT.h"
#include "TStyle.h"

/*
inline double evalHMETMassiveMt(double m0, double pt0, double phi0, double m1, double pt1, double phi1)
{

  return TMath::Sqrt(m0*m0 + m1*m1 + 2.0 * (TMath::Sqrt((m0*m0+pt0*pt0)*(m1*m1+pt1*pt1)) - pt0*pt1*TMath::Cos(deltaPhi(phi0, phi1)) ));

}


#include "TLorentzVector.h"
inline double evalHMETPt(double m0, double pt0, double phi0, double eta0, double m1, double pt1, double phi1, double eta1)

{

  TLorentzVector v1,v2;
  v1.SetPtEtaPhiM(pt0,eta0, phi0, m0);
  v2.SetPtEtaPhiM(pt1,eta1, phi1, m1);
  return (v1+v2).Pt();

}
*/




// To run, do "root -l plotHistos_Znn.C+"

void plotHistos_Znn_13TeV() {
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
    //    TCut cutmc_all   = "Vtype==4 && Sum$(hJets_btagCSV>0.9)>0 && H_mass<250 && Jet_pt[hJidx[1]] > 30";
   // change this to your channel
    //    TCut cutdata_all = "Vtype==4 && Sum$(hJets_btagCSV>0.9)>0 && H_mass<250 && Jet_pt[hJidx[1]] > 30";   // change this to your channel

    //     TCut cutmc_all = " Vtype==4 && Sum$(hJets_btagCSV>0.7)>0 && H_mass<250 && Jet_pt[hJidx[1]] > 30 && met_pt>200 && H_pt>150 && min(hJets_btagCSV[0], hJets_btagCSV[1])>0.7 && Jet_pt[hJidx[0]]>100 && Jet_pt[hJidx[1]]>80   && deltaR_jj<1.3 && Sum$(Jet_pt>30 && abs(Jet_eta)<2.5)== 2";

    //    TCut cutmc_all = "Vtype==4 &&  Jet_btagCSV[hJCidx[0]]>0.9 && Jet_btagCSV[hJCidx[1]]>0.8  && HCSV_mass<250 &&  met_pt>200 && HCSV_pt>0 && Jet_pt[hJCidx[0]]>80 && Jet_pt[hJCidx[1]]>30   && deltaR( Jet_eta[hJCidx[0]], Jet_eta[hJCidx[1]],  Jet_phi[hJCidx[0]], Jet_phi[hJCidx[1]] )<1.3 &&  Sum$(Jet_pt>30 && abs(Jet_eta)<2.5)<=3" ;



      TCut cutmc_all = "Vtype==4 &&  Jet_btagCSV[hJCidx[0]]>0.814 && Jet_btagCSV[hJCidx[1]]>0.423   && HCSV_mass<250 &&  met_pt>170 && HCSV_pt>130 && ( (Jet_pt[hJCidx[0]]>80 && Jet_pt[hJCidx[1]]>30) || (Jet_pt[hJCidx[1]]>80 && Jet_pt[hJCidx[0]]>30 ))    && deltaR( Jet_eta[hJCidx[0]], Jet_phi[hJCidx[0]], Jet_eta[hJCidx[1]],  Jet_phi[hJCidx[1]] )<3.2 &&  Sum$(aLeptons_pt>5  && aLeptons_relIso03<1. &&  ( deltaR(Jet_eta[hJCidx[0]], aLeptons_eta, Jet_phi[hJCidx[0]], aLeptons_phi)>0.4 &&  deltaR( Jet_eta[hJCidx[1]], aLeptons_eta, Jet_phi[hJCidx[1]], aLeptons_phi)>0.4))==0 && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId==1)<=99 && MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5 && abs(deltaPhi(HCSV_phi, met_phi))>2.5";

      //TCut cutmc_all = "Vtype==4 &&  Jet_btagCSV[hJCidx[0]]>0.9 && Jet_btagCSV[hJCidx[1]]>0.68   && HCSV_mass<250 &&  met_pt>200 && HCSV_pt>130 && ( (Jet_pt[hJCidx[0]]>80 && Jet_pt[hJCidx[1]]>30) || (Jet_pt[hJCidx[1]]>80 && Jet_pt[hJCidx[0]]>30 ))    && deltaR( Jet_eta[hJCidx[0]], Jet_phi[hJCidx[0]], Jet_eta[hJCidx[1]],  Jet_phi[hJCidx[1]] )<1.3 &&  Sum$(aLeptons_pt>5  && aLeptons_relIso03<1. &&  ( deltaR(Jet_eta[hJCidx[0]], aLeptons_eta, Jet_phi[hJCidx[0]], aLeptons_phi)>0.4 &&  deltaR( Jet_eta[hJCidx[1]], aLeptons_eta, Jet_phi[hJCidx[1]], aLeptons_phi)>0.4))==0 && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId==1)>=2  && abs(deltaPhi(HCSV_phi, met_phi))>0.90";

    //    TCut cutmc_all = "";

    // TCut cutdata_all = "Vtype==4 && Sum$(hJets_btagCSV>0.7)>0 && H_mass<250 && Jet_pt[hJidx[1]] > 30 && met_pt>200 && H_pt>150 && min(hJets_btagCSV[0], hJets_btagCSV[1])>0.7 && Jet_pt[hJidx[0]]>100 && Jet_pt[hJidx[1]]>80   && deltaR_jj<1.3 && Sum$(Jet_pt>30 && abs(Jet_eta)<2.5)== 2";

    TCut cutdata_all = cutmc_all;

    ev->read(cutmc_all, cutdata_all);




    TString var      = "HCSV_mass";                // the variable to plot
    TCut cut         = "HCSV_mass<150 && HCSV_mass>90";              // the selection cut
    TString title    = ";Hmass [GeV]";          // the title of the histogram
    TString plotname = "All_hmass_tightMass";      // the name of the image file
    int nbinsx       = 6;                      // number of bins
    double xlow      = 90;                    // the low edge of x-axis
    double xup       = 150;                   // the upper edge of x-axis
    TString options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";    // use "!plotLog" to plot on log-y scale,
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);



    var      = "HCSV_mass";                // the variable to plot
    cut         = "HCSV_mass>70 && HCSV_mass<170";              // the selection cut
    title    = ";Hmass [GeV]";          // the title of the histogram
    plotname = "All_hmass_tightMassNorm";      // the name of the image file
    nbinsx       = 5;                      // number of bins
    xlow      = 70;                    // the low edge of x-axis
    xup       = 170;                   // the upper edge of x-axis
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";    // use "!plotLog" to plot on log-y scale,
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);




    nbinsx       = 17;                      // number of bins
    xlow      = 0;                    // the low edge of x-axis
    xup       = 255.0;                   // the upper edge of x-axis 
    cut = "";
    var      = "HCSV_mass";                // the variable to plot
    title    = ";Hmass [GeV]";          // the title of the histogram
    plotname = "All_hmass_tightcuts";      // the name of the image file
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);




    nbinsx       = 5;                      // number of bins
    xlow      = 70.0;                    // the low edge of x-axis
    xup       = 170.0;                   // the upper edge of x-axis 
    var      = "HCSV_mass";                // the variable to plot
    title    = ";Hmass [GeV]";          // the title of the histogram
    plotname = "All_hmass_tightMass";      // the name of the image file
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);







    var      = "met_pt";                // the variable to plot
    cut         = "";              // the selection cut
    title    = ";MET [GeV]";          // the title of the histogram
    plotname = "metpt";      // the name of the image file
    nbinsx       = 7;                      // number of bins
    xlow      = 200.0;                    // the low edge of x-axis
    xup       = 550.0;                   // the upper edge of x-axis
    options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";    // use "!plotLog" to plot on log-y scale,





    plotname = "All_metpt_norm";
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);


    plotname = "All_metpt";
    options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);
    


    
    plotname = "All_metpt_tightcuts";
    options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";
    //    cut         = "Jet_pt[hJCidx[0]]>50  && Jet_pt[hJCidx[1]] > 50  && abs(Jet_eta[hJCidx[0]])<2.4 &&  abs(Jet_eta[hJCidx[1]])<2.4  &&  (Jet_btagCSV[hJCidx[0]]>0.7 &&  Jet_btagCSV[hJCidx[1]]>0.7)  && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5)<=2 && met_pt>200";              // the selection cut
    //    cut         = "Jet_pt[hJCidx[0]]>50    && abs(Jet_eta[hJCidx[0]])<2.4 &&  abs(Jet_eta[hJCidx[1]])<2.4  && met_pt>200";              // the selection cut
    cut         = "";
    nbinsx       = 6;                      // number of bins
    xlow      = 200.0;                    // the low edge of x-axis
    xup       = 500.0;                   // the upper edge of x-axis
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);



//&& Sum$(Jet_pt>30 && abs(Jet_eta)<2.5)<= (2 + %d)












    var      = "HCSV_pt";                // the variable to plot
    title    = ";Hpt [GeV]";          // the title of the histogram
    plotname = "All_hpt_tightcuts";      // the name of the image file
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);


    nbinsx       = 17;                      // number of bins
    xlow      = 0.0;                    // the low edge of x-axis
    xup       = 255.0;                   // the upper edge of x-axis 
    var      = "HCSV_mass";                // the variable to plot
    title    = ";Hmass [GeV]";          // the title of the histogram
    plotname = "All_hmass_tightcuts";      // the name of the image file
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);



    nbinsx       = 8;                      // number of bins
    xlow      = 0.0;                    // the low edge of x-axis
    xup       = 3.2;                   // the upper edge of x-axis 
    plotname = "All_dPhi_tightcuts_norm";      // the name of the image file
    var      = "abs(deltaPhi(HCSV_phi, met_phi))";                // the variable to plot
    title    = "#delta#Phi(MET,bb)";          // the title of the histogram
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);

    plotname = "All_dPhi_tightcuts";      // the name of the image file
    options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);


    nbinsx       = 6;                      // number of bins
    xlow      = 50.0;                    // the low edge of x-axis
    xup       = 350.;                   // the upper edge of x-axis 
    plotname = "All_Jet0_tightcuts_norm";      // the name of the image file
    var      = "max(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])";                // the variable to plot
    title    = "leading jet pt [GeV]";          // the title of the histogram
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);




    nbinsx       = 6;                      // number of bins
    xlow      = 50.0;                    // the low edge of x-axis
    xup       = 350;                   // the upper edge of x-axis 
    plotname = "All_Jet1_tightcuts_norm";      // the name of the image file
    var      = "min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])";                // the variable to plot
    title    = "trailing jet pt [GeV]";          // the title of the histogram
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);




    nbinsx       = 10;                      // number of bins
    xlow      = 50.0;                    // the low edge of x-axis
    xup       = 550;                   // the upper edge of x-axis 
    plotname = "All_Hpt_tightcuts_norm";      // the name of the image file
    var      = "HCSV_pt";                // the variable to plot
    title    = "higgs pt [GeV]";          // the title of the histogram
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);





    nbinsx       = 8;                      // number of bins
    xlow      = 0.0;                    // the low edge of x-axis
    xup       = 3.2;                   // the upper edge of x-axis 
    plotname = "All_dRbb_tightcuts_norm";      // the name of the image file
    var      = "deltaR_jj";                // the variable to plot
    title    = "#delta R(b,b)";          // the title of the histogram
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);



    // additionla jets    Jet_pt[aJidx[0]]
    nbinsx    = 7;                      // number of bins
    xlow      = 0.0;                    // the low edge of x-axis
    xup       = 210;                   // the upper edge of x-axis 
    plotname = "All_AddJetpt_tightcuts_norm";      // the name of the image file
    var      = "Jet_pt[aJidx[0]]";                // the variable to plot
    title    = "additional jet pt [GeV]";          // the title of the histogram
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);


    // additionla jets    Jet_pt[aJidx[0]]
    nbinsx       = 7;                      // number of bins
    xlow      = 0.0;                    // the low edge of x-axis
    xup       = 210;                   // the upper edge of x-axis 
    plotname = "All_AddJetpt_tightcuts";      // the name of the image file
    var      = "Jet_pt[aJidx[0]]";                // the variable to plot
    title    = "additional jet pt [GeV]";          // the title of the histogram
    options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);




    //    evalHMETMassiveMt(double m0, double pt0, double phi0, double m1, double pt1, double phi1)
    nbinsx       = 10;                      // number of bins                                                                                                                             
    xlow      = 50.0;                    // the low edge of x-axis                                                                                                                                  
    xup       = 1050;                   // the upper edge of x-axis                                                                                                                                 
    plotname = "All_MassiveMt_tightcuts_norm";      // the name of the image file                                                                                                                  
    var      = "evalHMETMassiveMt(1, met_pt, met_phi, 125, HCSV_pt,HCSV_phi) ";                // the variable to plot   
    title    = "MT(higgs+MET) [GeV]";          // the title of the histogram                                                                                                                     
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);


    //    evalHMETPt(double m0, double pt0, double phi0, double eta0,  double m1, double pt1, double phi1, double eta1)
    nbinsx       = 10;                      // number of bins                                                                                                                                       
    xlow      = 0.0;                    // the low edge of x-axis                                                                                                                                  
    xup       = 500;                   // the upper edge of x-axis                                                                                                                                 
    plotname = "All_MassivePt_tightcuts_norm";      // the name of the image file                                                                                                                  
    var      = "evalHMETPt(1, met_pt, met_phi, met_eta, HCSV_mass, HCSV_pt,HCSV_phi, HCSV_eta)";                // the variable to plot   
    title    = "PT(higgs+MET) [GeV]";          // the title of the histogram                                                                                                                     
    options  = "printStat:plotSig:!plotData:!plotLog:plotNorm";
    MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);

    

    /*

    ev->read(cutmc_all, cutdata_all);

    double vpt=150;
    double hpt=150;
    double mincsv=0.8;
    double jet1pt=80;
    double jet2pt=70;
    double   dPhiBB =1.3;
    int   najets =0;



    //	TCut  cut = "met_pt>150 && HCSV_pt>150 && min(hJets_btagCSV[0], hJets_btagCSV[1])>0.8 && Jet_pt[hJCidx[0]]>80 && Jet_pt[hJCidx[1]]>80  && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5)<=2";
		
		int  nbinsx       = 8;                      // number of bins
		double xlow      = 80.0;                    // the low edge of x-axis
		double xup       = 160.0;                   // the upper edge of x-axis 
		TString var      = "H_mass";                // the variable to plot
		TString title    = ";Hmass [GeV]";          // the title of the histogram
		TString plotname = "All_hmass_tightcuts";      // the name of the image file
		TString options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";
		MakePlots(ev, var, cutmc_all, cutdata_all , title, nbinsx, xlow, xup, plotname, plotdir, options);




     for (vpt=150; vpt<400; vpt+=50 ) {
       //   for (hpt=150 ; hpt<400 ; hpt+=50 ) {
       //	for (mincsv=0.8 ; mincsv<=1.0 ; mincsv+=0.1) {
	    for (jet2pt=80 ; jet2pt<120 ; jet2pt+=20) {
	      for (jet1pt=90 ; jet1pt<140 ; jet1pt+=10) {


         //	      for (dPhiBB=1.3 ; dPhiBB<1.5 ; dPhiBB+=0.2) {
	      for (najets=0 ; najets<=2 ; najets+=1) {

	TCut  cut = Form("met_pt>%.2f && H_pt>%.2f && min(hJets_btagCSV[0], hJets_btagCSV[1])>%.2f && Jet_pt[hJCidx[0]]>%.2f && Jet_pt[hJCidx[1]]>%.2f  && Sum$(Jet_pt>30 && abs(Jet_eta)<2.5)<= (2 + %d) && deltaR_jj<1.3", vpt , 150., 0.8, jet1pt, jet2pt , najets);
		
		 nbinsx       = 8;                      // number of bins
		 xlow      = 80.0;                    // the low edge of x-axis
		 xup       = 160.0;                   // the upper edge of x-axis 
		 var      = "H_mass";                // the variable to plot
		 title    = ";Hmass [GeV]";          // the title of the histogram
		 plotname = "All_hmass_tightcuts";      // the name of the image file
		 options  = "printStat:plotSig:!plotData:!plotLog:!plotNorm";
		MakePlots(ev, var, cut, cut , title, nbinsx, xlow, xup, plotname, plotdir, options);

	      }
	    }
	  }
	}
     //}
     // }

     */


    // You can put in the parameters directly as in the following commented out line:
    //MakePlot(ev->ZH, "H.pt", cut, "; p_{T}(jj) [GeV]", 16, 0, 240., process + "_Hpt", plotdir, options);




    ////////////////////////////////////////////////////////////////////////////
    // Task 1 (b)                                                             //
    // - please comment out other tasks                                       //
    ////////////////////////////////////////////////////////////////////////////
/*
    // Read from ntuples
    Events * ev = new Events();
    TCut cutmc_all   = "Vtype==4";   // change this to your channel
    TCut cutdata_all = "Vtype==4";   // change this to your channel
    ev->read(cutmc_all, cutdata_all, "VH:ZJ");  // read both VH and ZJ processes

    TString var      = "H.mass";
    //TCut cut         = "V.pt>120";
    TCut cut         = "V.pt>120 && hJet_pt[0]>30 && hJet_pt[1]>30";  // for Wln, Znn, change to tighter cut
    TString title    = ";m(jj) [GeV]";
    TString plotname = "ZH_vs_ZJ_Hmass";
    int nbinsx       = 15;
    double xlow      = 30.0;
    double xup       = 255.0;
    TString options  = "!!plotLog:plotNorm";

    // Using "ev->ZH" for ZH and "ev->ZjHF" for Z+HF
    MakePlot2(ev->ZH, ev->ZjHF, var, cut, title, nbinsx, xlow, xup, plotname, plotdir, options);
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
    //TString channel  = "Wen";

    // These are loose cuts for all plots in this particular channel
    //TCut cutmc_all   = "Vtype==3 && H.pt>80 && hJet_pt[0]>30 && hJet_pt[1]>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && vLepton_pt[0]>30 && abs(vLepton_eta[0])<2.5 && METtype1corr.et >45 && nalep==0 && Sum$(aJet_pt>20 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)==0 && hbhe==1";
    //cutmc_all       += "min(hJet_csv_nominal[0], hJet_csv_nominal[1])>0.4 && abs(HVdPhi)>2.0";  // tighter cut
    //TCut cutdata_all = cutmc_all;
    //cutmc_all       *= "weightTrig2012";  // apply trigger weight for MC

    // Scale factors in order of: WjLF, WjHF, ZjLF, ZjHF, TT
    // NOTE: ZjLF, ZjHF are not needed for Wln
    //double scalefactors[5] = {1.00, 1.00, 1.00, 1.00, 1.00};


    // Znn______________________________________________________________________
    TString channel  = "Znn";

    // These are loose cuts for all plots in this particular channel
    TCut cutmc_all   = "Vtype==4 && H.pt>130 && hJet_pt[0]>80 && hJet_pt[1]>30 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && nalep==0 && Sum$(aJet_pt>25 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)==0 && min(Min$(abs(deltaPhi(METtype1corr.phi,hJet_phi))),Min$(abs(deltaPhiMETjets(METtype1corr.phi,aJet_phi,aJet_pt,aJet_eta)))+999*(Sum$(aJet_pt>25 && abs(aJet_eta)<4.5 && aJet_id==1 && aJet_puJetIdL>0)==0) )>0.5 && hbhe==1";
    cutmc_all       += "min(hJet_csv_nominal[0], hJet_csv_nominal[1])>0.4 && abs(HVdPhi)>2.0";  // tighter cut
    TCut cutdata_all = cutmc_all;
    cutmc_all       *= "(triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1) && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent"; // apply trigger bits and MET cleaning for MC (they are already applied on data)
    cutmc_all       *= "triggercorrMET(METtype1corr.et)";

    // Scale factors in order of: WjLF, WjHF, ZjLF, ZjHF, TT
    double scalefactors[5] = {1.00, 1.00, 1.00, 1.00, 1.00};


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
    //double vpt    = 150.;
    //double hpt    = 100.;
    //double maxcsv = 0.898;
    //double mincsv = 0.5;
    //double dPhi   = 2.95;

    // Optimize these five variables (default: recommendations for Znn)
    double vpt    = 170.;  // for Znn, pT(V) = MET
    double hpt    = 170.;
    double maxcsv = 0.898;
    double mincsv = 0.5;
    double dPhi   = 2.95;

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
    TString options  = "printStat:plotSig:!plotData:!!plotLog";

    MakePlots(ev, var, cutmc, cutdata, title, nbinsx, xlow, xup, plotname, plotdir, options);

    // Or, just put in them directly as in the following commented out line:
    //MakePlots(ev, "H.mass", cutmc, cutdata, "m(jj) [GeV]", 15, 30.0, 255.0, channel+"_Hmass", plotdir, "printStat:plotSig:!plotData:!!plotLog");
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
