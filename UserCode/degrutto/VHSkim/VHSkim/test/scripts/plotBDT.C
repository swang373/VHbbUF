void plotBDT()
{

bool doData = true;
float xmax = 0.5;
float xmin = -1.2;
float deltaBin = 0.025;


double lumi = 0.186*2;
if(!doData) lumi = 10;
const double lzh = 7638.0.;
 
const double lwh = 1321.; 

const double lzjmad = 0.74456; 
const double lwj = 0.463; 
const double lwz = 42. ;
const double lww = 48.;
const double lzz = 357.; 

const double ltt = 5.05; 

const double lts = 499.; 
const double lTt = 23.; 
const double ltw= 46.; 

double lz0         = 0.586085; 
double lz10         = 30.858 * (2.03/33.3);  
double lz11         = 1080.613 * (0.78 /33.3); 
double lz13         = 60588.240 * (0.90 / 33.3);
double lz18         = 2879753.581 * ( 2.15 / 33.3);
double lz20         = 10.115 * (6.35/33.3); 
double lz21         = 104.075 * (4.54/33.3);
double lz23         = 4482.161 * (6.35/33.3);
double lz28         = 162135.973 * (8.71 / 33.3);
double lz30         =   (8.50/33.3) *       11.084;   
double lz31         =   (11.0/33.3) *       58.528;   
double lz33         =   (10.2/33.3) *     1729.905;   
double lz38         =   (23.1/33.3) *    77584.126;   
double lz40         =   (31.9/33.3) *       18.665;   
double lz41         =   (17.4/33.3) *       73.680;   
double lz43         =   (9.76/33.3) *      915.589;   
double lz48         =   (29.6/33.3) *    73920.903;   
double lz50         =   (15.0/33.3) *       21.582;   
double lz51         =   (13.1/33.3) *       50.731;   
double lz53         =   (19.4/33.3) *     1006.313;   
double lz58         =   (13.2/33.3) *   346530.969;  
double lz0b= 119.400  ; // 204000. ;
double lz1b= 123.100 ; // 211000. ;
double lz2b= 17.400 ; // 30000.;
double lz3b= 39.900 ; // 30000.;
double lz0c= 150.100 ;//* 4./ 5. ; //256000. * 4./5.;
double lz1c= 113.200 ;
double lz2c= 17.160 ;
double lz3c= 36.200 ;
float ll = lz0 +lz10+lz11+lz13+lz18+lz20+lz21+lz23+lz28+lz30+lz31+lz33+lz38+lz40+lz41+lz43+lz48+lz50+lz51+lz53+lz58;
float lh = lz0b+lz1b+lz2b+lz3b+lz0c+lz1c+lz2c+lz3c;

const double eelz0 = 0.586085;  //49530.  * 1.3;
const double eelz10 = 30.858 * (95.9/33.3);  //49530.  * 1.3;
const double eelz11 = 1080.613 * (98.4 /33.3); //1375000. ;
const double eelz13 = 60588.240 * (98.2 / 33.3); // 76945000.;
const double eelz18 = 2879753.581 * (95.7 / 33.3);// 3660000000. ;
const double eelz20 = 10.115 * (87.3/33.3); //12850. ;
const double eelz21 = 104.075 * (90.9/33.3); //132000. ;
const double eelz23 = 4482.161 * (87.2/33.3); //5693000. ;
const double eelz28 = 162135.973 * (82.8 / 33.3); //05291000. ;
const double eelz30  =   (83.2/33.33) *     11.084;
const double eelz31  =   (77.9/33.33) *     58.528;
const double eelz33  =   (79.6/33.33) *   1729.905;
const double eelz38  =   (54.2/33.33) *  77584.126;
const double eelz40  =   (36.3/33.33) *     18.665;
const double eelz41  =   (65.2/33.33) *     73.680;
const double eelz43  =   (80.3/33.33) *    915.589;
const double eelz48  =   (41.3/33.33) *  73920.903;
const double eelz50  =   (69.5/33.33) *     21.582;
const double eelz51  =   (74.4/33.33) *     50.731;
const double eelz53  =   (61.3/33.33) *   1006.313;
const double eelz58  =   (73.4/33.33) * 346530.969;
const double eelz0b= 119.400.  ; // 204000. ;
const double eelz1b= 123.100. ; // 211000. ;
const double eelz2b= 17.400. ; // 30000.;
const double eelz3b= 39.900. ; // 30000.;
const double eelz0c= 150.100. ;//* 4./ 5. ; //256000. * 4./5.;
const double eelz1c= 113.200. ;
const double eelz2c= 17.160. ;
const double eelz3c= 36.200. ;



TString *title  = new TString ("TMVA Output: BDT");


gROOT->SetStyle("Plain");
const Color_t hLineColor = kBlack;
const Color_t hFillColor = 0;
const Color_t vjLineColor = kAzure+2;
const Color_t vjFillColor = kAzure+6;
const Color_t vvLineColor = kViolet+3;
const Color_t vvFillColor = kViolet-5;
const Color_t ttLineColor =  kRed+3;
const Color_t ttFillColor = kRed+1;




TFile *fTT = new TFile("TMVA-TT.root","read");
TFile *fZZ = new TFile("TMVA-ZZ.root","read");
//TFile *fZJ = new TFile("TMVA-ZJ.root","read");
TFile *fZH = new TFile("TMVA-ZH.root","read");
TFile *fDT = new TFile("TMVA-DATA.root","read");
TFile *fZ0 = new TFile("TMVA-Z0Jets.root");
TFile *f10= new TFile("TMVA-Z1Jets0-100.root");
TFile *f11 = new TFile("TMVA-Z1Jets100-300.root");
TFile *f13 = new TFile("TMVA-Z1Jets300-800.root");
TFile *f18 = new TFile("TMVA-Z1Jets800-1600.root");
TFile *f20 = new TFile("TMVA-Z2Jets0-100.root");
TFile *f21 = new TFile("TMVA-Z2Jets100-300.root");
TFile *f23 = new TFile("TMVA-Z2Jets300-800.root");
TFile *f28 = new TFile("TMVA-Z2Jets800-1600.root");
TFile *f30 = new TFile("TMVA-Z3Jets0-100.root");
TFile *f31 = new TFile("TMVA-Z3Jets100-300.root");
TFile *f33 = new TFile("TMVA-Z3Jets300-800.root");
TFile *f38 = new TFile("TMVA-Z3Jets800-1600.root");
TFile *f40 = new TFile("TMVA-Z4Jets0-100.root");
TFile *f41 = new TFile("TMVA-Z4Jets100-300.root");
TFile *f43 = new TFile("TMVA-Z4Jets300-800.root");
TFile *f48 = new TFile("TMVA-Z4Jets800-1600.root");
TFile *f50 = new TFile("TMVA-Z5Jets0-100.root");
TFile *f51 = new TFile("TMVA-Z5Jets100-300.root");
TFile *f53 = new TFile("TMVA-Z5Jets300-800.root");
TFile *f58 = new TFile("TMVA-Z5Jets800-1600.root");
TFile *fb0 = new TFile("TMVA-Zbb0Jets.root");;
TFile *fb1 = new TFile("TMVA-Zbb1Jets.root");
TFile *fb2 = new TFile("TMVA-Zbb2Jets.root");
TFile *fb3 = new TFile("TMVA-Zbb3Jets.root");
TFile *fc0 = new TFile("TMVA-Zcc0Jets.root");
TFile *fc1 = new TFile("TMVA-Zcc1Jets.root");
TFile *fc2 = new TFile("TMVA-Zcc2Jets.root");
TFile *fc3 = new TFile("TMVA-Zcc3Jets.root");
TFile *fts = new TFile("TMVA-Ts.root");
TFile *fTt = new TFile("TMVA-Tt.root");
TFile *ftw = new TFile("TMVA-TtW.root");
TFile *fww = new TFile("TMVA-WW.root");
TFile *fwz = new TFile("TMVA-WZ.root");
TFile *fwj = new TFile("TMVA-WJ.root");




TTree *tTT = (TTree*) fTT->Get("tree");
TTree *tZH = (TTree*) fZH->Get("tree");
TTree *tZZ = (TTree*) fZZ->Get("tree");
TTree *tDT = (TTree*) fDT->Get("tree");
 TTree *tZ0 = (TTree*)  fZ0->Get("tree");
 TTree *t10 = (TTree*)  f10->Get("tree");
 TTree *t11 = (TTree*)  f11->Get("tree");
 TTree *t13 = (TTree*)  f13->Get("tree");
 TTree *t18 = (TTree*)  f18->Get("tree");
 TTree *t20 = (TTree*)  f20->Get("tree");
 TTree *t21 = (TTree*)  f21->Get("tree");
 TTree *t23 = (TTree*)  f23->Get("tree");
 TTree *t28 = (TTree*)  f28->Get("tree");
 TTree *t30 = (TTree*)  f30->Get("tree");
 TTree *t31 = (TTree*)  f31->Get("tree");
 TTree *t33 = (TTree*)  f33->Get("tree");
 TTree *t38 = (TTree*)  f38->Get("tree");
 TTree *t40 = (TTree*)  f40->Get("tree");
 TTree *t41 = (TTree*)  f41->Get("tree");
 TTree *t43 = (TTree*)  f43->Get("tree");
 TTree *t48 = (TTree*)  f48->Get("tree");
 TTree *t50 = (TTree*)  f50->Get("tree");
 TTree *t51 = (TTree*)  f51->Get("tree");
 TTree *t53 = (TTree*)  f53->Get("tree");
 TTree *t58 = (TTree*)  f58->Get("tree");
 TTree *tb0 = (TTree*)  fb0->Get("tree");
 TTree *tb1 = (TTree*)  fb1->Get("tree");
 TTree *tb2 = (TTree*)  fb2->Get("tree");
 TTree *tb3 = (TTree*)  fb3->Get("tree");
 TTree *tc0 = (TTree*)  fc0->Get("tree");
 TTree *tc1 = (TTree*)  fc1->Get("tree");
 TTree *tc2 = (TTree*)  fc2->Get("tree");
 TTree *tc3 = (TTree*)  fc3->Get("tree");
 TTree *tts = (TTree*)  fts->Get("tree");
 TTree *tTt = (TTree*)  fTt->Get("tree");
 TTree *ttw = (TTree*)  ftw->Get("tree");
 TTree *tww = (TTree*)  fww->Get("tree");
 TTree *twz = (TTree*)  fwz->Get("tree");
 TTree *twj = (TTree*)  fwj->Get("tree");






TH1F *hTT = new TH1F("hTT","hTT",(xmax-xmin)/deltaBin,xmin,xmax);
hTT->SetLineColor(ttLineColor);
hTT->SetFillColor(ttFillColor);

TH1F *hZZ = new TH1F("hZZ","hZZ",(xmax-xmin)/deltaBin,xmin,xmax);
hZZ->SetLineColor(vvLineColor);
hZZ->SetFillColor(vvFillColor);


TH1F *lZJ = new TH1F("lZJ","lZJ",(xmax-xmin)/deltaBin,xmin,xmax);
lZJ->SetLineColor(vjLineColor);
lZJ->SetFillColor(vjFillColor);

TH1F *hZJ = new TH1F("hZJ","hZJ",(xmax-xmin)/deltaBin,xmin,xmax);
hZJ->SetLineColor(kGreen + 3);
hZJ->SetFillColor(kGreen + 4);



TH1F *hDT = new TH1F("hDT","hDT",(xmax-xmin)/deltaBin,xmin,xmax);
hDT->SetLineColor(ttLineColor);
hDT->SetFillColor(ttFillColor);

TH1F *hZH = new TH1F("hZH","hZH",(xmax-xmin)/deltaBin,xmin,xmax);
hZH->SetLineColor(kBlack);
hZH->SetLineWidth(3);
hZH->SetFillColor(0);

TH1F *h0 = new TH1F("h0","h0",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h10 = new TH1F("h10","h10",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h11 = new TH1F("h11","h11",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h13 = new TH1F("h13","h13",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h18 = new TH1F("h18","h18",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h20 = new TH1F("h20","h20",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h21 = new TH1F("h21","h21",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h23 = new TH1F("h23","h23",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h28 = new TH1F("h28","h28",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h30 = new TH1F("h30","h30",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h31 = new TH1F("h31","h31",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h33 = new TH1F("h33","h33",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h38 = new TH1F("h38","h38",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h40 = new TH1F("h40","h40",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h41 = new TH1F("h41","h41",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h43 = new TH1F("h43","h43",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h48 = new TH1F("h48","h48",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h50 = new TH1F("h50","h50",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h51 = new TH1F("h51","h51",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h53 = new TH1F("h53","h53",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *h58 = new TH1F("h58","h58",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hb0 = new TH1F("hb0","hb0",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hb1 = new TH1F("hb1","hb1",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hb2 = new TH1F("hb2","hb2",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hb3 = new TH1F("hb3","hb3",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hc0 = new TH1F("hc0","hc0",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hc1 = new TH1F("hc1","hc1",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hc2 = new TH1F("hc2","hc2",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hc3 = new TH1F("hc3","hc3",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hts = new TH1F("hts","hts",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hTt = new TH1F("hTt","hTt",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *htw = new TH1F("htw","htw",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hww = new TH1F("hww","hww",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hwz = new TH1F("hwz","hwz",(xmax-xmin)/deltaBin,xmin,xmax);
TH1F *hwj = new TH1F("hwj","hwj",(xmax-xmin)/deltaBin,xmin,xmax);

  TH1F *l0  = new TH1F ("l0", "l0",   (xmax-xmin)/deltaBin, xmin, xmax);			
  TH1F *l10 = new TH1F ("l10","l10", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l11 = new TH1F ("l11","l11", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l13 = new TH1F ("l13","l13", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l18 = new TH1F ("l18","l18", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l20 = new TH1F ("l20","l20", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l21 = new TH1F ("l21","l21", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l23 = new TH1F ("l23","l23", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l28 = new TH1F ("l28","l28", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l30 = new TH1F ("l30","l30", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l31 = new TH1F ("l31","l31", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l33 = new TH1F ("l33","l33", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l38 = new TH1F ("l38","l38", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l40 = new TH1F ("l40","l40", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l41 = new TH1F ("l41","l41", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l43 = new TH1F ("l43","l43", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l48 = new TH1F ("l48","l48", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l50 = new TH1F ("l50","l50", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l51 = new TH1F ("l51","l51", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l53 = new TH1F ("l53","l53", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l58 = new TH1F ("l58","l58", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l0b = new TH1F ("lb0","l0b", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l1b = new TH1F ("lb1","l1b", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l2b = new TH1F ("lb2","l2b", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l3b = new TH1F ("lb3","l3b", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l0c = new TH1F ("lc0","l0c", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l1c = new TH1F ("lc1","l1c", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l2c = new TH1F ("lc2","l2c", (xmax-xmin)/deltaBin, xmin, xmax);
  TH1F *l3c = new TH1F ("lc3","l3c", (xmax-xmin)/deltaBin, xmin, xmax);




tTT->Project("hTT","BDT");
tZZ->Project("hZZ","BDT");

tZ0->Project("lZ0","BDT","heavy < 0");
t10->Project("l10","BDT","heavy < 0");
t11->Project("l11","BDT","heavy < 0");
t13->Project("l13","BDT","heavy < 0");
t18->Project("l18","BDT","heavy < 0");
t20->Project("l20","BDT","heavy < 0");
t21->Project("l21","BDT","heavy < 0");
t23->Project("l23","BDT","heavy < 0");
t28->Project("l28","BDT","heavy < 0");
t30->Project("l30","BDT","heavy < 0");
t31->Project("l31","BDT","heavy < 0");
t33->Project("l33","BDT","heavy < 0");
t38->Project("l38","BDT","heavy < 0");
t40->Project("l40","BDT","heavy < 0");
t41->Project("l41","BDT","heavy < 0");
t43->Project("l43","BDT","heavy < 0");
t48->Project("l48","BDT","heavy < 0");
t50->Project("l50","BDT","heavy < 0");
t51->Project("l51","BDT","heavy < 0");
t53->Project("l53","BDT","heavy < 0");
t58->Project("l58","BDT","heavy < 0");
tb0->Project("lb0","BDT","heavy < 0");
tb1->Project("lb1","BDT","heavy < 0");
tb2->Project("lb2","BDT","heavy < 0");
tb3->Project("lb3","BDT","heavy < 0");
tc0->Project("lc0","BDT","heavy < 0");
tc1->Project("lc1","BDT","heavy < 0");
tc2->Project("lc2","BDT","heavy < 0");
tc3->Project("lc3","BDT","heavy < 0");

tZ0->Project("hZ0","BDT","heavy > 0");
t10->Project("h10","BDT","heavy > 0");
t11->Project("h11","BDT","heavy > 0");
t13->Project("h13","BDT","heavy > 0");
t18->Project("h18","BDT","heavy > 0");
t20->Project("h20","BDT","heavy > 0");
t21->Project("h21","BDT","heavy > 0");
t23->Project("h23","BDT","heavy > 0");
t28->Project("h28","BDT","heavy > 0");
t30->Project("h30","BDT","heavy > 0");
t31->Project("h31","BDT","heavy > 0");
t33->Project("h33","BDT","heavy > 0");
t38->Project("h38","BDT","heavy > 0");
t40->Project("h40","BDT","heavy > 0");
t41->Project("h41","BDT","heavy > 0");
t43->Project("h43","BDT","heavy > 0");
t48->Project("h48","BDT","heavy > 0");
t50->Project("h50","BDT","heavy > 0");
t51->Project("h51","BDT","heavy > 0");
t53->Project("h53","BDT","heavy > 0");
t58->Project("h58","BDT","heavy > 0");
tb0->Project("hb0","BDT","heavy > 0");
tb1->Project("hb1","BDT","heavy > 0");
tb2->Project("hb2","BDT","heavy > 0");
tb3->Project("hb3","BDT","heavy > 0");
tc0->Project("hc0","BDT","heavy > 0");
tc1->Project("hc1","BDT","heavy > 0");
tc2->Project("hc2","BDT","heavy > 0");
tc3->Project("hc3","BDT","heavy > 0");


tts->Project("hts","BDT");
tTt->Project("hTt","BDT");
ttw->Project("htw","BDT");
tww->Project("hww","BDT");
twz->Project("hwz","BDT");
twj->Project("hwj","BDT");
if(doData) tDT->Project("hDT","BDT");
else tZH->Project("hZH","BDT");

h0->Scale(lumi/lz0);
h10->Scale(lumi/lz10);
h11->Scale(lumi/lz11);
h13->Scale(lumi/lz13);
h18->Scale(lumi/lz18);
h20->Scale(lumi/lz20);
h21->Scale(lumi/lz21);
h23->Scale(lumi/lz23);
h28->Scale(lumi/lz28);
h30->Scale(lumi/lz30);
h31->Scale(lumi/lz31);
h33->Scale(lumi/lz33);
h38->Scale(lumi/lz38);
h40->Scale(lumi/lz40);
h41->Scale(lumi/lz41);
h43->Scale(lumi/lz43);
h48->Scale(lumi/lz48);
h50->Scale(lumi/lz50);
h51->Scale(lumi/lz51);
h53->Scale(lumi/lz53);
h58->Scale(lumi/lz58);
hb0->Scale(lumi/lz0b);
hb1->Scale(lumi/lz1b);
hb2->Scale(lumi/lz2b);
hb3->Scale(lumi/lz3b);
hc0->Scale(lumi/lz0c);
hc1->Scale(lumi/lz1c);
hc2->Scale(lumi/lz2c); 
hc3->Scale(lumi/lz3c); 

l0->Scale(lumi/lz0);
l10->Scale(lumi/lz10);
l11->Scale(lumi/lz11);
l13->Scale(lumi/lz13);
l18->Scale(lumi/lz18);
l20->Scale(lumi/lz20);
l21->Scale(lumi/lz21);
l23->Scale(lumi/lz23);
l28->Scale(lumi/lz28);
l30->Scale(lumi/lz30);
l31->Scale(lumi/lz31);
l33->Scale(lumi/lz33);
l38->Scale(lumi/lz38);
l40->Scale(lumi/lz40);
l41->Scale(lumi/lz41);
l43->Scale(lumi/lz43);
l48->Scale(lumi/lz48);
l50->Scale(lumi/lz50);
l51->Scale(lumi/lz51);
l53->Scale(lumi/lz53);
l58->Scale(lumi/lz58);
lb0->Scale(lumi/lz0b);
lb1->Scale(lumi/lz1b);
lb2->Scale(lumi/lz2b);
lb3->Scale(lumi/lz3b);
lc0->Scale(lumi/lz0c);
lc1->Scale(lumi/lz1c);
lc2->Scale(lumi/lz2c); 
lc3->Scale(lumi/lz3c); 


hZZ->Scale(lumi/lzz);



hTT->Scale(lumi/ltt);

hts->Scale(lumi/lts);
hTt->Scale(lumi/lTt);
htw->Scale(lumi/ltw);
hww->Scale(lumi/lww);
hwz->Scale(lumi/lwz);
hwj->Scale(lumi/lwj);
hZH ->Scale(lumi/lzh);

lZJ->Add(l0);  lZJ->Add(l10); lZJ->Add(l11); lZJ->Add(l13); lZJ->Add(l18); lZJ->Add(l20); lZJ->Add(l21); lZJ->Add(l23); lZJ->Add(l28); lZJ->Add(l30); lZJ->Add(l31); lZJ->Add(l33); lZJ->Add(l38); lZJ->Add(l40); lZJ->Add(l41); lZJ->Add(l43); lZJ->Add(l48); lZJ->Add(l50); lZJ->Add(l51); lZJ->Add(l53); lZJ->Add(l58); lZJ->Add(lb0); lZJ->Add(lb1); lZJ->Add(lb2); lZJ->Add(lb3); lZJ->Add(lc0); lZJ->Add(lc1); lZJ->Add(lc2); lZJ->Add(lc3); 

hZJ->Add(h0);  hZJ->Add(h10); hZJ->Add(h11); hZJ->Add(h13); hZJ->Add(h18); hZJ->Add(h20); hZJ->Add(h21); hZJ->Add(h23); hZJ->Add(h28); hZJ->Add(h30); hZJ->Add(h31); hZJ->Add(h33); hZJ->Add(h38); hZJ->Add(h40); hZJ->Add(h41); hZJ->Add(h43); hZJ->Add(h48); hZJ->Add(h50); hZJ->Add(h51); hZJ->Add(h53); hZJ->Add(h58); hZJ->Add(hb0); hZJ->Add(hb1); hZJ->Add(hb2); hZJ->Add(hb3); hZJ->Add(hc0); hZJ->Add(hc1); hZJ->Add(hc2); hZJ->Add(hc3);   hZJ->Add(hwj);

hZZ->Add(hww);hZZ->Add(hwz);
hTT->Add(hTt);hTT->Add(htw);hTT->Add(hts);


   if(doData) {
    hDT->SetMarkerStyle(kFullCircle);
    hDT->SetMarkerSize(1);
    hDT->SetMarkerColor(kBlack);
    hDT->SetLineWidth(1);
    hDT->SetLineColor(kBlack);
    gStyle->SetEndErrorSize(1);
    hDT->GetXaxis()->SetLabelSize(0);
    hDT->GetYaxis()->SetLabelSize(0);
   }  

   THStack * hs = new THStack("hs","");

  float iZH = hZH->Integral();
  float eZH = sqrt(lumi/lzh*iZH);
  float iTT = hTT->Integral();
  float eTT = sqrt(lumi/ltt*iTT);
  float iZJl = lZJ->Integral();
  float eZJl = sqrt(lumi/ll*iZJl);
  float iZJh = hZJ->Integral();
  float eZJh = sqrt(lumi/lh*iZJh);
  float iZZ = hZZ->Integral();
  float eZZ = sqrt(lumi/lzz*iZZ);
  float iDT = hDT->Integral();
  float eDT = sqrt(iDT);
  float totMC =  iZZ+iTT+iZJl+iZJh+iZZ;
//  float totMCErr =   ((eTT*eTT)/(iTT*iTT) +  (eZZ*eZZ)/(iZZ*iZZ) + (eZJl*eZJl)/(iZJl*iZJl) + (eZJh*eZJh)/(iZJh*iZJh) ));

  float ratio = iDT/(iTT+iZJh+iZJl+iZZ);
    std::cout.setf(0,ios::floatfield);
   std::cout.setf(ios::fixed,ios::floatfield);
std::cout.precision(1);
std::cout << "ZH = " <<  iZH << " +/- " << eZH  << " </br> "<<std::endl;
std::cout << "ZJ lf =" <<  iZJl<< " +/- "  << eZJl << " </br> "<<std::endl;
std::cout << "ZJ hf = " <<  iZJh<< " +/- " << eZJh << " </br> "<<std::endl;
std::cout << "ZZ = " <<  iZZ<<  " +/- " << eZZ <<" </br> "<<std::endl;
std::cout << "TT = " <<  iTT<< " +/- " << eTT << " </br> "<<std::endl;
std::cout << "Total MC = " <<  totMC  << " +/- " << eZJl + eZJh + eZZ + eTT << " </br> "<<std::endl;

std::cout << "Data = " <<  iDT << " </br> "<<std::endl;
std::cout << "Data/MC = " << iDT/totMC << " </br> "<<std::endl;
//  hTT->Scale(ratio);
//  hZJ->Scale(ratio);
//  hZZ->Scale(ratio);
 

   hs->Add(hTT);
   hs->Add(hZZ);
   hs->Add(lZJ);
   hs->Add(hZJ);

if(doData) title->Append("            #font[12]{L} = 0.186 fb^{-1}");
else title->Append("             #font[12]{L} = 10 fb^{-1}");

TCanvas *c1= new TCanvas();


if(doData && hDT->GetMaximum() > hs->GetMaximum()) hs->SetMaximum(hDT->GetMaximum());
else if(doData) hDT->Draw("PE1SAME"); 

hs->Draw("HIST");
if(doData) hDT->Draw("PE1SAME"); 
hs->SetMinimum(0.1);
hs->GetXaxis()->SetTitle("BDT");
if(!doData) hZH->Draw("SAME");

c1->SetLogy();
TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.04);
latex.SetTextAlign(31); // align right
latex.DrawLatex(0.55,0.92,title->Data());


legend=new TLegend(0.69,0.72,0.85,0.85);
legend->SetLineColor(0);  legend->SetFillColor(0);
if (doData) legend->AddEntry("hDT","data", "pl");
if (!doData) legend->AddEntry("hZH","ZH","l");
legend->AddEntry(hZZ,"ZZ+WW+WZ","f");
legend->AddEntry(hZJ,"Z + Jets hf","f");
legend->AddEntry(lZJ,"Z + Jets lf","f");
legend->AddEntry(hTT,"t#bar{t} + st","f"); 
legend->SetShadowColor(kWhite);
legend->Draw();


float maxSig=0; float tempSigU,tempSigD;
for(int i=1;i<hZH->GetNbinsX();i++)
{
   tempSigU = hZH->Integral(i,9999) / (1.5 + sqrt(hZZ->Integral(i,9999) + hTT->Integral(i,9999) + hZJ->Integral(i,9999)));

   tempSigD = hZH->Integral(0,i) / (1.5 + sqrt(hZZ->Integral(0,i) + hTT->Integral(0,i) + hZJ->Integral(0,i)));
   
   if(tempSigU > maxSig) {maxSig = tempSigU; GT = true; goodbin = i;}
   if(tempSigD > maxSig) {maxSig = tempSigD; GT = false;goodbin = i;}

}

   std::cout.setf(0,ios::floatfield);
   std::cout.setf(ios::fixed,ios::floatfield);
   std::cout.precision(3);
 
  std::cout << "Cut At " << hZH ->GetBinCenter(goodbin) << " GT: " << GT << std::endl;
  std::cout << "Sig: " << maxSig << std::endl;

  TLine *ass = new TLine(hZH ->GetBinCenter(goodbin),-0.2  ,hZH ->GetBinCenter(goodbin),10000);
  ass->SetLineWidth(3);
  ass->SetLineColor(kGreen);
 /*here*/ if(doData) {c1->SaveAs("BDT-CCC.png");c1->SaveAs("BDT-CCC.eps");}

  else  c1->SaveAs("BDTMC.png");
  ass->Draw();
}



