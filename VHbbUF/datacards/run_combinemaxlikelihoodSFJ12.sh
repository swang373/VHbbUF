DIR=mlfit_SF_J12
mkdir -p $DIR

combineCards.py ZnunuHighPt_WjLF=vhbb_Znn_SF_J12_ZnunuHighPt_WjLF_8TeV.txt ZnunuHighPt_WjHF=vhbb_Znn_SF_J12_ZnunuHighPt_WjHF_8TeV.txt ZnunuHighPt_ZjLF=vhbb_Znn_SF_J12_ZnunuHighPt_ZjLF_8TeV.txt ZnunuHighPt_ZjHF=vhbb_Znn_SF_J12_ZnunuHighPt_ZjHF_8TeV.txt ZnunuHighPt_TT=vhbb_Znn_SF_J12_ZnunuHighPt_TT_8TeV.txt > vhbb_Znn_SF_J12_ZnunuHighPt_8TeV.txt
combineCards.py ZnunuMedPt_WjLF=vhbb_Znn_SF_J12_ZnunuMedPt_WjLF_8TeV.txt ZnunuMedPt_WjHF=vhbb_Znn_SF_J12_ZnunuMedPt_WjHF_8TeV.txt ZnunuMedPt_ZjLF=vhbb_Znn_SF_J12_ZnunuMedPt_ZjLF_8TeV.txt ZnunuMedPt_ZjHF=vhbb_Znn_SF_J12_ZnunuMedPt_ZjHF_8TeV.txt ZnunuMedPt_TT=vhbb_Znn_SF_J12_ZnunuMedPt_TT_8TeV.txt > vhbb_Znn_SF_J12_ZnunuMedPt_8TeV.txt
combineCards.py ZnunuLowPt_WjLF=vhbb_Znn_SF_J12_ZnunuLowPt_WjLF_8TeV.txt ZnunuLowPt_WjHF=vhbb_Znn_SF_J12_ZnunuLowPt_WjHF_8TeV.txt ZnunuLowPt_ZjLF=vhbb_Znn_SF_J12_ZnunuLowPt_ZjLF_8TeV.txt ZnunuLowPt_ZjHF=vhbb_Znn_SF_J12_ZnunuLowPt_ZjHF_8TeV.txt ZnunuLowPt_TT=vhbb_Znn_SF_J12_ZnunuLowPt_TT_8TeV.txt > vhbb_Znn_SF_J12_ZnunuLowPt_8TeV.txt
#combineCards.py ZnunuHighPt_WjLF=vhbb_Znn_SF_J12_ZnunuHighPt_WjLF_8TeV.txt ZnunuHighPt_WjHF=vhbb_Znn_SF_J12_ZnunuHighPt_WjHF_8TeV.txt ZnunuHighPt_ZjLF=vhbb_Znn_SF_J12_ZnunuHighPt_ZjLF_8TeV.txt ZnunuHighPt_ZjHF=vhbb_Znn_SF_J12_ZnunuHighPt_ZjHF_8TeV.txt ZnunuHighPt_TT=vhbb_Znn_SF_J12_ZnunuHighPt_TT_8TeV.txt ZnunuLowCSV_WjHF=vhbb_Znn_SF_J12_ZnunuLowCSV_WjHF_8TeV.txt ZnunuLowCSV_ZjHF=vhbb_Znn_SF_J12_ZnunuLowCSV_ZjHF_8TeV.txt > vhbb_Znn_SF_J12_ZnunuLowCSV_8TeV.txt

combineCards.py ZnunuHighPt_WjLF=vhbb_Znn_SF_J12_ZnunuHighPt_WjLF_8TeV.txt ZnunuHighPt_WjHF=vhbb_Znn_SF_J12_ZnunuHighPt_WjHF_8TeV.txt ZnunuHighPt_TT=vhbb_Znn_SF_J12_ZnunuHighPt_TT_8TeV.txt > vhbb_Znn_SF_J12_ZnunuLepHighPt_8TeV.txt
combineCards.py ZnunuMedPt_WjLF=vhbb_Znn_SF_J12_ZnunuMedPt_WjLF_8TeV.txt ZnunuMedPt_WjHF=vhbb_Znn_SF_J12_ZnunuMedPt_WjHF_8TeV.txt ZnunuMedPt_TT=vhbb_Znn_SF_J12_ZnunuMedPt_TT_8TeV.txt > vhbb_Znn_SF_J12_ZnunuLepMedPt_8TeV.txt
combineCards.py ZnunuLowPt_WjLF=vhbb_Znn_SF_J12_ZnunuLowPt_WjLF_8TeV.txt ZnunuLowPt_WjHF=vhbb_Znn_SF_J12_ZnunuLowPt_WjHF_8TeV.txt ZnunuLowPt_TT=vhbb_Znn_SF_J12_ZnunuLowPt_TT_8TeV.txt > vhbb_Znn_SF_J12_ZnunuLepLowPt_8TeV.txt

sed -i 's/kmax .*/kmax */' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
sed -i 's/CMS_vhbb_eff_b/#CMS_vhbb_eff_b/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
sed -i 's/CMS_vhbb_fake_b/#CMS_vhbb_fake_b/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
sed -i 's/CMS_vhbb_res_j/#CMS_vhbb_res_j/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
sed -i 's/CMS_vhbb_scale_j/#CMS_vhbb_scale_j/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
sed -i 's/CMS_vhbb_trigger_MET_Znunu_8TeV/#CMS_vhbb_trigger_MET_Znunu_8TeV/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
sed -i 's/CMS_vhbb_WJSlope/#CMS_vhbb_WJSlope/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
sed -i 's/CMS_vhbb_ZJSlope/#CMS_vhbb_ZJSlope/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/UEPS/#UEPS/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/CMS_vhbb_WJModel/#CMS_vhbb_WJModel/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/CMS_vhbb_ZJModel/#CMS_vhbb_ZJModel/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/CMS_vhbb_TTModel/#CMS_vhbb_TTModel/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt

#sed -i 's/kmax 24/kmax */' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/kmax 22/kmax */' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_eff_b/CMS_vhbb_eff_b/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_fake_b/CMS_vhbb_fake_b/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_res_j/CMS_vhbb_res_j/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_scale_j/CMS_vhbb_scale_j/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_trigger_MET_Znunu_8TeV/CMS_vhbb_trigger_MET_Znunu_8TeV/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_WJSlope/CMS_vhbb_WJSlope/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_ZJSlope/CMS_vhbb_ZJSlope/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_WJModel/CMS_vhbb_WJModel/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_ZJModel/CMS_vhbb_ZJModel/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt
#sed -i 's/#CMS_vhbb_TTModel/CMS_vhbb_TTModel/' vhbb_Znn_SF_J12_Znunu*Pt_8TeV.txt

combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J12_ZnunuHighPt_8TeV.txt
mv mlfit.root $DIR/mlfit_ZnunuHighPt.root
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J12_ZnunuLepHighPt_8TeV.txt
mv mlfit.root $DIR/mlfit_ZnunuLepHighPt.root

combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J12_ZnunuMedPt_8TeV.txt
mv mlfit.root $DIR/mlfit_ZnunuMedPt.root
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J12_ZnunuLepMedPt_8TeV.txt
mv mlfit.root $DIR/mlfit_ZnunuLepMedPt.root

combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J12_ZnunuLowPt_8TeV.txt
mv mlfit.root $DIR/mlfit_ZnunuLowPt.root
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J12_ZnunuLepLowPt_8TeV.txt
mv mlfit.root $DIR/mlfit_ZnunuLepLowPt.root

#combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J12_ZnunuLowCSV_8TeV.txt
#mv mlfit.root $DIR/mlfit_ZnunuLowCSV.root

#python diffNuisances.py -f html -a  $DIR/mlfit_ZnunuHighPt.root
#python diffNuisances.py -f html -a  $DIR/mlfit_ZnunuMedPt.root
#python diffNuisances.py -f html -a  $DIR/mlfit_ZnunuLowPt.root

echo "*** ZnunuHighPt ***"
python printNuisances.py $DIR/mlfit_ZnunuHighPt.root
echo "*** ZnunuMedPt ***"
python printNuisances.py $DIR/mlfit_ZnunuMedPt.root
echo "*** ZnunuLowPt ***"
python printNuisances.py $DIR/mlfit_ZnunuLowPt.root

echo "*** ZnunuLepHighPt ***"
python printNuisances.py $DIR/mlfit_ZnunuLepHighPt.root
echo "*** ZnunuLepMedPt ***"
python printNuisances.py $DIR/mlfit_ZnunuLepMedPt.root
echo "*** ZnunuLepLowPt ***"
python printNuisances.py $DIR/mlfit_ZnunuLepLowPt.root

