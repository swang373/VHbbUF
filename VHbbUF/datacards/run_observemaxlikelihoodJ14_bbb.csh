echo "*** ZnunuHighPt ***" >! M1.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm zhinv_Zbb_J14_bbb_ZnunuHighPt_8TeV.txt >> M1.log
cp mlfit.root mlfit_ZnunuHighPt.root
echo "*** ZnunuMedPt  ***" >! M2.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm zhinv_Zbb_J14_bbb_ZnunuMedPt_8TeV.txt >> M2.log
cp mlfit.root mlfit_ZnunuMedPt.root
echo "*** ZnunuLowPt  ***" >! M3.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm zhinv_Zbb_J14_bbb_ZnunuLowPt_8TeV.txt >> M3.log
cp mlfit.root mlfit_ZnunuLowPt.root
#echo "*** ZnunuLowCSV ***" >! M4.log
#combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm zhinv_Zbb_J14_bbb_ZnunuLowCSV_8TeV.txt >> M4.log
#cp mlfit.root mlfit_ZnunuLowCSV.root
echo "*** combo       ***" >! M5.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm zhinv_Zbb_J14_bbb_combo_8TeV.txt >> M5.log
cp mlfit.root mlfit_combo.root
echo "*** combo2      ***" >! M6.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm zhinv_Zbb_J14_bbb_combo2_8TeV.txt >> M6.log
cp mlfit.root mlfit_combo2.root

cat M1.log M2.log M3.log M5.log M6.log


echo "*** ZnunuHighPt ***" >! N1.html
python diffNuisances.py mlfit_ZnunuHighPt.root -f html >> N1.html
echo "*** ZnunuMedPt ***" >! N2.html
python diffNuisances.py mlfit_ZnunuMedPt.root -f html >> N2.html
echo "*** ZnunuLowPt ***" >! N3.html
python diffNuisances.py mlfit_ZnunuLowPt.root -f html >> N3.html
#echo "*** ZnunuLowCSV ***" >! N4.html
#python diffNuisances.py mlfit_ZnunuLowCSV.root -f html >> N4.html
echo "*** combo ***" >! N5.html
python diffNuisances.py mlfit_combo.root -f html >> N5.html
echo "*** combo2 ***" >! N6.html
python diffNuisances.py mlfit_combo2.root -f html >> N6.html

# To save post-fit shapes
#combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveShapes --saveWithUncertainties zhinv_Zbb_J14_bbb_combo_8TeV.txt

