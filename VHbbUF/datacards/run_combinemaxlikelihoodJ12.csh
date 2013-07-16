echo "*** ZnunuHighPt ***" >! m1.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuHighPt_8TeV.txt >> m1.log
cp mlfit.root mlfit_toy_ZnunuHighPt.root
echo "*** ZnunuMedPt  ***" >! m2.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuMedPt_8TeV.txt >> m2.log
cp mlfit.root mlfit_toy_ZnunuMedPt.root
echo "*** ZnunuLowPt  ***" >! m3.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuLowPt_8TeV.txt >> m3.log
cp mlfit.root mlfit_toy_ZnunuLowPt.root
#echo "*** ZnunuLowCSV ***" >! m4.log
#combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_ZnunuLowCSV_8TeV.txt >> m4.log
#cp mlfit.root mlfit_toy_ZnunuLowCSV.root
echo "*** combo       ***" >! m5.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_combo_8TeV.txt >> m5.log
cp mlfit.root mlfit_toy_combo.root
echo "*** combo2      ***" >! m6.log
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J12_combo2_8TeV.txt >> m6.log
cp mlfit.root mlfit_toy_combo2.root

cat m1.log m2.log m3.log m5.log m6.log


echo "*** ZnunuHighPt ***" >! n1.html
python diffNuisances.py mlfit_toy_ZnunuHighPt.root -f html >> n1.html
echo "*** ZnunuMedPt ***" >! n2.html
python diffNuisances.py mlfit_toy_ZnunuMedPt.root -f html >> n2.html
echo "*** ZnunuLowPt ***" >! n3.html
python diffNuisances.py mlfit_toy_ZnunuLowPt.root -f html >> n3.html
#echo "*** ZnunuLowCSV ***" >! n4.html
#python diffNuisances.py mlfit_toy_ZnunuLowCSV.root -f html >> n4.html
echo "*** combo ***" >! n5.html
python diffNuisances.py mlfit_toy_combo.root -f html >> n5.html
echo "*** combo2 ***" >! n6.html
python diffNuisances.py mlfit_toy_combo2.root -f html >> n6.html

