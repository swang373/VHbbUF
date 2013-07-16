echo "*** ZnunuHighPt ***" >! A1.log
combine -M Asymptotic vhbb_Znn_J12_ZnunuHighPt_8TeV.txt >> A1.log
echo "*** ZnunuMedPt  ***" >! A2.log
combine -M Asymptotic vhbb_Znn_J12_ZnunuMedPt_8TeV.txt >> A2.log
echo "*** ZnunuLowPt  ***" >! A3.log
combine -M Asymptotic vhbb_Znn_J12_ZnunuLowPt_8TeV.txt >> A3.log
#echo "*** ZnunuLowCSV ***" >! A4.log
#combine -M Asymptotic vhbb_Znn_J12_ZnunuLowCSV_8TeV.txt >> A4.log
echo "*** combo       ***" >! A5.log
combine -M Asymptotic vhbb_Znn_J12_combo_8TeV.txt >> A5.log
echo "*** combo2      ***" >! A6.log
combine -M Asymptotic vhbb_Znn_J12_combo2_8TeV.txt >> A6.log

cat A1.log A2.log A3.log A5.log A6.log
