echo "*** ZnunuHighPt ***" >! P1.log
combine -M ProfileLikelihood -m 125 --signif --pvalue --minimizerTolerance=1 vhbb_Znn_J12_bbb_ZnunuHighPt_8TeV.txt >> P1.log
echo "*** ZnunuMedPt  ***" >! P2.log
combine -M ProfileLikelihood -m 125 --signif --pvalue --minimizerTolerance=1 vhbb_Znn_J12_bbb_ZnunuMedPt_8TeV.txt >> P2.log
echo "*** ZnunuLowPt  ***" >! P3.log
combine -M ProfileLikelihood -m 125 --signif --pvalue --minimizerTolerance=1 vhbb_Znn_J12_bbb_ZnunuLowPt_8TeV.txt >> P3.log
#echo "*** ZnunuLowCSV ***" >! P4.log
#combine -M ProfileLikelihood -m 125 --signif --pvalue --minimizerTolerance=1 vhbb_Znn_J12_bbb_ZnunuLowCSV_8TeV.txt >> P4.log
echo "*** combo       ***" >! P5.log
combine -M ProfileLikelihood -m 125 --signif --pvalue --minimizerTolerance=1 vhbb_Znn_J12_bbb_combo_8TeV.txt >> P5.log
echo "*** combo2      ***" >! P6.log
combine -M ProfileLikelihood -m 125 --signif --pvalue --minimizerTolerance=1 vhbb_Znn_J12_bbb_combo2_8TeV.txt >> P6.log

cat P1.log P2.log P3.log P5.log P6.log
