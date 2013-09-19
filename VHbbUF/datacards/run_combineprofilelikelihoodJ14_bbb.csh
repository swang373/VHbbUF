echo "*** ZnunuHighPt ***" >! p1.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zbb_J14_bbb_ZnunuHighPt_8TeV.txt >> p1.log
echo "*** ZnunuMedPt  ***" >! p2.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zbb_J14_bbb_ZnunuMedPt_8TeV.txt >> p2.log
echo "*** ZnunuLowPt  ***" >! p3.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zbb_J14_bbb_ZnunuLowPt_8TeV.txt >> p3.log
#echo "*** ZnunuLowCSV ***" >! p4.log
#combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zbb_J14_bbb_ZnunuLowCSV_8TeV.txt >> p4.log
echo "*** combo       ***" >! p5.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zbb_J14_bbb_combo_8TeV.txt >> p5.log
echo "*** combo2      ***" >! p6.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --toysFreq --expectSignal=1 zhinv_Zbb_J14_bbb_combo2_8TeV.txt >> p6.log

cat p1.log p2.log p3.log p5.log p6.log
