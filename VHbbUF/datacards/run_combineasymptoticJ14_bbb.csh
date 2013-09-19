echo "*** ZnunuHighPt ***" >! a1.log
combine -M Asymptotic -t -1 zhinv_Zbb_J14_bbb_ZnunuHighPt_8TeV.txt >> a1.log
echo "*** ZnunuMedPt  ***" >! a2.log
combine -M Asymptotic -t -1 zhinv_Zbb_J14_bbb_ZnunuMedPt_8TeV.txt >> a2.log
echo "*** ZnunuLowPt  ***" >! a3.log
combine -M Asymptotic -t -1 zhinv_Zbb_J14_bbb_ZnunuLowPt_8TeV.txt >> a3.log
#echo "*** ZnunuLowCSV ***" >! a4.log
#combine -M Asymptotic -t -1 zhinv_Zbb_J14_bbb_ZnunuLowCSV_8TeV.txt >> a4.log
echo "*** combo       ***" >! a5.log
combine -M Asymptotic -t -1 zhinv_Zbb_J14_bbb_combo_8TeV.txt >> a5.log
echo "*** combo2      ***" >! a6.log
combine -M Asymptotic -t -1 zhinv_Zbb_J14_bbb_combo2_8TeV.txt >> a6.log

cat a1.log a2.log a3.log a5.log a6.log
