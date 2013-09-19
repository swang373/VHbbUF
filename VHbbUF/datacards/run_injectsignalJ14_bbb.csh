echo "*** combo       ***" >! SI1.log
combine -M Asymptotic zhinv_Zbb_SI_J14_bbb_combo_8TeV.txt >> SI1.log

echo "*** combo       ***" >! SI2.log
combine -M ProfileLikelihood -m 125 --signif --pvalue zhinv_Zbb_SI_J14_bbb_combo_8TeV.txt >> SI2.log
