echo "*** combo       ***" >! SI1.log
combine -M Asymptotic vhbb_Znn_SI_J12_combo_8TeV.txt >> SI1.log

echo "*** combo       ***" >! SI2.log
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_Znn_SI_J12_combo_8TeV.txt >> SI2.log
