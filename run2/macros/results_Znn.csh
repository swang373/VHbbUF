 


combine -M Asymptotic -t -1  znnhbb_HighPt_13TeV.txt > CLs_13_Znn.log
combine -M MaxLikelihoodFit -m 125 -t -1 --expectSignal=1 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm znnhbb_HighPt_13TeV.txt > Mu_13_Exp_Znn.log
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --expectSignal=1 znnhbb_HighPt_13TeV.txt > SigPre_13_Znn.log





