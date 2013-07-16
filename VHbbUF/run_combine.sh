combineCards.py ZnunuHighPt_8TeV=vhbb_Znn_ZnunuHighPt_8TeV.txt ZnunuLowPt_8TeV=vhbb_Znn_ZnunuLowPt_8TeV.txt > ! vhbb_Znn_combo2_8TeV.txt
combineCards.py ZnunuHighPt_8TeV=vhbb_Znn_ZnunuHighPt_8TeV.txt ZnunuMedPt_8TeV=vhbb_Znn_ZnunuMedPt_8TeV.txt ZnunuLowPt_8TeV=vhbb_Znn_ZnunuLowPt_8TeV.txt > ! vhbb_Znn_combo_8TeV.txt

#EXP
combine -M Asymptotic -t -1 vhbb_Znn_ZnunuHighPt_8TeV.txt
#OBS
combine -M Asymptotic vhbb_Znn_ZnunuHighPt_8TeV.txt

#EXP:
combine -m 125 -M ProfileLikelihood --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_Znn_ZnunuHighPt_8TeV.txt
#OBS:
combine -m 125 -M ProfileLikelihood --signif --pvalue vhbb_Znn_ZnunuHighPt_8TeV.txt
