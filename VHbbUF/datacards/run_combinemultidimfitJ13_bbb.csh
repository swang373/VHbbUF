#text2workspace.py vhbb_Znn_J13_bbb_combo_8TeV.txt -m 125 -P HiggsAnalysis.CombinedLimit.Diboson:floatingXSVZHF --PO modes=WZHF,ZZHF
#combine -M MultiDimFit vhbb_Znn_J13_bbb_combo_8TeV.root --robustFit=1 --algo=singles --cl=0.68

text2workspace.py vhbb_Znn_J13_bbb_combo_8TeV.txt -m 125 -P HiggsAnalysis.CombinedLimit.Diboson:floatingXSDiboson --PO modes=WZHF,ZZHF
combine -M MultiDimFit vhbb_Znn_J13_bbb_combo_8TeV.root --robustFit=1 --algo=singles --cl=0.68
