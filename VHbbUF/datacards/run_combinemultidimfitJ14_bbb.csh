# REMEMBER: inject signal
# REMEMBER: set ZH & WH as signal
cp vhbb_Znn_J14_bbb_ZnunuHighPt_8TeV.txt vhbb_Znn_2D_J14_bbb_ZnunuHighPt_8TeV.txt
cp vhbb_Znn_J14_bbb_ZnunuMedPt_8TeV.txt vhbb_Znn_2D_J14_bbb_ZnunuMedPt_8TeV.txt
cp vhbb_Znn_J14_bbb_ZnunuLowPt_8TeV.txt vhbb_Znn_2D_J14_bbb_ZnunuLowPt_8TeV.txt
sed -i 's/process     13         14         0/process     -2         -1         0/' vhbb_Znn_2D_J14_bbb_ZnunuHighPt_8TeV.txt
sed -i 's/process     13         14         0/process     -2         -1         0/' vhbb_Znn_2D_J14_bbb_ZnunuMedPt_8TeV.txt
sed -i 's/process     13         14         0/process     -2         -1         0/' vhbb_Znn_2D_J14_bbb_ZnunuLowPt_8TeV.txt

combineCards.py ZnunuHighPt_8TeV=vhbb_Znn_2D_J14_bbb_ZnunuHighPt_8TeV.txt ZnunuMedPt_8TeV=vhbb_Znn_2D_J14_bbb_ZnunuMedPt_8TeV.txt ZnunuLowPt_8TeV=vhbb_Znn_2D_J14_bbb_ZnunuLowPt_8TeV.txt > ! vhbb_Znn_2D_J14_bbb_combo_8TeV.txt
text2workspace.py vhbb_Znn_2D_J14_bbb_combo_8TeV.txt -m 125 -P HiggsAnalysis.CombinedLimit.InvisibleHiggs:floatingXSInvisBR
#combine -M MultiDimFit vhbb_Znn_2D_J14_bbb_combo_8TeV.root -m 125 --robustFit=1 --rMin=-5 --rMax=5 --algo=singles --cl=0.68
#combine -M MultiDimFit vhbb_Znn_2D_J14_bbb_combo_8TeV.root -m 125 --robustFit=1 --rMin=-5 --rMax=5 --algo=contour2d --points=30 --cl=0.68
#combine -M MultiDimFit vhbb_Znn_2D_J14_bbb_combo_8TeV.root -m 125 --robustFit=1 --rMin=-5 --rMax=5 --algo=grid --points=10000
