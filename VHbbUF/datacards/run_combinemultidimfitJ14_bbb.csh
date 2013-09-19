# REMEMBER: inject signal
# REMEMBER: set ZH & WH as signal
cp zhinv_Zbb_J14_bbb_ZnunuHighPt_8TeV.txt zhinv_Zbb_2D_J14_bbb_ZnunuHighPt_8TeV.txt
cp zhinv_Zbb_J14_bbb_ZnunuMedPt_8TeV.txt zhinv_Zbb_2D_J14_bbb_ZnunuMedPt_8TeV.txt
cp zhinv_Zbb_J14_bbb_ZnunuLowPt_8TeV.txt zhinv_Zbb_2D_J14_bbb_ZnunuLowPt_8TeV.txt
sed -i 's/process     13         14         0/process     -2         -1         0/' zhinv_Zbb_2D_J14_bbb_ZnunuHighPt_8TeV.txt
sed -i 's/process     13         14         0/process     -2         -1         0/' zhinv_Zbb_2D_J14_bbb_ZnunuMedPt_8TeV.txt
sed -i 's/process     13         14         0/process     -2         -1         0/' zhinv_Zbb_2D_J14_bbb_ZnunuLowPt_8TeV.txt

combineCards.py ZnunuHighPt_8TeV=zhinv_Zbb_2D_J14_bbb_ZnunuHighPt_8TeV.txt ZnunuMedPt_8TeV=zhinv_Zbb_2D_J14_bbb_ZnunuMedPt_8TeV.txt ZnunuLowPt_8TeV=zhinv_Zbb_2D_J14_bbb_ZnunuLowPt_8TeV.txt > ! zhinv_Zbb_2D_J14_bbb_combo_8TeV.txt
text2workspace.py zhinv_Zbb_2D_J14_bbb_combo_8TeV.txt -m 125 -P HiggsAnalysis.CombinedLimit.InvisibleHiggs:floatingXSInvisBR
#combine -M MultiDimFit zhinv_Zbb_2D_J14_bbb_combo_8TeV.root -m 125 --robustFit=1 --rMin=-5 --rMax=5 --algo=singles --cl=0.68
#combine -M MultiDimFit zhinv_Zbb_2D_J14_bbb_combo_8TeV.root -m 125 --robustFit=1 --rMin=-5 --rMax=5 --algo=contour2d --points=30 --cl=0.68
#combine -M MultiDimFit zhinv_Zbb_2D_J14_bbb_combo_8TeV.root -m 125 --robustFit=1 --rMin=-5 --rMax=5 --algo=grid --points=10000
