Z(nn)H(bb) Analysis, Step by Step
NOTE: All the scripts must be run from the base directory.

1. Use maketier2list.py to transfer Step 2 ntuples from Pisa to FNAL.
   Edit the directories in the script before running.
   '''
   python maketier2list.py
   '''

2. Use `Skim.C` or `Skim_backup.C` to skim the Step 2 ntuples with baseline selection. 

   Update `HelperNtuples.h`. In `inputstep2.ini`,
   - Section [Skim], edit tagMC, tagData, baseline, mettrigger, metfilter.
   - Section [Stitch], edit xxxLHECUTs.
   In `pyhelper.py`, enable reader.write_HelperNtuples() and disable the rest.
   '''
   python pyhelper.py > HelperNtuples.h
   '''
   To run the skim,
   '''
   source run_Skim.sh
   '''

3. Use `SkimRegression.C` to skim the Step 2 ntuples for BDTG regression training. 

   Update `HelperNtuples.h`. In `inputstep2.ini`,
   - Section [Skim], edit regression, fjregression.
   In `pyhelper.py`, enable reader.write_HelperNtuples() and disable the rest.
   '''
   python pyhelper.py > HelperNtuples.h
   '''
   To run the skim for ak5 jet regression,
   '''
   source run_SkimRegression.sh
   '''
   To run the skim for filter jet regression,
   '''
   source run_SkimRegressionFJ.sh
   '''

4. Use `TrainRegression.C` and `TrainRegressionFJ.C` to produce the BDT regression .xml weight files. 

   Update `HelperTMVA.h`. In `inputstep2.ini`,
   - Section [BDT Regression Variable], edit the variables to use.
   - Section [BDT Regression FJ Variable], edit the variables to use.
   In `pyhelper.py`, enable reader.write_HelperTMVA() and disable the rest.
   '''
   python pyhelper.py > HelperTMVA.h
   '''
   To run the BDT regression and copy the TMVA output and weight files,
   '''
   python run_TrainRegression.py && cp weights/TMVARegression_BDTG.weights.xml weights/TMVARegression_BDTG.testweights.xml && cp TMVAReg.root testTMVAReg.root
   '''
   To run the BDT regression for FJ,
   '''
   python run_TrainRegressionFJ.py && cp weights/TMVARegressionFJ_BDTG.weights.xml weights/TMVARegressionFJ_BDTG.testweights.xml && cp TMVARegFJ.root testTMVARegFJ.root
   '''
   Run the following scripts to assess regression performance,
   '''
   python ComparePtResolution.py
   python ComparePtOffset.py
   python CompareMass_sig.py
   '''

5. Update `HelperNtuples.h`
   In `skimmer.py`, enable skimmer.process() and disable the rest.
   '''
   python skimmer.py inputstep2.ini
   '''
   Copy the output to Section [Process] in `inputstep2.ini`.
   In `skimmer.py`, enable skimmer.stitch() and disable the rest.
   '''
   python skimmer.py inputstep2.ini
   '''
   Copy the output to Section [Stitch] in `inputstep2.ini`.
   In `pyhelper.py`, enable reader.write_HelperNtuples() and disable the rest.
   '''
   python pyhelper.py > HelperNtuples.h
   '''

6. Use `GrowTree.C` to create the Step 3 ntuples.
   '''
   python run_GrowTree.py
   '''
   Currently, I issue instead the following command.
   '''
   source runGrowTree.sh
   '''

7. Use `SkimClassification.C` to skim the Step 3 ntuples for BDT classification.
   This is currently bugged. For now, remove the string "ZnunuHighPt_VH:= " in `HelperTMVA.h`.
   '''
   source run_SkimClassification.sh
   '''

8. Use `TrainBDT.C` to produce the BDT classification .xml weight files.
   '''
   source run_TrainBDT.sh
   '''

9. Use `TrimTree.C` to create the Step 4 ntuples.
   
   In `inputstep2.ini`,
   - Section [Weight], edit the MC event weights.
   - Section [Trigger], edit the Data trigger.
   - Section [Selection], edit the signal and control regions.
   In `pyhelper.py`, enable reader.write_HelperTMVA() and disable the rest.
   '''
   python pyhelper.py > HelperTMVA.h
   '''

   In `TrimTree.C`, edit the BDT .xml weight files to load.
   '''
   python run_TrimTree.py
   '''

 stitch the Step 4's (remember to edit $DIR in the script)
source run_stitch.sh

10. run for now
cd macros
root -l -b -q plotHistos_Znn_13TeV_BDT.C++
source results_Znn.csh


later we should do
10. Use `ScaleFactorJ11.C` and `ScaleFactorWorkspaceJ11.C` to fit scale factors.
.. code:: bash
 set the variable indir in the scripts, reset the fitresults arrays, reset the iteration number in calc_scalefactors(...)
root -l -b -q ScaleFactorJ11.C++O
root -l -b -q ScaleFactorWorkspaceJ11.C++O
 use MaxLikelihoodFit to get the nuisances
combineCards.py ZnunuHighPt_WjLF=vhbb_Znn_SF_J11_ZnunuHighPt_WjLF_8TeV.txt ZnunuHighPt_WjHF=vhbb_Znn_SF_J11_ZnunuHighPt_WjHF_8TeV.txt ZnunuHighPt_ZjLF=vhbb_Znn_SF_J11_ZnunuHighPt_ZjLF_8TeV.txt ZnunuHighPt_ZjHF=vhbb_Znn_SF_J11_ZnunuHighPt_ZjHF_8TeV.txt ZnunuHighPt_TT=vhbb_Znn_SF_J11_ZnunuHighPt_TT_8TeV.txt > vhbb_Znn_SF_J11_ZnunuHighPt_8TeV.txt
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_SF_J11_ZnunuHighPt_8TeV.txt
python printNuisances.py mlfit.root

11. Use `BDTShapeJ12.C` and `BDTShapeWorkspaceJ12.C` to get final limit and significance.
 set the variable indir in the scripts, adjust binning, then run
root -l -b -q BDTShapeJ12.C++O
root -l -b -q BDTShapeWorkspaceJ12.C++O
 calculate expected limit and significance
combine -M Asymptotic -t -1 vhbb_Znn_J11_ZnunuHighPt_8TeV.txt
combine -M ProfileLikelihood -m 125 --signif --pvalue -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J11_ZnunuHighPt_8TeV.txt
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm -t -1 --toysFreq --expectSignal=1 vhbb_Znn_J11_ZnunuHighPt_8TeV.txt
python diffNuisances.py -f html mlfit.root
 calculate observed limit and significance
combine -M Asymptotic vhbb_Znn_J11_ZnunuHighPt_8TeV.txt
combine -M ProfileLikelihood -m 125 --signif --pvalue vhbb_Znn_J11_ZnunuHighPt_8TeV.txt
combine -M MaxLikelihoodFit -m 125 --robustFit=1 --stepSize=0.05 --rMin=-5 --rMax=5 --saveNorm vhbb_Znn_J11_ZnunuHighPt_8TeV.txt
python diffNuisances.py -f html mlfit.root

Tuning

1. In directory `tuneBDT`, use `prepare_tuneBDT.sh` to create symbolic links. This step is only needed at the first time.
 edit the links if necessary
source prepare_tuneBDT.sh

2. Use `BDTShapeJ12_reload.C.diff` to patch `BDTShapeJ12_reload.C`.

 the patch was created with `diff -u BDTShapeJ12_reload.C.orig BDTShapeJ12_reload.C`
 it might only work with a specific version of BDTShapeJ12.C
patch -p0 < BDTShapeJ12_reload.C.diff

3. Use `run_tune.py` to get commands to train and reload.
 edit r, n, w if necessary
python run_tune.py
 first train, then reload

4. Use `run_combine_tune.py` to calculate final limit and significance.
.. code:: bash
 edit r, n, w, logs if necessary
python run_combine_tune.py

5. Use `retrieve_combine_tune.py` to store the previous results into TTrees.
.. code:: bash
 edit r, w, logs if necessary
 it only does one channel at a time
python run_combine_tune.py
 use `hadd` to combine the TTrees

6. Open the combined TTree (e.g. TMVA_ZnunuHighPt_new.root), scan for the best significance.
.. code:: bash
root [0] fomtree->Scan("NTrees:nEventsMin:MaxDepth:COMB_limit:COMB_signif:TMVA_kolS:TMVA_kolB:USER_psig","COMB_limit<3.6 && COMB_signif>0.68")
98. BDT binning optimization instructions pending...
99. Find yield uncertainties (only for ZbbHinv cards).
.. code:: bash
 edit nuisances, processes, soverb if necessary
python systematica_ZbbHinv.py

