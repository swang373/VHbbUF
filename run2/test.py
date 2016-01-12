import ROOT


# Setup
signal = ROOT.TFile.Open('/afs/cern.ch/work/s/swang373/private/V14/Step3/ZnnH125.root', 'read')
sig_tree = signal.Get('tree')
background = ROOT.TFile.Open('/afs/cern.ch/work/s/swang373/private/V14/Step3/QCD.root', 'read')
bkg_tree = background.Get('tree')

variables = {
    # M(jj)
    'H mass': {
        'expr': 'HCSV_reg_mass',
        'unit': 'GeV',
        'type': 'F'
    },
    # p_T(jj)
    'H p_{T}': {
        'expr': 'HCSV_reg_pt',
        'unit': 'GeV',
        'type': 'F'
    },
    # p_T(j1)
    'H j1 p_{T}': {
        'expr': 'max(Jet_pt_reg[hJCidx[0]], Jet_pt_reg[hJCidx[1]])',
        'unit': 'GeV',
        'type': 'F'
    },
    # p_T(j2)
    'H j2 p_{T}': {
        'expr': 'min(Jet_pt_reg[hJCidx[0]], Jet_pt_reg[hJCidx[1]])',
        'unit': 'GeV',
        'type': 'F'
    },
    # p_T(V) (same as MET for Znn)
    'Type1 Corr. PF #slash{E}_{T}': {
        'expr': 'met_pt',
        'unit': 'GeV',
        'type': 'F'
    },
    # H Jet Max CSV
    'CSV_{max}(j1,j2)': {
        'expr': 'maxCSV := max(0, max(Jet_btagCSV[hJCidx[0]], Jet_btagCSV[hJCidx[1]]))',
        'unit': '',
        'type': 'F'
    },
    # H Jet Min CSV
    'CSV_{min}(j1,j2)': {
        'expr': 'minCSV := max(0, min(Jet_btagCSV[hJCidx[0]], Jet_btagCSV[hJCidx[1]]))',
        'unit': '',
        'type': 'F'
    },
    # dPhi(H,V)
    '#||{#Delta #varphi(H,PF #slash{E}_{T}}': {
        'expr': 'abs(HVdPhi)',
        'unit': '',
        'type': 'F'
    },
    # dEta(jj)
    '#||{#Delta #eta(j1,j2)}': {
        'expr': 'dEta_jj := abs(Jet_eta[hJCidx[0]] - Jet_eta[hJCidx[1]])',
        'unit': '',
        'type': 'F'
    },
    # dR(jj)
    '#Delta R(j1,j2)': {
        'expr': 'deltaR_jj',
        'unit': '',
        'type': 'F'
    },
    # N_aj
    '# Add. Jets {p_{T}>25}': {
        'expr': 'naJets_Znn := max(0, Sum$(Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)-2)',
        'unit': '',
        'type': 'I'
    },


    #deltaPhi(pfMET,J) ADD SUPPORT
    #mindPhiMETJet_dPhi :=  MinIf$(abs(deltaPhi(met_phi,Jet_phi)), Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)
    #min #||{#Delta #varphi(pfMET,j25)}
    #''
    #F

    # Add. Jet Max CSV
    'CSV_{max}(Add. CJ 20)': {
        'expr': 'maxAddCSV := MaxIf$(max(Jet_btagCSV[aJCidx],0), Jet_pt[aJCidx]>20 && abs(Jet_eta[aJCidx])<2.5 && Jet_puId[aJCidx]==1)',
        'unit': '',
        'type': 'F'
    },

    #mindeltaR(H,aj) ADD SUPPORT
    #mindRAddJetH := Min$(deltaR(Jet_eta[hJCidx], Jet_phi[hJCidx], Jet_eta[aJCidx], Jet_phi[aJCidx]))
    #min #DeltaR(H, add. j25)
    #''
    #F

}


#for var, prop in vardict.iteritems():
#     print '{},{title},{units},{type}'.format(var, **properties)




# Run TMVA
ROOT.TMVA.Tools.Instance()

tmva_outfile = ROOT.TFile('TMVA_BDT.root', 'recreate')

factory = ROOT.TMVA.Factory('TMVAClassification', 
                            tmva_outfile,
                            #'!V:!Silent:Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification'
                            '!Silent:Transformations=I:AnalysisType=Classification')


for title, var in variables.iteritems():
    factory.AddVariable(var['expr'], title, var['unit'], var['type'])

factory.AddSignalTree(sig_tree, 2.0, 'train')
factory.AddBackgroundTree(bkg_tree, 2.0, 'train')
factory.SetSignalWeightExpression('sample_lumi')
factory.SetBackgroundWeightExpression('sample_lumi')

sigCut = ROOT.TCut('')
bkgCut = ROOT.TCut('')
factory.PrepareTrainingAndTestTree(sigCut, bkgCut, 'nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None')


