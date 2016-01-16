"""
VHbb Analysis Cuts

See also settings.py
"""

def op(func):
    def check_args(self, other):
        other = Cut.convert(other)
        if not self:
            return other
        if not other:
            return self
        return func(self, other)
    return check_args

class Cut(str):
    """
    A subclass of string with operators overridden 
    to emulate the behaviour of ROOT's TCut class.

    This code borrows heavily from RootPy:
    github.com/rootpy/rootpy/blob/master/rootpy/tree/cut.py
    """

    def __new__(cls, cut = ''):
        return super(Cut, cls).__new__(cls, cut)

    @staticmethod
    def convert(other):
        if isinstance(other, Cut):
            return other
        elif isinstance(other, basestring):
            return Cut(other)
        elif other is None:
            return Cut()
        return Cut(str(other))

    # Comparison Operators

    @op
    def __eq__(self, other):
        """
        Equal to. Ex: a == b
        """
        return Cut('({0!s})==({1!s})'.format(self, other))

    @op
    def __ne__(self, other):
        """
        Not equal to. Ex: a != b
        """
        return Cut('({0!s})!=({1!s})'.format(self, other))

    @op
    def __lt__(self, other):
        """
        Less than. Ex: a < b
        """
        return Cut('({0!s})<({1!s})'.format(self, other))

    @op
    def __le__(self, other):
        """
        Less than or equal to. Ex: a <= b
        """
        return Cut('({0!s})<=({1!s})'.format(self, other))
    
    @op
    def __gt__(self, other):
        """
        Greater than. Ex: a > b
        """
        return Cut('({0!s})>({1!s})'.format(self, other))

    @op
    def __ge__(self, other):
        """
        Greater than or equal to. Ex: a >= b
        """
        return Cut('({0!s})>=({1!s})'.format(self, other))

    # Logical Operators

    def __neg__(self):
        """
        Logical negation. Ex: -a
        """
        if not self:
            return Cut()
        return Cut('!({0!s})'.format(self))

    @op
    def __and__(self, other):
        """
        Logical AND. Ex: a & b
        """
        return Cut('({0!s})&&({1!s})'.format(self, other))

    @op
    def __rand__(self, other):
        return other & self

    @op
    def __or__(self, other):
        """
        Logical OR. Ex: a | b
        """
        return Cut('({0!s})||({1!s})'.format(self, other))

    @op
    def __ror__(self, other):
        return other | self

    # Arithmetic Operators

    @op
    def __add__(self, other):
        """
        Addition. Ex: a + b
        """
        return Cut('({0!s})+({1!s})'.format(self, other))

    @op
    def __radd__(self, other):
        return other + self

    @op
    def __mul__(self, other):
        """
        Multiplication. Ex: a * b
        """
        return Cut('({0!s})*({1!s})'.format(self, other))

    @op
    def __rmul__(self, other):
        return other * self


#####################
#-- Register Cuts --#
#####################

"""
This module is * imported into the central configuration module settings.py.
While any number of cut strings may be defined here, register only those cuts
which are considered analysis-wide constants for safety.
"""

__all__ = ['SKIM', 'Vbb', 'Vb', 'Vcc', 'Vudsg', 'TARGET_LUMI', 'DATA_WEIGHT', 'MC_WEIGHT']

#########################
#-- Step1 Ntuple Skim --#
#########################

Preselection = (
    Cut('Vtype>=0') &
    'min(met_pt, mhtJet30)>150' &
    'HCSV_pt>100' &
    'Jet_btagCSV[hJCidx[1]]>0.605' &
    'HCSV_mass<500 && HCSV_mass>0' &
    'abs(TVector2::Phi_mpi_pi(HCSV_phi-met_phi))>0.7' &
    'max(abs(Jet_eta[hJCidx[0]]), abs(Jet_eta[hJCidx[1]]))<2.6' &
    'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>20' &
    'abs(Jet_eta[0])<2.6' &
    'Jet_id[0]>=3' &
    'min(Jet_id[hJCidx[0]], Jet_id[hJCidx[1]])>=3'
)

JSON_Triggers = (
    (Cut('Alt$(HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v, 0)') |
    'Alt$(HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v, 0)' | 
    'HLT_BIT_HLT_PFMET170_NoiseCleaned_v') &
    'json'
)

SKIM = Preselection & JSON_Triggers

#################################################
#-- Generator Level Higgs Jets Classification --#
#################################################

Vbb = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)>=2')
Vb = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)==1')
Vcc = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numCHadrons>0)>=1') & -Vb & -Vbb
Vudsg = -(Vbb | Vb | Vcc)

##################################
#-- Control Region Definitions --#
##################################

# Silvio's Cuts
FlagsMET = Cut('Flag_hbheFilterNew && Flag_eeBadScFilter && Flag_CSCTightHaloFilter && Flag_hbheIsoFilter && Flag_goodVertices')

NoQCD1 = 'MinIf$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi)), Jet_pt>20 && Jet_puId && abs(Jet_eta)<5.2)>0.5'

NoQCD2 = 'Jet_chHEF[0]>0.1'

NoQCD3 = 'MinIf$(abs(TVector2::Phi_mpi_pi(DiscardedJet_phi-met_phi))-3.1415, DiscardedJet_pt>20 && DiscardedJet_puId && abs(DiscardedJet_eta)<5.2)>(0.5-3.1415)'


VetoMuon = Cut('(Sum$(aLeptons_looseIdPOG>0 && aLeptons_pfRelIso04<0.50 && abs(aLeptons_pdgId)==13 && aLeptons_pt>5) + Sum$(vLeptons_looseIdPOG>0 && vLeptons_pfRelIso04<0.50 && abs(vLeptons_pdgId)==13 && vLeptons_pt>5))')

VetoElectron = Cut('(Sum$(aLeptons_eleCutIdSpring15_25ns_v1>0 && aLeptons_pfRelIso04<0.50 && abs(aLeptons_pdgId)==11 && aLeptons_pt>5) + Sum$(vLeptons_eleCutIdSpring15_25ns_v1>0 && vLeptons_pfRelIso04<0.50 && abs(vLeptons_pdgId)==11 && vLeptons_pt>5))')






# Michele's Cuts
# Minimal Preselection
minimal = (
    Cut('Vtype==2 || Vtype==3 || Vtype==4') & 
    'HCSV_pt>150' & 
    'Jet_btagCSV[hJCidx[1]]>0.3' &
    'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30' &
    'HCSV_mass<300'
)

# AntiQCD Cuts
antiQCD = (
    Cut('MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)>0.7') &
    'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))<0.7' &
    'tkMet_pt>30'
)

# Additional Jet Cuts
addCenJet30m0 = Cut('(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)>0')
addCenJet30e0 = Cut('(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0')
addCenJet30e1 = Cut('(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)<=1')

# Additional Lepton Cuts
naddGoodLeptons10e0 = Cut('(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 || aLeptons_jetDR>0.3)) + Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3)))==0')
naddGoodTaus20e0 = Cut('Sum$(TauGood_idDecayMode>=1 && TauGood_idCI3hit>=1 && TauGood_pt>20 && abs(TauGood_eta)<2.3)==0')

# Old Jet Flavor Cuts
light_flavor = Cut('abs(Jet_mcFlavour[hJCidx[0]])!=5 && abs(Jet_mcFlavour[hJCidx[1]])!=5')
heavy_flavor = Cut('abs(Jet_mcFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5')

# Regional Cuts
signal_loose = Cut('Vtype==4') & 'Jet_btagCSV[hJCidx[1]]>0.605' & 'HCSV_mass<100 || HCSV_mass>140' & naddGoodLeptons10e0 & naddGoodTaus20e0
signal_tight = Cut('Vtype==4') & 'Jet_btagCSV[hJCidx[1]]>0.8' & 'HCSV_mass<100 || HCSV_mass>140' & naddGoodLeptons10e0 & naddGoodTaus20e0 & addCenJet30e1
tt = Cut('Vtype==2 || Vtype==3') & 'vLeptons_pt>30' & addCenJet30m0 & 'Jet_btagCSV[hJCidx[0]]>0.97' & 'Jet_btagCSV[hJCidx[1]]<0.97'
z_light = Cut('Vtype==4') & addCenJet30e0 & naddGoodLeptons10e0 & 'Jet_btagCSV[hJCidx[0]]<0.97'
z_bb = Cut('Vtype==4') & 'HCSV_mass<100 || HCSV_mass>140' & addCenJet30e0 & naddGoodLeptons10e0 & 'Jet_btagCSV[hJCidx[1]]>0.8'
w_light = Cut('Vtype==2 || Vtype==3') & 'vLeptons_pt>30' & addCenJet30e0 & 'Jet_btagCSV[hJCidx[0]]<0.97'
w_bb = Cut('Vtype==2 || Vtype==3') & 'vLeptons_pt>30' & addCenJet30e0 & 'Jet_btagCSV[hJCidx[1]]>0.8'
qcd = (
    Cut('MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)<0.7') &
    'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))>0.7' & 
    'tkMet_pt<30'
)

#######################
#-- Event Weighting --#
#######################

# Target luminosity of the data in inverse picobarns (pb-1).
TARGET_LUMI = 2190

# 'puWeight' is currently broken in the Heppy ntuples. Use the custom reweighting below.
DATA_WEIGHT = '1'

MC_WEIGHT = Cut('sign(genWeight)') * TARGET_LUMI * '1./sample_lumi' * 'weight2(nTrueInt)'
 
