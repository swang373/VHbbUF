"""
VHbb Analysis Cuts

See also settings.py
"""

__all__ = [
    # The Cut string subclass. Do not remove.
    'Cut',

    # This module is * imported into settings.py. Any number of cut strings may
    # be defined within this module, but only those cuts listed below can be 
    # referenced directly in the central configuration module.
    'Vbb', 'Vb', 'Vcc', 'Vudsg', 
    'CSV_Loose', 'CSV_Medium', 'CSV_Tight', 
    'VetoLeptons', 'addCenJet30', 'FlagsMET', 'NoQCD',
    'Preselection', 'JSON_Triggers',
]


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

"""
twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Muon_Isolation
twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2
"""

# CSV Working Points
CSV_Loose = 0.605
CSV_Medium = 0.890
CSV_Tight = 0.970

# Lepton Veto
VetoMuon = Cut('(Sum$(aLeptons_looseIdPOG>0 && aLeptons_pfRelIso04<0.50 && abs(aLeptons_pdgId)==13 && aLeptons_pt>5) + Sum$(vLeptons_looseIdPOG>0 && vLeptons_pfRelIso04<0.50 && abs(vLeptons_pdgId)==13 && vLeptons_pt>5))')

VetoElectron = Cut('(Sum$(aLeptons_eleCutIdSpring15_25ns_v1>0 && aLeptons_pfRelIso04<0.50 && abs(aLeptons_pdgId)==11 && aLeptons_pt>5) + Sum$(vLeptons_eleCutIdSpring15_25ns_v1>0 && vLeptons_pfRelIso04<0.50 && abs(vLeptons_pdgId)==11 && vLeptons_pt>5))')

VetoLeptons = VetoMuon + VetoElectron

# Additional Central Jets
addCenJet30 = Cut('Sum$(Jet_pt>30 && abs(Jet_eta)<5.2)-2')

# MET Flags
FlagsMET = Cut('Flag_hbheFilterNew && Flag_eeBadScFilter && Flag_CSCTightHaloFilter && Flag_hbheIsoFilter && Flag_goodVertices')

# QCD Reduction 
NoQCD = (
    FlagsMET &
    Cut('HCSV_reg_pt>150') &
    'MinIf$(abs(TVector2::Phi_mpi_pi(Jet_phi-met_phi)), Jet_pt>20 && abs(Jet_eta)<5.2)>0.5' &
    'Jet_chHEF[0]>0.1' &
    'MinIf$(abs(TVector2::Phi_mpi_pi(DiscardedJet_phi-met_phi))-3.1415, DiscardedJet_pt>20 && abs(DiscardedJet_eta)<5.2)>(0.5-3.1415)'
)

