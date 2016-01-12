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

    def __neg__(self):
        """
        Logical negation. Ex: -a
        """
        if not self:
            return Cut()
        return Cut('!({})'.format(self))

    @op
    def __and__(self, other):
        """
        Logical AND. Ex: a & b
        """
        return Cut('({})&&({})'.format(self, other))

    @op
    def __rand__(self, other):
        return self & other

    @op
    def __or__(self, other):
        """
        Logical OR. Ex: a | b
        """
        return Cut('({})||({})'.format(self, other))

    @op
    def __ror__(self, other):
        return self | other

################
#-- Znn Cuts --#
################

#-- Generator Level Higgs Jets Classification --#

Vbb = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)>=2')
Vb = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)==1')
Vcc = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numCHadrons>0)>=1') & -Vb & -Vbb
Vudsg = -(Vbb | Vb | Vcc)

#-- Control Region Cuts --#

Minimal = (
    Cut('Vtype==2 || Vtype==3 || Vtype==4') & 
    'HCSV_pt>150' & 
    'Jet_btagCSV[hJCidx[1]]>0.3' &
    'min(Jet_pt[hJCidx[0]], Jet_pt[hJCidx[1]])>30' &
    'HCSV_mass<300'
)

noQCD = (
    Cut('MinIf$(abs(TVector2::Phi_mpi_pi(met_phi - Jet_phi)), Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)>0.7') &
    'abs(TVector2::Phi_mpi_pi(met_phi - tkMet_phi))<0.7' &
    'tkMet_pt>30'
)

# Additional Jet Cuts
addCenJet30m0 = Cut('(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)>0')
addCenJet30e0 = Cut('(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0')
addCenJet30e1 = Cut('(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)<=1')

# Lepton Cuts
naddGoodLeptons10e0 = Cut('(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 || aLeptons_jetDR>0.3)) + Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3)))==0')
naddGoodTaus20e0 = Cut('Sum$(TauGood_idDecayMode>=1 && TauGood_idCI3hit>=1 && TauGood_pt>20 && abs(TauGood_eta)<2.3)==0')

signal_loose = Cut('Vtype==4') & 'Jet_btagCSV[hJCidx[1]]>0.605' & 'HCSV_mass<100 || HCSV_mass>140' & naddGoodLeptons10e0 & naddGoodTaus20e0



