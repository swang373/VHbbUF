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

# Generator Level Higgs Jets Classification
Vbb = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)>=2')
Vb = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)==1')
Vcc = Cut('Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numCHadrons>0)>=1') & -Vb & -Vbb
Vudsg = -(Vbb | Vb | Vcc) 





