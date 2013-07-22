from HiggsAnalysis.CombinedLimit.PhysicsModel import *
#from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
#import ROOT, os

### This base class implements signal yields by production and decay mode
### Specific models can be obtained redefining getHiggsSignalYieldScale
def getHiggsProdDecMode_inv(bin,process,options):
    """Return a triple of (production, decay, energy)"""
    processSource = process
    decaySource   = options.fileName+":"+bin # by default, decay comes from the datacard name or bin label
    if "_" in process: (processSource, decaySource) = process.split("_")
    if processSource not in ["ggH", "qqH", "VH", "WH", "ZH", "ttH", "ZbbHinv"]:
        raise RuntimeError, "Validation Error: signal process %s not among the allowed ones." % processSource
    #
    foundDecay = None
    for D in [ "hww", "hzz", "hgg", "htt", "hbb", "hzg", "hmm", "hinv" ]:
        if D in decaySource:
            if foundDecay: raise RuntimeError, "Validation Error: decay string %s contains multiple known decay names" % decaySource
            foundDecay = D
    if not foundDecay: raise RuntimeError, "Validation Error: decay string %s does not contain any known decay name" % decaySource
    #
    foundEnergy = None
    for D in [ '7TeV', '8TeV', '14TeV' ]:
        if D in decaySource:
            if foundEnergy: raise RuntimeError, "Validation Error: decay string %s contains multiple known energies" % decaySource
            foundEnergy = D
    if not foundEnergy:
        foundEnergy = '7TeV' ## To ensure backward compatibility
        print "Warning: decay string %s does not contain any known energy, assuming %s" % (decaySource, foundEnergy)
    #
    return (processSource, foundDecay, foundEnergy)


class FloatingXSInvisBR(SMLikeHiggsModel):
    "Float independently Higgs cross section and invisible branching ratio"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverted order. Second must be larger the first"
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        
        self.modelBuilder.doVar("r_VH[1,0,10]")
        self.modelBuilder.doVar("BRInvUndet[0,0.0,1.0]")
        self.modelBuilder.factory_("expr::r_VH_SM(\"@0 * (1-@1)\",r_VH,BRInvUndet)")
        self.modelBuilder.factory_("expr::r_VH_BSM(\"@0 * @1\",r_VH,BRInvUndet)")
        
        poi = "r_VH,BRInvUndet"
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
            poi+=',MH'
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
        self.modelBuilder.doSet("POI",poi)
    def getHiggsSignalYieldScale(self,production,decay,energy):
        if production in [ "WH", "ZH", "VH" ]: return "r_VH_SM"
        if production == "ZbbHinv": return "r_VH_BSM"
        raise RuntimeError, "Unknown production mode '%s'" % production
    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        if not self.DC.isSignal[process]: return 1
        (processSource, foundDecay, foundEnergy) = getHiggsProdDecMode_inv(bin,process,self.options)
        return self.getHiggsSignalYieldScale(processSource, foundDecay, foundEnergy)


floatingXSInvisBR = FloatingXSInvisBR()
