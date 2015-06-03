#!/usr/bin/env python

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))

def rgb_to_hex(r, g, b):
    return '#%02X%02X%02X' % (r, g, b)

def uniquify(seq):
    # Order preserving
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

def keysorted(seq, i=0):
    tuples = []
    from operator import itemgetter
    return sorted(seq, key=itemgetter(i))


class ProcessBox:
    """A container of xsec, count, countwithpu, color"""
    
    def __init__(self, value):
        num = 4
        if not value.count(" , ") == num-1:
            raise ValueError("%s: Expected %i values: \"%s\"" % (self.__class__.__name__, num, value))
        values = [v.strip() for v in value.split(" , ")]
        self.xsec = eval(values[0])
        self.count = eval(values[1])
        self.countwithpu = eval(values[2])
        self.color = eval(values[3])
    
    def __repr__(self):
        return self.get().__repr__()
    
    def get(self):
        return (self.xsec, self.count, self.countwithpu, self.color)


class HistoBox:
    """A container of xlabel, unit, type, nbinsx, xlow, xup"""
    
    def __init__(self, value):
        num = 6
        if not value.count(" , ") == num-1:
            raise ValueError("%s: Expected %i values: \"%s\"" % (self.__class__.__name__, num, value))
        values = [v.strip() for v in value.split(" , ")]
        self.xlabel = eval(values[0])
        self.unit = eval(values[1])
        self.type = eval(values[2])
        self.nbinsx = eval(values[3])
        self.xlow = eval(values[4])
        self.xup = eval(values[5])
        self.xtitle = self.xlabel
        if self.unit:
            self.xtitle = "%s [%s]" % (self.xlabel, self.unit)
    
    def __repr__(self):
        return self.get().__repr__()
    
    def get(self):
        return (self.xlabel, self.unit, self.type, self.nbinsx, self.xlow, self.xup)


import warnings
from collections import OrderedDict
#from ordereddict import OrderedDict 
import ConfigParser
import re

class Reader:
    """
    I read the .ini config
    """

    def __init__(self):
        self.cfgparser = ConfigParser.ConfigParser(dict_type=OrderedDict)
        self.cfgparser.optionxform = str  #< want to preserve case
        self.has_read = False
        self.sections = [
            "Main",
            "Process",
            "Stitch",
            "Skim",
            "Weight",
            "Trigger",
            "Selection",
            "Mjj Selection",
            "Variable",
            "Histogramming",
            "BDT Variable",
            "BDT Regression Variable",
            "BDT Regression FJ Variable",
            "Systematics",
            ]

    def is_not_a_list(self, value):
        if not value.count(" , ") == 0:
            raise ValueError("%s: Expected only one value: \"%s\"" % (self.__class__.__name__, value))
        return 0

    def is_a_list(self, value, num):
        if not value.count(" , ") == num-1:
            raise ValueError("%s: Expected %i values: \"%s\"" % (self.__class__.__name__, num, value))
        return 0

    def read(self, ini):
        self.cfgparser.read(ini)
        
        if len(self.cfgparser.sections()) != len(self.sections):
            raise ValueError("%s: Expected number of sections: %i; found: %i" % (self.__class__.__name__, len(self.sections), len(self.cfgparser.sections())))

        for sec in self.cfgparser.sections():
            if sec not in self.sections:
                raise ValueError("%s: Unexpected section: [%s]" % (self.__class__.__name__, sec))

        # Set attributes from Section "Main"
        for (k, v) in self.cfgparser.items("Main"):
            setattr(self, k, eval(v))
        
        self.has_read = True
        return 0

    def get_processes(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Process"):
            self.is_a_list(v, 4)
            d[k] = ProcessBox(v)
        return d

    def get_skims(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Skim"):
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_stitch(self):
        assert(self.has_read)
        d = OrderedDict()
        n = 0
        name = ""
        for (k, v) in self.cfgparser.items("Stitch"):
            if n == 0:
                n = len(v.split(" , "))
                name = k[0:k.find("LHECUT")]
            elif not k.startswith(name):
                n = len(v.split(" , "))
                name = k[0:k.find("LHECUT")]
            self.is_a_list(v, n)
            d[k] = eval(v)
        return d

    def get_weights(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Weight"):
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_triggers(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Trigger"):
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_selections(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Selection"):
            if k.endswith("BIN") or k.endswith("VETO"):  # avoid text substitutions
                continue
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_mjj_selections(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Mjj Selection"):
            if k.endswith("BIN") or k.endswith("VETO"):  # avoid text substitutions
                continue
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_variables(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Variable"):
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_histogrammings(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Histogramming"):
            self.is_a_list(v, 6)
            d[k] = HistoBox(v)
        return d

    def get_bdt_variables(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("BDT Variable"):
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_bdt_regression_variables(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("BDT Regression Variable"):
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_bdt_regression_fj_variables(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("BDT Regression FJ Variable"):
            self.is_not_a_list(v)
            d[k] = eval(v)
        return d

    def get_systematics(self):
        assert(self.has_read)
        d = OrderedDict()
        for (k, v) in self.cfgparser.items("Systematics"):
            vv = eval(v)
            assert len(vv) == 2
            vv0 = [a.strip() for a in vv[0].split(" , ")]
            vv1 = [a.strip() for a in vv[1].split(" , ")]
            assert len(vv0) == len(vv1)
            d[k] = zip(vv0, vv1)
        return d

    def get_regression_substitutions(self):
        substitutions = {}
        attr = "apply_bdt_regression"
        if not hasattr(self, attr):
            raise AttributeError("%s: Failed to find attribute: %s" % (self.__class__.__name__, attr))
        if getattr(self, attr):
            substitutions["Jet_pt"] = "Jet_ptReg"
            substitutions["H\.pt"] = "HptReg"
            substitutions["H\.mass"] = "HmassReg"
        return substitutions

    def write_HelperNtuples(self):
        self.write_skim()
        self.write_lumi()
        self.write_stitch()


    def write_HelperTMVA(self):
        self.write_regression()
        self.write_regression_fj()
        self.write_classification()
        self.write_systematics()

    def write_skim(self):
        strings = []
        skims = self.get_skims()
        strings += skims.iteritems()

        for (k, v) in strings:
            print "const std::string {skim:12} = \"{skimstr}\";".format(skim=k, skimstr=v)
        print

    def write_lumi(self):
        processes = self.get_processes()
        attr = "lumi"
        if not hasattr(self, attr):
            raise AttributeError("%s: Failed to find attribute: %s" % (self.__class__.__name__, attr))
        strings = []
        for (k, v) in processes.iteritems():
            if k == "Data":
                string = "{lumi:.1f}".format(lumi=1.0)
            else:
                string = "lumi * {xsec:>14f} / {count:>13.4f}".format(xsec=v.xsec, count=v.countwithpu)
            strings.append((k, string))
        
        print "const double lumi = {lumi:.1f};".format(lumi=getattr(self, attr))
        self.write_C_function("GetLumis", "mapF", strings)

    def write_stitch(self):
        stitches = self.get_stitch()
        strings = []
        for (k, v) in stitches.iteritems():
            for vv in v:
                if "LHECUT" in k:
                    string = "\"{0} := {1}\"".format(k.replace("LHECUT",""), vv)
                    strings.append((len(strings), string))
        self.write_C_function_if("GetLHECuts", "vector", strings)
        
        del strings[:]
        processes = self.get_processes()
        lhecut, lhenorm = [], []
        for (k, v) in stitches.iteritems():
            if "LHECUT" in k:
                lhecut = v
            
            elif "LHENORM" in k:
                lhenorm = ["%.1f * %f / %f" % (getattr(self, "lumi"), processes[p].xsec, processes[p].countwithpu) for p in v]
            
            elif "LHECOUNT" in k:
                string = "\""
                for ic, c in enumerate(v):
                    string += ("((%s) * %s * %f)" % (lhecut[ic], lhenorm[ic], c))
                    if ic != len(v)-1:  string += " + "
                string += "\""
                strings.append((k.replace("LHECOUNT",""),string))
        self.write_C_function("GetLHEWeights", "map", strings)
    
    def write_regression(self):
        bdt_regression_variables = self.get_bdt_regression_variables()
        variables = self.get_variables()
        histogrammings = self.get_histogrammings()
        
        vstrings, vlabelstrings, vstrings0, vstrings1 = [], [], [], []
        
        # Setup variables
        for (k, v) in bdt_regression_variables.iteritems():
            if not k in variables or not k in histogrammings:
                raise ValueError("%s: Cannot find variable: %s" % (self.__class__.__name__, k))
            if not v == "V":
                raise ValueError("%s: Expect only \"V\" in [BDT Regression Variable]." % (self.__class__.__name__))
            
            expr = variables[k]
            h = histogrammings[k]
            string = "\"{0} := {1}\"".format(k, expr)
            vstrings.append((len(vstrings), string))
            string = "\"{0};{1};{2}\"".format(h.xlabel, h.unit, h.type)
            vlabelstrings.append((len(vlabelstrings), string))
            
            expr0 , expr1 = expr, expr
            if not k.startswith("breg_evt"):
                expr0 = re.sub(r"hJCidx", r"hJCidx[0]", expr)
                expr1 = re.sub(r"hJCidx", r"hJCidx[1]", expr)
            string = "\"{0}\"".format(expr0)
            vstrings0.append((len(vstrings0), string))
            string = "\"{0}\"".format(expr1)
            vstrings1.append((len(vstrings1), string))
        
        print "/// TMVA Regression"
        self.write_C_function("GetInputExpressionsReg", "vector", vstrings)
        self.write_C_function("GetInputExpressionLabelsReg", "vector", vlabelstrings)
        self.write_C_function("GetInputExpressionsReg0", "vector", vstrings0)
        self.write_C_function("GetInputExpressionsReg1", "vector", vstrings1)

    def write_regression_fj(self):
        bdt_regression_variables = self.get_bdt_regression_fj_variables()
        variables = self.get_variables()
        histogrammings = self.get_histogrammings()
        
        vstrings, vlabelstrings, vstrings0, vstrings1, vstrings2 = [], [], [], [], []
        
        # Setup variables
        for (k, v) in bdt_regression_variables.iteritems():
            if not k in variables or not k in histogrammings:
                raise ValueError("%s: Cannot find variable: %s" % (self.__class__.__name__, k))
            if not v == "V":
                raise ValueError("%s: Expect only \"V\" in [BDT Regression FJ Variable]." % (self.__class__.__name__))
            
            expr = variables[k]
            h = histogrammings[k]
            string = "\"{0} := {1}\"".format(k, expr)
            vstrings.append((len(vstrings), string))
            string = "\"{0};{1};{2}\"".format(h.xlabel, h.unit, h.type)
            vlabelstrings.append((len(vlabelstrings), string))
            
            expr0, expr1, expr2 = expr, expr, expr
            if not k.startswith("breg_evt"):
                expr0 = re.sub(r"fathFilterJets_(\w+)", r"fathFilterJets_\1[0]", expr)
                expr1 = re.sub(r"fathFilterJets_(\w+)", r"fathFilterJets_\1[1]", expr)
                expr2 = re.sub(r"fathFilterJets_(\w+)", r"fathFilterJets_\1[2]", expr)
            string = "\"{0}\"".format(expr0)
            vstrings0.append((len(vstrings0), string))
            string = "\"{0}\"".format(expr1)
            vstrings1.append((len(vstrings1), string))
            string = "\"{0}\"".format(expr2)
            vstrings2.append((len(vstrings2), string))
        
        print "/// TMVA Regression (filter jets)"
        self.write_C_function("GetInputExpressionsFJReg", "vector", vstrings)
        self.write_C_function("GetInputExpressionLabelsFJReg", "vector", vlabelstrings)
        self.write_C_function("GetInputExpressionsFJReg0", "vector", vstrings0)
        self.write_C_function("GetInputExpressionsFJReg1", "vector", vstrings1)
        self.write_C_function("GetInputExpressionsFJReg2", "vector", vstrings2)

    def write_classification(self):
        bdtvariables = self.get_bdt_variables()
        variables = self.get_variables()
        histogrammings = self.get_histogrammings()
        weights = self.get_weights()
        triggers = self.get_triggers()
        selections = self.get_selections()
        mjjselections = self.get_mjj_selections()
        
        substitutions = self.get_regression_substitutions()
        
        vstrings, vlabelstrings, sstrings, slabelstrings = [], [], [], []
        histstrings, weigstrings, trigstrings, selstrings, mjjselstrings = [], [], [], [], []
        
        # Setup variables and spectators
        for (k, v) in bdtvariables.iteritems():
            if not k in variables or not k in histogrammings:
                raise ValueError("%s: Cannot find variable: %s" % (self.__class__.__name__, k))

            expr = variables[k]
            h = histogrammings[k]
            for old, new in substitutions.iteritems():
                if k.startswith("breg"):  continue
                expr = re.sub(r"(\A|\W)%s(\Z|\W)" % old, r"\1%s\2" % new, expr)
                #expr = expr.replace(old, new)

            if v == "V":
                string = "\"{0} := {1}\"".format(k, expr)
                vstrings.append((len(vstrings), string))
                string = "\"{0};{1};{2}\"".format(h.xlabel, h.unit, h.type)
                vlabelstrings.append((len(vlabelstrings), string))
            elif v == "S":
                string = "\"{0} := {1}\"".format(k, expr)
                sstrings.append((len(sstrings), string))
                string = "\"{0};{1};{2}\"".format(h.xlabel, h.unit, h.type)
                slabelstrings.append((len(slabelstrings), string))
            string = "\"{0};{1};{2};{3:.2f};{4:.2f}\"".format(k, h.xtitle, h.nbinsx, h.xlow, h.xup)
            histstrings.append((len(histstrings), string))
        
        # Setup weights
        for (k, v) in weights.iteritems():
            weigstrings.append((k, "\"%s\"" % v))
        
        # Setup triggers
        for (k, v) in triggers.iteritems():
            trigstrings.append((k, "\"%s\"" % v))
        
        # Setup selections
        for (k, v) in selections.iteritems():
            for old, new in substitutions.iteritems():
                v = re.sub(r"(\A|\W)%s(\Z|\W)" % old, r"\1%s\2" % new, v)
                #v = v.replace(old, new)
            string = "\"{0} := {1}\"".format(k, v)
            selstrings.append((len(selstrings), string))
        
        # Setup Mjj selections
        for (k, v) in mjjselections.iteritems():
            for old, new in substitutions.iteritems():
                v = re.sub(r"(\A|\W)%s(\Z|\W)" % old, r"\1%s\2" % new, v)
                #v = v.replace(old, new)
            string = "\"{0} := {1}\"".format(k, v)
            mjjselstrings.append((len(mjjselstrings), string))
        
        print "/// TMVA Classification"
        if getattr(self, "apply_bdt_regression"):
            print "//! BDT regression is applied"
        #self.write_C_function("GetInputExpressions", "vector", vstrings)
        #self.write_C_function("GetInputExpressionLabels", "vector", vlabelstrings)
        #self.write_C_function("GetSpecExpressions", "vector", sstrings)
        #self.write_C_function("GetSpecExpressionLabels", "vector", slabelstrings)
        #self.write_C_function("GetHistogrammings", "vector", histstrings)
        self.write_C_function("GetWeightExpressions", "map", weigstrings)
        self.write_C_function("GetTriggerExpressions", "map", trigstrings)
        self.write_C_function_if("GetSelExpressions", "vector", selstrings)
        self.write_C_function_if("GetSelMjjExpressions", "vector", mjjselstrings)
        
        ik = 0
        old_kk = ""
        print "TObjString selectFlags_id = \"",
        for (k, v) in  selections.iteritems():
            kk = k.split("_")[0]
            if old_kk == "":
                print "\\n//%s" % kk,
                old_kk = kk
            elif old_kk != kk:
                print "\\n//%s" % kk,
                ik = 0
                old_kk = kk
            print "\\nselectFlags[%i][isyst] = \\\"%s\\\";" % (ik, k), 
            ik += 1
        print "\";"
        print
        

    def write_systematics(self):
        systematics = self.get_systematics()  # this is a {[(a,b)]}

        type_ = "std::vector<std::vector<std::pair<std::string, std::string> > >"
        func_ = "GetSystExpressions"
        print "%s %s() {" % (type_, func_)
        print "    %s values;" % (type_)
        print "    values.resize(%i);" % (len(systematics))
        ik = 0
        for (k, v) in systematics.iteritems():
            print "    values[%i].resize(%i); // %s" % (ik, len(v), k)
            ipair = 0
            for pair in v:
                print "    values[%i][%i] = make_pair(\"%s\", \"%s\");" % (ik, ipair, pair[0], pair[1])
                ipair += 1
            ik += 1
        print "    return values;"
        print "}"
        print
        
        ik = 0
        print "TObjString systFlags_id = \"",
        #print "\\nsystFlags.resize(%i,\\\"\\\");" % (len(systematics)), 
        for (k, v) in systematics.iteritems():
            print "\\nsystFlags[%i] = \\\"%s\\\";" % (ik, k), 
            ik += 1
        print "\";"
        print

    def write_C_function(self, func, returntype, strings):
        if not strings:
            raise ValueError("%s: No strings" % (self.__class__.__name__))
        types = {
            "vector"  : "std::vector<std::string>",
            "vectorF" : "std::vector<float>",
            "map"     : "std::map<std::string, std::string>",
            "mapF"    : "std::map<std::string, float>",
            }
        print "%s %s() {" % (types[returntype], func)
        print "    %s values;" % (types[returntype])
        if returntype == "vector":
            print "    values.resize(%i, \"\");" % (len(strings))
        elif returntype == "vectorF":
            print "    values.resize(%i, -999.);" % (len(strings))
        for (k, v) in strings:
            if returntype in ["map", "mapF"]:
                print "    values[{k:18s}] = {v};".format(k=("\""+k+"\""), v=v)
            elif returntype in ["vector", "vectorF"]:
                print "    values[{k}] = {v};".format(k=k, v=v)
        print "    return values;"
        print "}"
        print

    def write_C_function_if(self, func, returntype, strings):
        if not strings:
            raise ValueError("%s: No strings" % (self.__class__.__name__))
        types = {
            "vector"  : "std::vector<std::string>",
            "vectorF" : "std::vector<float>",
            "map"     : "std::map<std::string, std::string>",
            "mapF"    : "std::map<std::string, float>",
            }
        print "%s %s(const std::string id) {" % (types[returntype], func)
        print "    %s values;" % (types[returntype])
        strings_if = OrderedDict()
        for (k, v) in strings:
            if not " := " in v:
                raise ValueError("%s: A label is expected: %s" % (self.__class__.__name__, k))
            kk = v.strip("\"").split(" := ")[0].split("_")[0]
            if kk in strings_if:
                strings_if[kk].append((len(strings_if[kk]), v))
            else:
                strings_if[kk] = [(0, v)]
        print_if = True
        for (kk, strings) in strings_if.iteritems():
            if print_if:
                print "    if (id == \"%s\") {" % (kk)
                print_if = False
            else:
                print "    else if (id == \"%s\") {" % (kk)
            if returntype == "vector":
                print "        values.resize(%i, \"\");" % (len(strings))
            elif returntype == "vectorF":
                print "        values.resize(%i, -999.);" % (len(strings))
            for (k, v) in strings:
                if v.startswith("\"" + kk + " := "):
                    v = v.replace(kk + " := ", "")
                if returntype in ["map", "mapF"]:
                    print "        values[{k:18s}] = {v};".format(k=("\""+k+"\""), v=v)
                elif returntype in ["vector", "vectorF"]:
                    print "        values[{k}] = {v};".format(k=k, v=v)
            print "        return values;"
            print "    }"
        print "    return values;"
        print "}"
        print


## Main
## ----------------------------------------------------------------------------

if __name__ == "__main__":

    ini = "inputstep2.ini"

    reader = Reader()
    reader.read(ini)

    # Write codes
    #reader.write_HelperNtuples()
    reader.write_HelperTMVA()

