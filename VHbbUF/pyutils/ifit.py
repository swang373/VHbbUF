#   Example of a fitting program
#   ============================
#
#   The fitting function fcn is a simple chisquare function
#   The data consists of 5 data points (arrays x,y,z) + the errors in errorsz
#   More details on the various functions or parameters for these functions
#   can be obtained in an interactive ROOT session with:
#    Root > TMinuit *minuit = new TMinuit(10);
#    Root > minuit->mnhelp("*",0)  to see the list of possible keywords
#    Root > minuit->mnhelp("SET",0) explains most parameters
#
#   Adapted from $ROOTSYS/tutorials/fit/Ifit.C

from ROOT import gROOT, TMinuit, TMath, TFile, TH1, TH1F, Long, Double
from array import array
from math import sqrt


# Global options
ncount = 0
tolerance = 1e-2
error_on_data = True
xbins = array('d', [130,150,170,190,210,240,270,300,350,400,500])  # variable-sized bins
nxbins = len(xbins)-1

# Global fit options
processes = ["ZjLF", "ZjHF", "WjLF", "WjHF", "TT", "s_Top", "VV"]  # not including data
procerrors = [0.07, 0.15, 0.07, 0.15, 0.08, 0.25, 0.30]
data_obs = "data_obs"
categories = ["ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5"]

npars = 5  # number of processes to fit
nconsts = 2  # number of processes set as constant
assert(len(processes) == npars + nconsts)

# Some arrays
mc = (npars+nconsts)*[array('f')]  # MC histogram contents in python arrays
mcerr = (npars+nconsts)*[array('f')]  # MC histogram errors in python arrays
data = [array('f')]  # data histogram contents
dataerr = [array('f')]  # data histogram errors


##______________________________________________________________________________
def rebin_xbins(h, ngroup, xbins):
    # Do variable-sized rebinning
    hrebin = h.Rebin(ngroup, "%s_1" % h.GetName(), xbins)
    hrebin.Scale(1.0, "width")
    return hrebin

##______________________________________________________________________________
def fill_arrays(h1, scale=1.0):
    # Convert TH1 to 2 python arrays
    x, err = [], []
    for i in xrange(h1.GetNbinsX()+2):
        x.append(h1.GetBinContent(i) * scale)
        err.append(h1.GetBinError(i) * scale)
    ax = array('f', x)
    aerr = array('f', err)
    return (ax, aerr)

def join_arrays(x, err):
    # Join 2 lists of arrays into 2 arrays respectively
    ax = array('f')
    for xx in x:
        ax.extend(xx)
    aerr = array('f')
    for ex in err:
        aerr.extend(ex)
    return (ax, aerr)


##______________________________________________________________________________
def fill(infile, key):
    # Fill histograms and join them
    tfile = TFile.Open(infile)
    
    for ip, p in enumerate(processes + [data_obs]):  # including data 
        lax, laerr = [], []  # lists of arrays
        for c in categories:
            h = tfile.Get(key % (c, p))
            hrebin = rebin_xbins(h, nxbins, xbins)
            ax, aerr = fill_arrays(hrebin)
            lax.append(ax)
            laerr.append(aerr)
        if p == data_obs:
            data[0], dataerr[0] = join_arrays(lax, laerr)
        else:
            mc[ip], mcerr[ip] = join_arrays(lax, laerr)
    tfile.Close()
    return


##______________________________________________________________________________
def func(x, par):
    # Calculate MC sum
    assert(len(x)==len(processes))
    value = (
        par[0]*x[0] +
        par[1]*x[1] +
        par[2]*x[2] +
        par[3]*x[3] +
        par[4]*x[4] +
        1.0   *x[5] +
        1.0   *x[6]
        )
    return value

def error(x, err, systerr):
    # Calculate MC errors
    # err = stat. errors, systerr = syst. errors
    assert(len(x)==len(processes) and len(err)==len(x) and len(systerr)==len(x))
    value = 0.
    for ex in err:
        value += pow(ex, 2)
    for xx, sx in zip(x, systerr):
        value += pow(sx * xx, 2)
    return sqrt(value)

def fcn(npar, gin, f, par, iflag):
    # Calculate chisquare
    global ncount

    chi, chisq = 0., 0.
    for i in xrange(len(data[0])):  # loop over all bins
        mc_ = [m[i] for m in mc]  # [MC bin contents]
        mcerr_ = [m[i] for m in mcerr]  # [MC bin errors]
        # If this bin has too few events, ignore it
        if sum(mc_) < tolerance:
            continue
        if data[0][i] < tolerance:
            continue

        value = func(mc_, par)
        err = error(mc_, mcerr_, procerrors)
        # Use only poisson errors from data
        if error_on_data:
            err = sqrt(data[0][i])

        chi = (data[0][i] - value)/err
        chisq += chi * chi
    f[0] = chisq
    ncount += 1


##______________________________________________________________________________
def ifit():   
    gMinuit = TMinuit(npars)
    gMinuit.SetFCN(fcn)

    arglist = array('d', 10*[0.])
    ierflg = Long(1982)

    arglist[0] = 1  # 1 for chisquared fits, 0.5 for negative log likelihood
    gMinuit.mnexcm("SET ERR", arglist, 1, ierflg)

    arglist[0] = 0.00001  # floating point precision
    gMinuit.mnexcm("SET EPS", arglist, 1, ierflg)

    arglist[0] = 2 # 1 mean fewer function calls, 2 mean more reliable minimization
    gMinuit.mnexcm("SET STR", arglist, 1, ierflg)

    # Set starting values and step sizes for parameters
    vstart = array('d', npars*[1.0])
    #vstart = array('d', [0.909, 1.341, 0.942, 1.300, 1.003])
    step = array('d', npars*[0.001])
    gMinuit.mnparm(0, "ZjLF", vstart[0], step[0], 0, 0, ierflg)
    gMinuit.mnparm(1, "ZjHF", vstart[1], step[1], 0, 0, ierflg)
    gMinuit.mnparm(2, "WjLF", vstart[2], step[2], 0, 0, ierflg)
    gMinuit.mnparm(3, "WjHF", vstart[3], step[3], 1.3, 2.0, ierflg)
    gMinuit.mnparm(4, "TT"  , vstart[4], step[4], 0, 0, ierflg)

    # Fix parameter
    #arglist[0] = 1
    #gMinuit.mnexcm("FIX", arglist, 1, ierflg)
    #arglist[0] = 2
    #gMinuit.mnexcm("FIX", arglist, 1, ierflg)
    #arglist[0] = 3
    #gMinuit.mnexcm("FIX", arglist, 1, ierflg)
    #arglist[0] = 4
    #gMinuit.mnexcm("FIX", arglist, 1, ierflg)
    #arglist[0] = 5
    #gMinuit.mnexcm("FIX", arglist, 1, ierflg)

    # Scan for best parameter values
    #arglist[0] = 1
    #gMinuit.mnexcm("SCAN", arglist, 0, ierflg)

    # Now ready for minimization step
    arglist[0] = 500
    arglist[1] = 1.  # default: 0.1
    gMinuit.mnexcm("SIMPLEX", arglist, 2, ierflg)
    #gMinuit.mnexcm("MIGRAD", arglist, 2, ierflg)

    # Search for additional distinct local minima
    #arglist[0] = 10000 
    #gMinuit.mnexcm("IMP", arglist, 1, ierflg)

    # Do error analysis
    gMinuit.mnexcm("HESSE", arglist, 0, ierflg)
    #arglist[0] = 10000
    #gMinuit.mnexcm("MINOS", arglist, 1, ierflg)

    # Print results
    amin, edm, errdef = Double(0.18), Double(0.19), Double(0.20)  # passing doubles by reference is allowed, as long as on the python side doubles are unique (init them with different values, for example).
    nvpar, nparx, icstat = Long(1986), Long(1987), Long(1988)
    gMinuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat)
    gMinuit.mnprin(5, amin)

    ndf = 0
    for i in xrange(len(data[0])):  # loop over all bins
        mc_ = [m[i] for m in mc]
        if sum(mc_) < tolerance:
            continue
        if data[0][i] < tolerance:
            continue
        ndf += 1

    print "chi^2 = ", amin, ", NDF = ", ndf, ", chi^2/NDF = ", amin/ndf, ", prob = ", TMath.Prob(amin,ndf)
    print "amin: ", amin, ", edm: ", edm, ", errdef: ", errdef, ", nvpar: ", nvpar, ", nparx: ", nparx, ", icstat: ", icstat
    print "c0: ", vstart[0], ", c1: ", vstart[1], ", c2: ", vstart[2], ", c3: ", vstart[3], ", c4: ", vstart[4]


##______________________________________________________________________________
if __name__ == '__main__':
    infile = "results/ifit_HptReg.root"
    key = "HptReg_%s_%s"
    
    gROOT.SetBatch(True)
    TH1.AddDirectory(0)  # don't put TH1Fs in the TFile directory
    TH1.SetDefaultSumw2(1)

    # Fill histograms and join them
    fill(infile, key)
    
    # Fit
    ifit()

