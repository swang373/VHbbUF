#!/bin/py


class ParseJSON():
	def __init__(self,jfile):
		print1 = "Opening up json file: %s"%jfile
		self.infile = open(jfile,"r")
		print print1


		self.runs= []
		self.getJsonList(self.infile)	
		self.printVersion()


	def getJsonList(self,infile):	
		json = infile.readlines()[0]
		lines = json.rsplit("{")[1].rsplit("}")[0].rsplit("%s"%'"')
	#	print len(lines)
		for line in lines: 
			test = line.rsplit(":")
			if len(test) == 2: 
				lumis = test[1].rsplit(",")
				eve = -1
				odd = -1
				for i in range(0,len(lumis)):
					lumi = lumis[i]
					aa = lumi.lstrip().lstrip("[").rstrip("]")
					if aa == "": continue
					if i%2 == 0: eve = int(aa)
					else: odd = int(aa)
					if eve >=0 and odd >= 0:
						self.runs.append((theRun,eve,odd))
						eve = -1
						odd = -1	
			elif len(test[0]) ==0:
				continue
			else:
				run = int(test[0])
				theRun = run
	

	def dumpRuns(self): 
		print "dumping out the run and lumi list..."
		for shit in self.runs: print "%d: %d-%d"%(shit[0],shit[1],shit[2])
	

	def isValidRunLumi(self,run,lumi):
		for shit in self.runs:
			if run != shit[0]: continue	
			if lumi < shit[1]: continue
			if lumi > shit[2]: continue
			return True			
		return False

	def printVersion(self):
		self.version =\
		'''
		JSON Parser v. 1.0
		T.N.K October 2010
		Congratulations on being sick of running over your data over and
		over again with JSON files because elevnteen people decided
		that it should be slightly changed. You've taken the first step
		into the rest of your amazing evening where you can easily manipulate
		your plots and pretend that it took you hours upon hours. 
		'''
		print self.version


	def getCMSSWFormat(self):
		retstring = ""
		retstring +="%s\n"%("import FWCore.ParameterSet.Config as cms")
		retstring +="%s\n"%("jsonlist = cms.untracked.VLuminosityBlockRange(")

		for shit in self.runs:
			retstring += "\'{0}:{1}-{0}:{2}\',\n".format(shit[0],shit[1],shit[2])
		retstring +=")\n"

		return retstring


#!/bin/py


class ParseJSON():
	def __init__(self,jfile):
		print1 = "Opening up json file: %s"%jfile
		self.infile = open(jfile,"r")
		print print1


		self.runs= []
		self.getJsonList(self.infile)	
		self.printVersion()


	def getJsonList(self,infile):	
		json = infile.readlines()[0]
		lines = json.rsplit("{")[1].rsplit("}")[0].rsplit("%s"%'"')
	#	print len(lines)
		for line in lines: 
			test = line.rsplit(":")
			if len(test) == 2: 
				lumis = test[1].rsplit(",")
				eve = -1
				odd = -1
				for i in range(0,len(lumis)):
					lumi = lumis[i]
					aa = lumi.lstrip().lstrip("[").rstrip("]")
					if aa == "": continue
					if i%2 == 0: eve = int(aa)
					else: odd = int(aa)
					if eve >=0 and odd >= 0:
						self.runs.append((theRun,eve,odd))
						eve = -1
						odd = -1	
			elif len(test[0]) ==0:
				continue
			else:
				run = int(test[0])
				theRun = run
	

	def dumpRuns(self): 
		print "dumping out the run and lumi list..."
		for shit in self.runs: print "%d: %d-%d"%(shit[0],shit[1],shit[2])
	

	def isValidRunLumi(self,run,lumi):
		for shit in self.runs:
			if run != shit[0]: continue	
			if lumi < shit[1]: continue
			if lumi > shit[2]: continue
			return True			
		return False

	def printVersion(self):
		self.version =\
		'''
		JSON Parser v. 1.0
		T.N.K October 2010
		Congratulations on being sick of running over your data over and
		over again with JSON files because elevnteen people decided
		that it should be slightly changed. You've taken the first step
		into the rest of your amazing evening where you can easily manipulate
		your plots and pretend that it took you hours upon hours. 
		'''
		print self.version

        def getCMSSWFormat(self):
		retstring = ""
		retstring +="%s\n"%("import FWCore.ParameterSet.Config as cms")
		retstring +="%s\n"%("jsonlist = cms.untracked.VLuminosityBlockRange(")

		for shit in self.runs:
			retstring += "\'{0}:{1}-{0}:{2}\',\n".format(shit[0],shit[1],shit[2])
		retstring +=")\n"

		return retstring


if __name__ == "__main__":
#	pj = ParseJSON("Cert_146645-99999_7TeV_StreamExpress_Collisions10_JSON_OnlyDCS.txt")
	pj = ParseJSON("2011_v2_JSON.txt")
	pj.dumpRuns()
	print pj.isValidRunLumi(147217,1)
	print pj.isValidRunLumi(147216,1)
	print pj.isValidRunLumi(147207,1)
        st = pj.getCMSSWFormat()
        print st
 
