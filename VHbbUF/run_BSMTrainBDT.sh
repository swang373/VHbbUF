# sig: ZHinv, bkg: regular
root -l -b -q BSMTrainBDT.C+\(125,\"\",\"ZbbHinv\",\"\",\"weightsInvRegular\"\)

# sig: ZHinv, bkg: regular+VV
root -l -b -q BSMTrainBDT.C+\(125,\"allVV\",\"ZbbHinv\",\"\",\"weightsInvRegularVV\"\)

# sig: ZHinv, bkg: regular+VV+VHbb
root -l -b -q BSMTrainBDT.C+\(125,\"allVH\",\"ZbbHinv\",\"\",\"weightsInvRegularVH\"\)

# sig: ZHinv+VHbb, bkg: regular
root -l -b -q BSMTrainBDT.C+\(125,\"\",\"ZHZH\",\"\",\"weightsInvRegular__sigZHZH\"\)

# sig: ZHinv+ZZHF, bkg: regular
root -l -b -q BSMTrainBDT.C+\(125,\"\",\"ZbbHinvZZHF\",\"\",\"weightsInvRegular__sigZbbHinvZZHF\"\)

# sig: ZHinv, bkg: ZZ
root -l -b -q BSMTrainBDT.C+\(125,\"ZZ\",\"ZbbHinv\",\"\",\"weightsInvZZ\"\)

# sig: ZHinv, bkg: ZZHF
root -l -b -q BSMTrainBDT.C+\(125,\"ZZHF\",\"ZbbHinv\",\"\",\"weightsInvZZHF\"\)


## sig: ZHinv, bkg: regular (MET>200)
#root -l -b -q BSMTrainBDT.C+\(125,\"\",\"ZbbHinv\",\"useMET200\",\"weightsInvRegularMET200\"\)

## sig: ZHinv, bkg: regular (MET>220)
#root -l -b -q BSMTrainBDT.C+\(125,\"\",\"ZbbHinv\",\"useMET220\",\"weightsInvRegularMET220\"\)
