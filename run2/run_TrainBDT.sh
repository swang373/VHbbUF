root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"\",\"weightsRegular\"\) >&! train0_0.log
root -l -b -q TrainBDT.C+\(125,\"TT\",\"\",\"\",\"weightsTT\"\) >&! train0_1.log
root -l -b -q TrainBDT.C+\(125,\"Vj\",\"\",\"\",\"weightsVj\"\) >&! train0_2.log
root -l -b -q TrainBDT.C+\(125,\"VjLF\",\"\",\"\",\"weightsVjLF\"\) >&! train0_3.log
root -l -b -q TrainBDT.C+\(125,\"ZjLF\",\"\",\"\",\"weightsZjLF\"\) >&! train0_4.log
root -l -b -q TrainBDT.C+\(125,\"VjHF\",\"\",\"\",\"weightsVjHF\"\) >&! train0_5.log
root -l -b -q TrainBDT.C+\(125,\"ZjHF\",\"\",\"\",\"weightsZjHF\"\) >&! train0_6.log
root -l -b -q TrainBDT.C+\(125,\"ST\",\"\",\"\",\"weightsST\"\) >&! train0_7.log
root -l -b -q TrainBDT.C+\(125,\"VV\",\"\",\"\",\"weightsVV\"\) >&! train0_8.log
root -l -b -q TrainBDT.C+\(125,\"ZZ\",\"\",\"\",\"weightsZZ\"\) >&! train0_9.log
root -l -b -q TrainBDT.C+\(125,\"ZZHF\",\"\",\"\",\"weightsZZHF\"\) >&! train0_10.log

root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"FJ\",\"weightsRegular_FJ\"\) >&! train1_0.log
root -l -b -q TrainBDT.C+\(125,\"TT\",\"\",\"FJ\",\"weightsTT_FJ\"\) >&! train1_1.log
root -l -b -q TrainBDT.C+\(125,\"VjLF\",\"\",\"FJ\",\"weightsVjLF_FJ\"\) >&! train1_2.log
root -l -b -q TrainBDT.C+\(125,\"ZZ\",\"\",\"FJ\",\"weightsZZ_FJ\"\) >&! train1_3.log
root -l -b -q TrainBDT.C+\(125,\"ZZHF\",\"\",\"FJ\",\"weightsZZHF_FJ\"\) >&! train1_4.log

root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"FJReg\",\"weightsRegular_FJReg\"\) >&! train2_0.log
root -l -b -q TrainBDT.C+\(125,\"TT\",\"\",\"FJReg\",\"weightsTT_FJReg\"\) >&! train2_1.log
root -l -b -q TrainBDT.C+\(125,\"VjLF\",\"\",\"FJReg\",\"weightsVjLF_FJReg\"\) >&! train2_2.log
root -l -b -q TrainBDT.C+\(125,\"ZZ\",\"\",\"FJReg\",\"weightsZZ_FJReg\"\) >&! train2_3.log
root -l -b -q TrainBDT.C+\(125,\"ZZHF\",\"\",\"FJReg\",\"weightsZZHF_FJReg\"\) >&! train2_4.log

root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"Angle\",\"weightsRegular_Angle\"\) >&! train3_0.log

root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"KillTop\",\"weightsRegular_KillTop\"\) >&! train4_0.log

root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"MET\",\"weightsRegular_MET\"\) >&! train5_0.log

root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"NoMjj\",\"weightsRegular_NoMjj\"\) >&! train6_0.log

root -l -b -q TrainBDT.C+\(125,\"\",\"\",\"Super\",\"weightsSuper\"\) >&! train9_0.log

root -l -b -q TrainBDT.C+\(125,\"\",\"ZZHF\",\"\",\"weightsRegular__sigZZHF\"\) >&! train11_0.log
root -l -b -q TrainBDT.C+\(125,\"TT\",\"ZZHF\",\"\",\"weightsTT__sigZZHF\"\) >&! train11_1.log
root -l -b -q TrainBDT.C+\(125,\"VjLF\",\"ZZHF\",\"\",\"weightsVjLF__sigZZHF\"\) >&! train11_2.log
root -l -b -q TrainBDT.C+\(125,\"VV\",\"ZZHF\",\"\",\"weightsVV__sigZZHF\"\) >&! train11_3.log

#root -l -b -q TrainBDT_splitsig.C+\(125,\"\",\"\",\"\",\"weightsRegularSplitSig\"\) >&! train21_0.log
#root -l -b -q TrainBDT_splitsig.C+\(125,\"TT\",\"\",\"\",\"weightsTTSplitSig\"\) >&! train21_1.log
#root -l -b -q TrainBDT_splitsig.C+\(125,\"VjLF\",\"\",\"\",\"weightsVjLFSplitSig\"\) >&! train21_2.log
#root -l -b -q TrainBDT_splitsig.C+\(125,\"ZZ\",\"\",\"\",\"weightsZZSplitSig\"\) >&! train21_3.log
#root -l -b -q TrainBDT_splitsig.C+\(125,\"ZZHF\",\"\",\"\",\"weightsZZHFSplitSig\"\) >&! train21_4.log

#root -l -b -q TrainBDT_fullbkg.C+\(125,\"\",\"\",\"\",\"weightsRegularFullBkg\"\) >&! train31_0.log

#root -l -b -q TrainBDT_hybrid.C+\(125,\"\",\"\",\"\",\"weightsRegularHybrid\"\) >&! train41_0.log
#root -l -b -q TrainBDT_hybrid.C+\(125,\"\",\"\",\"FJ\",\"weightsRegularHybrid_FJ\"\) >&! train42_0.log
#root -l -b -q TrainBDT_hybrid.C+\(125,\"\",\"\",\"FJReg\",\"weightsRegularHybrid_FJReg\"\) >&! train43_0.log
