cd tuneBDT
ln -s ../HelperFunctions.h
ln -s ../HelperTMVA.h
ln -s ../tdrstyle.C
ln -s ../TrainBDT.C
ln -s ../TrimTree.C

ln -s ../datacards/125/HelperBDTShape.C
ln -s ../datacards/125/BDTShapeWorkspaceJ12.C
cp    ../datacards/125/BDTShapeJ12.C BDTShapeJ12_reload.C.orig
cp    ../datacards/125/BDTShapeJ12.C BDTShapeJ12_reload.C
