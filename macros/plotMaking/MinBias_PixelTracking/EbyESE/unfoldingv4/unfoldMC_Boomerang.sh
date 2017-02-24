eta="eta1.0"
NORDER="2"
endDir="dataResp"
fName="data$NORDER.root"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`
cd -
echo "Job started on $(date)"
echo "cp /rfs/jcastle/PbPb2015/skewnessMC/CastleEbyE.root data/PbPb_2015/MC/./"
cp /rfs/jcastle/PbPb2015/skewnessMC/CastleEbyE.root data/PbPb_2015/MC/./
echo "cp /rfs/jcastle/PbPb2015/skewnessMC/dataDrivenResponseAndPriors.root DDResp/./"
cp /rfs/jcastle/PbPb2015/skewnessMC/dataDrivenResponseAndPriors.root DDResp/./
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldMCBoomerang_JamesUnfold.C++
UnfoldMCBoomerang_JamesUnfold($NORDER)
EOF
echo "txt/PbPb_2015/MC/data$NORDER.root /rfs/jcastle/PbPb2015/skewnessMC/$fName"
cp txt/PbPb_2015/MC/data$NORDER.root /rfs/jcastle/PbPb2015/skewnessMC/$fName
echo "Job finished on $(date)"