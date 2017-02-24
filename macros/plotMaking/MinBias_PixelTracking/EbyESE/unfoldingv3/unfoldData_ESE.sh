eta="eta2.4"
NORDER="2"
endDir="dataResp"
fName="data$NORDER.root"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`
cd -

echo "Job started on $(date)"
echo "cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/v$NORDER/$eta/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./"
cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/v$NORDER/$eta/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
echo "cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/v$NORDER/$eta/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./"
cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/v$NORDER/$eta/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldDataEbyESE_JamesUnfold.C++
UnfoldDataEbyESE_JamesUnfold($NORDER)
EOF
echo "cp txt/PbPb_2015/data/data$NORDER.root /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/v$NORDER/$eta/UnfoldResults/$endDir/$fName"
cp txt/PbPb_2015/data/data$NORDER.root /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/v$NORDER/$eta/UnfoldResults/$endDir/$fName
echo "Job finished on $(date)"