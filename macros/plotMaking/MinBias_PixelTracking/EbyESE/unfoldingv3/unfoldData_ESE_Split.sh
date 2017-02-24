eta="eta2.4"
NORDER="2"
endDir="dataResp"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`
cd -
echo "Job started on $(date)"

for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
do
echo "Processing split iteration $SPLIT"
echo "cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/statErrHandle/v$NORDER/$eta/AnalyzerResults/CastleEbyE_Split$SPLIT.root data/PbPb_2015/data/CastleEbyE.root"
cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/statErrHandle/v$NORDER/$eta/AnalyzerResults/CastleEbyE_Split$SPLIT.root data/PbPb_2015/data/CastleEbyE.root
echo "cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/statErrHandle/v$NORDER/$eta/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split$SPLIT.root DDResp/dataDrivenResponseAndPriors.root"
cp /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/statErrHandle/v$NORDER/$eta/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split$SPLIT.root DDResp/dataDrivenResponseAndPriors.root
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldDataEbyESE_JamesUnfold.C++
UnfoldDataEbyESE_JamesUnfold($NORDER) 
EOF
echo "cp txt/PbPb_2015/data/data$NORDER.root /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/statErrHandle/v$NORDER/$eta/UnfoldResults/$endDir/data${NORDER}_Split${SPLIT}.root"
cp txt/PbPb_2015/data/data$NORDER.root /home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE/ESE/statErrHandle/v$NORDER/$eta/UnfoldResults/$endDir/data${NORDER}_Split${SPLIT}.root
done
echo "Job finished on $(date)"
