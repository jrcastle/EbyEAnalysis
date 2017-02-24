N="2"
ETA="eta1.0_WideBin_TruncPriorResp2"
EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"

cd $EBYESE/unfoldingv2
cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/CastleEbyE.root
cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/dataDrivenResponseAndPriors.root
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldDataBoomerang_SVD_Trunc2.C++
UnfoldDataBoomerang_SVD_Trunc2($N)
EOF
cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/UnfoldResults/dataResp/data${N}_svd.root

echo "DONE!"