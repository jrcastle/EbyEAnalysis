N="2"
EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"

cd $EBYESE/unfoldingv3
cp $EBYESE/boomerangPlots/v$N/eta2.4/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/CastleEbyE.root
cp $EBYESE/boomerangPlots/v$N/eta2.4/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/dataDrivenResponseAndPriors.root
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldDataBoomerangGaussResp_JamesUnfold.C++
UnfoldDataBoomerangGaussResp_JamesUnfold($N)
EOF
cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta2.4/UnfoldResults/dataResp/data${N}Gauss.root

echo "DONE!"