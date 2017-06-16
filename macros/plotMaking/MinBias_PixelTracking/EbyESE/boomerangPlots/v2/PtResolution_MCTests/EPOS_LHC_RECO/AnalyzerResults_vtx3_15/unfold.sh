N="2"
ETA="PtResolution_MCTests/EPOS_LHC_RECO"
EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`

cd $EBYESE/unfoldingv4
cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults_vtx3_15/CastleEbyE.root data/PbPb_2015/data/CastleEbyE.root
cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults_vtx3_15/DDResp/dataDrivenResponseAndPriors.root DDResp/dataDrivenResponseAndPriors.root
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldDataBoomerang_DoSys.C++
UnfoldDataBoomerang_DoSys($N)
EOF
cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/UnfoldResults/dataResp/data${N}vtx3_15.root

echo "DONE!"