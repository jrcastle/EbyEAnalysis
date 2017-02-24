N="4"
ETA="eta2.4"
EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"
OUT="!!! PROGRESS:"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`

echo "$OUT Job started at $(date)"

### ##-- statErrorHandle/runAnalysis.sh
### echo "$OUT statErrorHandle/runAnalysis.sh (1/22)"
### cd $EBYESE/ESE/statErrorHandle/v$N/$ETA/AnalyzerResults
### bash runAnalysis.sh
### 
### ##-- statErrorHandle/Unfolding
### echo "$OUT statErrorHandle/Unfolding (2/22)"
### cd $EBYESE/unfoldingv$N
### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
### do
### echo "$OUT Processing split iteration $SPLIT"
### cp $EBYESE/ESE/statErrorHandle/v$N/$ETA/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
### cp $EBYESE/ESE/statErrorHandle/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
### root -l -b<<EOF
### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold"); 
### .L UnfoldDataEbyESE_JamesUnfold.C++
### UnfoldDataEbyESE_JamesUnfold($N)
### EOF
### cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/statErrorHandle/v$N/$ETA/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
### done
### 
### ##-- statErrorHandle/statUncertAssess.C
### echo "$OUT statErrorHandle/statUncertAssess.C (3/22)"
### cd $EBYESE/ESE/statErrorHandle/v$N/$ETA
### root -l -b<<EOF
### .x statUncertAssess.C
### EOF
### 
##-- v$N/runAnalysis.sh
echo "$OUT v$N/runAnalysis.sh (4/22)"
cd $EBYESE/ESE/v$N/$ETA/AnalyzerResults
bash runAnalysis.sh

##-- v$N/Unfolding
echo "$OUT v$N/Unfolding (5/22)"
cd $EBYESE/unfoldingv$N
cp $EBYESE/ESE/v$N/$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
cp $EBYESE/ESE/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldDataEbyESE_JamesUnfold.C++
UnfoldDataEbyESE_JamesUnfold($N)
EOF
cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/v$N/$ETA/UnfoldResults/dataResp/./

##-- v$N/unfoldAssess.C
echo "$OUT v$N/unfoldAssess.C (6/22)"
cd $EBYESE/ESE/v$N/$ETA
root -l -b<<EOF
.x unfoldAssess.C
EOF

### ##-- v$N/systematics/vtx_leq_3/runAnalysis.sh
### echo "$OUT v$N/systematics/vtx_leq_3/runAnalysis.sh (7/22)"
### cd $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults
### bash runAnalysis.sh
### 
### ##-- v$N/systematics/vtx_leq_3/Unfolding
### echo "$OUT v$N/systematics/vtx_leq_3/Unfolding (8/22)"
### cd $EBYESE/unfoldingv$N
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
### root -l -b<<EOF
### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
### .L UnfoldDataEbyESE_JamesUnfold.C++
### UnfoldDataEbyESE_JamesUnfold($N)
### EOF
### cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/UnfoldResults/dataResp/./
### 
### ##-- v$N/systematics/vtx3_15/runAnalysis.sh
### echo "$OUT v$N/systematics/vtx3_15/runAnalysis (9/22)"
### cd $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults
### bash runAnalysis.sh
### 
### ##-- v$N/systematics/vtx3_15/Unfolding
### echo "$OUT v$N/systematics/vtx3_15/Unfolding (10/22)"
### cd $EBYESE/unfoldingv$N
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
### root -l -b<<EOF
### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
### .L UnfoldDataEbyESE_JamesUnfold.C++
### UnfoldDataEbyESE_JamesUnfold($N)
### EOF
### cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/UnfoldResults/dataResp/./
### 
### ##-- v$N/systematics/sysVtxCut.C
echo "$OUT v$N/systematics/sysVtxCut.C (11/22)"
cd $EBYESE/ESE/v$N/$ETA/systematicStudies/vtxCut
root -l -b<<EOF
.x sysVtxCut.C
EOF
### 
### ##-- v$N/systematics/tkEff/runAnalysis.sh
### echo "$OUT v$N/systematics/tkEff/runAnalysis.sh (12/22)"
### cd $EBYESE/ESE/v$N/$ETA/systematicStudies/tkEff/AnalyzerResults
### bash runAnalysis.sh
### 
### ##-- v$N/systematics/tkEff/Unfolding
### echo "$OUT v$N/systematics/tkEff/Unfolding (13/22)"
### cd $EBYESE/unfoldingv$N
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/tkEff/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/tkEff/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
### root -l -b<<EOF
### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
### .L UnfoldDataEbyESE_JamesUnfold.C++
### UnfoldDataEbyESE_JamesUnfold($N)
### EOF
### cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/v$N/$ETA/systematicStudies/tkEff/UnfoldResults/dataResp/./
### 
##-- v$N/systematics/tkEff/sysTkEff.C
echo "$OUT v$N/systematics/tkEff/sysTkEff.C (14/22)"
cd $EBYESE/ESE/v$N/$ETA/systematicStudies/tkEff
root -l -b<<EOF
.x sysTkEff.C
EOF

### ##-- v$N/systematics/tkQuality/loose/runAnalysis.sh
### echo "$OUT v$N/systematics/tkQuality/loose/runAnalysis.sh (15/22)"
### cd $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults
### bash runAnalysis.sh
### 
### ##-- v$N/systematics/tkQuality/loose/Unfolding 
### echo "$OUT v$N/systematics/tkQuality/loose/Unfolding (16/22)"
### cd $EBYESE/unfoldingv$N
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
### root -l -b<<EOF
### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
### .L UnfoldDataEbyESE_JamesUnfold.C++
### UnfoldDataEbyESE_JamesUnfold($N)
### EOF
### cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/loose/UnfoldResults/dataResp/./
### 
### ##-- v$N/systematics/tkQuality/tight/runAnalysis.sh
### echo "$OUT v$N/systematics/tkQuality/tight/runAnalysis.sh (17/22)"
### cd $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults
### bash runAnalysis.sh
### 
### ##-- v$N/systematics/tkQuality/tight/Unfolding
### echo "$OUT v$N/systematics/tkQuality/tight/Unfolding (18/22)"
### cd $EBYESE/unfoldingv$N
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
### cp $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
### root -l -b<<EOF
### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
### .L UnfoldDataEbyESE_JamesUnfold.C++
### UnfoldDataEbyESE_JamesUnfold($N)
### EOF
### cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality/tight/UnfoldResults/dataResp/./
### 
##-- v$N/systematics/tkQuality/sysTkQuality.C
echo "$OUT v$N/systematics/tkQuality/sysTkQuality.C (19/22)"
cd $EBYESE/ESE/v$N/$ETA/systematicStudies/tkQuality
root -l -b<<EOF
.x sysTkQuality.C
EOF

##-- v$N/systematics/sysChi2Cutoff.C
echo "$OUT systematics/sysChi2Cutoff.C (20/22)"
cd $EBYESE/ESE/v$N/$ETA/systematicStudies/chi2Cutoff
root -l -b<<EOF
.x sysChi2Cutoff.C
EOF

### ##-- v$N/Unfolding_dosys
### echo "$OUT v$N/Unfolding_dosys (21/22)"
### cd $EBYESE/unfoldingv$N
### cp $EBYESE/ESE/v$N/$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
### cp $EBYESE/ESE/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
### root -l -b<<EOF
### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
### .L UnfoldDataEbyESE_DoSys.C++
### UnfoldDataEbyESE_DoSys($N)
### EOF
### cp txt/PbPb_2015/data/data$N.root $EBYESE/ESE/v$N/$ETA/UnfoldResults/dataResp/data${N}_dosys.root
### 
### ##-- v$N/systematics/sysRespEl.C
### echo "$OUT v$N/systematics/sysRespEl.C (22/22)"
### cd $EBYESE/ESE/v$N/$ETA/systematicStudies/responseElements
### root -l -b<<EOF
### .x sysRespEl.C
### EOF
### 
### echo "$OUT Job finished at $(date)"