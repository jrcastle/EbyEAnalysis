N="3"
ETA="eta1.0"
EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"
OUT="!!! PROGRESS:"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`

echo "$OUT Job started at $(date)"

#### ##-- statErrorHandle/runAnalysis.sh
#### echo "$OUT statErrorHandle/runAnalysis.sh (1/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- statErrorHandle/Unfolding
#### echo "$OUT statErrorHandle/Unfolding (2/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold"); 
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/statUncertAssess.C
#### echo "$OUT statErrorHandle/statUncertAssess.C (3/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA
#### root -l -b<<EOF
#### .x statUncertAssess.C
#### EOF
#### 
#### ##-- statErrorHandle/Unfolding_Gauss
#### echo "$OUT statErrorHandle/Unfolding_Gauss (4/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerangGaussResp_JamesUnfold.C++
#### UnfoldDataBoomerangGaussResp_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/UnfoldResults/dataResp/data${N}Gauss_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/systematics/clusCompatTune/runAnalysis.sh
#### echo "$OUT statErrorHandle/systematics/clusCompatTune/newCCTune2pct/runAnalysis.sh (5/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults
#### bash runAnalysis.sh
####  
#### ##-- statErrorHandle/systematics/clusCompatTune/Unfolding
#### echo "$OUT statErrorHandle/systematics/clusCompatTune/newCCTune2pct/Unfolding (6/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/systematics/tkQuality/loose/runAnalysis.sh
#### echo "$OUT statErrorHandle/systematics/tkQuality/loose/runAnalysis.sh (7/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- statErrorHandle/systematics/tkQuality/loose/Unfolding
#### echo "$OUT statErrorHandle/systematics/tkQuality/loose/Unfolding (8/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/loose/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/systematics/tkQuality/tight/runAnalysis.sh
#### echo "$OUT statErrorHandle/systematics/tkQuality/tight/runAnalysis.sh (9/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- statErrorHandle/systematics/tkQuality/tight/Unfolding 
#### echo "$OUT statErrorHandle/systematics/tkQuality/tight/Unfolding (10/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/tight/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/systematics/vtxCut/vtx3_15/runAnalysis.sh
#### echo "$OUT statErrorHandle/systematics/vtxCut/vtx3_15/runAnalysis.sh (11/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- statErrorHandle/systematics/vtxCut/vtx3_15/Unfolding 
#### echo "$OUT statErrorHandle/systematics/vtxCut/vtx3_15/Unfolding (12/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/systematics/vtxCut/vtx_leq_3/runAnalysis.sh
#### echo "$OUT statErrorHandle/systematics/vtxCut/vtx_leq_3/runAnalysis.sh (13/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- statErrorHandle/systematics/vtxCut/vtx_leq_3/Unfolding 
#### echo "$OUT statErrorHandle/systematics/vtxCut/vtx_leq_3/Unfolding (14/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++ 
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/statUncertChi2Cutoff.C
#### echo "$OUT statErrorHandle/statUncertChi2Cutoff.C (15/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/chi2Cutoff
#### root -l -b<<EOF
#### .x statChi2Cutoff.C
#### EOF
#### 
#### ##-- statErrorHandle/statUncertRespEl.C
#### echo "$OUT statErrorHandle/statUncertRespEl.C (16/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/responseElements
#### root -l -b<<EOF
#### .x statRespEl.C
#### EOF
#### 
#### ##-- statErrorHandle/statUncertNewCC.C
#### echo "$OUT statErrorHandle/statUncertNewCC.C (17/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct
#### root -l -b<<EOF
#### .x statUncertAssess.C
#### EOF
#### 
#### ##-- statErrorHandle/statUncertTkQualityLoose.C
#### echo "$OUT statErrorHandle/statUncertTkQualityLoose.C (18/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/loose
#### root -l -b<<EOF
#### .x statUncertTkQualityLoose.C
#### EOF
#### 
#### ##-- statErrorHandle/statUncertTkQualityTight.C
#### echo "$OUT statErrorHandle/statUncertTkQualityTight.C (19/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/tkQuality/tight
#### root -l -b<<EOF
#### .x statUncertTkQualityTight.C
#### EOF
#### 
#### ##-- statErrorHandle/statUncertVtx3_15.C
#### echo "$OUT statErrorHandle/statUncertVtx3_15.C (20/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/
#### root -l -b<<EOF
#### .x statUncertVtx3_15.C
#### EOF
#### 
#### ##-- statErrorHandle/statUncertVtx_leq_3.C
#### echo "$OUT statErrorHandle/statUncertVtx_leq_3.C (21/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/
#### root -l -b<<EOF
#### .x statUncertVtx_leq_3.C
#### EOF
#### 
#### ####################################
#### 
#### ##-- v$N/runAnalysis.sh
#### echo "$OUT v$N/runAnalysis.sh (22/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- v$N/Unfolding
#### echo "$OUT v$N/Unfolding (23/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/UnfoldResults/dataResp/./
#### 
#### ##-- v$N/unfoldAssess.C
#### echo "$OUT v$N/unfoldAssess.C (24/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA
#### root -l -b<<EOF
#### .x unfoldAssess.C
#### EOF
#### 
#### ##-- v$N/systematics/vtx_leq_3/runAnalysis.sh
#### echo "$OUT v$N/systematics/vtx_leq_3/runAnalysis.sh (25/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- v$N/systematics/vtx_leq_3/Unfolding
#### echo "$OUT v$N/systematics/vtx_leq_3/Unfolding (26/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx_leq_3/UnfoldResults/dataResp/./
#### 
#### ##-- v$N/systematics/vtx3_15/runAnalysis.sh
#### echo "$OUT v$N/systematics/vtx3_15/runAnalysis.sh  (27/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- v$N/systematics/vtx3_15/Unfolding
#### echo "$OUT v$N/systematics/vtx3_15/Unfolding (28/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut/vtx3_15/UnfoldResults/dataResp/./
#### 
#### ##-- v$N/systematics/sysVtxCut.C
#### echo "$OUT v$N/systematics/sysVtxCut.C (29/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/vtxCut
#### root -l -b<<EOF
#### .x sysVtxCut.C
#### EOF
#### 
#### ##-- v$N/systematics/clusCompatTune/runAnalysis.sh
#### echo "$OUT v$N/systematics/clusCompatTune/newCCTune2pct/runAnalysis.sh (30/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- v$N/systematics/clusCompatTune/Unfolding
#### echo "$OUT v$N/systematics/clusCompatTune/newCCTune2pct/Unfolding (31/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/clusCompatTune/newCCTune2pct/UnfoldResults/dataResp/./
#### 
#### ##-- v$N/systematics/clusCompatTune/sysNewCC.C
#### echo "$OUT v$N/systematics/clusCompatTune/newCCTune2pct/sysClusCompatTune.C (32/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/clusCompatTune
#### root -l -b<<EOF
#### .x sysClusCompatTune.C
#### EOF
#### 
#### ##-- v$N/systematics/tkQuality/loose/runAnalysis.sh
#### echo "$OUT v$N/systematics/tkQuality/loose/runAnalysis.sh (33/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- v$N/systematics/tkQuality/loose/Unfolding 
#### echo "$OUT v$N/systematics/tkQuality/loose/Unfolding (34/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/loose/UnfoldResults/dataResp/./
#### 
#### ##-- v$N/systematics/tkQuality/tight/runAnalysis.sh
#### echo "$OUT v$N/systematics/tkQuality/tight/runAnalysis.sh (35/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults
#### bash runAnalysis.sh
#### 
#### ##-- v$N/systematics/tkQuality/tight/Unfolding
#### echo "$OUT v$N/systematics/tkQuality/tight/Unfolding (36/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_JamesUnfold.C++
#### UnfoldDataBoomerang_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality/tight/UnfoldResults/dataResp/./
#### 
#### ##-- v$N/systematics/tkQuality/sysTkQuality.C
#### echo "$OUT v$N/systematics/tkQuality/sysTkQuality.C (37/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/tkQuality
#### root -l -b<<EOF
#### .x sysTkQuality.C
#### EOF
#### 
#### ##-- v$N/systematics/sysChi2Cutoff.C
#### echo "$OUT systematics/sysChi2Cutoff.C (38/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/chi2Cutoff
#### root -l -b<<EOF
#### .x sysChi2Cutoff.C
#### EOF
#### 
#### ##-- v$N/Unfolding_Gauss
#### echo "$OUT v$N/Unfolding_Gauss (39/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerangGaussResp_JamesUnfold.C++
#### UnfoldDataBoomerangGaussResp_JamesUnfold($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/UnfoldResults/dataResp/data${N}Gauss.root

#### DOSYS
##-- v$N/Unfolding_dosys
echo "$OUT v$N/Unfolding_dosys (40/44)"
cd $EBYESE/unfoldingv$N
cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
root -l -b<<EOF
gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
.L UnfoldDataBoomerang_DoSys.C++
UnfoldDataBoomerang_DoSys($N)
EOF
cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/UnfoldResults/dataResp/data${N}_dosys.root
#### 
#### ##-- v$N/systematics/sysRespEl.C
#### echo "$OUT v$N/systematics/sysRespEl.C (41/44)"
#### cd $EBYESE/boomerangPlots/v$N/$ETA/systematicStudies/responseElements
#### root -l -b<<EOF
#### .x sysRespEl.C
#### EOF
#### 
#### ##-- statErrorHandle/SVDUnfold
#### echo "$OUT statErrorHandle/SVDUnfolding (42/44)"
#### cd $EBYESE/unfoldingv$N
#### for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
#### do
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
#### cp $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_SVD.C++
#### UnfoldDataBoomerang_SVD($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA/UnfoldResults/dataResp/data${N}_svd_Split${SPLIT}.root
#### done
#### 
#### ##-- statErrorHandle/SVDUnfold/statuncertassess 
#### echo "$OUT statErrorHandle/statUncertAssessSVD.C (43/44)"
#### cd $EBYESE/boomerangPlots/statErrorHandle/v$N/$ETA
#### root -l -b<<EOF
#### .x statUncertAssessSVD.C
#### EOF
#### 
#### ##-- v$N/SVDUnfolding
#### echo "$OUT v$N/SVDUnfolding (44/44)"
#### cd $EBYESE/unfoldingv$N
#### cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
#### cp $EBYESE/boomerangPlots/v$N/$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
#### root -l -b<<EOF
#### gSystem->Load("/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/libRooUnfold");
#### .L UnfoldDataBoomerang_SVD.C++
#### UnfoldDataBoomerang_SVD($N)
#### EOF
#### cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/$ETA/UnfoldResults/dataResp/data${N}_svd.root

echo "$OUT Job finished at $(date)"