#!/bin/bash
N="3"     # Flow Order
ETA="2.4" # Track eta max, |eta| < $ETA
PMN="0.3" # Track pT min
PMX="3.0" # Track pT max
TEST="0"  # Test run flag

EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"
ROOUNFOLD="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold"
OUT="!!! PROGRESS:"

FIRSTPASS="0"
################### ATTENTION! ####################
# Do not set FIRSTPASS="0" until you first set    #
# _dosys = 1 in $ROOUNFOLD/src/RooUnfoldBayes.cxx #
# and then recompile:                             #
# cd $ROOUNFOLD                                   #
# make                                            #
###################################################

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`

echo "$OUT Job started at $(date)"

##-- Check to see if statistical ananlysis directory already exists.  If not, copy from TEMPLATE
if [ ! -d "$EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA" ]
then 
    echo "Statistical analysis statErrorHandle/v$N/eta$ETA directory does not exist!"
    echo "Copying from template now...."
    mkdir $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA
    cp -r $EBYESE/boomerangPlots/TEMPLATE/statErrorResampling/* $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/./
fi

##-- Check to see if analysis directory already exists.  If not, copy from TEMPLATE 
if [ ! -d "$EBYESE/boomerangPlots/v$N/eta$ETA" ]
then
    echo "Analysis v$N/eta$ETA directory does not exist!"
    echo "Copying from template now...."
    mkdir $EBYESE/boomerangPlots/v$N/eta$ETA
    cp -r $EBYESE/boomerangPlots/TEMPLATE/analysis/* $EBYESE/boomerangPlots/v$N/eta$ETA/./
fi

##-- First pass of data, does a majority of the work.
if [ $FIRSTPASS == 1 ]
then

    ##-- ================================== statErrorHandle/runAnalysis.sh ==================================
    echo "$OUT statErrorHandle/runAnalysis.sh (1/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== statErrorHandle/Unfolding ==================================
    echo "$OUT statErrorHandle/Unfolding (2/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/statUncertAssess.C ==================================
    echo "$OUT statErrorHandle/statUncertAssess.C (3/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA
    root -l -b<<EOF
.L statUncertAssess.C 
statUncertAssess($N);
EOF

    ##-- ================================== statErrorHandle/Unfolding_Gauss ==================================
    echo "$OUT statErrorHandle/Unfolding_Gauss (4/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerangGaussResp_RooUnfold.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/UnfoldResults/dataResp/data${N}Gauss_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/systematics/clusCompatTune/runAnalysis.sh ==================================
    echo "$OUT statErrorHandle/systematics/clusCompatTune/newCCTune2pct/runAnalysis.sh (5/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST
 
    ##-- ================================== statErrorHandle/systematics/clusCompatTune/Unfolding ==================================
    echo "$OUT statErrorHandle/systematics/clusCompatTune/newCCTune2pct/Unfolding (6/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/systematics/tkQuality/loose/runAnalysis.sh ==================================
    echo "$OUT statErrorHandle/systematics/tkQuality/loose/runAnalysis.sh (7/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/loose/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== statErrorHandle/systematics/tkQuality/loose/Unfolding ==================================
    echo "$OUT statErrorHandle/systematics/tkQuality/loose/Unfolding (8/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/loose/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/systematics/tkQuality/tight/runAnalysis.sh ==================================
    echo "$OUT statErrorHandle/systematics/tkQuality/tight/runAnalysis.sh (9/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/tight/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== statErrorHandle/systematics/tkQuality/tight/Unfolding ==================================
    echo "$OUT statErrorHandle/systematics/tkQuality/tight/Unfolding (10/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/tight/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/systematics/vtxCut/vtx3_15/runAnalysis.sh ==================================
    echo "$OUT statErrorHandle/systematics/vtxCut/vtx3_15/runAnalysis.sh (11/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== statErrorHandle/systematics/vtxCut/vtx3_15/Unfolding ================================== 
    echo "$OUT statErrorHandle/systematics/vtxCut/vtx3_15/Unfolding (12/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/systematics/vtxCut/vtx_leq_3/runAnalysis.sh ==================================
    echo "$OUT statErrorHandle/systematics/vtxCut/vtx_leq_3/runAnalysis.sh (13/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== statErrorHandle/systematics/vtxCut/vtx_leq_3/Unfolding ==================================
    echo "$OUT statErrorHandle/systematics/vtxCut/vtx_leq_3/Unfolding (14/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/statUncertChi2Cutoff.C ==================================
    echo "$OUT statErrorHandle/statUncertChi2Cutoff.C (15/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/chi2Cutoff
    root -l -b <<EOF
.L statChi2Cutoff.C
statChi2Cutoff($N)
EOF

    ##-- ================================== statErrorHandle/statUncertRespEl.C ==================================
    echo "$OUT statErrorHandle/statUncertRespEl.C (16/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/responseElements
    root -l -b<<EOF 
.L statRespEl.C
statRespEl($N)
EOF

    ##-- ================================== statErrorHandle/statUncertNewCC.C ==================================
    echo "$OUT statErrorHandle/statUncertNewCC.C (17/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct
    root -l -b<<EOF 
.L statUncertAssess.C
statUncertAssess($N)
EOF

    ##-- ================================== statErrorHandle/statUncertTkQualityLoose.C ==================================
    echo "$OUT statErrorHandle/statUncertTkQualityLoose.C (18/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/loose
    root -l -b<<EOF
.L statUncertTkQualityLoose.C
statUncertTkQualityLoose($N)
EOF

    ##-- ================================== statErrorHandle/statUncertTkQualityTight.C ==================================
    echo "$OUT statErrorHandle/statUncertTkQualityTight.C (19/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/tkQuality/tight
    root -l -b<<EOF 
.L statUncertTkQualityTight.C
statUncertTkQualityTight($N)
EOF

    ##-- ================================== statErrorHandle/statUncertVtx3_15.C ==================================
    echo "$OUT statErrorHandle/statUncertVtx3_15.C (20/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/
    root -l -b<<EOF
.L statUncertVtx3_15.C
statUncertVtx3_15($N)
EOF

    ##-- ================================== statErrorHandle/statUncertVtx_leq_3.C ==================================
    echo "$OUT statErrorHandle/statUncertVtx_leq_3.C (21/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/
    root -l -b<<EOF
.L statUncertVtx_leq_3.C
statUncertVtx_leq_3($N)
EOF

    ##-- ================================== v$N/runAnalysis.sh ==================================
    echo "$OUT v$N/runAnalysis.sh (22/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== v$N/Unfolding ==================================
    echo "$OUT v$N/Unfolding (23/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/UnfoldResults/dataResp/./

    ##-- ================================== v$N/systematics/vtx_leq_3/runAnalysis.sh ==================================
    echo "$OUT v$N/systematics/vtx_leq_3/runAnalysis.sh (24/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== v$N/systematics/vtx_leq_3/Unfolding ==================================
    echo "$OUT v$N/systematics/vtx_leq_3/Unfolding (25/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx_leq_3/UnfoldResults/dataResp/./

    ##-- ================================== v$N/systematics/vtx3_15/runAnalysis.sh ==================================
    echo "$OUT v$N/systematics/vtx3_15/runAnalysis.sh  (26/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== v$N/systematics/vtx3_15/Unfolding ==================================
    echo "$OUT v$N/systematics/vtx3_15/Unfolding (27/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut/vtx3_15/UnfoldResults/dataResp/./

    ##-- ================================== v$N/systematics/sysVtxCut.C ==================================
    echo "$OUT v$N/systematics/sysVtxCut.C (28/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/vtxCut
    root -l -b<<EOF
.L sysVtxCut.C
sysVtxCut($N)
EOF

    ##-- ================================== v$N/systematics/clusCompatTune/runAnalysis.sh ==================================
    echo "$OUT v$N/systematics/clusCompatTune/newCCTune2pct/runAnalysis.sh (29/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== v$N/systematics/clusCompatTune/Unfolding ==================================
    echo "$OUT v$N/systematics/clusCompatTune/newCCTune2pct/Unfolding (30/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/clusCompatTune/newCCTune2pct/UnfoldResults/dataResp/./

    ##-- ================================== v$N/systematics/clusCompatTune/sysNewCC.C ==================================
    echo "$OUT v$N/systematics/clusCompatTune/newCCTune2pct/sysClusCompatTune.C (31/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/clusCompatTune
    root -l -b<<EOF
.L sysClusCompatTune.C
sysClusCompatTune($N)
EOF

    ##-- ================================== v$N/systematics/tkQuality/loose/runAnalysis.sh ==================================
    echo "$OUT v$N/systematics/tkQuality/loose/runAnalysis.sh (32/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/loose/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== v$N/systematics/tkQuality/loose/Unfolding ==================================
    echo "$OUT v$N/systematics/tkQuality/loose/Unfolding (33/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/loose/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/loose/UnfoldResults/dataResp/./

    ##-- ================================== v$N/systematics/tkQuality/tight/runAnalysis.sh ==================================
    echo "$OUT v$N/systematics/tkQuality/tight/runAnalysis.sh (34/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/tight/AnalyzerResults
    bash runAnalysis.sh $N $ETA $PMN $PMX $TEST

    ##-- ================================== v$N/systematics/tkQuality/tight/Unfolding ==================================
    echo "$OUT v$N/systematics/tkQuality/tight/Unfolding (35/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/tight/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality/tight/UnfoldResults/dataResp/./

    ##-- ================================== v$N/systematics/tkQuality/sysTkQuality.C ==================================
    echo "$OUT v$N/systematics/tkQuality/sysTkQuality.C (36/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/tkQuality
    root -l -b<<EOF
.L sysTkQuality.C
sysTkQuality($N)
EOF

    ##-- ================================== v$N/systematics/sysChi2Cutoff.C ==================================
    echo "$OUT systematics/sysChi2Cutoff.C (37/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/chi2Cutoff
    root -l -b<<EOF
.L sysChi2Cutoff.C
sysChi2Cutoff($N)
EOF

    ##-- ================================== v$N/Unfolding_Gauss ==================================
    echo "$OUT v$N/Unfolding_Gauss (38/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerangGaussResp_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/UnfoldResults/dataResp/data${N}Gauss.root

    ##-- ================================== v$N/systematics/sysRespEl.C ==================================
    echo "$OUT v$N/systematics/sysRespEl.C (39/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies/responseElements
    root -l -b<<EOF
.L sysRespEl.C
sysRespEl($N)
EOF

    ##-- ================================== statErrorHandle/SVDUnfold ==================================
    echo "$OUT statErrorHandle/SVDUnfolding (40/51)"
    cd $EBYESE/unfoldingv$N
    for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
    do
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/AnalyzerResults/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
	cp $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
	bash unfoldDataBoomerang_RooUnfoldSVD.sh $N $ROOUNFOLD
	cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA/UnfoldResults/dataResp/data${N}_svd_Split${SPLIT}.root
    done

    ##-- ================================== statErrorHandle/SVDUnfold/statuncertassess ==================================
    echo "$OUT statErrorHandle/statUncertAssessSVD.C (41/51)"
    cd $EBYESE/boomerangPlots/statErrorHandle/v$N/eta$ETA
    root -l -b<<EOF
.L statUncertAssessSVD.C
statUncertAssessSVD($N)
EOF

    ##-- ================================== v$N/SVDUnfolding ================================== 
    echo "$OUT v$N/SVDUnfolding (42/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfoldSVD.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/UnfoldResults/dataResp/data${N}_svd.root

else 
    ##-- ================================== v$N/Unfolding_dosys ==================================
    echo "$OUT v$N/Unfolding_dosys (43/51)"
    cd $EBYESE/unfoldingv$N
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/CastleEbyE.root data/PbPb_2015/data/./
    cp $EBYESE/boomerangPlots/v$N/eta$ETA/AnalyzerResults/DDResp/dataDrivenResponseAndPriors.root DDResp/./
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $EBYESE/boomerangPlots/v$N/eta$ETA/UnfoldResults/dataResp/data${N}_dosys.root

    ##-- ================================== v$N/sysUnfoldDistns.C  ==================================
    echo "$OUT v$N/sysUnfoldDistns.C (44/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies
    root -l -b<<EOF
.L sysUnfoldDistns.C
sysUnfoldDistns($N, $ETA)
EOF

    ##-- ================================== v$N/SmoothSys.C  ==================================
    echo "$OUT v$N/SmoothSys.C (45/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies
    root -l -b<<EOF
.L SmoothSys.C
SmoothSys($N)
EOF

    ##-- ================================== v$N/sysResultsPlots.C  ==================================
    echo "$OUT v$N/sysResultsPlots.C (46/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies
    root -l -b<<EOF
.L sysResultsPlots.C
sysResultsPlots($N, $ETA)
EOF

    ##-- ================================== v$N/svdUnfoldResults.C  ==================================
    echo "$OUT v$N/svdUnfoldResults.C (47/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/systematicStudies
    root -l -b<<EOF
.L svdUnfoldResults.C
svdUnfoldResults($N, $ETA)
EOF

    ##-- ================================== v$N/FitPvn.C  ==================================
    echo "$OUT v$N/FitPvn.C (48/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/
    root -l -b<<EOF
.L FitPvn.C
FitPvn($N, $ETA)
EOF

    ##-- ================================== v$N/BottomLineTest.C  ==================================
    echo "$OUT v$N/BottomLineTest.C (49/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/
    root -l -b<<EOF
.L BottomLineTest.C
BottomLineTest($N)
EOF

    ##-- ================================== v$N/FitPvn.C  ==================================
    echo "$OUT v$N/FitPvn.C (50/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/
    root -l -b<<EOF
.L FitPvn.C
FitPvn($N, $ETA)
EOF

    ##-- ================================== v$N/svdComp.C  ==================================
    echo "$OUT v$N/svdComp.C (51/51)"
    cd $EBYESE/boomerangPlots/v$N/eta$ETA/crossChecks/svdUnfold
    root -l -b<<EOF
.L svdComp.C
svdComp($N)
EOF

fi

echo "$OUT Job finished at $(date)"