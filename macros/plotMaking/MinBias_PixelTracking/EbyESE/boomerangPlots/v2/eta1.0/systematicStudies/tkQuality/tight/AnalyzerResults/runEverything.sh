N=2
WORK=$PWD
ROOUNFOLD="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold"
EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`
cd $WORK

echo "Job started on $(date)"

root -l -b <<EOF
.x makeVNDet.C++
EOF
root -l -b <<EOF
.x ReadTree_normDet.C++
EOF
cd DDResp
root -l -b <<EOF
.x makeDDResp.C++
EOF

cd $EBYESE/unfoldingv$N
cp $WORK/CastleEbyE.root data/PbPb_2015/data/./
cp $WORK/DDResp/dataDrivenResponseAndPriors.root DDResp/./
bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
cp txt/PbPb_2015/data/data$N.root $WORK/../UnfoldResults/dataResp/./

echo "Analysis complete!"
echo "Job ended on $(date)"