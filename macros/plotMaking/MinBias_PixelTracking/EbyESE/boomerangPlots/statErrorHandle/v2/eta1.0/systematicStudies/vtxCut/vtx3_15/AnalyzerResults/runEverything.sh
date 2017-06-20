N=2
WORK=$PWD
ROOUNFOLD="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold"
EBYESE="/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/plotMaking/MinBias_PixelTracking/EbyESE"

cd /home/j550c590/CMSSW_7_5_8_patch2/src/
eval `scram runtime -sh`
cd $WORK

echo "Job started on $(date)"

#root -l -b <<EOF
#.x makeVNDet.C++
#EOF
#root -l -b <<EOF
#.x ReadTree_normDet.C++
#EOF
#cd DDResp
#root -l -b <<EOF
#.x makeDDResp.C++
#EOF

cd $EBYESE/unfoldingv3
for (( SPLIT=0; SPLIT < 10; SPLIT++ ))
do
    cp $WORK/CastleEbyE_Split${SPLIT}.root data/PbPb_2015/data/CastleEbyE.root
    cp $WORK/DDResp/dataDrivenResponseAndPriors_Split${SPLIT}.root DDResp/dataDrivenResponseAndPriors.root
    bash unfoldDataBoomerang_RooUnfold.sh $N $ROOUNFOLD
    cp txt/PbPb_2015/data/data$N.root $WORK/../UnfoldResults/dataResp/data${N}_Split${SPLIT}.root
done

echo "Analysis complete!"
echo "Job ended on $(date)"