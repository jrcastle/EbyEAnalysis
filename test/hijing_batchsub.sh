# Lxplus Batch Job Script
CMSSW_PROJECT_SRC="/afs/cern.ch/work/j/jcastle/CMSSW_7_5_8_patch5/src/"
CFG_FILE="qw_Hijing_afterburner.py"
WORK="/afs/cern.ch/work/j/jcastle/CMSSW_7_5_8_patch5/src/HeavyIonsAnalysis/EbyEAnalysis/test"

cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
cd $WORK
cmsRun $CFG_FILE