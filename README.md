============================================= ANALYZER =============================================
This setup is designed to run over PbPb minimum bias data/MC.  It is capable of running over 
GEN, hiGeneralTracks, and hiGeneralAndPixelTracks.  To see which variables you can set for the 
analyzer, look in python/ebyeana_cfi.py. The analyzer is located in src/EbyeAnalyzer.cc. The 
test/ directory contains some example python configurations for cmsRun and CRAB3 as well as the 
golden JSON for 2015 PbPb data. The example cmsRun configuration file has several flags for various 
systematic studies, be sure you know which are set before you run the analyzer.  These settings 
are conveniently commented with "PAY ATTENTION" to help guide your eye. After the analyzer is 
finished running, head over to the macros directory to run all the 
unfolding procedures, systematic studies, and plotting macros.

========================================== POST ANALYZER ==========================================
This is the offline software used to perform the offline EbyE analysis in CMS.  It is assumed that 
the user is using CMSSW, ROOT, and ROOUnfold for this. For the most part the software is location 
independent aside from three essential files.  The three files in question are header files that 
contain essential namespaces/classes for virtually every macro in this repository:

1. EbyESEBinning.h 
2. HiEvtPlaneList.h
3. EbyECumu.h

Check these files and make sure the pt/eta/centrality binning matches those in the analyzer. In all 
macros, these header files are included assuming the software is part of a CMSSW release.  For 
example:

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"

is stored in $CMSSW_BASE/src/HeavyIonsAnalysis/EbyEAnalysis/interface/. If anyone were to use this 
code, they will be faced with two options.  The first option  is to run this code in CMSSW and 
mimic the path to the header files and the second option is to go in an manually change every .C 
file and correct the include statements. There are other places to with location-dependent paths 
that need to be fixed.  I'll touch on those as we go along.

Okay, now how do we use this code? The boomerangPlots and ESE directories contain most of the 
important code to be run.  In each you will find bash scripts titled "runEverything*.sh" open and 
read these.  This is the basic outline of what will be run in each case.  There are two 
location-dependent paths that need to be fixed in these file.  The first is the shell variable 
$EBYESE, this points to the directory that this file is located on your machine.  The second is in 
the $ROOUNFOLD variable in the load statements for all unfolding steps.  This variable needs to 
point to the .so file for the ROOUnfold library.

The runEverything.sh script has two "passes." The first pass does a majority of the work, but in 
order to complete the response matrix systematic shape study, the ROOUnfold software needs to be 
hacked.  Head over to the ROOUnfold source code and edit the RooUnfoldBayes.cxx file.  Here, the 
variable _dosys needs to be set to true.  This will include the uncertainties on the response 
matrix elements in error propagation. When finished navigate up one directory from src and 
recompile the ROOUnfold code by typing "make." Now, head over to the runEverything.sh script and 
set FIRSTPASS="0" and rerun. Once this has finished running, all the pretty physics plots and 
several performance plots can be found in the v$N/eta$ETA/plots/ directory. Chances are likely you 
will need to go into the various plotting macros and manually change the settings.  All the default 
settings are based on v2, |eta|<1.0, 0.3<pT<3.0 GeV. 

If I haven't emphasized it enough, the runEverything.sh file is key in understanding the workflow. 
Each plotting macro has flags to catch if prerequisit procedures have been run prior to it and will 
exit with an error message if they haven't. To start a new analysis set the $N, $ETA, $PMN, and 
$PMX variables in the script and the script will then check to see if there is already a directory 
containing software for this.  If not, it will copy from the TEMPLATE directory and begin working. 
The software in the ESE directory is very outdated and needs further development. 

This is a very rough version of the README file, if you have further questions, please feel free 
to contact me and I'll use your feedback to improve the user friendliness of this repository. 


