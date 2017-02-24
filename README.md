This is the offline software used to perform the offline EbyE analysis in 
CMS.  It is assumed that the user is using CMSSW, ROOT, and ROOUnfold for 
this. 

For the most part the software is location independent aside from 
three essential files (forseen fix).  The three files in question are header 
files that contain essential namespaces for virtually every macro in this 
repository:

1. EbyESEBinning.h 
2. HiEvtPlaneList.h
3. EbyECumu.h
4. DAgostiniUnfold.h (See personal_projects repository for this file)

In all macros, these header files are included assuming the software is part of a 
CMSSW release.  For example:

"HeavyIonsAnalysis/EbyEAnalysis/interface/HiEvtPlaneList.h"

In this case the file is stored in $CMSSW_BASE/src/HeavyIonsAnalysis/EbyEAnalysis/interface/.
If anyone were to use this code, they will be faced with two options.  The first option  
is to run this code in CMSSW and mimic the path to the header files and the second 
option is to go in an manually change every .C file and correct the include statements.

There are other places to with location-dependent paths that need to be fixed.  I'll 
touch on those as we go along.

Okay, now how do we use this code? The boomeranPlots and ESE directories contain most 
of the important code to be run.  In each you will find bash scripts titled 
"runEverything_v*.sh" open and read these.  This is the basic outline of what will 
be run in each case.  There are two location-dependent paths that need to be fixed in 
these file.  The first is the shell variable $EBYESE, this points to the directory that 
this file is located on your machine.  The second is in the load statements for all 
unfolding steps.  This points to the .so file for the ROOUnfold library.

Once these files have finished running (the will take on the order of 1.5 days to run) 
you will need to go into the respective v*/eta2.4 directories and run the additional 
macros not included in the bash scripts.  These will make all the pretty physics plots 
and several performance plots.

This is a very rough version of the README file, if you have further questions, please 
feel free to contact me and I'll use your feedback to improve the user friendliness 
of this repository. 


