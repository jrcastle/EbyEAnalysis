#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"

#include <iostream>
using namespace std;
using namespace ebyese;

void makeEvtSplitTree(){

  TRandom3 * ran;

  //-- Analyzer Tree
  TFile * fAna;
  TTree * treeAna;

  //-- Ouput Tree
  TFile * fOut;
  TTree * treeSplit;

  //-- Event splitter
  int iSplit;

  //
  // MAIN
  //

  //-- Grab the analyzer tree
  fAna = new TFile(fAnaTreeName);
  treeAna = (TTree*) fAna->Get("ebyeana/tree");

  //-- Set up the output tree
  fOut = new TFile("/rfs/jcastle/PbPb2015/PixelTracking_MB2/SplitTree.root", "recreate");
  treeSplit = new TTree("SplitTree", "SplitTree");
  treeSplit->Branch("iSplit", &iSplit, "iSplit/I");

  //-- Set the RNG
  ran = new TRandom3(0);

  int N = treeAna->GetEntries();

  for(int ievent = 0; ievent < N; ievent++){

    if((ievent+1)% 500000 == 0) std::cout << "Processing Event " << ievent+1 << "\t" << (100.*(ievent+1)/N) << "% Completed" << std::endl;
    iSplit = ran->Uniform(0, NSPLIT);
    treeSplit->Fill();

  }

  fOut->Write();

}




