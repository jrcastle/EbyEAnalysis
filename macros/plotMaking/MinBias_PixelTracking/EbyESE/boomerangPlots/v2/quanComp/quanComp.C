#include "HeavyIonsAnalysis/EbyEAnalysis/interface/EbyESEBinning.h"
using namespace ebyese;

void quanComp(){

  const int N = 11;
  const int cmin[N] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  const int cmax[N] = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  const int ccol[N] = {};


  bool pixel_hiGeneralAndPixel = 0;
  bool pixel_hiGeneral         = 1;


  TFile * fQuan;
  TH1D * hMultQuan[N];

  TFile * fJames;
  TH2D hMultCentJames;
  TH1D * hMultJames;


  //
  // MAIN
  //

  fQuan = 0;
  if(pixel_hiGeneralAndPixel) fQuan = new TFile("HIMB2Pixel_pixel.root");
  if(pixel_hiGeneral)         fQuan = new TFile("HIMB2Pixel_general.root");
  if(pixel_hiGeneralAndPixel && pixel_hiGeneral){
    std::cout << "Pick one input Quan file and try again..." << std::endl;
    exit(0);
  }










}
