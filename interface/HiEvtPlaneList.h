#ifndef __HiEvtPlaneList__
#define __HiEvtPlaneList__
/*
Index     Name   Detector Order hmin1 hmax1 hmin2 hmax2 minpt maxpt nsub mcw    rmate1    rmate2
    0      HFm2        HF     2 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp2 trackmid2
    1      HFp2        HF     2  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2 trackmid2
    2     HFm2a        HF     2 -5.00 -2.00  0.00  0.00  0.01 30.00 3sub  no      HFp2   trackp2
    3     HFm2b        HF     2 -5.00 -2.50  0.00  0.00  0.01 30.00 3sub  no      HFp2   trackp2
    4     HFm2c        HF     2 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp2   trackp2
    5     HFm2d        HF     2 -5.00 -3.50  0.00  0.00  0.01 30.00 3sub  no      HFp2   trackp2
    6     HFm2e        HF     2 -5.00 -4.00  0.00  0.00  0.01 30.00 3sub  no      HFp2   trackp2
    7     HFm2f        HF     2 -5.00 -4.50  0.00  0.00  0.01 30.00 3sub  no      HFp2   trackp2
    8     HFp2a        HF     2  2.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2   trackm2
    9     HFp2b        HF     2  2.50  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2   trackm2
   10     HFp2c        HF     2  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2   trackm2
   11     HFp2d        HF     2  3.50  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2   trackm2
   12     HFp2e        HF     2  4.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2   trackm2
   13     HFp2f        HF     2  4.50  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2   trackm2
*/
#include <string>
using namespace std;
namespace hi{

  const std::string  EPNames[] = {
        "HFm2",      "HFp2",
       "HFm2a",     "HFm2b",     "HFm2c",     "HFm2d",     "HFm2e",
       "HFm2f",     "HFp2a",     "HFp2b",     "HFp2c",     "HFp2d",
       "HFp2e",     "HFp2f" 
  };

  const double EPEtaMin1[] = {
         -5.00,        3.00,
         -5.00,       -5.00,       -5.00,       -5.00,       -5.00,
         -5.00,        2.00,        2.50,        3.00,        3.50,
          4.00,        4.50 
  };

  const double EPEtaMax1[] = {
         -3.00,        5.00,
         -2.00,       -2.50,       -3.00,       -3.50,       -4.00,
         -4.50,        5.00,        5.00,        5.00,        5.00,
          5.00,        5.00 
  };

  const double minTransverse[] = {
          0.01,        0.01,
          0.01,        0.01,        0.01,        0.01,        0.01,
          0.01,        0.01,        0.01,        0.01,        0.01,
          0.01,        0.01 
  };

  const double maxTransverse[] = {
         30.00,       30.00,
         30.00,       30.00,       30.00,       30.00,       30.00,
         30.00,       30.00,       30.00,       30.00,       30.00,
         30.00,       30.00 
  };

  const std::string  EPSymmNames[] = {
     "HFpm2",    "HFpm2a",    "HFpm2b",    "HFpm2c",    "HFpm2d",
    "HFpm2e",    "HFpm2f"
  };

  const int    EPSymmPartnerBin[]  = {
    0,    0,    1,    2,    3,
    4,    5,    6,    1,    2,
    3,    4,    5,    6
  };

  const double EPSymmAbsEtaMin[] = {
    3.00,    2.00,    2.50,    3.00,    3.50,
    4.00,    4.50
  };

  const double EPSymmAbsEtaMax[] = {
    5.00,    5.00,    5.00,    5.00,    5.00,
    5.00,    5.00
  };


  static const int NumEPNames  = 14;
  static const int NEP         = 14;
  //const double epbinsDefault[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

  const int NEPSymm   = 7;
  const int EPSymmBin = 0;
  const double epbinsDefault[] = {0, 1, 2, 3, 4, 5, 6, 7};
}
#endif
