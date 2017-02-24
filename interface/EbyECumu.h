#ifndef __EbyECumu__
#define __EbyECumu__

#include "TMath.h"
#include "TH1D.h"

class EbyECumu{

 public:

  //-- Constructor/Destructor
  EbyECumu(TH1D * h);
  ~EbyECumu();

  //-- Methods
  void CalculateCumulants();
  void CalculateMoments();
  void CalculateGamma1Exp();

  double GetCumu_vn2();
  double GetCumu_vn4();
  double GetCumu_vn6();
  double GetCumu_vn8();

  double GetCumu_vn2_Error();
  double GetCumu_vn4_Error();
  double GetCumu_vn6_Error();
  double GetCumu_vn8_Error();

  double GetGamma1Exp();
  double GetGamma1ExpError();

  double GetMoment_vn2();
  double GetMoment_vn4();
  double GetMoment_vn6();
  double GetMoment_vn8();

  double GetMoment_vn2_Error();
  double GetMoment_vn4_Error();
  double GetMoment_vn6_Error();
  double GetMoment_vn8_Error();

  double PropagateError_Cumu_vn2();
  double PropagateError_Cumu_vn4();
  double PropagateError_Cumu_vn6();
  double PropagateError_Cumu_vn8();

  double PropagateSysError_Gamma1Exp(double sysErr_Vn2, double sysErr_Vn4, double sysErr_Vn6);

 private:

  double cumu_vn2_;
  double cumu_vn4_;
  double cumu_vn6_;
  double cumu_vn8_;

  double cumu_vn2_Error_;
  double cumu_vn4_Error_;
  double cumu_vn6_Error_;
  double cumu_vn8_Error_;

  bool cumu_vn2_valid;
  bool cumu_vn4_valid;
  bool cumu_vn6_valid;
  bool cumu_vn8_valid;

  double gamma1exp_;
  double gamma1expError_;

  TH1D * hpVn_;

  double moment_vn2_;
  double moment_vn4_;
  double moment_vn6_;
  double moment_vn8_;

  double moment_vn2_Error_;
  double moment_vn4_Error_;
  double moment_vn6_Error_;
  double moment_vn8_Error_;

  int nbins_;

};

//================== Constructor ==================
EbyECumu::EbyECumu(TH1D * h)
{

  //-- Get p(vn) and normalize it
  hpVn_ = (TH1D*) h->Clone("hpVn_");
  hpVn_->Scale( 1./hpVn_->Integral() );
  nbins_ = hpVn_->GetNbinsX();

  //-- Reset all cumulants, moments and their respective errors
  cumu_vn2_         = 0.;
  cumu_vn4_         = 0.;
  cumu_vn6_         = 0.;
  cumu_vn8_         = 0.;
  cumu_vn2_Error_   = 0.;
  cumu_vn4_Error_   = 0.;
  cumu_vn6_Error_   = 0.;
  cumu_vn8_Error_   = 0.;
  moment_vn2_       = 0.;
  moment_vn4_       = 0.;
  moment_vn6_       = 0.;
  moment_vn8_       = 0.;
  moment_vn2_Error_ = 0.;
  moment_vn4_Error_ = 0.;
  moment_vn6_Error_ = 0.;
  moment_vn8_Error_ = 0.;

  gamma1exp_      = 0.;
  gamma1expError_ = 0.;

  //-- Set cumu quality flags to true
  cumu_vn2_valid = 1;
  cumu_vn4_valid = 1;
  cumu_vn6_valid = 1;
  cumu_vn8_valid = 1;

  //-- Calculate the moments of p(vn) and their errors
  CalculateMoments();

  //-- Calculate the cumulants and their errors from the moments of p(vn)
  CalculateCumulants();

  //-- Calculate gamma1exp
  CalculateGamma1Exp();  
}

//================== Destructor ==================
EbyECumu::~EbyECumu()
{
  if( hpVn_ ) delete hpVn_;
}

//================== CalculateCumulants() ==================
void EbyECumu::CalculateCumulants()
{

  double m2  = moment_vn2_;
  double m4  = moment_vn4_;
  double m6  = moment_vn6_;
  double m8  = moment_vn8_;

  double vn22 = m2;
  double vn44 = -m4 + 2.*pow(m2, 2);
  double vn66 = ( m6 - 9.*m4*m2 + 12.*pow(m2, 3) ) / 4.;
  double vn88 = -( m8 - 16.*m6*m2 - 18.*pow(m4, 2) + 144.*m4*pow(m2, 2) - 144.*pow(m2, 4) ) / 33.;

  if( vn22 < 0 ){
    std::cout << "WARNING! vn22 < 0, which means vn{2} is nan!" << std::endl;
    cumu_vn2_valid  = 0;
    cumu_vn2_       = 0.;
    cumu_vn2_Error_ = 100000;
  }
  else{
    cumu_vn2_ = pow(vn22, 1./2.);
    cumu_vn2_Error_= PropagateError_Cumu_vn2();
  }

  if( vn44 < 0 ){
    std::cout << "WARNING! vn44 < 0, which means vn{4} is nan!" <<std::endl;
    cumu_vn4_valid  = 0;
    cumu_vn4_       = 0.;
    cumu_vn4_Error_ = 100000;
  }
  else{
    cumu_vn4_ = pow(vn44, 1./4.);
    cumu_vn4_Error_= PropagateError_Cumu_vn4();
  }

  if( vn66 < 0 ){
    std::cout << "WARNING! vn66 < 0, which means vn{6} is nan!" <<std::endl;
    cumu_vn6_valid  = 0;
    cumu_vn6_       = 0.;
    cumu_vn6_Error_ = 100000;
  }
  else{
    cumu_vn6_ = pow(vn66, 1./6.);
    cumu_vn6_Error_= PropagateError_Cumu_vn6();
  }

  if( vn88 < 0 ){
    std::cout << "WARNING! vn88 < 0, which means vn{8} is nan!" <<std::endl;
    cumu_vn8_valid  = 0;
    cumu_vn8_       = 0.;
    cumu_vn8_Error_ = 100000;
  }
  else{
    cumu_vn8_ = pow(vn88, 1./8.);
    cumu_vn8_Error_ = PropagateError_Cumu_vn8();
  }

}


//================== CalculateGamma1Exp() ==================
void EbyECumu::CalculateGamma1Exp()
{

  double vn2  = cumu_vn2_;
  double vn4  = cumu_vn4_;
  double vn6  = cumu_vn6_;
  double vn8  = cumu_vn8_;

  double vn2e = cumu_vn2_Error_;
  double vn4e = cumu_vn4_Error_;
  double vn6e = cumu_vn6_Error_;
  double vn8e = cumu_vn8_Error_;

  double g1  = -6.*sqrt(2.)*pow(vn4, 2) * (vn4 - vn6) / pow( pow(vn2, 2) - pow(vn4, 2), 3./2. );

  //-- Propagate errors
  double g1e = 6*sqrt(2)*sqrt((pow(vn4,2)*(pow(vn4,4)*pow(vn4e,2)*pow(vn6,2) + pow(vn4,6)*pow(vn6e,2) + pow(vn2,2)*pow(vn4,2)*(9*pow(vn2e,2)*pow(vn4 - vn6,2) - 6*vn4*pow(vn4e,2)*vn6 + 4*pow(vn4e,2)*pow(vn6,2) - 2*pow(vn4,2)*pow(vn6e,2)) + pow(vn2,4)*(-12*vn4*pow(vn4e,2)*vn6 + 4*pow(vn4e,2)*pow(vn6,2) + pow(vn4,2)*(9*pow(vn4e,2) + pow(vn6e,2)))))/pow(pow(vn2,2) - pow(vn4,2),5));

  //-- Check cumu qualities firts
  if( !cumu_vn2_valid || !cumu_vn4_valid || !cumu_vn6_valid ){
    gamma1exp_      = -10000.;
    gamma1expError_ = 100000.;
  }
  else{
    gamma1exp_      = g1;
    gamma1expError_ = g1e;
  }

}

//================== CalculateMoments() ==================
void EbyECumu::CalculateMoments()
{

  double sum_m2 = 0.;
  double sum_m4 = 0.;
  double sum_m6 = 0.;
  double sum_m8 = 0.;

  double sumw2_m2 = 0.;
  double sumw2_m4 = 0.;
  double sumw2_m6 = 0.;
  double sumw2_m8 = 0.;

  for(int i = 1; i<= nbins_; i++){
   
    //-- Moments
    double x   = hpVn_->GetBinCenter(i);
    double px  = hpVn_->GetBinContent(i);
    double pxe = hpVn_->GetBinError(i);

    double x2 = pow(x, 2);
    double x4 = pow(x, 4);
    double x6 = pow(x, 6); 
    double x8 = pow(x, 8);

    sum_m2 += x2 * px;
    sum_m4 += x4 * px;
    sum_m6 += x6 * px;
    sum_m8 += x8 * px;

    //-- Moment Errors
    sumw2_m2 += pow( x2 * pxe, 2);
    sumw2_m2 += pow( x4 * pxe, 2);
    sumw2_m2 += pow( x6 * pxe, 2);
    sumw2_m2 += pow( x8 * pxe, 2);

  }

  //-- Set Moments
  moment_vn2_ = sum_m2;
  moment_vn4_ = sum_m4;
  moment_vn6_ = sum_m6;
  moment_vn8_ = sum_m8;

  //-- Set Moment Errors
  moment_vn2_Error_ = sqrt( sumw2_m2 );
  moment_vn4_Error_ = sqrt( sumw2_m4 );
  moment_vn6_Error_ = sqrt( sumw2_m6 );
  moment_vn8_Error_ = sqrt( sumw2_m8 );

}

//================== GetCumu_*() ==================
double EbyECumu::GetCumu_vn2(){ return cumu_vn2_; }
double EbyECumu::GetCumu_vn4(){ return cumu_vn4_; }
double EbyECumu::GetCumu_vn6(){ return cumu_vn6_; }
double EbyECumu::GetCumu_vn8(){ return cumu_vn8_; }

double EbyECumu::GetCumu_vn2_Error(){ return cumu_vn2_Error_; }
double EbyECumu::GetCumu_vn4_Error(){ return cumu_vn4_Error_; }
double EbyECumu::GetCumu_vn6_Error(){ return cumu_vn6_Error_; }
double EbyECumu::GetCumu_vn8_Error(){ return cumu_vn8_Error_; }

//================== GetGamma1Exp*() ==================
double EbyECumu::GetGamma1Exp(){ return gamma1exp_; }
double EbyECumu::GetGamma1ExpError(){ return gamma1expError_; }

//================== GetMoment_*() ==================
double EbyECumu::GetMoment_vn2() { return moment_vn2_; }
double EbyECumu::GetMoment_vn4() { return moment_vn4_; }
double EbyECumu::GetMoment_vn6() { return moment_vn6_; }
double EbyECumu::GetMoment_vn8() { return moment_vn8_; }

double EbyECumu::GetMoment_vn2_Error() { return moment_vn2_Error_; }
double EbyECumu::GetMoment_vn4_Error() { return moment_vn4_Error_; }
double EbyECumu::GetMoment_vn6_Error() { return moment_vn6_Error_; }
double EbyECumu::GetMoment_vn8_Error() { return moment_vn8_Error_; }

//================== PropagateError_Cumu_vn2() ==================
double EbyECumu::PropagateError_Cumu_vn2()
{
  double m2  = moment_vn2_;
  double m2e = moment_vn2_Error_;
  double vn2 = cumu_vn2_;

  double vn22e = sqrt(pow(m2e,2));
  double vn2e  = 0.5 * fabs( vn22e / vn2 );

  return vn2e;
}

//================== PropagateError_Cumu_vn4() ==================
double EbyECumu::PropagateError_Cumu_vn4()
{
  double m2  = moment_vn2_;
  double m2e = moment_vn2_Error_;
  double vn2 = cumu_vn2_;

  double m4  = moment_vn4_;
  double m4e = moment_vn4_Error_;
  double vn4 = cumu_vn4_;

  double vn44e = sqrt(16*pow(m2,2)*pow(m2e,2) + pow(m4e,2));
  double vn4e  = sqrt( pow(vn44e, 2) / pow(vn4, 6) ) / 4.;

  return vn4e;
}
//================== PropagateError_Cumu_vn6() ==================
double EbyECumu::PropagateError_Cumu_vn6()
{
  double m2  = moment_vn2_;
  double m2e = moment_vn2_Error_;
  double vn2 = cumu_vn2_;

  double m4  = moment_vn4_;
  double m4e = moment_vn4_Error_;
  double vn4 = cumu_vn4_;

  double m6  = moment_vn6_;
  double m6e = moment_vn6_Error_;
  double vn6 = cumu_vn6_;

  double vn66e = sqrt(81*pow(m2e,2)*pow(-4*pow(m2,2) + m4,2) + 81*pow(m2,2)*pow(m4e,2) + pow(m6e,2))/4.;
  double vn6e  = sqrt( pow(vn66e, 2) / pow(vn6, 10) ) / 6.;

  return vn6e;
}
//================== PropagateError_Cumu_vn8() ==================
double EbyECumu::PropagateError_Cumu_vn8()
{
  double m2  = moment_vn2_;
  double m2e = moment_vn2_Error_;
  double vn2 = cumu_vn2_;

  double m4  = moment_vn4_;
  double m4e = moment_vn4_Error_;
  double vn4 = cumu_vn4_;

  double m6  = moment_vn6_;
  double m6e = moment_vn6_Error_;
  double vn6 = cumu_vn6_;

  double m8  = moment_vn8_;
  double m8e = moment_vn8_Error_;
  double vn8 = cumu_vn8_;

  double vn88e = sqrt(1296*pow(-4*pow(m2,2) + m4,2)*pow(m4e,2) + 256*pow(m2e,2)*pow(36*pow(m2,3) - 18*m2*m4 + m6,2) + 256*pow(m2,2)*pow(m6e,2) + pow(m8e,2))/33.;
  double vn8e  = sqrt( pow(vn88e, 2) / pow(vn8, 14) ) / 8.;

  return vn8e;
}

//================== PropagateSysError_Gamma1Exp() ==================
double EbyECumu::PropagateSysError_Gamma1Exp(double sysErr_Vn2, double sysErr_Vn4, double sysErr_Vn6)
{

  double vn2  = cumu_vn2_;
  double vn4  = cumu_vn4_;
  double vn6  = cumu_vn6_;

  double vn2e = sysErr_Vn2;
  double vn4e = sysErr_Vn4;
  double vn6e = sysErr_Vn6;

  double g1  = -6.*sqrt(2.)*pow(vn4, 2) * (vn4 - vn6) / pow( pow(vn2, 2) - pow(vn4, 2), 3./2. );

  //-- Propagate errors
  double g1e = 6*sqrt(2)*sqrt((pow(vn4,2)*(pow(vn4,4)*pow(vn4e,2)*pow(vn6,2) + pow(vn4,6)*pow(vn6e,2) + pow(vn2,2)*pow(vn4,2)*(9*pow(vn2e,2)*pow(vn4 - vn6,2) - 6*vn4*pow(vn4e,2)*vn6 + 4*pow(vn4e,2)*pow(vn6,2) - 2*pow(vn4,2)*pow(vn6e,2)) + pow(vn2,4)*(-12*vn4*pow(vn4e,2)*vn6 + 4*pow(vn4e,2)*pow(vn6,2) + pow(vn4,2)*(9*pow(vn4e,2) + pow(vn6e,2)))))/pow(pow(vn2,2) - pow(vn4,2),5));

  return g1e;

}
#endif
