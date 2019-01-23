//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/pdfinterface.cpp $
//$LastChangedDate: 2014-02-12 19:12:33 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1619 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include "pdfinterface.h"

#include "configBAMPS.h"

#ifdef LHAPDF_FOUND
#include "LHAPDF/LHAPDF.h"
#endif



//                                                         
// parameters of the Glück-Reya-Vogt structure function    
// for a proton     Z.Phys C 67, 433-447(1995)             
//                                                         

#define param(a,s) (a[0]+a[1]*s+a[2]*s*s)
#define param1(a,s) (a[0]+a[1]*s+a[2]*s*s+a[3]*s*s*s)
#define param2(a,s) (a[0]+a[1]*sqrt(s)+a[2]*s)
#define B(a,x) a[0]*pow(double(x),double(a[1]))*(1.0+a[2]*pow(double(x),double(a[3]))+a[4]*x+a[5]*pow(double(x),1.5))*pow(double(1.0-x),double(a[6]))
#define C(a,s,x) (pow(double(x),double(a[0]))*(a[1]+a[2]*x+a[3]*x*x)*pow(double(log(1.0/x)),double(a[4]))+pow(double(s),double(a[5]))*exp(-a[6]+sqrt(fabs(a[7]*pow(double(s),double(a[8]))*log(1.0/x)))))*pow(double(1.0-x),double(a[9]))
#define D(a,s,x) pow(double(s),double(a[0]))/pow(double(log(1.0/x)),double(a[1]))*(1.0+a[6]*sqrt(x)+a[7]*x)*pow(double(1.0-x),double(a[2]))*exp(-a[3]+sqrt(fabs(a[4]*pow(double(s),double(a[5]))*log(1.0/x))))


void interfacePDF_GRV::eval(double Q2, double x1, double x2, int nuc1, int nuc2, double F1[13], double F2[13])
{
// nuc1 and nuc2 indicate the nucleons are either proton(1) or neutron(0).

  const double mu2=0.23;
  const double lambda2=0.232*0.232;//GeV2
  double s,tmp;
  double xpp[6][10];
  static double set_u[7][3]={{ 2.284, 0.802, 0.055},
                             { 0.590,-0.024, 0.000},
                             {-0.499,-0.138,-0.076},
                             { 0.131, 0.063, 0.000},
                             { 0.213, 2.669,-0.728},
                             { 8.854,-9.135, 1.979},
                             { 2.997, 0.753,-0.076}},
                set_d[7][3]={{ 0.371, 0.083, 0.039},
                             { 0.376, 0.000, 0.000},
                             {-0.509, 3.310,-1.248},
                             { 0.486, 0.062, 0.000},
                             { 12.41,-10.52, 2.267},
                             { 6.373,-6.208, 1.418},
                             { 3.691, 0.799,-0.071}},
               set_du[7][3]={{ 0.082, 0.014, 0.008},
                             { 0.409,-0.005, 0.000},
                             {-38.07, 36.13,-0.656},
                             { 0.799, 0.071, 0.000},
                             { 90.31,-74.15, 7.645},
                             { 0.000, 0.000, 0.000},
                             { 7.486, 1.217,-0.159}},
                set_G[9][3]={{ 1.742,-0.930, 0.000},
                             { 7.486,-2.185, 0.000},
                             { 16.69,-22.74, 5.779},
                             {-25.59, 29.71,-7.296},
                             { 0.000, 0.000,-0.399},
                             { 0.524, 0.000, 0.000},
                             { 0.807, 2.005, 0.000},
                             { 3.841, 0.316, 0.000},
                             { 1.088, 0.000, 0.000}},
                  set_GD[4]= { 2.792, 2.215, 0.422,-0.104},     
              set_ud[10][3]={{ 0.410,-0.232, 0.000},
                             { 0.890,-0.140, 0.000},
                             {-0.981, 0.000, 0.000},
                             { 0.320, 0.683, 0.000},
                             { 0.534,-0.457, 0.000},
                             { 1.451, 0.000, 0.000},
                             { 4.119, 1.713, 0.000},
                             { 0.682, 2.978, 0.000},
                             { 0.271, 0.000, 0.000},
                             { 4.752, 1.164, 0.286}},
                set_s[6][3]={{ 0.914, 0.000, 0.000},
                             { 1.798,-0.596, 0.000},
                             { 6.379,-0.350, 0.142},
                             { 3.981, 1.638, 0.000},
                             { 6.402, 0.000, 0.000},
                             { 0.577, 0.000, 0.000}},
              set_sAB[2][3]={{-5.548, 3.669,-0.616},
                             { 18.92,-16.73, 5.168}};

  s=log(log(Q2/lambda2)/log(mu2/lambda2));

  //parameters for xu_v, xd_v, x(dbar-ubar)
  for (int i=0; i<7; i++)
  {
    xpp[0][i]=param(set_u[i],s);
    xpp[1][i]=param(set_d[i],s);
    xpp[2][i]=param(set_du[i],s);
  }

  //parameters for xG
  for (int j=0; j<9; j++)
  {
    xpp[3][j]=param(set_G[j],s);
  }
  xpp[3][9]=param1(set_GD,s);

  //parameters for x(ubar+dbar)
  for (int i=0; i<10; i++)
  {
    xpp[4][i]=param(set_ud[i],s);
  }

  //parameters for xs
  for (int j=0; j<6; j++)
  {
    xpp[5][j]=param(set_s[j],s);
  }
  for (int i=0; i<2; i++)
  {
    xpp[5][i+6]=param2(set_sAB[i],s);
  }

  // [0]g,[1]u,[2]ubar,[3]d,[4]dbar,[5]s,[6]sbar,[7]c,[8]cbar 

  F1[2]=(-B(xpp[2],x1)+C(xpp[4],s,x1))/2.0;
  F1[4]=(B(xpp[2],x1)+C(xpp[4],s,x1))/2.0;
  F1[1]=B(xpp[0],x1);
  F1[3]=B(xpp[1],x1);
  if(nuc1==0)
  {
    tmp=F1[1];
    F1[1]=F1[3];
    F1[3]=tmp;
  }
  F1[1]=F1[1]+F1[2];
  F1[3]=F1[3]+F1[4];
  F1[5]=F1[6]=D(xpp[5],s,x1);
  F1[7]=F1[8]=0.0;
  F1[0]=C(xpp[3],s,x1);
  for (int i=8; i<13; i++) F1[i] = 0.0;

  F2[2]=(-B(xpp[2],x2)+C(xpp[4],s,x2))/2.0;
  F2[4]=(B(xpp[2],x2)+C(xpp[4],s,x2))/2.0;
  F2[1]=B(xpp[0],x2);
  F2[3]=B(xpp[1],x2);
  if(nuc2==0)
  {
    tmp=F2[1];
    F2[1]=F2[3];
    F2[3]=tmp;
  }
  F2[1]=F2[1]+F2[2];
  F2[3]=F2[3]+F2[4];
  F2[5]=F2[6]=D(xpp[5],s,x2);
  F2[7]=F2[8]=0.0;
  F2[0]=C(xpp[3],s,x2);
  for (int i=8; i<13; i++) F2[i] = 0.0;
}


interfacePDF_LHAPDF::interfacePDF_LHAPDF ( const std::string& name, const unsigned int member, const bool useGrid, const bool _nPDF, const std::string _nPDFname, int _A1, int _A2 ) :
  A1( _A1 ),
  A2( _A2 ),
  useNuclearPDFs(_nPDF)
{
#ifdef LHAPDF_FOUND
  //  LHAPDF::setVerbosity( LHAPDF::SILENT );

  // initialize use of nuclear PDFs if requested
  if ( useNuclearPDFs && ( ( A1 > 1 ) || ( A2 > 1 ) ) )
  {
    LHAPDF::setParameter( _nPDFname.c_str() );
    std::cout << "Using nuclear PDFs from data set: " << _nPDFname << std::endl;
  }
  
  if ( useGrid )
  {
    LHAPDF::initPDFSet( name, LHAPDF::LHGRID );
  }
  else
  {
    LHAPDF::initPDFSet( name, LHAPDF::LHPDF );
  }
  LHAPDF::usePDFMember( member );
  
#else
  std::string errMsg = "LHAPDF not available";
  throw ePDF_error( errMsg );
#endif
}


void interfacePDF_LHAPDF::eval(double Q2, double x1, double x2, int nuc1, int nuc2, double F1[13], double F2[13])
{
#ifdef LHAPDF_FOUND
  double F[13];
  double Q = sqrt(Q2);

  if ( useNuclearPDFs && A1 > 1) 
  {
    LHAPDF::xfxa(x1, Q, A1, F);
  }
  else
  {
    LHAPDF::xfx(x1, Q, F);
  }
  Map(nuc1, F, F1);

  if ( useNuclearPDFs && A2 > 1) 
  {
    LHAPDF::xfxa(x2, Q, A2, F);
  }
  else
  {
    LHAPDF::xfx(x2, Q, F);
  }
  Map(nuc2, F, F2);
#endif
}


void interfacePDF_LHAPDF::Map(const int charge, const double F_LHAPDF[13], double F[13])
{
  switch (charge)
  {
    case 1: // proton
    {
      F[ 1] = F_LHAPDF[ 8]; // u
      F[ 2] = F_LHAPDF[ 4]; // ubar
      F[ 3] = F_LHAPDF[ 7]; // d
      F[ 4] = F_LHAPDF[ 5]; // dbar
    } break;
    case 0: // neutron
    {
      F[ 1] = F_LHAPDF[ 7]; // u !!
      F[ 2] = F_LHAPDF[ 5]; // ubar !!
      F[ 3] = F_LHAPDF[ 8]; // d !!
      F[ 4] = F_LHAPDF[ 4]; // dbar !!
    } break;
    default:
    {
      std::string errMsg = "wrong nucleon charge";
      throw ePDF_error( errMsg );
    }
  }
  F[ 0] = F_LHAPDF[ 6]; // g
  F[ 5] = F_LHAPDF[ 9]; // s
  F[ 6] = F_LHAPDF[ 3]; // sbar
  F[ 7] = F_LHAPDF[10]; // c
  F[ 8] = F_LHAPDF[ 2]; // cbar
  F[ 9] = F_LHAPDF[11]; // b
  F[10] = F_LHAPDF[ 1]; // bbar
  F[11] = F_LHAPDF[12]; // t
  F[12] = F_LHAPDF[ 0]; // tbar
}
