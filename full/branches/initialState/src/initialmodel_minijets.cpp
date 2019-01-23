//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_minijets.cpp $
//$LastChangedDate: 2014-02-12 19:12:33 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1619 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>

#include "configuration.h"
#include "initialmodel_minijets.h"
#include "pdfinterface.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_minijets::initialModel_minijets( const config& _config ) :
  initialModelWS(_config),
  nTestparticles( _config.getTestparticles() ),
  nParticlesToGenerate( 0 )
{
  distrPT = 0;  // assign a null-pointer
  maxIntegrandPT = 0; // assign a null-pointer
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  P0 = _config.getPtCutoff();
  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw eMiniJet_error( errMsg );
  }

  switch ( _config.getPDFsource() )
  {
    case builtInGRV:
      PDF = new interfacePDF_GRV;
      break;
    case LHAPDF:
      PDF = new interfacePDF_LHAPDF ( _config.getLHAPDFdatasetName(), _config.getLHAPDFmember(), _config.getLHAPDFuseGrid(), 
        _config.useNuclearPDFs(), _config.getNuclearPDFdatasetName(), _config.getA(), _config.getB() );
      break;
    default:
      throw eMiniJet_error( "error in initialModel_minijets: invalid specification of PDFsource" );
  }
    
  generateSamplingDataSets();
}

initialModel_minijets::~initialModel_minijets()
{
  delete distrPT;
  delete maxIntegrandPT;
}



void initialModel_minijets::generateSamplingDataSets()
{
  double Tab;
  double sigmaJet;

  cout << "======= Generating data sets for sampling of initial state =======" << endl;

  generateTimeDistributionWS(Tab);
  cout << "++++  Tab = " << Tab << "1/mb" << endl;

  generatePtDistribution(sigmaJet);
  cout << "++++  sigma_jet = " << sigmaJet << " mb" << endl;

  nParticlesToGenerate = 2 * static_cast<int>( Tab * sigmaJet ) * nTestparticles;
  cout << "++++ N = " << nParticlesToGenerate / nTestparticles 
       << " * " << nTestparticles
       << " = " << nParticlesToGenerate << " testparticles" << endl;
  cout << "==================================================================" << endl;
}



void initialModel_minijets::populateParticleVector( std::vector< Particle >& _particles )
{
  /*
   * Reserve memory for the Particle vector. Due to possible particle
   * creation this size should be chosen somewhat larger than the
   * initial number of particles. Doing this, push_back operations to
   * add new particles won't lead to internal memory re-allocation in
   * the vector, which could possibly be time consuming. 
   * However push_back operations are always possible, even when the
   * reserved size is used up (as long as there is physical memory
   * available of course). And modern day compilers probably optimize
   * to the extent where it doesn't really matter - but still, it's
   * good practice. 
   */
  _particles.reserve( static_cast<int>( nParticlesToGenerate * 1.2 ) );
  /*
   * Now the particle vector is re-sized to hold the initial number of
   * particles (this is NOT done when reserving memory!). The particles
   * are initialised with the standard constructor, actual attributes
   * such as momentum etc. MUST be set later! 
   */
  _particles.resize( nParticlesToGenerate );

  sample_PXYZE_FLAV( _particles );
  sample_TXYZ( _particles );
  
  // the routines are not suitable for heavy quarks. Therefore delete all of them
  int number = _particles.size();
  for(unsigned int j = 0; j < _particles.size(); j++ )
  {
    if( _particles[j].FLAVOR > 2 * Particle::max_N_light_flavor )
    {
      while( _particles.back().FLAVOR > 2 * Particle::max_N_light_flavor )
      {
        _particles.pop_back();
      }
      _particles[j] = _particles.back();
      _particles.pop_back();
    }
  }
  cout << "Minijet sampling is not written for heavy quarks. " << number - _particles.size() << " heavy quarks out of " << number << " particles have been deleted." << endl;
}



void initialModel_minijets::sample_TXYZ( std::vector< Particle >& _particles )
{
  for ( unsigned int index = 0; index < _particles.size(); index += 2 )
  {
    sample_TXYZ_singleParticle( _particles[index] );
    _particles[index + 1].Pos = _particles[index].Pos;
  }
}



void initialModel_minijets::sample_PXYZE_FLAV( std::vector< Particle >& _particles ) const
{
  double Y1, Y2, Y1max, Y2max, Y2min, xT;
  double XS, maxXS, tryXS;
  int f1, f2;
  static double arr[13][13];

  for ( unsigned int index = 0; index < _particles.size(); index += 2 )
  {
    // (1): Sample pT:

    double pT = distrPT->GenerateNumber();

    // (1a): set pX and pY:

    double theta = 2.0 * M_PI * ran2();
    _particles[index].Mom.Px() = pT * cos( theta );
    _particles[index].Mom.Py() = pT * sin( theta );
    _particles[index+1].Mom.Px() = - _particles[index].Mom.Px();
    _particles[index+1].Mom.Py() = - _particles[index].Mom.Py();

    // (2): Sample y1 and y2:
    //sampling via rejection method: comparison function is
    //f(y)=cs(y)=const=max, i.e. maximum value of cross section
    //for given PT when varying Y1 and Y2. Uniformly distributed
    //point is choosen: [Y,f(Y)] = [Ymin..Ymax , 0..max]. If
    //csjet(Y) < f(Y) the value for Y is accepted, otherwise it is rejected.

    // (2a): get maximum value  of integrand:
    maxXS = maxIntegrandPT->eval(pT);

    // (2b): set kinematic boundaries
    integrand_distPT fdistPT( pT, sqrtS_perNN, PDF);
    xT    = fdistPT.getXT();
    Y1max = fdistPT.getY1max();


    do
    {
      Y1 = 2 * Y1max * ran2() - Y1max;
      Y2 = 2 * Y1max * ran2() - Y1max;

      Y2max = 2 / xT - exp(  Y1 );
      Y2min = 2 / xT - exp( -Y1 );

      if (( Y2max <= 0 ) || ( Y2min <= 0 ) )
      {
        throw eMiniJet_error( "problem2 in sample_PXYZE_FLAV" );
      }

      Y2max =  log( Y2max );
      Y2min = -log( Y2min );

      XS = (( Y2>Y2max ) || ( Y2<Y2min )) ? -1.0 : fdistPT.calculate(Y1,Y2, 1,1, arr);

      if ( XS > maxXS )
      {
        std::stringstream tempStr;
        tempStr << "problem3 in sample_PXYZE_FLAV  pt = " << pT << "  XS = " << XS << "  maxXS = " << maxXS;
        std::string errMsg = tempStr.str();
        throw eMiniJet_error( errMsg );
      }

      tryXS = maxXS * ran2();
    }
    while ( XS < tryXS );

    _particles[index].Mom.Pz()   = 0.5 * pT * ( exp( Y1 ) - exp( -Y1 ) );
    _particles[index+1].Mom.Pz() = 0.5 * pT * ( exp( Y2 ) - exp( -Y2 ) );
    _particles[index].Mom.E()    = sqrt( _particles[index].Mom.vec2() );
    _particles[index+1].Mom.E()  = sqrt( _particles[index+1].Mom.vec2() );

    // (3): Sample the flavours:

    // (3a): The nucleons are sampled to being either proton(1) or
    // neutron(0).

    unsigned int n1 = ( A*ran2() < Aatomic ) ? 1 : 0;
    unsigned int n2 = ( B*ran2() < Batomic ) ? 1 : 0;

    // (3b): Calculate the isospin corrected XS
    maxXS = fdistPT.calculate(Y1,Y2, n1,n2, arr);
    if (maxXS == 0.0) 
    {
      throw eMiniJet_error( "problem4 in sample_PXYZE_FLAV" );
    }

    // (3c): get some random number and loop over all final states
    tryXS = maxXS * ran2();

    f1 = f2 = -1;
    XS = 0.0;
    for (int i=0;i<13;i++)
    {
      for (int j=0;j<13;j++)
      {
        XS += arr[i][j];
        if (XS >= tryXS) // we have found the solution!
        {
          f1 = i;
          f2 = j;
          i = j = 99; // --> quit the loops
        }
      }
    }

    if (f1 < 0) 
    {
      throw eMiniJet_error( "problem5 in sample_PXYZE_FLAV" );
    }

    _particles[index].FLAVOR   = static_cast<FLAVOR_TYPE>( f1 );
    _particles[index+1].FLAVOR = static_cast<FLAVOR_TYPE>( f2 );

    _particles[index].m   = Particle::getMass( _particles[index].FLAVOR );
    _particles[index+1].m = Particle::getMass( _particles[index+1].FLAVOR );
  }

}



void initialModel_minijets::generatePtDistribution( double& sigma_jet )
{
  double PT;
  double factor, end, dpt, old, s, y;
  double tgral, sd, chi2a, jetdist;
  double mm;

  std::vector<double> valPT;
  std::vector<double> valIntegral;
  std::vector<double> valMax;

  cout << "---- generatePtDistribution:" << endl;

  // reserve some minimal size for the arrays:
  valPT.reserve(200);
  valIntegral.reserve(200);
  valMax.reserve(200);
  

  integrand_distPT fdistPT( 0.0, sqrtS_perNN, PDF);

  factor = 1. / 2.5682;// transform dimension from 1/GeV^2 to mb
  factor = factor * K;

  end = sqrtS_perNN / 5;   // large and somewhat arbitrary PT cut < sqrtS/2  
  dpt = 0.01;

  old = infinity;
  s = 0.0;

  PT = P0;
  do
  {
    fdistPT.setPT( PT );
    vegas( 2, fdistPT, &tgral, &sd, &chi2a );  //integrate dsigma/(dPT^2*dY1*dY2) over Y1 and Y2
    jetdist = double( tgral ) * factor * 2.0 * PT;
    mm = fdistPT.maxIntegrand( y );

    valPT.push_back( PT );
    valIntegral.push_back( jetdist ); // we first fill the array with the function values
    valMax.push_back( mm );

    if (( old - jetdist ) < 0.001 )
    {
      dpt = dpt * 10.;
      dpt = ( dpt < 1.0 ) ? dpt : 1.0;
    }
    old = jetdist;

    PT += dpt;
  }
  while ( PT <= end );

  // we add some 'zero' at the end:
  valPT.push_back( PT );
  valIntegral.push_back( 0.0 );
  valMax.push_back( 0.0 );

  maxIntegrandPT = new interpolationGSL(valPT,valMax, interp_cspline);

  // now we calculate the integral from the end:

  std::reverse(valPT.begin(),valPT.end());
  std::reverse(valIntegral.begin(),valIntegral.end());
  
  old = 0.0;
  s = 0.0;
  for (unsigned int i=1; i<valPT.size(); i++)
  {
    dpt = valPT[i]-valPT[i-1];

    s += 0.5 * (old + valIntegral[i]) * dpt;
    old = valIntegral[i];
    valIntegral[i] = s;
  }

  for (unsigned int i=0; i<valPT.size(); i++)
  {
    valIntegral[i] /= s; 
  }

  distrPT = new ranGen_Distr(valPT,valIntegral, interp_cspline);

  sigma_jet = - s / 2.0; // minus sign because of std::reverse(...)

  cout << "---- generatePtDistribution: finished" << endl;
}


void initialModel_minijets::Plot(void)
{
  double PT;
  double factor, end, dpt;
  double tgral, sd, chi2a, jetdist;

  integrand_distPT fdistPT( 0.0, sqrtS_perNN, PDF);

  factor = 1. / 2.5682;// transform dimension from 1/GeV^2 to mb
  factor = factor * K;

  end = 35.0;//large PT cut < sqrtS/2
  dpt = 0.01;

  cout << endl << endl;

  PT = P0;
  do
  {
    fdistPT.setPT( PT );
    vegas( 2, fdistPT, &tgral, &sd, &chi2a );  //integrate dsigma/(dPT^2*dY1*dY2) over Y1 and Y2
    jetdist = double( tgral ) * factor * 2.0 * PT;

    cout << PT << " " << jetdist << endl;
    PT += dpt;
  }
  while ( PT <= end );

  std::string errMsg = "stop the code";
  throw eMiniJet_error( errMsg );
  
}




void integrand_distPT::operator()( const int *ndim, const double xx[], const int *ncomp, double ff[] ) const
{
  double wgt = 0;
  ff[0] = this->operator()( xx, wgt );
}

/**
 * Y1 and Y2 run from min[Yi] to max[Yi] (i=1,2) as x[1] and x[2] run
 * from 0 to 1
 **/
double integrand_distPT::operator()( const double x[], double wgt ) const
{
  double Y1max, Y2max, Y2min, Y1, Y2, V;
  static double arr[13][13];

  Y1max = 1.0 / xT + sqrt( 1.0 / (xT*xT) - 1.0 );
  if ( Y1max <= 0. )
  {
    cout << "problem 0 in fdistPT" << endl;
    return 0.0;
  }
  Y1max = log( Y1max );
  Y1 = Y1max * double( x[1] );

  Y2max = 2.0 / xT - exp( Y1 );
  Y2min = 2.0 / xT - exp( -Y1 );
  if (( Y2max <= 0. ) || ( Y2min <= 0. ) )
  {
    cout << "problem 1 in fdistPT" << endl;
    return 0.0;
  }
  Y2max = log( Y2max );
  Y2min = -log( Y2min );

  Y2 = ( Y2max - Y2min ) * double( x[2] ) + Y2min;

  V = 2.0 * Y1max * ( Y2max - Y2min );//factor 2.0 due to the symmetry in Y1

  return double( V*calculate(Y1,Y2, 1,1, arr) );
}

/**
 * we only give y1 and y2 as a function argument here, because sqrtS
 * and pT have been already set as the class parameters.
 **/
double integrand_distPT::calculate(const double Y1, const double Y2, const int n1, const int n2, double arr[13][13]) const
{
  double pT2, x1, x2;
  double s,t,u;
  double F1[13], F2[13];
  double total;


  // set default return values:
  total = 0.0;
  for (int i=0;i<13;i++)
    for (int j=0;j<13;j++)
      arr[i][j] = 0.0;


  x1 = 0.5 * xT * ( exp(  Y1 ) + exp(  Y2 ) );
  x2 = 0.5 * xT * ( exp( -Y1 ) + exp( -Y2 ) );

  if ((( 1.0 - x1 ) < 0.0001 ) || (( 1.0 - x2 ) < 0.0001 ) ) return 0.0;
  
  pT2 = pT*pT;
  s = x1 * x2 * sqrtS*sqrtS;
  t = -pT2 * ( 1. + exp( Y2 - Y1 ) );
  u = -pT2 * ( 1. + exp( Y1 - Y2 ) );
  t /= s; // in the massless case, we rescale t and u with s
  u /= s;
  
  PDF->eval( pT2, x1,x2, n1,n2, F1,F2 );

  // [0] preparational work:
  // --- sum of all quark-antiquark products:
  double ff0 = 0.0;
  for (int i=1;i<13;i+=2) ff0 += F1[i]*F2[i+1] + F2[i]*F1[i+1];
  //  for (int i=1;i<9;i+=2) ff0 += F1[i]*F2[i+1] + F2[i]*F1[i+1];
  // --- some abbreviations:
  double c3 = qaqd(t,u) + qaqd(u,t);
  double c4 = F1[0]*F2[0] * ( ggqq(t,u) + ggqq(u,t) );

  // [1] gluon-gluon in the output channel: ( qqbar -> gg, gg -> gg )
  arr[0][0] = ff0 * qqgg(t,u) + F1[0]*F2[0] * gggg(t,u);

  // [2] gluon-quark in the output channel:
  for (int i=1;i<13;i++)
  {
    arr[0][i] = F1[0]*F2[i] * gqgq(t,u) + F2[0]*F1[i] * gqgq(u,t);
    arr[i][0] = F1[0]*F2[i] * gqgq(u,t) + F2[0]*F1[i] * gqgq(t,u);
  }
  
  // [3] quark-antiquark in the output channel: 
  for (int i=1;i<13;i+=2)
  {
    double d1 = F1[i]*F2[i+1];
    double d2 = F2[i]*F1[i+1];
    double h = (ff0-d1-d2) * c3 + c4;

    arr[i][i+1] = d1 * qaqs(t,u) + d2 * qaqs(u,t) + 0.5 * h;
    arr[i+1][i] = d1 * qaqs(u,t) + d2 * qaqs(t,u) + 0.5 * h;
  }

  // [4] quark-quark in the output channel: 
  // ... the symmetric case: (i1==i2)
  for (int i=1;i<13;i++)
  {
    arr[i][i] = F1[i]*F2[i] * qqqq(t,u);
  }

  // ... the asymmetric case: (i1 != i2)
  // (in addition, we have to avoid the quark-antiquark case
  // (i2==i1+1), since this is already covered above)

  for (int i=1;i<13;i+=2)
  {
    for (int j=i+2;j<13;j++)
    {
      arr[i][j] = F1[i]*F2[j] * qq12(t,u) + F1[j]*F2[i] * qq12(u,t);
      arr[j][i] = F1[i]*F2[j] * qq12(u,t) + F1[j]*F2[i] * qq12(t,u);

      arr[i+1][j] = F1[i+1]*F2[j] * qq12(t,u) + F1[j]*F2[i+1] * qq12(u,t);
      arr[j][i+1] = F1[i+1]*F2[j] * qq12(u,t) + F1[j]*F2[i+1] * qq12(t,u);
    }
  }

  // the following three lines are just for 'historical' reasons. 
  // in the future, alpha_s should be consistent with the used PDF.

  const int flav = 4;//u,d,s,c Quark
  double alpha_s = 12.*M_PI / ( 33. - 2.*flav ) / log( pT2 / lambda2 );
  double factor = M_PI * alpha_s * alpha_s / ( s * s );

  // calculate the sum of all entries:
  for (int i=0;i<13;i++)
  {
    for (int j=0;j<13;j++)
    {
      arr[i][j] *= factor;
      total += arr[i][j];
    }
  }

  //...only gg:
  // total = arr[0][0];

  //...only gq:
  // total = 0.0;
  // for (int i=1;i<13;i++) total += arr[0][i]+arr[i][0];

  //..only qq:
  // total = 0.0;
  // for (int i=1;i<13;i++)
  //   for (int j=1;j<13;j++)
  //     total += arr[i][j];

  return total;

}

/**
 * give the maximum value of G~d(sigma)/d(pt2)d(y1)d(y2)
 * for a fixed PT for preparing the sampling in y1-y2 plane
 **/

double integrand_distPT::maxIntegrand( double& y ) const
{
  const int nm = 1000;
  double Ymax, Y1, Y2, tmp, maxx;
  static double arr[13][13];

  Ymax = log( 1.0 / xT + sqrt( 1.0 / (xT*xT) - 1.0 ) );

  maxx = calculate(0.0,0.0, 1,1, arr);
  y = 0.0;

  for ( int i = 0;i < nm;i++ )
  {
    // for symmetric systems:
    //    Y1 = Ymax * ran2();
    //    Y2 = -Y1;
    // for asymmetric systems
    Y1 = Ymax * (2.0*ran2() - 1.0);
    Y2 = Ymax * (2.0*ran2() - 1.0);

    tmp = calculate(Y1,Y2, 1,1, arr);
    if ( tmp > maxx )
    {
      maxx = tmp;
      y = Y1;
    }
  }
  return 1.03*maxx;// 3% uncertainty
}

double integrand_distPT::getY1max(void) const
{
  if (xT==0.0) 
  {
    std::string errMsg = "getY1max: xT==0!";
    throw eMiniJet_error( errMsg );
  }

  double Y1max = 1 / xT + sqrt( 1 / (xT*xT) - 1.0 );
  if ( Y1max <= 0 )
  {
    std::string errMsg = "getY1max: problem 1";
    throw eMiniJet_error( errMsg );
  }
  return log( Y1max );
}



// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
