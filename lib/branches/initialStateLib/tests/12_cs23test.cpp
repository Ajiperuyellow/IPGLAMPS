//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/12_cs23test.cpp $
//$LastChangedDate: 2014-12-11 22:00:00 +0100 (Do, 11. Dez 2014) $
//$LastChangedRevision: 2015 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the sampling of t in 2->2 processes.
*/

// Vegas
#define NSTART 50
#define NINCREASE 100
#define NBATCH 1000

//-----------------------------------------------------------


#include <iostream>
using namespace std;

#include "FPT_compare.h"
#include "binning.h"
#include "binary_cross_sections.h"
#include <fstream>
#include "interpolation23.h"
#include "integrand23.h"
#include "integrate.h"
#include <time.h>
#include "lorentz.h"

/* 
 * Number of clock ticks per second. A clock tick is the unit by which 
 * processor time is measured and is returned by 'clock'. 
 */ 
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC  1000l
#endif


// Some constants:
const bool running_coupling = true;
const int Nf_light = 3;
const int Nf_heavy = 0;
const double maxRunningCoupling = 1.5;

void totalcs_output(void)
{
  const double T = 0.2;
  // Boltzmann definition
  const double md2g_wo_as = ( ns_casc::Ncolor + ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
  const double md2q_wo_as = 16 / ( 3 * M_PI ) * pow( T, 2.0 ); //GeV^2
//   const double md2g_wo_as = 0; //GeV^2
//   const double md2q_wo_as = 0; //GeV^2
  
  stringstream filename_sigma;
  filename_sigma << "output/sigma23_runningCoupling" << coupling::isRunning() << "asMax" << maxRunningCoupling << "T" << T << ".f1";
  ofstream file_sigma( filename_sigma.str().c_str() );

  const int n_s = 10000;
  const int n_lambda = 6;
  const double s_min = 0.1;
  const double s_max = 1000.0;
  const double ds = (s_max - s_min) / n_s;

  const double lambda_min = exp(-3.0);
  const double lambda_max = exp(3.0);
  const double dlambda = ( lambda_max - lambda_min ) / n_lambda;
  
  lorentz LL_CM;
  VectorEPxPyPz P1cell,P2cell,P1cm,P2cm;
  
  P1cell = VectorEPxPyPz( sqrt(2.0), 1.0, 1.0, 0.0 );
  P2cell = VectorEPxPyPz( sqrt(2.0), -1.0, 1.0, 0.0 );
  
  LL_CM.setBetaCM(P1cell,P2cell);
  LL_CM.boost(P1cell,P2cell, P1cm,P2cm);

  const double beta = LL_CM.betaVal();
  const double cos_theta = CosTheta(LL_CM.beta(), P1cm);
  cout << beta << "\t" << cos_theta << endl;
  
//   cout << "T" << "  " << "md2g_wo_as" << "  " << "md2q_wo_as" << "  " << "s[si]" << endl;
  cout << T << "  " << md2g_wo_as << "  " << md2q_wo_as << endl;
  
//   file_dsigma_dt.precision( 8 );

  vector<binning> binning_sigma ( 8 );
  for ( int i = 0; i < 8 ; i++ )
  {
    binning_sigma[i].setMinMaxN( s_min, s_max, 3000 );
  }
  
  interpolation23 theI23( false );
  theI23.configure( false, 1, 0.0, 1.0, "bamps_org_extended", "GBimproved", false, 1.0, log_interpol, false );

  clock_t start = clock();
  
  for( double s = s_min; s <= s_max; s += ds )
  {
    file_sigma << s << "\t";
    
    for ( double lambda = lambda_min; lambda <= lambda_max; lambda += dlambda )
    {
      bool I23_result_trustable = true;
      file_sigma << 1.0 / s * ns_casc::Ncolor * theI23.getI23( log( md2g_wo_as / s ), log( lambda * sqrt( s ) ), beta, fabs( cos_theta ), log( s ), I23_result_trustable ) << "\t";
    }
    file_sigma << endl;
  }
  file_sigma << endl << endl;


  for( double s = s_min; s <= s_max; s += ds )
  {
    file_sigma << s << "\t";
    
    for ( double lambda = lambda_min; lambda <= lambda_max; lambda += dlambda )
    {
      integrand23 theIntegrand;
      theIntegrand.set_md2_wo_as( md2g_wo_as / s );
      theIntegrand.set_lambda( lambda * sqrt( s ) );
      theIntegrand.set_cos_theta( cos_theta );
      theIntegrand.set_beta( beta );
      theIntegrand.set_m1( 0.0 );
      theIntegrand.set_m2( 0.0 );
      theIntegrand.set_s( s );
      theIntegrand.set_kappa( 1.0 );
      theIntegrand.set_formationTimeTyp( "bamps_org_extended" );
      theIntegrand.set_matrixElement23( "GBimproved" );
      theIntegrand.set_matrixElement23_22qt( false );
      theIntegrand.set_md2_counter_term_in_I23( false );
      theIntegrand.set_fudge_factor_lpm( 1.0 );

      //--------------------------------------------------------------
      // parameters needed for calls to the integration routines
      int neval;  // actual number of integrand evaluations needed
      int fail;   // 0 = desired accuracy was reached, 1 = accuracy not reached, -1 = dimension out of range
      double intResult[NCOMP_23]; // result of the integration over the unit hypercube, NCOMP_23 = #components, 1 in our case
      double error[NCOMP_23];     // presumed absolute error of integral
      double prob[NCOMP_23];      // xi^2 probability that error is NOT a reliable estimate of the true integration error
      //--------------------------------------------------------------

      //--------------------------------------------------------------
      // create the functionoid that handles the integration
      // will be called later with: integrate( theIntegrand, neval, fail, intResult, error, prob );
      // the integration routine is chosen via precompiler defined switches
      integrate_nr_vegas integrate( NDIM_23 );
      //--------------------------------------------------------------

      integrate( theIntegrand, neval, fail, intResult, error, prob );

      double I23_gg_ggg = ( 9.0 / ( 2.0 * M_PI ) * intResult[0] );

      if ( I23_gg_ggg < 0.0 )
        file_sigma << 0.0 << "\t";
      else
        file_sigma << 1.0 / s * ns_casc::Ncolor * I23_gg_ggg << "\t";
    }
    file_sigma << endl;
  }
  file_sigma << endl << endl;

  clock_t end = clock();
  cout << "time = " << (end-start)/CLOCKS_PER_SEC << " ms \t time / sampling = " << (end-start)/CLOCKS_PER_SEC/static_cast<double>(n_s * n_lambda) << " ms" << endl;  
}


int main(int argc, char **argv)
{

  try
  {
    coupling::set_isRunning( running_coupling );
    coupling::set_maxRunningCoupling( maxRunningCoupling );  
    
    ParticlePrototype::set_N_light_flavor( Nf_light );
    ParticlePrototype::set_N_heavy_flavor( Nf_heavy );
  
    ParticlePrototype::setCharmMass( 1.3 );
    ParticlePrototype::setBottomMass( 4.6 );
   
    totalcs_output();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
