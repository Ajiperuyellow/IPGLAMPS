//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/04_tsampling_HQ.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the sampling of t in 2->2 processes.
*/

#include <iostream>
using namespace std;

#include "FPT_compare.h"
#include "binning.h"
#include "binary_cross_sections.h"
#include <fstream>
// #include "interpolation22.h"
#include <time.h>

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
const int Nf_heavy = 1;

void t_sampling(void)
{
  const double T = 0.4;
  // Boltzmann definition
  const double md2g_wo_as = ( ns_casc::Ncolor + ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
  const double md2q_wo_as = 16 / ( 3 * M_PI ) * pow( T, 2.0 ); //GeV^2
//   const double md2g_wo_as = 0; //GeV^2
//   const double md2q_wo_as = 0; //GeV^2
//   const double s = 18.0 * T * T;
  double s[3];
  s[0] = 0.7;
  s[1] = 3.0;
  s[2] = 240.0;
  
  xsection_generic* crossObj[2];

  stringstream filename_dsigma_dt_sampled;
  filename_dsigma_dt_sampled << "output/dsigma_dt_sampled_HQ_runningCoupling" << coupling::isRunning() << ".f2";
  stringstream filename_dsigma_dt;
  filename_dsigma_dt << "output/dsigma_dt_HQ_runningCoupling" << coupling::isRunning() << ".f1";
  ofstream file_dsigma_dt( filename_dsigma_dt.str().c_str() );
  ofstream file_dsigma_dt_sampled( filename_dsigma_dt_sampled.str().c_str() );

  for ( int si = 0; si < 3; si++ )
  {
    const double c = pow( s[si] - pow( ParticlePrototype::getMass(charm) , 2.0 ) , 2.0 );
  
  // the range in which the variable t needs to be sampled
    const double t_min = -c/s[si];
    const double t_max = 0.0;
    const int n = 10000;
    const double dt = (t_max - t_min) / n;

   cout << "T" << "  " << "md2g_wo_as" << "  " << "md2q_wo_as" << "  " << "s[si]" << endl;
    cout << T << "  " << md2g_wo_as << "  " << md2q_wo_as << "  " << s[si] << endl;
  
//   file_dsigma_dt.precision( 8 );

    vector<binning> binning_dsigma_dt ( 2 );
    for ( int i = 0; i < binning_dsigma_dt.size() ; i++ )
    {
      binning_dsigma_dt[i].setMinMaxN( t_min, t_max, 500 );
    }
  
//   interpolation22 theI22;
//   theI22.configure( running_coupling, Nf_light, Nf_heavy );
  
    crossObj[0] = new xsection_cg_cg ( s[si], md2g_wo_as, md2q_wo_as );
    crossObj[1] = new xsection_cq_cq ( s[si], md2g_wo_as, md2q_wo_as );
//     crossObj[2] = new xsection_gg_ccbar ( s[si], md2g_wo_as, md2q_wo_as );
//     crossObj[3] = new xsection_ccbar_gg ( s[si], md2g_wo_as, md2q_wo_as );
//     crossObj[4] = new xsection_ccbar_qqbar ( s[si], md2g_wo_as, md2q_wo_as );
//     crossObj[5] = new xsection_qqbar_ccbar ( s[si], md2g_wo_as, md2q_wo_as );
    
    cout << ".......... 1 file_dsigma_dt" << endl;

    for( double t = t_min; t <= t_max; t += dt )
    {
      file_dsigma_dt << t << "\t";
      for ( int i = 0; i < 2; i++ )
      {
        file_dsigma_dt << crossObj[i]->differentialCrossSection( t ) << "\t";
      }
      file_dsigma_dt << endl;
    }
    file_dsigma_dt << endl << endl;
  
    //    int nSamples = 1000000;
    int nSamples = 10000;
    double t;
    clock_t start, end;

    cout << ".......... 2 binning_dsigma_dt: " << nSamples << endl;
  
    start = clock();
    for ( int i = 0; i < nSamples; i++ )
    {
      for ( int j = 0; j < 2; j++ )
      {
        t = crossObj[j]->get_mandelstam_t();
        binning_dsigma_dt[j].add( t );
      }
    }
    end = clock();
    cout << "time = " << (end-start)/CLOCKS_PER_SEC << " ms \t time / sampling = " << (end-start)/CLOCKS_PER_SEC/static_cast<double>(nSamples)/8. << " ms" << endl;
  
    cout << ".......... 3 file_dsigma_dt_sampled: " <<  binning_dsigma_dt[0].getNBins() << endl;

    for ( int i = 0; i < binning_dsigma_dt[0].getNBins(); i++ )
    {
      file_dsigma_dt_sampled << binning_dsigma_dt[0].getBinLabel( i ) << "\t";
      for ( int j = 0; j < 2; j++ )
      {
//       file_dsigma_dt_sampled << binning_dsigma_dt[j].getBinRelative( i )/normalization_factor[j] << "\t";
        file_dsigma_dt_sampled << binning_dsigma_dt[j].getBinRelative( i ) << "\t" << binning_dsigma_dt[j].getBinRelative( i ) / crossObj[j]->differentialCrossSection( binning_dsigma_dt[0].getBinLabel( i ) ) << "\t";
      }
      file_dsigma_dt_sampled << endl;
    }
    file_dsigma_dt_sampled << endl << endl;
    
    for ( int i = 0; i < 2; i++ )
    {
      delete crossObj[i];
    }

    cout << ".......... 4 clean up" << endl << endl;
  }
}


int main(int argc, char **argv)
{

  try
  {
    coupling::set_isRunning( running_coupling );
  
    ParticlePrototype::set_N_light_flavor( Nf_light );
    ParticlePrototype::set_N_heavy_flavor( Nf_heavy );
  
    ParticlePrototype::setCharmMass( 1.3 );
    ParticlePrototype::setBottomMass( 4.6 );
   
    t_sampling();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
