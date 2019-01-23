//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: svn+ssh://senzel@th.physik.uni-frankfurt.de/home/bamps/svn/lib/trunk/tests/03_tsampling.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10 Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the sampling of t in 2->2 processes based on the old outdated asConst routines.
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
const int Nf_light = 3;
const int Nf_heavy = 0;
const double T = 0.4;

void t_sampling( double s )
{
//   Boltzmann definition
  const double md2g_wo_as = ( ns_casc::Ncolor + ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
  const double md2q_wo_as = 16 / ( 3 * M_PI ) * pow( T, 2.0 ); //GeV^2

  const int n = 10000;
  const double t_min = -s;
  const double t_max = 0.0;
  const double dt = (t_max - t_min) / n;

  cout << "T" << "  " << "md2g_wo_as" << "  " << "md2q_wo_as" << "  " << "s" << endl;
  cout << T << "  " << md2g_wo_as << "  " << md2q_wo_as << "  " << s << endl;


  xsection_qg_qg_asConst crossObj1( s, md2g_wo_as, md2q_wo_as );
  xsection_qqdash_qqdash_asConst crossObj2( s, md2g_wo_as, md2q_wo_as );
  xsection_qqbar_qqbarDash_asConst crossObj3( s, md2g_wo_as, md2q_wo_as );
  xsection_qq_qq_asConst crossObj4( s, md2g_wo_as, md2q_wo_as );
  xsection_gg_gg_asConst crossObj5( s, md2g_wo_as, md2q_wo_as );
  xsection_qqbar_qqbar_asConst crossObj6( s, md2g_wo_as, md2q_wo_as );
  xsection_gg_qqbar_asConst crossObj7( s, md2g_wo_as, md2q_wo_as );
  xsection_qqbar_gg_asConst crossObj8( s, md2g_wo_as, md2q_wo_as );
  
  cout << ".......... 1 file_dsigma_dt" << endl;

  stringstream filename_dsigma_dt;
  filename_dsigma_dt << "output/dsigma_dt_asConst_s" << s << ".f1";
  ofstream file_dsigma_dt( filename_dsigma_dt.str().c_str() );
  file_dsigma_dt.precision( 8 );

  for( double t = t_min; t <= t_max; t += dt )
  {
    file_dsigma_dt << t << "\t";
    file_dsigma_dt << crossObj1.differentialCrossSection( t ) / crossObj1.totalCrossSection() << "\t";
    file_dsigma_dt << crossObj2.differentialCrossSection( t ) / crossObj2.totalCrossSection() << "\t";
    file_dsigma_dt << crossObj3.differentialCrossSection( t ) / crossObj3.totalCrossSection() << "\t";
    file_dsigma_dt << crossObj4.differentialCrossSection( t ) / crossObj4.totalCrossSection() << "\t";
    file_dsigma_dt << crossObj5.differentialCrossSection( t ) / crossObj5.totalCrossSection() << "\t";
    file_dsigma_dt << crossObj6.differentialCrossSection( t ) / crossObj6.totalCrossSection() << "\t";
    file_dsigma_dt << crossObj7.differentialCrossSection( t ) / crossObj7.totalCrossSection() << "\t";
    file_dsigma_dt << crossObj8.differentialCrossSection( t ) / crossObj8.totalCrossSection() << "\t";
    file_dsigma_dt << endl;
  }
  file_dsigma_dt << endl << endl;
  
  cout << ".......... 2 binning_dsigma_dt" << endl;

  vector<binning> binning_dsigma_dt ( 8 );
  for ( int i = 0; i < 8 ; i++ )
  {
    binning_dsigma_dt[i].setMinMaxN( t_min, t_max, 100 );
  }

  int nSamples = 1000000;
  clock_t start, end;
  start = clock();

  double t;
  for ( int i = 0; i < nSamples; i++ )
  {
    t = sample_t_hat( crossObj1 );
    binning_dsigma_dt[0].add( t );
    t = sample_t_hat( crossObj2 );
    binning_dsigma_dt[1].add( t );
    t = sample_t_hat( crossObj3 );
    binning_dsigma_dt[2].add( t );
    t = sample_t_hat( crossObj4 );
    binning_dsigma_dt[3].add( t );
    t = sample_t_hat( crossObj5 );
    binning_dsigma_dt[4].add( t );
    t = sample_t_hat( crossObj6 );
    binning_dsigma_dt[5].add( t );
    t = sample_t_hat( crossObj7 );
    binning_dsigma_dt[6].add( t );
    t = sample_t_hat( crossObj8 );
    binning_dsigma_dt[7].add( t );
  }
  end = clock();
  cout << "time = " << (end-start)/CLOCKS_PER_SEC << " ms \t time / sampling = " << (end-start)/CLOCKS_PER_SEC/static_cast<double>(nSamples)/8. << " ms" << endl;

  cout << ".......... 3 file_dsigma_dt_sampled: " <<  binning_dsigma_dt[0].getNBins() << endl;

  stringstream filename_dsigma_dt_sampled;
  filename_dsigma_dt_sampled << "output/dsigma_dt_sampled_runningCoupling" << coupling::isRunning() << "_s" << s << ".f2";
  ofstream file_dsigma_dt_sampled( filename_dsigma_dt_sampled.str().c_str() );

  for ( int i = 0; i < binning_dsigma_dt[0].getNBins(); i++ )
  {
    file_dsigma_dt_sampled << binning_dsigma_dt[0].getBinLabel( i ) << "\t";
    for ( int j = 0; j < 8; j++ )
    {
      //       file_dsigma_dt_sampled << binning_dsigma_dt[j].getBinRelative( i )/normalization_factor[j] << "\t";
      file_dsigma_dt_sampled << binning_dsigma_dt[j].getBinRelative( i ) << "\t";
    }
    file_dsigma_dt_sampled << endl;
  }
  file_dsigma_dt_sampled << endl << endl;

  cout << ".......... 4 clean up" << endl << endl;
}


int main(int argc, char **argv)
{

  try
  {
    ParticlePrototype::set_N_light_flavor( Nf_light );
    ParticlePrototype::set_N_heavy_flavor( Nf_heavy );
  
    ParticlePrototype::setCharmMass( 1.3 );
    ParticlePrototype::setBottomMass( 4.6 );

    double s[3];
    s[0] = 0.3;
    s[1] = 3.0;
    s[2] = 240.0;

    for( int i = 0; i < 3; i++ )
    {
      t_sampling( s[i] );
    }
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;
}
