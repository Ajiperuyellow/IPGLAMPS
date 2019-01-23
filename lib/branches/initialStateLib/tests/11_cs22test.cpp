//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/11_cs22test.cpp $
//$LastChangedDate: 2017-04-06 15:21:43 +0200 (Do, 06. Apr 2017) $
//$LastChangedRevision: 2558 $
//$LastChangedBy: greif $
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
const double maxRunningCoupling = 1.0;

void totalcs_output( const double T )
{
  // Boltzmann definition
  const double md2g_wo_as = ( ns_casc::Ncolor + ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
  const double md2q_wo_as = 16 / ( 3 * M_PI ) * pow( T, 2.0 ); //GeV^2
  
  stringstream filename_sigma;
  filename_sigma << "output/sigma_md2g" << md2g_wo_as * coupling::get_constant_coupling() << "_md2q" << md2q_wo_as * coupling::get_constant_coupling() << "_runCoup" << coupling::isRunning() << "_asMax" << maxRunningCoupling << ".f1";
  ofstream file_sigma( filename_sigma.str().c_str() );

  int n = 1000;
  double s_min = 0.1;
  double s_max = 1000.0;
  double ds = (s_max - s_min) / n;

    cout << "T" << "  " << "md2g_wo_as" << "  " << "md2q_wo_as" << endl;
  cout << T << "  " << md2g_wo_as << "  " << md2q_wo_as << endl;
  
  interpolation22 theI22;
  theI22.configure( running_coupling, Nf_light, Nf_heavy, ParticlePrototype::Mcharm, ParticlePrototype::Mbottom, maxRunningCoupling );

  //  clock_t start = clock();
  for( double s = s_min; s <= s_max; s += ds )
  {
    xsection_qg_qg crossObj1( s, md2g_wo_as, md2q_wo_as, &theI22 );
    xsection_qqdash_qqdash crossObj2( s, md2g_wo_as, md2q_wo_as, &theI22 );
    xsection_qqbar_qqbarDash crossObj3( s, md2g_wo_as, md2q_wo_as, &theI22 );
    xsection_qq_qq crossObj4( s, md2g_wo_as, md2q_wo_as, &theI22 );
    xsection_gg_gg crossObj5( s, md2g_wo_as, md2q_wo_as, &theI22 );
    xsection_qqbar_qqbar crossObj6( s, md2g_wo_as, md2q_wo_as, &theI22 );
    xsection_gg_qqbar crossObj7( s, md2g_wo_as, md2q_wo_as, &theI22 );
    xsection_qqbar_gg crossObj8( s, md2g_wo_as, md2q_wo_as, &theI22 );
    
    file_sigma << s << "\t";
    file_sigma << crossObj1.totalCrossSection() << "\t";
    file_sigma << crossObj2.totalCrossSection() << "\t";
    file_sigma << crossObj3.totalCrossSection() << "\t";
    file_sigma << crossObj4.totalCrossSection() << "\t";
    file_sigma << crossObj5.totalCrossSection() << "\t";
    file_sigma << crossObj6.totalCrossSection() << "\t";
    file_sigma << crossObj7.totalCrossSection() << "\t";
    file_sigma << crossObj8.totalCrossSection() << "\t";
    file_sigma << endl;
  }
  file_sigma << endl << endl;
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

    for( double T = 0.2; T <= 1.0; T += 0.2 )
    {
      totalcs_output( T );
    }
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
