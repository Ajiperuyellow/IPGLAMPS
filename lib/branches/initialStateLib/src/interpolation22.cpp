//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interpolation22.cpp $
//$LastChangedDate: 2014-11-20 11:25:15 +0100 (Do, 20. Nov 2014) $
//$LastChangedRevision: 1955 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <boost/lexical_cast.hpp>

#include "interpolation22.h"
#include "globalsettings.h"
#include "FPT_compare.h"

using std::string;

/**
 * Configure interpolation routine. Load the tables for all processes.
 * @param[in] asRun_arg whether coupling is running or not
 * @param[in] N_light_flav number of light flavors
 * @param[in] N_heavy_flav number of light flavors
 * @param[in] M_charm charm mass
 * @param[in] M_bottom bottom mass
 */
void interpolation22::configure( const bool asRun, const int N_light_flav, const int N_heavy_flav, const double M_charm, const double M_bottom, const double maximum_coupling, const double fixed_coupling_value )
{
  if( !asRun && !FPT_COMP_E( fixed_coupling_value, 0.3 ) )
    throw eI22_read_error( "Value for fixed coupling for I22 is unequal to 0.3.");
  
  if( N_light_flav >= 0 ) // N_light_flav = -1 if no gluons should be present, only heavy quarks
  {
    const bool is_heavy_quark_involved = false;
    I_gg_gg.loadTable( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling );
    
    if( N_light_flav > 0 )
    {
      I_gg_qqbar.loadTable( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling );
      I_qg_qg.loadTable( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling );
      I_qqbar_qqbar.loadTable( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling );
      I_qqbar_qqbarDash.loadTable( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling );
      I_qq_qq.loadTable( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling );
      I_qqdash_qqdash.loadTable( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling );
    }
  }
  
  if( N_heavy_flav > 0 )
  {
    if( asRun && !FPT_COMP_E( maximum_coupling, 1.0 ) )
    {
      throw eI22_read_error( "Maximum coupling for I22 unequal to 1.");
    }
    
    I_gg_ccbar.setParametersAndLoadTable( asRun, N_light_flav, M_charm );
    I_gc_gc.setParametersAndLoadTable( asRun, N_light_flav, M_charm );
    
    if( N_light_flav > 0 )
    {
      I_qqbar_ccbar.setParametersAndLoadTable( asRun, N_light_flav, M_charm );
      I_qc_qc.setParametersAndLoadTable( asRun, N_light_flav, M_charm );
    }

    if( N_heavy_flav > 1 )
    {
      I_gg_bbbar.setParametersAndLoadTable( asRun, N_light_flav, M_bottom );
      I_gb_gb.setParametersAndLoadTable( asRun, N_light_flav, M_bottom );
      
      if( N_light_flav > 0 )
      {
        I_qqbar_bbbar.setParametersAndLoadTable( asRun, N_light_flav, M_bottom );
        I_qb_qb.setParametersAndLoadTable( asRun, N_light_flav, M_bottom );
      }
    }
  }
}


/**
 * Returns the total cross section of a given process
 * @param[in] interactionType interaction type of the process according to interactiontype.h
 * @param[in] s Mandelstam s
 * @param[in] md2g_wo_as_arg gluon debye mass squared divided by alpha_s
 * @param[in] md2q_wo_as_arg quark debye mass squared divided by alpha_s
 * @return total cross section of process
 */
double interpolation22::get_cs( const int interactionType, const double _s, const double md2g_wo_as, const double md2q_wo_as ) const
{
  double cs;
  double s = _s;
  
  // maximum value of s in tables
  const double s_max = 2000.0;
  
  // do not interpolate for large s, but use last point. This is okay since the cross section is a flat constant in this regime. Interpolation from the last two data points to large s can lead to large deviations if these last to points have some statistical fluctuations
  if( s > s_max )
    s = s_max;
  
  // set parameters for the different processes
  switch( interactionType )
  {
    // light partons
    case 221 : //g+g -> g+g
      cs = I_gg_gg.get_cs( s, md2g_wo_as );
      break;
    case 222 : //g+g -> q+qbar
      cs = I_gg_qqbar.get_cs( s, md2g_wo_as, md2q_wo_as );
      break;
    case 223 : //g+q -> g+q
      cs = I_qg_qg.get_cs( s, md2g_wo_as, md2q_wo_as );
      break;
    case 224 : //q+qbar -> q+qbar
      cs = I_qqbar_qqbar.get_cs( s, md2g_wo_as );
      break;
    case 225 : //q+qbar -> q'+qbar'
      cs = I_qqbar_qqbarDash.get_cs( s, md2g_wo_as );
      break;
    case 227 : //q+q -> q+q
      cs = I_qq_qq.get_cs( s, md2g_wo_as );
      break;
    case 228 : //q+q'-> q+q'
      cs = I_qqdash_qqdash.get_cs( s, md2g_wo_as );
      break;

    // heavy partons
    case 2210 : //g+g -> c+cbar
      cs = I_gg_ccbar.get_cs( s, md2g_wo_as, md2q_wo_as );
      break;
    case 2211 : //q+qbar -> c+cbar
      cs = I_qqbar_ccbar.get_cs( s, md2g_wo_as );
      break;
    case 2214 : //g+c -> g+c
      cs = I_gc_gc.get_cs( s, md2g_wo_as, md2q_wo_as );
      break;
    case 2215 : //q+c-> q+c
      cs = I_qc_qc.get_cs( s, md2g_wo_as );
      break;

    case 2220 : //g+g -> b+bbar
      cs = I_gg_bbbar.get_cs( s, md2g_wo_as, md2q_wo_as );
      break;
    case 2221 : //q+qbar -> b+bbar
      cs = I_qqbar_bbbar.get_cs( s, md2g_wo_as );
      break;
    case 2224 : //g+b -> g+b
      cs = I_gb_gb.get_cs( s, md2g_wo_as, md2q_wo_as );
      break;
    case 2225 : //q+b-> q+b
      cs = I_qb_qb.get_cs( s, md2g_wo_as );
      break;
    default :
      throw eI22_read_error( "Wrong interaction type given in interpolation22::get_cs().");
      break;
  }
  
  // in case of numerical rounding errors prevent that cross section could be negative
  if( cs < 0 )
  {
    std::cout << "cs22 interpolation < 0: " << cs << "  " << interactionType << "  " << s << "  " << md2g_wo_as << "  " << md2q_wo_as << std::endl;
    cs = 0;
  }

  return cs;
}


double interpolation22_process_1d::get_cs( const double s, const double md2g_wo_as, const double md2q_wo_as ) const
{
  double cs;
  double ln_s = log( s );
  
  if( ln_s < ln_s_start )
    cs = 0.0;
  else
    cs = I22.getInterpolatedData( ln_s );
//     cs = exp( I22.getInterpolatedData( ln_s ) );

  return cs;
}


double interpolation22_process_2d::get_cs( const double s, const double md2g_wo_as, const double md2q_wo_as ) const
{
  double cs;
  double ln_s = log( s );
  double ln_md2g_wo_as = log( md2g_wo_as );
  
  // If Debye is very small and not tabularized anymore set to smallest value of the tabularized value
  if( ln_md2g_wo_as < ln_md2g_wo_as_start )
    ln_md2g_wo_as = ln_md2g_wo_as_start;
  
  if( ln_s < ln_s_start )
    cs = 0.0;
  else
    cs = I22.getInterpolatedData( ln_md2g_wo_as, ln_s );
//     cs = exp( I22.getInterpolatedData( ln_md2g_wo_as, ln_s ) );
    
  return cs;
}


double interpolation22_process_3d::get_cs( const double s, const double md2g_wo_as, const double md2q_wo_as ) const
{
  double cs;
  double ln_s = log( s );
  double ln_md2g_wo_as = log( md2g_wo_as );
  double ln_md2q_wo_as = log( md2q_wo_as );
  
  // If Debye is very small and not tabularized anymore set to smallest value of the tabularized value
  if( ln_md2g_wo_as < ln_md2g_wo_as_start )
    ln_md2g_wo_as = ln_md2g_wo_as_start;
  if( ln_md2q_wo_as < ln_md2q_wo_as_start )
    ln_md2q_wo_as = ln_md2q_wo_as_start;
  
  if( ln_s < ln_s_start )
    cs = 0.0;
  else
    cs = I22.getInterpolatedData( ln_md2q_wo_as, ln_md2g_wo_as, ln_s );
  //     cs = exp( I22.getInterpolatedData( ln_md2q_wo_as, ln_md2g_wo_as, ln_s ) );
    
  return cs;
}




void interpolation22_gg_ccbar::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_charm )
{
  ln_s_start = std::max( log( pow( 2.0 * M_charm ,2.0) ), log( 1.1 * ns_casc::lambda2 ) );
  delta_ln_s = 0.02;
  
  process_name = "gg_ccbar";

  if( M_charm == 1.3 )
  {
    n_ln_s = 285;
    process_name = process_name + "_13";
  }
  else if( M_charm == 1.5 )
  {
    n_ln_s = 271;
    process_name = process_name + "_15";
  }
  else if( M_charm == 1.87 )
  {
    n_ln_s = 249;
    process_name = process_name + "_187";
  }
  else if( M_charm == 0.1 )
  {
    n_ln_s = 541;
    process_name = process_name + "_01";
  }
  else if( M_charm == 0.0 )
  {
    n_ln_s = 537;
    process_name = process_name + "_0";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_gg_ccbar::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}


void interpolation22_gg_bbbar::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_bottom )
{
  ln_s_start = std::max( log( pow( 2.0 * M_bottom ,2.0) ), log( 1.1 * ns_casc::lambda2 ) );
  delta_ln_s = 0.01;
  
  process_name = "gg_bbbar";

  if( M_bottom == 4.6 )
  {
    n_ln_s = 317;
    process_name = process_name + "_46";
  }
  else if( M_bottom == 4.8 )
  {
    n_ln_s = 308;
    process_name = process_name + "_48";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_gg_bbbar::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}


void interpolation22_qqbar_ccbar::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_charm )
{
  ln_s_start = std::max( log( pow( 2.0 * M_charm ,2.0) ), log( 1.1 * ns_casc::lambda2 ) );
  delta_ln_s = 0.02;
  
  process_name = "qqbar_ccbar";

  if( M_charm == 1.3 )
  {
    n_ln_s = 285;
    process_name = process_name + "_13";
  }
  else if( M_charm == 1.5 )
  {
    n_ln_s = 271;
    process_name = process_name + "_15";
  }
  else if( M_charm == 1.87 )
  {
    n_ln_s = 249;
    process_name = process_name + "_187";
  }
  else if( M_charm == 0.1 )
  {
    n_ln_s = 541;
    process_name = process_name + "_01";
  }
  else if( M_charm == 0.0 )
  {
    n_ln_s = 537;
    process_name = process_name + "_0";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_qqbar_ccbar::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}


void interpolation22_qqbar_bbbar::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_bottom )
{
  ln_s_start = std::max( log( pow( 2.0 * M_bottom ,2.0) ), log( 1.1 * ns_casc::lambda2 ) ); 
  delta_ln_s = 0.01;
  
  process_name = "qqbar_bbbar";

  if( M_bottom == 4.6 )
  {
    n_ln_s = 317;
    process_name = process_name + "_46";
  }
  else if( M_bottom == 4.8 )
  {
    n_ln_s = 308;
    process_name = process_name + "_48";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_qqbar_bbbar::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}



void interpolation22_gc_gc::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_charm )
{
  ln_s_start = std::max( log( pow( M_charm , 2.0 ) + 0.001 ), log( 1.1 * ns_casc::lambda2 ) ); 
  delta_ln_s = 0.02;
  
  process_name = "gc_gc";

  if( M_charm == 1.3 )
  {
    n_ln_s = 354;
    process_name = process_name + "_13";
  }
  else if( M_charm == 1.5 )
  {
    n_ln_s = 340;
    process_name = process_name + "_15";
  }
  else if( M_charm == 1.87 )
  {
    n_ln_s = 318;
    process_name = process_name + "_187";
  }
  else if( M_charm == 0.1 )
  {
    n_ln_s = 606;
    process_name = process_name + "_01";
  }
  else if( M_charm == 0.0 )
  {
    n_ln_s = 537;
    process_name = process_name + "_0";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_gc_gc::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}


void interpolation22_gb_gb::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_bottom )
{
  ln_s_start = std::max( log( pow( M_bottom , 2.0 ) + 0.001 ), log( 1.1 * ns_casc::lambda2 ) ); 
  delta_ln_s = 0.01;
  
  process_name = "gb_gb";

  if( M_bottom == 4.6 )
  {
    n_ln_s = 455;
    process_name = process_name + "_46";
  }
  else if( M_bottom == 4.8 )
  {
    n_ln_s = 447;
    process_name = process_name + "_48";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_gb_gb::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}


void interpolation22_qc_qc::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_charm )
{
  ln_s_start = std::max( log( pow( M_charm , 2.0 ) + 0.001 ), log( 1.1 * ns_casc::lambda2 ) ); 
  delta_ln_s = 0.02;
  
  process_name = "qc_qc";

  if( M_charm == 1.3 )
  {
    n_ln_s = 354;
    process_name = process_name + "_13";
  }
  else if( M_charm == 1.5 )
  {
    n_ln_s = 340;
    process_name = process_name + "_15";
  }
  else if( M_charm == 1.87 )
  {
    n_ln_s = 318;
    process_name = process_name + "_187";
  }
  else if( M_charm == 0.1 )
  {
    n_ln_s = 606;
    process_name = process_name + "_01";
  }
  else if( M_charm == 0.0 )
  {
    n_ln_s = 537;
    process_name = process_name + "_0";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_qc_qc::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}


void interpolation22_qb_qb::setParametersAndLoadTable( const bool asRun, const int N_light_flav, const double M_bottom )
{
  ln_s_start = std::max( log( pow( M_bottom , 2.0 ) + 0.001 ), log( 1.1 * ns_casc::lambda2 ) ); 
  delta_ln_s = 0.01;
  
  process_name = "qb_qb";

  if( M_bottom == 4.6 )
  {
    n_ln_s = 455;
    process_name = process_name + "_46";
  }
  else if( M_bottom == 4.8 )
  {
    n_ln_s = 447;
    process_name = process_name + "_48";
  }
  else
  {
    throw eI22_read_error( "error in interpolation22_qb_qb::setParametersAndLoadTable.");
  }
  
  loadTable( asRun, N_light_flav, true );
}

