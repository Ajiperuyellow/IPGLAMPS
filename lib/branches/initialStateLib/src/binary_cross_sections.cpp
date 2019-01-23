//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binary_cross_sections.cpp $
//$LastChangedDate: 2016-09-21 12:23:47 +0200 (Mi, 21. Sep 2016) $
//$LastChangedRevision: 2435 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#include <math.h>
#include "binary_cross_sections.h"
#include "globalsettings.h"
#include "particleprototype.h"
#include "interpolation22.h"
#include "FPT_compare.h"
#include "coupling.h"

using std::cout;
using std::endl;


/**
 * Get the cross section for the process g+g -> Q+Qbar
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_gg_qqbar::totalCrossSection() const
{
  double cs;
  
  // heavy quarks can only be produced if their flavor is switched on
  if( ( ParticlePrototype::N_heavy_flavor < 2 && flavor == bottom ) || ( ParticlePrototype::N_heavy_flavor < 1 && flavor == charm ) || ( ParticlePrototype::N_light_flavor == 0 && flavor == light_quark ) )
  {
    cs = 0.0;
  }
  else if( s <= 4.0 * pow( M_Q , 2.0 ) ) // 2 times quark mass required to produce QQbar -> squared
  {
    cs = 0.0; 
  }
  else
  { 

    if ( flavor == light_quark )
      cs = ParticlePrototype::N_light_flavor * theI22->get_cs( 222, s, md2g_wo_as, md2q_wo_as );
    else if( flavor == charm )
      cs = theI22->get_cs( 2210, s, md2g_wo_as, md2q_wo_as );
    else if( flavor == bottom )
      cs = theI22->get_cs( 2220, s, md2g_wo_as, md2q_wo_as );
    else
    {
      std::string errMsg = "Wrong heavy flavor in xsection_gg_QQbar.";
      throw eBinCS_error( errMsg );
    }
      
    
    if( backreaction )
    {
      double chi = pow((1.0-4.0*pow(M_Q,2.0)/s),0.5);
      cs = cs * 32.0/9.0/pow(chi,2.0); // cs_ccb->gg = 1/2 * 64/9 * 1/chi^2 cs_gg->ccb
      
      if( flavor == light_quark )
        cs = cs / ParticlePrototype::N_light_flavor;
    }
    
    cs = Kfactor * cs; // K factor
  }
   
  return cs;
}


/**
 * Get the mandelstam t for gg->QQb via the metropolis algorithm
 *
 * @return Mandelstam t
 */
double xsection_gg_qqbar::get_mandelstam_t() const
{ 
  double t;
  double t_new;
  double g = 0;
  double r;
  
  const double M = M_Q; // heavy quark mass in GeV
  const double md2q = md2q_wo_as * coupling::get_constant_coupling(); // quark debye mass with constant coupling
  
  const double chi = sqrt(1.0-4.0*pow(M,2.0)/s);
  // the range in which the variable t needs to be sampled
  // range is different for massless quarks because sampling is done only for half of t range and afterwards mirrored
  double t_min;
  double t_max;
  if ( M != 0.0 )
  {
    t_min = pow(M,2.0)-s/2.0*(1.0+chi);
    t_max = pow(M,2.0)-s/2.0*(1.0-chi);
  }
  else
  {
    t_min = -s/2;
    t_max = 0.0;
  }
  
  // select initial values of t
  do
  {
    r = ran2();
    // different comparison function for massive and massless partons. For heavy quarks uniform distribution is okay since their production rate is small, but for light quarks some more concrete function is necessary.
    if ( M != 0.0 )
    {
      t = t_min + ran2()*(t_max - t_min);
    }
    else
    {
      t = ( md2q * s * ( r - 1 ) )/(2*md2q + s * r); // inverse function of integrated comparison function based on 1/t^2-term.
    }
    
    g = getMatrixElement(t);
  } while(FPT_COMP_E(g,0.0));
  
  
  // do markov_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for(int i=0; i < markov_steps; i++)
  {
    do
    {
      r = ran2();
      
      if ( M != 0.0 )
      {
        t_new = t_min + ran2()*(t_max - t_min);          // propose new u using a uniform distribution over the entire range
      }
      else
      {
        t_new = ( md2q * s * ( r - 1 ) )/(2*md2q + s * r); // inverse function of integrated comparison function based on 1/t^2-term.
      }
    } while( t_new < t_min || t_new > t_max);
    
    double g_new = getMatrixElement(t_new);              // calculate the matrix element at the proposed point
    
    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    if ( M == 0.0 )
    {
      ratio = ratio * ( 1.0/pow(t-md2q,2.0) ) / ( 1.0/pow(t_new-md2q,2.0) ); // necessary if one does not use a symmetric propose fct. for t_new
    }
        
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      t = t_new;
      g = g_new;
    }
  }

  // sampled value for massless partons is mirrored around -s/2 because sampling is done for half of the t range.
  if ( M == 0.0 )
  {
    bool reflect = false;
    if ( ran2() < 0.5 )
    {
      reflect = true;
    }
  
    t = t - reflect * ( 2 * t + s );
  }

  return t;
}


/**
 * Get matrix element or d sigma / d t (without the prefactor) of gg->QQb for given t, cf. Combridge 79
 *
 * @param[in] t Mandelstam t
 * @return matrix element of gg->QQb
 */
double xsection_gg_qqbar::getMatrixElement(const double t) const
{
  const double M2 = M_Q * M_Q;
  
  //  const double sM2 = s-M2;
  const double u = 2.0 * M2 - s - t;
  const double M2u = s+t-M2; // M2 - u
  const double M2t = M2 - t; // M2 - t
  const double smd2g = s + md2g_wo_as * coupling::get_coupling( s );
  const double M2umd2q = M2u + md2q_wo_as * coupling::get_coupling( -( M2 - u ) );
  const double M2tmd2q = M2t + md2q_wo_as * coupling::get_coupling( -( M2 - t ) );
  
  const double ss_channel = 12.0 / pow( smd2g, 2.0 ) * M2t * M2u;
  const double tt_channel = 8.0 / 3.0 * ( M2t * M2u - 2.0 * M2 * (M2 + t ) ) / pow( M2tmd2q , 2.0 );
  const double uu_channel = 8.0 / 3.0 * ( M2t * M2u - 2.0 * M2 * (M2 + u ) ) / pow( M2umd2q , 2.0 );
  const double tu_channel = - 2.0 / 3.0 * M2 * ( s - 4.0 * M2 ) / M2tmd2q / M2umd2q;
  const double st_channel = - 6.0 * ( M2t * M2u + M2 * ( u - t ) ) / smd2g / M2tmd2q;
  const double su_channel = - 6.0 * ( M2t * M2u + M2 * ( t - u ) ) / smd2g / M2umd2q;
  
  const double as_t = coupling::get_coupling( -( M2 - t ) ); // t channel
  const double as_s = coupling::get_coupling( s ); // s channel
  const double as_u = coupling::get_coupling( -( M2 - u ) ); // u channel
  
  const double M = pow( M_PI , 2.0 ) * ( pow( as_s , 2.0) * ss_channel + pow( as_t , 2.0) * tt_channel + pow( as_u , 2.0) * uu_channel + as_t * as_u * tu_channel + as_s * as_t * st_channel + as_s * as_u * su_channel );

  return M;
}





/**
 * Get the cross section for the process q+qbar -> Q+Qbar
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_qqbar_qqbarDash::totalCrossSection() const
{
  double cs;
  
  // heavy quarks can only be produced if their flavor is switched on
  if( ( ParticlePrototype::N_heavy_flavor < 2 && flavor == bottom ) || ( ParticlePrototype::N_heavy_flavor < 1 && flavor == charm ) || ParticlePrototype::N_light_flavor == 0 )
  {
    cs = 0.0;
  }
  else if( s <= 4.0 * pow( M_Q , 2.0 ) ) // 2 times charm mass required to produce QQbar -> squared
  {
    cs = 0.0; 
  }
  else
  { 
    if ( flavor == light_quark )
      cs = (ParticlePrototype::N_light_flavor - 1) * theI22->get_cs( 225, s, md2g_wo_as );
    else if( flavor == charm )
      cs = theI22->get_cs( 2211, s, md2g_wo_as );
    else if( flavor == bottom )
      cs = theI22->get_cs( 2221, s, md2g_wo_as );
    else
    {
      std::string errMsg = "Wrong heavy flavor in xsection_qqbar_QQbar.";
      throw eBinCS_error( errMsg );
    }
      
    
    if( backreaction )
    {
      double chi = pow((1.0-4.0*pow(M_Q,2.0)/s),0.5);
      cs = cs / pow(chi,2.0); // cs_ccb->qqb = N_f/chi^2 cs_qqb->ccb
      
      if( flavor == charm || flavor == bottom )
        cs = cs * ParticlePrototype::N_light_flavor; // more number of final states for backreaction
    }
    
    cs = Kfactor * cs; // K factor
  }
   
  return cs;
}


/**
 * Get the mandelstam t for q+qbar -> Q+Qbar via the metropolis algorithm
 *
 * @return Mandelstam t
 */
double xsection_qqbar_qqbarDash::get_mandelstam_t() const
{ 
  double t;
  double t_new;
  double g = 0;
  //  double r;
  
  const double M = M_Q; // heavy quark mass in GeV

  const double chi = sqrt(1.0-4.0*pow(M,2.0)/s);
  // the range in which the variable t needs to be sampled
  const double t_min = pow(M,2.0)-s/2.0*(1.0+chi);
  const double t_max = pow(M,2.0)-s/2.0*(1.0-chi);
  
  // select initial values of t
  do
  {
    t = t_min + ran2()*(t_max - t_min);
    //     r = ran2();
    //     t = md2g_kappa*c*(r-1.0)/(r*c+md2g_kappa*s); // inverse function of integrated comparison function 1/t^2. Thus, samples t according to 1/t^2 (first term in matrix element)
  
    g = getMatrixElement(t);
  } while(FPT_COMP_E(g,0.0));
  
  
  // do markov_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for(int i=0; i<markov_steps; i++)
  {
    do
    {
      t_new = t_min + ran2()*(t_max - t_min);          // propose new u using a uniform distribution over the entire range
      //       r = ran2();
      //       t_new = md2g_kappa*c*(r-1.0)/(r*c+md2g_kappa*s); // inverse function of integrated comparison function 1/t^2
      
    } while( t_new < t_min || t_new > t_max);
    
    double g_new = getMatrixElement(t_new);              // calculate the matrix element at the proposed point
    
    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    
    //     ratio = ratio * ( 1.0/pow(t-md2g_kappa,2.0) ) / ( 1.0/pow(t_new-md2g_kappa,2.0) ); // necessary if one does not use a symmetric propose fct. for t_new
    
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      t = t_new;
      g = g_new;
    }
  }
 
  return t;
}


/**
 * Get matrix element or d sigma / d t (without the prefactor) of q+qbar -> Q+Qbar for given t, cf. Combridge 79
 *
 * @param[in] t Mandelstam t
 * @return matrix element of q+qbar -> Q+Qbar
 */
double xsection_qqbar_qqbarDash::getMatrixElement(const double t) const
{
  const double M2 = M_Q * M_Q;
  
  const double M2u = s+t-M2; // M2 - u
  const double M2t = M2 - t; // M2 - t
  
  const double smd2g = s + md2g_wo_as * coupling::get_coupling( s );
  const double as_s = coupling::get_coupling( s ); // s channel
  
  const double M = pow( M_PI , 2.0 ) * pow( as_s , 2.0) * 64.0 / 9.0 / pow( smd2g , 2.0 ) * ( pow( M2t , 2.0 ) + pow( M2u , 2.0 ) + 2.0 * M2 * s );

  return M;
}








/**
 * Get the cross section for running or constant alpha_s and matrix element for gQ -> gQ, gQbar -> gQbar
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_qg_qg::totalCrossSection() const
{
  // const double cs_gQ_max = 1000.0 * 2.568; // 1000 mb 
  const double cs_gQ_max = 500.0 * 2.568; //  !!
  //   const double cs_gQ_max = 200.0 * 2.568;

  double cs;
  
  const double md2g_wo_as_kappa = kappa_gQgQ * md2g_wo_as;
  const double md2q_wo_as_kappa = kappa_gQgQ * md2q_wo_as;

  if(isConstantCrossSecGQ)
    return constantCrossSecValueGQ * 2.568; // 1/GeV^2
  
  if ( flavor == light_quark )
    cs = theI22->get_cs( 223 , s, md2g_wo_as_kappa, md2q_wo_as_kappa );
  else if( flavor == charm )
    cs = theI22->get_cs( 2214 , s, md2g_wo_as_kappa, md2q_wo_as_kappa );
  else if( flavor == bottom )
    cs = theI22->get_cs( 2224 , s, md2g_wo_as_kappa, md2q_wo_as_kappa );
  else
  {
    std::string errMsg = "Wrong heavy flavor in xsection_gg_QQbar.";
    throw eBinCS_error( errMsg );
  }

  cs = Kfactor * cs; // K factor
  
  if( ( flavor == charm || flavor == bottom ) && cs >= cs_gQ_max)
    return cs_gQ_max;
  else
    return cs;
}


/**
 * Get the mandelstam t for gQ->gQ with the metropolis algoritm according to function getMatrixElement_GQToGQ(t)
 *
 * @return Mandelstam t
 */
double xsection_qg_qg::get_mandelstam_t() const
{
  double t;
  double t_new;
  double g = 0;
  double r;
  
  // kappa * md2g with a constant coupling to get an estimate
  const double md2g_kappa = kappa_gQgQ * md2g_wo_as * coupling::get_constant_coupling();
    
  const double c = pow( s - pow( M_Q , 2.0 ) , 2.0 );
  
  // the range in which the variable t needs to be sampled
  const double t_min = -c/s;
  const double t_max = 0.0;
  
  if( isotropicCrossSecGQ )
  {
    double costheta = 2.0 * ran2() - 1.0;
    t = c / 2.0 / s * ( costheta - 1.0 );
    return t;
  }

  // select initial values of t
  do
  {
    //     t = t_min + ran2()*(t_max - t_min);
    r = ran2();
    t = md2g_kappa*c*(r-1.0)/(r*c+md2g_kappa*s); // inverse function of integrated comparison function 1/t^2. Thus, samples t according to 1/t^2 (first term in matrix element)
  
    g = getMatrixElement(t);
  } while(FPT_COMP_E(g,0.0));
  
  
  // do markov_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for(int i=0; i<markov_steps; i++)
  {
    do
    {
      //       t_new = t_min + ran2()*(t_max - t_min);          // propose new u using a uniform distribution over the entire range
      r = ran2();
      t_new = md2g_kappa*c*(r-1.0)/(r*c+md2g_kappa*s); // inverse function of integrated comparison function 1/t^2
      
    } while( t_new < t_min || t_new > t_max);
    
    double g_new = getMatrixElement(t_new);              // calculate the matrix element at the proposed point
    
    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    
    ratio = ratio * ( 1.0/pow(t-md2g_kappa,2.0) ) / ( 1.0/pow(t_new-md2g_kappa,2.0) ); // necessary if one does not use a symmetric propose fct. for t_new
    
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      t = t_new;
      g = g_new;
    }
  }
 
  return t;
}


/**
 * Get matrix element or d sigma / d t (without the prefactor 1 / (16.0 * M_PI * pow(s-pow(Mcharm,2.0),2.0) ) ) of gQ->gQ for given t, cf. Combridge 79
 * no approximations, just Debye screened in t channel
 *
 * @param[in] t Mandelstam t
 * @return matrix element of gQ->gQ
 */
double xsection_qg_qg::getMatrixElement(const double t) const
{
  const double M2 = M_Q * M_Q;
  
  double md2g_wo_as_kappa = kappa_gQgQ * md2g_wo_as;
  double md2q_wo_as_kappa = kappa_gQgQ * md2q_wo_as;

  // screen 1/t in t channel, then write down matrix element, do not cancel any t -> according to Andre Peshier the best approach
  const double sM2 = s-M2;
  const double sM2md2q = sM2 + md2q_wo_as_kappa * coupling::get_coupling( s-M2 ); // s - M2
  const double stM2 = s+t-M2; // M2 - u
  const double stM2md2q = stM2 + md2q_wo_as_kappa * coupling::get_coupling( -(s+t-M2) ); // M2 - u
  const double tmd2g = t - md2g_wo_as_kappa * coupling::get_coupling( t );
  const double sM2stM2 = sM2*stM2;
  const double a = 32.0 * sM2stM2 / pow(tmd2g,2.0);
  const double b = 64.0/9.0 * (sM2stM2+2.0*M2*(s+M2)) / pow( sM2md2q, 2.0 );
  const double c = 64.0/9.0 * (sM2stM2+2.0*M2*(3.0*M2-s-t)) / pow(stM2md2q,2.0);
  const double d = 16.0/9.0 * M2*(4.0*M2-t) / sM2md2q / stM2md2q;
  const double e = 16.0 * (sM2stM2+M2*(2.0*s+t-2.0*M2)) / tmd2g / sM2md2q;
  const double f = -16.0 * (sM2stM2-M2*(2.0*s+t-2.0*M2)) / tmd2g / stM2md2q;

  const double as_t = coupling::get_coupling( t ); // t channel
  const double as_s = coupling::get_coupling( s-M2 ); // s channel
  const double as_u = coupling::get_coupling( -(s+t-M2) ); // u channel
  
  const double M = pow(as_t,2.0) * a + pow(as_s,2.0) * b + pow(as_u,2.0) * c
    + as_s * as_u * d + as_t * as_s * e + as_t * as_u * f;

  return pow( M_PI , 2.0 ) * M;
}



/**
 * Get the cross section for running or constant alpha_s and matrix element for qQ -> qQ, qQbar -> qQbar, etc.
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_qqdash_qqdash::totalCrossSection() const
{
  // const double cs_qQ_max = 1000.0 * 2.568; // 1000 mb 
  const double cs_qQ_max = 500.0 * 2.568; //  !!
  //   const double cs_qQ_max = 200.0 * 2.568;

  double cs;
  
  const double md2g_wo_as_kappa = kappa_qQqQ * md2g_wo_as;

  if(isConstantCrossSecqQ)
    return constantCrossSecValueqQ * 2.568; // 1/GeV^2
  
  if ( flavor == light_quark )
    cs = theI22->get_cs( 228 , s, md2g_wo_as_kappa );
  else if( flavor == charm )
    cs = theI22->get_cs( 2215 , s, md2g_wo_as_kappa );
  else if( flavor == bottom )
    cs = theI22->get_cs( 2225 , s, md2g_wo_as_kappa );
  else
  {
    std::string errMsg = "Wrong heavy flavor in xsection_gg_QQbar.";
    throw eBinCS_error( errMsg );
  }

  cs = Kfactor * cs; // K factor
  
  if( ( flavor == charm || flavor == bottom ) && cs >= cs_qQ_max)
    return cs_qQ_max;
  else
    return cs;
}


/**
 * Get the mandelstam t for qQ->qQ with the metropolis algoritm according to function getMatrixElement_qQToqQ(t)
 *
 * @return Mandelstam t
 */
double xsection_qqdash_qqdash::get_mandelstam_t() const
{
  double t;
  double t_new;
  double g = 0;
  double r;
  
  // kappa * md2g with a constant coupling to get an estimate
  const double md2g_kappa = kappa_qQqQ * md2g_wo_as * coupling::get_constant_coupling();
    
  const double c = pow( s - pow( M_Q , 2.0 ) , 2.0 );
  
  // the range in which the variable t needs to be sampled
  const double t_min = -c/s;
  const double t_max = 0.0;
  
  if( isotropicCrossSecqQ )
  {
    double costheta = 2.0 * ran2() - 1.0;
    t = c / 2.0 / s * ( costheta - 1.0 );
    return t;
  }

  // select initial values of t
  do
  {
    //     t = t_min + ran2()*(t_max - t_min);
    r = ran2();
    t = md2g_kappa*c*(r-1.0)/(r*c+md2g_kappa*s); // inverse function of integrated comparison function 1/t^2. Thus, samples t according to 1/t^2 (first term in matrix element)
  
    g = getMatrixElement(t);
  } while(FPT_COMP_E(g,0.0));
  
  
  // do markov_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for(int i=0; i<markov_steps; i++)
  {
    do
    {
      //       t_new = t_min + ran2()*(t_max - t_min);          // propose new u using a uniform distribution over the entire range
      r = ran2();
      t_new = md2g_kappa*c*(r-1.0)/(r*c+md2g_kappa*s); // inverse function of integrated comparison function 1/t^2
      
    } while( t_new < t_min || t_new > t_max);
    
    double g_new = getMatrixElement(t_new);              // calculate the matrix element at the proposed point
    
    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    
    ratio = ratio * ( 1.0/pow(t-md2g_kappa,2.0) ) / ( 1.0/pow(t_new-md2g_kappa,2.0) ); // necessary if one does not use a symmetric propose fct. for t_new
    
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      t = t_new;
      g = g_new;
    }
  }
 
  return t;
}


/**
 * Get matrix element or d sigma / d t (without the prefactor 1 / (16.0 * M_PI * pow(s-pow(Mcharm,2.0),2.0) ) ) of qQ->qQ for given t, cf. Combridge 79
 * no approximations, just Debye screened in t channel
 *
 * @param[in] t Mandelstam t
 * @return matrix element of qQ->qQ
 */
double xsection_qqdash_qqdash::getMatrixElement(const double t ) const
{
  const double M2 = M_Q * M_Q;
  
  double md2g_wo_as_kappa = kappa_qQqQ * md2g_wo_as;
  double md2g_kappa= md2g_wo_as_kappa * coupling::get_coupling( t ); // running Debye mass

  // screen 1/t in t channel, then write down matrix element, do not cancel any t -> according to Andre Peshier the best approach
  const double sM2 = s-M2;
  const double stM2 = s+t-M2; // M2 - u
  const double tmd2g = t-md2g_kappa;
  const double as_t = coupling::get_coupling( t ); // t channel
  const double M = 64.0 / 9.0 * pow( M_PI , 2.0 ) * pow( as_t , 2.0 ) * ( pow( stM2 , 2.0 ) + pow( sM2 , 2.0 ) + 2.0 * M2 * t ) / pow( tmd2g , 2.0 );

  return M;
}




/**
 * Get the cross section for the process g+Psi -> Q+Qbar
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_gPsi_QQbar::totalCrossSection() const
{
  double cs;
  
  // heavy quarks can only be produced if their flavor is switched on
  if( ( ParticlePrototype::N_heavy_flavor < 2 && heavy_quark_flavor == bottom ) || ( ParticlePrototype::N_heavy_flavor < 1 && heavy_quark_flavor == charm )  || ( ParticlePrototype::N_psi_states < 1 && psi_flavor == jpsi ) )
  {
    cs = 0.0;
  }
  else if( s < 4.0  *pow( M_Q , 2.0 ) ) // 2 times charm mass required to produce QQbar -> squared
  {
    cs = 0.0; 
  }
  else
  {
    const double w_0 = ( 4.0 * pow( M_Q , 2.0 ) - pow( M_psi , 2.0 ) ) / 2.0 / M_psi;
    const double A_0 = pow( 2.0 , 11.0 ) * M_PI / 27.0 / sqrt( pow( M_Q , 3.0 ) * e_psi );
    const double w = ( s - pow( M_psi , 2.0 ) ) / 2.0 / M_psi;

    cs = A_0 * pow( w/w_0 - 1.0 , 1.5 ) / pow( w/w_0 , 5.0 );

    if( isConstantCrossSec )
      cs = constantCrossSecValue * 2.568; // 1/GeV^2
      
    if( backreaction )
    {
      // If the temperature of the medium is larger than the dissociation temperature Td, Jpsi cannot be produced
      if( temperature >= dissociation_temperature ) 
      {
        return 0.0;
      }
      cs = cs * 4.0/3.0 * pow( s - pow( M_psi , 2.0 ) , 2.0 ) / s / ( s - 4.0 * pow( M_Q , 2.0 ) ); // detailed balance
    }
  }

  return cs;
}


/**
 * Get the mandelstam t for g+Psi -> Q+Qbar which gives an isotropical distribution
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_gPsi_QQbar::get_mandelstam_t() const
{
  // masses which are actually used in the simulation, the heavy quark mass M_Q for the cross section is an artificial value which is usable only for the cross section, but not for the kinematics
  const double M_psi_simulation = ParticlePrototype::getMass( psi_flavor );
  const double M_Q_simulation = ParticlePrototype::getMass( heavy_quark_flavor );
  
  // angle does not matter, just interested in yield -> isotrop
  double costheta = 2.0 * ran2() - 1.0;
  const double chi = sqrt( 1.0 - 4.0 * pow( M_Q_simulation , 2.0 ) / s );
  
  return ( pow( M_Q_simulation , 2.0 ) - ( s - pow( M_psi_simulation , 2.0 ) ) / 2.0 * ( 1.0 - chi * costheta ) );
};


/**
 * Get the cross section for the process g+g -> g+g
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_gg_gg::totalCrossSection() const
{
  double cs;
  
  const double md2g_wo_as_kappa = kappa * md2g_wo_as;

  if(isConstantCrossSec)
    return constantCrossSecValue * 2.568; // 1/GeV^2
  
  cs = theI22->get_cs( 221, s, md2g_wo_as_kappa );
  
  cs = Kfactor * cs; // K factor
  
  return cs;
}


/**
 * Get the mandelstam t for gg->gg via the metropolis algorithm
 *
 * @return Mandelstam t
 */
double xsection_gg_gg::get_mandelstam_t() const
{ 
  double t;
  double t_new;
  double g = 0;
  double r;
  
  double md2g = md2g_wo_as * coupling::get_constant_coupling(); // gluon debye mass with constant coupling
  
  // the range in which the variable t needs to be sampled
  // because distribution is symmetric around -s/2 it is possible to sample in range from -s/2 to 0 and mirror afterwards with probability 0.5 value around -s/2
  const double t_min = -s/2;
  const double t_max = 0;
  
  // select initial values of t
  do
  {
    r = ran2();
    t = ( md2g * s * ( r - 1 ) ) / (2 * md2g + s * r); // inverse function of integrated comparison function based on 1/t^2-term. Thus, samples t according to 1/t^2 (first term in matrix element)
  
    g = getMatrixElement(t);
  } while(FPT_COMP_E(g,0.0));
  
  
  // do markov_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for(int i=0; i<markov_steps; i++)
  {
    do
    {
      r = ran2();
      t_new = ( md2g * s * ( r - 1 ) ) / (2 * md2g + s * r);
     
    } while( t_new < t_min || t_new > t_max);
    
    double g_new = getMatrixElement(t_new);              // calculate the matrix element at the proposed point
    
    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    
    ratio = ratio * ( 1.0/pow(t-md2g,2.0) ) / ( 1.0/pow(t_new-md2g,2.0) ); // necessary if one does not use a symmetric propose fct. for t_new
    
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      t = t_new;
      g = g_new;
    }
  }
 
  // sampled value is mirrored around -s/2 because sampling is done for half of the t range.
  bool reflect = false;
  if ( ran2() < 0.5 )
  {
    reflect = true;
  }
  
  t = t - reflect * ( 2 * t + s );

  return t;
}


/**
 * Get matrix element or d sigma / d t (without the prefactor) of gg->gg for given t, cf. Cutlers,Sivers 1978
 *
 * @param[in] t Mandelstam t
 * @return matrix element of gg->gg
 */
double xsection_gg_gg::getMatrixElement(const double t) const
{
  const double md2g_wo_as_kappa = kappa * md2g_wo_as;
  
  const double u = - s - t;
  const double tmd2g = t - md2g_wo_as_kappa * coupling::get_coupling( t );
  const double smd2g = s + md2g_wo_as_kappa * coupling::get_coupling( s );
  const double umd2g = u - md2g_wo_as_kappa * coupling::get_coupling( u );
  
  //   const double tt_channel = 18.0 * ( 17.0 / 2.0 - 4.0 * u * s / pow( tmd2g, 2.0 ) );
  //   const double uu_channel = 18.0 * ( 17.0 / 2.0 - 4.0 * t * s / pow( umd2g, 2.0 ) );
  //   const double ss_channel = 18.0 * ( 17.0 / 2.0 - 4.0 * u * t / pow( smd2g, 2.0 ) );
  //   const double fourgluon_channel = 18.0 * 27.0;
  //   
  //   const double tu_channel = 9.0 * ( 15.0 - pow( s , 2.0 ) / tmd2g / umd2g );
  //   const double ts_channel = 9.0 * ( 15.0 - pow( u , 2.0 ) / tmd2g / smd2g );
  //   const double us_channel = 9.0 * ( 15.0 - pow( t , 2.0 ) / smd2g / umd2g );
  //   const double tus_fourgluon_channel = 54.0 * ( -81.0 / 4.0 );
  //   const double error_in_result = 18.0 * ( - 3.0 / 4.0 ); // The result of Cutlers,Sivers 1978 differs by this term from the correct result, see for instance Peskin Schröder. However, according to the note in their paper shortly before the appendix this term should be 18.0 * ( - 3.0 / 2.0 )
  
  //   const double tt_channel = 18.0 * ( 17.0 / 2.0 * pow( t, 2.0) / pow( tmd2g, 2.0 ) - 4.0 * u * s / pow( tmd2g, 2.0 ) );
  //   const double uu_channel = 18.0 * ( 17.0 / 2.0 * pow( u, 2.0) / pow( umd2g, 2.0 ) - 4.0 * t * s / pow( umd2g, 2.0 ) );
  //   const double ss_channel = 18.0 * ( 17.0 / 2.0 * pow( s, 2.0) / pow( smd2g, 2.0 ) - 4.0 * u * t / pow( smd2g, 2.0 ) );
  //   const double fourgluon_channel = 18.0 * 27.0;
  //   
  //   const double tu_channel = 9.0 * ( 15.0 * t * u / tmd2g / umd2g - pow( s , 2.0 ) / tmd2g / umd2g );
  //   const double ts_channel = 9.0 * ( 15.0 * t * s / tmd2g / smd2g - pow( u , 2.0 ) / tmd2g / smd2g );
  //   const double us_channel = 9.0 * ( 15.0 * s * u / smd2g / umd2g - pow( t , 2.0 ) / smd2g / umd2g );
  //   const double tus_fourgluon_channel = 54.0 * ( -81.0 / 4.0 ) * ( s / smd2g + t / tmd2g + u / umd2g ) / 3.0;
  //   const double error_in_result = 18.0 * ( - 3.0 / 4.0 ) * ( s / smd2g + t / tmd2g + u / umd2g ) / 3.0; // The result of Cutlers,Sivers 1978 differs by this term from the correct result, see for instance Peskin Schröder. However, according to the note in their paper shortly before the appendix this term should be 18.0 * ( - 3.0 / 2.0 )
  
  //   const double tt_channel = 18.0 * ( 17.0 / 2.0 * pow( t, 2.0) / pow( tmd2g, 2.0 ) - 4.0 * u * s / pow( tmd2g, 2.0 ) );
  //   const double uu_channel = 18.0 * ( 17.0 / 2.0 * pow( u, 2.0) / pow( umd2g, 2.0 ) - 4.0 * t * s / pow( umd2g, 2.0 ) );
  //   const double ss_channel = 18.0 * ( 17.0 / 2.0 * pow( s, 2.0) / pow( smd2g, 2.0 ) - 4.0 * u * t / pow( smd2g, 2.0 ) );
  //   const double fourgluon_channel = 0;
  //   
  //   const double tu_channel = 9.0 * ( 15.0 * t * u / tmd2g / umd2g - pow( s , 2.0 ) / tmd2g / umd2g );
  //   const double ts_channel = 9.0 * ( 15.0 * t * s / tmd2g / smd2g - pow( u , 2.0 ) / tmd2g / smd2g );
  //   const double us_channel = 9.0 * ( 15.0 * s * u / smd2g / umd2g - pow( t , 2.0 ) / smd2g / umd2g );
  //   const double tus_fourgluon_channel = 0;
  //   const double error_in_result = 0; // The result of Cutlers,Sivers 1978 differs by this term from the correct result, see for instance Peskin Schröder. However, according to the note in their paper shortly before the appendix this term should be 18.0 * ( - 3.0 / 2.0 )
  
  //   const double tt_channel = 18.0 * ( - 4.0 * u * s / pow( tmd2g, 2.0 ) );
  //   const double uu_channel = 18.0 * ( - 4.0 * t * s / pow( umd2g, 2.0 ) );
  //   const double ss_channel = 18.0 * ( - 4.0 * u * t / pow( smd2g, 2.0 ) );
  //   const double fourgluon_channel = 0;
  //   
  //   const double tu_channel = 9.0 * ( - pow( s , 2.0 ) / tmd2g / umd2g );
  //   const double ts_channel = 9.0 * ( - pow( u , 2.0 ) / tmd2g / smd2g );
  //   const double us_channel = 9.0 * ( - pow( t , 2.0 ) / smd2g / umd2g );
  //   const double tus_fourgluon_channel = 0;
  //   const double error_in_result = 0; // The result of Cutlers,Sivers 1978 differs by this term from the correct result, see for instance Peskin Schröder. However, according to the note in their paper shortly before the appendix this term should be 18.0 * ( - 3.0 / 2.0 )
  
  
  
  
  const double tt_channel = 18.0 * ( - 4.0 * u * s / pow( tmd2g, 2.0 ) );
  const double uu_channel = 18.0 * ( - 4.0 * t * s / pow( umd2g, 2.0 ) );
  const double ss_channel = 18.0 * ( - 4.0 * u * t / pow( smd2g, 2.0 ) );
  const double fourgluon_channel = 18.0 * (4.0 * 3.0 );
  
  const double tu_channel = 0;
  const double ts_channel = 0;
  const double us_channel = 0;
  const double tus_fourgluon_channel = 0;
  const double error_in_result = 0; // The result of Cutlers,Sivers 1978 differs by this term from the correct result, see for instance Peskin Schröder. However, according to the note in their paper shortly before the appendix this term should be 18.0 * ( - 3.0 / 2.0 )
  
  
  
  
  
  
  
  const double as_t = coupling::get_coupling( t ); // t channel
  const double as_s = coupling::get_coupling( s ); // s channel
  const double as_u = coupling::get_coupling( u ); // u channel
  const double as_fourgluon = coupling::get_coupling( s ); // what scale? Just take s
  
  const double M = pow( M_PI, 2.0 ) * ( pow( as_s , 2.0 ) * ss_channel + pow( as_t , 2.0 ) * tt_channel + pow( as_u , 2.0 ) * uu_channel
					+ as_s * as_t * ts_channel + as_s * as_u * us_channel + as_t * as_u * tu_channel 
					+ pow( as_fourgluon , 2.0 ) * fourgluon_channel + as_fourgluon * ( ( as_t + as_u + as_s ) / 3.0 ) * tus_fourgluon_channel + as_fourgluon * ( ( as_t + as_u + as_s ) / 3.0 ) * error_in_result );
  return M;
}

/**
 * Get the cross section for the process q+q -> q+q
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_qq_qq::totalCrossSection() const
{
  double cs;
  
  const double md2g_wo_as_kappa = kappa * md2g_wo_as;

  if(isConstantCrossSec)
    return constantCrossSecValue * 2.568; // 1/GeV^2
  
  cs = theI22->get_cs( 227, s, md2g_wo_as_kappa );
  
  cs = Kfactor * cs; // K factor
  
  return cs;
}


/**
 * Get the mandelstam t for qq->qq via the metropolis algorithm
 *
 * @return Mandelstam t
 */
double xsection_qq_qq::get_mandelstam_t() const
{ 
  double t;
  double t_new;
  double g = 0;
  double r;
  
  const double M = 0.0; // heavy quark mass in GeV

  const double chi = sqrt(1.0-4.0*pow(M,2.0)/s);
  // the range in which the variable t needs to be sampled
  // because distribution is symmetric around -s/2 it is possible to sample in range from -s/2 to 0 and mirror afterwards with probability 0.5 value around -s/2
  const double t_min = ( pow(M,2.0)-s/2.0*(1.0+chi) ) / 2;
  const double t_max = pow(M,2.0)-s/2.0*(1.0-chi) ;
  
  const double md2q = md2q_wo_as * coupling::get_constant_coupling(); // quark debye mass with constant coupling
  // select initial values of t
  do
  {
    r = ran2();
    t = ( md2q * s * ( r - 1 ) ) / (2 * md2q + s * r); // inverse function of integrated comparison function based on 1/t^2-term. Thus, samples t according to 1/t^2 (first term in matrix element)
  
    g = getMatrixElement(t);
  } while(FPT_COMP_E(g,0.0));
  
  
  // do markov_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for(int i=0; i<markov_steps; i++)
  {
    do
    {
      r = ran2();
      t_new = ( md2q * s * ( r - 1 ) ) / (2 * md2q + s * r); // inverse function of integrated comparison function based on 1/t^2-term. Thus, samples t according to 1/t^2 (first term in matrix element)
      
    } while( t_new < t_min || t_new > t_max);
    
    double g_new = getMatrixElement(t_new);              // calculate the matrix element at the proposed point
    
    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    
    ratio = ratio * ( 1.0/pow(t-md2q,2.0) ) / ( 1.0/pow(t_new-md2q,2.0) ); // necessary if one does not use a symmetric propose fct. for t_new
    
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      t = t_new;
      g = g_new;
    }
  }
 
  // sampled value is mirrored around -s/2 because sampling is done for half of the t range.
  bool reflect = false;
  if ( ran2() < 0.5 )
  {
    reflect = true;
  }
  
  t = t - reflect * ( 2 * t + s );

  return t;
}


/**
 * Get matrix element or d sigma / d t (without the prefactor) of qq->qq for given t, cf. Peskin/Schröder
 *
 * @param[in] t Mandelstam t
 * @return matrix element of qq->qq
 */
double xsection_qq_qq::getMatrixElement(const double t) const
{
  const double md2g_wo_as_kappa = kappa * md2g_wo_as;
  
  const double u = - s - t;
  const double tmd2g = t - md2g_wo_as_kappa * coupling::get_coupling( t );
  const double umd2g = u - md2g_wo_as_kappa * coupling::get_coupling( u );
  const double s2 = pow( s , 2.0 );
  const double t2 = pow( t , 2.0 );
  const double u2 = pow( u , 2.0 );
  
  const double tt_channel = 64.0 / 9.0 * ( s2 + u2 ) / pow( tmd2g, 2.0 );
  const double uu_channel = 64.0 / 9.0 * ( s2 + t2 ) / pow( umd2g, 2.0 );
  const double tu_channel = 64.0 / 9.0 * ( - 2.0 / 3.0 ) * s2 / tmd2g / umd2g;
  
  const double as_t = coupling::get_coupling( t ); // t channel
  const double as_u = coupling::get_coupling( u ); // u channel
   
  const double M = pow( M_PI , 2.0 ) * ( pow( as_t , 2.0) * tt_channel + pow( as_u , 2.0) * uu_channel + as_t * as_u * tu_channel );

  return M;
}

/**
 * Get the cross section for the process q+qbar -> q+qbar
 *
 * @return cross section in units 1/GeV^2
 */
double xsection_qqbar_qqbar::totalCrossSection() const
{
  double cs;
  
  const double md2g_wo_as_kappa = kappa * md2g_wo_as;

  if(isConstantCrossSec)
    return constantCrossSecValue * 2.568; // 1/GeV^2
  
  cs = theI22->get_cs( 224, s, md2g_wo_as_kappa );
  
  cs = Kfactor * cs; // K factor
  
  return cs;
}


/**
 * Get the mandelstam t for qqbar->qqbar via the metropolis algorithm
 *
 * @return Mandelstam t
 */
double xsection_qqbar_qqbar::get_mandelstam_t() const
{
  double t;
  double t_new;
  double g = 0;
  double r;
  
  const double M = 0.0; // heavy quark mass in GeV

  const double chi = sqrt(1.0-4.0*pow(M,2.0)/s);
  // the range in which the variable t needs to be sampled
  const double t_min = pow(M,2.0)-s/2.0*(1.0+chi);
  const double t_max = pow(M,2.0)-s/2.0*(1.0-chi);
  
  const double md2g = md2g_wo_as * coupling::get_constant_coupling(); // gluon debye mass with constant coupling
    
  const double c = pow( s - pow( M , 2.0 ) , 2.0 );
  
  // select initial values of t
  do
  {
    r = ran2();
    t = md2g*c*(r-1.0)/(r*c+md2g*s); // inverse function of integrated comparison function 1/t^2. Thus, samples t according to 1/t^2 (first term in matrix element)
  
    g = getMatrixElement(t);
  } while(FPT_COMP_E(g,0.0));
  
  
  // do markov_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  // remark: number of steps in the Markov chain has to be scaled up in qqbar->qqbar process
  for(int i=0; i < markov_steps; i++)
  {
    do
    {
      r = ran2();
      t_new = md2g*c*(r-1.0)/(r*c+md2g*s); // inverse function of integrated comparison function 1/t^2
      
    } while( t_new < t_min || t_new > t_max);
    
    double g_new = getMatrixElement(t_new);              // calculate the matrix element at the proposed point
    
    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    
    ratio = ratio * ( 1.0/pow(t-md2g,2.0) ) / ( 1.0/pow(t_new-md2g,2.0) ); // necessary if one does not use a symmetric propose fct. for t_new
    
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      t = t_new;
      g = g_new;
    }
  }
 
  return t;
}


/**
 * Get matrix element or d sigma / d t (without the prefactor) of qqbar->qqbar for given t, cf. Peskin / Schröder
 *
 * @param[in] t Mandelstam t
 * @return matrix element of qqbar->qqbar
 */
double xsection_qqbar_qqbar::getMatrixElement(const double t) const
{
  const double md2g_wo_as_kappa = kappa * md2g_wo_as;
  
  const double u = -s - t;
  
  const double s2 = pow( s , 2.0 );
  const double t2 = pow( t , 2.0 );
  const double u2 = pow( u , 2.0 );
  
  const double smd2g = s + md2g_wo_as_kappa * coupling::get_coupling( s );
  const double tmd2g = t - md2g_wo_as_kappa * coupling::get_coupling( t );
  
  const double ss_channel = 64.0 / 9.0 * (t2 + u2) / pow( smd2g, 2.0 );
  const double tt_channel = 64.0 / 9.0 * (s2 + u2) / pow( tmd2g, 2.0 );
  const double st_channel = 64.0 / 9.0 * ( - 2.0 / 3.0 ) * u2 / smd2g / tmd2g;

  const double as_t = coupling::get_coupling( t ); // t channel
  const double as_s = coupling::get_coupling( s ); // s channel
  
  const double M = pow( M_PI , 2.0 ) * ( pow( as_s , 2.0) * ss_channel + pow( as_t , 2.0) * tt_channel + as_t * as_s * st_channel );
  
  return M;
}
