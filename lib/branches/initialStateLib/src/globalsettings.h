//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/globalsettings.h $
//$LastChangedDate: 2014-11-27 20:37:27 +0100 (Do, 27. Nov 2014) $
//$LastChangedRevision: 1976 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Global settings that need to be shared across various routines and translation units
 */


#ifndef GLOBALSETTINGS_H
#define GLOBALSETTINGS_H


/** @brief Namespace for global BAMPS paramters, objects etc. */
namespace ns_casc
{
  //-- global parameters -----------------------------//
  /** @brief (lambda_QCD)^2 in GeV^2, global parameter*/
  const double lambda2 = 0.04;
  
  /** @brief number of colors, global parameter */
  const unsigned int Ncolor = 3;

  /** @brief degeneracy of gluons */
  const unsigned int gG = 2 * ( pow( Ncolor, 2 ) - 1 );

  /** @brief degeneracy of quarks */
  const unsigned int gQ = 2 * Ncolor;
  
  /** @brief K-factor for minijet sampling */
  const double K = 2.0;                

  /** @brief lower cutoff for some t-chanel processes */
  const double tcut = -0.1;
  
  /**
  * @brief auxiliary variable, denotes "a large number"
  *
  * yeah, I know that is's not really infinity, but what the heck
  */
  const double infinity = 300000;
  //--------------------------------------------------//
}
#endif
