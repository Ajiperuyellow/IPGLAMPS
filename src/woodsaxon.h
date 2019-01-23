//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/woodsaxon.h $
//$LastChangedDate: 2012-01-27 15:25:24 +0100 (金, 27  1月 2012) $
//$LastChangedRevision: 351 $
//$LastChangedBy: fochler $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



#ifndef WOODSAXON_H
#define WOODSAXON_H

/**
 * @brief Class to implement the Wood-Saxon geometry
 */
class WoodSaxon
{
public:
  double gamma; ///< boost parameter 
  double velocity; ///< the velocity
  double RA0; ///< radius, where density dropped to treshold value
  double RA; ///< radius parameter
  double dA; ///< skin parameter
  double n0A;
  
  
  WoodSaxon() : gamma( 0 ), velocity( 0 ), RA0( 0 ), RA( 0 ), dA( 0 ), n0A( 0 ) {};
  ~WoodSaxon() {};
  
  /**
   * @brief Calculate the parameters.
   *
   * This routine calculates the parameters concerning radii
   * etc. according the input parameters
   * @param[in] A number of nucleons
   * @param[in] b impact parameter (in fm)
   * @param[in] sqrtS The collision energy (in GeV)
   * @return FALSE, if impact parameter too large
   **/
  bool Calculate(double A, double b, double sqrtS);
    
  /**
   * @brief nuclear density (Woods-Saxon distribution) 
   *
   * Calculate the density at transversal distance b and
   * longitudinal distance z. 
   *
   * @param[in] b transversal distance (in fm)
   * @param[in] z longitudinal coordinate (in fm)
   * @return the function value
   **/
  double densityA(double b, double z) const;
    
};

#endif // WOODSAXON_H
