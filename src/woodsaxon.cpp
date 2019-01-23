//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/woodsaxon.cpp $
//$LastChangedDate: 2014-02-12 19:12:33 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1619 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include <string>
#include <iostream>
#include "woodsaxon.h"

using std::cout;
using std::endl;


double WoodSaxon::densityA( double b, double z) const
{
  return gamma*n0A / ( exp(( sqrt( b*b + pow( double( gamma*z ), 2.0 ) ) - RA ) / dA ) + 1.0 );
}

bool WoodSaxon::Calculate(const double A, const double b, const double sqrtS)
{
  const double mproton = 0.938272;
  const double neps = 1.0e-3;
  bool flagOK = true;

  gamma = sqrtS / ( 2 * mproton );
  velocity = sqrt( 1.0 - 1 / pow( gamma, 2 ) );

  RA = 1.12 * pow( A, 1.0 / 3.0 ) - 0.86 * pow( A, -1.0 / 3.0 );
  dA = 0.54;
  n0A = 3 * A / ( 4 * M_PI * pow( RA, 3.0 ) * ( 1.0 + pow(( M_PI * dA / RA ), 2 ) ) );

  RA0 = dA * log( gamma * n0A / neps - 1.0 ) + RA;  //radius where density has dropped to neps = 10^-3

  flagOK = ( b < ( 2 * RA0 ) );

  cout << "WoodSaxon parameters:" << endl;
  cout << "  long. extension = " << 2 * RA0 / gamma << " fm" << endl;
  cout << "  trans. radius = " << RA << " fm" << endl;
  cout << "  impact parameter = " << b << " fm" << endl;
  cout << "  overlap time = " << sqrt( pow( 2 * RA0, 2 ) - pow( b, 2 ) ) / gamma / velocity << " fm/c" << endl;

  return flagOK;
}
