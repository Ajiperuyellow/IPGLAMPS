//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/rangen_distr.cpp $
//$LastChangedDate: 2014-11-16 22:00:01 +0100 (So, 16. Nov 2014) $
//$LastChangedRevision: 1938 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>

#include "rangen_distr.h"

using std::cout;
using std::endl;


ranGen_Distr::ranGen_Distr( const double x[], const double y[], const int N, INTERP_TYPE interp_type) :
  distr(y,x,N,interp_type) // we have to interchange the role of x and y 
{
  CheckMinMax();
}

ranGen_Distr::ranGen_Distr(const std::vector<double>& x, const std::vector<double>& y, INTERP_TYPE interp_type) :
  distr(y,x,interp_type) // we have to interchange the role of x and y 
{
  CheckMinMax();
}

/**
 * We use here an exact comparion of the bounds at 0 and 1, since they
 * may be defined exactly.
 **/
void ranGen_Distr::CheckMinMax(void)
{
  double x1 = distr.minX();
  double x2 = distr.maxX();

  if (x1 < x2)
  {
    if (x1 != 0.0) // exact comparison!
    {
      cout << "CDF(1): min y = " << x1 << endl 
           << "        max y = " << x2 
           << "  " << 1.0-x2 << endl;
      std::string errMsg = "x1 should be 0";
      throw ranGen_error( errMsg );
    }
    if (x2 != 1.0) // exact comparison!
    {
      cout << "CDF(2): min y = " << x1 << endl 
           << "        max y = " << x2
           << "  " << 1.0-x2 << endl;
      std::string errMsg = "x2 should be 1";
      throw ranGen_error( errMsg );
    }
  }
  else
  {
    std::string errMsg = "Error: x2 < x1. monotonically increasing!";
    throw ranGen_error( errMsg );
  }

}
