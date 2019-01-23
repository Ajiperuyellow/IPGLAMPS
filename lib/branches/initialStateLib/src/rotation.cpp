//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/rotation.cpp $
//$LastChangedDate: 2015-01-28 18:23:23 +0100 (Mi, 28. Jan 2015) $
//$LastChangedRevision: 2061 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



#include <iostream>
#include <math.h>

#include "rotation.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


///// ROTATE ONE VECTOR /////

// at the moment, we only provide scalar versions...

vector4D rotation::rotate( const vector4D & x ) const
{
  
  return vector4D( Arr[0+4*0]*x[0]-Arr[1+4*0]*x[1]-Arr[2+4*0]*x[2]-Arr[3+4*0]*x[3],
                  -Arr[0+4*3]*x[0]+Arr[1+4*3]*x[1]+Arr[2+4*3]*x[2]+Arr[3+4*3]*x[3],
                  -Arr[0+4*2]*x[0]+Arr[1+4*2]*x[1]+Arr[2+4*2]*x[2]+Arr[3+4*2]*x[3],
                  -Arr[0+4*1]*x[0]+Arr[1+4*1]*x[1]+Arr[2+4*1]*x[2]+Arr[3+4*1]*x[3]);
}

vector4D rotation::rotateInv( const vector4D & x ) const
{
  return vector4D( Arr[0+4*0]*x[0]+Arr[1+4*0]*x[1]+Arr[2+4*0]*x[2]+Arr[3+4*0]*x[3],
                   Arr[0+4*3]*x[0]+Arr[1+4*3]*x[1]+Arr[2+4*3]*x[2]+Arr[3+4*3]*x[3],
                   Arr[0+4*2]*x[0]+Arr[1+4*2]*x[1]+Arr[2+4*2]*x[2]+Arr[3+4*2]*x[3],
                   Arr[0+4*1]*x[0]+Arr[1+4*1]*x[1]+Arr[2+4*1]*x[2]+Arr[3+4*1]*x[3]);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
