//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/FPT_compare.h $
//$LastChangedDate: 2015-03-03 17:31:53 +0100 (Di, 03. Mär 2015) $
//$LastChangedRevision: 2097 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** @file
 * @brief Comparison functions for floating point types (FPT)
 **/

#ifndef FPT_COMPARE_H
#define FPT_COMPARE_H

#include <math.h>

/**
 * @def FPT_COMP_PRECISION
 *
 * The accuracy of the FPT_... comparisons (absolute and relative)
 **/
#define FPT_COMP_PRECISION 1.0e-9

/** 
 * @brief Compute the maximum of three doubles
 * @return max(a,b,c)
 **/
inline double FPT_MAX( const double& a, const double& b, const double& c )
{
  double res = (a > b ? a : b );
  return ( res > c ? res : c );
}

/** 
 * @brief Check two double numbers for "equality"
 * @return a == b
 *
 * @note 
 * please check, whether the actual implementation is okay or whether
 * we should only use the maximum of *absolute* values.!!! 
 * The latter would be more systematic:
 * * a,b>0: the maximal deviation is **max(1,a,b)** times ::FPT_COMP_PRECISION
 * * a,b<0: the maximal deviation is **1** times ::FPT_COMP_PRECISION. 
 * This is different from above! 
 **/
inline bool FPT_COMP_E(const double a, const double b)
{
  return ( fabs( a - b ) <= FPT_COMP_PRECISION * FPT_MAX( 1.0, a, b ) );
  //  return ( fabs( a - b ) <= FPT_COMP_PRECISION * FPT_MAX( 1.0, fabs(a), fabs(b) ) );
}  

/** 
 * @brief Check whether one double is truly smaller than another double
 * @return a < b
 **/
inline bool FPT_COMP_L(const double a, const double b)
{
  return( (a < b) && !FPT_COMP_E(a,b) );
}

/** 
 * @brief Check whether one double is truly larger than another double
 * @return a > b
 **/
inline bool FPT_COMP_G(const double a, const double b)
{
  return( (a > b) && !FPT_COMP_E(a,b) );
}

/** 
 * @brief Check whether one double is smaller than or "equal" to another double
 * @return a <= b
 **/
inline bool FPT_COMP_LE(const double a, const double b)
{
  return( (a < b) || FPT_COMP_E(a,b) );
}

/** @brief Check whether one double is larger than or "equal" to another double
 * @return a >= b
 **/
inline bool FPT_COMP_GE(const double a, const double b) // check for a >= b
{
  return( (a > b) || FPT_COMP_E(a,b) );
}

/** @brief Check double numbers against "zero"
 * @return a == 0.0
 **/
inline bool FPT_COMP_Z(const double a)
{
  return ( fabs( a ) <= FPT_COMP_PRECISION );
}  

/** @brief Check double numbers against "zero"
 * @return a != 0.0
 **/
inline bool FPT_COMP_NZ(const double a)
{
  return ( fabs( a ) > FPT_COMP_PRECISION );
}  

#endif  //#ifndef FPT_COMPARE
