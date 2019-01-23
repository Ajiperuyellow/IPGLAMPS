//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/lorentz.h $
//$LastChangedDate: 2019-01-05 18:04:29 +0100 (Sa, 05. Jan 2019) $
//$LastChangedRevision: 2918 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declaration for routines related to Lorentz transformations
 */

#ifndef LORENTZ_H
#define LORENTZ_H

#include "bampsvector.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>

#ifndef __INTEL_COMPILER
#define _MM_PREFETCH(A,B) _mm_prefetch(A,B)
#else
#define _MM_PREFETCH(A,B) {}
#endif
//#define _MM_PREFETCH(A,B) {}


/** 
 * @brief class to implement Lorentz boosts for vector4D
 * 
 * The internal representation is independent of the user interface.
 *
 * Some routines exists as a 'Scalar' or a vectorized 'SSE3' version. 
 * The scalar versions are implemented as a fallback option and exist
 * as an independent iplementation mainly for testing purposes.
 * The overall behaviour is steared by the compiler define
 * "VECTOR_IMPL_SSE"
 * 
 */
class lorentz
{
public: 
  /** 
   * @brief The type of the numerical representation
   **/
  typedef double Scalar;
  
protected:
  /** 
   * @brief The internal array to store the values
   **/
  SSE_ALIGNED(Scalar) Arr[16];
  
public:
  
  /** 
   * @brief Constructor
   * 
   * Default constructor of an empty boost
   */
  lorentz( ) 
  { 
    SSE_ALIGNED(static const Scalar) Arr1[16] = 
      { 1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0 };
    for (unsigned int i=0;i<16;i++) Arr[i] = Arr1[i];
  };

  /** 
   * @brief Constructor
   * 
   * Default constructor with given beta vector
   */
  lorentz( const vector4D & beta )
  {
    setBeta( beta );
  }

  /** 
   * @brief Set the beta vector
   * 
   * This routine sets the beta vector and calculates the 
   * transformation matrix.
   */
  void setBeta( const vector4D & beta );

  /** 
   * @brief Set the beta vector for transformation to CM
   *
   * This routine calculates the beta vector for a boost to the CM
   * of the given momentum. It sets the beta vector and the
   * transformation matrix.
   **/
  void setBetaCM( const vector4D & x1 );

  /** 
   * @brief Set the beta vector for transformation to CM
   *
   * This routine calculates the beta vector for a boost to the CM
   * of the given momenta. It sets the beta vector and the
   * transformation matrix.
   **/
  void setBetaCM( const vector4D & x1, const vector4D & x2 );

  /** 
   * @brief Set the beta vector for transformation to CM
   *
   * This routine calculates the beta vector for a boost to the CM
   * of the given momenta. It sets the beta vector and the
   * transformation matrix.
   **/
  void setBetaCM( const vector4D & x1, const vector4D & x2, const vector4D & x3 );
  
   /** 
   * @brief get the spatial distance between two vectors
   *
   * 
   **/ 
  double getSpatialDistance(const vector4D & x1, const vector4D & x2);
  
   /** 
   * @brief get the distance of closest approach DCA
   *
   * 
   **/ 
  double getDCAsquared(const vector4D & x1, const vector4D & x2,const vector4D & p1, const vector4D & p2);
  
   
  
  
#if defined (VECTOR_IMPL_SSE) || defined (VECTOR_IMPL_AVX)
protected:
  /** 
   * @brief Some helper routine for set the beta vector
   * 
   */
  inline void setBeta(__m128d & Lo, __m128d & Hi);
#endif

public:
  /** 
   * @brief Return the spatial length of the beta vector
   *
   * This is a shortcut to sqrt(beta().vec2()).
   * It is questionable, how often it is used.
   **/
  Scalar betaVal(void)
  { 
    return sqrt( (Arr[1]*Arr[1]+Arr[2]*Arr[2]+Arr[3]*Arr[3]) )/Arr[0]; 
  };

  /** 
   * @brief Return the boost gamma
   **/
  Scalar gammaVal(void)
  { 
    return 1./(1.-pow ( sqrt( (Arr[1]*Arr[1]+Arr[2]*Arr[2]+Arr[3]*Arr[3]) )/Arr[0],2.0 ) ); 
  };  
  
  
  
  /** 
   * @brief Return the beta vector
   **/
  vector4D beta(void) const
  {
    return vector4D(1.0, Arr[3]/Arr[0], Arr[2]/Arr[0], Arr[1]/Arr[0]); 
  }

  ///// BOOST ONE VECTOR /////

  /** 
   * @brief Return a boosted vector
   **/
  vector4D boost(const vector4D & x) const;

  /** 
   * @brief Return a boosted-back vector
   **/
  vector4D boostInv(const vector4D & x) const;

  /** 
   * @brief Operator overload: Return a boosted vector
   *
   * The multiplication of a boost with a vector returns a boosted vector.
   **/
  vector4D operator* (const vector4D & x) const 
  { return boost(x); }

  /** 
   * @brief Boost one vector
   **/
  void boost(const vector4D & x, vector4D & xNew) const 
  { xNew = boost(x); }; 

  /** 
   * @brief Boost-back one vector
   **/
  void boostInv(const vector4D & x, vector4D & xNew) const 
  { xNew = boostInv(x); }; 


  ///// BOOST TWO VECTORS /////

  /** 
   * @brief Boost two vectors simultanously
   **/
  void boost(const vector4D & x1, const vector4D & x2, 
             vector4D & x1New, vector4D & x2New) const;

  /** 
   * @brief Boost-back two vectors simultanously
   **/
  void boostInv(const vector4D & x1, const vector4D & x2, 
                vector4D & x1New, vector4D & x2New) const;


  ///// BOOST THREE VECTORS /////

  /** 
   * @brief Boost three vectors simultanously
   **/
  void boost(const vector4D & x1, const vector4D & x2, const vector4D & x3, 
             vector4D & x1New, vector4D & x2New, vector4D & x3New) const;

  /** 
   * @brief Boost-back three vectors simultanously
   **/
  void boostInv(const vector4D & x1, const vector4D & x2, const vector4D & x3, 
                vector4D & x1New, vector4D & x2New, vector4D & x3New) const;

  ///// MISC /////

  /** 
   * @brief Output routine
   *
   * Writes out the transformation matrix in human readable form.
   **/
  friend std::ostream& operator<<(std::ostream &os, const lorentz &obj)
  {
    os << std::setprecision(15) << " (" << obj.Arr[0] 
       << ","  << obj.Arr[3] 
       << ","  << obj.Arr[2] 
       << ","  << obj.Arr[1] << ") " << std::endl;
    for (unsigned int i=3;i>0;i--)
      os << std::setprecision(15) << " (" << obj.Arr[0+4*i] 
         << ","  << obj.Arr[3+4*i] 
         << ","  << obj.Arr[2+4*i] 
         << ","  << obj.Arr[1+4*i] << ") " << std::endl;

    return os;
  }

  ////// now the Scalar versions:
  //
  // For testing purposes we may declare the following routines also
  // as public. In the 'final' version, they should be at least
  // 'protected'. 

  //protected:


  /** 
   * @brief Set the beta vector (Scalar version)
   */
  void setBetaScalar( const vector4D & beta );

  /** 
   * @brief Set the beta vector for transformation to CM (Scalar version)
   **/
  void setBetaCMScalar( const vector4D & x1 );

  /** 
   * @brief Set the beta vector for transformation to CM (Scalar version)
   **/
  void setBetaCMScalar( const vector4D & x1, const vector4D & x2 );

  ///// BOOST ONE VECTOR /////

  /** 
   * @brief Return a boosted vector (Scalar version)
   **/
  vector4D boostScalar(const vector4D & x) const;

  /** 
   * @brief Boost one vector (Scalar version)
   **/
  void boostScalar(const vector4D & x, vector4D & xNew) const 
  { xNew = boostScalar(x); }; 

  /** 
   * @brief Return a boosted-back vector (Scalar version)
   **/
  vector4D boostInvScalar(const vector4D & x) const;

  /** 
   * @brief Boost-back one vector (Scalar version)
   **/
  void boostInvScalar(const vector4D & x, vector4D & xNew) const 
  { xNew = boostInvScalar(x); }; 


  ///// BOOST TWO VECTORS /////

  /** 
   * @brief Boost two vectors simultanously (Scalar version)
   **/
  void boostScalar(const vector4D & x1, const vector4D & x2, 
                   vector4D & x1New, vector4D & x2New) const;

  /** 
   * @brief Boost-back two vectors simultanously (Scalar version)
   **/
  void boostInvScalar(const vector4D & x1, const vector4D & x2, 
                      vector4D & x1New, vector4D & x2New) const;

  ///// BOOST THREE VECTORS /////

  void boostScalar( const vector4D & x1, const vector4D & x2, const vector4D & x3,
                    vector4D & x1New, vector4D & x2New, vector4D & x3New ) const;

  void boostInvScalar( const vector4D & x1, const vector4D & x2, const vector4D & x3, 
                       vector4D & x1New, vector4D & x2New, vector4D & x3New ) const;


  FREE_STORE_OPERATORS_ALIGNED
} __attribute__((aligned(VectorAlignment)));

#if defined (VECTOR_IMPL_SSE) || defined (VECTOR_IMPL_AVX)
inline void lorentz::setBeta( __m128d & xG_lo, __m128d & xG_hi)
{
  __m128d xOne = _mm_set1_pd( 1.0 );

  __m128d xBeta_hi = xG_hi; // save it for future usage
  __m128d xBeta_lo = xG_lo;

  xG_lo = _mm_mul_pd( xG_lo, xG_lo ); // square (1.0,Z)
  xG_hi = _mm_mul_pd( xG_hi, xG_hi ); // square (X,Y)

  __m128d xG = _mm_addsub_pd( xG_lo, xG_hi ); // x0=lo0-hi0, x1=lo1+hi1
  xG = _mm_hsub_pd( xG, xG );

  xG = _mm_sqrt_pd( xG );
  xG = _mm_div_pd( xOne, xG );

  // xGG = xG*xG/(xG+1)   == beta*gamma
  __m128d xGplusOne = _mm_add_pd( xOne, xG );
  __m128d xGG = _mm_mul_pd( xG, xG );
  xGG = _mm_div_pd( xGG, xGplusOne );

  // preparations:
  __m128d xZero = _mm_setzero_pd( ); // set it to zero

  _MM_PREFETCH(Arr, _MM_HINT_NTA);

  // now set the array to beta[i]*beta[j] with beta[0]=1

  __m128d xBB0 = _mm_unpacklo_pd( xBeta_lo, xBeta_lo );
  __m128d xBB1 = _mm_unpackhi_pd( xBeta_lo, xBeta_lo );
  __m128d xBB2 = _mm_unpacklo_pd( xBeta_hi, xBeta_hi );
  __m128d xBB3 = _mm_unpackhi_pd( xBeta_hi, xBeta_hi );
  
  __m128d ArrLo0 = _mm_mul_pd( xBB0, xBeta_lo );
  __m128d ArrHi0 = _mm_mul_pd( xBB0, xBeta_hi );

  __m128d ArrLo1 = _mm_mul_pd( xBB1, xBeta_lo );
  __m128d ArrHi1 = _mm_mul_pd( xBB1, xBeta_hi );

  __m128d ArrLo2 = _mm_mul_pd( xBB2, xBeta_lo );
  __m128d ArrHi2 = _mm_mul_pd( xBB2, xBeta_hi );

  __m128d ArrLo3 = _mm_mul_pd( xBB3, xBeta_lo );
  __m128d ArrHi3 = _mm_mul_pd( xBB3, xBeta_hi );

  // now multiply [0][i] with g and [i][j] with gg

  __m128d xGGG = _mm_unpacklo_pd( xG, xGG );

  ArrLo0 = _mm_mul_pd( ArrLo0, xG );
  ArrHi0 = _mm_mul_pd( ArrHi0, xG );

  ArrLo1 = _mm_mul_pd( ArrLo1, xGGG );
  ArrHi1 = _mm_mul_pd( ArrHi1, xGG );

  ArrLo2 = _mm_mul_pd( ArrLo2, xGGG );
  ArrHi2 = _mm_mul_pd( ArrHi2, xGG );

  ArrLo3 = _mm_mul_pd( ArrLo3, xGGG );
  ArrHi3 = _mm_mul_pd( ArrHi3, xGG );

  // now add 1 to the spatial diagonal eements

  __m128d xAdd01 = _mm_unpacklo_pd( xZero, xOne );
  __m128d xAdd10 = _mm_unpacklo_pd( xOne, xZero );
  
  ArrLo1 = _mm_add_pd( ArrLo1, xAdd01 );
  ArrHi2 = _mm_add_pd( ArrHi2, xAdd10 );
  ArrHi3 = _mm_add_pd( ArrHi3, xAdd01 );


  // store the matrix to memory

  _mm_store_pd( Arr+0, ArrLo0 );
  _mm_store_pd( Arr+2, ArrHi0 );
  _mm_store_pd( Arr+4, ArrLo1);
  _mm_store_pd( Arr+6, ArrHi1 );
  _mm_store_pd( Arr+8, ArrLo2 );
  _mm_store_pd( Arr+10, ArrHi2 );
  _mm_store_pd( Arr+12, ArrLo3 );
  _mm_store_pd( Arr+14, ArrHi3 );
} 
#endif


/** 
 * @brief exception class for handling unexpected critical behaviour
 * within lorentz class
 */
class eLorentz_error : public std::runtime_error
{
public:
  explicit eLorentz_error(const std::string& what) : std::runtime_error(what) 
  {};
    
  virtual ~eLorentz_error() throw() {};
};



#endif
