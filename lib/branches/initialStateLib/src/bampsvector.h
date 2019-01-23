//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/bampsvector.h $
//$LastChangedDate: 2017-04-06 14:55:30 +0200 (Do, 06. Apr 2017) $
//$LastChangedRevision: 2556 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Declaration of 3D and 4D vectors
 *
 * Here we implement the 4D vector classes vector4D, where the basic
 * parameters are (E, pX, pY, pZ).
 *
 * Special improvements are SSE or AVX vectorization implementations
 * of the class methods.
 *
 * There is no distinction between 3D and 4D vectors. There is some
 * subset of routines, which calculates 3D scalar and vector products.
 *
 * The output routines are enhanced by modifiers:
 *  - #plainvector
 *  - #bracketvector
 *
 * Example, how to use this:
 * ~~~
 * std::cout << plainvector << x << std::endl;
 * ~~~
 */

#ifndef BAMPSVECTOR_H
#define BAMPSVECTOR_H

#include "allocator.h"
#include "configBAMPS.h"

#include <iostream>
#include <math.h>
#include <stdint.h> // for uintptr_t

#include <cmath> // std::isinf, std::isnan

//--------------------------------
#if defined(VECTOR_IMPL_SSE)
//--------------------------------

#ifndef HAVE_VECTOR_SSE3
#error "SSE selected, but not provided by CPU"
#endif

#if defined(__GNUC__)
#define SSE_ALIGNED(type) type __attribute__((aligned(VectorAlignment)))
#else
#define SSE_ALIGNED(type) __declspec(align(VectorAlignment)) type
#endif

#include <immintrin.h>

//--------------------------------
#elif defined(VECTOR_IMPL_AVX)
//--------------------------------

#ifndef HAVE_VECTOR_AVX
#error "AVX selected, but not provided by CPU"
#endif

#define SSE_ALIGNED(type) type

#include <immintrin.h>


//--------------------------------
#else // Scalar
//--------------------------------

#define SSE_ALIGNED(type) type

//--------------------------------
#endif
//--------------------------------

#define is_alignedSSE(POINTER)                          \
  (((uintptr_t)(const void *)(POINTER)) % (VectorAlignment) == 0)




/**
 * @brief Access routine used internally for the vector4D output modifiers
 *
 * This is used to get access to std::ios_base::xalloc()
 **/
inline int vector4D_xalloc() { static int i = std::ios_base::xalloc(); return i; }

/** @brief Output modifier: vector without any 'decoration', just the plain numbers **/
inline std::ostream& plainvector(std::ostream& os) { os.iword(vector4D_xalloc()) = 1; return os; }

/** @brief Output modifier: vector with 'decorations' (default) **/
inline std::ostream& bracketvector(std::ostream& os) { os.iword(vector4D_xalloc()) = 0; return os; }




/**
 * @brief class to implement 4D (and implicit 3D) vectors.
 *
 * The internal representation is independent of the user interface.
 * Some routines exists as a 'Scalar' or a vectorized 'SSE3' version.
 * The overall behaviour is steared by the compiler define
 * "VECTOR_IMPL_SSE"
 *
 * The interface seperates between momenta (e.g. Px()) and coordinates
 * (e.g. X()). (At the moment, this distinction is artificial.)
 */
class vector4D
{
public:

  /** @brief The type of the numerical representation */
  typedef double Scalar;

protected:
  /** @brief The internal array to store the values */
  SSE_ALIGNED(Scalar) Mem[4];

public:

  // the class lorentz needs direct access to the stored values:
  friend class lorentz;
  // the class rotation needs direct access to the stored values:
  friend class rotation;

  /**
   * @brief Constructor
   *
   * Default constructor of an empty vector (Px = Py = Pz = E = 0 )
   */
  vector4D( ) : Mem()  { };

  /**
   * @brief Constructor
   *
   * Generic constructor from four scalar values.
   */
  vector4D(const Scalar & E,
           const Scalar & x,
           const Scalar & y,
           const Scalar & z)
  { Mem[0] = E; Mem[1] = z; Mem[2] = y; Mem[3] = x;};

  /**
   * @brief Constructor
   *
   * Generic constructor from only three scalar values
   * (the spatial components). The temporal (0th) component
   * is set to zero.
   */
  vector4D(const Scalar & x,
           const Scalar & y,
           const Scalar & z)
  { Mem[0] = 0.0; Mem[1] = z; Mem[2] = y; Mem[3] = x;};

#if defined(VECTOR_IMPL_SSE) || defined(VECTOR_IMPL_AVX)
  /**
   * @brief Constructor
   *
   * Constructor from two SSE registers
   */
  vector4D(const __m128d & Lo,
           const __m128d & Hi)
  {
    _mm_store_pd(Mem, Lo);
    _mm_store_pd(Mem+2, Hi);
  };

#endif

#if defined(VECTOR_IMPL_AVX)
  /**
   * @brief Constructor
   *
   * Constructor from one AVX register
   */
  vector4D(const __m256d & R)
  {
    _mm256_store_pd(Mem, R);
  };

#endif




  /**
   * @brief Set all entries
   *
   * 'Set routine' from four scalar values.
   */
  void SetTXYZ(const Scalar & t,
               const Scalar & x,
               const Scalar & y,
               const Scalar & z)
  { Mem[0] = t; Mem[1] = z; Mem[2] = y; Mem[3] = x; };

  /**
   * @brief Set all entries
   *
   * 'Set routine' from array of four scalar values. (no checking)
   */
  void SetTXYZ(const Scalar *Arr)
  { Mem[0]=Arr[0]; Mem[1]=Arr[3]; Mem[2]=Arr[2]; Mem[3]=Arr[1]; };

  /** @brief Get x component */
  Scalar Px() const { return Mem[3]; }
  /** @brief Get y component */
  Scalar Py() const { return Mem[2]; }
  /** @brief Get z component */
  Scalar Pz() const { return Mem[1]; }
  /** @brief Get 0th component */
  Scalar E() const { return Mem[0]; }

  /** @brief Set x component */
  Scalar & Px() { return Mem[3]; }
  /** @brief Set y component */
  Scalar & Py() { return Mem[2]; }
  /** @brief Set z component */
  Scalar & Pz() { return Mem[1]; }
  /** @brief Set 0th component */
  Scalar & E() { return Mem[0]; }

  /** @brief Get x component */
  Scalar X() const { return Mem[3]; }
  /** @brief Get y component */
  Scalar Y() const { return Mem[2]; }
  /** @brief Get z component */
  Scalar Z() const { return Mem[1]; }
  /** @brief Get 0th component */
  Scalar T() const { return Mem[0]; }

  /** @brief Set x component */
  Scalar & X() { return Mem[3]; }
  /** @brief Set y component */
  Scalar & Y() { return Mem[2]; }
  /** @brief Set z component */
  Scalar & Z() { return Mem[1]; }
  /** @brief Set 0th component */
  Scalar & T() { return Mem[0]; }


  /** @brief Get x component squared */
  Scalar Px2() const { return Mem[3]*Mem[3]; }
  /** @brief Get y component squared */
  Scalar Py2() const { return Mem[2]*Mem[2]; }
  /** @brief Get z component squared */
  Scalar Pz2() const { return Mem[1]*Mem[1]; }
  /** @brief Get 0th component squared */
  Scalar E2() const { return Mem[0]*Mem[0]; }

  /** @brief Get x component squared */
  Scalar X2() const { return Mem[3]*Mem[3]; }
  /** @brief Get y component squared */
  Scalar Y2() const { return Mem[2]*Mem[2]; }
  /** @brief Get z component squared */
  Scalar Z2() const { return Mem[1]*Mem[1]; }
  /** @brief Get 0th component squared */
  Scalar T2() const { return Mem[0]*Mem[0]; }

  /**
   * @brief Get entry
   *
   * return the ith entry: i=0..3 -> t,x,y,z
   */
  Scalar operator()(unsigned int i) const { return Mem[i==0?i:4-i]; }
  /**
   * @brief Set entry
   *
   * set the ith entry: i=0..3 -> t,x,y,z
   */
  Scalar & operator()(unsigned int i) { return Mem[i==0?i:4-i]; }

protected:
  /**
   * @brief Get entry
   *
   * return the ith entry, implementation dependent!
   */
  Scalar operator[](unsigned int i) const { return Mem[i]; }
  /**
   * @brief Set entry
   *
   * set the ith entry, implementation dependent!
   */
  Scalar & operator[](unsigned int i) { return Mem[i]; }

public:

  /**
   * @brief Get all entries
   *
   * get routine from array of four scalar values. (no checking!)
   */
  void GetTXYZ(Scalar *Arr) const
  { Arr[0]=T(); Arr[1]=X(); Arr[2]=Y(); Arr[3]=Z(); };

  /**
   * @brief Get all entries
   *
   * get routine from four scalar values.
   */
  void GetTXYZ(Scalar & t,
               Scalar & x,
               Scalar & y,
               Scalar & z) const
  { t=T(); x=X(); y=Y(); z=Z(); };

  /**
   * @brief Get all entries
   *
   * get routine from array of four scalar values. (no checking!)
   */
  void GetEPxPyPz(Scalar *Arr) const
  { Arr[0]=T(); Arr[1]=X(); Arr[2]=Y(); Arr[3]=Z(); };

  /**
   * @brief Check all entries against std::isinf
   */
  bool isInf() const
  {
    return std::isinf(T()) || std::isinf(X()) || std::isinf(Y()) || std::isinf(Z());
  }

  /**
   * @brief Check all entries against std::isnan
   */
  bool isNan() const
  {
    return std::isnan(T()) || std::isnan(X()) || std::isnan(Y()) || std::isnan(Z());
  }


  /**
   * @brief Set Energy as sqrt( p^2+m^2 )
   **/
  void SetEbyM(double mass)
  {
    E() = sqrt( vec2() + mass*mass );
  };

  /**
   * @brief Negate the 3D component
   **/
  void Minus3()
  {
    X() = -X(); Y() = -Y(); Z() = -Z();
  };

  /** @brief Unary minus. */
  inline vector4D operator - () const;

  /** @brief The Plus Operator */
  inline vector4D operator + (const vector4D &) const;

  /** @brief The Plus Operator with a scalar */
  inline vector4D operator + (const Scalar a) const;


  /** @brief The Increment Operator */
  inline vector4D & operator += (const vector4D &);

  /** @brief The Minus Operator */
  inline vector4D operator - (const vector4D &) const;

  /** @brief The Minus Operator with a scalar*/
  inline vector4D operator - (const Scalar a) const;

  /** @brief The Decrement Operator */
  inline vector4D & operator -= (const vector4D &);

  /**
   * @brief The Times Operator
   *
   * This is the multiplication with a scalar.
   * We intentionaly do not implement the division operator, because
   * it is better to use "*(1/a)" instead.
   **/
  inline vector4D operator * (const Scalar a) const;

  /**
   * @brief The Multiplicated Operator
   *
   * This is the multiplication with a scalar.
   * We intentionaly do not implement the division operator, because
   * it is better to use "*(1/a)" instead.
   **/
  inline vector4D & operator *= (const Scalar a);


  /**
   * @brief The Minimum
   *
   * returns a vector which holds the minimum of each component
   **/
  friend inline vector4D min(const vector4D &a, const vector4D &b );

  /**
   * @brief The Maximum
   *
   * returns a vector which holds the maximum of each component
   **/
  friend inline vector4D max(const vector4D &a, const vector4D &b );

  /**
   * @brief Normalize to the zeroth component
   *
   * i.e.: x[] *= 1/x[0]
   **/
  inline vector4D NormalizeToE(void);


  /** @brief The squared length of the 3D part */
  inline Scalar vec2() const;

  /** @brief The squared length of the vector (+--- metric) */
  inline Scalar M2() const;

  /** @brief The squared length of the vector (+--- metric) */
  inline Scalar M2_Scalar() const
  {
    return( Mem[0]*Mem[0] - Mem[1]*Mem[1] - Mem[2]*Mem[2] - Mem[3]*Mem[3] );
  }


  /** @brief The transverse momentum squared */
  Scalar Pt2() const { return X()*X()+Y()*Y(); };

  /** @brief The transverse momentum */
  Scalar Pt() const { return sqrt(Pt2()); };

  /** @brief The transverse component squared */
  Scalar Perp2() const { return Pt2(); };

  /** @brief The transverse component */
  Scalar Perp() const { return sqrt(Pt2()); };

  /** @brief The squared transverse mass
   *
   * We use here MT^2 = PT^2+ M^2 with M given as parameter.
   * Instead one can also use E^2-PZ^2, but there are some minor
   * numerical differences.
   */
  Scalar Mt2(const Scalar m) const { return m*m + Pt2(); };

  /** @brief The transverse mass
   *
   * We use here MT^2 = PT^2+ M^2 with M given as parameter.
   * Instead one can also use E^2-PZ^2, but there are some minor
   * numerical differences.
   */
  Scalar Mt(const Scalar m) const { return sqrt( Mt2(m) ); };

  /**
   * @brief The (longitudinal) eigentime tau squared
   **/
  Scalar Eigentime2() const { return T2()-Z2(); };

  /**
   * @brief The (longitudinal) eigentime tau
   **/
  Scalar Eigentime() const { return sqrt(Eigentime2()); };


  /**
   * @brief The (longitudinal) rapidity
   *
   * returns 0.5*log( (E+Pz)/(E-Pz) )
   **/
  Scalar Rapidity() const { return 0.5*log( (E()+Pz())/(E()-Pz()) ); };
  
  /**
   * @brief The (longitudinal) rapidity (with shift)
   *
   * returns 0.5*log( (E+Pz-shift)/(E-Pz-shift) )
   **/
  Scalar Rapidity(const Scalar shift) const { return 0.5*log( (E()+Pz()-shift)/(E()-Pz()-shift) ); };


  /**
   * @brief The (longitudinal) pseudorapidity
   *
   * @param[in] m The mass of the particle
   *
   * returns 0.5*log( (P+Pz)/(P-Pz) )
   *
   * One has to give the mass explicit as argument, since the length
   * of the vector is calculated vis sqrt(E^2-m^2), and not via #vec2.
   * This should be checked and improved!!!
   **/
  Scalar Pseudorapidity(const Scalar m) const;

  /** @brief The 4D scalar product of two vectors */
  friend Scalar Dot(const vector4D &x,const vector4D &y )
  {
    return ( (x[0]*y[0] - x[1]*y[1] - x[2]*y[2] - x[3]*y[3]) );
  }

  /** @brief The 3D scalar product of two vectors */
  friend Scalar Dot3(const vector4D &x,const vector4D &y )
  {
    return ( (x[1]*y[1] + x[2]*y[2] + x[3]*y[3]) );
  }

  /** @brief The angle between two vectors */
  friend Scalar CosTheta(const vector4D &x,const vector4D &y )
  {
    return ( (x[1]*y[1] + x[2]*y[2] + x[3]*y[3])/sqrt(x.vec2()*y.vec2()) );
  }

  /** @brief The angle between two vectors (assumed to be massless) */
  friend Scalar CosThetaMassless(const vector4D &x,const vector4D &y )
  {
    return ( (x[1]*y[1] + x[2]*y[2] + x[3]*y[3])/(x[0]*y[0]) );
  }

  /** @brief The angle between two vectors in transversal plane*/
  friend Scalar CosPhi(const vector4D &x,const vector4D &y )
  {
    return ( (x[2]*y[2] + x[3]*y[3])/sqrt(x.Perp2()*y.Perp2()) );
  }

  /** @brief The angle between the vector and the x-axis in transversal plane
   *
   * @return Principal arc tangent of y/x, in the interval [-pi,+pi] radians.
   */
  inline Scalar Phi_x( ) const
  {
    return ( atan2( Y(), X() ) );
  }

  /** @brief The angle between the vector and the y-axis in transversal plane
   *
   * @return Principal arc tangent of x/y, in the interval [-pi,+pi] radians.
   */

  inline Scalar Phi_y( ) const
  {
    return ( atan2( X(), Y() ) );
  }


  /** @brief The 3D cross product of two vectors */
  friend vector4D Cross(const vector4D &x,const vector4D &y )
  {
    return vector4D( 0.0,
                     x.Y()*y.Z()-x.Z()*y.Y(),
                     x.Z()*y.X()-x.X()*y.Z(),
                     x.X()*y.Y()-x.Y()*y.X() );
  }

  /** @brief Square every component */
  friend vector4D Square(const vector4D &x)
  {
    return vector4D( x[0]*x[0], x[1]*x[1], x[2]*x[2], x[3]*x[3] );
  }

  /** @brief The relative velocity of two particles
   *
   * @param[in] p1 4-momentum of particle 1
   * @param[in] p2 4-momentum of particle 2
   * @param[in] M1 mass of particle 1
   * @param[in] M2 mass of particle 2
   *
   * no checking is done, whether sqrt or division by energies is
   * possible.
   *
   * In principle, M1^2 and M2^2 can be calculated by P1^2 and P2^2,
   * but if we know them, giving them as parameter is faster.
   */
  friend Scalar VelRel(const vector4D &p1, const vector4D &p2, const Scalar M1, const Scalar M2)
  {
    return ( sqrt( pow( Dot(p1,p2), 2.0) - pow( M1*M2, 2.0))/( p1.E()*p2.E() ) );
  }

  /** @brief The v2 flow component in x-y-frame */
  Scalar FlowV2(void) const
  {
    return ( (Px2() - Py2())/Perp2() );
  }

  /** @brief The v4 flow component in x-y-frame */
  Scalar FlowV4(void) const
  {
    return ( pow(Px2(),2) - 6.0 * Px2() * Py2() + pow(Py2(),2) ) / pow(Perp2(),2);
  }

  // ------ Input/Output ------

  /**
   * @brief The standard output routine
   *
   * The output may be changed by the modifiers:
   *   #plainvector
   *   #bracketvector
   */
  friend std::ostream& operator<<(std::ostream &os, const vector4D &obj)
  {
    int osval = os.iword(vector4D_xalloc());
    switch (osval)
    {
      case 0:
      {
        os << "<" << obj.T() << "," << obj.X() << "," << obj.Y() << "," << obj.Z() << ">";
      } break;
      case 1:
      {
        os << "" << obj.T() << " " << obj.X() << " " << obj.Y() << " " << obj.Z() << "";
      } break;
      default:
      {
        os << "***** wrong value of output style = " << osval << "*****";
      }
    }
    return os;
  }

  FREE_STORE_OPERATORS_ALIGNED
} __attribute__((aligned(VectorAlignment)));



/////
///// Some inline definitions:
/////

inline vector4D vector4D::operator - () const
{
  return vector4D( -T(), -X(), -Y(), -Z() );
}



inline vector4D vector4D::operator + (const vector4D & q) const
{
#if defined(VECTOR_IMPL_SSE)

  __m128d qlo = _mm_load_pd( q.Mem );
  __m128d qhi = _mm_load_pd( q.Mem+2 );
  __m128d xlo = _mm_load_pd( Mem );
  __m128d xhi = _mm_load_pd( Mem+2 );
  xlo = _mm_add_pd( xlo, qlo );
  xhi = _mm_add_pd( xhi, qhi );

  return vector4D( xlo, xhi );

#elif defined(VECTOR_IMPL_AVX)

  __m256d qq = _mm256_load_pd( q.Mem );
  __m256d xx = _mm256_load_pd( Mem );
  xx = _mm256_add_pd( xx, qq );

  return vector4D( xx );

#else
  return vector4D( T()+q.T(), X()+q.X(), Y()+q.Y(), Z()+q.Z() );
#endif
}

inline vector4D vector4D::operator + (const Scalar a) const
{
#if defined(VECTOR_IMPL_SSE)

  __m128d aa = _mm_load1_pd( &a );
  __m128d xlo = _mm_load_pd( Mem );
  __m128d xhi = _mm_load_pd( Mem+2 );
  xlo = _mm_add_pd( xlo, aa );
  xhi = _mm_add_pd( xhi, aa );

  return vector4D( xlo, xhi );

#elif defined(VECTOR_IMPL_AVX)

  __m256d aa = _mm256_set1_pd( a );
  __m256d xx = _mm256_load_pd( Mem );
  xx = _mm256_add_pd( xx, aa );

  return vector4D( xx );

#else
  return vector4D( T()+a, X()+a, Y()+a, Z()+a );
#endif
}

inline vector4D &vector4D::operator += (const vector4D & q)
{
#if defined(VECTOR_IMPL_SSE)

  __m128d qlo = _mm_load_pd( q.Mem );
  __m128d qhi = _mm_load_pd( q.Mem+2 );
  __m128d xlo = _mm_load_pd( Mem );
  __m128d xhi = _mm_load_pd( Mem+2 );
  xlo = _mm_add_pd( xlo, qlo );
  xhi = _mm_add_pd( xhi, qhi );

  _mm_store_pd( Mem,   xlo );
  _mm_store_pd( Mem+2, xhi );

#elif defined(VECTOR_IMPL_AVX)

  __m256d qq = _mm256_load_pd( q.Mem );
  __m256d xx = _mm256_load_pd( Mem );
  xx = _mm256_add_pd( xx, qq );
  _mm256_store_pd( Mem, xx );

#else

  Mem[0] += q.Mem[0];
  Mem[1] += q.Mem[1];
  Mem[2] += q.Mem[2];
  Mem[3] += q.Mem[3];

#endif

  return *this;
}

inline vector4D vector4D::operator - (const vector4D & q) const
{
#if defined(VECTOR_IMPL_SSE)

  __m128d qlo = _mm_load_pd( q.Mem );
  __m128d qhi = _mm_load_pd( q.Mem+2 );
  __m128d xlo = _mm_load_pd( Mem );
  __m128d xhi = _mm_load_pd( Mem+2 );
  xlo = _mm_sub_pd( xlo, qlo );
  xhi = _mm_sub_pd( xhi, qhi );

  return vector4D( xlo, xhi );

#elif defined(VECTOR_IMPL_AVX)

  __m256d qq = _mm256_load_pd( q.Mem );
  __m256d xx = _mm256_load_pd( Mem );
  xx = _mm256_sub_pd( xx, qq );

  return vector4D( xx );

#else
  return vector4D( T()-q.T(), X()-q.X(), Y()-q.Y(), Z()-q.Z() );
#endif
}

inline vector4D vector4D::operator - (const Scalar a) const
{
#if defined(VECTOR_IMPL_SSE)

  __m128d aa = _mm_load1_pd( &a );
  __m128d xlo = _mm_load_pd( Mem );
  __m128d xhi = _mm_load_pd( Mem+2 );
  xlo = _mm_sub_pd( xlo, aa );
  xhi = _mm_sub_pd( xhi, aa );

  return vector4D( xlo, xhi );

#elif defined(VECTOR_IMPL_AVX)

  __m256d aa = _mm256_set1_pd( a );
  __m256d xx = _mm256_load_pd( Mem );
  xx = _mm256_sub_pd( xx, aa );

  return vector4D( xx );

#else
  return vector4D( T()-a, X()-a, Y()-a, Z()-a );
#endif
}

inline vector4D &vector4D::operator -= (const vector4D & q)
{
#if defined(VECTOR_IMPL_SSE)

  __m128d qlo = _mm_load_pd( q.Mem );
  __m128d qhi = _mm_load_pd( q.Mem+2 );
  __m128d xlo = _mm_load_pd( Mem );
  __m128d xhi = _mm_load_pd( Mem+2 );
  xlo = _mm_sub_pd( xlo, qlo );
  xhi = _mm_sub_pd( xhi, qhi );

  _mm_store_pd( Mem,   xlo );
  _mm_store_pd( Mem+2, xhi );

#elif defined(VECTOR_IMPL_AVX)

  __m256d qq = _mm256_load_pd( q.Mem );
  __m256d xx = _mm256_load_pd( Mem );
  xx = _mm256_sub_pd( xx, qq );

  _mm256_store_pd( Mem, xx );

#else

  Mem[0] -= q.Mem[0];
  Mem[1] -= q.Mem[1];
  Mem[2] -= q.Mem[2];
  Mem[3] -= q.Mem[3];

#endif

  return *this;
}

inline vector4D   vector4D::operator * (const vector4D::Scalar a) const
{
  return vector4D( a*T(), a*X(), a*Y(), a*Z() );
}

inline vector4D & vector4D::operator *= (const vector4D::Scalar a)
{
  Mem[0]*=a;
  Mem[1]*=a;
  Mem[2]*=a;
  Mem[3]*=a;
  return *this;
}

inline vector4D min(const vector4D & a, const vector4D & b)
{
#if defined(VECTOR_IMPL_SSE) || defined(VECTOR_IMPL_AVX)

  __m128d alo = _mm_load_pd( a.Mem );
  __m128d ahi = _mm_load_pd( a.Mem+2 );
  __m128d xlo = _mm_load_pd( b.Mem );
  __m128d xhi = _mm_load_pd( b.Mem+2 );
  xlo = _mm_min_pd( xlo, alo );
  xhi = _mm_min_pd( xhi, ahi );

  return vector4D( xlo, xhi );

#else
  return vector4D( std::min(a.T(),b.T()),
                   std::min(a.X(),b.X()),
                   std::min(a.Y(),b.Y()),
                   std::min(a.Z(),b.Z()) );
#endif
}

inline vector4D max(const vector4D & a, const vector4D & b)
{
#if defined(VECTOR_IMPL_SSE) || defined(VECTOR_IMPL_AVX)

  __m128d alo = _mm_load_pd( a.Mem );
  __m128d ahi = _mm_load_pd( a.Mem+2 );
  __m128d xlo = _mm_load_pd( b.Mem );
  __m128d xhi = _mm_load_pd( b.Mem+2 );
  xlo = _mm_max_pd( xlo, alo );
  xhi = _mm_max_pd( xhi, ahi );

  return vector4D( xlo, xhi );

#else
  return vector4D( std::max(a.T(),b.T()),
                   std::max(a.X(),b.X()),
                   std::max(a.Y(),b.Y()),
                   std::max(a.Z(),b.Z()) );
#endif
}

inline vector4D::Scalar vector4D::Pseudorapidity(const vector4D::Scalar m) const
{
  vector4D::Scalar pp = sqrt( E()*E() - m*m );
  return 0.5*log( (pp+Pz())/(pp-Pz()) );
}

inline vector4D vector4D::NormalizeToE(void)
{
#if defined(VECTOR_IMPL_SSE) || defined(VECTOR_IMPL_AVX)

  __m128d xlo = _mm_load_pd( Mem );
  __m128d xhi = _mm_load_pd( Mem+2 );
  __m128d xOne = _mm_set1_pd( 1.0 );
  __m128d xE = _mm_unpacklo_pd( xlo, xlo );
  xE = _mm_div_pd( xOne, xE );

  xlo = _mm_mul_pd( xlo, xE );
  xhi = _mm_mul_pd( xhi, xE );

  _mm_store_pd( Mem, xlo );
  _mm_store_pd( Mem+2, xhi );

#else
  double a = 1/Mem[0];
  Mem[0]*=a;
  Mem[1]*=a;
  Mem[2]*=a;
  Mem[3]*=a;
#endif
  return *this;
}

inline vector4D::Scalar vector4D::vec2() const
{
#if defined(VECTOR_IMPL_SSE) || defined(VECTOR_IMPL_AVX)

  __m128d x1 = _mm_load_sd(Mem+1); // only load 1 double (other==0)
  __m128d x2 = _mm_load_pd(Mem+2);

  x1 = _mm_mul_pd( x1, x1); // square (0,Z)
  x2 = _mm_mul_pd( x2, x2); // square (X,Y)
  x1 = _mm_hadd_pd( x1, x2); // horizontal add, step 1
  x1 = _mm_hadd_pd( x1, x1); // horizontal add, step 2

  SSE_ALIGNED(double) c;
  _mm_store_sd(&c, x1);
  return c;

#else
  return( Mem[1]*Mem[1] + Mem[2]*Mem[2] + Mem[3]*Mem[3] );
#endif
}

inline vector4D::Scalar vector4D::M2() const
{
#if defined(VECTOR_IMPL_SSE) || defined(VECTOR_IMPL_AVX)

  __m128d x1 = _mm_load_pd(Mem);
  __m128d x2 = _mm_load_pd(Mem+2);

  x1 = _mm_mul_pd( x1, x1); // square (E,Z)
  x2 = _mm_mul_pd( x2, x2); // square (X,Y)

  x1 = _mm_addsub_pd( x1, x2 ); // x0=lo0-hi0, x1=lo1+hi1
  x1 = _mm_hsub_pd( x1, x1 );

  SSE_ALIGNED(double) c;
  _mm_store_sd(&c, x1);
  return c;

#else
  return( Mem[0]*Mem[0] - Mem[1]*Mem[1] - Mem[2]*Mem[2] - Mem[3]*Mem[3] );
#endif
}


/**
 * @brief Define the type of a (T,X,Y,Z) - vector
 */
typedef vector4D VectorTXYZ;

/**
 * @brief Define the type of a (E,Px,Py,Pz) - vector
 */
typedef vector4D VectorEPxPyPz;

/**
 * @brief Define the type of a (X,Y,Z) - vector
 *
 * This is identical to a (T,X,Y,Z) - vector (but T is unused)
 */
typedef vector4D VectorXYZ;

#ifdef BAMPS_DECLARE_ALLOCATOR

BAMPS_DECLARE_ALLOCATOR(vector4D);
// BAMPS_DECLARE_ALLOCATOR(VectorTXYZ);
// BAMPS_DECLARE_ALLOCATOR(VectorEPxPyPz);
// BAMPS_DECLARE_ALLOCATOR(VectorXYZ);

#else
#warning "ALLOCATOR not defined"
#endif

#pragma omp declare                                     \
  reduction(+ : vector4D :                              \
            omp_out += omp_in )                         \
  initializer( omp_priv = omp_orig )

#endif
