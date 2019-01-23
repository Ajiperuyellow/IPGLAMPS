//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/lorentz.cpp $
//$LastChangedDate: 2019-01-05 18:04:29 +0100 (Sa, 05. Jan 2019) $
//$LastChangedRevision: 2918 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



#include <iostream>
#include <math.h>

#include "lorentz.h"


///// set Beta /////

void lorentz::setBeta( const vector4D & beta0 )
{
#if defined (VECTOR_IMPL_SSE) || defined (VECTOR_IMPL_AVX)

  SSE_ALIGNED(static const double) One = 1.0;

  __m128d xG_lo = _mm_load_pd( beta0.Mem );
  __m128d xG_hi = _mm_load_pd( beta0.Mem+2 );

  xG_lo = _mm_loadl_pd( xG_lo, &One ); // set the lower value to 1.0

  setBeta( xG_lo, xG_hi );

#else

  setBetaScalar( beta0 );

#endif

}

void lorentz::setBetaScalar( const vector4D & beta0 )
{

  vector4D beta = beta0;

  if (fabs(beta(0)-1.0)>1e-5)
  {
    // std::cout << "WARNING: beta(0) != 1 !!!" << std::endl;
    // std::cout << beta << std::endl << std::endl;

    // std::string errMsg = "WARNING: beta(0) != 1 !!!";
    // throw eLorentz_error( errMsg );

    beta(0) = 1.0;
  }

  double g = beta.M2_Scalar(); // = 1 - b^2
  if (g < 0)
  {
    std::cout << "ERROR: beta^2 > 1: g=" << g << std::endl;
    std::cout << beta << std::endl << std::endl;
    return;
  }
  g = 1/sqrt(g);
  double gg = g*g/(1+g);

  for (unsigned int i=1;i<4;i++)
    for (unsigned int j=1;j<4;j++)
      Arr[i+4*j] = gg * beta[i] * beta[j];
  for (unsigned int i=1;i<4;i++)
    Arr[i] = Arr[4*i] = g * beta[i];
  Arr[0] = g;
  for (unsigned int i=1;i<4;i++)
    Arr[i+4*i] += 1;
}

void lorentz::setBetaCM( const vector4D & x1 )
{
#if defined (VECTOR_IMPL_SSE) || defined (VECTOR_IMPL_AVX)

  // x1:

  __m128d xlo = _mm_load_pd( x1.Mem );
  __m128d xhi = _mm_load_pd( x1.Mem+2 );

  // normalize:

  __m128d xOne = _mm_set1_pd( 1.0 );
  __m128d xE = _mm_unpacklo_pd( xlo, xlo );
  xE = _mm_div_pd( xOne, xE );

  xlo = _mm_mul_pd( xlo, xE );
  xhi = _mm_mul_pd( xhi, xE );

  // set beta:

  setBeta( xlo, xhi );

#else
  setBetaCMScalar( x1 );
#endif
}

void lorentz::setBetaCMScalar( const vector4D & x1 )
{
  vector4D beta = x1;
  beta *= 1/beta[0];
  setBetaScalar( beta );
}

void lorentz::setBetaCM( const vector4D & x1, const vector4D & x2 )
{
#if defined (VECTOR_IMPL_SSE) || defined (VECTOR_IMPL_AVX)

  // x1 + x2:

  __m128d xlo = _mm_load_pd( x1.Mem );
  __m128d xhi = _mm_load_pd( x1.Mem+2 );
  __m128d xlo2 = _mm_load_pd( x2.Mem );
  __m128d xhi2 = _mm_load_pd( x2.Mem+2 );
  xlo = _mm_add_pd( xlo, xlo2 );
  xhi = _mm_add_pd( xhi, xhi2 );

  // normalize:

  __m128d xOne = _mm_set1_pd( 1.0 );
  __m128d xE = _mm_unpacklo_pd( xlo, xlo );
  xE = _mm_div_pd( xOne, xE );

  xlo = _mm_mul_pd( xlo, xE );
  xhi = _mm_mul_pd( xhi, xE );

  // set beta:

  setBeta( xlo, xhi );

#else
  setBetaCMScalar( x1, x2 );
#endif
}

void lorentz::setBetaCMScalar( const vector4D & x1, const vector4D & x2 )
{
  vector4D beta = x1 + x2;
  beta *= 1/beta[0];
  setBetaScalar( beta );
}

void lorentz::setBetaCM( const vector4D & x1, const vector4D & x2, const vector4D & x3)
{
  setBeta( (x1+x2+x3).NormalizeToE() );
}

///// BOOST ONE VECTOR /////

vector4D lorentz::boost( const vector4D & x ) const
{
#if defined (VECTOR_IMPL_SSE)
  // please note: we have to be careful with the signs of the matrix:
  // The components [0][i] and [i][0] have to be negated!

  // we also use _mm_prefetch in order to load some memory into a
  // cache closer to the processor.
  _MM_PREFETCH(x.Mem, _MM_HINT_NTA);

  SSE_ALIGNED(const unsigned long long) S = 0x8000000000000000ull ;

  // cppcheck-suppress invalidPointerCast
  __m128d xS10 = _mm_load_sd( (double *)&S ); // only lower part negates
  __m128d xS01 = _mm_setzero_pd( ); // set it to zero
  __m128d xS11 = _mm_unpacklo_pd( xS10, xS10 ); // both parts negate
  xS01 = _mm_unpacklo_pd( xS01, xS10 ); // only upper part negates

  _MM_PREFETCH(Arr, _MM_HINT_NTA);

  __m128d V1 = _mm_load_pd(x.Mem);
  __m128d V3 = _mm_load_pd(x.Mem+2);
  __m128d V0 = _mm_unpacklo_pd( V1, V1 );
  __m128d V2 = _mm_unpacklo_pd( V3, V3 );
  V1 = _mm_unpackhi_pd( V1, V1 );
  V3 = _mm_unpackhi_pd( V3, V3 );

  // load the matrix components into SSE registers:
  __m128d ArrLo0 = _mm_load_pd( Arr+0 );
  __m128d ArrHi0 = _mm_load_pd( Arr+2 );
  __m128d ArrLo1 = _mm_load_pd( Arr+4 );
  __m128d ArrHi1 = _mm_load_pd( Arr+6 );
  __m128d ArrLo2 = _mm_load_pd( Arr+8 );
  __m128d ArrHi2 = _mm_load_pd( Arr+10 );
  __m128d ArrLo3 = _mm_load_pd( Arr+12 );
  __m128d ArrHi3 = _mm_load_pd( Arr+14 );

  // correct the signs of the matrix components:
  ArrLo0 = _mm_xor_pd( ArrLo0, xS01 );
  ArrHi0 = _mm_xor_pd( ArrHi0, xS11 );
  ArrLo1 = _mm_xor_pd( ArrLo1, xS10 );
  ArrLo2 = _mm_xor_pd( ArrLo2, xS10 );
  ArrLo3 = _mm_xor_pd( ArrLo3, xS10 );


  __m128d resLo = _mm_mul_pd( V0, ArrLo0 );
  __m128d resHi = _mm_mul_pd( V0, ArrHi0 );

  resLo = _mm_add_pd( resLo, _mm_mul_pd( V1, ArrLo1 ));
  resHi = _mm_add_pd( resHi, _mm_mul_pd( V1, ArrHi1 ));

  resLo = _mm_add_pd( resLo, _mm_mul_pd( V2, ArrLo2 ));
  resHi = _mm_add_pd( resHi, _mm_mul_pd( V2, ArrHi2 ));

  resLo = _mm_add_pd( resLo, _mm_mul_pd( V3, ArrLo3 ));
  resHi = _mm_add_pd( resHi, _mm_mul_pd( V3, ArrHi3 ));

  return vector4D( resLo, resHi );

#elif defined (VECTOR_IMPL_AVX)

  // __m256d a1 = _mm256_set_pd(1.0,2.0,3.0,4.0);
  // __m256d a2 = _mm256_set_pd(10.0,20.0,30.0,40.0);

  // cout << (vector4D)a1  << endl;
  // cout << (vector4D)a2  << endl;
  // __m256d b1 = _mm256_addsub_pd(a1,a2);
  // cout << (vector4D)b1  << endl;


  // load the vector into registers:
  // (V0 consist of 4 copies of component 0, V1 has 4 copies of
  // component 1 and so on)
  __m256d V0 = _mm256_broadcast_sd(x.Mem);
  __m256d V1 = _mm256_broadcast_sd(x.Mem+1);
  __m256d V2 = _mm256_broadcast_sd(x.Mem+2);
  __m256d V3 = _mm256_broadcast_sd(x.Mem+3);

  // load the matrix components into registers:
  __m256d Arr0 = _mm256_load_pd( Arr+0 );
  __m256d Arr1 = _mm256_load_pd( Arr+4 );
  __m256d Arr2 = _mm256_load_pd( Arr+8 );
  __m256d Arr3 = _mm256_load_pd( Arr+12 );

  // correct the signs of the matrix components:
  // in Arr0, we have to flip signs .123
  // in Arri, we have to flip signs 0...

  static const __m256d S01 = _mm256_set_pd(-0.,-0.,-0., 0.);
  static const __m256d S10 = _mm256_set_pd( 0., 0., 0.,-0.);

  Arr0 = _mm256_xor_pd( Arr0, S01 );
  Arr1 = _mm256_xor_pd( Arr1, S10 );
  Arr2 = _mm256_xor_pd( Arr2, S10 );
  Arr3 = _mm256_xor_pd( Arr3, S10 );

  // cout << "changed Arr0:" << endl;
  // cout << (vector4D)Arr0  << endl;
  // cout << (vector4D)Arr1  << endl;
  // cout << (vector4D)Arr2  << endl;
  // cout << (vector4D)Arr3  << endl;

  // do the actual multiplication:
  __m256d res = _mm256_mul_pd( V0, Arr0 );
  res = _mm256_add_pd( res, _mm256_mul_pd( V1, Arr1 ));
  res = _mm256_add_pd( res, _mm256_mul_pd( V2, Arr2 ));
  res = _mm256_add_pd( res, _mm256_mul_pd( V3, Arr3 ));
  
  return vector4D( res );

#else

  return boostScalar( x );

#endif
}

vector4D lorentz::boostScalar( const vector4D & x ) const
{

  return vector4D( Arr[0+4*0]*x[0]-Arr[1+4*0]*x[1]-Arr[2+4*0]*x[2]-Arr[3+4*0]*x[3],
                  -Arr[0+4*3]*x[0]+Arr[1+4*3]*x[1]+Arr[2+4*3]*x[2]+Arr[3+4*3]*x[3],
                  -Arr[0+4*2]*x[0]+Arr[1+4*2]*x[1]+Arr[2+4*2]*x[2]+Arr[3+4*2]*x[3],
                  -Arr[0+4*1]*x[0]+Arr[1+4*1]*x[1]+Arr[2+4*1]*x[2]+Arr[3+4*1]*x[3]);
}


vector4D lorentz::boostInv( const vector4D & x ) const
{
#if defined (VECTOR_IMPL_SSE)

  __m128d V1 = _mm_load_pd(x.Mem);
  __m128d V3 = _mm_load_pd(x.Mem+2);
  __m128d V0 = _mm_unpacklo_pd( V1, V1 );
  __m128d V2 = _mm_unpacklo_pd( V3, V3 );
  V1 = _mm_unpackhi_pd( V1, V1 );
  V3 = _mm_unpackhi_pd( V3, V3 );

  // load the matrix components into SSE registers:
  __m128d ArrLo0 = _mm_load_pd( Arr+0 );
  __m128d ArrHi0 = _mm_load_pd( Arr+2 );
  __m128d ArrLo1 = _mm_load_pd( Arr+4 );
  __m128d ArrHi1 = _mm_load_pd( Arr+6 );
  __m128d ArrLo2 = _mm_load_pd( Arr+8 );
  __m128d ArrHi2 = _mm_load_pd( Arr+10 );
  __m128d ArrLo3 = _mm_load_pd( Arr+12 );
  __m128d ArrHi3 = _mm_load_pd( Arr+14 );


  __m128d resLo = _mm_mul_pd( V0, ArrLo0 );
  __m128d resHi = _mm_mul_pd( V0, ArrHi0 );

  resLo = _mm_add_pd( resLo, _mm_mul_pd( V1, ArrLo1 ));
  resHi = _mm_add_pd( resHi, _mm_mul_pd( V1, ArrHi1 ));

  resLo = _mm_add_pd( resLo, _mm_mul_pd( V2, ArrLo2 ));
  resHi = _mm_add_pd( resHi, _mm_mul_pd( V2, ArrHi2 ));

  resLo = _mm_add_pd( resLo, _mm_mul_pd( V3, ArrLo3 ));
  resHi = _mm_add_pd( resHi, _mm_mul_pd( V3, ArrHi3 ));

  return vector4D( resLo, resHi );

#elif defined (VECTOR_IMPL_AVX)

  // load the vector into registers:
  // (V0 consist of 4 copies of component 0, V1 has 4 copies of
  // component 1 and so on)
  __m256d V0 = _mm256_broadcast_sd(x.Mem);
  __m256d V1 = _mm256_broadcast_sd(x.Mem+1);
  __m256d V2 = _mm256_broadcast_sd(x.Mem+2);
  __m256d V3 = _mm256_broadcast_sd(x.Mem+3);

  // load the matrix components into registers:
  __m256d Arr0 = _mm256_load_pd( Arr+0 );
  __m256d Arr1 = _mm256_load_pd( Arr+4 );
  __m256d Arr2 = _mm256_load_pd( Arr+8 );
  __m256d Arr3 = _mm256_load_pd( Arr+12 );

  // do the actual multiplication:
  __m256d res = _mm256_mul_pd( V0, Arr0 );
  res = _mm256_add_pd( res, _mm256_mul_pd( V1, Arr1 ));
  res = _mm256_add_pd( res, _mm256_mul_pd( V2, Arr2 ));
  res = _mm256_add_pd( res, _mm256_mul_pd( V3, Arr3 ));
  
  return vector4D( res );

#else

  return boostInvScalar( x );

#endif
}

vector4D lorentz::boostInvScalar( const vector4D & x ) const
{
  return vector4D( Arr[0+4*0]*x[0]+Arr[1+4*0]*x[1]+Arr[2+4*0]*x[2]+Arr[3+4*0]*x[3],
                   Arr[0+4*3]*x[0]+Arr[1+4*3]*x[1]+Arr[2+4*3]*x[2]+Arr[3+4*3]*x[3],
                   Arr[0+4*2]*x[0]+Arr[1+4*2]*x[1]+Arr[2+4*2]*x[2]+Arr[3+4*2]*x[3],
                   Arr[0+4*1]*x[0]+Arr[1+4*1]*x[1]+Arr[2+4*1]*x[2]+Arr[3+4*1]*x[3]);
}


///// BOOST TWO VECTORS /////


void lorentz::boost( const vector4D & x1, const vector4D & x2, 
                     vector4D & x1New, vector4D & x2New ) const
{
#if defined (VECTOR_IMPL_SSE)
  // please note: we have to be careful with the signs of the matrix:
  // The components [0][i] and [i][0] have to be negated!

  // we also use _mm_prefetch in order to load some memory into a
  // cache closer to the processor.

  _MM_PREFETCH(x1.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2.Mem, _MM_HINT_NTA);

  SSE_ALIGNED(const unsigned long long) S = 0x8000000000000000ull ;

  // cppcheck-suppress invalidPointerCast
  __m128d xS10 = _mm_load_sd( (double *)&S ); // only lower part negates
  __m128d xS01 = _mm_setzero_pd( ); // set it to zero
  __m128d xS11 = _mm_unpacklo_pd( xS10, xS10 ); // both parts negate
  xS01 = _mm_unpacklo_pd( xS01, xS10 ); // only upper part negates

  _MM_PREFETCH(Arr, _MM_HINT_NTA);

  __m128d V1a = _mm_load_pd(x1.Mem);
  __m128d V1b = _mm_load_pd(x2.Mem);
  __m128d V3a = _mm_load_pd(x1.Mem+2);
  __m128d V3b = _mm_load_pd(x2.Mem+2);
  __m128d V0a = _mm_unpacklo_pd( V1a, V1a );
  __m128d V0b = _mm_unpacklo_pd( V1b, V1b );
  __m128d V2a = _mm_unpacklo_pd( V3a, V3a );
  __m128d V2b = _mm_unpacklo_pd( V3b, V3b );
  V1a = _mm_unpackhi_pd( V1a, V1a );
  V1b = _mm_unpackhi_pd( V1b, V1b );
  V3a = _mm_unpackhi_pd( V3a, V3a );
  V3b = _mm_unpackhi_pd( V3b, V3b );

  // load the matrix components into SSE registers:
  __m128d ArrLo0 = _mm_load_pd( Arr+0 );
  __m128d ArrHi0 = _mm_load_pd( Arr+2 );
  __m128d ArrLo1 = _mm_load_pd( Arr+4 );
  __m128d ArrHi1 = _mm_load_pd( Arr+6 );
  __m128d ArrLo2 = _mm_load_pd( Arr+8 );
  __m128d ArrHi2 = _mm_load_pd( Arr+10 );
  __m128d ArrLo3 = _mm_load_pd( Arr+12 );
  __m128d ArrHi3 = _mm_load_pd( Arr+14 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);

  // correct the signs of the matrix components:
  ArrLo0 = _mm_xor_pd( ArrLo0, xS01 );
  ArrHi0 = _mm_xor_pd( ArrHi0, xS11 );
  ArrLo1 = _mm_xor_pd( ArrLo1, xS10 );
  ArrLo2 = _mm_xor_pd( ArrLo2, xS10 );
  ArrLo3 = _mm_xor_pd( ArrLo3, xS10 );

  __m128d resLo1 = _mm_mul_pd( V0a, ArrLo0 );
  __m128d resHi1 = _mm_mul_pd( V0a, ArrHi0 );
  __m128d resLo2 = _mm_mul_pd( V0b, ArrLo0 );
  __m128d resHi2 = _mm_mul_pd( V0b, ArrHi0 );

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V1a, ArrLo1 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V1a, ArrHi1 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V1b, ArrLo1 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V1b, ArrHi1 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V2a, ArrLo2 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V2a, ArrHi2 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V2b, ArrLo2 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V2b, ArrHi2 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V3a, ArrLo3 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V3a, ArrHi3 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V3b, ArrLo3 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V3b, ArrHi3 ));

  x1New = vector4D( resLo1, resHi1 );
  x2New = vector4D( resLo2, resHi2 );


#elif defined (VECTOR_IMPL_AVX)

  _MM_PREFETCH(Arr, _MM_HINT_NTA);
  _MM_PREFETCH(x1.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2.Mem, _MM_HINT_NTA);

  // load the vector into registers:
  // (V0 consist of 4 copies of component 0, V1 has 4 copies of
  // component 1 and so on)
  __m256d V0a = _mm256_broadcast_sd(x1.Mem);
  __m256d V1a = _mm256_broadcast_sd(x1.Mem+1);
  __m256d V2a = _mm256_broadcast_sd(x1.Mem+2);
  __m256d V3a = _mm256_broadcast_sd(x1.Mem+3);

  __m256d V0b = _mm256_broadcast_sd(x2.Mem);
  __m256d V1b = _mm256_broadcast_sd(x2.Mem+1);
  __m256d V2b = _mm256_broadcast_sd(x2.Mem+2);
  __m256d V3b = _mm256_broadcast_sd(x2.Mem+3);

  // load the matrix components into registers:
  __m256d Arr0 = _mm256_load_pd( Arr+0 );
  __m256d Arr1 = _mm256_load_pd( Arr+4 );
  __m256d Arr2 = _mm256_load_pd( Arr+8 );
  __m256d Arr3 = _mm256_load_pd( Arr+12 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);

  // correct the signs of the matrix components:
  // in Arr0, we have to flip signs .123
  // in Arri, we have to flip signs 0...

  static const __m256d S01 = _mm256_set_pd(-0.,-0.,-0., 0.);
  static const __m256d S10 = _mm256_set_pd( 0., 0., 0.,-0.);

  Arr0 = _mm256_xor_pd( Arr0, S01 );
  Arr1 = _mm256_xor_pd( Arr1, S10 );
  Arr2 = _mm256_xor_pd( Arr2, S10 );
  Arr3 = _mm256_xor_pd( Arr3, S10 );


  // do the actual multiplication:
  __m256d res1 = _mm256_mul_pd( V0a, Arr0 );
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V1a, Arr1 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V2a, Arr2 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V3a, Arr3 ));
  
  __m256d res2 = _mm256_mul_pd( V0b, Arr0 );
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V1b, Arr1 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V2b, Arr2 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V3b, Arr3 ));
  
  x1New = vector4D( res1 );
  x2New = vector4D( res2 );

#else

  boostScalar( x1, x2, x1New, x2New );

#endif

}

void lorentz::boostScalar( const vector4D & x1, const vector4D & x2, 
                     vector4D & x1New, vector4D & x2New ) const
{
  x1New = vector4D( Arr[0+4*0]*x1[0]-Arr[1+4*0]*x1[1]-Arr[2+4*0]*x1[2]-Arr[3+4*0]*x1[3],
                   -Arr[0+4*3]*x1[0]+Arr[1+4*3]*x1[1]+Arr[2+4*3]*x1[2]+Arr[3+4*3]*x1[3],
                   -Arr[0+4*2]*x1[0]+Arr[1+4*2]*x1[1]+Arr[2+4*2]*x1[2]+Arr[3+4*2]*x1[3],
                   -Arr[0+4*1]*x1[0]+Arr[1+4*1]*x1[1]+Arr[2+4*1]*x1[2]+Arr[3+4*1]*x1[3]);
  
  x2New = vector4D( Arr[0+4*0]*x2[0]-Arr[1+4*0]*x2[1]-Arr[2+4*0]*x2[2]-Arr[3+4*0]*x2[3],
                   -Arr[0+4*3]*x2[0]+Arr[1+4*3]*x2[1]+Arr[2+4*3]*x2[2]+Arr[3+4*3]*x2[3],
                   -Arr[0+4*2]*x2[0]+Arr[1+4*2]*x2[1]+Arr[2+4*2]*x2[2]+Arr[3+4*2]*x2[3],
                   -Arr[0+4*1]*x2[0]+Arr[1+4*1]*x2[1]+Arr[2+4*1]*x2[2]+Arr[3+4*1]*x2[3]);

}



double lorentz::getSpatialDistance(const vector4D& x1, const vector4D& x2)
{
  double xsq = pow(x1[1]-x2[1],2.0);
  double ysq = pow(x1[2]-x2[2],2.0);
  double zsq = pow(x1[3]-x2[3],2.0);
  return sqrt(  xsq + ysq + zsq );
}


double lorentz::getDCAsquared(const vector4D& x1, const vector4D& x2, const vector4D& p1, const vector4D& p2)
{
  double dvx,dvy,dvz,dx,dy,dz;
  dvx = p1[1]/p1[0] - p2[1]/p2[0];
  dvy = p1[2]/p1[0] - p2[2]/p2[0];
  dvz = p1[3]/p1[0] - p2[3]/p2[0];
  
  double dv2 = dvx*dvx +dvy*dvy + dvz*dvz;
  
  
  dx = x1[1]-x2[1];
  dy = x1[2]-x2[2];
  dz = x1[3]-x2[3];
  
  double dx2 = dx*dx + dy*dy + dz*dz;
  double dvdx = dx*dvx + dy*dvy + dz*dvz;
  
  double t=-(dvdx)/dv2;
  return (dx2 - dvdx*dvdx/dv2);
}






void lorentz::boostInv( const vector4D & x1, const vector4D & x2, 
                        vector4D & x1New, vector4D & x2New ) const
{
#if defined (VECTOR_IMPL_SSE)

  __m128d V1a = _mm_load_pd(x1.Mem);
  __m128d V1b = _mm_load_pd(x2.Mem);
  __m128d V3a = _mm_load_pd(x1.Mem+2);
  __m128d V3b = _mm_load_pd(x2.Mem+2);
  __m128d V0a = _mm_unpacklo_pd( V1a, V1a );
  __m128d V0b = _mm_unpacklo_pd( V1b, V1b );
  __m128d V2a = _mm_unpacklo_pd( V3a, V3a );
  __m128d V2b = _mm_unpacklo_pd( V3b, V3b );
  V1a = _mm_unpackhi_pd( V1a, V1a );
  V1b = _mm_unpackhi_pd( V1b, V1b );
  V3a = _mm_unpackhi_pd( V3a, V3a );
  V3b = _mm_unpackhi_pd( V3b, V3b );

  // load the matrix components into SSE registers:
  __m128d ArrLo0 = _mm_load_pd( Arr+0 );
  __m128d ArrHi0 = _mm_load_pd( Arr+2 );
  __m128d ArrLo1 = _mm_load_pd( Arr+4 );
  __m128d ArrHi1 = _mm_load_pd( Arr+6 );
  __m128d ArrLo2 = _mm_load_pd( Arr+8 );
  __m128d ArrHi2 = _mm_load_pd( Arr+10 );
  __m128d ArrLo3 = _mm_load_pd( Arr+12 );
  __m128d ArrHi3 = _mm_load_pd( Arr+14 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);

  __m128d resLo1 = _mm_mul_pd( V0a, ArrLo0 );
  __m128d resHi1 = _mm_mul_pd( V0a, ArrHi0 );
  __m128d resLo2 = _mm_mul_pd( V0b, ArrLo0 );
  __m128d resHi2 = _mm_mul_pd( V0b, ArrHi0 );

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V1a, ArrLo1 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V1a, ArrHi1 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V1b, ArrLo1 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V1b, ArrHi1 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V2a, ArrLo2 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V2a, ArrHi2 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V2b, ArrLo2 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V2b, ArrHi2 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V3a, ArrLo3 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V3a, ArrHi3 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V3b, ArrLo3 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V3b, ArrHi3 ));

  x1New = vector4D( resLo1, resHi1 );
  x2New = vector4D( resLo2, resHi2 );


 #elif defined (VECTOR_IMPL_AVX)

  _MM_PREFETCH(Arr, _MM_HINT_NTA);
  _MM_PREFETCH(x1.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2.Mem, _MM_HINT_NTA);

  // load the vector into registers:
  // (V0 consist of 4 copies of component 0, V1 has 4 copies of
  // component 1 and so on)
  __m256d V0a = _mm256_broadcast_sd(x1.Mem);
  __m256d V1a = _mm256_broadcast_sd(x1.Mem+1);
  __m256d V2a = _mm256_broadcast_sd(x1.Mem+2);
  __m256d V3a = _mm256_broadcast_sd(x1.Mem+3);

  __m256d V0b = _mm256_broadcast_sd(x2.Mem);
  __m256d V1b = _mm256_broadcast_sd(x2.Mem+1);
  __m256d V2b = _mm256_broadcast_sd(x2.Mem+2);
  __m256d V3b = _mm256_broadcast_sd(x2.Mem+3);

  // load the matrix components into registers:
  __m256d Arr0 = _mm256_load_pd( Arr+0 );
  __m256d Arr1 = _mm256_load_pd( Arr+4 );
  __m256d Arr2 = _mm256_load_pd( Arr+8 );
  __m256d Arr3 = _mm256_load_pd( Arr+12 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);

  // do the actual multiplication:
  __m256d res1 = _mm256_mul_pd( V0a, Arr0 );
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V1a, Arr1 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V2a, Arr2 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V3a, Arr3 ));
  
  __m256d res2 = _mm256_mul_pd( V0b, Arr0 );
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V1b, Arr1 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V2b, Arr2 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V3b, Arr3 ));
  
  x1New = vector4D( res1 );
  x2New = vector4D( res2 );
 
#else

  boostInvScalar( x1, x2, x1New, x2New );

#endif
}

void lorentz::boostInvScalar( const vector4D & x1, const vector4D & x2, 
                              vector4D & x1New, vector4D & x2New ) const
{
  x1New = vector4D(Arr[0+4*0]*x1[0]+Arr[1+4*0]*x1[1]+Arr[2+4*0]*x1[2]+Arr[3+4*0]*x1[3],
                   Arr[0+4*3]*x1[0]+Arr[1+4*3]*x1[1]+Arr[2+4*3]*x1[2]+Arr[3+4*3]*x1[3],
                   Arr[0+4*2]*x1[0]+Arr[1+4*2]*x1[1]+Arr[2+4*2]*x1[2]+Arr[3+4*2]*x1[3],
                   Arr[0+4*1]*x1[0]+Arr[1+4*1]*x1[1]+Arr[2+4*1]*x1[2]+Arr[3+4*1]*x1[3]);
  
  x2New = vector4D(Arr[0+4*0]*x2[0]+Arr[1+4*0]*x2[1]+Arr[2+4*0]*x2[2]+Arr[3+4*0]*x2[3],
                   Arr[0+4*3]*x2[0]+Arr[1+4*3]*x2[1]+Arr[2+4*3]*x2[2]+Arr[3+4*3]*x2[3],
                   Arr[0+4*2]*x2[0]+Arr[1+4*2]*x2[1]+Arr[2+4*2]*x2[2]+Arr[3+4*2]*x2[3],
                   Arr[0+4*1]*x2[0]+Arr[1+4*1]*x2[1]+Arr[2+4*1]*x2[2]+Arr[3+4*1]*x2[3]);
}


///// BOOST THREE VECTORS /////



void lorentz::boost( const vector4D & x1, const vector4D & x2, const vector4D & x3,
                     vector4D & x1New, vector4D & x2New, vector4D & x3New ) const
{
#if defined (VECTOR_IMPL_SSE)
  // please note: we have to be careful with the signs of the matrix:
  // The components [0][i] and [i][0] have to be negated!

  // we also use _mm_prefetch in order to load some memory into a
  // cache closer to the processor.

  _MM_PREFETCH(x1.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3.Mem, _MM_HINT_NTA);

  SSE_ALIGNED(const unsigned long long) S = 0x8000000000000000ull ;

  // cppcheck-suppress invalidPointerCast
  __m128d xS10 = _mm_load_sd( (double *)&S ); // only lower part negates
  __m128d xS01 = _mm_setzero_pd( ); // set it to zero
  __m128d xS11 = _mm_unpacklo_pd( xS10, xS10 ); // both parts negate
  xS01 = _mm_unpacklo_pd( xS01, xS10 ); // only upper part negates

  _MM_PREFETCH(Arr, _MM_HINT_NTA);

  __m128d V1a = _mm_load_pd(x1.Mem);
  __m128d V1b = _mm_load_pd(x2.Mem);
  __m128d V1c = _mm_load_pd(x3.Mem);
  __m128d V3a = _mm_load_pd(x1.Mem+2);
  __m128d V3b = _mm_load_pd(x2.Mem+2);
  __m128d V3c = _mm_load_pd(x3.Mem+2);
  __m128d V0a = _mm_unpacklo_pd( V1a, V1a );
  __m128d V0b = _mm_unpacklo_pd( V1b, V1b );
  __m128d V0c = _mm_unpacklo_pd( V1c, V1c );
  __m128d V2a = _mm_unpacklo_pd( V3a, V3a );
  __m128d V2b = _mm_unpacklo_pd( V3b, V3b );
  __m128d V2c = _mm_unpacklo_pd( V3c, V3c );
  V1a = _mm_unpackhi_pd( V1a, V1a );
  V1b = _mm_unpackhi_pd( V1b, V1b );
  V1c = _mm_unpackhi_pd( V1c, V1c );
  V3a = _mm_unpackhi_pd( V3a, V3a );
  V3b = _mm_unpackhi_pd( V3b, V3b );
  V3c = _mm_unpackhi_pd( V3c, V3c );

  // load the matrix components into SSE registers:
  __m128d ArrLo0 = _mm_load_pd( Arr+0 );
  __m128d ArrHi0 = _mm_load_pd( Arr+2 );
  __m128d ArrLo1 = _mm_load_pd( Arr+4 );
  __m128d ArrHi1 = _mm_load_pd( Arr+6 );
  __m128d ArrLo2 = _mm_load_pd( Arr+8 );
  __m128d ArrHi2 = _mm_load_pd( Arr+10 );
  __m128d ArrLo3 = _mm_load_pd( Arr+12 );
  __m128d ArrHi3 = _mm_load_pd( Arr+14 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3New.Mem, _MM_HINT_NTA);

  // correct the signs of the matrix components:
  ArrLo0 = _mm_xor_pd( ArrLo0, xS01 );
  ArrHi0 = _mm_xor_pd( ArrHi0, xS11 );
  ArrLo1 = _mm_xor_pd( ArrLo1, xS10 );
  ArrLo2 = _mm_xor_pd( ArrLo2, xS10 );
  ArrLo3 = _mm_xor_pd( ArrLo3, xS10 );

  __m128d resLo1 = _mm_mul_pd( V0a, ArrLo0 );
  __m128d resHi1 = _mm_mul_pd( V0a, ArrHi0 );
  __m128d resLo2 = _mm_mul_pd( V0b, ArrLo0 );
  __m128d resHi2 = _mm_mul_pd( V0b, ArrHi0 );
  __m128d resLo3 = _mm_mul_pd( V0c, ArrLo0 );
  __m128d resHi3 = _mm_mul_pd( V0c, ArrHi0 );

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V1a, ArrLo1 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V1a, ArrHi1 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V1b, ArrLo1 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V1b, ArrHi1 ));
  resLo3 = _mm_add_pd( resLo3, _mm_mul_pd( V1c, ArrLo1 ));
  resHi3 = _mm_add_pd( resHi3, _mm_mul_pd( V1c, ArrHi1 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V2a, ArrLo2 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V2a, ArrHi2 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V2b, ArrLo2 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V2b, ArrHi2 ));
  resLo3 = _mm_add_pd( resLo3, _mm_mul_pd( V2c, ArrLo2 ));
  resHi3 = _mm_add_pd( resHi3, _mm_mul_pd( V2c, ArrHi2 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V3a, ArrLo3 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V3a, ArrHi3 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V3b, ArrLo3 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V3b, ArrHi3 ));
  resLo3 = _mm_add_pd( resLo3, _mm_mul_pd( V3c, ArrLo3 ));
  resHi3 = _mm_add_pd( resHi3, _mm_mul_pd( V3c, ArrHi3 ));

  x1New = vector4D( resLo1, resHi1 );
  x2New = vector4D( resLo2, resHi2 );
  x3New = vector4D( resLo3, resHi3 );


#elif defined (VECTOR_IMPL_AVX)

  _MM_PREFETCH(Arr, _MM_HINT_NTA);
  _MM_PREFETCH(x1.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3.Mem, _MM_HINT_NTA);

  // load the vector into registers:
  // (V0 consist of 4 copies of component 0, V1 has 4 copies of
  // component 1 and so on)
  __m256d V0a = _mm256_broadcast_sd(x1.Mem);
  __m256d V1a = _mm256_broadcast_sd(x1.Mem+1);
  __m256d V2a = _mm256_broadcast_sd(x1.Mem+2);
  __m256d V3a = _mm256_broadcast_sd(x1.Mem+3);

  __m256d V0b = _mm256_broadcast_sd(x2.Mem);
  __m256d V1b = _mm256_broadcast_sd(x2.Mem+1);
  __m256d V2b = _mm256_broadcast_sd(x2.Mem+2);
  __m256d V3b = _mm256_broadcast_sd(x2.Mem+3);

  __m256d V0c = _mm256_broadcast_sd(x3.Mem);
  __m256d V1c = _mm256_broadcast_sd(x3.Mem+1);
  __m256d V2c = _mm256_broadcast_sd(x3.Mem+2);
  __m256d V3c = _mm256_broadcast_sd(x3.Mem+3);

  // load the matrix components into registers:
  __m256d Arr0 = _mm256_load_pd( Arr+0 );
  __m256d Arr1 = _mm256_load_pd( Arr+4 );
  __m256d Arr2 = _mm256_load_pd( Arr+8 );
  __m256d Arr3 = _mm256_load_pd( Arr+12 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3New.Mem, _MM_HINT_NTA);

  // correct the signs of the matrix components:
  // in Arr0, we have to flip signs .123
  // in Arri, we have to flip signs 0...

  static const __m256d S01 = _mm256_set_pd(-0.,-0.,-0., 0.);
  static const __m256d S10 = _mm256_set_pd( 0., 0., 0.,-0.);

  Arr0 = _mm256_xor_pd( Arr0, S01 );
  Arr1 = _mm256_xor_pd( Arr1, S10 );
  Arr2 = _mm256_xor_pd( Arr2, S10 );
  Arr3 = _mm256_xor_pd( Arr3, S10 );


  // do the actual multiplication:
  __m256d res1 = _mm256_mul_pd( V0a, Arr0 );
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V1a, Arr1 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V2a, Arr2 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V3a, Arr3 ));
  
  __m256d res2 = _mm256_mul_pd( V0b, Arr0 );
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V1b, Arr1 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V2b, Arr2 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V3b, Arr3 ));

  __m256d res3 = _mm256_mul_pd( V0c, Arr0 );
  res3 = _mm256_add_pd( res3, _mm256_mul_pd( V1c, Arr1 ));
  res3 = _mm256_add_pd( res3, _mm256_mul_pd( V2c, Arr2 ));
  res3 = _mm256_add_pd( res3, _mm256_mul_pd( V3c, Arr3 ));
  
  x1New = vector4D( res1 );
  x2New = vector4D( res2 );
  x3New = vector4D( res3 );
  
#else

  boostScalar( x1, x2, x3, x1New, x2New, x3New );

#endif
}

void lorentz::boostScalar( const vector4D & x1, const vector4D & x2, const vector4D & x3,
                           vector4D & x1New, vector4D & x2New, vector4D & x3New ) const
{
  x1New = vector4D( Arr[0+4*0]*x1[0]-Arr[1+4*0]*x1[1]-Arr[2+4*0]*x1[2]-Arr[3+4*0]*x1[3],
                   -Arr[0+4*3]*x1[0]+Arr[1+4*3]*x1[1]+Arr[2+4*3]*x1[2]+Arr[3+4*3]*x1[3],
                   -Arr[0+4*2]*x1[0]+Arr[1+4*2]*x1[1]+Arr[2+4*2]*x1[2]+Arr[3+4*2]*x1[3],
                   -Arr[0+4*1]*x1[0]+Arr[1+4*1]*x1[1]+Arr[2+4*1]*x1[2]+Arr[3+4*1]*x1[3]);
  
  x2New = vector4D( Arr[0+4*0]*x2[0]-Arr[1+4*0]*x2[1]-Arr[2+4*0]*x2[2]-Arr[3+4*0]*x2[3],
                   -Arr[0+4*3]*x2[0]+Arr[1+4*3]*x2[1]+Arr[2+4*3]*x2[2]+Arr[3+4*3]*x2[3],
                   -Arr[0+4*2]*x2[0]+Arr[1+4*2]*x2[1]+Arr[2+4*2]*x2[2]+Arr[3+4*2]*x2[3],
                   -Arr[0+4*1]*x2[0]+Arr[1+4*1]*x2[1]+Arr[2+4*1]*x2[2]+Arr[3+4*1]*x2[3]);

  x3New = vector4D( Arr[0+4*0]*x3[0]-Arr[1+4*0]*x3[1]-Arr[2+4*0]*x3[2]-Arr[3+4*0]*x3[3],
                   -Arr[0+4*3]*x3[0]+Arr[1+4*3]*x3[1]+Arr[2+4*3]*x3[2]+Arr[3+4*3]*x3[3],
                   -Arr[0+4*2]*x3[0]+Arr[1+4*2]*x3[1]+Arr[2+4*2]*x3[2]+Arr[3+4*2]*x3[3],
                   -Arr[0+4*1]*x3[0]+Arr[1+4*1]*x3[1]+Arr[2+4*1]*x3[2]+Arr[3+4*1]*x3[3]);

}

void lorentz::boostInv( const vector4D & x1, const vector4D & x2, const vector4D & x3, 
                        vector4D & x1New, vector4D & x2New, vector4D & x3New ) const
{
#if defined (VECTOR_IMPL_SSE)
  // we also use _mm_prefetch in order to load some memory into a
  // cache closer to the processor.

  _MM_PREFETCH(x1.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(Arr, _MM_HINT_NTA);

  __m128d V1a = _mm_load_pd(x1.Mem);
  __m128d V1b = _mm_load_pd(x2.Mem);
  __m128d V1c = _mm_load_pd(x3.Mem);
  __m128d V3a = _mm_load_pd(x1.Mem+2);
  __m128d V3b = _mm_load_pd(x2.Mem+2);
  __m128d V3c = _mm_load_pd(x3.Mem+2);
  __m128d V0a = _mm_unpacklo_pd( V1a, V1a );
  __m128d V0b = _mm_unpacklo_pd( V1b, V1b );
  __m128d V0c = _mm_unpacklo_pd( V1c, V1c );
  __m128d V2a = _mm_unpacklo_pd( V3a, V3a );
  __m128d V2b = _mm_unpacklo_pd( V3b, V3b );
  __m128d V2c = _mm_unpacklo_pd( V3c, V3c );
  V1a = _mm_unpackhi_pd( V1a, V1a );
  V1b = _mm_unpackhi_pd( V1b, V1b );
  V1c = _mm_unpackhi_pd( V1c, V1c );
  V3a = _mm_unpackhi_pd( V3a, V3a );
  V3b = _mm_unpackhi_pd( V3b, V3b );
  V3c = _mm_unpackhi_pd( V3c, V3c );

  // load the matrix components into SSE registers:
  __m128d ArrLo0 = _mm_load_pd( Arr+0 );
  __m128d ArrHi0 = _mm_load_pd( Arr+2 );
  __m128d ArrLo1 = _mm_load_pd( Arr+4 );
  __m128d ArrHi1 = _mm_load_pd( Arr+6 );
  __m128d ArrLo2 = _mm_load_pd( Arr+8 );
  __m128d ArrHi2 = _mm_load_pd( Arr+10 );
  __m128d ArrLo3 = _mm_load_pd( Arr+12 );
  __m128d ArrHi3 = _mm_load_pd( Arr+14 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3New.Mem, _MM_HINT_NTA);

  __m128d resLo1 = _mm_mul_pd( V0a, ArrLo0 );
  __m128d resHi1 = _mm_mul_pd( V0a, ArrHi0 );
  __m128d resLo2 = _mm_mul_pd( V0b, ArrLo0 );
  __m128d resHi2 = _mm_mul_pd( V0b, ArrHi0 );
  __m128d resLo3 = _mm_mul_pd( V0c, ArrLo0 );
  __m128d resHi3 = _mm_mul_pd( V0c, ArrHi0 );

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V1a, ArrLo1 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V1a, ArrHi1 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V1b, ArrLo1 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V1b, ArrHi1 ));
  resLo3 = _mm_add_pd( resLo3, _mm_mul_pd( V1c, ArrLo1 ));
  resHi3 = _mm_add_pd( resHi3, _mm_mul_pd( V1c, ArrHi1 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V2a, ArrLo2 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V2a, ArrHi2 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V2b, ArrLo2 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V2b, ArrHi2 ));
  resLo3 = _mm_add_pd( resLo3, _mm_mul_pd( V2c, ArrLo2 ));
  resHi3 = _mm_add_pd( resHi3, _mm_mul_pd( V2c, ArrHi2 ));

  resLo1 = _mm_add_pd( resLo1, _mm_mul_pd( V3a, ArrLo3 ));
  resHi1 = _mm_add_pd( resHi1, _mm_mul_pd( V3a, ArrHi3 ));
  resLo2 = _mm_add_pd( resLo2, _mm_mul_pd( V3b, ArrLo3 ));
  resHi2 = _mm_add_pd( resHi2, _mm_mul_pd( V3b, ArrHi3 ));
  resLo3 = _mm_add_pd( resLo3, _mm_mul_pd( V3c, ArrLo3 ));
  resHi3 = _mm_add_pd( resHi3, _mm_mul_pd( V3c, ArrHi3 ));

  x1New = vector4D( resLo1, resHi1 );
  x2New = vector4D( resLo2, resHi2 );
  x3New = vector4D( resLo3, resHi3 );
  
#elif defined (VECTOR_IMPL_AVX)

  _MM_PREFETCH(Arr, _MM_HINT_NTA);
  _MM_PREFETCH(x1.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3.Mem, _MM_HINT_NTA);

  // load the vector into registers:
  // (V0 consist of 4 copies of component 0, V1 has 4 copies of
  // component 1 and so on)
  __m256d V0a = _mm256_broadcast_sd(x1.Mem);
  __m256d V1a = _mm256_broadcast_sd(x1.Mem+1);
  __m256d V2a = _mm256_broadcast_sd(x1.Mem+2);
  __m256d V3a = _mm256_broadcast_sd(x1.Mem+3);

  __m256d V0b = _mm256_broadcast_sd(x2.Mem);
  __m256d V1b = _mm256_broadcast_sd(x2.Mem+1);
  __m256d V2b = _mm256_broadcast_sd(x2.Mem+2);
  __m256d V3b = _mm256_broadcast_sd(x2.Mem+3);

  __m256d V0c = _mm256_broadcast_sd(x3.Mem);
  __m256d V1c = _mm256_broadcast_sd(x3.Mem+1);
  __m256d V2c = _mm256_broadcast_sd(x3.Mem+2);
  __m256d V3c = _mm256_broadcast_sd(x3.Mem+3);

  // load the matrix components into registers:
  __m256d Arr0 = _mm256_load_pd( Arr+0 );
  __m256d Arr1 = _mm256_load_pd( Arr+4 );
  __m256d Arr2 = _mm256_load_pd( Arr+8 );
  __m256d Arr3 = _mm256_load_pd( Arr+12 );

  // already load the output memory into cache
  _MM_PREFETCH(x1New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x2New.Mem, _MM_HINT_NTA);
  _MM_PREFETCH(x3New.Mem, _MM_HINT_NTA);

  // do the actual multiplication:
  __m256d res1 = _mm256_mul_pd( V0a, Arr0 );
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V1a, Arr1 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V2a, Arr2 ));
  res1 = _mm256_add_pd( res1, _mm256_mul_pd( V3a, Arr3 ));
  
  __m256d res2 = _mm256_mul_pd( V0b, Arr0 );
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V1b, Arr1 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V2b, Arr2 ));
  res2 = _mm256_add_pd( res2, _mm256_mul_pd( V3b, Arr3 ));

  __m256d res3 = _mm256_mul_pd( V0c, Arr0 );
  res3 = _mm256_add_pd( res3, _mm256_mul_pd( V1c, Arr1 ));
  res3 = _mm256_add_pd( res3, _mm256_mul_pd( V2c, Arr2 ));
  res3 = _mm256_add_pd( res3, _mm256_mul_pd( V3c, Arr3 ));
  
  x1New = vector4D( res1 );
  x2New = vector4D( res2 );
  x3New = vector4D( res3 );

#else

  boostInvScalar( x1, x2, x3, x1New, x2New, x3New );

#endif
}

void lorentz::boostInvScalar( const vector4D & x1, const vector4D & x2, const vector4D & x3, 
                              vector4D & x1New, vector4D & x2New, vector4D & x3New ) const
{
  x1New = vector4D(Arr[0+4*0]*x1[0]+Arr[1+4*0]*x1[1]+Arr[2+4*0]*x1[2]+Arr[3+4*0]*x1[3],
                   Arr[0+4*3]*x1[0]+Arr[1+4*3]*x1[1]+Arr[2+4*3]*x1[2]+Arr[3+4*3]*x1[3],
                   Arr[0+4*2]*x1[0]+Arr[1+4*2]*x1[1]+Arr[2+4*2]*x1[2]+Arr[3+4*2]*x1[3],
                   Arr[0+4*1]*x1[0]+Arr[1+4*1]*x1[1]+Arr[2+4*1]*x1[2]+Arr[3+4*1]*x1[3]);
  
  x2New = vector4D(Arr[0+4*0]*x2[0]+Arr[1+4*0]*x2[1]+Arr[2+4*0]*x2[2]+Arr[3+4*0]*x2[3],
                   Arr[0+4*3]*x2[0]+Arr[1+4*3]*x2[1]+Arr[2+4*3]*x2[2]+Arr[3+4*3]*x2[3],
                   Arr[0+4*2]*x2[0]+Arr[1+4*2]*x2[1]+Arr[2+4*2]*x2[2]+Arr[3+4*2]*x2[3],
                   Arr[0+4*1]*x2[0]+Arr[1+4*1]*x2[1]+Arr[2+4*1]*x2[2]+Arr[3+4*1]*x2[3]);

  x3New = vector4D(Arr[0+4*0]*x3[0]+Arr[1+4*0]*x3[1]+Arr[2+4*0]*x3[2]+Arr[3+4*0]*x3[3],
                   Arr[0+4*3]*x3[0]+Arr[1+4*3]*x3[1]+Arr[2+4*3]*x3[2]+Arr[3+4*3]*x3[3],
                   Arr[0+4*2]*x3[0]+Arr[1+4*2]*x3[1]+Arr[2+4*2]*x3[2]+Arr[3+4*2]*x3[3],
                   Arr[0+4*1]*x3[0]+Arr[1+4*1]*x3[1]+Arr[2+4*1]*x3[2]+Arr[3+4*1]*x3[3]);

}


// #ifdef FOR_HISTORICAL_REASONS


// /**
//  * Given a velocity vector beta, this routine boosts the vector A with this velocity.
//  *
//  * @param[in] beta[] Boost velocity vector.
//  * @param[in] A[] 4-Vector to be boosted.
//  * @param[out] LA[] The result - A boosted with beta
//  */
// void lorentz(const double beta[4], const double A[4], double LA[4])
// {
//   double beta2,cc,c0,c1,c2,c3,gama;

//   beta2=0.0;
//   for(int i=1;i<=3;i++)
//     beta2 += beta[i]*beta[i];

//   if(beta2 < 1.0e-10)
//   {
//     for(int k=0;k<=3;k++) 
//       LA[k]=A[k];
//     return;
//   }

//   if(beta2 > (1.0-1.0e-10)) 
//   {
//     std::cout << beta2  <<"  error in lorentz()" << std::endl;
//     std::cout << A[0] << " " << A[1] << " " << A[2] << " " << A[3] << std::endl;
//     std::cout << beta[1] << " " << beta[2] << " "  << beta[3] << "  " << std::endl;
//     beta2 = (1.0-1.0e-10);
//     std::cout << "Set beta2 to  " << beta2 << std::endl;
//   }
//   else
//   {
//     gama=1.0/sqrt(1.0-beta2);
//       //cout << beta2 << endl;
//       //cout << A[0] << " " << A[1] << " " << A[2] << " " << A[3] << endl;
//   }

//   LA[0]=gama*(A[0]-beta[1]*A[1]-beta[2]*A[2]-beta[3]*A[3]);
//   for(int j=1;j<=3;j++){
//     cc=(gama-1.0)*beta[j]/beta2;
//     c0=-gama*beta[j];
//     c1=cc*beta[1];
//     c2=cc*beta[2];
//     c3=cc*beta[3];
//     LA[j]=c0*A[0]+c1*A[1]+c2*A[2]+c3*A[3]+A[j];
//   }
// }

// #endif
