//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/dSFMT.hpp $
//$LastChangedDate: 2014-05-13 19:11:47 +0200 (Di, 13. Mai 2014) $
//$LastChangedRevision: 1683 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef DSFMT_HPP
#define DSFMT_HPP

#include <iostream>
#include <iomanip>


#if defined(HAVE_SSE2)
#include <emmintrin.h>

/** 128-bit data structure */
union W128_T 
{
    __m128i si;
    __m128d sd;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};

union X128I_T 
{
    uint64_t u[2];
    __m128i  i128;
    X128I_T(uint64_t u1,uint64_t u2) {u[0]=u1; u[1]=u2;} 
};
union X128D_T 
{
    double d[2];
    __m128d d128;
};

#define SSE2_SHUFF 0x1b

#else  /* standard C */
/** 128-bit data structure */
union W128_T 
{
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
    W128_T(uint64_t u1,uint64_t u2) {u[0]=u1; u[1]=u2;};
};
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

#define DSFMT_LOW_MASK  UINT64_C(0x000FFFFFFFFFFFFF)
#define DSFMT_HIGH_CONST UINT64_C(0x3FF0000000000000)
#define DSFMT_SR	12



template<unsigned int MEXP, unsigned int POS1, int SL1, unsigned long long MSK1, unsigned long long MSK2, unsigned long long FIX1, unsigned long long FIX2, unsigned long long PCV1, unsigned long long PCV2>
class dsfmt_generator
{
  public:
    static const unsigned int N = ((MEXP - 128) / 104 + 1);
    static const unsigned int N32 = (N * 4);
    static const unsigned int N64 = (N * 2);
  
#if defined(HAVE_SSE2)
    //  static const union X128I_T sse2_param_mask = {{MSK1, MSK2}};
    //  const union X128I_T sse2_param_mask = {{MSK1, MSK2}};
    const union X128I_T sse2_param_mask;
#else
    const union W128_T sse2_param_mask;
#endif

    dsfmt_generator() :
      sse2_param_mask( MSK1, MSK2 ),
      idx( )
    {
      if(posix_memalign((void**)&dsfmt, 16, sizeof(w128_t) * (N+1))) 
      {
        throw std::bad_alloc();
      }
      // sse2_param_mask.u[0] = MSK1;
      // sse2_param_mask.u[1] = MSK2;
    
      // std::cout << MEXP << " "
      // 	      << POS1 << " "
      // 	      << SL1 << std::endl
      // 	      << std::hex << MSK1 << " "
      // 	      << std::hex << MSK2 << std::endl
      // 	      << std::hex << FIX1 << " "
      // 	      << std::hex << FIX2 << std::endl
      // 	      << std::hex << PCV1 << " "
      // 	      << std::hex << PCV2
      // 	      << std::dec 
      // 	      << std::endl;
    }

    ~dsfmt_generator()
    {
      free((void*)dsfmt);
    }
  
    /**
     * This function initializes the internal state array with a 32-bit
     * integer seed.
     * @param[in] seed a 32-bit integer used as the seed.
     */
    void seed(uint32_t seed) 
    {
      unsigned int i;
      uint32_t *psfmt = (uint32_t *)dsfmt;

      psfmt[idxof(0)] = seed;
      for (i = 1; i < (N + 1) * 4; i++) 
      {
        psfmt[idxof(i)] = 1812433253UL
          * (psfmt[idxof(i - 1)] ^ (psfmt[idxof(i - 1)] >> 30)) + i;
      }
      initial_mask();
      period_certification();
      idx = N64;
    }

    /**
     * This function generates and returns unsigned 32-bit integer.
     * This is slower than SFMT, only for convenience usage.
     * @return double precision floating point pseudorandom number
     */
    uint32_t uint32() 
    {
      uint32_t r;
      uint64_t *psfmt64 = (uint64_t *)dsfmt;

      if (idx >= N64) 
      {
        gen_rand_all();
        idx = 0;
      }
      r = psfmt64[idx++] & 0xffffffffU;
      return r;
    }

    /**
     * This function generates and returns double precision pseudorandom
     * number which distributes uniformly in the range [1, 2).  This is
     * the primitive and faster than generating numbers in other ranges.
     * @return double precision floating point pseudorandom number
     */
    double rand_close1_open2(void) 
    {
      double r;
      double *psfmt64 = (double *)dsfmt;

      if (idx >= N64) 
      {
        gen_rand_all();
        idx = 0;
      }
      r = psfmt64[idx++];
      return r;
    }
  
    /**
     * This function generates and returns double precision pseudorandom
     * number which distributes uniformly in the range [0, 1).
     * @return double precision floating point pseudorandom number
     */
    double rand_close_open(void) 
    {
      return rand_close1_open2() - 1.0;
    }

    /**
     * This function generates and returns double precision pseudorandom
     * number which distributes uniformly in the range (0, 1].
     * @return double precision floating point pseudorandom number
     */
    double rand_open_close(void) 
    {
      return 2.0 - rand_close1_open2();
    }

    /**
     * This function generates and returns double precision pseudorandom
     * number which distributes uniformly in the range (0, 1).
     * @return double precision floating point pseudorandom number
     */
    double rand_open_open(void) 
    {
      double *dsfmt64 = (double *)dsfmt;
      union {
          double d;
          uint64_t u;
      } r;

      if (idx >= N64) 
      {
        gen_rand_all();
        idx = 0;
      }
      r.d = dsfmt64[idx++];
      r.u |= 1;
      return r.d - 1.0;
    }


  private:

    /**
     * This function fills the internal state array with double precision
     * floating point pseudorandom numbers of the IEEE 754 format.
     */
    void gen_rand_all() 
    {
      unsigned int i;
      w128_t lung( dsfmt[N] );

      do_recursion(&dsfmt[0], &dsfmt[0],
                   &dsfmt[POS1], &lung);
      for (i = 1; i < N - POS1; i++) 
      {
        do_recursion(&dsfmt[i], &dsfmt[i],
                     &dsfmt[i + POS1], &lung);
      }
      for (; i < N; i++) 
      {
        do_recursion(&dsfmt[i], &dsfmt[i],
                     &dsfmt[i + POS1 - N], &lung);
      }
      dsfmt[N] = lung;
    }

#if defined(HAVE_SSE2)
    /**
     * This function represents the recursion formula.
     * @param[out] r output 128-bit
     * @param[in] a a 128-bit part of the internal state array
     * @param[in] b a 128-bit part of the internal state array
     * @param[in,out] d a 128-bit part of the internal state array (I/O)
     */
    void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *u) 
    {
      __m128i v, w, x, y, z;
    
      x = a->si;
      z = _mm_slli_epi64(x, SL1);
      y = _mm_shuffle_epi32(u->si, SSE2_SHUFF);
      z = _mm_xor_si128(z, b->si);
      y = _mm_xor_si128(y, z);
    
      v = _mm_srli_epi64(y, DSFMT_SR);
      w = _mm_and_si128(y, sse2_param_mask.i128);
      v = _mm_xor_si128(v, x);
      v = _mm_xor_si128(v, w);
      r->si = v;
      u->si = y;
    }
#else
    /**
     * This function represents the recursion formula.
     * @param[out] r output 128-bit
     * @param[in] a a 128-bit part of the internal state array
     * @param[in] b a 128-bit part of the internal state array
     * @param[in,out] lung a 128-bit part of the internal state array
     */
    void do_recursion(w128_t *r, w128_t *a, w128_t * b, w128_t *lung) 
    {
      uint64_t t0, t1, L0, L1;
    
      t0 = a->u[0];
      t1 = a->u[1];
      L0 = lung->u[0];
      L1 = lung->u[1];
      lung->u[0] = (t0 << SL1) ^ (L1 >> 32) ^ (L1 << 32) ^ b->u[0];
      lung->u[1] = (t1 << SL1) ^ (L0 >> 32) ^ (L0 << 32) ^ b->u[1];
      r->u[0] = (lung->u[0] >> DSFMT_SR) ^ (lung->u[0] & MSK1) ^ t0;
      r->u[1] = (lung->u[1] >> DSFMT_SR) ^ (lung->u[1] & MSK2) ^ t1;
    }
#endif
  
    /**
     * This function initializes the internal state array to fit the IEEE
     * 754 format.
     */
    void initial_mask() 
    {
      unsigned int i;
      uint64_t *psfmt = (uint64_t *)dsfmt;

      for (i = 0; i < N * 2; i++) 
      {
        psfmt[i] = (psfmt[i] & DSFMT_LOW_MASK) | DSFMT_HIGH_CONST;
      }
    }

    /**
     * This function simulate a 32-bit array index overlapped to 64-bit
     * array of LITTLE ENDIAN in BIG ENDIAN machine.
     */
    inline const int idxof(int i) { return i; }

    /**
     * This function certificate the period of 2^{MEXP}-1.
     */
    void period_certification() 
    {
      uint64_t pcv[2] = {PCV1, PCV2};
      uint64_t tmp[2];
      uint64_t inner;
      int i;

      tmp[0] = (dsfmt[N].u[0] ^ FIX1);
      tmp[1] = (dsfmt[N].u[1] ^ FIX2);

      inner = tmp[0] & pcv[0];
      inner ^= tmp[1] & pcv[1];
      for (i = 32; i > 0; i >>= 1) 
      {
        inner ^= inner >> i;
      }
      inner &= 1;
      /* check OK */
      if (inner == 1) 
      {
        return;
      }
      /* check NG, and modification */

      if ( (PCV2 & 1) == 1 )
      {
        dsfmt[N].u[1] ^= 1;
      }
      else
      {
        for (i = 1; i >= 0; i--) 
        {
          uint64_t work = 1;
          for (int j = 0; j < 64; j++) 
          {
            if ((work & pcv[i]) != 0) 
            {
              dsfmt[N].u[i] ^= work;
              return;
            }
            work = work << 1;
          }
        }
      }
      return;
    }


    w128_t *dsfmt; ///< The internal state

    unsigned int idx; ///< counter pointing into internal state
};


typedef dsfmt_generator<19937,117,19,0x000ffafffffffb3fULL,0x000ffdfffc90fffdULL,0x90014964b32f4329ULL,0x3b8d12ac548a7c7aULL,0x3d84e1ac0dc82880ULL,0x0000000000000001ULL> dsfmt_19937_generator;

#endif
