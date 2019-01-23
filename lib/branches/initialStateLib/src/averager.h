//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/averager.h $
//$LastChangedDate: 2014-11-12 14:31:26 +0100 (Mi, 12. Nov 2014) $
//$LastChangedRevision: 1924 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declaration of calculator of averages
 *
 */

#ifndef AVERAGER_H
#define AVERAGER_H



/**
 * @brief class to implement a templated statistical 'averager'
 *
 * The class T has to provide:
 * - T += T
 * - T *= T
 * - T = T + T
 * - T = T * T
 * - T *= a
 *
 * Instead of T/a we use T * (1.0/a) -- due to speed
 *
 * Instead of T-T we use T+(T*(-1)) -- due to minimal writing effort
 * (you do not have to implement -= and - operator). Since we always
 * use it as T-T*a, thus T+(T*(-a)), this is not really a blocker.
 *
 **/
template<typename T>
class tMultiAverager
{  
protected:
  long double _count; ///< store number of values
  T _sum;   ///< store sum of values
  T _sum2;  ///< store sum of values squared
  
public:
  /** @brief Default constructor */
  tMultiAverager() :  _count(0), _sum(), _sum2() {};
  
  /** @brief reset all values */
  void reset(void)
  {
    _count = 0;
    _sum = T();
    _sum2 = T();
  }
  
  /** @brief add some value */
  void add(const T val) 
  {
    ++_count;
    _sum += val;
    _sum2 += val*val;
  }
  
  /** @brief add some value */
  void add(const T val, long double weight) 
  {
    _count += weight;
    _sum += val*weight;
    _sum2 += val*val*weight;
  }
  
  /** @brief return the number of entries */
  long double Count() const
  {
    return _count;
  }
  
  /** @brief return the mean value */
  T Mean() const
  {
    return (_count==0.0 ? T() : _sum * (1.0/_count));
  }
  
  /** @brief return the mean value squared */
  T Mean2()
  {
    return (_count==0.0 ? T() : _sum*_sum * (1.0/(_count*_count)));
  }
  
  /** @brief return the variance */
  T Var()
  {
    return (_count<=1.0 ? T() : _sum2 * (1.0/(_count-1)) + (_sum*_sum) * (-1.0/(_count*(_count-1))));
  }
  
  /** @brief return the variance in the large number limit */
  T VarA()
  {
    return (_count==0.0 ? T() : _sum2 * (1.0/_count) + (_sum*_sum) * (-1.0/(_count*_count)));
  }
  
  /** 
   * @brief The Increment Operator 
   *
   * every component is added
   **/
  inline tMultiAverager & operator += (const tMultiAverager & a)
  {
    _count += a._count;
    _sum += a._sum;
    _sum2 += a._sum2;
    return *this;
  }
  
};



/**
 * @brief for convenience: this is the simplest case
 **/
typedef tMultiAverager<long double> tAverager;

#pragma omp declare                            \
  reduction(+ : tAverager :                    \
            omp_out += omp_in )                \
  initializer( omp_priv = tAverager() )


#endif
