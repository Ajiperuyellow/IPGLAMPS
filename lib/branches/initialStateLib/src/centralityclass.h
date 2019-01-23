//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/centralityclass.h $
//$LastChangedDate: 2014-12-11 00:12:00 +0100 (Do, 11. Dez 2014) $
//$LastChangedRevision: 2014 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief This file contains classes to implement centrality classes
 **/

#ifndef CENTRALITYCLASS_H
#define CENTRALITYCLASS_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** 
 * @brief base class to implement centrality classes/bins
 *
 * Please note: The function of these classes is not really bound to
 * centrality classes. In contrary, it could be used as a basis for
 * complicated binning or histogram routines.
 **/
class tCentralityBase
{
public:
  /** 
   * @brief Constructor
   **/
  tCentralityBase( ) : 
    isUptodate( false ),
    binleft() 
  {};

  /**
   * @brief get the number of classes
   **/
  unsigned int getClassMax(void);

  /**
   * @brief get the class number
   *
   * The return value is in the range 0 ... #getClassMax.
   * The value 0 stands for 'below the first bin'/'not in this
   * binning scheme'/'invalid'
   *
   * You have to adjust your array indexing as e.g.:
   * ~~~
   * unsigned int iclass = getClass( x );
   * if (iclass > 0) bin[iclass-1] += value;
   * ~~~
   *
   * The open/close logic of lower/upper boundary is:
   * \f$
   *   \left[\,b_i, b_{i+1} \right) 
   *   = \left\{ x : b_i \leq x < b_{i+1} \right\}
   * \f$
   * i.e. if the value coincides with the lower boundary, it is
   * included in the corresponding bin.
   **/
  unsigned int getClass(const double val);

  /**
   * @brief Add some new class/bin
   *
   * This routine pushes the given value to the list of lower/left
   * boundaries of the bins.
   **/
  void addLowerBound(const double lbound);

  /**
   * @brief output operator
   **/
  friend std::ostream& operator<<(std::ostream &os, const tCentralityBase &obj);

  /**
   * @brief input operator
   **/
  friend std::istream& operator>>(std::istream &is, tCentralityBase &obj);

  /**
   * @brief debug output
   **/
  virtual void print(std::ostream &os, const std::string & text = "") const;

  /**
   * @brief force the call of the routine #construct
   **/
  void update(void) { construct(); };

  /**
   * @brief reset lowest boundary
   **/
  void resetLowestBoundary( const double lbound );

protected:
  /**
   * @brief Routine to calculate the bins
   **/
  virtual void construct(void) 
  {
    isUptodate = true;
  };
  

protected:

  bool isUptodate; ///< indicates whether #construct has to be called

  std::vector<double> binleft; ///< the left border of the bins

};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// now the real fun pops up: we define different methods to define the
// bins:

/** 
 * @brief class to implement centrality classes/bins by giving the
 * weights and a list of values
 *
 * In this class, we continously add the values, which will be
 * internally histogramed, and the final bin boundaries are
 * calculated according pre-given integral ratios.
 *
 * We start by generating the list of weights, e.g.
 * ~~~
 * addWeight( 0.5 );
 * addWeight( 0.3 );
 * ~~~
 * This will prepare the calculations to a setup, where the first 50%
 * make up the first bin, the next 30% the second bin and the final
 * 20% of entries the final bin.
 * 
 * The list of values will be searched for finding the given
 * boundaries. Internally, the list of values will be sorted, and
 * then, in the given example, the lowest boundary will be the value
 * at position given by half of the size of the vector. 
 * Correspondingly, the lower boundary of bin 3 will be the 
 * value of the entry at position 80% of the vector's size (and reach
 * till infinity).
 *
 * If the weight of some bin is too small, then the left boundaries of
 * that class and the one of the next class will be identical, thus
 * making no harm; this bin will just be ignored.
 *
 **/
class tCentrality_Values : public tCentralityBase
{
public:
  /** 
   * @brief Constructor
   **/
  tCentrality_Values( ) : 
    tCentralityBase(),
    sumWeights(),
    weights(),
    values()
  {};

  /**
   * @brief add the weight of a new bin
   **/
  void addWeight(const double w);

  /**
   * @brief add the weight of a new bin
   **/
  void addValue(const double val);

protected:
  /**
   * @brief Routine to calculate the bins
   **/
  virtual void construct(void);

protected:
  /**
   * @brief for debugging: print weight  list
   **/
  void printWeights(void);

protected:

  double sumWeights; ///< the sum of all weights
  std::vector<double> weights; ///< the list of weights

  std::vector<double> values; ///< the list of values
};


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** 
 * @brief exception class for handling unexpected critical behaviour within the centrality class  
 * */
class eCentrality_error : public std::runtime_error
{
public:
  explicit eCentrality_error ( const std::string& what ) : std::runtime_error ( what ) {};
  
  virtual ~eCentrality_error() throw() {};
};


#endif
