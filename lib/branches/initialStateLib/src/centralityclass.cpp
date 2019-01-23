//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/centralityclass.cpp $
//$LastChangedDate: 2014-12-11 22:00:00 +0100 (Do, 11. Dez 2014) $
//$LastChangedRevision: 2015 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "centralityclass.h"

#include <algorithm>

unsigned int tCentralityBase::getClassMax(void)
{
  if (!isUptodate) construct();
  return binleft.size(); 
}

unsigned int tCentralityBase::getClass(const double val)
{
  if (!isUptodate) construct();

  if (binleft.size() == 0) return 0;

  unsigned int nBin = 0;
  for (auto it=begin(binleft); it!=end(binleft); ++it)
  {
    if (val < *it) break;
    nBin++;
  }
  return nBin;

}

/**
 * After a call to this routine, the flag #isUptodate is true.
 **/
void tCentralityBase::addLowerBound(const double lbound)
{
  if (!isUptodate) construct();

  if (binleft.size() > 0)
  {
    if (lbound <= binleft.back())
    {
      throw eCentrality_error("boundaries not increasing!");
    }
  }
  binleft.push_back( lbound );
}

std::ostream& operator<<(std::ostream &os, const tCentralityBase &obj)
{
  if (!obj.isUptodate)
  {
    throw eCentrality_error("object not uptodate when trying to print.");
  }
  
  os << obj.binleft.size() << "  ";
  for (auto it=begin(obj.binleft); it!=end(obj.binleft); ++it)
  {
    os << *it << "  ";
  }
  return os;
}

std::istream& operator>>(std::istream &is, tCentralityBase &obj)
{
  obj.binleft.clear();
  
  unsigned int n;
  is >> n;
  
  for( unsigned int i=0; i<n; i++)
  {
    double x;
    is >> x;
    
    obj.binleft.push_back( x );
  }
  
  obj.isUptodate = true;
  
  return is;
}

void tCentralityBase::print(std::ostream &os, const std::string & text) const
{
  if (!isUptodate)
  {
    throw eCentrality_error("object not uptodate when trying to print.");
  }
  
  unsigned int n = binleft.size();
  os << "=============================" << std::endl
     << "# " << text << std::endl
     << "# number of bins = " << n << std::endl
     << "# bin = left boundary ... right boundary :" << std::endl;

  if (n>0)
  {

    for( unsigned int i=0; i<n-1; i++)
    {
      os << "Bin " << i << " = " 
         << binleft[i] << " ... "
         << binleft[i+1] << " : " 
         << binleft[i+1]-binleft[i]
         << std::endl;
    }
    os << "Bin " << n-1 << " = " 
       << binleft[n-1] << " ... "
       << "inf" << " : " 
       << std::endl; 
  }
  
}

void tCentralityBase::resetLowestBoundary( const double lbound )
{
  if (!isUptodate) construct();
  if (binleft.size() > 0)
  {
    binleft[0] = lbound;
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void tCentrality_Values::addWeight(const double w)
{
  // remove last bin, because it was added automatically.
  if (isUptodate) 
  {
    weights.pop_back();
  }

  sumWeights += w;
  if (sumWeights > 1.0)
  {
    std::cout << "error while adding weight " << w << std::endl;
    throw eCentrality_error("sum of weights > 1");
  }
  weights.push_back( w );
  isUptodate = false;
}

void tCentrality_Values::addValue(const double val)
{
  values.push_back( val );
  isUptodate = false;
}

/**
 * It is not necessary, that the given weights sum up to 1. (At least,
 * they must not overshoot 1!) The last bin is given by 
 * ''boundary ... infinity'', 
 * which automatically takes care of the remaining bin weights.
 **/
void tCentrality_Values::construct(void)
{
  const bool verbose = false;

  if (isUptodate) return; // ministry of redundant redundancy

  isUptodate = false; // default return value: failure
  
  if (values.size() == 0) return; // failure!!!

  // okay, now we are doing the real work:

  //  printWeights();

  binleft.clear(); // remove all previous information

  std::sort(values.begin(),values.end()); // sort the values !!!

  binleft.push_back( values[0] ); // add the lowest boundary

  double sWeight = 0.0;
  for (auto it=begin(weights); it!=end(weights); ++it)
  {
    sWeight += (*it);
    // this is the first guess:
    unsigned int iZero = sWeight * values.size() - 1;


    // look for different values:

    unsigned int iMinus = iZero;
    for ( ; iMinus+1 > 0; iMinus--) // 'unsigned int' down-loop!
    {
      if (values[iMinus] != values[iZero]) 
      {
        iMinus++;
        break;
      }
    }
    unsigned int iPlus = iZero;
    for ( ; iPlus<values.size(); iPlus++)
    {
      if (values[iPlus] != values[iZero]) 
      {
        break;
      }
    }
    
    // now all values[iMinus...iZero...iPlus-1] have the same value.
    // we select the i-value closest to iZero: 

    unsigned int iDo = (iPlus-iZero)<(iZero-iMinus)? iPlus : iMinus;

    if (verbose)
    {
      std::cout << "doing cut at entry = " 
                << iDo << " "
                << iZero << " "
                << iMinus << " "
                << iPlus << " "
                << values.size() << " "
                << values[iZero] << " "
                << (*it)  << " "
                << std::endl;
    }
    binleft.push_back( values[iDo] );

  }
  
  isUptodate = true;

  
  // unsigned int i = 0;
  // for( auto it = values.begin(); it != values.end(); ++it,++i )
  // {
  //   std::cout << i << "  " 
  //             << getClass( *it ) << "  "
  //             << *it << std::endl;
  // }
  
}


void tCentrality_Values::printWeights(void)
{
  std::cout <<  __PRETTY_FUNCTION__ << std::endl
            << weights.size() << std::endl;
  for(unsigned int i=0; i<weights.size(); i++)
  {
    std::cout << weights[i] << std::endl;
  }
  std::cout << "sum = " << sumWeights << std::endl << std::endl;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
