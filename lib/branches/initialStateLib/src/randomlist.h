//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/randomlist.h $
//$LastChangedDate: 2014-07-04 15:09:30 +0200 (Fr, 04. Jul 2014) $
//$LastChangedRevision: 1782 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef RANDOMLIST_H
#define RANDOMLIST_H

#include <cstdlib>
#include <iostream>
#include <vector>

/** 
 * @brief This class implements a randomized list
 *
 * This class shuffles the numbers 1...N into random order.
 *
 * In addition, it allows for accessing slices. This may be used,
 * if you e.g. want to seperate the particle vector into some number M
 * of sublists, perhaps with M=number of tesparticles.
 *
 * It uses internally the Fisher–Yates shuffle, as implemented by 
 * Durstenfeld, as an in-place shuffle.
 *
 * Example
 * \code{.cpp}
    randomlist<int> rList(10,3);
    std::cout << "# slices: " << rList.numberSlice() << endl;
    for (unsigned int iSlice=0; iSlice < rList.numberSlice(); iSlice++)
    {
      std::cout << "#" << iSlice << " : ";
      for (auto it=rList.beginSlice(iSlice); it!=rList.endSlice(iSlice); ++it) 
      {
        std::cout << (*it) << "  ";
      }
      std::cout << std::endl;
    }
\endcode
 * will produce for example follwing output
 * \code
# slices: 3
   #0 : 9  1  5  6  
   #1 : 2  0  7  4  
   #2 : 8  3 
\endcode
 **/

template<class T>
class randomlist
{
  protected:
    std::vector<T> TheList;
    size_t nSlice; ///< number of the slices
    size_t sSlice; ///< size of the slices

  public:
    /**
     * @brief constructor
     *
     * The list is set to the given size and initialiazed in random
     * order via the Durstenfeld in-place shuffle.
     **/
    randomlist( size_t size, size_t numberSlice = 1 ) :
      TheList(),
      nSlice( numberSlice )
    {
      sSlice = (size-1)/nSlice + 1;


      TheList.reserve( size );

      for (unsigned int i = 0; i < size; i++)
      {
        unsigned int j = static_cast<int>( (i+1)*ran2() );
        if ( j == i )
        {
          TheList.push_back( i );
        }
        else
        {
          TheList.push_back( TheList[j] );
          TheList[j] = i;
        }
      }
    };

    /** @brief return the slice size (setter routine) */
    size_t &numberSlice(void) { return nSlice; };
    
    /** @brief return the slice size (getter routine) */
    size_t numberSlice(void) const { return nSlice; };

    /** @brief return the number of slices */
    size_t sizeSlice(void) const { return sSlice; }
    
    
    /** @brief Defining an iterator **/
    typedef typename std::vector<T>::iterator iterator;
    
    /** @brief Defining an iterator (const version) **/
    typedef typename std::vector<T>::const_iterator const_iterator;

    /** @brief return begin iterator */
    inline iterator begin(void) { return TheList.begin(); }

    /** @brief return end iterator */
    inline iterator end(void) { return TheList.end(); }
    
    /** @brief return begin iterator (const) */
    inline const_iterator begin(void) const { return TheList.begin(); }

    /** @brief return end iterator (const) */
    inline const_iterator end(void) const { return TheList.end(); }

    /** @brief return begin iterator of slice i */
    inline iterator beginSlice(size_t i) 
    { 
      return i*sSlice < TheList.size() ? TheList.begin() + i*sSlice : TheList.end(); 
    }

    /** @brief return end iterator of slice i */
    inline iterator endSlice(size_t i) 
    { 
      return (i+1)*sSlice < TheList.size() ? TheList.begin() + (i+1)*sSlice : TheList.end(); 
    }

    /** @brief return begin iterator of slice i (const)*/
    inline const_iterator beginSlice(size_t i) const 
    { 
      return i*sSlice < TheList.size() ? TheList.begin() + i*sSlice : TheList.end(); 
    }

    /** @brief return end iterator of slice i (const)*/
    inline const_iterator endSlice(size_t i) const 
    { 
      return (i+1)*sSlice < TheList.size() ? TheList.begin() + (i+1)*sSlice : TheList.end(); 
    }

};

#endif
