//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/global.h $
//$LastChangedDate: 2014-02-12 20:15:39 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1620 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef GLOBAL_H
#define GLOBAL_H

/** @file
 * @brief  This file defines some global functions in the namespace nc_casc
 **/

#include <iostream>
#include <ostream>

namespace ns_casc
{
  /**
   * @brief Writing out some short SVN revision message
   *
   * actually, it just writes the content of SVN_REVISION.
   *
   * This is made as some own compiled function in order to release
   * dependencies on 'revision.h'
   **/
  void printSVNshort(std::ostream & f);

  /**
   * @brief Writing out some long SVN revision message
   **/
  void printSVN(std::ostream & f);

  
}




#endif
