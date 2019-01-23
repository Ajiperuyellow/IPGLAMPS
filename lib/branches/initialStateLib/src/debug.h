//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/debug.h $
//$LastChangedDate: 2014-12-03 22:40:43 +0100 (Mi, 03. Dez 2014) $
//$LastChangedRevision: 2002 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef DEBUG_H
#define DEBUG_H

#if _OPENMP
#warning "using OpenMP"
#include <omp.h>
#else
#warning "no OpenMP"
#endif


/**
 * @brief This is a poor mans way of debugging: print out line numbers
 **/
//#define DDD {cout << "DDD: " << __FILE__ << " at " << __LINE__ << endl;};
//#define DDD {cout << "DDD: " << __FILE__ << " at " << __LINE__ << "  " << omp_get_thread_num() << endl;};



#ifndef DDD
#define DDD
#endif

#endif
