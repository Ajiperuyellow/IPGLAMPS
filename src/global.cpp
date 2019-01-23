//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/global.cpp $
//$LastChangedDate: 2014-02-12 20:15:39 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1620 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "global.h"
#include "revision.h"


void ns_casc::printSVNshort(std::ostream & f)
{
  f << SVN_REVISION;
}

void ns_casc::printSVN(std::ostream & f)
{
  f << "# Version: SVN revision " << SVN_REVISION << std::endl
    << "# Version: SVN status -q (SRC)" << std::endl 
    << SVN_STATUS_SRC << std::endl
    << "# Version: SVN status -q (LIB)" << std::endl 
    << SVN_STATUS_LIB << std::endl
    << "#" << std::endl << std::endl;

}
