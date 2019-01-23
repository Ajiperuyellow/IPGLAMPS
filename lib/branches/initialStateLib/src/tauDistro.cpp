//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/tauDistro.cpp $
//$LastChangedDate: 2016-03-23 21:04:01 +0100 (Mi, 23. Mär 2016) $
//$LastChangedRevision: 2319 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "tauDistro.h"

#include <fstream>
#include <iomanip>      // std::setw
#include <sstream>

tTauDistro::tTauDistro(const std::vector<double> & tstep,
                       const std::string & fName0,
                       const double vMax, const double vDelta) :
  taustep( tstep ),
  fNamePre( fName0 ),
  toPrint( true )
{
  int nStep = tstep.size();
  int nBins = 2*vMax/vDelta;

  dNdeta.resize(nStep, binningValues("",-vMax,vMax, nBins) );
  dEdeta.resize(nStep, binningValues("",-vMax,vMax, nBins) );
  dNdy.resize(nStep, binningValues("",-vMax,vMax, nBins) );
  dEdy.resize(nStep, binningValues("",-vMax,vMax, nBins) );
}

void tTauDistro::print(void)
{
  if (!toPrint) return;
  toPrint = false;

  for (int iStep=0; iStep<taustep.size(); ++iStep)
  {
    const std::string fName = filename_step(iStep) + "_tauDistro.rep";
    std::fstream file(fName.c_str(),std::ios::out);

    file << "# dN/deta_s:" << std::endl;
    dNdeta[iStep].print(file);

    file << std::endl << std::endl << "# dE_T/deta_s:" << std::endl;
    dEdeta[iStep].print(file);

    file << std::endl << std::endl << "# dN/dy:" << std::endl;
    dNdy[iStep].print(file);

    file << std::endl << std::endl << "# dE_T/dy:" << std::endl;
    dEdy[iStep].print(file);

    file.close();
  }

}

std::string tTauDistro::filename_step( const int step ) const
{
  std::string name;

  if ( step == 0 )
    name = "initial";
  else if ( step == taustep.size()+1 )
    name = "final";
  else
  {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(3) << step;
    name = "step" + ss.str();
  }
  return fNamePre + "_" + name;
}
