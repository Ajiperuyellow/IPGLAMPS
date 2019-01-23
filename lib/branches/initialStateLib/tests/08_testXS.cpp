//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/08_testXS.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "hadroncrosssection.h"

#include <iostream>

using namespace std;

void test1(void)
{
  //  if (!FPT_COMP_E(Res, -28.)) throw 101;
  XSgenerator_pp_tot_PDG sigma_pp_tot;
  cout << sigma_pp_tot(30.0) << endl;
  cout << sigma_pp_tot(200.0) << endl;
  cout << sigma_pp_tot(2760.0) << endl;

  cout << endl;
}

void test2(void)
{
  //  if (!FPT_COMP_E(Res, -28.)) throw 101;
  XSgenerator_pp_tot_PDG sigma_tot0;
  XSgenerator_pp_tot_Pythia sigma_tot;
  XSgenerator_pp_inelast_Pythia sigma_inelast;

  double sqrts = 30.0;

  cout << sqrts << " "
       << sigma_tot0(sqrts) << " " 
       << sigma_tot(sqrts) << " " 
       << sigma_inelast(sqrts) << " " 
       << endl;

  sqrts = 200.0;

  cout << sqrts << " "
       << sigma_tot0(sqrts) << " " 
       << sigma_tot(sqrts) << " " 
       << sigma_inelast(sqrts) << " " 
       << endl;
  
  sqrts = 2760.0;

  cout << sqrts << " "
       << sigma_tot0(sqrts) << " " 
       << sigma_tot(sqrts) << " " 
       << sigma_inelast(sqrts) << " " 
       << endl;

}

int main(int argc, char **argv)
{

  try
  {
    test1();
    test2();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
