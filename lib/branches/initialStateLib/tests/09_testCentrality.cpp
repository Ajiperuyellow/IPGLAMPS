//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/09_testCentrality.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "centralityclass.h"


#include <iostream>

using namespace std;

void test1(void)
{
  tCentralityBase CClass;

  if (CClass.getClassMax() != 0) throw 101;
  if (CClass.getClass(99.9) != 0)throw 102;

  CClass.addLowerBound( 0.0 );
  if (CClass.getClass(-99.9) != 0) throw 103;
  if (CClass.getClass(99.9) != 1) throw 104;

  CClass.addLowerBound( 1.0 );
  CClass.addLowerBound( 2.0 );
  CClass.addLowerBound( 3.0 );

  if (CClass.getClassMax() != 4) throw 105;
  if (CClass.getClass(2.7) != 3) throw 106;
  if (CClass.getClass(3.7) != 4) throw 107;

  cout << CClass << endl;

  CClass.print( cout, "test1" );

}

void test2(void)
{
  tCentrality_Values CClass;

  if (CClass.getClassMax() != 0) throw 201;
  if (CClass.getClass(99.9) != 0)throw 202;

  CClass.addWeight( 0.5 ); // defining bin 1
  CClass.addWeight( 0.3 ); // defining bin 2
                           // bin 3 defined automatically

  CClass.addWeight( 1e-4 );



  for (int i=0; i<10; i++)
  {
    for (int j=0; j<10; j++)
    {
      CClass.addValue( 3.0*i + 2.0*j + 0.1 );
    }
  }
  CClass.update();

  CClass.print( cout, "test2" );

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
