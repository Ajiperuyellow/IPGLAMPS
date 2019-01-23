//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/hydroTensor.cpp $
//$LastChangedDate: 2016-05-05 14:20:29 +0200 (Do, 05. Mai 2016) $
//$LastChangedRevision: 2339 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <cmath>
#include <complex>
#include <iostream>
#include <string>

#include "hydroTensor.h"
#include "FPT_compare.h"
//#include "tools.h"

using namespace std;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void tTmunuNmu_Base::add(const VectorEPxPyPz & Mom, const double fak)
{
  double fakE = fak / Mom.E();

  _T00 += Mom.E() * fak;
  _T10 += Mom.Px() * fak;
  _T20 += Mom.Py() * fak;
  _T30 += Mom.Pz() * fak;

  _T11 += Mom.Px2() * fakE;
  _T22 += Mom.Py2() * fakE;
  _T33 += Mom.Pz2() * fakE;
  
  _T21 += Mom.Py()*Mom.Px() *  fakE;
  _T31 += Mom.Pz()*Mom.Px() *  fakE;
  _T32 += Mom.Pz()*Mom.Py() *  fakE;
  
  _N0 += fak;
  _N1 += Mom.Px() * fakE;
  _N2 += Mom.Py() * fakE;
  _N3 += Mom.Pz() * fakE;
  
  _NN++;
};

std::string tTmunuNmu_Base::header(int & nr)
{
  std::string text = "";
  text += std::to_string( nr++ ) + ": T00 ";
  text += std::to_string( nr++ ) + ": T11 ";
  text += std::to_string( nr++ ) + ": T22 ";
  text += std::to_string( nr++ ) + ": T33 ";
  text += std::to_string( nr++ ) + ": T10 ";
  text += std::to_string( nr++ ) + ": T20 ";
  text += std::to_string( nr++ ) + ": T30 ";
  text += std::to_string( nr++ ) + ": T21 ";
  text += std::to_string( nr++ ) + ": T31 ";
  text += std::to_string( nr++ ) + ": T32 ";
  text += std::to_string( nr++ ) + ": N0 ";
  text += std::to_string( nr++ ) + ": N1 ";
  text += std::to_string( nr++ ) + ": N2 ";
  text += std::to_string( nr++ ) + ": N3 ";
  text += std::to_string( nr++ ) + ": NN ";
  return text;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

std::string tPimunu::header(int & nr)
{
  std::string text = "";
  text += std::to_string( nr++ ) + ": Pi00 ";
  text += std::to_string( nr++ ) + ": Pi11 ";
  text += std::to_string( nr++ ) + ": Pi22 ";
  text += std::to_string( nr++ ) + ": Pi33 ";
  text += std::to_string( nr++ ) + ": Pi10 ";
  text += std::to_string( nr++ ) + ": Pi20 ";
  text += std::to_string( nr++ ) + ": Pi30 ";
  text += std::to_string( nr++ ) + ": Pi21 ";
  text += std::to_string( nr++ ) + ": Pi31 ";
  text += std::to_string( nr++ ) + ": Pi32 ";
  return text;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


tdeltamunu::tdeltamunu( const VectorTXYZ & v )
{
  double gamma2 = 1.0/(fabs(1.0 - v.vec2()));

  D_0_0 = D$0_0 = D_0$0 = D$0$0 = 1.0 - gamma2;

  D_1_1 = D$1$1 = -1.0 - gamma2 * v.X2();
  D_2_2 = D$2$2 = -1.0 - gamma2 * v.Y2();
  D_3_3 = D$3$3 = -1.0 - gamma2 * v.Z2();

  D$1_1 = 1.0 + gamma2 * v.X2();
  D$2_2 = 1.0 + gamma2 * v.Y2();
  D$3_3 = 1.0 + gamma2 * v.Z2();

  D$1_0 = D$1$0 = D$0$1 = - gamma2 * v.X();
  D$2_0 = D$2$0 = D$0$2 = - gamma2 * v.Y();
  D$3_0 = D$3$0 = D$0$3 = - gamma2 * v.Z();

  D$0_1 = D_1_0 = D_0_1 = gamma2 * v.X();
  D$0_2 = D_2_0 = D_0_2 = gamma2 * v.Y();
  D$0_3 = D_3_0 = D_0_3 = gamma2 * v.Z();

  D_2_1 = D$2$1 = D_1_2 = D$1$2 = - gamma2 * v.X() * v.Y();
  D_3_1 = D$3$1 = D_1_3 = D$1$3 = - gamma2 * v.X() * v.Z();
  D_3_2 = D$3$2 = D_2_3 = D$2$3 = - gamma2 * v.Y() * v.Z();

  D$2_1 = D$1_2 = gamma2 * v.X() * v.Y();
  D$3_1 = D$1_3 = gamma2 * v.X() * v.Z();
  D$3_2 = D$2_3 = gamma2 * v.Y() * v.Z();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

VectorTXYZ tTmunuNmu::Eckart_velocity( ) const
{
  return vector4D(N0(),N1(),N2(),N3()).NormalizeToE();
}

VectorTXYZ tTmunuNmu::Landau_velocity( const double energyDensity ) const
{
  const double x11 = T11()+energyDensity;
  const double x22 = T22()+energyDensity;
  const double x33 = T33()+energyDensity;

  return vector4D( x11  *(x22  *x33 - T32()*T32()) - T21()*(T21()*x33   - T32()*T31()) + T31()*(T21()*T32() - T31()*x22  ),
                   T10()*(x22  *x33 - T32()*T32()) - T21()*(T20()*x33   - T32()*T30()) + T31()*(T20()*T32() - T30()*x22  ),
                   x11  *(T20()*x33 - T32()*T30()) - T10()*(T21()*x33   - T32()*T31()) + T31()*(T21()*T30() - T20()*T31()),
                   x11  *(T30()*x22 - T20()*T32()) - T21()*(T21()*T30() - T20()*T31()) + T10()*(T21()*T32() - T31()*x22  ) ).NormalizeToE();
}

/**
 * instead of making the transition v -> u = v*gamma, we use v and
 * multiply the whole result with gamma^2
 *
 * please note: we have changed the signs of the 0i components, since
 * v(i)*gamma corresponds to u$i, and not u_i, as in the original
 * formula
 **/
double tTmunuNmu::Eckart_energyDensity( const VectorTXYZ & v ) const
{
  const double gamma2 = 1.0/(fabs(1.0 - v.vec2()));

  return
    (   v(0)*T00()*v(0) - v(0)*T01()*v(1) - v(0)*T02()*v(2) - v(0)*T03()*v(3)
      - v(1)*T10()*v(0) + v(1)*T11()*v(1) + v(1)*T12()*v(2) + v(1)*T13()*v(3)
      - v(2)*T20()*v(0) + v(2)*T21()*v(1) + v(2)*T22()*v(2) + v(2)*T23()*v(3)
      - v(3)*T30()*v(0) + v(3)*T31()*v(1) + v(3)*T32()*v(2) + v(3)*T33()*v(3) )
    * gamma2;
}

double tTmunuNmu::Landau_energyDensity( bool & solution ) const
{
  solution = false; // set default value
  double energyDensity = 0.0; // set default value

  //-----------------------------------------
  //energy density in landau frame
  //this is valid for a common equation of 4. degree
  //still in beta-phase!!
  //-----------------------------------------
    
  std::complex<double> y1(0.0, 0.0 );
  std::complex<double> y2(0.0, 0.0 );
  std::complex<double> y3(0.0, 0.0 );
  std::complex<double> y4(0.0, 0.0 );
  std::complex<double> e1(0.0, 0.0 );
  std::complex<double> e2(0.0, 0.0 );
  std::complex<double> e3(0.0, 0.0 );
  std::complex<double> e4(0.0, 0.0 ); 
  std::complex<double> z1(0.0, 0.0 );
  std::complex<double> z2(0.0, 0.0 );
  std::complex<double> z3(0.0, 0.0 );
  std::complex<double> sqrt_z1(0.0, 0.0 );
  std::complex<double> sqrt_z2(0.0, 0.0 );
  std::complex<double> sqrt_z3(0.0, 0.0 );
    
  //-----------------------------------------
  //Calculating the energydensity with analytical method
  //-----------------------------------------
  //it has to be of 4th degree, otherwise the solution is not correct. Be careful!!

  double c_A = -1.0;
  double c_B = T00() - (T11() + T22() + T33());
    
  double c_C = (T00()*(T11()+T22()+T33()) - T11()*T22() - T11()*T33() - T22()*T33() 
         + pow(T32(),2) + pow(T21(),2) + pow(T31(),2) 
         - pow(T10(),2) - pow(T20(),2) - pow(T30(),2));

  double c_D = (T00()*(T11()*T22() + T11()*T33() + T22()*T33() 
                  - pow(T32(),2) - pow(T21(),2) - pow(T31(),2))
	 - T11()*T22()*T33() - 2.0*T21()*T32()*T31() 
         + T11()*pow(T32(),2) + T33()*pow(T21(),2) + T22()*pow(T31(),2)
	 - pow(T10(),2)*(T22() + T33()) + 2.0*T21()*T10()*T20() + 2.0*T31()*T10()*T30()
	 - pow(T20(),2)*(T11() + T33()) + 2.0*T20()*T32()*T30()
	 - pow(T30(),2)*(T11() + T22()));

  double c_E = (T00()*T11()*T22()*T33() + 2.0*T00()*T21()*T32()*T31()
         + pow(T10(),2)*pow(T32(),2) + pow(T31(),2)*pow(T20(),2) + pow(T21(),2)*pow(T30(),2)
	 + 2.0*T21()*T10()*T20()*T33() + 2.0*T20()*T11()*T32()*T30() 
         + 2*T31()*T10()*T30()*T22()
	 - 2.0*T21()*T10()*T32()*T30() - 2.0*T31()*T10()*T20()*T32() 
         - 2.0*T31()*T21()*T30()*T20()
	 - T00()*T11()*pow(T32(),2) - T00()*T33()*pow(T21(),2) 
         - T00()*T22()*pow(T31(),2) - T22()*T33()*pow(T10(),2)
         - T11()*T33()*pow(T20(),2) - T11()*T22()*pow(T30(),2));

  //switch to general form
  //x^4 + d_B*y^3 + d_C*y^2 + d_D*y + d_E = 0
  double d_B = c_B/c_A;
  double d_C = c_C/c_A;
  double d_D = c_D/c_A;
  double d_E = c_E/c_A;
    
  //substitution: x = y -d_A; still 4 degree, but without y^3
  //y^4 + e_C*y^2 + e_D*y + e_E = 0
  double e_C = -3.0/8.0*pow( d_B,2 ) + d_C;
  double e_D = 1.0/8.0*pow( d_B,3 ) - d_B*d_C/2.0 + d_D;
  double e_E = -3.0*pow( d_B,4 )/256.0 + pow( d_B,2 )*d_C/16.0 - d_B*d_D/4.0 + d_E;
    
  //transfrom to a cubic equation;
  //z^3 + e_R*z^2 + e_S*z + e_T = 0
  double e_R = -2.0 * e_C;
  double e_S = pow( e_C,2 ) - 4.0*e_E;
  double e_T = pow( e_D,2 );

  //-----------------------------------------
  //the solutions of the cubic equation
  //-----------------------------------------
  {
    double p = e_S - pow( e_R,2 ) / 3.0;
    double q = 2.0/27.0*pow( e_R,3 ) - e_R*e_S/3.0 + e_T;
    double R = pow(q/2.0,2) + pow(p/3.0,3);
      
    if(R < 0.0)
    {
      double a = q/2.0;
      double b = - sqrt(-R);
      double r = fabs(sqrt( pow(a,2) - R));
      double phi = atan2(b,a);
      if (phi <= 0.0) phi += M_PI;
      z1 = + 2.0*pow(r,1.0/3.0)*cos(phi/3.0);
      z2 = - 2.0*pow(r,1.0/3.0)*cos(phi/3.0 - M_PI/3.0);
      z3 = - 2.0*pow(r,1.0/3.0)*cos(phi/3.0 + M_PI/3.0);   
    }
    else
    {
      double v,u,u$3,v$3;
                        
      u$3 = -q/2.0 + sqrt(R);
      v$3 = -q/2.0 - sqrt(R);
                        
      if(v$3 < 0.0)
      {
	v$3 = -v$3;
	v = -pow(v$3,1.0/3.0);
      }
      else
      {
	v = pow(v$3,1.0/3.0);
      }

      if(u$3 < 0.0)
      {
	u$3 = -u$3;
	u = -pow(u$3,1.0/3.0);
      }
      else
      {
	u = pow(u$3,1.0/3.0);
      }
            
      std::complex<double> z2_temp(-(u+v)/2.0,-(u-v)/2.0*sqrt(3.0));
      std::complex<double> z3_temp(-(u+v)/2.0, (u-v)/2.0*sqrt(3.0));   
          
      z1 = u + v;
      z2 = z2_temp;
      z3 = z3_temp;
    }
  }

  //Substituing back with z = z - e_R/3
  z1 -= e_R/3.0;
  z2 -= e_R/3.0;
  z3 -= e_R/3.0;

    
  //define the square roots, they have to be negativ except
  //for sqrt(z2), where sign depends on e_D
  sqrt_z1 = - sqrt( -z1 );
  sqrt_z2 = ( -e_D >= 0.0 ) ? sqrt( -z2 ) : - sqrt( -z2 );
  sqrt_z3 = - sqrt( -z3 );  

  //find solution for y by adding all z
  y1 = (  sqrt_z1 + sqrt_z2 + sqrt_z3 )/2.0; 
  y2 = (  sqrt_z1 - sqrt_z2 - sqrt_z3 )/2.0;
  y3 = ( -sqrt_z1 + sqrt_z2 - sqrt_z3 )/2.0; 
  y4 = ( -sqrt_z1 - sqrt_z2 + sqrt_z3 )/2.0;    
    
  //substituing back x = y-d_B/4
  e1 = y1 - d_B/4.0;
  e2 = y2 - d_B/4.0;
  e3 = y3 - d_B/4.0;
  e4 = y4 - d_B/4.0;
    
  //in general there is one positive solution, which represents the energy density
  //and three negative one, which are the pressure of the energy momentum tensor;
  //having only 2 particles, there is a possibility for two positive solution, so
  //the largest one has to be chosen;
  //in the case there is only one particle, the solution is definitve not exact!!
  //then one should set it to T00, as approximation!!
    
  int nn_e = 0;
  if(e1.real()> 0.0){energyDensity = e1.real(); nn_e++;}
  if(e2.real()> 0.0){energyDensity = e2.real(); nn_e++;}
  if(e3.real()> 0.0){energyDensity = e3.real(); nn_e++;}
  if(e4.real()> 0.0){energyDensity = e4.real(); nn_e++;}
    

  switch (nn_e)
  {
    case 0:
    {
      //cout << "No positive solution in eDensity Landau frame! " << endl;
      energyDensity = T00();
      solution = false;
    } break;
    case 1:
    {
      solution = true;
    } break;
    case 2:
    {
      //cout << "2 positive solution in eDensity Landau frame! " << endl;      
      if( e1.real() > energyDensity  ) { energyDensity = e1.real(); }
      if( e2.real() > energyDensity  ) { energyDensity = e2.real(); }
      if( e3.real() > energyDensity  ) { energyDensity = e3.real(); }
      if( e4.real() > energyDensity  ) { energyDensity = e4.real(); }
      solution = true; 
    } break;
    default:
    {
      //cout << "3 or more positive solution in eDensity Landau frame! " << endl;
      energyDensity = T00();
      solution = false;    
    }
  }
    
  /*
    if(!(nn_e == 1))
    {
    //------------------------------------------------------------------
    //Energiedichte im 1-dimensionalen Fall, dient als Abschätzung
    double root_F;
    double p_OneD,q_OneD,energyDensityOneD;

    p_OneD = T33 - T00;
    q_OneD = - T33*T00 + pow(T30,2);
    energyDensityOneD = -(p_OneD/2.0) + sqrt( pow((p_OneD/2.0),2) - q_OneD); // GeV / fm^3
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    //Newton-Verfahren //übergangsweise
    double initial = energyDensityOneD + 1.0;
    //------------------------
              
    cout << initial << endl;
              
    for(int j = 1; j<100000; j++)
    {
    double x_nv,F,dF;
    x_nv = initial;

    F = pow(x_nv,4) + e_C*pow(x_nv,2) + e_D*x_nv + e_E;
    dF = 4.0*pow(x_nv,3) + 2.0*e_C*x_nv + e_D;

    root_F = initial - F/dF;
    initial = root_F;
    }

    energyDensity = root_F;

    cout << "Error in energy density calculations" << endl;
    cout << "Starting Newton-Verfahren for energy density:" << endl;
    cout << "Solutions for e: " << e1.real() << " " << e2.real() << " " << e3.real() << " " << e4.real() << endl;
    cout << "root_F Newton:   " <<  root_F << endl;
    cout << "eDensity old:    " <<  energyDensity << endl;
    cout << "difference:      " <<  energyDensity - root_F << endl;
    }
  */

  return energyDensity;
}

double tTmunuNmu::LanEck_particleDensity( const VectorTXYZ & v ) const
{
  const double gamma = 1.0/(sqrt(1.0 - v.vec2()));
  return Dot(v, VectorTXYZ(N0(),N1(),N2(),N3()))*gamma;
}

tPimunu tTmunuNmu::LanEck_shearStress( const VectorTXYZ & v ) const
{
  tdeltamunu D( v );
  tPimunu Pimunu;

  Pimunu._PI00 = (D.D$0_0*D.D$0_0 - 1.0/3.0*D.D$0$0*D.D_0_0)*T00()
    + (D.D$0_1*D.D$0_1 - 1.0/3.0*D.D$0$0*D.D_1_1)*T11() + (D.D$0_2*D.D$0_2 - 1.0/3.0*D.D$0$0*D.D_2_2)*T22() + (D.D$0_3*D.D$0_3 - 1.0/3.0*D.D$0$0*D.D_3_3)*T33()
    + 2.0*((D.D$0_1*D.D$0_0 - 1.0/3.0*D.D$0$0*D.D_1_0)*T10() + (D.D$0_2*D.D$0_0 - 1.0/3.0*D.D$0$0*D.D_2_0)*T20() + (D.D$0_3*D.D$0_0 - 1.0/3.0*D.D$0$0*D.D_3_0)*T30())
    + 2.0*((D.D$0_2*D.D$0_1 - 1.0/3.0*D.D$0$0*D.D_2_1)*T21() + (D.D$0_3*D.D$0_1 - 1.0/3.0*D.D$0$0*D.D_3_1)*T31() + (D.D$0_3*D.D$0_2 - 1.0/3.0*D.D$0$0*D.D_3_2)*T32());

  Pimunu._PI11 = (D.D$1_0*D.D$1_0 - 1.0/3.0*D.D$1$1*D.D_0_0)*T00()
    + (D.D$1_1*D.D$1_1 - 1.0/3.0*D.D$1$1*D.D_1_1)*T11() + (D.D$1_2*D.D$1_2 - 1.0/3.0*D.D$1$1*D.D_2_2)*T22() + (D.D$1_3*D.D$1_3 - 1.0/3.0*D.D$1$1*D.D_3_3)*T33()
    + 2.0*((D.D$1_1*D.D$1_0 - 1.0/3.0*D.D$1$1*D.D_1_0)*T10() + (D.D$1_2*D.D$1_0 - 1.0/3.0*D.D$1$1*D.D_2_0)*T20() + (D.D$1_3*D.D$1_0 - 1.0/3.0*D.D$1$1*D.D_3_0)*T30())
    + 2.0*((D.D$1_2*D.D$1_1 - 1.0/3.0*D.D$1$1*D.D_2_1)*T21() + (D.D$1_3*D.D$1_1 - 1.0/3.0*D.D$1$1*D.D_3_1)*T31() + (D.D$1_3*D.D$1_2 - 1.0/3.0*D.D$1$1*D.D_3_2)*T32());

  Pimunu._PI22 = (D.D$2_0*D.D$2_0 - 1.0/3.0*D.D$2$2*D.D_0_0)*T00()
    + (D.D$2_1*D.D$2_1 - 1.0/3.0*D.D$2$2*D.D_1_1)*T11() + (D.D$2_2*D.D$2_2 - 1.0/3.0*D.D$2$2*D.D_2_2)*T22() + (D.D$2_3*D.D$2_3 - 1.0/3.0*D.D$2$2*D.D_3_3)*T33()
    + 2.0*((D.D$2_1*D.D$2_0 - 1.0/3.0*D.D$2$2*D.D_1_0)*T10() + (D.D$2_2*D.D$2_0 - 1.0/3.0*D.D$2$2*D.D_2_0)*T20() + (D.D$2_3*D.D$2_0 - 1.0/3.0*D.D$2$2*D.D_3_0)*T30())
    + 2.0*((D.D$2_2*D.D$2_1 - 1.0/3.0*D.D$2$2*D.D_2_1)*T21() + (D.D$2_3*D.D$2_1 - 1.0/3.0*D.D$2$2*D.D_3_1)*T31() + (D.D$2_3*D.D$2_2 - 1.0/3.0*D.D$2$2*D.D_3_2)*T32());

  Pimunu._PI33 = (D.D$3_0*D.D$3_0 - 1.0/3.0*D.D$3$3*D.D_0_0)*T00()
    + (D.D$3_1*D.D$3_1 - 1.0/3.0*D.D$3$3*D.D_1_1)*T11() + (D.D$3_2*D.D$3_2 - 1.0/3.0*D.D$3$3*D.D_2_2)*T22() + (D.D$3_3*D.D$3_3 - 1.0/3.0*D.D$3$3*D.D_3_3)*T33()
    + 2.0*((D.D$3_1*D.D$3_0 - 1.0/3.0*D.D$3$3*D.D_1_0)*T10() + (D.D$3_2*D.D$3_0 - 1.0/3.0*D.D$3$3*D.D_2_0)*T20() + (D.D$3_3*D.D$3_0 - 1.0/3.0*D.D$3$3*D.D_3_0)*T30())
    + 2.0*((D.D$3_2*D.D$3_1 - 1.0/3.0*D.D$3$3*D.D_2_1)*T21() + (D.D$3_3*D.D$3_1 - 1.0/3.0*D.D$3$3*D.D_3_1)*T31() + (D.D$3_3*D.D$3_2 - 1.0/3.0*D.D$3$3*D.D_3_2)*T32());

  Pimunu._PI10 = (D.D$1_0*D.D$0_0 - 1.0/3.0*D.D$1$0*D.D_0_0)*T00()
    + (D.D$1_1*D.D$0_1 - 1.0/3.0*D.D$1$0*D.D_1_1)*T11() + (D.D$1_2*D.D$0_2 - 1.0/3.0*D.D$1$0*D.D_2_2)*T22() + (D.D$1_3*D.D$0_3 - 1.0/3.0*D.D$1$0*D.D_3_3)*T33()
    + (D.D$1_0*D.D$0_1 - 1.0/3.0*D.D$1$0*D.D_0_1)*T10() + (D.D$1_0*D.D$0_2 - 1.0/3.0*D.D$1$0*D.D_0_2)*T20() + (D.D$1_0*D.D$0_3 - 1.0/3.0*D.D$1$0*D.D_0_3)*T30()
    + (D.D$1_1*D.D$0_0 - 1.0/3.0*D.D$1$0*D.D_1_0)*T10() + (D.D$1_2*D.D$0_0 - 1.0/3.0*D.D$1$0*D.D_2_0)*T20() + (D.D$1_3*D.D$0_0 - 1.0/3.0*D.D$1$0*D.D_3_0)*T30()
    + (D.D$1_2*D.D$0_1 - 1.0/3.0*D.D$1$0*D.D_2_1)*T21() + (D.D$1_3*D.D$0_1 - 1.0/3.0*D.D$1$0*D.D_3_1)*T31() + (D.D$1_3*D.D$0_2 - 1.0/3.0*D.D$1$0*D.D_3_2)*T32()
    + (D.D$1_1*D.D$0_2 - 1.0/3.0*D.D$1$0*D.D_1_2)*T21() + (D.D$1_1*D.D$0_3 - 1.0/3.0*D.D$1$0*D.D_1_3)*T31() + (D.D$1_2*D.D$0_3 - 1.0/3.0*D.D$1$0*D.D_2_3)*T32();

  Pimunu._PI20 = (D.D$2_0*D.D$0_0 - 1.0/3.0*D.D$2$0*D.D_0_0)*T00()
    + (D.D$2_1*D.D$0_1 - 1.0/3.0*D.D$2$0*D.D_1_1)*T11() + (D.D$2_2*D.D$0_2 - 1.0/3.0*D.D$2$0*D.D_2_2)*T22() + (D.D$2_3*D.D$0_3 - 1.0/3.0*D.D$2$0*D.D_3_3)*T33()
    + (D.D$2_0*D.D$0_1 - 1.0/3.0*D.D$2$0*D.D_0_1)*T10() + (D.D$2_0*D.D$0_2 - 1.0/3.0*D.D$2$0*D.D_0_2)*T20() + (D.D$2_0*D.D$0_3 - 1.0/3.0*D.D$2$0*D.D_0_3)*T30()
    + (D.D$2_1*D.D$0_0 - 1.0/3.0*D.D$2$0*D.D_1_0)*T10() + (D.D$2_2*D.D$0_0 - 1.0/3.0*D.D$2$0*D.D_2_0)*T20() + (D.D$2_3*D.D$0_0 - 1.0/3.0*D.D$2$0*D.D_3_0)*T30()
    + (D.D$2_2*D.D$0_1 - 1.0/3.0*D.D$2$0*D.D_2_1)*T21() + (D.D$2_3*D.D$0_1 - 1.0/3.0*D.D$2$0*D.D_3_1)*T31() + (D.D$2_3*D.D$0_2 - 1.0/3.0*D.D$2$0*D.D_3_2)*T32()
    + (D.D$2_1*D.D$0_2 - 1.0/3.0*D.D$2$0*D.D_1_2)*T21() + (D.D$2_1*D.D$0_3 - 1.0/3.0*D.D$2$0*D.D_1_3)*T31() + (D.D$2_2*D.D$0_3 - 1.0/3.0*D.D$2$0*D.D_2_3)*T32();

  Pimunu._PI30 = (D.D$3_0*D.D$0_0 - 1.0/3.0*D.D$3$0*D.D_0_0)*T00()
    + (D.D$3_1*D.D$0_1 - 1.0/3.0*D.D$3$0*D.D_1_1)*T11() + (D.D$3_2*D.D$0_2 - 1.0/3.0*D.D$3$0*D.D_2_2)*T22() + (D.D$3_3*D.D$0_3 - 1.0/3.0*D.D$3$0*D.D_3_3)*T33()
    + (D.D$3_0*D.D$0_1 - 1.0/3.0*D.D$3$0*D.D_0_1)*T10() + (D.D$3_0*D.D$0_2 - 1.0/3.0*D.D$3$0*D.D_0_2)*T20() + (D.D$3_0*D.D$0_3 - 1.0/3.0*D.D$3$0*D.D_0_3)*T30()
    + (D.D$3_1*D.D$0_0 - 1.0/3.0*D.D$3$0*D.D_1_0)*T10() + (D.D$3_2*D.D$0_0 - 1.0/3.0*D.D$3$0*D.D_2_0)*T20() + (D.D$3_3*D.D$0_0 - 1.0/3.0*D.D$3$0*D.D_3_0)*T30()
    + (D.D$3_2*D.D$0_1 - 1.0/3.0*D.D$3$0*D.D_2_1)*T21() + (D.D$3_3*D.D$0_1 - 1.0/3.0*D.D$3$0*D.D_3_1)*T31() + (D.D$3_3*D.D$0_2 - 1.0/3.0*D.D$3$0*D.D_3_2)*T32()
    + (D.D$3_1*D.D$0_2 - 1.0/3.0*D.D$3$0*D.D_1_2)*T21() + (D.D$3_1*D.D$0_3 - 1.0/3.0*D.D$3$0*D.D_1_3)*T31() + (D.D$3_2*D.D$0_3 - 1.0/3.0*D.D$3$0*D.D_2_3)*T32();

  Pimunu._PI21 = (D.D$2_0*D.D$1_0 - 1.0/3.0*D.D$2$1*D.D_0_0)*T00()
    + (D.D$2_1*D.D$1_1 - 1.0/3.0*D.D$2$1*D.D_1_1)*T11() + (D.D$2_2*D.D$1_2 - 1.0/3.0*D.D$2$1*D.D_2_2)*T22() + (D.D$2_3*D.D$1_3 - 1.0/3.0*D.D$2$1*D.D_3_3)*T33()
    + (D.D$2_0*D.D$1_1 - 1.0/3.0*D.D$2$1*D.D_0_1)*T10() + (D.D$2_0*D.D$1_2 - 1.0/3.0*D.D$2$1*D.D_0_2)*T20() + (D.D$2_0*D.D$1_3 - 1.0/3.0*D.D$2$1*D.D_0_3)*T30()
    + (D.D$2_1*D.D$1_0 - 1.0/3.0*D.D$2$1*D.D_1_0)*T10() + (D.D$2_2*D.D$1_0 - 1.0/3.0*D.D$2$1*D.D_2_0)*T20() + (D.D$2_3*D.D$1_0 - 1.0/3.0*D.D$2$1*D.D_3_0)*T30()
    + (D.D$2_2*D.D$1_1 - 1.0/3.0*D.D$2$1*D.D_2_1)*T21() + (D.D$2_3*D.D$1_1 - 1.0/3.0*D.D$2$1*D.D_3_1)*T31() + (D.D$2_3*D.D$1_2 - 1.0/3.0*D.D$2$1*D.D_3_2)*T32()
    + (D.D$2_1*D.D$1_2 - 1.0/3.0*D.D$2$1*D.D_1_2)*T21() + (D.D$2_1*D.D$1_3 - 1.0/3.0*D.D$2$1*D.D_1_3)*T31() + (D.D$2_2*D.D$1_3 - 1.0/3.0*D.D$2$1*D.D_2_3)*T32();

  Pimunu._PI31 = (D.D$3_0*D.D$1_0 - 1.0/3.0*D.D$3$1*D.D_0_0)*T00()
    + (D.D$3_1*D.D$1_1 - 1.0/3.0*D.D$3$1*D.D_1_1)*T11() + (D.D$3_2*D.D$1_2 - 1.0/3.0*D.D$3$1*D.D_2_2)*T22() + (D.D$3_3*D.D$1_3 - 1.0/3.0*D.D$3$1*D.D_3_3)*T33()
    + (D.D$3_0*D.D$1_1 - 1.0/3.0*D.D$3$1*D.D_0_1)*T10() + (D.D$3_0*D.D$1_2 - 1.0/3.0*D.D$3$1*D.D_0_2)*T20() + (D.D$3_0*D.D$1_3 - 1.0/3.0*D.D$3$1*D.D_0_3)*T30()
    + (D.D$3_1*D.D$1_0 - 1.0/3.0*D.D$3$1*D.D_1_0)*T10() + (D.D$3_2*D.D$1_0 - 1.0/3.0*D.D$3$1*D.D_2_0)*T20() + (D.D$3_3*D.D$1_0 - 1.0/3.0*D.D$3$1*D.D_3_0)*T30()
    + (D.D$3_2*D.D$1_1 - 1.0/3.0*D.D$3$1*D.D_2_1)*T21() + (D.D$3_3*D.D$1_1 - 1.0/3.0*D.D$3$1*D.D_3_1)*T31() + (D.D$3_3*D.D$1_2 - 1.0/3.0*D.D$3$1*D.D_3_2)*T32()
    + (D.D$3_1*D.D$1_2 - 1.0/3.0*D.D$3$1*D.D_1_2)*T21() + (D.D$3_1*D.D$1_3 - 1.0/3.0*D.D$3$1*D.D_1_3)*T31() + (D.D$3_2*D.D$1_3 - 1.0/3.0*D.D$3$1*D.D_2_3)*T32();

  Pimunu._PI32 = (D.D$3_0*D.D$2_0 - 1.0/3.0*D.D$3$2*D.D_0_0)*T00()
    + (D.D$3_1*D.D$2_1 - 1.0/3.0*D.D$3$2*D.D_1_1)*T11() + (D.D$3_2*D.D$2_2 - 1.0/3.0*D.D$3$2*D.D_2_2)*T22() + (D.D$3_3*D.D$2_3 - 1.0/3.0*D.D$3$2*D.D_3_3)*T33()
    + (D.D$3_0*D.D$2_1 - 1.0/3.0*D.D$3$2*D.D_0_1)*T10() + (D.D$3_0*D.D$2_2 - 1.0/3.0*D.D$3$2*D.D_0_2)*T20() + (D.D$3_0*D.D$2_3 - 1.0/3.0*D.D$3$2*D.D_0_3)*T30()
    + (D.D$3_1*D.D$2_0 - 1.0/3.0*D.D$3$2*D.D_1_0)*T10() + (D.D$3_2*D.D$2_0 - 1.0/3.0*D.D$3$2*D.D_2_0)*T20() + (D.D$3_3*D.D$2_0 - 1.0/3.0*D.D$3$2*D.D_3_0)*T30()
    + (D.D$3_2*D.D$2_1 - 1.0/3.0*D.D$3$2*D.D_2_1)*T21() + (D.D$3_3*D.D$2_1 - 1.0/3.0*D.D$3$2*D.D_3_1)*T31() + (D.D$3_3*D.D$2_2 - 1.0/3.0*D.D$3$2*D.D_3_2)*T32()
    + (D.D$3_1*D.D$2_2 - 1.0/3.0*D.D$3$2*D.D_1_2)*T21() + (D.D$3_1*D.D$2_3 - 1.0/3.0*D.D$3$2*D.D_1_3)*T31() + (D.D$3_2*D.D$2_3 - 1.0/3.0*D.D$3$2*D.D_2_3)*T32();

  return Pimunu;
}

double tTmunuNmu::LanEck_isotropicPressure( const VectorTXYZ & v) const
{
  tdeltamunu D( v );

  return - 1.0/3.0*(D.D_0_0*T00() + D.D_1_1*T11() + D.D_2_2*T22() + D.D_3_3*T33()) 
    - 2.0/3.0*( D.D_1_0*T10() + D.D_2_0*T20() + D.D_3_0*T30() + D.D_2_1*T21() + D.D_3_1*T31() + D.D_3_2*T32());
}

VectorTXYZ tTmunuNmu::LanEck_energyMomentumFlow( const VectorTXYZ & v) const
{
  tdeltamunu D( v );
  VectorTXYZ u(v);
  u *= 1.0/(sqrt(1.0 - v.vec2())); // multiply with gamma
  u.Minus3(); // change sign of spatial components, since we need u_i
              // and not u$i

  return VectorTXYZ( D.D$0_0*T00()*u(0) + D.D$0_0*T01()*u(1) + D.D$0_0*T02()*u(2) + D.D$0_0*T03()*u(3)
                    +D.D$0_1*T10()*u(0) + D.D$0_1*T11()*u(1) + D.D$0_1*T12()*u(2) + D.D$0_1*T13()*u(3)
                    +D.D$0_2*T20()*u(0) + D.D$0_2*T21()*u(1) + D.D$0_2*T22()*u(2) + D.D$0_2*T23()*u(3)
                    +D.D$0_3*T30()*u(0) + D.D$0_3*T31()*u(1) + D.D$0_3*T32()*u(2) + D.D$0_3*T33()*u(3),
                     D.D$1_0*T00()*u(0) + D.D$1_0*T01()*u(1) + D.D$1_0*T02()*u(2) + D.D$1_0*T03()*u(3)
                    +D.D$1_1*T10()*u(0) + D.D$1_1*T11()*u(1) + D.D$1_1*T12()*u(2) + D.D$1_1*T13()*u(3)
                    +D.D$1_2*T20()*u(0) + D.D$1_2*T21()*u(1) + D.D$1_2*T22()*u(2) + D.D$1_2*T23()*u(3)
                    +D.D$1_3*T30()*u(0) + D.D$1_3*T31()*u(1) + D.D$1_3*T32()*u(2) + D.D$1_3*T33()*u(3),
                     D.D$2_0*T00()*u(0) + D.D$2_0*T01()*u(1) + D.D$2_0*T02()*u(2) + D.D$2_0*T03()*u(3)
                    +D.D$2_1*T10()*u(0) + D.D$2_1*T11()*u(1) + D.D$2_1*T12()*u(2) + D.D$2_1*T13()*u(3)
                    +D.D$2_2*T20()*u(0) + D.D$2_2*T21()*u(1) + D.D$2_2*T22()*u(2) + D.D$2_2*T23()*u(3)
                    +D.D$2_3*T30()*u(0) + D.D$2_3*T31()*u(1) + D.D$2_3*T32()*u(2) + D.D$2_3*T33()*u(3),
                     D.D$3_0*T00()*u(0) + D.D$3_0*T01()*u(1) + D.D$3_0*T02()*u(2) + D.D$3_0*T03()*u(3)
                    +D.D$3_1*T10()*u(0) + D.D$3_1*T11()*u(1) + D.D$3_1*T12()*u(2) + D.D$3_1*T13()*u(3)
                    +D.D$3_2*T20()*u(0) + D.D$3_2*T21()*u(1) + D.D$3_2*T22()*u(2) + D.D$3_2*T23()*u(3)
                    +D.D$3_3*T30()*u(0) + D.D$3_3*T31()*u(1) + D.D$3_3*T32()*u(2) + D.D$3_3*T33()*u(3) );         
}

VectorTXYZ tTmunuNmu::LanEck_particleFlow( const VectorTXYZ & v) const
{
  tdeltamunu D( v );

  return VectorTXYZ( (N0()*D.D$0$0 - N1()*D.D$1$0 - N2()*D.D$2$0 - N3()*D.D$3$0),
                     (N0()*D.D$0$1 - N1()*D.D$1$1 - N2()*D.D$2$1 - N3()*D.D$3$1),
                     (N0()*D.D$0$2 - N1()*D.D$1$2 - N2()*D.D$2$2 - N3()*D.D$3$2),
                     (N0()*D.D$0$3 - N1()*D.D$1$3 - N2()*D.D$2$3 - N3()*D.D$3$3) );
}

double tTmunuNmu::Temperature( const double eDens, const double nDens, const double mass )
{
  if ( FPT_COMP_Z(mass) )
  {
    return eDens / ( 3.0 * nDens );
  }

  const double precision = 0.001;
  const double border = pow(10.0,10);
  const double T_max = eDens / ( 3.0 * nDens );
  const double T_min = std::max( 0.001, T_max - mass/3.0 );

  const double F_min = 3.0 * T_min + mass * gsl_sf_bessel_Kn(1,mass/T_min) / gsl_sf_bessel_Kn(2,mass/T_min) - eDens / nDens;
  double F;
  
  int counter = 0;
  double T = T_max;

  do
  {
    if ( (mass > T * border) || (mass < T * 1/border) )
    {
      return 0.0; // ==> failure
    }
    if (++counter > 1000)
    {
      return 0.0; // ==> failure
    }

    F = 3.0 * T + mass * gsl_sf_bessel_Kn(1,mass/T) / gsl_sf_bessel_Kn(2,mass/T) - eDens / nDens;
    T -= (T-T_min) / (F-F_min) * F;
    
  } while ( fabs(F) > precision );

  if ((T>T_max) || (T<T_min)) return 0.0; // ==> failure
  return T;
}

double tTmunuNmu::Boltzmann_nDens( const double T, const double mass, const double degen )
{
  if ( FPT_COMP_Z(mass) )
  {
    return degen / pow( M_PI, 2 ) * pow( T/0.197, 3 );
  }
  else
  {
    return degen * 0.5 / pow( M_PI, 2 ) * pow( mass, 2 ) * T * gsl_sf_bessel_Kn( 2, mass/T ) / pow( 0.197, 3 );
  }
}

double tTmunuNmu::Boltzmann_eDens( const double T, const double mass, const double degen )
{
  if ( FPT_COMP_Z(mass) )
  {
    return degen / pow( M_PI, 2 ) * pow( T/0.197, 3 ) * 3 * T;
  }
  else
  {
    return degen * 0.5 * T * pow( mass, 2 ) / (pow( M_PI, 2 )*pow( 0.197, 3 ))*
      ( gsl_sf_bessel_Kn( 2, mass/T ) * 3 * T
        + mass * gsl_sf_bessel_Kn( 1, mass/T ) );   
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
