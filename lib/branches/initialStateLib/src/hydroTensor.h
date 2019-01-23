//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/hydroTensor.h $
//$LastChangedDate: 2016-05-08 08:46:53 +0200 (So, 08. Mai 2016) $
//$LastChangedRevision: 2340 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef HYDROTENSOR_H
#define HYDROTENSOR_H

//#include <fstream>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <gsl/gsl_sf_bessel.h>

#include "bampsvector.h"
#include "FPT_compare.h"

 //----------------------------------------------------------------------------
 //----------------------------------------------------------------------------
// some forward definitions; for use before its actual definition.
class tTmunuNmu;


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * @brief Class to encapsulate Tmunu and Nmu
 *
 * This provides the basic functionality.
 * It has to be extended for calculating hydro properties.
 *
 * The calculated tensors are \f$T^{\mu\nu}\f$ and
 * \f$N^{\mu}\f$.
 **/
class tTmunuNmu_Base
{
public:
  /**
   * @brief constructor
   **/
  tTmunuNmu_Base( ) :
    _T00( ),
    _T10( ),
    _T20( ),
    _T30( ),
    _T11( ),
    _T22( ),
    _T33( ),
    _T21( ),
    _T31( ),
    _T32( ),
    _N0( ),
    _N1( ),
    _N2( ),
    _N3( ),
    _NN( )
  { };

  /**
   * @brief constructor from momentum
   **/
  tTmunuNmu_Base(const VectorEPxPyPz & Mom) :
    _T00( Mom.E() ),
    _T10( Mom.Px() ),
    _T20( Mom.Py() ),
    _T30( Mom.Pz() ),
    _T11( Mom.Px2()/Mom.E() ),
    _T22( Mom.Py2()/Mom.E() ),
    _T33( Mom.Pz2()/Mom.E() ),
    _T21( Mom.Py()*Mom.Px()/Mom.E() ),
    _T31( Mom.Pz()*Mom.Px()/Mom.E() ),
    _T32( Mom.Pz()*Mom.Py()/Mom.E() ),
    _N0( 1.0 ),
    _N1( Mom.Px()/Mom.E()),
    _N2( Mom.Py()/Mom.E()),
    _N3( Mom.Pz()/Mom.E()),
    _NN( 1 )
  { };

  /**
   * @brief add the momentum for a given particle type
   *
   * @param[in] Mom 4-momentum
   * @param[in] fak multiplicative factor to use, probably 1/(cell
   * volume * #testparticles)
   **/
  void add(const VectorEPxPyPz & Mom, const double fak);

  /** @brief get component */
  double T00() const { return _T00; };
  /** @brief get component */
  double T11() const { return _T11; };
  /** @brief get component */
  double T22() const { return _T22; };
  /** @brief get component */
  double T33() const { return _T33; };
  /** @brief get component */
  double T10() const { return _T10; };
  /** @brief get component (for completeness) */
  double T01() const { return _T10; };
  /** @brief get component */
  double T20() const { return _T20; };
  /** @brief get component (for completeness) */
  double T02() const { return _T20; };
  /** @brief get component */
  double T30() const { return _T30; };
  /** @brief get component (for completeness) */
  double T03() const { return _T30; };
  /** @brief get component */
  double T21() const { return _T21; };
  /** @brief get component (for completeness) */
  double T12() const { return _T21; };
  /** @brief get component */
  double T31() const { return _T31; };
  /** @brief get component (for completeness) */
  double T13() const { return _T31; };
  /** @brief get component */
  double T32() const { return _T32; };
  /** @brief get component (for completeness) */
  double T23() const { return _T32; };
  /** @brief get component */
  double N0() const { return _N0; };
  /** @brief get component */
  double N1() const { return _N1; };
  /** @brief get component */
  double N2() const { return _N2; };
  /** @brief get component */
  double N3() const { return _N3; };
  /** @brief get component */
  int NN() const { return _NN; };

  /** @brief The Increment Operator */
  inline tTmunuNmu_Base & operator += (const tTmunuNmu_Base &);

  /**
   * @brief The Multiplicated Operator
   *
   * This is the multiplication with a scalar.
   * We intentionaly do not implement the division operator, because
   * it is better to use "*(1/a)" instead.
   **/
  inline tTmunuNmu_Base & operator *= (const double a);

  /**
   * @brief The Times Operator
   *
   * This is the multiplication with a scalar.
   * We intentionaly do not implement the division operator, because
   * it is better to use "*(1/a)" instead.
   **/
  inline tTmunuNmu_Base operator * (const double a) const
  {
    return tTmunuNmu_Base(*this) *= a;
  }

  /**
   * @brief The standard output routine
   *
   * Prints all information in one single line
   **/
  friend std::ostream& operator<<(std::ostream &os, const tTmunuNmu_Base &obj)
  {
    os << obj.T00() << " "
       << obj.T11() << " "
       << obj.T22() << " "
       << obj.T33() << " "
       << obj.T10() << " "
       << obj.T20() << " "
       << obj.T30() << " "
       << obj.T21() << " "
       << obj.T31() << " "
       << obj.T32() << " "
       << obj.N0() << " "
       << obj.N1() << " "
       << obj.N2() << " "
       << obj.N3() << " "
       << obj.NN();
    return os;
  }

  /**
   * @brief The standard input routine
   **/
  friend std::istream& operator>>(std::istream &is, tTmunuNmu_Base &obj)
  {
    is >> obj._T00;
    is >> obj._T11;
    is >> obj._T22;
    is >> obj._T33;
    is >> obj._T10;
    is >> obj._T20;
    is >> obj._T30;
    is >> obj._T21;
    is >> obj._T31;
    is >> obj._T32;
    is >> obj._N0;
    is >> obj._N1;
    is >> obj._N2;
    is >> obj._N3;
    is >> obj._NN;

    return is;
  }

  /**
   * @brief return a string given all variable names
   *
   * the given integer is the starting number (and increased at exit)
   **/
  static std::string header(int & nr);

protected:
  double _T00;
  double _T10;
  double _T20;
  double _T30;
  double _T11;
  double _T22;
  double _T33;
  double _T21;
  double _T31;
  double _T32;
  double _N0;
  double _N1;
  double _N2;
  double _N3;
  int _NN;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * @brief Class to hold the shear stress tensor
 *
 * \f${\pi}^{\mu\nu} =  \left[\frac{1}{2}\left( \Delta^{\mu}_{\alpha}
 * \Delta^{\nu}_{\beta} + \Delta^{\nu}_{\alpha} \Delta^{\mu}_{\beta}
 * \right) - \frac{1}{3} \Delta^{\mu \nu} \Delta_{\alpha \beta} \right]
 * T^{\alpha \beta}\f$
 *
 **/
class tPimunu
{
public:
  /**
   * @brief constructor
   **/
  tPimunu( ) :
    _PI00( ),
    _PI10( ),
    _PI20( ),
    _PI30( ),
    _PI11( ),
    _PI22( ),
    _PI33( ),
    _PI21( ),
    _PI31( ),
    _PI32( )
  {};

  /**
   * @brief return the contraction PimunuPimunu
   **/
  double PImunuPImunu() const
  {
    return
      pow( _PI00,2 ) - pow( _PI10,2 ) - pow( _PI20,2 ) - pow( _PI30,2 )
      - pow( _PI10,2 ) + pow( _PI11,2 ) + pow( _PI21,2 ) + pow( _PI31,2 )
      - pow( _PI20,2 ) + pow( _PI21,2 ) + pow( _PI22,2 ) + pow( _PI32,2 )
      - pow( _PI30,2 ) + pow( _PI31,2 ) + pow( _PI32,2 ) + pow( _PI33,2 );
  }

  /**
   * @brief The standard output routine
   *
   * Prints all information in one single line
   **/
  friend std::ostream& operator<<(std::ostream &os, const tPimunu &obj)
  {
    os << obj._PI00 << " "
       << obj._PI11 << " "
       << obj._PI22 << " "
       << obj._PI33 << " "
       << obj._PI10 << " "
       << obj._PI20 << " "
       << obj._PI30 << " "
       << obj._PI21 << " "
       << obj._PI31 << " "
       << obj._PI32;
    return os;
  }

  /**
   * @brief return a string given all variable names
   *
   * the given integer is the starting number (and increased at exit)
   **/
  static std::string header(int & nr);

  /**
   * @brief declare tTmunuNmu as friend
   *
   * now tTmunuNmu can access the data directly.
   **/
  friend class tTmunuNmu;

protected:
  double _PI00;
  double _PI10;
  double _PI20;
  double _PI30;
  double _PI11;
  double _PI22;
  double _PI33;
  double _PI21;
  double _PI31;
  double _PI32;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * @brief Class to hold the transversal projection tensor, connected
 * with a velocity vector
 *
 * \f$\Delta^{\mu\nu}=g^{\mu\nu}-u^{\mu}u^{\nu}\f$
 **/
class tdeltamunu
{
public:
  /**
   * @brief constructor
   **/
  tdeltamunu( const VectorTXYZ & v );

  /**
   * @brief declare tTmunuNmu as friend
   *
   * now tTmunuNmu can access the data directly.
   **/
  friend class tTmunuNmu;

protected:
  double D_0_0,D$0_0,D_0$0,D$0$0;
  double D_1_1,D$1$1,D_2_2,D$2$2,D_3_3,D$3$3;
  double D$1_1,D$2_2,D$3_3;
  double D$1_0,D$1$0,D$0$1,D$2_0,D$2$0,D$0$2,D$3_0,D$3$0,D$0$3;
  double D$0_1,D_1_0,D_0_1,D$0_2,D_2_0,D_0_2,D$0_3,D_3_0,D_0_3;
  double D_2_1,D$2$1,D_1_2,D$1$2,D_3_1,D$3$1,D_1_3,D$1$3,D_3_2,D$3$2,D_2_3,D$2$3;
  double D$2_1,D$1_2,D$3_1,D$1_3,D$3_2,D$2_3;
};


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/**
 * @brief Class to extend Tmunu with hydro calculations
 *
 * Some calculations are still not possible, since we do not know,
 * which particle class is accumulated in the stored vales.
 **/
class tTmunuNmu: public tTmunuNmu_Base
{
public:
  /**
   * @brief constructor
   **/
  tTmunuNmu( ) :
    tTmunuNmu_Base( )
  { };


  /**
   * @brief Gives the energy density in the Landau frame.
   *
   * @return velocity vector beta, v[0]==1
   */
  VectorTXYZ Eckart_velocity( ) const;

  /**
   * @brief Gives the velocities in the Landau frame.
   *
   * @return velocity vector beta, v[0]==1
   */
  VectorTXYZ Landau_velocity( const double energyDensity ) const;

  /**
   * @brief Gives the energy density in the Eckart frame.
   *
   * @param[in] v velocity vector beta, v[0]==1
   */
  double Eckart_energyDensity( const VectorTXYZ & v ) const;

  /**
   * @brief Gives the energy density in the Landau frame.
   *
   * @param[out] solution
   */
  double Landau_energyDensity( bool & solution ) const;

  /**
   * @brief Gives the particle density in Eckart or Landau frame.
   *
   * @return particle density
   */
  double LanEck_particleDensity( const VectorTXYZ & v ) const;

  /**
   * @brief Gives the components of the shear stress tensor in
   * Eckart or Landau frame.
   */
  tPimunu LanEck_shearStress( const VectorTXYZ & v ) const;

  /**
   * @brief Gives the isotropic pressure in Eckart or Landau frame.
   *
   * \f$P = p + \Pi = -1/3 \Delta_{\mu\nu} T^{\mu\nu}\f$
   */
  double LanEck_isotropicPressure( const VectorTXYZ & v ) const;

  /**
   * @brief Gives the energy momentum flow in Eckart or Landau frame.
   *
   * \f$W^\mu = \Delta^{\mu\alpha} T_{\alpha\beta} u^{\beta}\f$
   */
  VectorTXYZ LanEck_energyMomentumFlow( const VectorTXYZ & v ) const;

  /**
   * @brief Gives the particle flow in Eckart or Landau frame.
   */
  VectorTXYZ LanEck_particleFlow( const VectorTXYZ & v ) const;

  /**
   * @brief Calculate the temperature for given eDens, nDens and mass
   */
  static double Temperature( const double eDens, const double nDens, const double mass );

  /**
   * @brief Calculate the particle density according Boltzmann distribution
   *
   * return value given in 1/fm^3
   */
  static double Boltzmann_nDens( const double T, const double mass, const double degen );

  /**
   * @brief Calculate the energy density according Boltzmann distribution
   *
   * return value given in GeV/fm^3
   */
  static double Boltzmann_eDens( const double T, const double mass, const double degen );



};



//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @brief exception class for handling unexpected critical behaviour
    within hydro tensors  */
class eHydroTensors : public std::runtime_error
{
public:
  explicit eHydroTensors(const std::string& what) : std::runtime_error(what) {};

  virtual ~eHydroTensors() throw() {};
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

inline tTmunuNmu_Base & tTmunuNmu_Base::operator += (const tTmunuNmu_Base & q)
{
  _T00 += q._T00;
  _T10 += q._T10;
  _T20 += q._T20;
  _T30 += q._T30;
  _T11 += q._T11;
  _T22 += q._T22;
  _T33 += q._T33;
  _T21 += q._T21;
  _T31 += q._T31;
  _T32 += q._T32;
  _N0 += q._N0;
  _N1 += q._N1;
  _N2 += q._N2;
  _N3 += q._N3;
  _NN += q._NN;
  return *this;
}

inline tTmunuNmu_Base & tTmunuNmu_Base::operator *= (const double a)
{
  _T00 *= a;
  _T10 *= a;
  _T20 *= a;
  _T30 *= a;
  _T11 *= a;
  _T22 *= a;
  _T33 *= a;
  _T21 *= a;
  _T31 *= a;
  _T32 *= a;
  _N0 *= a;
  _N1 *= a;
  _N2 *= a;
  _N3 *= a;
  _NN *= a;
  return *this;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#endif
