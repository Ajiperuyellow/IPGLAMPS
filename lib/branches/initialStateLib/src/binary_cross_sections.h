//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binary_cross_sections.h $
//$LastChangedDate: 2014-12-11 22:00:00 +0100 (Do, 11. Dez 2014) $
//$LastChangedRevision: 2015 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef BINARY_XSECTIONS_H
#define BINARY_XSECTIONS_H

#include <math.h>
#include <iostream>
#include "particleprototype.h"
#include "globalsettings.h"
#include "random.h"
#include "coupling.h"
#include "interpolation22.h"


/** 
 * @brief Prototype for the computation of binary cross sections and
 * the sampling of momentum transfers 
 **/
class xsection_generic
{
public:
  /**
   * @brief Constructor
   * @param[in] _s Mandelstam s
   * @param[in] _md2g_wo_as Gluon debye mass squared without alpha_s
   * @param[in] _md2q_wo_as Quark debye mass squared without alpha_s
   * @param{in] _Kfactor
   */
  xsection_generic(const double _s, 
                   const double _md2g_wo_as, 
                   const double _md2q_wo_as, 
                   const double _Kfactor = 1.0 )
    : s(_s), 
      md2g_wo_as(_md2g_wo_as), 
      md2q_wo_as(_md2q_wo_as), 
      markov_steps( 50 ), 
      Kfactor( _Kfactor ) {};

  /**
   * @brief Destructor
   *
   * we define it as virtual, since we are using polymorph pointers to
   * derived objects.
   **/
  virtual ~xsection_generic() {};
  
  /**
   * @brief Compute the total cross section
   * @return Total cross section in 1/GeV^2
   */
  virtual double totalCrossSection() const { return 0; };
  
  /**
   * @brief Compute the total cross section without prefactors -
   * needed for sampling 
   * @return Total cross section without prefactors
   */
  double totalCrossSection_noPreFactors() const { return 0; };
    
  /**
   * @brief Differential cross section d(sigma)/dt
   * @return Differential cross section d(sigma)/dt in 1/GeV^2
   */
  virtual double differentialCrossSection( const double t ) const { return 0; };
    
  /**
   * @brief Inverse of the integral over the differential cross
   * section without prefactors - needed for sampling 
   *
   * Let the integral over the differential cross section yield
   * Integrate[f(x'),{x',0,x}] = A(x) then this function returns x(A)
   * (without any prefactors as in totalCrossSection_noPreFactors()).
   * A uniformly sampled A (in [0, A_max]) then immediately yields x
   * sampled according to f(x). Here A_max corresponds to the total
   * cross section (without prefactors). 
   *
   * @return Inverse integral over the differential cross section
   */
  double inverseIntegral_noPreFactors( const double _A ) const { return 0; };
    
  /**
   * @brief Returns sampled Mandelstam t. This routine is only used
   * for heavy quarks and qqbar -> qqbarDash 
   * @return sampled Mandelstam t
   */
  virtual double get_mandelstam_t() const { return 0; };
  
  double s; ///<  Mandelstam s
  
protected:

  double md2g_wo_as; ///< Gluon debye mass squared without alpha_s
  double md2q_wo_as; ///< Quark debye mass squared without alpha_s

  const int markov_steps; ///< Number of steps in Markov chain in
                          ///sampling of Mandelstam t (burn-in) 

  double Kfactor; ///< K factor for this process. The total cross
                  ///section is multiplied by this factor. Is used for
                  ///both light and heavy parton processes. 
};

/** 
 * @brief Enhanced version of the generic cross section. 
 *
 * It also stores I22, some additional parameters, etc. This is used
 * in the new light parton scatterings with running coupling. But has
 * no new features, is only less copy paste. 
 */ 
class xsection_generic_enhanced : public xsection_generic
{
public:
  xsection_generic_enhanced(const double _s, 
                            const double _md2g_wo_as = 0.0, 
                            const double _md2q_wo_as = 0.0, 
                            const interpolation22 * const _theI22 = 0, 
                            const double _Kfactor = 1.0, 
                            const bool _isConstantCrossSec = false, 
                            const double _constantCrossSecValue = 0.0, 
                            const double _kappa = 1.0, 
                            const bool _isotropicCrossSec = false )
    : xsection_generic( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ),
      kappa(_kappa),
      isConstantCrossSec(_isConstantCrossSec), 
      constantCrossSecValue(_constantCrossSecValue), 
      isotropicCrossSec(_isotropicCrossSec), 
    theI22(_theI22) {};
  
protected:

  double kappa; ///< Kappa for Debye screening for the process
  
  bool isConstantCrossSec; ///< Whether a constant cross section is
                           ///employed for process 

  double constantCrossSecValue; ///< Value of constant cross section
                                ///for process 

  bool isotropicCrossSec; ///< Whether an isotropic momentum sampling
                          ///is employed for process 
  
  const interpolation22 * const theI22; ///< Pointer to interpolation
                                        ///routine for tabularized
                                        ///total cross section. 
  
  /** @brief Matrix element of the process without prefactor */
  double getMatrixElement(const double) const;
};


/** 
 * @brief Cross section generic class for constant coupling. 
 * 
 * In the new version of this file xsection_generic is intended mainly
 * to serve as a base class for processes with a running
 * coupling. Hence the Debye masses are not multiplied with alpha_s
 * here. Therefore this class here is provided for the old processes
 * with a fixed coupling. This is the old form and will be obsolet in
 * the near future. In this class the Debye mass without coupling is
 * just multiplied with the value for constant coupling. 
 **/ 
class xsection_generic_asConst : public xsection_generic
{
public:
  xsection_generic_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor )
  {
    as_const = coupling::get_constant_coupling();
    md2g = as_const * _md2g_wo_as;
    md2q = as_const * _md2q_wo_as;
  };
  
protected:

  double as_const; ///< Value for constant coupling alpha_s
  double md2g; ///< Gluon debye mass squared (includes alpha_s)
  double md2q; ///< Quark debye mass squared (includes alpha_s)
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gg -> gg processes 
 **/ 
class xsection_gg_gg_asConst : public xsection_generic_asConst
{
public:
  xsection_gg_gg_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * 9 * M_PI * pow( as_const , 2 ) * s / ( 2 * md2g * ( 4*md2g + s) ) );
    //       return 0;
  }
  
  double totalCrossSection_noPreFactors() const  { return ( s / ( 4*md2g + s) ); }
  double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
    
  double differentialCrossSection( const double t ) const
  {
    //       // gives back dsigma/dt (t)
    //       const double qt2 = -( pow(t,2)/ s + t );
    
    // gives back dsigma/dqt2 (qt2)
    double qt2 = -t;
    if( qt2 > s*3.0/4.0 ) // u channel contribution
      qt2 = s - qt2;
    else if( qt2 > s/4.0 ) // no contribution for intermediate t, u (no small angle)
      return 0;
    
    return 9.0 / 2.0 * M_PI * pow( as_const , 2.0 ) / ( pow( qt2 + md2g , 2.0 ) ); // /2.0 because due to assignment of qt2 both channels are considered, t and u
  }
  
private:
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqbar -> gg processes 
 **/ 
class xsection_qqbar_gg_asConst : public xsection_generic_asConst
{
public:
  xsection_qqbar_gg_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * 32 * M_PI * pow( as_const , 2 ) / ( 27 * s ) * log( 1 + s / (4*md2q) ) ) ;
  }
  
  double totalCrossSection_noPreFactors() const  { return ( log( 1 + s / (4*md2q) ) ) ; }
  double inverseIntegral_noPreFactors( const double _A ) const  { return ( md2q * (exp(_A) - 1) ); }
  
  double differentialCrossSection( const double t ) const
  {
    //       // gives back dsigma/dt (t)
    //       const double qt2 = -( pow(t,2)/ s + t );
    
    // gives back dsigma/dqt2 (qt2)
    double qt2 = -t;
    if( qt2 > s*3.0/4.0 ) // u channel contribution
      qt2 = s - qt2;
    else if( qt2 > s/4.0 ) // no contribution for intermediate t, u (no small angle)
      return 0;
    
    return 64.0 / 27.0 / 2.0 * M_PI * pow( as_const , 2.0 ) / s / ( qt2 + md2g ); // /2.0 because due to assignment of qt2 both channels are considered, t and u
  }
  
private:
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqbar -> qqbar processes 
 **/
class xsection_qqbar_qqbar_asConst : public xsection_generic_asConst
{
public:
  xsection_qqbar_qqbar_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * 8 * M_PI * pow( as_const , 2 ) * s / ( 9 * md2g * (4 * md2g + s) ) ) ;
  }
  
  double totalCrossSection_noPreFactors() const  { return  ( s / (4 * md2g + s) ) ; }
  double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
  
  double differentialCrossSection( const double t ) const
  {
    //       const double qt2 = -( pow(t,2)/ s + t );
    const double qt2 = -t;
    if( qt2 > s/4.0 )
      return 0;
    
    return 8.0 / 9.0 * M_PI * pow( as_const , 2.0 ) / ( pow( qt2 + md2g , 2.0 ) );
  }
  
private:
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqbar -> q'qbar' processes 
 **/ 
class xsection_qqbar_qqbarDash_asConst : public xsection_generic_asConst
{
public:
  xsection_qqbar_qqbarDash_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * (ParticlePrototype::N_light_flavor - 1) * 8 * M_PI * pow( as_const , 2 ) / ( 27 * s ) ) ;
  }
  
  double totalCrossSection_noPreFactors() const  { return ( 2 * s / 3); }
  double inverseIntegral_noPreFactors( const double _A ) const
  {
    return (-( s/2 + pow(s,2)/(2.*pow(-6*_A*pow(s,2) + 2*pow(s,3) + sqrt(36*pow(_A,2)*pow(s,4) - 24*_A*pow(s,5) + 5*pow(s,6)), 0.3333333))
               - pow(-6*_A*pow(s,2) + 2*pow(s,3) + sqrt(36*pow(_A,2)*pow(s,4) - 24*_A*pow(s,5) + 5*pow(s,6)), 0.33333333)/2. ));
  }
  double get_mandelstam_t() const {return inverseIntegral_noPreFactors( ran2() * totalCrossSection_noPreFactors() ); }
  
  double differentialCrossSection( const double t ) const
  {
    return 4.0 / 9.0 * M_PI * pow( as_const , 2.0 ) / pow( s , 2.0 ) * ( ( pow( t , 2.0 ) + pow( -s-t , 2.0 ) )  / pow( s , 2.0 ));
  }
  
private:
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qg -> qg processes 
 **/
class xsection_qg_qg_asConst : public xsection_generic_asConst
{
public:
  xsection_qg_qg_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * 2 * M_PI * pow( as_const , 2 ) * s / ( md2g * ( 4*md2g + s) )  ) ;
  }
  
  double totalCrossSection_noPreFactors() const  { return ( s /  (4*md2g + s) ); }
  double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A)  ); }
  
  //     double differentialCrossSection( const double t ) const
  //     {
  //       if( t < - s / 2.0 )
  //         return 0;
  //       
  //       const double qt2 = -( pow(t,2)/ s + t );
  //       if( qt2 > s/4.0 )
  //         return 0;
  //       
  //       double substitution_factor = 1.0 / ( - 4.0 * t / s + 1.0 );
  //       
  //       return substitution_factor * 2.0 * M_PI * pow( as_const , 2.0 ) / ( pow( qt2 + md2g , 2.0 ) );
  //     }
  
  double differentialCrossSection( const double t ) const
  {
    //       const double qt2 = -( pow(t,2)/ s + t );
    const double qt2 = -t;
    if( qt2 > s/4.0 )
      return 0;
    
    return 2.0 * M_PI * pow( as_const , 2.0 ) / ( pow( qt2 + md2g , 2.0 ) );
  }
  
private:
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gg -> qqbar processes 
 **/
class xsection_gg_qqbar_asConst : public xsection_generic_asConst
{
public:
  xsection_gg_qqbar_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * ParticlePrototype::N_light_flavor * M_PI * pow( as_const , 2 ) / ( 3 * s ) * log( 1 + s / (4*md2q) ) ) ;
    //       return 0;
  }
  
  double totalCrossSection_noPreFactors() const  { return ( log( 1 + s / (4*md2q) ) ); }
  double inverseIntegral_noPreFactors( const double _A ) const  { return ( md2q * (exp(_A) - 1) ); }
  
  double differentialCrossSection( const double t ) const
  {
    //       // gives back dsigma/dt (t)
    //       const double qt2 = -( pow(t,2)/ s + t );
    
    // gives back dsigma/dqt2 (qt2)
    double qt2 = -t;
    if( qt2 > s*3.0/4.0 ) // u channel contribution
      qt2 = s - qt2;
    else if( qt2 > s/4.0 ) // no contribution for intermediate t, u (no small angle)
      return 0;
    
    return 1.0 / 3.0 / 2.0 * M_PI * pow( as_const , 2.0 ) / s / ( qt2 + md2q ); // /2.0 because due to assignment of qt2 both channels are considered, t and u
  }
  
private:
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qq -> qq (qbarqbar -> qbarqbar) processes 
 **/
class xsection_qq_qq_asConst : public xsection_generic_asConst
{
public:
  xsection_qq_qq_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * 8 * M_PI * pow( as_const , 2 ) * s / ( 9 * md2g * ( 4*md2g + s) ) ) ;
  }
  
  double totalCrossSection_noPreFactors() const  { return ( s /  (4*md2g + s) ); }
  double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
  
  double differentialCrossSection( const double t ) const
  {
    //       // gives back dsigma/dt (t)
    //       const double qt2 = -( pow(t,2)/ s + t );
    
    // gives back dsigma/dqt2 (qt2)
    double qt2 = -t;
    if( qt2 > s*3.0/4.0 ) // u channel contribution
      qt2 = s - qt2;
    else if( qt2 > s/4.0 ) // no contribution for intermediate t, u (no small angle)
      return 0;
    
    return 16.0 / 9.0 / 2.0 * M_PI * pow( as_const , 2.0 ) / ( pow( qt2 + md2g , 2.0 ) ); // /2.0 because due to assignment of qt2 both channels are considered, t and u
  }
  
private:
};



/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qq' -> qq' processes
 **/
class xsection_qqdash_qqdash_asConst : public xsection_generic_asConst
{
public:
  xsection_qqdash_qqdash_asConst(const double _s, const double _md2g_wo_as, const double _md2q_wo_as, const double _Kfactor = 1.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ) {};
  
  double totalCrossSection() const
  {
    return ( Kfactor * 8 * M_PI * pow( as_const , 2 ) * s / ( 9 * md2g * (4 * md2g + s) ) ) ;
  }
  
  double totalCrossSection_noPreFactors() const  { return  ( s / (4 * md2g + s) ) ; }
  double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
  
  double differentialCrossSection( const double t ) const
  {
    //       const double qt2 = -( pow(t,2)/ s + t );
    const double qt2 = -t;
    if( qt2 > s/4.0 )
      return 0;
    
    return 8.0 / 9.0 * M_PI * pow( as_const , 2.0 ) / ( pow( qt2 + md2g , 2.0 ) );
  }
    
private:
};






//============================================================//
//=====================  heavy quarks  =======================//
//============================================================//

// production
// light or heavy quark production from gluons
/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gg -> QQb processes, Q being a light/heavy quark 
 **/
class xsection_gg_qqbar : public xsection_generic
{
public:
  xsection_gg_qqbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0 )
    : xsection_generic( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ), 
      M_Q( 0.0 ), 
      backreaction( false ),
      flavor( light_quark ), 
      theI22(_theI22) {}; 
    
  /** 
   * @brief Returns the total cross section with the correct summation
   * over the number of final states flavors. 
   **/ 
  double totalCrossSection() const;

  /** 
   * @brief Returns differential cross section, but just for one
   * active flavor in the final state. To get total cross section one
   * must also multiply with the number of allowed final state quark
   * flavors. 
   **/
  double differentialCrossSection( const double t ) const
  {
    double backreaction_factor = 1.0; 
    if( backreaction )
      backreaction_factor = 64.0/9.0; // cs_ccb->gg = (1/2) * 64/9 * 1/chi^2 cs_gg->ccb, the first 1/2 (=1/nu) not needed since just differential and not total cross section (it is just a definition to not multiply 1/nu for differential cross section)
      
    return backreaction_factor * Kfactor * 1.0 / ( 16.0 * M_PI * pow( s , 2.0 ) ) * getMatrixElement( t );
  }
  double get_mandelstam_t() const;
    
protected:
  double M_Q; ///< Heavy quark mass
  bool backreaction; ///< For back reaction a different prefactor is
                     ///needed for detailed balance. The matrix
                     ///element is the same 
  FLAVOR_TYPE flavor; ///< flavor
    
  const interpolation22 * const theI22; ///< Pointer to interpolation
                                        ///routine for tabularized
                                        ///total cross section. 
    
  /** @brief Matrix element of the process */
  double getMatrixElement(const double) const;
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqb -> gg processes
 **/
class xsection_qqbar_gg : public xsection_gg_qqbar
{
public:
  xsection_qqbar_gg(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0 )
    : xsection_gg_qqbar( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor ) { backreaction = true; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gg -> ccb processes  
 **/
class xsection_gg_ccbar : public xsection_gg_qqbar
{
public:
  xsection_gg_ccbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0 )
    : xsection_gg_qqbar( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor ) { M_Q = ParticlePrototype::Mcharm; flavor = charm; backreaction = false; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gg -> bbb processes 
 **/
class xsection_gg_bbbar : public xsection_gg_qqbar
{
public:
  xsection_gg_bbbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0 )
    : xsection_gg_qqbar( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor ) { M_Q = ParticlePrototype::Mbottom; flavor = bottom; backreaction = false; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * ccb -> gg  processes 
 **/
class xsection_ccbar_gg : public xsection_gg_qqbar
{
public:
  xsection_ccbar_gg(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0 )
    : xsection_gg_qqbar( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor ) { M_Q = ParticlePrototype::Mcharm; flavor = charm; backreaction = true; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * bbb -> gg processes 
 **/
class xsection_bbbar_gg : public xsection_gg_qqbar
{
public:
  xsection_bbbar_gg(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0 )
    : xsection_gg_qqbar( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor ) { M_Q = ParticlePrototype::Mbottom; flavor = bottom; backreaction = true; };
};


// light quark (with other flavor) or heavy quark production from light quarks

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqbar -> QQb processes, Q being another species of light/heavy
 * quark  
 **/
class xsection_qqbar_qqbarDash : public xsection_generic
{
public:
  xsection_qqbar_qqbarDash(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0 )
    : xsection_generic( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ), 
      M_Q( 0.0 ), 
      backreaction( false ),
      flavor( light_quark ), 
      theI22(_theI22) {};
    
  /** 
   * @brief Returns the total cross section with the correct summation
   * over the number of final states flavors. 
   **/
  double totalCrossSection() const;

  /** 
   * @brief Returns differential cross section, but just for one
   * active flavor in the final state. 
   *
   * To get total cross section one must also multiply with the number
   * of allowed final state quark flavors.
   **/ 
  double differentialCrossSection( const double t ) const
  {
    return 1.0 / ( 16.0 * M_PI * pow( s , 2.0 ) ) * getMatrixElement( t );
  }
  double get_mandelstam_t() const;
    
protected:
  double M_Q; ///< Heavy quark mass

  bool backreaction; ///< For back reaction a different prefactor is
                     ///needed for detailed balance. The matrix
                     ///element is the same 

  FLAVOR_TYPE flavor; ///< flavor
    
  const interpolation22 * const theI22; ///< Pointer to interpolation
                                        ///routine for tabularized
                                        ///total cross section. 
    
  /** @brief Matrix element of the process */
  double getMatrixElement(const double) const;
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqbar -> ccb processes
 **/
class xsection_qqbar_ccbar : public xsection_qqbar_qqbarDash
{
public:
  xsection_qqbar_ccbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0)
    : xsection_qqbar_qqbarDash( _s, _md2g_wo_as, _md2q_wo_as, _theI22 ) { M_Q = ParticlePrototype::Mcharm; flavor = charm; backreaction = false; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqbar -> bbb processes 
 **/
class xsection_qqbar_bbbar : public xsection_qqbar_qqbarDash
{
public:
  xsection_qqbar_bbbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0)
    : xsection_qqbar_qqbarDash( _s, _md2g_wo_as, _md2q_wo_as, _theI22 ) { M_Q = ParticlePrototype::Mbottom; flavor = bottom; backreaction = false; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * ccb -> qqbar  processes 
 **/
class xsection_ccbar_qqbar : public xsection_qqbar_qqbarDash
{
public:
  xsection_ccbar_qqbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0)
    : xsection_qqbar_qqbarDash( _s, _md2g_wo_as, _md2q_wo_as, _theI22 ) { M_Q = ParticlePrototype::Mcharm; flavor = charm; backreaction = true; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * bbb -> qqbar processes 
 **/
class xsection_bbbar_qqbar : public xsection_qqbar_qqbarDash
{
public:
  xsection_bbbar_qqbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0)
    : xsection_qqbar_qqbarDash( _s, _md2g_wo_as, _md2q_wo_as, _theI22 ) { M_Q = ParticlePrototype::Mbottom; flavor = bottom; backreaction = true; };
};



// elastic scattering
// elastic scattering of heavy or light quarks with gluons

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gQ -> gQ processes, Q being a light/heavy quark 
 **/
class xsection_qg_qg : public xsection_generic
{
public:
  xsection_qg_qg(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSecGQ = false , const double _constantCrossSecValueGQ = 0.0, const double _kappa_gQgQ = 1.0, const bool _isotropicCrossSecGQ = false )
    : xsection_generic( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ),
      M_Q( 0.0 ),
      flavor( light_quark ),
      isConstantCrossSecGQ(_isConstantCrossSecGQ), 
      constantCrossSecValueGQ(_constantCrossSecValueGQ), 
    kappa_gQgQ(_kappa_gQgQ),
    isotropicCrossSecGQ(_isotropicCrossSecGQ), 
    theI22(_theI22) {};
    
  double totalCrossSection() const;
  double differentialCrossSection( const double t ) const
  {
    return Kfactor * 1.0 / ( 16.0 * M_PI * pow( s - pow( M_Q , 2.0 ) , 2.0 ) ) * getMatrixElement( t );
  }
  double get_mandelstam_t() const;
    
protected:

  double M_Q; ///< Heavy quark mass

  FLAVOR_TYPE flavor; ///< flavor (charm or bottom)
    
  bool isConstantCrossSecGQ; ///< Whether a constant cross section is
                             ///employed for process g + Q -> g + Q 

  double constantCrossSecValueGQ; ///< Value of constant cross section
                                  ///for process g + Q -> g + Q 

  double kappa_gQgQ; ///< Kappa for Debye screening for process g + Q
                     ///-> g + Q, usually 0.2 (Peshier,Gossiaux) 

  bool isotropicCrossSecGQ; ///< Whether an isotropic momentum
                            ///sampling is employed for process g + Q
                            ///-> g + Q 
    
  const interpolation22 * const theI22; ///< Pointer to interpolation
                                        ///routine for tabularized
                                        ///total cross section. 
    
  /** @brief Matrix element of the process */
  double getMatrixElement(const double) const;
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gc -> gc processes 
 **/
class xsection_cg_cg : public xsection_qg_qg
{
public:
  xsection_cg_cg(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSecGQ = false , const double _constantCrossSecValueGQ = 0.0, const double _kappa_gQgQ = 1.0, const bool _isotropicCrossSecGQ = false )
    : xsection_qg_qg( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor, _isConstantCrossSecGQ, _constantCrossSecValueGQ, _kappa_gQgQ, _isotropicCrossSecGQ ) { M_Q = ParticlePrototype::Mcharm; flavor = charm; };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gb -> gb processes 
 **/
class xsection_bg_bg : public xsection_qg_qg
{
public:
  xsection_bg_bg(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSecGQ = false , const double _constantCrossSecValueGQ = 0.0, const double _kappa_gQgQ = 1.0, const bool _isotropicCrossSecGQ = false )
    : xsection_qg_qg( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor, _isConstantCrossSecGQ, _constantCrossSecValueGQ, _kappa_gQgQ, _isotropicCrossSecGQ ) { M_Q = ParticlePrototype::Mbottom; flavor = bottom; };
};


// elastic scattering of light quarks with another flavor of light or
// heavy quarks 

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qQ -> qQ processes, Q being a light/heavy quark
 **/
class xsection_qqdash_qqdash : public xsection_generic
{
public:
  xsection_qqdash_qqdash(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSecqQ = false , const double _constantCrossSecValueqQ = 0.0, const double _kappa_qQqQ = 1.0, const bool _isotropicCrossSecqQ = false )
    : xsection_generic( _s, _md2g_wo_as, _md2q_wo_as, _Kfactor ), 
      M_Q( 0.0 ),
      flavor( light_quark ), 
      isConstantCrossSecqQ(_isConstantCrossSecqQ),
      constantCrossSecValueqQ(_constantCrossSecValueqQ), 
    kappa_qQqQ(_kappa_qQqQ), 
    isotropicCrossSecqQ(_isotropicCrossSecqQ),
    theI22(_theI22) {};
    
  double totalCrossSection() const;
  double differentialCrossSection( const double t ) const
  {
    return Kfactor * 1.0 / ( 16.0 * M_PI * pow( s - pow( M_Q , 2.0 ) , 2.0 ) ) * getMatrixElement( t );
  }
  double get_mandelstam_t() const;
    
protected:

  double M_Q; ///< Heavy quark mass

  FLAVOR_TYPE flavor; ///< flavor (charm or bottom)
    
  bool isConstantCrossSecqQ; ///< Whether a constant cross section is
                             ///employed for process q + Q -> q + Q 

  double constantCrossSecValueqQ; ///< Value of constant cross section
                                  ///for process q + Q -> q + Q 

  double kappa_qQqQ; ///< Kappa for Debye screening for process q + Q
                     ///-> q + Q, usually 0.2 (Peshier,Gossiaux) 

  bool isotropicCrossSecqQ; ///< Whether an isotropic momentum
                            ///sampling is employed for process q + Q
                            ///-> q + Q 
    
  const interpolation22 * const theI22; ///< Pointer to interpolation
                                        ///routine for tabularized
                                        ///total cross section. 
    
  /** @brief Matrix element of the process without prefactor */
  double getMatrixElement(const double) const;
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qc -> qc processes 
 **/
class xsection_cq_cq : public xsection_qqdash_qqdash
{
public:
  xsection_cq_cq(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSecqQ = false , const double _constantCrossSecValueqQ = 0.0, const double _kappa_qQqQ = 1.0, const bool _isotropicCrossSecqQ = false )
    : xsection_qqdash_qqdash( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor, _isConstantCrossSecqQ, _constantCrossSecValueqQ, _kappa_qQqQ, _isotropicCrossSecqQ ) { M_Q = ParticlePrototype::Mcharm; flavor = charm; };
};


/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qb -> qb processes 
 **/
class xsection_bq_bq : public xsection_qqdash_qqdash
{
public:
  xsection_bq_bq(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSecqQ = false , const double _constantCrossSecValueqQ = 0.0, const double _kappa_qQqQ = 1.0, const bool _isotropicCrossSecqQ = false )
    : xsection_qqdash_qqdash( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor, _isConstantCrossSecqQ, _constantCrossSecValueqQ, _kappa_qQqQ, _isotropicCrossSecqQ ) { M_Q = ParticlePrototype::Mbottom; flavor = bottom; };
};


// Jpsi production and dissociation

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * g+Psi -> QQb processes, Q being a heavy quark and Psi a heavy meson 
 **/
class xsection_gPsi_QQbar : public xsection_generic_asConst
{
public:
  xsection_gPsi_QQbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const double _temperature = 0.0, const double _dissociation_temperature = 0.0, const bool _isConstantCrossSec = false, const double _constantCrossSecValue = 0.0 )
    : xsection_generic_asConst( _s, _md2g_wo_as, _md2q_wo_as ), 
      M_Q(),
      M_psi(),
      backreaction(),
      heavy_quark_flavor(),
      psi_flavor(),
      temperature(_temperature), 
      dissociation_temperature(_dissociation_temperature), 
    isConstantCrossSec(_isConstantCrossSec), 
    constantCrossSecValue(_constantCrossSecValue),
    e_psi()
  {};
    
  double totalCrossSection() const;
  double get_mandelstam_t() const;
    
protected:

  double M_Q; ///< Heavy quark mass
  double M_psi; ///< Heavy meson mass

  bool backreaction; ///< For back reaction a different prefactor is
                     ///needed for detailed balance. The matrix
                     ///element is the same 

  FLAVOR_TYPE heavy_quark_flavor; ///< flavor of Q

  FLAVOR_TYPE psi_flavor; ///< "flavor" of psi, eg. 100 for Jpsi
    
  double temperature; ///< Temperature of the medium. At higher
                      ///temperature as dissociation temperatur psi
                      ///production from QQbar is not possible. 

  double dissociation_temperature; ///< Dissociation temperature. At
                                   ///higher temperature psi
                                   ///production from QQbar is not
                                   ///possible. 
    
  bool isConstantCrossSec; ///< Whether a constant cross section is
                           ///employed for process. 

  double constantCrossSecValue; ///< Value of constant cross section
                                ///for process. 
    
  double e_psi; ///< Binding energy of psi.
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * g+Psi -> QQb processes, Q being a heavy quark and Psi a heavy meson 
 **/
class xsection_gJpsi_ccbar : public xsection_gPsi_QQbar
{
public:
  xsection_gJpsi_ccbar(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const double _temperature = 0.0, const double _dissociation_temperature = 0.32, const bool _isConstantCrossSec = false, const double _constantCrossSecValue = 0.0 )
    : xsection_gPsi_QQbar( _s, _md2g_wo_as, _md2q_wo_as, _temperature , _dissociation_temperature, _isConstantCrossSec, _constantCrossSecValue ) 
  { 
    //       // values from Kai for J/psi production via c + cb -> J/psi + g
    //       M_Q = 1.87; // GeV, a rather large value, larger than our typical charm mass!!
    //       M_psi = 3.6; // GeV
    //       e_psi = 0.15; // GeV
    //       // with these values above, w_0 = 0.1427
      
    // new values
    M_psi = 3.1; // GeV
    e_psi = 0.1; // GeV -> equivalent to Kai's values
    M_Q = ( e_psi + M_psi ) / 2.0; // GeV
      
    // with these values above, w_0 = 0.1427
    heavy_quark_flavor = charm;
    psi_flavor = jpsi;
    backreaction = false;
  };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * QQb -> g+Psi  processes, Q being a heavy quark and Psi a heavy
 * meson 
 **/ 
class xsection_ccbar_gJpsi : public xsection_gPsi_QQbar
{
public:
  xsection_ccbar_gJpsi(const double _s, const double _md2g_wo_as = 0.0, const double _md2q_wo_as = 0.0, const double _temperature = 0.0, const double _dissociation_temperature = 0.32, const bool _isConstantCrossSec = false, const double _constantCrossSecValue = 0.0 )
    : xsection_gPsi_QQbar( _s, _md2g_wo_as, _md2q_wo_as, _temperature , _dissociation_temperature, _isConstantCrossSec, _constantCrossSecValue ) 
  {
    //       // values from Kai for J/psi production via c + cb -> J/psi + g
    //       M_Q = 1.87; // GeV, a rather large value, larger than our typical charm mass!!
    //       M_psi = 3.6; // GeV
    //       e_psi = 0.15; // GeV
    //       // with these values above, w_0 = 0.1427
      
    // new values
    M_psi = 3.1; // GeV
    e_psi = 0.1; // GeV -> equivalent to Kai's values
    M_Q = ( e_psi + M_psi ) / 2.0; // GeV
      
    heavy_quark_flavor = charm;
    psi_flavor = jpsi;
    backreaction = true;
  };
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * gg -> gg processes 
 **/
class xsection_gg_gg : public xsection_generic_enhanced
{
public:
  xsection_gg_gg(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSec = false , const double _constantCrossSecValue = 0.0, const double _kappa = 1.0, const bool _isotropicCrossSec = false )
    : xsection_generic_enhanced( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor, _isConstantCrossSec, _constantCrossSecValue, _kappa, _isotropicCrossSec ) {};
    
  double totalCrossSection() const;
  double differentialCrossSection( const double t ) const
  {
    return 1.0 / ( 16.0 * M_PI * pow( s , 2.0 ) ) * getMatrixElement( t );
  }
  double get_mandelstam_t() const;
    
protected:
  /** @brief Matrix element of the process without prefactor */
  double getMatrixElement(const double) const;
};

/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qq -> qq processes 
 **/
class xsection_qq_qq : public xsection_generic_enhanced
{
public:
  xsection_qq_qq(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSec = false , const double _constantCrossSecValue = 0.0, const double _kappa = 1.0, const bool _isotropicCrossSec = false )
    : xsection_generic_enhanced( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor, _isConstantCrossSec, _constantCrossSecValue, _kappa, _isotropicCrossSec ) {};
    
  double totalCrossSection() const;
  double differentialCrossSection( const double t ) const
  {
    return 1.0 / ( 16.0 * M_PI * pow( s , 2.0 ) ) * getMatrixElement( t );
  }
  double get_mandelstam_t() const;
    
protected:
  /** @brief Matrix element of the process without prefactor */
  double getMatrixElement(const double) const;
};


/** 
 * @brief Binary cross section and sampling of momentum transfers for
 * qqbar -> qqbar processes 
 **/
class xsection_qqbar_qqbar : public xsection_generic_enhanced
{
public:
  xsection_qqbar_qqbar(const double _s, const double _md2g_wo_as, const double _md2q_wo_as = 0.0, const interpolation22 * const _theI22 = 0, const double _Kfactor = 1.0, const bool _isConstantCrossSec = false , const double _constantCrossSecValue = 0.0, const double _kappa = 1.0, const bool _isotropicCrossSec = false )
    : xsection_generic_enhanced( _s, _md2g_wo_as, _md2q_wo_as, _theI22, _Kfactor, _isConstantCrossSec, _constantCrossSecValue, _kappa, _isotropicCrossSec ) {};
    
  double totalCrossSection() const;
  double differentialCrossSection( const double t ) const
  {
    return 1.0 / ( 16.0 * M_PI * pow( s , 2.0 ) ) * getMatrixElement( t );
  }
  double get_mandelstam_t() const;
    
protected:
  /** @brief Matrix element of the process without prefactor */
  double getMatrixElement(const double) const;
};





/** 
 * @brief Default type trait for sampling in mandelstam t (instead of qt^2) 
 *
 * This template serves the purpose to determine whether a class
 * derived from xsection_generic actually samples in mandelstam t (the
 * default) or directly in the transverse momentum transfer qt^2. 
 * In order to change the default value for a certain process type, a
 * specialized version of this template needs to be created.
 */
template<typename T>
struct returns_mandelstam_t { static const bool value = true; };

/** 
 * @brief Specialized type trait for sampling in mandelstam qt2
 * instead of t 
 **/ 
template<>
struct returns_mandelstam_t<xsection_gg_gg_asConst> { static const bool value = false; };
template<>
struct returns_mandelstam_t<xsection_qqbar_gg_asConst> { static const bool value = false; };
template<>
struct returns_mandelstam_t<xsection_qqbar_qqbar_asConst> { static const bool value = false; };
template<>
struct returns_mandelstam_t<xsection_qg_qg_asConst> { static const bool value = false; };
template<>
struct returns_mandelstam_t<xsection_gg_qqbar_asConst> { static const bool value = false; };
template<>
struct returns_mandelstam_t<xsection_qq_qq_asConst> { static const bool value = false; };
template<>
struct returns_mandelstam_t<xsection_qqdash_qqdash_asConst> { static const bool value = false; };


// /** @brief Sample the transverse momentum transfer qt^2*/
// template<class T>
// inline double sampleBinaryPT2( T& _xsectionObj )
// {
//   // if the cross section object for this process returns mandelstam t instead of qt^2, we need to convert it
//   if ( returns_mandelstam_t<T>::value )
//   {
//     double t = _xsectionObj.inverseIntegral_noPreFactors( ran2() * _xsectionObj.totalCrossSection_noPreFactors() );
//     return ( -( pow(t,2)/ _xsectionObj.s + t ) );
//   }
//   else
//   {
//     return _xsectionObj.inverseIntegral_noPreFactors( ran2() * _xsectionObj.totalCrossSection_noPreFactors() );
//   }
// }

/** 
 * @brief Sample Mandelstam t 
 **/
template<class T>
inline double sample_t_hat( T& _xsectionObj )
{
  // if the cross section object for this process returns mandelstam t instead of qt^2, we need to convert it
  if ( returns_mandelstam_t<T>::value )
  {
    double t = _xsectionObj.get_mandelstam_t();
    return t;
  }
  else
  {
    double PT2 = _xsectionObj.inverseIntegral_noPreFactors( ran2() * _xsectionObj.totalCrossSection_noPreFactors() );
    double s = _xsectionObj.s;
    return ( sqrt( s * ( s/4.0 - PT2 ) ) - s/2.0 ); // this is only true for massless particles
  }
}


inline void sampleFlavor(FLAVOR_TYPE & F1, FLAVOR_TYPE & F2, FLAVOR_TYPE exclude = gluon )
{
  do {
    int flav = static_cast<int>( ParticlePrototype::N_light_flavor * ran2() ) + 1;
    if ( flav > ParticlePrototype::N_light_flavor ) 
      flav = ParticlePrototype::N_light_flavor;
    F2 = static_cast<FLAVOR_TYPE>( 2 * flav );
    F1 = static_cast<FLAVOR_TYPE>( 2 * flav - 1 );
  } while ( F1 == exclude || F2 == exclude );
}


/** 
 * @brief exception class for handling unexpected critical behaviour
 * within binary cross section routines
 **/
class eBinCS_error : public std::runtime_error
{
public:
  explicit eBinCS_error(const std::string& what) : std::runtime_error(what) {};

  virtual ~eBinCS_error() throw() {};
};

#endif
