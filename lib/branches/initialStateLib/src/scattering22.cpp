//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering22.cpp $
//$LastChangedDate: 2017-09-11 10:52:20 +0200 (Mo, 11. Sep 2017) $
//$LastChangedRevision: 2630 $
//$LastChangedBy: greif $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <algorithm>

#include "scattering22.h"
#include "lorentz.h"
#include "globalsettings.h"
#include "FPT_compare.h"
#include "random.h"
#include "binary_cross_sections.h"
#include "coupling.h"
#include "interpolation22.h"

using std::cout;
using std::endl;


/** @brief Standard constructor for only light parton scattering. When used, setParameter needs to be called prior to other methods! */
scattering22::scattering22() :
  P1(), P2(), P1cm(), P2cm(), F1(gluon), F2(gluon), 
  M1(), M2(), M3(), M4(),
  md2_gluon_wo_as(0), md2_quark_wo_as(0), 
  s(0), temperature(0), jpsi_dissociation_temperature(),
  isConstantCrossSecJpsi(), constantCrossSecValueJpsi(),
  Kfactor_light( 1.0 ), KggQQbar(0), KgQgQ(0), kappa_gQgQ(0), 
  isotropicCrossSecGQ(false), isConstantCrossSecGQ(false), constantCrossSecValueGQ(0), 
  theI22( NULL )
{
}


/** @brief Standard constructor for heavy quark scattering. When used, setParameter needs to be called prior to other methods! */
scattering22::scattering22( const interpolation22 * const theI22_arg ) :
  P1(), P2(), P1cm(), P2cm(), F1(gluon), F2(gluon), 
  M1(), M2(), M3(), M4(),
  md2_gluon_wo_as(0), md2_quark_wo_as(0), 
  s(0), temperature(0), jpsi_dissociation_temperature(),
  isConstantCrossSecJpsi(), constantCrossSecValueJpsi(),
  Kfactor_light( 1.0 ), KggQQbar(0), KgQgQ(0), kappa_gQgQ(0), 
  isotropicCrossSecGQ(false), isConstantCrossSecGQ(false), constantCrossSecValueGQ(0), 
  theI22( theI22_arg )
{
}


/**
 * This works only for light parton interactions!
 * Call this constructor when creating a scattering22 object with already know properties of a particle pair.
 * Alternatively the parameters can be set by calling #scattering22::setParameter. 
 *
 * @param[in] P1_arg 4-momentum vector of ingoing particle 1
 * @param[in] P2_arg 4-momentum vector of ingoing particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 1
 * @param[in] s_arg mandelstam s
 * @param[in] md2g_arg gluon debye mass squared (not divided by alpha_s)
 * @param[in] md2q_arg quark debye mass squared (not divided by alpha_s)
 */
scattering22::scattering22( const VectorEPxPyPz & P1_arg, 
                            const VectorEPxPyPz & P2_arg,
                            const FLAVOR_TYPE F1_arg, 
                            const FLAVOR_TYPE F2_arg,
                            const double s_arg, 
                            const double md2g_arg, 
                            const double md2q_arg, 
			    const double Kfactor_light_arg ) :
  P1(P1_arg), P2(P2_arg),
  F1(F1_arg), F2(F2_arg), 
  M1(0.0), M2(0.0), M3(), M4(),
  s(s_arg), temperature(0), jpsi_dissociation_temperature(),
  isConstantCrossSecJpsi(), constantCrossSecValueJpsi(), 
  Kfactor_light( Kfactor_light_arg ), KggQQbar(0.0), KgQgQ(0.0), kappa_gQgQ(1.0), 
  isotropicCrossSecGQ (false), isConstantCrossSecGQ(false), constantCrossSecValueGQ(0.0), 
  theI22( NULL )
{
  if( F1_arg > 6 || F2_arg > 6 )
  {
    std::string errMsg = "Heavy quark scattering in scattering22 although heavy quark parameters are not set.";
    throw eScatt22_error( errMsg );
  }
  
  double as = coupling::get_constant_coupling();
  md2_gluon_wo_as = md2g_arg / as;
  md2_quark_wo_as = md2q_arg / as;  
  
  LL_CM.setBetaCM(P1,P2);
  LL_CM.boost(P1,P2, P1cm,P2cm);

}


/**
 * Call this constructor when creating a scattering22 object with already know properties of a particle pair.
 * Alternatively the parameters can be set by calling #scattering22::setParameter.
 *
 * @param[in] P1_arg 4-momentum vector of ingoing particle 1
 * @param[in] P2_arg 4-momentum vector of ingoing particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 1
 * @param[in] s_arg mandelstam s
 * @param[in] md2g_wo_as_arg gluon debye mass squared divided by alpha_s
 * @param[in] md2q_wo_as_arg quark debye mass squared divided by alpha_s
 * @param[in] KggQQbar_arg K factor for process g + g -> Q + Qbar
 * @param[in] KgQgQ_arg K factor for process g + Q -> g + Q
 * @param[in] kappa_gQgQ_arg Kappa for Debye screening for process g + Q -> g + Q, usually 0.2 (Peshier,Gossiaux)
 * @param[in] isConstantCrossSecGQ_arg Whether a constant cross section is employed for process g + Q -> g + Q
 * @param[in] constantCrossSecValueGQ_arg Value of constant cross section for process g + Q -> g + Q
 * @param[in] isotropicCrossSecGQ_arg Value of constant cross section for process g + Q -> g + Q
 * @param[in] temperature_arg (optional, default = 0) temperature of the cell of the particles
 * @param[in] jpsi_dissociation_temperature_arg (optional, default = 0) dissociation temperature of Jpsi
 * @param[in] isConstantCrossSecJpsi_arg (optional, default = false) Whether a constant cross section is employed for process g+Jpsi -> ccb
 * @param[in] constantCrossSecValueJpsi_arg (optional, default = 0) Value of constant cross section for process g+Jpsi -> ccb
 */
scattering22::scattering22( const interpolation22 * const theI22_arg,
                            const VectorEPxPyPz & P1_arg, 
                            const VectorEPxPyPz & P2_arg,
                            const FLAVOR_TYPE F1_arg, 
                            const FLAVOR_TYPE F2_arg, 
                            const double M1_arg, 
                            const double M2_arg, 
                            const double s_arg, 
                            const double md2g_wo_as_arg, 
                            const double md2q_wo_as_arg, 
                            const double KggQQbar_arg, 
                            const double KgQgQ_arg, 
                            const double kappa_gQgQ_arg, 
                            const bool isConstantCrossSecGQ_arg, 
                            const double constantCrossSecValueGQ_arg, 
                            const bool isotropicCrossSecGQ_arg, 
                            const double Kfactor_light_arg,
                            const double temperature_arg, 
                            const double jpsi_dissociation_temperature_arg, 
                            const bool isConstantCrossSecJpsi_arg, 
                            const double constantCrossSecValueJpsi_arg ) :
  P1(P1_arg), P2(P2_arg),
  F1(F1_arg), F2(F2_arg), 
  M1(M1_arg), M2(M2_arg), M3(), M4(),
  md2_gluon_wo_as(md2g_wo_as_arg), md2_quark_wo_as(md2q_wo_as_arg), 
  s(s_arg), 
  temperature(temperature_arg), jpsi_dissociation_temperature(jpsi_dissociation_temperature_arg), 
  isConstantCrossSecJpsi(isConstantCrossSecJpsi_arg), constantCrossSecValueJpsi(constantCrossSecValueJpsi_arg), 
  Kfactor_light( Kfactor_light_arg ), KggQQbar(KggQQbar_arg), KgQgQ(KgQgQ_arg), 
  kappa_gQgQ(kappa_gQgQ_arg), 
  isotropicCrossSecGQ (isotropicCrossSecGQ_arg), 
  isConstantCrossSecGQ(isConstantCrossSecGQ_arg), constantCrossSecValueGQ(constantCrossSecValueGQ_arg), 
  theI22( theI22_arg )
{
  LL_CM.setBetaCM(P1,P2);
  LL_CM.boost(P1,P2, P1cm,P2cm);
}  



scattering22::~scattering22()
{

}


/**
 * This works only for light parton interactions!
 * Used to set the parameters. With this an existing scattering22 object can be re-used for a new particle pair.
 *
 * @param[in] P1_arg 4-momentum vector of ingoing particle 1
 * @param[in] P2_arg 4-momentum vector of ingoing particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 1
 * @param[in] s_arg mandelstam s
 * @param[in] md2g_arg gluon debye mass squared (not divided by alpha_s)
 * @param[in] md2q_arg quark debye mass squared (not divided by alpha_s)
 */
void scattering22::setParameter( const VectorEPxPyPz & P1_arg, 
                                 const VectorEPxPyPz & P2_arg,
                                 const FLAVOR_TYPE F1_arg, 
                                 const FLAVOR_TYPE F2_arg, 
                                 const double s_arg, 
                                 const double vrel_arg, 
                                 const double md2g_arg, 
                                 const double md2q_arg, 
                                 const double Kfactor_light_arg )
{
  if( F1_arg > 6 || F2_arg > 6 )
  {
    std::string errMsg = "Heavy quark scattering in scattering22 although heavy quark parameters are not set.";
    throw eScatt22_error( errMsg );
  }
  
  const double as = coupling::get_constant_coupling(); // alpha_s for gluons and debye mass
  double md2g_wo_as_arg = md2g_arg / as;
  double md2q_wo_as_arg = md2q_arg / as;
  
  setParameter(P1_arg, P2_arg, F1_arg, F2_arg, 0.0, 0.0, s_arg, vrel_arg, md2g_wo_as_arg, md2q_wo_as_arg, 0.0, 0.0, 1.0, false, 0.0, false, Kfactor_light_arg );
}


/**
 * Used to set the parameters. With this an existing scattering22 object can be re-used for a new particle pair.
 *
 * @param[in] P1_arg 4-momentum vector of ingoing particle 1
 * @param[in] P2_arg 4-momentum vector of ingoing particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 2
 * @param[in] M1_arg mass of ingoing particle 1
 * @param[in] M2_arg mass of ingoing particle 2
 * @param[in] s_arg mandelstam s
 * @param[in] md2g_wo_as_arg gluon debye mass squared divided by alpha_s
 * @param[in] md2q_wo_as_arg quark debye mass squared divided by alpha_s
 * @param[in] KggQQbar_arg K factor for process g + g -> Q + Qbar
 * @param[in] KgQgQ_arg K factor for process g + Q -> g + Q
 * @param[in] kappa_gQgQ_arg Kappa for Debye screening for process g + Q -> g + Q, usually 0.2 (Peshier,Gossiaux)
 * @param[in] isConstantCrossSecGQ_arg Whether a constant cross section is employed for process g + Q -> g + Q
 * @param[in] constantCrossSecValueGQ_arg Value of constant cross section for process g + Q -> g + Q
 * @param[in] isotropicCrossSecGQ_arg Value of constant cross section for process g + Q -> g + Q
 * @param[in] temperature_arg (optional, default = 0) temperature of the cell of the particles
 * @param[in] jpsi_dissociation_temperature_arg (optional, default = 0) dissociation temperature of Jpsi
 * @param[in] isConstantCrossSecJpsi_arg (optional, default = false) Whether a constant cross section is employed for process g+Jpsi -> ccb
 * @param[in] constantCrossSecValueJpsi_arg (optional, default = 0) Value of constant cross section for process g+Jpsi -> ccb
 */
void scattering22::setParameter( const VectorEPxPyPz & P1_arg, 
                                 const VectorEPxPyPz & P2_arg,
                                 const FLAVOR_TYPE F1_arg, 
                                 const FLAVOR_TYPE F2_arg, 
                                 const double M1_arg, 
                                 const double M2_arg, 
                                 const double s_arg, 
                                 const double vrel_arg,
                                 const double md2g_wo_as_arg, 
                                 const double md2q_wo_as_arg, 
                                 const double KggQQbar_arg, 
                                 const double KgQgQ_arg, 
                                 const double kappa_gQgQ_arg,
                                 const bool isConstantCrossSecGQ_arg, 
                                 const double constantCrossSecValueGQ_arg, 
                                 const bool isotropicCrossSecGQ_arg, 
                                 const double Kfactor_light_arg,
                                 const double temperature_arg, 
                                 const double jpsi_dissociation_temperature_arg,
                                 const bool isConstantCrossSecJpsi_arg, 
                                 const double constantCrossSecValueJpsi_arg)
{
  P1 = P1_arg;
  P2 = P2_arg;
  F1 = F1_arg;
  F2 = F2_arg;
  s = s_arg;
  vrel = vrel_arg;
  M1 = M1_arg;
  M2 = M2_arg;
  Kfactor_light = Kfactor_light_arg;
  KggQQbar = KggQQbar_arg;
  KgQgQ = KgQgQ_arg;
  kappa_gQgQ = kappa_gQgQ_arg;
  isConstantCrossSecGQ = isConstantCrossSecGQ_arg;
  constantCrossSecValueGQ = constantCrossSecValueGQ_arg;
  isotropicCrossSecGQ = isotropicCrossSecGQ_arg;
  md2_gluon_wo_as = md2g_wo_as_arg;
  md2_quark_wo_as = md2q_wo_as_arg;
  temperature = temperature_arg;
  jpsi_dissociation_temperature = jpsi_dissociation_temperature_arg;
  isConstantCrossSecJpsi = isConstantCrossSecJpsi_arg;
  constantCrossSecValueJpsi = constantCrossSecValueJpsi_arg;
  
  LL_CM.setBetaCM(P1,P2);
  LL_CM.boost(P1,P2, P1cm,P2cm);
}


/**
 * Calculates the cross section for the two incoming particles
 *
 * @param[out] initialStateIndex Integer flag set according to the initial state that has been processed.
 * @return cross section in units 1/GeV^2
 */
double scattering22::getXSection22( int& initialStateIndex ) const
{
  if ( s < ( 1.1 * ns_casc::lambda2 ) )
  {
    return 0;
  }

  // sort F1 and F2 such that comparisons below are easier
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  double cs = 0;

  if ( _F2 > 10 ) // One of the scattering partners is not a parton. These cross sections here are only for parton scatterings except for Jpsi.
  {
    if( _F1 == 0 && static_cast<FLAVOR_TYPE>(_F2) == jpsi ) // Jpsi + g -> ccbar
    {
      xsection_gJpsi_ccbar csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, temperature, jpsi_dissociation_temperature, isConstantCrossSecJpsi, constantCrossSecValueJpsi );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 20;
    }
    else
    {
      cs = 0.0;
      initialStateIndex = -1;
    }
  }
  else if ( (_F1 + _F2) == 0 )  // gg -> gg, gg -> qqbar, gg -> ccbar, gg -> bbbar 
  {
    xsection_gg_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    xsection_gg_qqbar csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    xsection_gg_ccbar csObj3( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );
    xsection_gg_bbbar csObj4( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );

    cs = csObj1.totalCrossSection() + csObj2.totalCrossSection() + csObj3.totalCrossSection() + csObj4.totalCrossSection();
    initialStateIndex = 0;
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar  (only for light quarks)
  {
    if( _F1 >= 7 ) // heavy quarks do not interact with each other
    {
      cs = 0.0;
    }
    else
    {
      xsection_qq_qq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 1;
    }
  }
  else if ( (_F1 * _F2) == 0 ) // gq -> gq, gqbar -> gqbar, gc -> gc,  gcbar -> gcbar, gb -> gb,  gbbar -> gbbar
  {
    if( _F2 == 7 || _F2 == 8 ) // gc -> gc,  gcbar -> gcbar,
    {
      xsection_cg_cg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 5;
    }
    else if( _F2==9 || _F2==10 ) // gb -> gb,  gbbar -> gbbar
    {
      xsection_bg_bg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ  );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 6;
    }
    else // gq -> gq, gqbar -> gqbar
    {
      xsection_qg_qg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 2;
    }
  }
  else if ( (_F2 - _F1) == 1 &&  (_F2 % 2) == 0 )  // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg, qqbar -> ccbar, qqbar -> bbbar, ccbar -> gg, ccbar -> qqbar, bbbar -> gg, bbbar -> qqbar
  {
    if( _F1 == 7 ) // ccbar -> gg, ccbar -> qqbar, ccbar -> g Jpsi
    {
      xsection_ccbar_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );
      xsection_ccbar_qqbar csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
      xsection_ccbar_gJpsi csObj3( s, md2_gluon_wo_as, md2_quark_wo_as, temperature, jpsi_dissociation_temperature, isConstantCrossSecJpsi, constantCrossSecValueJpsi );
      
      cs = csObj1.totalCrossSection() + csObj2.totalCrossSection() + csObj3.totalCrossSection();
      initialStateIndex = 7;
    }
    else if( _F1 == 9 ) // bbbar -> gg, bbbar -> qqbar
    {
      xsection_bbbar_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );
      xsection_bbbar_qqbar csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
      
      cs = csObj1.totalCrossSection() + csObj2.totalCrossSection();
      initialStateIndex = 8;
    }
    else  // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg, qqbar -> ccbar, qqbar -> bbbar
    {
      xsection_qqbar_qqbar csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      xsection_qqbar_qqbarDash csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      xsection_qqbar_gg csObj3( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      xsection_qqbar_ccbar csObj4( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
      xsection_qqbar_bbbar csObj5( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
    
      cs = csObj1.totalCrossSection() + csObj2.totalCrossSection() + csObj3.totalCrossSection() + csObj4.totalCrossSection() + csObj5.totalCrossSection();
      initialStateIndex = 3;
    }
  }
  else // qq' -> qq', qqbar' -> qqbar', qc -> qc, qb -> qb and all anti particles
  {
    if( _F1 < 7 && _F2 < 7 ) // light quarks
    {
      xsection_qqdash_qqdash csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 4;
    }
    else if( _F1 < 7 && ( _F2 == 7 || _F2 == 8 ) )
    {
      xsection_cq_cq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 9;
    }
    else if( _F1 < 7 && ( _F2 == 9 || _F2 == 10 ) )
    {
      xsection_bq_bq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 10;
    }
    else
      cs = 0.0;
  }

  return cs;
}



/**
 * Get the cross section for ELASTIC 2->2 processes
 *
 * @param[out] initialStateIndex Integer flag set according to the initial state that has been processed.
 * @return cross section in units 1/GeV^2
 */
double scattering22::getXSectionElastic( int& initialStateIndex ) const
{
  if ( s < ( 1.1 * ns_casc::lambda2 ) )
  {
    return 0;
  }

  // sort F1 and F2 such that comparisons below are easier
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  double cs = 0;

  if ( _F2 > 10 ) // One of the scattering partners is not a parton. These cross sections here are only for parton scatterings.
  {
    cs = 0.0;
    initialStateIndex = -1;
  }
  else if (( _F1 + _F2 ) == 0 ) // gg -> gg
  {
    xsection_gg_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    initialStateIndex = 0;
    cs = csObj1.totalCrossSection();
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar  (only for light quarks)
  {
    if( _F1 >= 7 ) // heavy quarks do not interact with each other
    {
      cs = 0.0;
    }
    else
    {
      xsection_qq_qq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 1;
    }
  }
  else if ( (_F1 * _F2) == 0 ) // gq -> gq, gqbar -> gqbar, gc -> gc,  gcbar -> gcbar, gb -> gb,  gbbar -> gbbar
  {
    if( _F2 == 7 || _F2 == 8 ) // gc -> gc,  gcbar -> gcbar,
    {
      xsection_cg_cg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 5;
    }
    else if( _F2==9 || _F2==10 ) // gb -> gb,  gbbar -> gbbar
    {
      xsection_bg_bg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 6;
    }
    else // gq -> gq, gqbar -> gqbar
    {
      xsection_qg_qg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 2;
    }
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg
  {
    if( _F1 == 7 ) // no elastic heavy-heavy quark scattering
    {
      cs = 0.0;
      initialStateIndex = 7;
    }
    else if( _F1 == 9 ) // no elastic heavy-heavy quark scattering
    {
      cs = 0.0;
      initialStateIndex = 8;
    }
    else
    {
      xsection_qqbar_qqbar csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      initialStateIndex = 3;
      cs = csObj1.totalCrossSection();
    }
  }
  else // qq' -> qq', qqbar' -> qqbar', qc -> qc, qb -> qb and all anti particles
  {
    if( _F1 < 7 && _F2 < 7 ) // light quarks
    {
      xsection_qqdash_qqdash csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 4;
    }
    else if( _F1 < 7 && ( _F2 == 7 || _F2 == 8 ) )
    {
      xsection_cq_cq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 9;
    }
    else if( _F1 < 7 && ( _F2 == 9 || _F2 == 10 ) )
    {
      xsection_bq_bq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      cs = csObj1.totalCrossSection();
      initialStateIndex = 10;
    }
    else
      cs = 0.0;
  }

  return cs;
}

/**
 * This routine samples the transverse momentum transfer of a given collision isotropically.
 * @param[out] F1_out Sampled flavor of outgoing particle 1
 * @param[out] F2_out Sampled flavor of outgoing particle 2
 * @param[out] M1_out mass of outgoing particle 1
 * @param[out] M2_out mass of outgoing particle 2
 * @param[out] t_hat Mandelstam t of the process
 */
void scattering22::getMomentaAndMasses22_isotropic( FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg, double& M1_out, double& M2_out, double& t_hat )
{
  // First assume elastic scatterings, that is that flavor will not change.
  F1arg = F1;
  F2arg = F2;
  
  M1_out = ParticlePrototype::getMass( F1arg );
  M2_out = ParticlePrototype::getMass( F2arg );
  
  // set masses of outgoing particles in scatt22_object, needed in setNewMomenta22()
  M3 = M1_out;
  M4 = M2_out;
  
  double Pa2 = pow((s+M2_out*M2_out-M1_out*M1_out),2.0)/4.0/s-M2_out*M2_out;
  t_hat = -2.0*Pa2*2.0*ran2();
}

/**
 * This routine samples the transverse momentum transfer of a given collision according to the differential
 * cross section.
 * @param[out] F1_out Sampled flavor of outgoing particle 1
 * @param[out] F2_out Sampled flavor of outgoing particle 2
 * @param[out] M1_out mass of outgoing particle 1
 * @param[out] M2_out mass of outgoing particle 2
 * @param[out] t_hat Mandelstam t of the process
 * @param[out] typ the type of collision, used for analysis purposes (e.g. 221 for gg->gg)
 */
void scattering22::getMomentaAndMasses22(FLAVOR_TYPE& F1_out, FLAVOR_TYPE& F2_out, double& M1_out, double& M2_out, double& t_hat, int& typ)
{
  // First assume elastic scatterings, that is that flavor will not change. For all other cases the flavor will be changed below.
  F1_out = F1;
  F2_out = F2;
  
  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

  if ( _F2 > 10 ) // One of the scattering partners is not a parton. These cross sections here are only for parton scatterings.
  {
    if( _F1 == 0 && static_cast<FLAVOR_TYPE>(_F2) == jpsi ) // Jpsi + g -> ccbar
    {
      xsection_gJpsi_ccbar csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, temperature, jpsi_dissociation_temperature, isConstantCrossSecJpsi, constantCrossSecValueJpsi );
      typ = 2240; // g+Jpsi -> c+cbar
      t_hat = sample_t_hat( csObj1 );
      F1_out = charm;
      F2_out = anti_charm;
    }
    else
    {
      std::string errMsg = "Error in getMomentaAndMasses22().";
      throw eScatt22_error( errMsg );
    }
  }
  else if ( (_F1 + _F2) == 0 )  // gg -> gg, gg -> qqbar, gg -> ccbar, gg -> bbbar 
  {
    xsection_gg_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    xsection_gg_qqbar csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    xsection_gg_ccbar csObj3( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );
    xsection_gg_bbbar csObj4( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );
    double cs1 = csObj1.totalCrossSection();
    double cs2 = csObj2.totalCrossSection();
    double cs3 = csObj3.totalCrossSection();
    double cs4 = csObj4.totalCrossSection();
    
    double select = ran2() * ( cs1 + cs2 + cs3 + cs4 );
    
    if ( select < cs1 )
    {
      typ = 221; // gg -> gg
      t_hat = sample_t_hat( csObj1 );
    }
    else if ( select < ( cs1 + cs2 ) )
    {
      typ = 222; // gg -> qqbar
      t_hat = sample_t_hat( csObj2 );
      sampleFlavor( F1_out, F2_out );
    }
    else if ( select < ( cs1 + cs2 + cs3 ) )
    {
      typ = 2210; // gg -> ccbar
      t_hat = sample_t_hat( csObj3 );
      F1_out = charm;
      F2_out = anti_charm;
    }
    else
    {
      typ = 2220; // gg -> bbbar 
      t_hat = sample_t_hat( csObj4 );
      F1_out = bottom;
      F2_out = anti_bottom;
    }
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar (only for light quarks)
  {
    xsection_qq_qq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    typ = 227;
    t_hat = sample_t_hat( csObj1 );
  }
  else if ( (_F1 * _F2) == 0 ) // gq -> gq, gqbar -> gqbar, gc -> gc,  gcbar -> gcbar, gb -> gb,  gbbar -> gbbar
  {
    if( _F2==7 || _F2==8 ) // gc -> gc,  gcbar -> gcbar,
    {
      xsection_cg_cg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      typ = 2214;
      t_hat = sample_t_hat( csObj1 );
    }
    else if( _F2==9 || _F2==10 ) // gb -> gb,  gbbar -> gbbar
    {
      xsection_bg_bg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      typ = 2224;
      t_hat = sample_t_hat( csObj1 );
    }
    else // gq -> gq, gqbar -> gqbar
    {
      xsection_qg_qg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      typ = 223;
      t_hat = sample_t_hat( csObj1 );
    }
  }
  else if ( (_F2 - _F1) == 1 &&  (_F2 % 2) == 0 ) // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg, qqbar -> ccbar, qqbar -> bbbar, ccbar -> gg, ccbar -> qqbar, bbbar -> gg, bbbar -> qqbar
  {
    if( _F1 == 7 ) // ccbar -> gg, ccbar -> qqbar
    {
      xsection_ccbar_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );
      xsection_ccbar_qqbar csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
      xsection_ccbar_gJpsi csObj3( s, md2_gluon_wo_as, md2_quark_wo_as, temperature, jpsi_dissociation_temperature, isConstantCrossSecJpsi, constantCrossSecValueJpsi );

      double cs1 = csObj1.totalCrossSection();
      double cs2 = csObj2.totalCrossSection();
      double cs3 = csObj3.totalCrossSection();
      
      double select = ran2() * ( cs1 + cs2 + cs3 );
    
      if ( select < cs1 )
      {
        typ = 2212; // ccbar -> gg
        t_hat = sample_t_hat( csObj1 );
        F1_out = gluon;
        F2_out = gluon;
      }
      else if ( select < ( cs1 + cs2 ) )
      {
        typ = 2213; // ccbar -> qqbar
        t_hat = sample_t_hat( csObj2 );
        sampleFlavor( F1_out, F2_out );
      }
      else
      {
        typ = 2241; // ccbar -> g Jpsi
        t_hat = sample_t_hat( csObj3 );
        F1_out = gluon;
        F2_out = jpsi;
      } 
    }
    else if( _F1 == 9 ) // bbbar -> gg, bbbar -> qqbar
    {
      xsection_bbbar_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KggQQbar );
      xsection_bbbar_qqbar csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
      
      double cs1 = csObj1.totalCrossSection();
      double cs2 = csObj2.totalCrossSection();
      
      if ( ran2() < cs1 / (cs1 + cs2) )
      {
        typ = 2222; // bbbar -> gg
        t_hat = sample_t_hat( csObj1 );
        F1_out = gluon;
        F2_out = gluon;
      }
      else
      {
        typ = 2223; // bbbar -> qqbar
        t_hat = sample_t_hat( csObj2 );
        sampleFlavor( F1_out, F2_out );
      }
    }
    else  // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg, qqbar -> ccbar, qqbar -> bbbar
    {
      xsection_qqbar_qqbar csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      xsection_qqbar_qqbarDash csObj2( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      xsection_qqbar_gg csObj3( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      xsection_qqbar_ccbar csObj4( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
      xsection_qqbar_bbbar csObj5( s, md2_gluon_wo_as, md2_quark_wo_as, theI22 );
      
      double cs1 = csObj1.totalCrossSection();
      double cs2 = csObj2.totalCrossSection();
      double cs3 = csObj3.totalCrossSection();
      double cs4 = csObj4.totalCrossSection();
      double cs5 = csObj5.totalCrossSection();
      
      double select = ran2() * ( cs1 + cs2 + cs3 + cs4 + cs5 );
      
      if ( select < cs1 )  // qqbar -> qqbar
      {
        typ = 224;
        t_hat = sample_t_hat( csObj1 );
      }
      else if ( select < ( cs1 + cs2 ) )  // qqbar -> q'qbar'
      {
        typ = 225;
        t_hat = sample_t_hat( csObj2 );
        sampleFlavor( F1_out, F2_out, F1 );  // sample flavor excluding F1 (or anti-F1 = F2)
      }
      else if ( select < ( cs1 + cs2 + cs3 ) )   // qqbar -> gg
      {
        typ = 226;
        t_hat = sample_t_hat( csObj3 );
        F1_out = gluon;
        F2_out = gluon;
      }
      else if ( select < ( cs1 + cs2 + cs3 + cs4 ) )   // qqbar -> ccbar
      {
        typ = 2211; // qqbar -> ccbar
        t_hat = sample_t_hat( csObj4 );
        F1_out = charm;
        F2_out = anti_charm;
      }
      else // qqbar -> bbbar
      {
        typ = 2221; // qqbar -> bbbar
        t_hat = sample_t_hat( csObj5 );
        F1_out = bottom;
        F2_out = anti_bottom;
      }
    }
  }  
  else // qq' -> qq', qqbar' -> qqbar', qc -> qc, qb -> qb and all anti particles
  {
    if( _F1 < 7 && _F2 < 7 ) // light quarks
    {
      typ = 228; // qq' -> qq', qqbar' -> qqbar'
      xsection_qqdash_qqdash csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      t_hat = sample_t_hat( csObj1 );
    }
    else if( _F1 < 7 && ( _F2 == 7 || _F2 == 8 ) )
    {
      typ = 2215; // qc -> qc
      xsection_cq_cq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      t_hat = sample_t_hat( csObj1 );
    }
    else if( _F1 < 7 && ( _F2 == 9 || _F2 == 10 ) )
    {
      typ = 2225; // qb -> qb
      xsection_bq_bq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      t_hat = sample_t_hat( csObj1 );
    }
  }
  
  M1_out = ParticlePrototype::getMass( F1_out );
  M2_out = ParticlePrototype::getMass( F2_out );
  
  // set masses of outgoing particles in scatt22_object, needed in setNewMomenta22()
  M3 = M1_out;
  M4 = M2_out;
}

/**
 * This routine gives the masses for the particles. M3 and M4
 * @param[out] F1_out Sampled flavor of outgoing particle 1
 * @param[out] F2_out Sampled flavor of outgoing particle 2
 * @param[out] M1_out mass of outgoing particle 1
 * @param[out] M2_out mass of outgoing particle 2
 */
void scattering22::getOnlyMasses22(FLAVOR_TYPE& F1_out, FLAVOR_TYPE& F2_out, double& M1_out, double& M2_out)
{
  // First assume elastic scatterings, that is that flavor will not change. For all other cases the flavor will be changed below.
  F1_out = F1;
  F2_out = F2;
  
  M1_out = ParticlePrototype::getMass( F1_out );
  M2_out = ParticlePrototype::getMass( F2_out );
  
  // set masses of outgoing particles in scatt22_object, needed in setNewMomenta22()
  M3 = M1_out;
  M4 = M2_out;
}


/**
 * This works only for light parton interactions!
 * This routine samples the transverse momentum transfer of a given collision according to the differential
 * cross section.
 *
 * @param[out] PT2 transverse momentum transfer squared
 * @param[out] typ the type of collision, used for analysis purposes (e.g. 221 for gg->gg)
 * @param[out] F1arg Sampled flavor of outgoing particle 1
 * @param[out] F2arg Sampled flavor of outgoing particle 2
 */
void scattering22::getMomenta22( double& t_hat, int& typ, FLAVOR_TYPE & F1arg, FLAVOR_TYPE & F2arg )
{
  double M1_out, M2_out;
  
  getMomentaAndMasses22( F1arg, F2arg, M1_out, M2_out, t_hat, typ);
  
  if( M1_out != 0.0 || M2_out != 0.0 )
  {
    std::string errMsg = "Production of massive quarks in scattering22 although heavy quark parameters are not set.";
    throw eScatt22_error( errMsg );
  }
    
}


/**
 * This routine samples the transverse momentum transfer of a given ELASTIC collision according to the differential
 * cross section.
 *
 * @param[out] PT2 transverse momentum transfer squared
 * @param[out] typ the type of collision, used for analysis purposes (e.g. 221 for gg->gg)
 */
void scattering22::getMomentaElastic( double& t_hat, int& typ )
{
  // Only elastic scatterings, that is that flavor and mass will not change.
  M3 = M1;
  M4 = M2;

  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

  if ( (_F1 + _F2) == 0 )  // gg -> gg, gg -> qqbar, gg -> ccbar, gg -> bbbar 
  {
    xsection_gg_gg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    typ = 221; // gg -> gg
    t_hat = sample_t_hat( csObj1 );
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar (only for light quarks)
  {
    xsection_qq_qq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    typ = 227;
    t_hat = sample_t_hat( csObj1 );
  }
  else if ( (_F1 * _F2) == 0 ) // gq -> gq, gqbar -> gqbar, gc -> gc,  gcbar -> gcbar, gb -> gb,  gbbar -> gbbar
  {
    if( _F2==7 || _F2==8 ) // gc -> gc,  gcbar -> gcbar,
    {
      xsection_cg_cg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      typ = 2212;
      t_hat = sample_t_hat( csObj1 );
    }
    else if( _F2==9 || _F2==10 ) // gb -> gb,  gbbar -> gbbar
    {
      xsection_bg_bg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      typ = 2218;
      t_hat = sample_t_hat( csObj1 );
    }
    else // gq -> gq, gqbar -> gqbar
    {
      xsection_qg_qg csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      typ = 223;
      t_hat = sample_t_hat( csObj1 );
    }
  }
  else if ( (_F2 - _F1) == 1 &&  (_F2 % 2) == 0 ) // qqbar -> qqbar
  {
    xsection_qqbar_qqbar csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
    typ = 224;
    t_hat = sample_t_hat( csObj1 );
  }  
  else // qq' -> qq', qqbar' -> qqbar', qc -> qc, qb -> qb and all anti particles
  {
    if( _F1 < 7 && _F2 < 7 ) // light quarks
    {
      typ = 228; // qq' -> qq', qqbar' -> qqbar'
      xsection_qqdash_qqdash csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, Kfactor_light );
      t_hat = sample_t_hat( csObj1 );
    }
    else if( _F1 < 7 && ( _F2 == 7 || _F2 == 8 ) )
    {
      typ = 2213; // qc -> qc
      xsection_cq_cq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      t_hat = sample_t_hat( csObj1 );
    }
    else if( _F1 < 7 && ( _F2 == 9 || _F2 == 10 ) )
    {
      typ = 2219; // qb -> qb
      xsection_bq_bq csObj1( s, md2_gluon_wo_as, md2_quark_wo_as, theI22, KgQgQ, isConstantCrossSecGQ, constantCrossSecValueGQ, kappa_gQgQ, isotropicCrossSecGQ );
      t_hat = sample_t_hat( csObj1 );
    }
  }
}




// void scattering22::getMomenta22_isotropic( double& PT2, int& typ )
// {
//   double u_min = 0.0;
//   double u_max = 1.0;
//   double u;
//
//   
//   // Only elastic scatterings, that is that flavor and mass will not change.
//   M3 = M1;
//   M4 = M2;
// 
//   if (( F1 + F2 ) == 0 ) // gluon,gluon -> gluon,gluon
//   {
//      typ  = 221;// gg->gg
//              u   = u_min + ran2() * (u_max - u_min);
//              PT2 = s / 4.0 *  (1.0 - u*u );
//   }
//   else
//     cout << "Purley gluonic plasma expected when using the isotropic 2->2 routines." << endl;
// }





/**
 * Samples new momentum vectors for the outgoing particles, using a given PT2.
 *
 * @param[out] P1 Momentum vector of outgoing particle 1
 * @param[out] P2 Momentum vector of outgoing particle 2
 * @param[in] R1 Space-time vector of ingoing particle 1
 * @param[in] R2 Space-time vector of ingoing particle 2
 * @param[in] t_hat Mandelstam t of the process
 */
void scattering22::setNewMomenta22( VectorEPxPyPz & P1, 
                                    VectorEPxPyPz & P2,
                                    const VectorTXYZ & R1,
                                    const VectorTXYZ & R2,
                                    const double t_hat)
{
  double PT, PZ, c;
  VectorEPxPyPz PP, TT;
  double M1s,M2s,M3s,M4s; // 1,2 -> 3,4
  VectorTXYZ R1cm, R2cm;

  // squared masses of incoming particles
  M1s=M1*M1;
  M2s=M2*M2;
  
  // squared masses of outgoing particles
  M3s=M3*M3;
  M4s=M4*M4;
  
  // determine PT and PZ of outgoing particle 1
//   // first (old) approach
//   double PT2,c1;
//   c1 = t_hat + 0.5*(s - (M1s+M2s+M3s+M4s) + (M1s-M2s)*(M3s-M4s)/s);
//   PT2 = 0.25*(pow(s-M3s-M4s,2.0)/s - 4.0*M3s*M4s/s - pow(c1,2.0)/(pow(P1cm[0],2.0)-M1s));
//   PT = sqrt(PT2);
//   
// //     PZ = sqrt(pow(P1cm[0],2.0)-M1s) - sqrt(-(t_hat+PT2));
//   PZ = sqrt( ( pow( s - M3s - M4s ,2.0) - 4.0 * M3s * M4s ) / 4.0 / s - PT2 );
  
  // second (new) approach
  double p1, p3, p3z, p3xt, costheta;
  p1 = sqrt( pow( s -M1s -M2s , 2.0) - 4.0 * M1s * M2s ) / 2.0 / sqrt(s);
  p3 = sqrt( pow( s -M3s -M4s , 2.0) - 4.0 * M3s * M4s ) / 2.0 / sqrt(s);
  costheta = ( t_hat + ( pow(s,2.0) - s * ( M1s + M2s + M3s + M4s ) + ( M1s - M2s ) * ( M3s - M4s ) ) / 2.0 / s ) / 2.0 / p1 / p3;
  
  p3z = p3 * costheta;
  p3xt = sqrt( p3*p3 - p3z*p3z );
  
//   // check if both approaches give the same, two approaches are useful for error checking...
//   if( !FPT_COMP_E(fabs(p3z),PZ) )
//     cout << "error pz: " << p3z << "  " << PZ << "  "  << M1s << "  " << M2s  << "  "  << M3s << "  " << M4s << "  " << t_hat << endl;
//   if( !FPT_COMP_E(p3xt,PT) && fabs(t_hat) > 1E-3 ) // numerical problems if t_hat = 0: Due to rounding errors PT can get nan which causes also this error but is not actually an error.
//     cout << "error pt: " << p3xt << "  " << PT << "  "  << M1s << "  " << M2s  << "  "  << M3s << "  " << M4s << "  " << t_hat << endl;

  PZ = p3z;
  PT = p3xt;
  
  //<<---------------------------------------------
  // problem with numerical precision: very small t causes an error: PT is nan
  if( ( std::isnan(PT) || std::isnan(PZ) ) && fabs(t_hat) < 1E-3 )
  {
//     cout << "t very small: t=" << t_hat << ". Set t=0." << endl;
    PT = 0.0;
    PZ = sqrt( ( pow( s - M3s - M4s ,2.0) - 4.0 * M3s * M4s ) / 4.0 / s );
  }
  //----------------------------------------------->>
  
  LL_CM.boost( R1, R2, R1cm, R2cm );

  // set coordinate system

  rotation( P1cm, R1cm, R2cm, PP, TT );
  

  // set new momenta
  P1cm = TT * PT + PP * PZ;
  P2cm = -P1cm;
  c = P1cm.vec2();

  P1cm(0) = sqrt(c+M3s);
  P2cm(0) = sqrt(c+M4s);

  LL_CM.boostInv(P1cm, P2cm, P1, P2);

  
  if(P1(0) != P1(0))
  {
  cout << p1 << " " <<  p3  << " " << p3z  << " " << p3xt  << " " << costheta << endl;
  cout << M3 << " " << M4 << endl;
  //cout << P1cm << " " << R1cm << " " << R2cm << " " << PP << " " << TT << endl;   
  }
}


/**
 * Used by #scattering22::setNewMomenta22
 *
 * @param[in] P Momentum vector of ingoing particle 1
 * @param[in] R1 Space-time vector of ingoing particle 1
 * @param[in] R2 Space-time vector of ingoing particle 2
 * @param[out] PP Direction of momentum (i.e. normalized P)
 * @param[out] TT Direction in the azimuthal plane between the two particles.
 */
void scattering22::rotation( const VectorEPxPyPz & P, const VectorTXYZ & R1, const VectorTXYZ & R2, VectorEPxPyPz & PP, VectorEPxPyPz & TT ) const
{
  double phi, sinus, cosinus, c1, c2, c3;

  c1 = P.vec2();

  VectorTXYZ dR = R2 - R1;
  c2 = Dot3(P, dR);


  PP = P * (1.0/sqrt( c1 ));
  TT = P * c2 - dR * c1;

  c3 = sqrt( TT.vec2() );

  // r2-r1 // p1 => TT[] = 0
  if ( c3 < 1.0e-8 )
  {
    double min = 1.0;
    int mini = -1;
    for ( int i = 1;i <= 3;i++ )
    {
      if ( fabs( PP(i) ) < min )
      {
        min = fabs( PP(i) );
        mini = i;
      }
    }
    switch ( mini )
    {
    case 1:
      TT(1) = 0.0;
      TT(2) =  PP(3);
      TT(3) = -PP(2);
      break;
    case 2:
      TT(1) = -PP(3);
      TT(2) = 0.0;
      TT(3) =  PP(1);
      break;
    case 3:
      TT(1) =  PP(2);
      TT(2) = -PP(1);
      TT(3) = 0.0;
      break;
    default:
      cout << "Error in rotation()" << endl;
    }
  }
  //----------------------

  phi = 2.0 * M_PI * ran2();
  sinus = sin( phi );
  cosinus = cos( phi );

  VectorEPxPyPz transv = TT * cosinus + Cross(PP, TT) * sinus;

  c3 = transv.vec2();
  TT = transv * (1.0/sqrt( c3 ));

}


/**
 * Computes the transport cross section in the case of gluon gluon scattering.
 *
 * @return transport cross section in units 1/GeV^2
 */
double scattering22::getTransportXSection22() const
{
  double as = coupling::get_constant_coupling();
  double prefactor = 36 * M_PI * pow( as, 2 ) / 2;
  double A = log(( s + 4 * md2_gluon_wo_as*as ) / ( 4 * md2_gluon_wo_as*as ) );
  double B = s / ( s + 4 * md2_gluon_wo_as*as );

  return prefactor / s * ( A - B );
}


