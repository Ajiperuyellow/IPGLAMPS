//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/particleprototype.h $
//$LastChangedDate: 2019-01-05 18:04:29 +0100 (Sa, 05. Jan 2019) $
//$LastChangedRevision: 2918 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the ParticlePrototype class from which the actual particle classes can be derived
 */

#ifndef PARTICLEPROTOTYPE_H
#define PARTICLEPROTOTYPE_H

#include <iostream>
#include <math.h>
#include <string>
#include "coupling.h"
#include "bampsvector.h"

/**
 * @brief Enumeration type for the flavor of a particle 
 *
 * With this, the names given in the enum list can be used instead of
 * integers (though their use is perfectly legal still).
 */
enum FLAVOR_TYPE {
  gluon,                  ///< 0: g (gluon)
  up,                     ///< 1: u (up)
  anti_up,                ///< 2: ub (anti-up)
  down,                   ///< 3: d (down)
  anti_down,              ///< 4: db (anti-down)
  strange,                ///< 5: s (strange)
  anti_strange,           ///< 6: sb (anti-strange)
  charm,                  ///< 7: c (charm)
  anti_charm,             ///< 8: cb (anti-charm)
  bottom,                 ///< 9: b (bottom)
  anti_bottom,            ///< 10: bb (anti-bottom)
  photon          =   13, ///< 13: p (photon)
  jpsi            =   50, ///< 50: Jpsi
  jpsi_ini        =   51, ///< 51: Jpsi (ini)
  jpsi_sec        =   52, ///< 52: Jpsi (sec)
  electron        = 1000, ///< 1000: electron
  positron        = 1001, ///< 1001: positron
  dmeson_plus     =  411, ///< 411: D+
  dmeson_minus    = -411, ///< -411: D-
  dmeson_zero     =  421, ///< 421: D0
  dmeson_zero_bar = -421, ///< -421: D0bar
  bmeson_plus     =  521, ///< 521: B+
  bmeson_minus    = -521, ///< -521: B-
  bmeson_zero     =  511, ///< 511: B0
  bmeson_zero_bar = -511, ///< -511: B0bar
  // generalized flavor type which should not assigned to an individual
  // particle, but is used for analysis purposes
  quark           =   77, ///< 77: quark = light quark or anti-light-quark
  light_quark     =   88, ///< 88: light quark (generic)
  anti_light_quark=   99, ///< 99: light anti-quark (generic)
  allFlavors      =  111, ///< 111: any flavor (generic)
  dmeson_gen      =  222, ///< 222: D-meson (generic)
  bmeson_gen      =  333, ///< 333: B-meson (generic)
  heavy_quark,            ///< 334: heavy quarks (generic)
  electron_gen,           ///< 335:  all electrons/positrons (generic)
  c_electron,             ///< 336: electrons from charm
  b_electron,             ///< 337: electrons from bottom
  massA           =   13, ///< 13: massA
  massB           =   14, ///< 14: massB
  masslessC       =   15, ///< 15: masslessC
  masslessD       =   16  ///< 16: masslessD
};  


/**
 * @brief Access routine used internally for the #ParticlePrototype output modifiers
 *
 * This is used to get access to std::ios_base::xalloc()
 **/
inline int particle_xalloc() { static int i = std::ios_base::xalloc(); return i; }

/** @brief Output modifier for #ParticlePrototype: only basic info (default) **/
inline std::ostream& particlestyle0(std::ostream& os) { os.iword(particle_xalloc()) = 0; return os; }

/** @brief Output modifier for #ParticlePrototype: also info about rates **/
inline std::ostream& particlestyle1(std::ostream& os) { os.iword(particle_xalloc()) = 1; return os; }

/**
 * @brief Provides generic properties of a particle, can be extended by inheritance
 *
 * This class encapsulates properties of particles that are needed for
 * the simulation process, such as position and momentum variables.
 * Derive from this class if more properties are needed.
 */
class ParticlePrototype
{
public:
  /** @brief Provide standard constructor (for completeness) */
  ParticlePrototype() : 
    Pos( 0,0,0,0 ), 
    Mom( 0,0,0,0 ), 
    m( 0 ), 
    FLAVOR( gluon ),
    cell_id( -1 ), 
    dead( false ), 
    unique_id( -1 ), 
    rate22( 0 ), 
    rate23( 0 ), 
    rate32( 0 ) 
  { };
    
  /** @brief time [fm/c] and spatial [fm] coordinates */
  VectorTXYZ Pos;

  /** @brief energy [GeV] and momentum [GeV/c] */
  VectorEPxPyPz Mom;

  /** @brief mass [GeV] */
  double m;
  
  /** @brief flavor of the particle mapped to int according to FLAVOR_TYPE */
  FLAVOR_TYPE FLAVOR;
  
  /** @brief the cell the particle belongs to */
  int cell_id;
  
  /** @brief switch that indicates whether this particle has been annihilated within this timestep */
  bool dead;
    
  /** @brief unique particle ID */
  int unique_id;

  /** @brief counter for unique particle IDs (static) */
  static int unique_id_counter;
  
  /** @brief 2->2 rate currently associated with this particle [GeV] */
  double rate22;
  /** @brief 2->3 rate currently associated with this particle [GeV] */
  double rate23;
  /** @brief 3->2 rate currently associated with this particle [GeV] */
  double rate32;


 /** 
   * @brief Map a flavor to a generic particle type (gluon, quark, anti-quark, etc.). 
   * 
   * Do not group heavy quarks since we are interested in rates for
   * charm and bottom individually  
   */ 
  static int mapToPDG( const FLAVOR_TYPE _flav )
  {
    switch ( _flav )
    {
      case up:
        return 2;
      case down:
        return 1;
      case strange:
        return 3;

      case anti_up:
        return -2;
      case anti_down:
        return -1;
      case anti_strange:
        return -3;

      case gluon:
        return 21;
      
      default:
        return 0;
    }
  }
  


  
  /** 
   * @brief Map a flavor to a generic particle type (gluon, quark, anti-quark, etc.). 
   * 
   * Do not group heavy quarks since we are interested in rates for
   * charm and bottom individually  
   */ 
  static FLAVOR_TYPE mapToGenericFlavorType( const FLAVOR_TYPE _flav )
  {
    switch ( _flav )
    {
      case up:
      case down:
      case strange:
        return light_quark;

      case anti_up:
      case anti_down:
      case anti_strange:
        return anti_light_quark;

      case dmeson_minus:
      case dmeson_plus:
      case dmeson_zero:
      case dmeson_zero_bar:
        return dmeson_gen;

      case bmeson_minus:
      case bmeson_plus:
      case bmeson_zero:
      case bmeson_zero_bar:
        return bmeson_gen;

      case electron:
      case positron:
      case c_electron:
      case b_electron:
        return electron_gen;

      case jpsi_ini:
      case jpsi_sec:
        return jpsi;

      default:
        return _flav;
    }
  }
  
  /** 
   * @brief Map a flavor to generic particle type quark with no difference quark-antiquark.
   */   
  static FLAVOR_TYPE mapQuarkGluon( const FLAVOR_TYPE _flav )
  {
    switch ( _flav )
    {
      case up:
      case down:
      case strange:
      case anti_up:
      case anti_down:
      case anti_strange:
        return quark;  
        
      default:
        return _flav;        
    }
  }
  /** @brief Return the name of the flavor as a string */
  static std::string getName( const FLAVOR_TYPE _flav );
    

  // /** @brief assignment operator */
  //Particle& operator=(Particle const &rhs);
    
    
  /** @brief Whether particle is a heavy quark (charm or bottom) */
  inline bool isHeavyQuark() const 
  {
    return ( FLAVOR >= 7 && FLAVOR <= 10 );
  };  

  /** @brief Whether particle is a heavy quark (charm or bottom) */
  static bool isHeavyQuark( FLAVOR_TYPE _flavor )
  {
    return ( _flavor >= 7 && _flavor <= 10 );
  };  

  /** @brief Whether particle is gluon or light quark */
  inline bool isLightFlavor() const
  {
    return (FLAVOR <= 2 * N_light_flavor);
  }

  /** @brief Whether particle is gluon */
  inline bool isGluon() const
  {
    return (FLAVOR == 0);
  }

  /** @brief maximum number of flavors that are considered to be light */
  static int max_N_light_flavor;

  /** @brief number of active light flavors 
   * ( only gluons (0), up (1), down (2), strange (3))
   */
  static int N_light_flavor;

  /** @brief number of active heavy flavors 
   * ( 0: no charm and bottom, 1: only charm, 2: charm and bottom)
   */
  static int N_heavy_flavor;

  /** @brief number of active psi states
   * ( 0: no psi states, 1: only jpsi ), can be extended to psi prime, chi_c
   */
  static int N_psi_states;

  /** @brief Set  number of active light flavors */
  static void set_N_light_flavor( const int _N ) { N_light_flavor = _N;  coupling::set_Nflavor( _N ); };

  /** @brief Set number of active heavy flavors */
  static void set_N_heavy_flavor( const int _N ) { N_heavy_flavor = _N; };

  /** @brief Set number of active flavors (both light and heavy) */
  static void set_N_flavor( const int _N );

  /** @brief Set number of active psi states */
  static void set_N_psi_states( const int _N ) { N_psi_states = _N; };
    
  /** @brief returns mass of given flavor */
  static double getMass( const FLAVOR_TYPE _flav );
    
  /** @brief Mass of charm quark */
  static double Mcharm; // GeV

  /** @brief Mass of bottom quark */
  static double Mbottom; // GeV

  /** @brief Mass of charged D meson */
  static double MchargedD;

  /** @brief Mass of neutral D meson */
  static double MneutralD;

  /** @brief Mass of charged B meson */
  static double MchargedB;

  /** @brief Mass of neutral B meson */
  static double MneutralB;

  /** @brief Mass of J/psi */
  static double Mjpsi; // GeV

  /** @brief Set mass of charm quark */
  static void setCharmMass( const double _M ) { Mcharm = _M; };

  /** @brief Set mass of bottom quark */
  static void setBottomMass( const double _M ) { Mbottom = _M; };

  /** @brief Set mass of J/psi */
  static void setJpsiMass( const double _M ) { Mjpsi = _M; };
    
  /** @brief Mass of massA */
  static double MmassA; // GeV

  /** @brief Mass of massB */
  static double MmassB; // GeV

  /** @brief Mass of masslessC */
  static double MmasslessC;

  /** @brief Mass of masslessD */
  static double MmasslessD;

  /**
   * @brief Propagte this particle till given time
   *
   * This is a shortcut for
   *   Pos += Mom * ( T - Pos.T )/ Mom.E
   */
  void Propagate( const double time);

  /**
   * @brief Propagte this particle till given time
   *
   * This is a shortcut for
   *   Pos += Mom * ( T - Pos.T )/ Mom.E
   *
   * returns also the traveled distance
   */
  void Propagate( const double time, double &distance);

  /** 
   * @brief The standard output routine 
   *
   * The output may be changed by the modifiers:
   *  - #particlestyle0 (default): short version, only basic info
   *  - #particlestyle1 : longer version, also info about rates
   *
   * Example:
   * ~~~
   * std::cout << particlestyle1 << part << std::endl;
   * ~~~
   */
  friend std::ostream& operator<<(std::ostream &os, const ParticlePrototype &obj);


private:
  
};


#endif
