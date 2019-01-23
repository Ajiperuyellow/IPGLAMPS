//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/particleprototype.cpp $
//$LastChangedDate: 2017-09-11 10:52:20 +0200 (Mo, 11. Sep 2017) $
//$LastChangedRevision: 2630 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Definitions for the ParticeProtoype class
 */


#include <iostream>
#include "particleprototype.h"

/** @brief maximum number of flavors that are considered to be light*/
int ParticlePrototype::max_N_light_flavor = 3;
/** @brief number of active light flavors */
int ParticlePrototype::N_light_flavor = 3;
/** @brief number of active heavy flavors */
int ParticlePrototype::N_heavy_flavor = 0;
/** @brief number of active psi states like J/psi et al. */
int ParticlePrototype::N_psi_states = 0;

/** @brief Initialize static member for mass of charm quark */
double ParticlePrototype::Mcharm = 1.5; // GeV
/** @brief Initialize static member for mass of bottom quark */
double ParticlePrototype::Mbottom = 4.8; // GeV
/** @brief Initialize static member for mass of charged D meson */
double ParticlePrototype::MchargedD = 1.870; // GeV
/** @brief Initialize static member for mass of neutral D meson */
double ParticlePrototype::MneutralD = 1.865; // GeV
/** @brief Initialize static member for mass of charged B meson */
double ParticlePrototype::MchargedB = 5.280; // GeV
/** @brief Initialize static member for mass of neutral B meson */
double ParticlePrototype::MneutralB = 5.279; // GeV
/** @brief Initialize static member for mass of J/psi */
double ParticlePrototype::Mjpsi = 3.1; // GeV

/** @brief Initialize static member for mass of massA */
double ParticlePrototype::MmassA = 1.0; // GeV - also should be always not 0
/** @brief Initialize static member for mass of massB */
double ParticlePrototype::MmassB = 0.2; // GeV - also should be always not 0
/** @brief Initialize static member for mass of masslessC */
double ParticlePrototype::MmasslessC = 0.0; // GeV - also should be always 0
/** @brief Initialize static member for mass of masslessD */
double ParticlePrototype::MmasslessD = 0.0; // GeV - also should be always 0

/**
* Initialize static member for unique numbering of particles
*/
int ParticlePrototype::unique_id_counter = 0;


/**
 * Return mass of particle with flavor _flav
 *
 * @param[in] _flav flavor of particle 
 * @return mass of particle
 */
double ParticlePrototype::getMass( const FLAVOR_TYPE _flav )
{
  if( _flav <= 6 || _flav == electron || _flav == positron || _flav == photon )
    return 0.0;
  else if( _flav == charm || _flav == anti_charm )
    return Mcharm;
  else if( _flav == bottom || _flav == anti_bottom )
    return Mbottom;
  else if( _flav == jpsi )
    return Mjpsi;
  else if( _flav == dmeson_plus || _flav == dmeson_minus )
    return MchargedD;
  else if( _flav == dmeson_zero || _flav == dmeson_zero_bar )
    return MneutralD;
  else if( _flav == bmeson_plus || _flav == bmeson_minus )
    return MchargedB;
  else if( _flav == bmeson_zero || _flav == bmeson_zero_bar )
    return MneutralB;
  else
    return -1.0;
}



/**
 * Return mass of particle with flavor _flav
 *
 * @param[in] _flav flavor of particle 
 * @return mass of particle
 */
std::string ParticlePrototype::getName( const FLAVOR_TYPE _flav )
{
  std::string name;
  
  switch ( _flav )
  {
    case gluon:
      name = "gluon";
      break;
    case up:
      name = "up";
      break;
    case down:
      name = "down";
      break;
    case strange:
      name = "strange";
      break;
    case charm:
      name = "charm";
      break;
    case bottom:
      name = "bottom";
      break;
    case anti_up:
      name = "anti_up";
      break;
    case anti_down:
      name = "anti_down";
      break;
    case anti_strange:
      name = "anti_strange";
      break;
    case anti_charm:
      name = "anti_charm";
      break;
    case anti_bottom:
      name = "anti_bottom";
      break;
    case photon:
      name = "photon";
      break; 
    case jpsi:
      name = "jpsi";
      break;
    case jpsi_ini:
      name = "jpsi_ini";
      break;
    case jpsi_sec:
      name = "jpsi_sec";
      break;
    case electron:
      name = "electron";
      break;
    case positron:
      name = "positron";
      break;
    case dmeson_plus:
      name = "dmeson_plus";
      break;
    case dmeson_minus:
      name = "dmeson_minus";
      break;
    case dmeson_zero:
      name = "dmeson_zero";
      break;
    case dmeson_zero_bar:
      name = "dmeson_zero_bar";
      break;
    case bmeson_plus:
      name = "bmeson_plus";
      break;
    case bmeson_minus:
      name = "bmeson_minus";
      break;
    case bmeson_zero:
      name = "bmeson_zero";
      break;
    case bmeson_zero_bar:
      name = "bmeson_zero_bar";
      break;
    case quark:
      name = "quark";
      break;     
    case light_quark:
      name = "light_quark";
      break;
    case anti_light_quark:
      name = "anti_light_quark";
      break;
    case allFlavors:
      name = "allFlavors";
      break;
    case dmeson_gen:
      name = "dmeson";
      break;
    case bmeson_gen:
      name = "bmeson";
      break;
    case heavy_quark:
      name = "heavy_quark";
      break;
    case electron_gen:
      name = "all_electron";
      break;
    case c_electron:
      name = "c_electron";
      break;
    case b_electron:
      name = "b_electron";
      break;
    default:
      name = "unknown";
      break;
  }
  

  
  return name;
}


/**
 * Set number of active flavor
 *
 * @param[in] _N number of flavor
 */
void ParticlePrototype::set_N_flavor( const int _N ) 
{ 
  if( _N <= 3 )
  {
    N_light_flavor = _N;
    N_heavy_flavor = 0;
  }
  else if( _N <= 5 )
  {
    N_light_flavor = 3;
    N_heavy_flavor = _N - 3;
  }
  else
  {
    std::cout << "Error! Too many flavor set: " << _N << ". Set maximum number: 5 = 3 + 2." << std::endl;
    N_light_flavor = 3;
    N_heavy_flavor = 2;
  }
};

/**
 * This is a shortcut for
 *   Pos += Mom * ( T - Pos.T )/ Mom.E
 *
 * @param[in] time final time
 */
void ParticlePrototype::Propagate( const double time)
{
  Pos += Mom * (( time - Pos.T() ) / Mom.E());
  Pos.T() = time; // necessary due to rounding errors
};


/**
 * This is a shortcut for
 *   Pos += Mom * ( T - Pos.T )/ Mom.E
 *
 * @param[in] time final time
 * @param[out] distance The given parameter is increased by the
 *   traveled distance
 */
void ParticlePrototype::Propagate( const double time, double & distance)
{
  VectorTXYZ Dist;
  Dist = Mom * (( time - Pos.T() ) / Mom.E());
  Pos += Dist;
  Pos.T() = time; // necessary due to rounding errors
  distance += sqrt( Dist.vec2() );
};

std::ostream& operator<<(std::ostream &os, const ParticlePrototype &obj)
{
  int osval = os.iword(particle_xalloc());
  switch (osval)
  {
    case 0:
    {
      os << obj.unique_id << " " << obj.FLAVOR << " " << obj.m << " " 
         << obj.Mom << " " << obj.Pos;
    } break;
    case 1:
    {
      os << obj.unique_id << " " << obj.FLAVOR << " " << obj.m << " " 
         << obj.Mom << " " << obj.Pos << " " 
         << obj.cell_id << " " << obj.dead << " " 
         << obj.rate22 << " " << obj.rate23 << " " << obj.rate32 << " " ;
    } break;
    default:
    {
      os << "***** wrong value of output style = " << osval << "*****";
    }
    
  }
  return os;
}
