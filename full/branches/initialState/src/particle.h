//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/particle.h $
//$LastChangedDate: 2017-04-06 15:22:13 +0200 (Do, 06. Apr 2017) $
//$LastChangedRevision: 2559 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Declarations for the Particle class
 */



#ifndef PARTICLE_H
#define PARTICLE_H

#include "particleprototype.h"

#include "allocator.h"
#include "bampsvector.h"

/**
* @brief Provides properties of a particle needed in the simulations.
*
* This class encapsulates properties of particles that are needed for the simulation process, such as position and momentum variables.
* It is derived from ParticlePrototype (and extends it) that provides general properties of particles.
*/
class Particle : public ParticlePrototype
{
  public:
    /** @brief Provide standard constructor (for completeness) */
    Particle() : 
      ParticlePrototype(), 
      eta( 0 ), 
      md2g( 0 ), md2q( 0 ), 
      N_EVENT_pp( 0 ), HARD( true ), N_EVENT_AA( 0 ),
      edge( -1 ), coll_id( -1 ),
      free( true ), init( true ),
      collisionTime( 0 ), collisionPartner( -1 ), 
      Old( 0,0,0,0 ), 
      as22( 0 ), as23( 0 ), 
      rate23v( 0 ), rate32v( 0 ), rate22v( 0), 
      cs22( 0 ), cs23( 0 ), 
      lambda_scaled( 0 ), 
      md2g_scaled_22( 0 ),
      md2q_scaled_22( 0 ), 
      md2g_scaled_23( 0 ), 
      md2q_scaled_23( 0 ), 
      step( 0 ), tstep( 0 ), taustep( 0 ) 
{
};
    
    /** @brief space time rapidity \eta */
    double eta;
    
    /** @brief gluonic screening mass currently associated with the particle object, scaled by alpha_s as md2g / alpha_s [GeV^2] */
    double md2g;
    /** @brief quark screening mass currently associated with the particle object, scaled by alpha_s as md2g / alpha_s [GeV^2] */
    double md2q;
    
    /** @brief Pythia event number */
    int N_EVENT_pp;

    /** @brief Pythia hard or soft scattering */
    bool HARD; // true/1 if parton comes from hard scattering

    /**
     * @brief Event number of heavy ion collision to which particle belongs. 
     *
     * Necessary if the number of added particles is much larger than the 
     * number of particles which would be present in a event according
     * to the test particles number of offline particles. 
     * This is only important if one considers scatterings among the 
     * added particles.  
     **/
    int N_EVENT_AA;
    
    /** @brief index of edge cell the particle belongs to, edge = -1 corresponds to no edge cell */
    short int edge;
    
    /** @brief Unique ID of the last collision of the particle */
    long coll_id;
    
    /** @brief Flag for free streaming, true for particles in regions with energy density < config::freezeOutEnergyDensity or for particles outside the grid */
    bool free;
    /** @brief Flag for discerning particles that are still within their initial formation time */
    bool init;
    
    /** @brief collision ordering time [fm] (geometric collisions) */
    double collisionTime;
    /** @brief collision partner (geometric collisions) */
    int collisionPartner;

    /** @brief Momentum and Energy prior to geometric collision, needed for particles in "edge cell"*/
    VectorEPxPyPz Old;

    /** @brief Mean alpha_s for 2->2 interactions associated with this particle, averaged over cell in previous time step */ 
    double as22;
    /** @brief Mean alpha_s for 2->3 interactions associated with this particle, averaged over cell in previous time step */
    double as23;
    /** @brief Previous mean 2->3 rate (in GeV) associated with this particle, averaged over cell in the second to last time step */ 
    double rate23v;
    /** @brief Previous mean 3->2 rate (in GeV) associated with this particle, averaged over cell in the second to last time step */ 
    double rate32v;
    /** @brief Previous mean 2->2 rate (in GeV) associated with this particle, averaged over cell in the second to last time step */ 
    double rate22v;
    /** @brief Mean 2->2 cross section (1/GeV^2) associated with this particle, averaged over cell in previous time step */
    double cs22;
    /** @brief Mean 2->3 cross section (1/GeV^2) associated with this particle, averaged over cell in previous time step */
    double cs23; 
    /** @brief Mean lambda_scaled associated with this particle, averaged over cell in previous time step */
    double lambda_scaled;
    /** @brief Mean md2g (scaled with s) from 2->2 interactions, averaged over cell in previous time step */
    double md2g_scaled_22;
    /** @brief Mean md2q (scaled with s) from 2->2 interactions, averaged over cell in previous time step */
    double md2q_scaled_22;
    /** @brief Mean md2g (scaled with s) from 2->3 interactions, averaged over cell in previous time step */
    double md2g_scaled_23;
    /** @brief Mean md2q (scaled with s) from 2->3 interactions, averaged over cell in previous time step */
    double md2q_scaled_23;
    
    int step,tstep,taustep;//fm
};

#ifdef BAMPS_DECLARE_ALLOCATOR
BAMPS_DECLARE_ALLOCATOR(Particle);
#else
#warning "ALLOCATOR not defined"
#endif

#endif
