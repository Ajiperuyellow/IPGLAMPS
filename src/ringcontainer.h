//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/ringcontainer.h $
//$LastChangedDate: 2018-04-28 19:57:48 +0200 (Sa, 28. Apr 2018) $
//$LastChangedRevision: 2742 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// at revision 912, this is identical to full/branches/vector4D/src/ringcontainer.h

#ifndef RINGCONTAINER_H
#define RINGCONTAINER_H

#include "bampsvector.h"
#include "ratesmanager.h"
#include "particle.h"


/**
 * @brief class to store information of one ring
 **/

/**
 * Calculation of the energy density
 * =================================
 *
 * We want to calculate the energy density in the rest frame of the
 * cell, thus we have to start with the general boost formula of the
 * energy momentum tensor, given in any frame (here: in the lab frame)
 * to the cell rest frame, 
 * \f{eqnarray*}{  
 *   \hat{T}^{\mu\nu} 
 *   = \Lambda^\mu_{\mu'}\Lambda^\nu_{\nu'} T^{\mu'\nu'}\quad,
 * \f}
 * with \f$ \hat{T} \f$ standing for values in the rest frame of the
 * cell and \f$ T \f$ being in the lab frame. The energy
 * density is defined to be the 00-component in the rest frame, 
 * \f{eqnarray*}{
 * \varepsilon = {\hat{T}}^{00}\quad.
 * \f}
 *
 * The most general formulation of the lorentz boost is given by
 * \f{eqnarray*}{
 * \Lambda = \left(\begin{matrix}
 *  \gamma     & -\gamma v_x & -\gamma v_y & -\gamma v_z \\ 
 * -\gamma v_x & 1+(\gamma-1)\frac{v_x^2}{v^2} & (\gamma-1)\frac{v_xv_y}{v^2} & (\gamma-1)\frac{v_xv_z}{v^2} \\
 * -\gamma v_x & (\gamma-1)\frac{v_xv_y}{v^2} & 1+(\gamma-1)\frac{v_y^2}{v^2} & (\gamma-1)\frac{v_yv_z}{v^2} \\
 * -\gamma v_z & (\gamma-1)\frac{v_xv_z}{v^2} & (\gamma-1)\frac{v_yv_z}{v^2} & 1+(\gamma-1)\frac{v_z^2}{v^2} \\
 * \end{matrix}\right) \quad.
 * \f}
 *
 * Thus we have 
 * \f{eqnarray*}{
 * \varepsilon &=&
 * \left(\Lambda^0_0\right)^2T^{00} + \left(\Lambda^0_1\right)^2T^{11}
 * + \left(\Lambda^0_2\right)^2T^{22} + \left(\Lambda^0_3\right)^2T^{33}\\
 * && \quad 
 * +2\Lambda^0_1\Lambda^0_0T^{10} +2\Lambda^0_2\Lambda^0_0T^{20}
 * +2\Lambda^0_3\Lambda^0_1T^{30} +2\Lambda^0_2\Lambda^0_1T^{21}    
 * +2\Lambda^0_3\Lambda^0_1T^{31} +2\Lambda^0_3\Lambda^0_2T^{32}\\
 * &=&
 * \gamma^2T^{00} + \gamma^2v_x^2T^{11}
 * + \gamma^2v_y^2T^{22} + \gamma^2v_z^2T^{33}\\
 * && \quad 
 * -2\gamma^2v_xT^{10} -2\gamma^2v_yT^{20} -2\gamma^2v_zT^{30}
 * +2\gamma^2v_xv_yT^{21} +2\gamma^2v_xv_zT^{31} +2\gamma^2v_yv_zT^{32}
 * \f}
 *
 * Now we are defininig radial coordinates
 * \f$ \vec v=(v_r\cos\theta,v_r\sin\theta,v_z) \f$
 * with the tensor components
 * \f{eqnarray*}{
 * T^{\mu r}&=&T^{\mu 1}\cos\theta+T^{\mu 2}\sin\theta\\
 * T^{rr}&=&T^{11}\cos^2\theta +2T^{12}\sin\theta\cos\theta + T^{22}\sin^2\theta
 * \f}
 * yielding for the energy density
 * \f{eqnarray*}{
 * \varepsilon &=&
 * \gamma^2T^{00} -2\gamma^2v_rT^{0r} -2\gamma^2v_zT^{0z}
 * +\gamma^2 v_r^2T^{rr} +\gamma^2 v_z^2T^{zz}
 * +2\gamma^2v_rv_zT^{zr}\quad.
 * \f}
 *
 * The tensor components per cell (in the lab frame) are defined as
 * \f{eqnarray*}{
 * T^{\mu\nu} &=&\int\frac{d^3p}{(2\pi)^3}\frac{1}{E}\,p^\mu\,p^\nu\,f(x,p)\quad,
 * \f}
 * which are replaced via the MonteCarlo summation as
 * \f{eqnarray*}{
 * T^{\mu\nu} &\longrightarrow&
 * \frac{1}{N\,V}\sum_i^N \frac{1}{E_i}\,p_i^\mu\,p_i^\nu
 * \quad=\quad
 * \frac{1}{V}\left\langle \frac{p^\mu p^\nu}{p^0}\right\rangle
 * \f}
 * and therefore
 * \f{eqnarray*}{
 * T^{00}\longrightarrow\frac{1}{V}\left\langle E\right\rangle\quad,
 * \ T^{0X}\longrightarrow\frac{1}{V}\left\langle p^X\right\rangle\quad,
 * \ T^{XY}\longrightarrow\frac{1}{V}\left\langle \frac{p^Xp^Y}{p^0}\right\rangle\quad.
 * \f}
 **/

class ringContainer
{
public:
  /** 
   * @brief Constructor
   */
  ringContainer();

  /** 
   * @brief Constructor
   */
  ringContainer( const double _minR, const double _maxR );

  /** 
   * @brief Destructor
   */
  ~ringContainer() {};
  
  /**
   * @brief Reset all values
   **/
  void clear();

  /**
   * @brief set inner and outer radius
   **/
  void relocate( const double _minR, const double _maxR );
  
  /**
   * @brief add information of some particle
   **/
  void addParticle( const Particle& _particle );

  /**
   * @brief add information of some particle, which is in formation
   * after a geometrical collision
   **/
  void addParticleInFormGeom( const Particle& _particle, const double _time );

  /**
   * @brief add the rates given by the particle
   **/
  void addRates( const Particle& _particle );
  
  /** @brief return the averaged value **/
  double getAveraged_md2g() const; 

  /** @brief return the averaged value **/
  double getAveraged_md2q() const;
  
  /** @brief return the averaged value **/
  double getAveraged_v_x() const { return ( numberOfParticles > 0 )?( v.X() / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_v_y() const { return ( numberOfParticles > 0 )?( v.Y() / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_v_z() const { return ( numberOfParticles > 0 )?( v.Z() / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  VectorXYZ getAveraged_v() const
  {
    if ( numberOfParticles > 0 )
      return v * (1./ numberOfParticles);
    else
      return VectorXYZ(0,0,0);
  }


  /** @brief return the averaged value **/
  double getAveraged_v_r() const { return ( numberOfParticles > 0 )?( v_r / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getGamma() const;

  /** @brief return the averaged value **/
  double getAccumulated_E() const { return ( numberOfParticles > 0 )?( E ):0; }

  /** @brief return the averaged value **/
  double getAveraged_E() const { return ( numberOfParticles > 0 )?( E / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_p_z() const { return ( numberOfParticles > 0 )?( p_z / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_p_r() const { return ( numberOfParticles > 0 )?( p_r / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_p_t() const { return ( numberOfParticles > 0 )?( p_t / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_pr2_over_E() const { return ( numberOfParticles > 0 )?( pr2_over_E / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_pz2_over_E() const { return ( numberOfParticles > 0 )?( pz2_over_E / numberOfParticles ):0; }

  /** @brief return the averaged value **/
  double getAveraged_pr_pz_over_E() const { return ( numberOfParticles > 0 )?( pr_pz_over_E / numberOfParticles ):0; }
  
  /** @brief return the volume of the ring **/
  double getVolume( const double _dz) const;
  
  /** @brief Do the actual calculation of the averaged values */
  void prepareAverages( const double _dz, const int _Ntest );
  
  /** @brief return the averaged value **/
  double getParticleDensity() const;

  /** @brief return the averaged value **/
  double getGluonDensity() const;

  /** @brief return the averaged value **/
  double getQuarkDensity() const;
  
  /** @brief return the averaged value **/
  double getEnergyDensity() const;

  /** @brief return the averaged value, epsilon/(3*n) **/
  double getEffectiveTemperature() const;
  
  /** @brief return the averaged energy in a comoving frame **/
  double transformEnergyToComovingFrame( VectorEPxPyPz & P ) const;

  /** @brief return the value of #minRadius **/
  double getMinRadius(void) const { return minRadius; };

  /** @brief return the value of #maxRadius **/
  double getMaxRadius(void) const { return maxRadius; };
  
  /** @brief return value of averagesPrepared **/
  bool getAveragesPrepared(void) const { return averagesPrepared; };

  void setValuesAsWorkaround(const double Edens,const double Pdens,const double _gamma, const int _numberOfParticles,const double Gdens,const double Qdens);
  
private:
  double minRadius; ///< the minimal radius
  double maxRadius; ///< the maximal radius
  double deltaR; ///< the thickness of the ring
  
public:
  ratesManager rates; ///< the storage for the rates
    
private:
  bool averagesPrepared; ///< flag that indicates, whether the
                         ///averages are calculated
  
  int numberOfParticles; ///< number of particles
  int numberOfGluons; ///< number of gluons
  int numberOfQuarks; ///< number of quarks
  int numberOfActiveParticles; ///< number of active particles

  double md2g; ///< The averaged value of the gluon Debye mass
  double md2q; ///< The averaged value of the quark Debye mass
  
  VectorXYZ v;
  double v_r;
  
  double E;
  double inverseE_gluons;
  double inverseE_quarks;
  double p_z;
  double p_r;
  double p_t;
  double pr2_over_E;
  double pz2_over_E;
  double pr_pz_over_E;
  
  double gamma;
  double energyDensity;
  double particleDensity;
  double gluonDensity;
  double quarkDensity;
  
  double volume;

  int numberOfCollectedRateObjects;
};



/** @brief exception class */
class eRingContainer_error : public std::runtime_error
{
public:
  explicit eRingContainer_error(const std::string& what) : std::runtime_error(what) {};
  
  virtual ~eRingContainer_error() throw() {};
};


#endif // RINGCONTAINER_H
