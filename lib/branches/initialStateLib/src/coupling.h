//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/coupling.h $
//$LastChangedDate: 2014-07-07 15:55:59 +0200 (Mo, 07. Jul 2014) $
//$LastChangedRevision: 1784 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef COUPLING
#define COUPLING


/** @brief Prototype for the computation of binary cross sections and the sampling of momentum transfers */
class coupling
{
  public:
    
    /** @brief returns value of constant coupling */
    static double get_constant_coupling() { return fixedCoupling_global; };
    
    /**
    * @brief Get the running coupling
    * @param[in] Q2 scale for running coupling, e.g. momentum transfer of process
    * @return coupling strength at this scale
    */
    static double get_running_coupling( const double Q2 );
    
    /**
    * @brief returns coupling strength for given scale
    * @param[in] Q2 scale for running coupling, e.g. momentum transfer of process
    * @param[in] isRunning whether the coupling should run or be constant
    * @return coupling strength
    */
    static double get_coupling( const double Q2, const bool isRunning );
    
    /**
    * @brief returns coupling strength for given scale
    * @param[in] Q2 scale for running coupling, e.g. momentum transfer of process
    * @return coupling strength
    */
    static double get_coupling( const double Q2 );
    
    /**
    * @brief returns the maximum value of the coupling
    * @return maximum value of the coupling
    */
    static double get_maximum_coupling() 
    { 
      if( isRunning_global )
        return alpha_max_global;
      else
        return get_constant_coupling();
    };
    
    /**
    * @brief Set number of active flavors for running coupling determination. Usually number of light quarks
    */
    static void set_Nflavor( const int _N ) { Nflavor = _N; };
    
    /**
    * @brief Set whether coupling is running or not.
    */
    static void set_isRunning( const bool _x ) { isRunning_global = _x; };
    
    /**
    * @brief returns whether coupling is running or not.
    */
    static bool isRunning() { return isRunning_global; };

    /**
    * @brief Set value of fixed coupling.
    */
    static void set_fixedCoupling( const double _x ) { fixedCoupling_global = _x; };
    
    /**
    * @brief Set maximum value of running coupling.
    */
    static void set_maxRunningCoupling( const double _x ) { alpha_max_global = _x; };
    
  private:
    /** @brief Number of active flavor. Usually only light flavors are counted as active in running coupling determination. See for instance Bethke et al. */
    static int Nflavor;
    /** @brief Whether coupling is running or not. Set it in the beginning. get_coupling( const double Q2 ) returns then value for constant or running coupling according to this variable. */
    static bool isRunning_global;
    /** @brief Value of fixed coupling constant. */
    static double fixedCoupling_global;
    /** @brief Maximum value of alpha_s */
    static double alpha_max_global;
};

#endif
