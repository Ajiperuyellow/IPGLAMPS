//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/ratesmanager.cpp $
//$LastChangedDate: 2016-09-14 12:41:00 +0200 (Mi, 14. Sep 2016) $
//$LastChangedRevision: 2410 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include "ratesmanager.h"

std::vector< std::vector< double > > ratesManager::faks22;
std::vector< std::vector< double > > ratesManager::faks23;
std::vector< std::vector< double > > ratesManager::faks32;

bool ratesManager::faksInitialized = false;

void ratesManager::initializeFaks(void)
{
  faks22.assign( interactionType::indexProcessesInclusive22.size(),
                 std::vector<double>( 7, 0 ) );
  faks23.assign( interactionType::indexProcessesInclusive23.size(),
                 std::vector<double>( 7, 0 ) );
  faks32.assign( interactionType::indexProcessesInclusive32.size(),
                 std::vector<double>( 7, 0 ) );

  const std::vector<FLAVOR_TYPE> flav = 
    { gluon, light_quark, anti_light_quark,
      charm, anti_charm, bottom, anti_bottom };

  for ( unsigned int i = 0; i < rates22.size(); i++ )
  {
    int type = interactionType::getInclusiveProcessTypeFromIndex( i, c22 );
    for ( unsigned int j = 0; j < flav.size(); j++ )
    {
      faks22[i][j] = interactionType::getInvolvedInitialInclusive( flav[j], type );
    }
  }

  for ( unsigned int i = 0; i < rates23.size(); i++ )
  {
    int type = interactionType::getInclusiveProcessTypeFromIndex( i, c23 );
    for ( unsigned int j = 0; j < flav.size(); j++ )
    {
      faks23[i][j] = interactionType::getInvolvedInitialInclusive( flav[j], type );
    }
  }

  for ( unsigned int i = 0; i < rates32.size(); i++ )
  {
    int type = interactionType::getInclusiveProcessTypeFromIndex( i, c32 );
    for ( unsigned int j = 0; j < flav.size(); j++ )
    {
      faks32[i][j] = interactionType::getInvolvedInitialInclusive( flav[j], type );
    }
  }

  faksInitialized = true;
};


/**
 * Interface routine to store (collect) a probability for a given 2->X
 * interaction process 
 *
 * @param[in] _type Interaction process type id (221, 3280, etc., see class interactionType)
 * @param[in] _P Probability to be stored
 * @return 1 for successful execution
 */
int ratesManager::add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const double _P )
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::add called after having normalized the rates via ratesManager::normalizeRates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  
  // interactions with zero probability are usually not allowed and
  // have therefore often no process type in class interactionType
  // which will cause an error in getting the index. Thus abort the
  // function call here 
  if( _P == 0.0 )
    return 0;

  int index = interactionType::getIndexFromProcessType( interactionType::getInclusiveProcessType( _F1, _F2, genericType ), genericType );

  switch ( genericType )
  {
    case c22:
      rates22[index] += _P;
      ++simulatedNumber22[index];
      break;
    case c23:
      rates23[index] += _P;
      ++simulatedNumber23[index];
      break;
    default:
    {
      std::string errMsg = "Generic interaction type not found in ratesManager::add. Unrecoverable error.";
      throw eRatesManager_error( errMsg );
    }
  }

  return 1;
}




/**
* Interface routine to store (collect) a probability for a given 3->X interaction process
*
* @param[in] _type Interaction process type id (221, 3280, etc., see class interactionType)
* @param[in] _P Probability to be stored
* @return 1 for successful execution
*/
int ratesManager::add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, const double _P )
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::add called after having normalized the rates via ratesManager::normalizeRates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  
  // interactions with zero probability are usually not allowed and
  // have therefore often no process type in class interactionType
  // which will cause an error in getting the index. Thus abort the
  // function call here 
  if( _P == 0.0 )
    return 0;

  int index = interactionType::getIndexFromProcessType( interactionType::getInclusiveProcessType( _F1, _F2, _F3, genericType ), genericType );

  if ( genericType == c32 )
  {
    rates32[index] += _P;
    ++simulatedNumber32[index];
  }
  else
  {
    std::string errMsg = "Generic interaction type not found in ratesManager::add (3->2). Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }

  return 1;
}






/**
 * Normalize the accumulated probabilties such that one gets the
 * interaction rate per single particle and per time step dt. 
 * For this the accumulated probabilities for a given process is
 * weighted with the number of gluons (quarks, anti-quarks) involved
 * in the initial state of this interaction. 
 * e.g. the rate for a single gluons in ggg -> gg processes would be:
 *   3 * sum(probabilities) / (Ng * dt) 
 * See notes dating 27.05.2010 for more details.
 *
 * @param[in] Ng Number of gluons for which the interaction probabilities have been accumulated
 * @param[in] Nq Number of light quarks for which the interaction probabilities have been accumulated
 * @param[in] Nqbar Number of light anti-quarks for which the interaction probabilities have been accumulated
 * @param[in] Nc Number of charm quarks for which the interaction probabilities have been accumulated
 * @param[in] Ncbar Number of charm anti-quarks for which the interaction probabilities have been accumulated
 * @param[in] Nb Number of bottom quarks for which the interaction probabilities have been accumulated
 * @param[in] Nbbar Number of bottom anti-quarks for which the interaction probabilities have been accumulated
 * @param[in] dt Time step dt in fm/c
 * @param[in] dV Cell volume dV in fm^3/c
 * @return 1 for successful execution
 */
int ratesManager::normalizeRates( const int Ng, const int Nq, const int Nqbar, const int Nc, const int Ncbar, const int Nb, const int Nbbar, const double dt, const double dV )
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::normalizeRates called twice. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    ratesFlav22.assign( 7, 0 );
    ratesFlav23.assign( 7, 0 );
    ratesFlav32.assign( 7, 0 );
    
    for ( unsigned int i = 0; i < rates22.size(); i++ )
    {
      for ( unsigned int j = 0; j < ratesFlav22.size(); j++)
      {
        ratesFlav22[j] += faks22[i][j] * rates22[i];
      }

      rates22[i] /= dt * dV;  // this yields the rate for the given process per volume - not per particle as for the entries in ratesFlav etc.
    }

    for ( unsigned int i = 0; i < rates23.size(); i++ )
    {
      for ( unsigned int j = 0; j < ratesFlav23.size(); j++)
      {
        ratesFlav23[j] += faks23[i][j] * rates23[i];
      }

      rates23[i] /= dt * dV;  // this yields the rate for the given process per volume - not per particle as for the entries in ratesFlav etc.
    }

    for ( unsigned int i = 0; i < rates32.size(); i++ )
    {
      for ( unsigned int j = 0; j < ratesFlav32.size(); j++)
      {
        ratesFlav32[j] += faks32[i][j] * rates32[i];
      }

      rates32[i] /= dt * dV;  // this yields the rate for the given process per volume - not per particle as for the entries in ratesFlav etc.
    }


    // it can happen that initially one particle species (say light
    // quarks) were not present in the cell for which the rates are
    // calculated. However, it is possible that light quarks are
    // produced in this timestep and then accumulate probabilities
    // and hence rates. Avoid here deviding by zero for this case
    // and also set the rates to zero since they cannot be properly
    // normalized. Anyways, this should be a small effect. 
    
    const std::vector<int> givenN = { Ng, Nq, Nqbar, Nc, Ncbar, Nb, Nbbar };

    for (unsigned int i=0; i< givenN.size(); i++)
    {
      double fak = ( givenN[i] != 0 ) ? 1.0/( dt * givenN[i] ) : 0.0;
      ratesFlav22[i] *= fak;
      ratesFlav23[i] *= fak;
      ratesFlav32[i] *= fak;
    }

    for (unsigned int i=0; i < givenN.size(); i++)
    {
      ratesFlavTotal[i] = ratesFlav22[i] + ratesFlav23[i] + ratesFlav32[i];
    }

    isNormalized = true;
    return 1;
  }
}



int ratesManager::normalizeRates()
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::normalizeRates called twice. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    ratesFlav22.assign( 7, 0 );
    ratesFlav23.assign( 7, 0 );
    ratesFlav32.assign( 7, 0 );
    ratesFlavTotal.assign( 7, 0 );

    isNormalized = true;
    return 1;
  }
}



/**
 * Overloaded function for energyloss, create_mfp_tables,...
 * Normalize the accumulated probabilties such that one gets the
 * interaction rate per single particle and per time step dt. 
 * This routine is used in situations where no actual probabilities were stored, but only sigma*vrel. This is the case for energyloss, create_mfp_tables, ...
 * Replaces routines from "ratesManager_mfp"
 *
 * @param[in] ng gluon density
 * @param[in] nq quark density (summed over all flavors!)
 * @param[in] nqbar anti-quark density (summed over all flavors!)
 * @return 1 for successful execution
 */
int ratesManager::normalizeRates( const double ng, const double nq, const double nqbar )
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::normalizeRates called twice. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    ratesFlav22.assign( 7, 0 );
    ratesFlav23.assign( 7, 0 );
    ratesFlav32.assign( 7, 0 );

    std::vector<double> densFaks22;
    std::vector<double> densFaks23;
    std::vector<double> densFaks32;

    double inclusiveFaks22 = 0.0;
    double inclusiveFaks23 = 0.0;
    double inclusiveFaks32 = 0.0;

    int type = 0;

    for ( int i = 0; i < rates22.size(); i++ )
    {
      densFaks22.assign( 7, 0 );
      
      type = interactionType::getInclusiveProcessTypeFromIndex( i, c22 );

      switch ( type )
      {
        // light partons
        case ( 9221 ):      // 9221 : g+g -> X
          densFaks22[0] = ng;
          densFaks22[1] = 0;
          densFaks22[2] = 0;
          inclusiveFaks22 = pow( ng, 2.0 ) / 2.0;
          break;
        case ( 9222 ):      // 9222 : g+q -> X
          densFaks22[0] = nq;
          densFaks22[1] = ng;
          densFaks22[2] = 0;
          inclusiveFaks22 = ng * nq;
          break;
        case ( 9223 ):      // 9223 : g+qbar -> X
          densFaks22[0] = nqbar;
          densFaks22[1] = 0;
          densFaks22[2] = ng;
          inclusiveFaks22 = ng * nqbar;
          break;
        case ( 9224 ):      // 9224 : q+qbar -> X
          densFaks22[0] = 0;
          densFaks22[1] = nqbar / 3.0;
          densFaks22[2] = nq / 3.0;
          inclusiveFaks22 = 1.0 / 3.0 * nq * nqbar;
          break;
        case ( 9225 ):      // 9225 : q+q -> X
          densFaks22[0] = 0;
          densFaks22[1] = nq / 3.0;
          densFaks22[2] = 0;
          inclusiveFaks22 = 1.0 / 6.0 * pow( nq, 2.0 );
          break;
        case ( 9226 ):      // 9226 : qbar+qbar -> X
          densFaks22[0] = 0;
          densFaks22[1] = 0;
          densFaks22[2] = nqbar / 3.0;
          inclusiveFaks22 = 1.0 / 6.0 * pow( nqbar, 2.0 );
          break;
        case ( 9227 ):      // 9227 : q+q' -> X
          densFaks22[0] = 0;
          densFaks22[1] = 2.0 / 3.0 * nq;
          densFaks22[2] = 0;
          inclusiveFaks22 = 1.0 / 3.0 * pow( nq, 2.0 );
          break;
        case ( 9228 ):      // 9228 : q+qbar' -> X
          densFaks22[0] = 0;
          densFaks22[1] = 2.0 / 3.0 * nqbar;
          densFaks22[2] = 2.0 / 3.0 * nq;
          inclusiveFaks22 = 2.0 / 3.0 * nq * nqbar;
          break;
        case ( 9229 ):      // 9229 : qbar+qbar' -> X
          densFaks22[0] = 0;
          densFaks22[1] = 0;
          densFaks22[2] = 2.0 / 3.0 * nqbar;
          inclusiveFaks22 = 1.0 / 3.0 * pow( nqbar, 2.0 );
          break;
        
        // heavy quarks
        case ( 92210 ):      // g+c -> X
          densFaks22[0] = 0; // negligible since proportional to n_c
          densFaks22[3] = ng;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92211 ):      // g+cbar -> X
          densFaks22[0] = 0; // negligible since proportional to n_c
          densFaks22[4] = ng;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92212 ):      // q+c -> X
          densFaks22[1] = 0; // negligible since proportional to n_c
          densFaks22[3] = nq;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92213 ):      // q+cbar -> X
          densFaks22[1] = 0; // negligible since proportional to n_c
          densFaks22[4] = nq;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92214 ):      // qbar+c -> X
          densFaks22[2] = 0; // negligible since proportional to n_c
          densFaks22[3] = nqbar;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92215 ):      // qbar+cbar -> X
          densFaks22[2] = 0; // negligible since proportional to n_c
          densFaks22[4] = nqbar;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92216 ):      // c+cbar -> X
          densFaks22[3] = 0; // negligible since proportional to n_c
          densFaks22[4] = 0; // negligible since proportional to n_c
          inclusiveFaks22 = 0.0;
          break;
          
        case ( 92220 ):      // g+b -> X
          densFaks22[0] = 0; // negligible since proportional to n_b
          densFaks22[5] = ng;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92221 ):      // g+bbar -> X
          densFaks22[0] = 0; // negligible since proportional to n_b
          densFaks22[6] = ng;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92222 ):      // q+b -> X
          densFaks22[1] = 0; // negligible since proportional to n_b
          densFaks22[5] = nq;
          inclusiveFaks22 = 0.0;
         break;
        case ( 92223 ):      // q+bbar -> X
          densFaks22[1] = 0; // negligible since proportional to n_b
          densFaks22[6] = nq;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92224 ):      // qbar+b -> X
          densFaks22[2] = 0; // negligible since proportional to n_b
          densFaks22[5] = nqbar;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92225 ):      // qbar+bbar -> X
          densFaks22[2] = 0; // negligible since proportional to n_b
          densFaks22[6] = nqbar;
          inclusiveFaks22 = 0.0;
          break;
        case ( 92226 ):      // b+bbar -> X
          densFaks22[3] = 0; // negligible since proportional to n_b
          densFaks22[4] = 0; // negligible since proportional to n_b
          inclusiveFaks22 = 0.0;
          break;
        default:
//           std::string errMsg = "ratesManager_mfp::normalizeRates 22 type not found";
//           throw eRatesManager_error( errMsg );
          break;
      }
      
      if ( simulatedNumber22[i] > 0 )
      {
        for( int j = 0; j < 7; j++ )
        {
          ratesFlav22[j] += pow( 0.197, 2.0 ) * densFaks22[j] * rates22[i] / simulatedNumber22[i];
        }
        
        rates22[i] = inclusiveFaks22 * pow( 0.197, 2.0 ) * rates22[i] / simulatedNumber22[i];
      }
    }

    for ( int i = 0; i < rates23.size(); i++ )
    {
      densFaks23.assign( 7, 0 );
      
      type = interactionType::getInclusiveProcessTypeFromIndex( i, c23 );

      switch ( type )
      {
        // light partons
        case ( 9231 ):      // 9231 : g+g -> X
          densFaks23[0] = ng;
          densFaks23[1] = 0;
          densFaks23[2] = 0;
          inclusiveFaks23 = pow( ng, 2.0 ) / 2.0;
          break;
        case ( 9232 ):      // 9232 : g+q -> X
          densFaks23[0] = nq;
          densFaks23[1] = ng;
          densFaks23[2] = 0;
          inclusiveFaks23 = ng * nq;
          break;
        case ( 9233 ):      // 9233 : g+qbar -> X
          densFaks23[0] = nqbar;
          densFaks23[1] = 0;
          densFaks23[2] = ng;
          inclusiveFaks23 = ng * nqbar;
          break;
        case ( 9234 ):      // 9234 : q+qbar -> X
          densFaks23[0] = 0;
          densFaks23[1] = nqbar / 3.0;
          densFaks23[2] = nq / 3.0;
          inclusiveFaks23 = 1.0 / 3.0 * nq * nqbar;
          break;
        case ( 9235 ):      // 9235 : q+q -> X
          densFaks23[0] = 0;
          densFaks23[1] = nq / 3.0;
          densFaks23[2] = 0;
          inclusiveFaks23 = 1.0 / 6.0 * pow( nq, 2.0 );
          break;
        case ( 9236 ):      // 9236 : qbar+qbar -> X
          densFaks23[0] = 0;
          densFaks23[1] = 0;
          densFaks23[2] = nqbar / 3.0;
          inclusiveFaks23 = 1.0 / 6.0 * pow( nqbar, 2.0 );
          break;
        case ( 9237 ):      // 9237 : q+q' -> X
          densFaks23[0] = 0;
          densFaks23[1] = 2.0 / 3.0 * nq;
          densFaks23[2] = 0;
          inclusiveFaks23 = 1.0 / 3.0 * pow( nq, 2.0 );
          break;
        case ( 9238 ):      // 9238 : q+qbar' -> X
          densFaks23[0] = 0;
          densFaks23[1] = 2.0 / 3.0 * nqbar;
          densFaks23[2] = 2.0 / 3.0 * nq;
          inclusiveFaks23 = 2.0 / 3.0 * nq * nqbar;
          break;
        case ( 9239 ):      // 9239 : qbar+qbar' -> X
          densFaks23[0] = 0;
          densFaks23[1] = 0;
          densFaks23[2] = 2.0 / 3.0 * nqbar;
          inclusiveFaks23 = 1.0 / 3.0 * pow( nqbar, 2.0 );
          break;
        
        // heavy quarks
        case ( 92310 ):      // g+c -> X
          densFaks23[0] = 0; // negligible since proportional to n_c
          densFaks23[3] = ng;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92311 ):      // g+cbar -> X
          densFaks23[0] = 0; // negligible since proportional to n_c
          densFaks23[4] = ng;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92312 ):      // q+c -> X
          densFaks23[1] = 0; // negligible since proportional to n_c
          densFaks23[3] = nq;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92313 ):      // q+cbar -> X
          densFaks23[1] = 0; // negligible since proportional to n_c
          densFaks23[4] = nq;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92314 ):      // qbar+c -> X
          densFaks23[2] = 0; // negligible since proportional to n_c
          densFaks23[3] = nqbar;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92315 ):      // qbar+cbar -> X
          densFaks23[2] = 0; // negligible since proportional to n_c
          densFaks23[4] = nqbar;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92316 ):      // c+cbar -> X
          densFaks23[3] = 0; // negligible since proportional to n_c
          densFaks23[4] = 0; // negligible since proportional to n_c
          inclusiveFaks23 = 0.0;
          break;
          
        case ( 92320 ):      // g+b -> X
          densFaks23[0] = 0; // negligible since proportional to n_b
          densFaks23[5] = ng;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92321 ):      // g+bbar -> X
          densFaks23[0] = 0; // negligible since proportional to n_b
          densFaks23[6] = ng;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92322 ):      // q+b -> X
          densFaks23[1] = 0; // negligible since proportional to n_b
          densFaks23[5] = nq;
          inclusiveFaks23 = 0.0;
         break;
        case ( 92323 ):      // q+bbar -> X
          densFaks23[1] = 0; // negligible since proportional to n_b
          densFaks23[6] = nq;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92324 ):      // qbar+b -> X
          densFaks23[2] = 0; // negligible since proportional to n_b
          densFaks23[5] = nqbar;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92325 ):      // qbar+bbar -> X
          densFaks23[2] = 0; // negligible since proportional to n_b
          densFaks23[6] = nqbar;
          inclusiveFaks23 = 0.0;
          break;
        case ( 92326 ):      // b+bbar -> X
          densFaks23[5] = 0; // negligible since proportional to n_b
          densFaks23[6] = 0; // negligible since proportional to n_b
          inclusiveFaks23 = 0.0;
          break;
        default:
//           std::string errMsg = "ratesManager_mfp::normalizeRates 23 type not found";
//           throw eRatesManager_error( errMsg );
          break;
      }
      
      if ( simulatedNumber23[i] > 0 )
      {
        for( int j = 0; j < 7; j++ )
        {
          ratesFlav23[j] += pow( 0.197, 2.0 ) * densFaks23[j] * rates23[i] / simulatedNumber23[i];
        }
        
        rates23[i] = inclusiveFaks23 * pow( 0.197, 2.0 ) * rates23[i] / simulatedNumber23[i];
      }
    }


    for ( int i = 0; i < rates32.size(); i++ )
    {
      densFaks32.assign( 7, 0 );
      
      type = interactionType::getInclusiveProcessTypeFromIndex( i, c32 );
      
      switch ( type )
      {
        case ( 9321 ):      // 9321 : g+g+g -> X
          densFaks32[0] = pow( ng, 2.0 ) / 2.0;
          inclusiveFaks32 = pow( ng, 3.0 ) / 6.0;
          break;
        case ( 9322 ):      // 9322 : g+q+g -> X
          densFaks32[0] = nq * ng;
          densFaks32[1] = pow( ng, 2.0 ) / 2.0;
          inclusiveFaks32 = pow( ng, 2.0 ) * nq / 2.0;
          break;
        case ( 9323 ):      // 9323 : g+qbar+g -> X
          densFaks32[0] = nqbar * ng;
          densFaks32[2] = pow( ng, 2.0 ) / 2.0;
          inclusiveFaks32 = pow( ng, 2.0 ) * nqbar / 2.0;
          break;
        case ( 9324 ):      // 9324 : q+qbar+g -> X
          densFaks32[0] = nq * nqbar / 3.0;
          densFaks32[1] = nqbar * ng / 3.0;
          densFaks32[2] = nq * ng / 3.0;
          inclusiveFaks32 = ng * nq * nqbar / 3.0;
          break;
        case ( 9325 ):      // 9325 : q+q+g -> X
          densFaks32[0] = pow( nq, 2.0 ) / 6.0; // this differs from previous revisions, cf. notes from 13.09.2016.
          densFaks32[1] = ng * nq / 3.0;
          inclusiveFaks32 = ng * pow( nq, 2.0 ) / 6.0;
          break;
        case ( 9326 ):      // 9326 : qbar+qbar+g -> X
          densFaks32[0] = pow( nqbar, 2.0 ) / 6.0; // this differs from previous revisions, cf. notes from 13.09.2016.
          densFaks32[2] = ng * nqbar / 3.0;
          inclusiveFaks32 = ng * pow( nqbar, 2.0 ) / 6.0;
          break;
        case ( 9327 ):      // 9327 : q+q'+g -> X
          densFaks32[0] = 1.0 / 3.0 * pow( nq, 2.0 );
          densFaks32[1] = 2.0 / 3.0 * nq * ng;
          inclusiveFaks32 = ng * pow( nq, 2.0 ) * 1.0 / 3.0;
          break;
        case ( 9328 ):      // 9328 : q+qbar'+g -> X
          densFaks32[0] = 2.0 / 3.0 * nq * nqbar; // remark: this was different in Oli's thesis (p.163), compare notes from 13.09.2016. However in previous revisions the right formular was used.
          densFaks32[1] = 2.0 / 3.0 * nqbar * ng;
          densFaks32[2] = 2.0 / 3.0 * nq * ng;
          inclusiveFaks32 = ng * nq * nqbar * 2.0 / 3.0;
          break;
        case ( 9329 ):      // 9329 : qbar+qbar'+g -> X
          densFaks32[0] = 1.0 / 3.0 * pow( nqbar, 2 );
          densFaks32[2] = 2.0 / 3.0 * nqbar * ng;
          inclusiveFaks32 = ng * pow( nqbar, 2.0 ) * 1.0 / 3.0;
          break;
        default:
//           std::string errMsg = "ratesManager_mfp::normalizeRates 32 type not found";
//           throw eRatesManager_error( errMsg );
          break;
      }
      
      std::vector< double > averagedRate;
      averagedRate.assign( 7, 0 );
      
      if ( simulatedNumber32[i] > 0 )
      {
        for( int j = 0; j < 7; j++ )
        {
          ratesFlav32[j] += pow( 0.197, 5.0 ) * densFaks32[j] * rates32[i] / simulatedNumber32[i];
        }
        
        rates32[i] = inclusiveFaks32 * pow( 0.197, 5.0 ) * rates32[i] / simulatedNumber32[i];
      }
    }

    for( int j = 0; j < 7; j++ )
    {
      ratesFlavTotal[j] = ratesFlav22[j] + ratesFlav23[j] + ratesFlav32[j];
    }

    isNormalized = true;

    return 1;
  }
}



/**
 * Get the rate for a certain particle type
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, light_quark, anti_light_quark, charm, anit_charm, bottom, anti_bottom)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getRate( const FLAVOR_TYPE particleType, const UNIT_TYPE unit ) const
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::getRate called without prior normalization of rates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double conversionFactor = 1;
    if ( unit == GeV )
    {
      conversionFactor = 0.197;
    }

    unsigned int index = flavIndex( ParticlePrototype::mapToGenericFlavorType( particleType ) );
    return ratesFlavTotal[index] * conversionFactor;
  }
}


/**
 * Get the rate for a certain particle type and a certain generic collision type (22, 23, 32)
 * e.g. to get the rate of quarks in 2->2 collisions, call ratesManager::getRate( quark, c22, fm )
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, light_quark, anti_light_quark, charm, anit_charm, bottom, anti_bottom)
 * @param[in] collType Generic collision type for which the rate is requested (c22, c23, c32)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getRate( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit ) const
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::getRate called without prior normalization of rates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double conversionFactor = 1;
    if ( unit == GeV )
    {
      conversionFactor = 0.197;
    }

    unsigned int index = flavIndex( ParticlePrototype::mapToGenericFlavorType( particleType ) );

    switch( collType )
    {
      case c22:
        return ratesFlav22[index] * conversionFactor;
      case c23:
        return ratesFlav23[index] * conversionFactor;
      case c32:
        return ratesFlav32[index] * conversionFactor;
    }
  }
  return 0.0;
}



/**
 * Get the rate for a certain inclusive process (9231 etc.).
 * NB: this rate is NOT per particle as all other rates but per process
 *
 * @param[in] inclCollType Inclusive collision type (9221, 9326 etc)
 * @param[in] genCollType Generic collision type for which the rate is requested (c22, c23, c32)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return Rate for the inclusive process. NOT per particle but per gq -> X etc..
 */
double ratesManager::getRateInclusive( const int inclCollType, const GENERIC_COLL_TYPE genCollType, const UNIT_TYPE unit ) const
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::getRateInclusive called without prior normalization of rates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else if( interactionType::getInvolvedInitialInclusive( charm, inclCollType ) != 0 || interactionType::getInvolvedInitialInclusive( anti_charm, inclCollType ) != 0 
    || interactionType::getInvolvedInitialInclusive( bottom, inclCollType ) != 0 || interactionType::getInvolvedInitialInclusive( anti_bottom, inclCollType ) != 0 )
  {
    std::string errMsg = "ratesManager::getRateInclusive called for inclusive process where heavy quark is present. This is not supported at the moment. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double conversionFactor = 1;
    if ( unit == GeV )
    {
      conversionFactor = 0.197;
    }

    int index = interactionType::getIndexFromProcessType( inclCollType, genCollType );

    switch ( genCollType )
    {
      case c22:
        return rates22[index] * conversionFactor;
      case c23:
        return rates23[index] * conversionFactor;
      case c32:
        return rates32[index] * conversionFactor;
      default:
        std::string errMsg = "Generic type not found in ratesManager::getRateInclusive. Unrecoverable error.";
        throw eRatesManager_error( errMsg );
        break;
    }
  }
}

/**
 * Get the mean free path for a certain particle type
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, light_quark, anti_light_quark, charm, anit_charm, bottom, anti_bottom)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getLambda( const FLAVOR_TYPE particleType, const UNIT_TYPE unit ) const
{
  double lambda_out;
  
  if( particleType > 2*ParticlePrototype::max_N_light_flavor )
  {
    std::cout << "For heavy particles with flavor " << particleType << " the velocity of the particle must be specified." << std::endl;
    std::string errMsg = "For heavy particles the velocity of the particle must be specified.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double velocity = 1.0;
    lambda_out = getLambda( particleType, velocity, unit );
  }
  
  return lambda_out;
}

/**
 * Get the mean free path for a certain particle type with specifying the velocity
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, quark, anti-quark)
 * @param[in] velocity Velocity of the particle
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getLambda( const FLAVOR_TYPE particleType, const double velocity, const UNIT_TYPE unit ) const
{
  double R = getRate( particleType, unit );

  if ( R > 0 )
  {
    return ( velocity / R );
  }
  else
  {
    return 0;
  }
}


/**
 * Get the mean free path for a certain particle type and a certain generic collision type (22, 23, 32)
 * e.g. to get the mean free path of quarks in 2->2 collisions, call ratesManager::getRate( quark, c22, fm )
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, quark, anti-quark)
 * @param[in] collType Generic collision type for which the rate is requested (c22, c23, c32)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getLambda( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit ) const
{
  double lambda_out;
  
  if( particleType > 2*ParticlePrototype::max_N_light_flavor )
  {
    std::cout << "For heavy particles with flavor " << particleType << " the velocity of the particle must be specified." << std::endl;
    std::string errMsg = "For heavy particles the velocity of the particle must be specified.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double velocity = 1.0;
    lambda_out = getLambda( particleType, collType, velocity, unit );
  }
  
  return lambda_out;
}


/**
 * Get the mean free path for a certain particle type and a certain generic collision type (22, 23, 32)
 * e.g. to get the mean free path of quarks in 2->2 collisions, call ratesManager::getRate( quark, c22, fm )
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, quark, anti-quark)
 * @param[in] collType Generic collision type for which the rate is requested (c22, c23, c32)
 * @param[in] velocity Velocity of the particle
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getLambda( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const double velocity, const UNIT_TYPE unit ) const
{
  double R = getRate( particleType, collType, unit );
  return ( R > 0 ) ? ( velocity / R ) : 0;
}



void ratesManager::clear()
{
  isNormalized = false;
  rates22.assign( interactionType::indexProcessesInclusive22.size(), 0 );
  rates23.assign( interactionType::indexProcessesInclusive23.size(), 0 );
  rates32.assign( interactionType::indexProcessesInclusive32.size(), 0 );
  simulatedNumber22.assign( interactionType::indexProcessesInclusive22.size(), 0 );
  simulatedNumber23.assign( interactionType::indexProcessesInclusive23.size(), 0 );
  simulatedNumber32.assign( interactionType::indexProcessesInclusive32.size(), 0 );
  
  ratesFlav22.assign( 7, 0 );
  ratesFlav23.assign( 7, 0 );
  ratesFlav32.assign( 7, 0 );
  ratesFlavTotal.assign( 7, 0 );
  nCollectedParticleBasedRates.assign( 7, 0 );
}



void ratesManager::addParticleBasedRates( const ParticlePrototype& _particle, const UNIT_TYPE unit )
{
  double conversionFactor = 1;
  if ( unit == GeV )
  {
    conversionFactor = 1 / 0.197;
  }
  
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::addParticleBasedRates only works for previously finalized instances";
    throw eRatesManager_error( errMsg );
  }

  unsigned int index = flavIndex( ParticlePrototype::mapToGenericFlavorType( _particle.FLAVOR ) );

  ratesFlav22[index] += _particle.rate22 * conversionFactor;
  ratesFlav23[index] += _particle.rate23 * conversionFactor;
  ratesFlav32[index] += _particle.rate32 * conversionFactor;
  ratesFlavTotal[index] += ( _particle.rate22 + _particle.rate23 + _particle.rate32 ) * conversionFactor;
  ++nCollectedParticleBasedRates[index];
}


void ratesManager::prepareParticleBasedAverages()
{
  for (unsigned int i = 0; i < 7; i++)
  {
    double fak = ( nCollectedParticleBasedRates[i] > 0 ) ? 1.0/nCollectedParticleBasedRates[i] : 0.0;
    ratesFlav22[i] *= fak;
    ratesFlav23[i] *= fak;
    ratesFlav32[i] *= fak;
    ratesFlavTotal[i] *= fak;
  }
}




ratesManager& ratesManager::operator+=( const ratesManager & rhs )
{
  if ( !( rhs.isNormalized && ( *this ).isNormalized ) )
  {
    std::string errMsg = "ratesManager::operator+= only works for previously finalized instances";
    throw eRatesManager_error( errMsg );
  }

  // four independent loops would probably be faster (cache lines!)
  for ( int j = 0; j < 7; j++ )
  {
    ratesFlav22[j] += rhs.ratesFlav22[j];
    ratesFlav23[j] += rhs.ratesFlav23[j];
    ratesFlav32[j] += rhs.ratesFlav32[j];
    ratesFlavTotal[j] += rhs.ratesFlavTotal[j];
  }

  for ( int j = 0; j < rates22.size(); j++ )
  {
    rates22[j] += rhs.rates22[j];
  }

  for ( int j = 0; j < rates23.size(); j++ )
  {
    rates23[j] += rhs.rates23[j];
  }

  for ( int j = 0; j < rates32.size(); j++ )
  {
    rates32[j] += rhs.rates32[j];
  }
  
  return ( *this );
}




ratesManager& ratesManager::operator/=( const double & arg )
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::operator/= only works for previously finalized instances";
    throw eRatesManager_error( errMsg );
  }

  // four independent loops would probably be faster (cache lines!)
  for ( int j = 0; j < 7; j++ )
  {
    ratesFlav22[j] /= arg;
    ratesFlav23[j] /= arg;
    ratesFlav32[j] /= arg;
    ratesFlavTotal[j] /= arg;
  }

  for ( int j = 0; j < rates22.size(); j++ )
  {
    rates22[j] /= arg;
  }

  for ( int j = 0; j < rates23.size(); j++ )
  {
    rates23[j] /= arg;
  }

  for ( int j = 0; j < rates32.size(); j++ )
  {
    rates32[j] /= arg;
  }
  
  return ( *this );
}



ratesManager& ratesManager::operator*=( const double & arg )
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::operator*= only works for previously finalized instances";
    throw eRatesManager_error( errMsg );
  }

  // four independent loops would probably be faster (cache lines!)
  for ( int j = 0; j < 7; j++ )
  {
    ratesFlav22[j] *= arg;
    ratesFlav23[j] *= arg;
    ratesFlav32[j] *= arg;
    ratesFlavTotal[j] *= arg;
  }

  for ( int j = 0; j < rates22.size(); j++ )
  {
    rates22[j] *= arg;
  }

  for ( int j = 0; j < rates23.size(); j++ )
  {
    rates23[j] *= arg;
  }

  for ( int j = 0; j < rates32.size(); j++ )
  {
    rates32[j] *= arg;
  }
  
  return ( *this );
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
