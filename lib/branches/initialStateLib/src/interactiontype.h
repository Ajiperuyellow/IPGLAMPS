//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interactiontype.h $
//$LastChangedDate: 2014-11-18 16:05:52 +0100 (Di, 18. Nov 2014) $
//$LastChangedRevision: 1947 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



#ifndef INTERACTIONTYPE_H
#define INTERACTIONTYPE_H

#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include "particleprototype.h"


/**
 * @brief enum for generic 2->2, 2->3 and 3->2 interactions
 */
enum GENERIC_COLL_TYPE { 
  c22, ///< 0 = 2 -> 2
  c23, ///< 1 = 2 -> 3
  c32  ///< 2 = 3 -> 2
};

/**
 * @brief Provide classification of scattering processes
 *
 * @note
 * Many features and routines aren't used (yet). But may be helpful in
 * the future.  
 *
 *
 * Interactions types are mapped according to the following scheme:
 *
 * @par 2 -> 2
 * light partons:
 * @arg @c 9221 : g+g -> X
 * @arg @c 9222 : g+q -> X
 * @arg @c 9223 : g+qbar -> X
 * @arg @c 9224 : q+qbar -> X
 * @arg @c 9225 : q+q -> X
 * @arg @c 9226 : qbar+qbar -> X
 * @arg @c 9227 : q+q' -> X
 * @arg @c 9228 : q+qbar' -> X
 * @arg @c 9229 : qbar+qbar' -> X
 *
 * heavy partons:
 * @arg @c 92210 : g+c -> X
 * @arg @c 92211 : g+cbar -> X
 * @arg @c 92212 : q+c -> X
 * @arg @c 92213 : q+cbar -> X
 * @arg @c 92214 : qbar+c -> X
 * @arg @c 92215 : qbar+cbar -> X
 * @arg @c 92216 : c+cbar -> X
 * @arg @c 92220 : g+b -> X
 * @arg @c 92221 : g+bbar -> X
 * @arg @c 92222 : q+b -> X
 * @arg @c 92223 : q+bbar -> X
 * @arg @c 92224 : qbar+b -> X
 * @arg @c 92225 : qbar+bbar -> X
 * @arg @c 92226 : b+bbar -> X
 *
 * @arg @c 92240 : g+Jpsi -> X
 *
 *
 * light partons:
 * @arg @c 221 : g+g -> g+g
 * @arg @c 222 : g+g -> q+qbar
 * @arg @c 223 : g+q -> g+q                
 * @arg @c 2230 : g+qbar -> g+qbar
 * @arg @c 224 : q+qbar -> q+qbar
 * @arg @c 225 : q+qbar -> q'+qbar'
 * @arg @c 226 : q+qbar -> g+g
 * @arg @c 227 : q+q -> q+q                
 * @arg @c 2270 : qbar+qbar -> qbar+qbar
 * @arg @c 228 : q+q'-> q+q'              
 * @arg @c 2280 : q+qbar'-> q+qbar'           
 * @arg @c 22800 : qbar+qbar' -> qbar+qbar'
 *
 * heavy partons:
 * @arg @c 2210 : g+g -> c+cbar
 * @arg @c 2211 : q+qbar -> c+cbar
 * @arg @c 2212 : c+cbar -> g+g
 * @arg @c 2213 : c+cbar -> q+qbar
 * @arg @c 2214 : g+c -> g+c               
 * @arg @c 22140 : g+cbar -> g+cbar
 * @arg @c 2215 : q+c-> q+c                
 * @arg @c 22150 : q+cbar-> q+cbar            
 * @arg @c 221500 : qbar+c -> qbar+c            
 * @arg @c 2215000 : qbar+cbar -> qbar+cbar
 * @arg @c 2220 : g+g -> b+bbar
 * @arg @c 2221 : q+qbar -> b+bbar
 * @arg @c 2222 : b+bbar -> g+g
 * @arg @c 2223 : b+bbar -> q+qbar
 * @arg @c 2224 : g+b -> g+b               
 * @arg @c 22240 : g+bbar -> g+bbar
 * @arg @c 2225 : q+b-> q+b                
 * @arg @c 22250 : q+bbar-> q+bbar            
 * @arg @c 222500 : qbar+b -> qbar+b            
 * @arg @c 2225000 : qbar+bbar -> qbar+bbar
 *
 * @arg @c 2240 : g+Jpsi -> c+cbar
 * @arg @c 2241 : c+cbar -> g+Jpsi
 *
 *
 *
 * @par 2 -> 3
 * light partons:
 * @arg @c 9231 : g+g -> X+g
 * @arg @c 9232 : g+q -> X+g
 * @arg @c 9233 : g+qbar -> X+g
 * @arg @c 9234 : q+qbar -> X+g
 * @arg @c 9235 : q+q -> X+g
 * @arg @c 9236 : qbar+qbar -> X+g
 * @arg @c 9237 : q+q' -> X+g
 * @arg @c 9238 : q+qbar' -> X+g
 * @arg @c 9239 : qbar+qbar' -> X+g
 *
 * heavy partons:
 * @arg @c 92310 : g+c -> X+g
 * @arg @c 92311 : g+cbar -> X+g
 * @arg @c 92312 : q+c -> X+g
 * @arg @c 92313 : q+cbar -> X+g
 * @arg @c 92314 : qbar+c -> X+g
 * @arg @c 92315 : qbar+cbar -> X+g
 * @arg @c 92316 : c+cbar -> X+g
 * @arg @c 92320 : g+b -> X+g
 * @arg @c 92321 : g+bbar -> X+g
 * @arg @c 92322 : q+b -> X+g
 * @arg @c 92323 : q+bbar -> X+g
 * @arg @c 92324 : qbar+b -> X+g
 * @arg @c 92325 : qbar+bbar -> X+g
 * @arg @c 92326 : b+bbar -> X+g
 *
 * light partons:
 * @arg @c 231 : g+g -> g+g+g
 * @arg @c 232 : g+g -> q+qbar+g
 * @arg @c 233 : g+q -> g+q+g              
 * @arg @c 2330 : g+qbar -> g+qbar+g
 * @arg @c 234 : q+qbar -> q+qbar+g
 * @arg @c 235 : q+qbar -> q'+qbar'+g
 * @arg @c 236 : q+qbar -> g+g+g
 * @arg @c 237 : q+q -> q+q+g              
 * @arg @c 2370 : qbar+qbar -> qbar+qbar+g
 * @arg @c 238 : q+q'-> q+q'+g             
 * @arg @c 2380 : q+qbar'-> q+qbar'+g           
 * @arg @c 23800 : qbar+qbar' -> qbar+qbar'+g
 *
 * heavy partons:
 * @arg @c 2310 : g+g -> c+cbar+g
 * @arg @c 2311 : q+qbar -> c+cbar+g
 * @arg @c 2312 : c+cbar -> g+g+g
 * @arg @c 2313 : c+cbar -> q+qbar+g
 * @arg @c 2314 : g+c -> g+c+g             
 * @arg @c 23140 : g+cbar -> g+cbar+g
 * @arg @c 2315 : q+c-> q+c+g              
 * @arg @c 23150 : q+cbar-> q+cbar+g          
 * @arg @c 231500 : qbar+c -> qbar+c+g          
 * @arg @c 2315000 : qbar+cbar -> qbar+cbar+g
 * @arg @c 2320 : g+g -> b+bbar+g
 * @arg @c 2321 : q+qbar -> b+bbar+g
 * @arg @c 2322 : b+bbar -> g+g+g
 * @arg @c 2323 : b+bbar -> q+qbar+g
 * @arg @c 2324 : g+b -> g+b+g             
 * @arg @c 23240 : g+bbar -> g+bbar+g
 * @arg @c 2325 : q+b-> q+b+g              
 * @arg @c 23250 : q+bbar-> q+bbar+g          
 * @arg @c 232500 : qbar+b -> qbar+b+g          
 * @arg @c 2325000 : qbar+bbar -> qbar+bbar+g
 *
 *
 * @par 3 -> 2
 * light partons:
 * @arg @c 9321 : g+g+g -> X
 * @arg @c 9322 : g+q+g -> X
 * @arg @c 9323 : g+qbar+g -> X
 * @arg @c 9324 : q+qbar+g -> X
 * @arg @c 9325 : q+q+g -> X
 * @arg @c 9326 : qbar+qbar+g -> X
 * @arg @c 9327 : q+q'+g -> X
 * @arg @c 9328 : q+qbar'+g -> X
 * @arg @c 9329 : qbar+qbar'+g -> X
 *
 * heavy partons:
 * @arg @c 93210 : g+c+g -> X
 * @arg @c 93211 : g+cbar+g -> X
 * @arg @c 93212 : q+c+g -> X
 * @arg @c 93213 : q+cbar+g -> X
 * @arg @c 93214 : qbar+c+g -> X
 * @arg @c 93215 : qbar+cbar+g -> X
 * @arg @c 93216 : c+cbar+g -> X
 * @arg @c 93220 : g+b+g -> X
 * @arg @c 93221 : g+bbar+g -> X
 * @arg @c 93222 : q+b+g -> X
 * @arg @c 93223 : q+bbar+g -> X
 * @arg @c 93224 : qbar+b+g -> X
 * @arg @c 93225 : qbar+bbar+g -> X
 * @arg @c 93226 : b+bbar+g -> X
 *
 * light partons:
 * @arg @c 321 : g+g+g -> g+g
 * @arg @c 322 : g+g+g -> q+qbar
 * @arg @c 323 : g+q+g -> g+q              
 * @arg @c 3230 : g+qbar+g -> g+qbar
 * @arg @c 324 : q+qbar+g -> q+qbar
 * @arg @c 325 : q+qbar+g -> q'+qbar'
 * @arg @c 326 : q+qbar+g -> g+g
 * @arg @c 327 : q+q+g -> q+q              
 * @arg @c 3270 : qbar+qbar+g -> qbar+qbar
 * @arg @c 328 : q+q'+g-> q+q'             
 * @arg @c 3280 : q+qbar'+g-> q+qbar'         
 * @arg @c 32800 : qbar+qbar'+g -> qbar+qbar'
 *
 * heavy partons:
 * @arg @c 3210 : g+g+g -> c+cbar
 * @arg @c 3211 : q+qbar+g -> c+cbar
 * @arg @c 3212 : c+cbar+g -> g+g
 * @arg @c 3213 : c+cbar+g -> q+qbar
 * @arg @c 3214 : g+c+g -> g+c             
 * @arg @c 32140 : g+cbar+g -> g+cbar
 * @arg @c 3215 : q+c+g -> q+c             
 * @arg @c 32150 : q+cbar+g -> q+cbar         
 * @arg @c 321500 : qbar+c+g -> qbar+c          
 * @arg @c 3215000 : qbar+cbar+g -> qbar+cbar
 * @arg @c 3220 : g+g+g -> b+bbar
 * @arg @c 3221 : q+qbar+g -> b+bbar
 * @arg @c 3222 : b+bbar+g -> g+g
 * @arg @c 3223 : b+bbar+g -> q+qbar
 * @arg @c 3224 : g+b+g -> g+b             
 * @arg @c 32240 : g+bbar+g -> g+bbar
 * @arg @c 3225 : q+b+g -> q+b             
 * @arg @c 32250 : q+bbar+g -> q+bbar         
 * @arg @c 322500 : qbar+b+g -> qbar+b          
 * @arg @c 3225000 : qbar+bbar+g -> qbar+bbar
 *
 */
class interactionType
{
public:
  /** @brief Standard constructor */
  interactionType() : id( 0 ) {};
  /** @brief Standard destructor */
  ~interactionType() {};

  /** @brief id of the interaction type */
  int id;

  /** @brief Returns the number of particles of type _F involved in the initial state of interaction type id */
  double getInvolvedInitial( const FLAVOR_TYPE _F ) const
  {
    return getInvolvedInitial( _F, id );
  }
  /** @brief (static version) Returns the number of particles of type _F involved in the initial state of interaction type id */
  static double getInvolvedInitial( const FLAVOR_TYPE _F, const int _id );

  
  /** @brief Returns the number of particles of type _F involved in the initial state of inclusive interaction type id */
  double getInvolvedInitialInclusive( const FLAVOR_TYPE _F ) const
  {
    return getInvolvedInitialInclusive( _F, id );
  }
  /** @brief (static version) Returns the number of particles of type _F involved in the initial state of inclusive interaction type id */
  static double getInvolvedInitialInclusive( const FLAVOR_TYPE _F, const int _id );
  
  /** @brief Get scaling factor for identical particles in the initial state for given inclusive process of type _id*/
  static int getScalingFactorForIdenticalParticles( const int _id );
    
  /** @brief Returns a vector of interaction type ids that involve particles with flavor _F */
  static std::vector<int> getInteractionTypes( const FLAVOR_TYPE _F );

  /** @brief Generic collision type (2->2 etc) from specific interaction type */
  static GENERIC_COLL_TYPE getCollType( const int _id );
  
  /** @brief Get the index as denoted in processIndices22 (23, 32) for the given process */
  int getIndexFromProcessType( const GENERIC_COLL_TYPE _genType ) const { return getIndexFromProcessType(id, _genType) ;}
  /** @brief (static version) Get the index as denoted in processIndices22 (23, 32) for the given process */
  static int getIndexFromProcessType( const int _id, const GENERIC_COLL_TYPE _genType );
  
  /** @brief (static version) Get the process id at a certain index for a given generic type (22, 23, 32) */
  static int getProcessTypeFromIndex( const unsigned int _index, const GENERIC_COLL_TYPE _genType );
  /** @brief (static version) Get the inclusive process id at a certain index for a given generic type (22, 23, 32) */
  static int getInclusiveProcessTypeFromIndex( const unsigned int _index, const GENERIC_COLL_TYPE _genType );
  
  /** @brief (static version) Get the inclusive process id from an exclusive process id for a given generic type (22, 23, 32) */
  static int getInclusiveProcessTypeFromExclusiveProcessType ( const unsigned int _exclProcess, const GENERIC_COLL_TYPE _genType );
  
  
  /** @brief Get the id for an inclusive process given the initial flavor (2->X) */
  static int getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const GENERIC_COLL_TYPE genericType);  
  /** @brief Get the id for an inclusive process given the initial flavor (3->2) */
  static int getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, const GENERIC_COLL_TYPE genericType );
  
  
  /** @brief Associates process id numbers for 2->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndices22;
  /** @brief Associates process id numbers for 2->3 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndices23;
 /** @brief Associates process id numbers for 3->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndices32;
  
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for 2->2 processes  */
  static std::vector<int> indexProcesses22;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for 2->3 processes  */
  static std::vector<int> indexProcesses23;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for 3->2 processes  */
  static std::vector<int> indexProcesses32;
  
   /** @brief Associates process id numbers for inclusive (specified initial state, arbitrary final state) 2->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndicesInclusive22;
  /** @brief Associates process id numbers for inclusive (specified initial state, arbitrary final state) 2->3 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndicesInclusive23;
  /** @brief Associates process id numbers for inclusive (specified initial state, arbitrary final state) 3->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndicesInclusive32;
  
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for inclusive 2->2 processes  */
  static std::vector<int> indexProcessesInclusive22;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for inclusive 2->3 processes  */
  static std::vector<int> indexProcessesInclusive23;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for inclusive 3->2 processes  */
  static std::vector<int> indexProcessesInclusive32;
  
  /** @brief Associate inclusive process id numbers with every exclusive process number */
  static std::map<int, int> inclusiveProcessFromExclusiveProcess22;
  
  /** @brief Associates process id numbers for inclusive processes with the number of gluons in the initial state */
  static std::map<int, int> involvedInitialGluons;
  /** @brief Associates process id numbers for inclusive processes with the number of light quarks in the initial state */
  static std::map<int, int> involvedInitialLightQuarks;
  /** @brief Associates process id numbers for inclusive processes with the number of light anti-quarks in the initial state */
  static std::map<int, int> involvedInitialAntiLightQuarks;
  /** @brief Associates process id numbers for inclusive processes with the number of charm quarks in the initial state */
  static std::map<int, int> involvedInitialCharmQuarks;
  /** @brief Associates process id numbers for inclusive processes with the number of charm anti-quarks in the initial state */
  static std::map<int, int> involvedInitialAntiCharmQuarks;
  /** @brief Associates process id numbers for inclusive processes with the number of bottom quarks in the initial state */
  static std::map<int, int> involvedInitialBottomQuarks;
  /** @brief Associates process id numbers for inclusive processes with the number of bottom anti-quarks in the initial state */
  static std::map<int, int> involvedInitialAntiBottomQuarks;

  /** @brief Associates process id numbers for inclusive processes with the scaling factor for identical particles in the initial state */
  static std::map<int, int> scalingFactorForIdenticalParticles;
  
private:
  
};


/** @brief exception class for handling unexpected behaviour when dealing with interaction types */
class eInteraction_type_error : public std::runtime_error
{
  public:
    explicit eInteraction_type_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eInteraction_type_error() throw() {};
};


#endif // INTERACTIONTYPE_H
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
