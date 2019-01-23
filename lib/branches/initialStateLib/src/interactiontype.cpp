//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interactiontype.cpp $
//$LastChangedDate: 2015-03-12 15:23:55 +0100 (Do, 12. Mär 2015) $
//$LastChangedRevision: 2122 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "interactiontype.h"

#include "tools.h"

/**
 * Initialize map containers that associate interaction type id
 * numbers with index numbers for use in vector access etc. 
 * When modifications are necessary, also change indexProcessesXX (see below)!
 */
std::map<int, int> interactionType::processIndices22 = create_map<int, int >
( 221, 0 )( 222, 1 )( 223, 2 )( 2230, 3 )( 224, 4 )( 225, 5 )( 226, 6 )( 227, 7 )( 2270, 8 )( 228, 9 )( 2280, 10 )( 22800, 11 )
( 2210, 12 )( 2211, 13 )( 2212, 14 )( 2213, 15 )( 2214, 16 )( 22140, 17 )( 2215, 18 )( 22150, 19 )( 221500, 20 )( 2215000, 21 )
( 2220, 22 )( 2221, 23 )( 2222, 24 )( 2223, 25 )( 2224, 26 )( 22240, 27 )( 2225, 28 )( 22250, 29 )( 222500, 30 )( 2225000, 31 )
( 2240, 32 )( 2241, 33 );
std::map<int, int> interactionType::processIndices23 = create_map<int, int >
( 231, 0 )( 232, 1 )( 233, 2 )( 2330, 3 )( 234, 4 )( 235, 5 )( 236, 6 )( 237, 7 )( 2370, 8 )( 238, 9 )( 2380, 10 )( 23800, 11 )
( 2310, 12 )( 2311, 13 )( 2312, 14 )( 2313, 15 )( 2314, 16 )( 23140, 17 )( 2315, 18 )( 23150, 19 )( 231500, 20 )( 2315000, 21 )
( 2320, 22 )( 2321, 23 )( 2322, 24 )( 2323, 25 )( 2324, 26 )( 23240, 27 )( 2325, 28 )( 23250, 29 )( 232500, 30 )( 2325000, 31 );
std::map<int, int> interactionType::processIndices32 = create_map<int, int >
( 321, 0 )( 322, 1 )( 323, 2 )( 3230, 3 )( 324, 4 )( 325, 5 )( 326, 6 )( 327, 7 )( 3270, 8 )( 328, 9 )( 3280, 10 )( 32800, 11 )
( 3210, 12 )( 3211, 13 )( 3212, 14 )( 3213, 15 )( 3214, 16 )( 32140, 17 )( 3215, 18 )( 32150, 19 )( 321500, 20 )( 3215000, 21 )
( 3220, 22 )( 3221, 23 )( 3222, 24 )( 3223, 25 )( 3224, 26 )( 32240, 27 )( 3225, 28 )( 32250, 29 )( 322500, 30 )( 3225000, 31 );


/**
 * Initialize vector containers that associate index numbers (for use
 * in vector access etc.) with interaction type id numbers  
 * When modifications are necessary, also change processIndicesXX (see
 * above)! 
 */
std::vector<int> interactionType::indexProcesses22 = create_vector<int>
( 221 )( 222 )( 223 )( 2230 )( 224 )( 225 )( 226 )( 227 )( 2270 )( 228 )( 2280 )( 22800 )
( 2210 )( 2211 )( 2212 )( 2213 )( 2214 )( 22140 )( 2215 )( 22150 )( 221500 )( 2215000 )
( 2220 )( 2221 )( 2222 )( 2223 )( 2224 )( 22240 )( 2225 )( 22250 )( 222500 )( 2225000 )
( 2240 )( 2241 );
std::vector<int> interactionType::indexProcesses23 = create_vector<int>
( 231 )( 232 )( 233 )( 2330 )( 234 )( 235 )( 236 )( 237 )( 2370 )( 238 )( 2380 )( 23800 )
( 2310 )( 2311 )( 2312 )( 2313 )( 2314 )( 23140 )( 2315 )( 23150 )( 231500 )( 2315000 )
( 2320 )( 2321 )( 2322 )( 2323 )( 2324 )( 23240 )( 2325 )( 23250 )( 232500 )( 2325000 );
std::vector<int> interactionType::indexProcesses32 = create_vector<int>
( 321 )( 322 )( 323 )( 3230 )( 324 )( 325 )( 326 )( 327 )( 3270 )( 328 )( 3280 )( 32800 )
( 3210 )( 3211 )( 3212 )( 3213 )( 3214 )( 32140 )( 3215 )( 32150 )( 321500 )( 3215000 )
( 3220 )( 3221 )( 3222 )( 3223 )( 3224 )( 32240 )( 3225 )( 32250 )( 322500 )( 3225000 );


/**
 * Initialize map containers that associate interaction type id
 * numbers for inclusive processes (specified initial state, arbitrary
 * final state, IDs starting with 9) with index numbers for use in
 * vector access etc. 
 * When modifications are necessary, also change
 * indexProcessesInclusiveXX (see below)! 
 */
std::map<int, int> interactionType::processIndicesInclusive22 = create_map<int, int >
( 9221, 0 )( 9222, 1 )( 9223, 2 )( 9224, 3 )( 9225, 4 )( 9226, 5 )( 9227, 6 )( 9228, 7 )( 9229, 8 )
( 92210,  9 )( 92211, 10 )( 92212, 11 )( 92213, 12 )( 92214, 13 )( 92215, 14 )( 92216, 15 )
( 92220, 16 )( 92221, 17 )( 92222, 18 )( 92223, 19 )( 92224, 20 )( 92225, 21 )( 92226, 22 )
( 92240, 23 );
std::map<int, int> interactionType::processIndicesInclusive23 = create_map<int, int >
( 9231, 0 )( 9232, 1 )( 9233, 2 )( 9234, 3 )( 9235, 4 )( 9236, 5 )( 9237, 6 )( 9238, 7 )( 9239, 8 )
( 92310,  9 )( 92311, 10 )( 92312, 11 )( 92313, 12 )( 92314, 13 )( 92315, 14 )( 92316, 15 )
( 92320, 16 )( 92321, 17 )( 92322, 18 )( 92323, 19 )( 92324, 20 )( 92325, 21 )( 92326, 22 );
std::map<int, int> interactionType::processIndicesInclusive32 = create_map<int, int >
( 9321, 0 )( 9322, 1 )( 9323, 2 )( 9324, 3 )( 9325, 4 )( 9326, 5 )( 9327, 6 )( 9328, 7 )( 9329, 8 )
( 93210,  9 )( 93211, 10 )( 93212, 11 )( 93213, 12 )( 93214, 13 )( 93215, 14 )( 93216, 15 )
( 93220, 16 )( 93221, 17 )( 93222, 18 )( 93223, 19 )( 93224, 20 )( 93225, 21 )( 93226, 22 );


/**
 * Initialize vector containers that associate index numbers (for use
 * in vector access etc.) with inclusive interaction type id numbers
 * When modifications are necessary, also change
 * processIndicesInclusiveXX (see above)! 
 */
std::vector<int> interactionType::indexProcessesInclusive22 = create_vector<int>
( 9221 )( 9222 )( 9223 )( 9224 )( 9225 )( 9226 )( 9227 )( 9228 )( 9229 )
( 92210 )( 92211 )( 92212 )( 92213 )( 92214 )( 92215 )( 92216 )
( 92220 )( 92221 )( 92222 )( 92223 )( 92224 )( 92225 )( 92226 )
( 92240 );
std::vector<int> interactionType::indexProcessesInclusive23 = create_vector<int>
( 9231 )( 9232 )( 9233 )( 9234 )( 9235 )( 9236 )( 9237 )( 9238 )( 9239 )
( 92310 )( 92311 )( 92312 )( 92313 )( 92314 )( 92315 )( 92316 )
( 92320 )( 92321 )( 92322 )( 92323 )( 92324 )( 92325 )( 92326 );
std::vector<int> interactionType::indexProcessesInclusive32 = create_vector<int>
( 9321 )( 9322 )( 9323 )( 9324 )( 9325 )( 9326 )( 9327 )( 9328 )( 9329 )
( 93210 )( 93211 )( 93212 )( 93213 )( 93214 )( 93215 )( 93216 )
( 93220 )( 93221 )( 93222 )( 93223 )( 93224 )( 93225 )( 93226 );


/**
 * Initialize map containers that associate an interaction type id
 * numbers for inclusive processes with every exclusive process type
 * id 
 */
std::map< int, int > interactionType::inclusiveProcessFromExclusiveProcess22 = create_map< int, int >
( 221, 9221 )( 222, 9221 )( 223, 9222 )( 2230, 9223 )( 224, 9224 )( 225, 9224 )( 226, 9224 )( 227, 9225 )( 2270, 9226 )( 228, 9227 )( 2280, 9228 )( 22800, 9229 )
( 2210, 9221 )( 2211, 9224 )( 2212, 92216 )( 2213, 92216 )( 2214, 92210 )( 22140, 92211 )( 2215, 92212 )( 22150, 92213 )( 221500, 92214 )( 2215000, 92215 )
( 2220, 9221 )( 2221, 9224 )( 2222, 92226 )( 2223, 92226 )( 2224, 92220 )( 22240, 92221 )( 2225, 92222 )( 22250, 92223 )( 222500, 92224 )( 2225000, 92225 )
( 2240, 92240 )( 2241, 92216 );


/**
* Initialize map containers that associate interaction type id numbers
* with the number of involved particles in the initial state. 
* So far only IDs starting with 9, i.e. inclusive processes
*/
std::map<int, int> interactionType::involvedInitialGluons = create_map<int, int >
( 9221, 2 )( 9222, 1 )( 9223, 1 )( 9224, 0 )( 9225, 0 )( 9226, 0 )( 9227, 0 )( 9228, 0 )( 9229, 0 )
( 92210, 1 )( 92211, 1 )( 92212, 0 )( 92213, 0 )( 92214, 0 )( 92215, 0 )( 92216, 0 )
( 92220, 1 )( 92221, 1 )( 92222, 0 )( 92223, 0 )( 92224, 0 )( 92225, 0 )( 92226, 0 )
( 92240, 1 )
( 9231, 2 )( 9232, 1 )( 9233, 1 )( 9234, 0 )( 9235, 0 )( 9236, 0 )( 9237, 0 )( 9238, 0 )( 9239, 0 )
( 92310, 1 )( 92311, 1 )( 92312, 0 )( 92313, 0 )( 92314, 0 )( 92315, 0 )( 92316, 0 )
( 92320, 1 )( 92321, 1 )( 92322, 0 )( 92323, 0 )( 92324, 0 )( 92325, 0 )( 92326, 0 )
( 9321, 3 )( 9322, 2 )( 9323, 2 )( 9324, 1 )( 9325, 1 )( 9326, 1 )( 9327, 1 )( 9328, 1 )( 9329, 1 )
( 93210, 2 )( 93211, 2 )( 93212, 1 )( 93213, 1 )( 93214, 1 )( 93215, 1 )( 93216, 1 )
( 93220, 2 )( 93221, 2 )( 93222, 1 )( 93223, 1 )( 93224, 1 )( 93225, 1 )( 93226, 1 );
std::map<int, int> interactionType::involvedInitialLightQuarks = create_map<int, int >
( 9221, 0 )( 9222, 1 )( 9223, 0 )( 9224, 1 )( 9225, 2 )( 9226, 0 )( 9227, 2 )( 9228, 1 )( 9229, 0 )
( 92210, 0 )( 92211, 0 )( 92212, 1 )( 92213, 1 )( 92214, 0 )( 92215, 0 )( 92216, 0 )
( 92220, 0 )( 92221, 0 )( 92222, 1 )( 92223, 1 )( 92224, 0 )( 92225, 0 )( 92226, 0 )
( 92240, 0 )
( 9231, 0 )( 9232, 1 )( 9233, 0 )( 9234, 1 )( 9235, 2 )( 9236, 0 )( 9237, 2 )( 9238, 1 )( 9239, 0 )
( 92310, 0 )( 92311, 0 )( 92312, 1 )( 92313, 1 )( 92314, 0 )( 92315, 0 )( 92316, 0 )
( 92320, 0 )( 92321, 0 )( 92322, 1 )( 92323, 1 )( 92324, 0 )( 92325, 0 )( 92326, 0 )
( 9321, 0 )( 9322, 1 )( 9323, 0 )( 9324, 1 )( 9325, 2 )( 9326, 0 )( 9327, 2 )( 9328, 1 )( 9329, 0 )
( 93210, 0 )( 93211, 0 )( 93212, 1 )( 93213, 1 )( 93214, 0 )( 93215, 0 )( 93216, 0 )
( 93220, 0 )( 93221, 0 )( 93222, 1 )( 93223, 1 )( 93224, 0 )( 93225, 0 )( 93226, 0 );
std::map<int, int> interactionType::involvedInitialAntiLightQuarks = create_map<int, int >
( 9221, 0 )( 9222, 0 )( 9223, 1 )( 9224, 1 )( 9225, 0 )( 9226, 2 )( 9227, 0 )( 9228, 1 )( 9229, 2 )
( 92210, 0 )( 92211, 0 )( 92212, 0 )( 92213, 0 )( 92214, 1 )( 92215, 1 )( 92216, 0 )
( 92220, 0 )( 92221, 0 )( 92222, 0 )( 92223, 0 )( 92224, 1 )( 92225, 1 )( 92226, 0 )
( 92240, 0 )
( 9231, 0 )( 9232, 0 )( 9233, 1 )( 9234, 1 )( 9235, 0 )( 9236, 2 )( 9237, 0 )( 9238, 1 )( 9239, 2 )
( 92310, 0 )( 92311, 0 )( 92312, 0 )( 92313, 0 )( 92314, 1 )( 92315, 1 )( 92316, 0 )
( 92320, 0 )( 92321, 0 )( 92322, 0 )( 92323, 0 )( 92324, 1 )( 92325, 1 )( 92326, 0 )
( 9321, 0 )( 9322, 0 )( 9323, 1 )( 9324, 1 )( 9325, 0 )( 9326, 2 )( 9327, 0 )( 9328, 1 )( 9329, 2 )
( 93210, 0 )( 93211, 0 )( 93212, 0 )( 93213, 0 )( 93214, 1 )( 93215, 1 )( 93216, 0 )
( 93220, 0 )( 93221, 0 )( 93222, 0 )( 93223, 0 )( 93224, 1 )( 93225, 1 )( 93226, 0 );
std::map<int, int> interactionType::involvedInitialCharmQuarks = create_map<int, int >
( 9221, 0 )( 9222, 0 )( 9223, 0 )( 9224, 0 )( 9225, 0 )( 9226, 0 )( 9227, 0 )( 9228, 0 )( 9229, 0 )
( 92210, 1 )( 92211, 0 )( 92212, 1 )( 92213, 0 )( 92214, 1 )( 92215, 0 )( 92216, 1 )
( 92220, 0 )( 92221, 0 )( 92222, 0 )( 92223, 0 )( 92224, 0 )( 92225, 0 )( 92226, 0 )
( 92240, 0 )
( 9231, 0 )( 9232, 0 )( 9233, 0 )( 9234, 0 )( 9235, 0 )( 9236, 0 )( 9237, 0 )( 9238, 0 )( 9239, 0 )
( 92310, 1 )( 92311, 0 )( 92312, 1 )( 92313, 0 )( 92314, 1 )( 92315, 0 )( 92316, 1 )
( 92320, 0 )( 92321, 0 )( 92322, 0 )( 92323, 0 )( 92324, 0 )( 92325, 0 )( 92326, 0 )
( 9321, 0 )( 9322, 0 )( 9323, 0 )( 9324, 0 )( 9325, 0 )( 9326, 0 )( 9327, 0 )( 9328, 0 )( 9329, 0 )
( 93210, 1 )( 93211, 0 )( 93212, 1 )( 93213, 0 )( 93214, 1 )( 93215, 0 )( 93216, 1 )
( 93220, 0 )( 93221, 0 )( 93222, 0 )( 93223, 0 )( 93224, 0 )( 93225, 0 )( 93226, 0 );
std::map<int, int> interactionType::involvedInitialAntiCharmQuarks = create_map<int, int >
( 9221, 0 )( 9222, 0 )( 9223, 0 )( 9224, 0 )( 9225, 0 )( 9226, 0 )( 9227, 0 )( 9228, 0 )( 9229, 0 )
( 92210, 0 )( 92211, 1 )( 92212, 0 )( 92213, 1 )( 92214, 0 )( 92215, 1 )( 92216, 1 )
( 92220, 0 )( 92221, 0 )( 92222, 0 )( 92223, 0 )( 92224, 0 )( 92225, 0 )( 92226, 0 )
( 92240, 0 )
( 9231, 0 )( 9232, 0 )( 9233, 0 )( 9234, 0 )( 9235, 0 )( 9236, 0 )( 9237, 0 )( 9238, 0 )( 9239, 0 )
( 92310, 0 )( 92311, 1 )( 92312, 0 )( 92313, 1 )( 92314, 0 )( 92315, 1 )( 92316, 1 )
( 92320, 0 )( 92321, 0 )( 92322, 0 )( 92323, 0 )( 92324, 0 )( 92325, 0 )( 92326, 0 )
( 9321, 0 )( 9322, 0 )( 9323, 0 )( 9324, 0 )( 9325, 0 )( 9326, 0 )( 9327, 0 )( 9328, 0 )( 9329, 0 )
( 93210, 0 )( 93211, 1 )( 93212, 0 )( 93213, 1 )( 93214, 0 )( 93215, 1 )( 93216, 1 )
( 93220, 0 )( 93221, 0 )( 93222, 0 )( 93223, 0 )( 93224, 0 )( 93225, 0 )( 93226, 0 );
std::map<int, int> interactionType::involvedInitialBottomQuarks = create_map<int, int >
( 9221, 0 )( 9222, 0 )( 9223, 0 )( 9224, 0 )( 9225, 0 )( 9226, 0 )( 9227, 0 )( 9228, 0 )( 9229, 0 )
( 92210, 0 )( 92211, 0 )( 92212, 0 )( 92213, 0 )( 92214, 0 )( 92215, 0 )( 92216, 0 )
( 92220, 1 )( 92221, 0 )( 92222, 1 )( 92223, 0 )( 92224, 1 )( 92225, 0 )( 92226, 1 )
( 92240, 0 )
( 9231, 0 )( 9232, 0 )( 9233, 0 )( 9234, 0 )( 9235, 0 )( 9236, 0 )( 9237, 0 )( 9238, 0 )( 9239, 0 )
( 92310, 0 )( 92311, 0 )( 92312, 0 )( 92313, 0 )( 92314, 0 )( 92315, 0 )( 92316, 0 )
( 92320, 1 )( 92321, 0 )( 92322, 1 )( 92323, 0 )( 92324, 1 )( 92325, 0 )( 92326, 1 )
( 9321, 0 )( 9322, 0 )( 9323, 0 )( 9324, 0 )( 9325, 0 )( 9326, 0 )( 9327, 0 )( 9328, 0 )( 9329, 0 )
( 93210, 0 )( 93211, 0 )( 93212, 0 )( 93213, 0 )( 93214, 0 )( 93215, 0 )( 93216, 0 )
( 93220, 1 )( 93221, 0 )( 93222, 1 )( 93223, 0 )( 93224, 1 )( 93225, 0 )( 93226, 1 );
std::map<int, int> interactionType::involvedInitialAntiBottomQuarks = create_map<int, int >
( 9221, 0 )( 9222, 0 )( 9223, 0 )( 9224, 0 )( 9225, 0 )( 9226, 0 )( 9227, 0 )( 9228, 0 )( 9229, 0 )
( 92210, 0 )( 92211, 0 )( 92212, 0 )( 92213, 0 )( 92214, 0 )( 92215, 0 )( 92216, 0 )
( 92220, 0 )( 92221, 1 )( 92222, 0 )( 92223, 1 )( 92224, 0 )( 92225, 1 )( 92226, 1 )
( 92240, 0 )
( 9231, 0 )( 9232, 0 )( 9233, 0 )( 9234, 0 )( 9235, 0 )( 9236, 0 )( 9237, 0 )( 9238, 0 )( 9239, 0 )
( 92310, 0 )( 92311, 0 )( 92312, 0 )( 92313, 0 )( 92314, 0 )( 92315, 0 )( 92316, 0 )
( 92320, 0 )( 92321, 1 )( 92322, 0 )( 92323, 1 )( 92324, 0 )( 92325, 1 )( 92326, 1 )
( 9321, 0 )( 9322, 0 )( 9323, 0 )( 9324, 0 )( 9325, 0 )( 9326, 0 )( 9327, 0 )( 9328, 0 )( 9329, 0 )
( 93210, 0 )( 93211, 0 )( 93212, 0 )( 93213, 0 )( 93214, 0 )( 93215, 0 )( 93216, 0 )
( 93220, 0 )( 93221, 1 )( 93222, 0 )( 93223, 1 )( 93224, 0 )( 93225, 1 )( 93226, 1 );

/**
* Initialize map container that associate interaction type id numbers
* with the scaling factor for indentical particles in the initial
* state. 
* For example in g+g -> X this would be 2! = 2
*/
std::map<int, int> interactionType::scalingFactorForIdenticalParticles = create_map<int, int >
( 9221, 2 )( 9222, 1 )( 9223, 1 )( 9224, 1 )( 9225, 2 )( 9226, 2 )( 9227,1  )( 9228, 1 )( 9229, 1 )
( 92210, 1 )( 92211, 1 )( 92212, 1 )( 92213, 1 )( 92214, 1 )( 92215, 1 )( 92216, 1 )
( 92220, 1 )( 92221, 1 )( 92222, 1 )( 92223, 1 )( 92224, 1 )( 92225, 1 )( 92226, 1 )
( 92240, 1 )
( 9231, 2 )( 9232, 1 )( 9233, 1 )( 9234, 1 )( 9235, 2 )( 9236, 2 )( 9237, 1 )( 9238, 1 )( 9239, 1 )
( 92310, 1 )( 92311, 1 )( 92312, 1 )( 92313, 1 )( 92314, 1 )( 92315, 1 )( 92316, 1 )
( 92320, 1 )( 92321, 1 )( 92322, 1 )( 92323, 1 )( 92324, 1 )( 92325, 1 )( 92326, 1 )
( 9321, 6 )( 9322, 2 )( 9323, 2 )( 9324, 1 )( 9325, 2 )( 9326, 2 )( 9327, 1 )( 9328, 1 )( 9329, 1 )
( 93210, 2 )( 93211, 2 )( 93212, 1 )( 93213, 1 )( 93214, 1 )( 93215, 1 )( 93216, 1 )
( 93220, 2 )( 93221, 2 )( 93222, 1 )( 93223, 1 )( 93224, 1 )( 93225, 1 )( 93226, 1 );




/**
 * Look up the process denoted by _id in the std::maps and return the
 * appropriate index number. 
 *
 * @param[in] _id process id to be looked up (231, 32800 etc.)
 * @param[in] _genType generic collision type the specific process belongs to
 * @return the index number as denoted in the std:maps
 */
int interactionType::getIndexFromProcessType( const int _id, const GENERIC_COLL_TYPE _genType )
{
  std::map<int, int>::const_iterator find_iter;

  switch ( _genType )
  {
  case c22:
    find_iter = processIndicesInclusive22.find( _id );
    if ( find_iter != processIndicesInclusive22.end() )
    {
      return find_iter->second;
    }
    else
    {
      find_iter = processIndices22.find( _id );
      if ( find_iter != processIndices22.end() )
      {
        return find_iter->second;
      }
      else
      {
        std::string errMsg = "Determination of process index for 2->2 process failed. Unrecoverable error.";
        throw eInteraction_type_error( errMsg );
      }
    }
    break;
  case c23:
    find_iter = processIndicesInclusive23.find( _id );
    if ( find_iter != processIndicesInclusive23.end() )
    {
      return find_iter->second;
    }
    else
    {
      find_iter = processIndices23.find( _id );
      if ( find_iter != processIndices23.end() )
      {
        return find_iter->second;
      }
      else
      {
        std::string errMsg = "Determination of process index for 2->3 process failed. Unrecoverable error.";
        throw eInteraction_type_error( errMsg );
      }
    }
    break;
  case c32:
    find_iter = processIndicesInclusive32.find( _id );
    if ( find_iter != processIndicesInclusive32.end() )
    {
      return find_iter->second;
    }
    else
    {
      find_iter = processIndices32.find( _id );
      if ( find_iter != processIndices32.end() )
      {
        return find_iter->second;
      }
      else
      {
        std::string errMsg = "Determination of process index for 3->2 process failed. Unrecoverable error.";
        throw eInteraction_type_error( errMsg );
      }
    }
    break;
  default:
    std::string errMsg = "Determination of process index failed. Bad generic type. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}



/**
 * Look up the process denoted by _id in the std::maps and return the appropriate index number.
 *
 * @param[in] _index index that needs to be looked up
 * @param[in] _genType generic collision type the specific process belongs to
 * @return process id to be looked up (231, 32800 etc.)
 */
int interactionType::getProcessTypeFromIndex( const unsigned int _index, const GENERIC_COLL_TYPE _genType )
{
  switch ( _genType )
  {
  case c22:
    if ( _index < indexProcesses22.size() )
    {
      return indexProcesses22[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c23:
    if ( _index < indexProcesses23.size() )
    {
      return indexProcesses23[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->3 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c32:
    if ( _index < indexProcesses32.size() )
    {
      return indexProcesses32[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 3->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  default:
    std::string errMsg = "Resolution of index failed. Bad generic type. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}



/**
* Look up the process index and return the inclusive process id
*
* @param[in] _index index that needs to be looked up
* @param[in] _genType generic collision type the specific process belongs to
* @return process id to be looked up (9231, 9223 etc.)
*/
int interactionType::getInclusiveProcessTypeFromIndex( const unsigned int _index, const GENERIC_COLL_TYPE _genType )
{
  switch ( _genType )
  {
  case c22:
    if ( _index < indexProcessesInclusive22.size() )
    {
      return indexProcessesInclusive22[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c23:
    if ( _index < indexProcessesInclusive23.size() )
    {
      return indexProcessesInclusive23[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->3 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c32:
    if ( _index < indexProcessesInclusive32.size() )
    {
      return indexProcessesInclusive32[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 3->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  default:
    std::string errMsg = "Resolution of index failed. Bad generic type. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}


/**
 * Look up the exclusive process denoted by _exclProcess in the std::maps and return the associated inclusive process id.
 *
 * @param[in] _exclProcess id of the exclusive process that is to be looked up
 * @param[in,out] _genType generic collision type the specific process belongs to
 * @return the inclusive process id associated with the exclusive process id
 */
int interactionType::getInclusiveProcessTypeFromExclusiveProcessType( const unsigned int _exclProcess, const GENERIC_COLL_TYPE _genType )
{
  if ( _genType == c22 )
  {
    std::map<int, int>::const_iterator find_iter;    
    
    find_iter = inclusiveProcessFromExclusiveProcess22.find( _exclProcess );
    if ( find_iter != inclusiveProcessFromExclusiveProcess22.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "interactionType::getInclusiveProcessTypeFromExclusiveProcessType: determination of process id for 2->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
  }
  else
  {
    std::string errMsg = "interactionType::getInclusiveProcessTypeFromExclusiveProcessType not yet written for process types other than 2->2";
    throw eInteraction_type_error( errMsg );
  }
}



/**
 * Look up the scaling factor for identical particles in the intial state for the inclusiv process as given by the process id.
 *
 * @param[in] _id id of the inclusive process that is to be looked up
 * @return Scaling factor for identical particles in the initial state (e.g. 2! = 2 for gg -> X)
 */
int interactionType::getScalingFactorForIdenticalParticles( const int _id )
{
  std::map<int, int>::const_iterator find_iter;    
  
  find_iter = scalingFactorForIdenticalParticles.find( _id );
  if ( find_iter != scalingFactorForIdenticalParticles.end() )
  {
    return find_iter->second;
  }
  else
  {
    std::string errMsg = "interactionType::getScalingFactorForIdenticalParticles: determination of process id failed. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
  }
}



#include <iostream>
/**
 * Given the flavors of the incoming particles, return the process type ID associated with this inclusive process (for 2->X)
 *
 * @param[in] _F1 Flavor of incoming particle 1
 * @param[in] _F2 Flavor of incoming particle 2
 * @param[in] genericType generic collision type (c22, c23, c32) the specific process belongs to
 * @return inclusive process type
 */
int interactionType::getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const GENERIC_COLL_TYPE genericType )
{
  if ( genericType == c32 )
  {
    std::string errMsg = "Generic interaction type is 3->2, but only two initial flavors are given.";
    throw eInteraction_type_error( errMsg );
  }
  
  int type = 0;
  bool error = false;

  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int F1 = std::min( static_cast<unsigned int>( _F1 ), static_cast<unsigned int>( _F2 ) );
  unsigned int F2 = std::max( static_cast<unsigned int>( _F1 ), static_cast<unsigned int>( _F2 ) );

  if ( F2 > 10 ) // One of the scattering partners is not a parton. These cross sections here are only for parton scatterings.
  {
    if( F1 == 0 && static_cast<FLAVOR_TYPE>( F2 ) == jpsi ) // Jpsi + g -> ccbar
    {
      type = 92240;  // g+Jpsi -> X
    }
    else
    {
      error = true;
      std::cout << "1" << std::endl;
    }
  }
  else if (( F1 + F2 ) == 0 ) // gg -> X
  {
    type = 9221;  // gg -> X
  }
  else if (( F1 * F2 ) == 0 && ( F2 % 2 ) == 0 ) // gqbar -> X
  {
    if( F2 <= 6 )
      type = 9223;  // gqbar -> X
    else if( F2 == 8 )
      type = 92211;  // gcbar -> X
    else if( F2 == 10 )
      type = 92221;  // gbbar -> X
    else
      error = true;
  }
  else if (( F1 * F2 ) == 0 ) // gq -> X
  {
    if( F2 <= 6 )
      type = 9222;  // gq -> X
    else if( F2 == 7 )
      type = 92210;  // gc -> X
    else if( F2 == 9 )
      type = 92220;  // gb -> X
    else
      error = true;
  }
  else if (( F1 % 2 ) == 1 && ( F2 % 2 ) == 1 )  // qq -> X, qq' -> X
  {
    if ( F1 == F2 )  // qq -> X
    {
      if ( F1 > 6 )
        error = true; // no elastic interactions among heavy quarks
      else
        type = 9225;  // qq -> X
    }
    else  // qq' -> X
    {
      if( F1 <= 6 )
      {
        if( F2 <= 6 )
          type = 9227;  // qq' -> X
        else if( F2 == 7 )
          type = 92212; // qc -> X
        else if( F2 == 9 )
          type = 92222; // qb -> X
        else
          error = true;
      }
      else
        error = true; // no elastic interactions among heavy quarks
    }
  }
  else if (( F1 % 2 ) == 0 && ( F2 % 2 ) == 0 )  // qbarqbar -> X, qbar+qbar' -> X
  {
    if ( F1 == F2 )  // qbarqbar -> X
    {
      if ( F1 > 6 )
        error = true; // no elastic interactions among heavy quarks
      else
        type = 9226;  // qbarqbar -> X
    }
    else  // qbar+qbar' -> X
    {
      if( F1 <= 6 )
      {
        if( F2 <= 6 )
          type = 9229;  // qbar+qbar' -> X
        else if( F2 == 8 )
          type = 92215; // qbar+cbar -> X
        else if( F2 == 10 )
          type = 92225; // qbar+bbar -> X
        else
          error = true;
      }
      else
        error = true; // no elastic interactions among heavy quarks
    }
  }
  else if ((( F1 % 2 ) == 1 && ( F2 % 2 ) == 0 ) && ( (F2 - F1) == 1 ) )  // qqbar -> X
  {
    if( F1 <= 6 )
      type = 9224; // qqbar -> X
    else if( F1 == 7 )
      type = 92216; // ccbar -> X
    else if( F1 == 9 )
      type = 92226; // bbbar -> X
    else
      error = true;
  }
  else if ( (( F1 % 2 ) == 1 && ( F2 % 2 ) == 0 ) ||  (( F2 % 2 ) == 1 && ( F1 % 2 ) == 0 ) )  //qqbar' -> X
  {
    if( F1 <= 6 )
    {
      if( F2 <= 6 )
        type = 9228; //qqbar' -> X
      else if( F2 == 7 )
        type = 92214; // qbar+c -> X
      else if( F2 == 8 )
        type = 92213; // q+cbar -> X
      else if( F2 == 9 )
        type = 92224; // qbar+b -> X
      else if( F2 == 10 )
        type = 92223; // q+bbar -> X
    }
    else
      error = true;     
  }
  else
    error = true;
  
  
  if( error )
  {
    std::string errMsg = "Inclusive process type could not be determined. Unrecoverable error.";
    std::cout << "error: " << errMsg << std::endl;
    throw eInteraction_type_error( errMsg );
  }

  //fix me: ugly!
  if ( genericType == c23 )
  {
    if( F2 <= 6 )
      type += 10;  // 9221 becomes 9231 etc.
    else if( F2 <= 10 )
      type += 100;  // 92210 becomes 92310 etc.
    else
    {
      std::string errMsg = "No 2 -> 3 processes involving a non-partonic particle.";
      throw eInteraction_type_error( errMsg );
    }
  }

  return type;
}



/**
* Given the flavors of the incoming particles, return the process type ID associated with this inclusive process (for 3->2)
*
* @param[in] _F1 Flavor of incoming particle 1
* @param[in] _F2 Flavor of incoming particle 2
* @param[in] _F3 Flavor of incoming particle 3
* @param[in] genericType generic collision type (c22, c23, c32) the specific process belongs to
* @return inclusive process type
*/
int interactionType::getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, const GENERIC_COLL_TYPE genericType )
{
  if ( genericType == c23 || genericType == c22 )
  {
    std::string errMsg = "Generic interaction type is 2->X, but three initial flavors are given.";
    throw eInteraction_type_error( errMsg );
  }
  
  int type = 0;

  FLAVOR_TYPE __F1, __F2;
  if ( _F1 == gluon )
  {
    __F1 = _F2;
    __F2 = _F3;
  }
  else if ( _F2 == gluon )
  {
    __F1 = _F1;
    __F2 = _F3;
  }
  else if ( _F3 == gluon )
  {
    __F1 = _F1;
    __F2 = _F2;
  }
  else
  {
    std::string errMsg = "No gluon in 3->2 process. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
  }

  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  //  unsigned int F1 = std::min( static_cast<unsigned int>( __F1 ), static_cast<unsigned int>( __F2 ) );
  unsigned int F2 = std::max( static_cast<unsigned int>( __F1 ), static_cast<unsigned int>( __F2 ) );
  
  // interaction types of 3->2 are similar to 2->2, only with one more gluon in the initial state
  type = getInclusiveProcessType( __F1, __F2, c22 );

  //fix me: ugly!
  if( F2 <= 6 )
    type += 100;  // 9221 becomes 9321 etc.
  else if( F2 <= 10 )
    type += 1000;  // 92210 becomes 93210 etc.
  else
  {
    std::string errMsg = "No 3 -> 2 processes involving a non-partonic particle.";
    throw eInteraction_type_error( errMsg );
  }

  return type;
}



/**
 * Get generic collision type (c22, c23, c32) for a specific collision type
 *
 * @param[in] _id process id to be looked up (231, 32800 etc.)
 * @return generic collision type (c22, c23, c32)
 */
GENERIC_COLL_TYPE interactionType::getCollType( const int _id )
{
  int number = _id;
  std::vector<unsigned short int> digits;

  while ( number > 0 )
  {
    digits.push_back( number % 10 );
    number /= 10;
  }
  int _typ = digits[ digits.size() - 1 ] * 10 + digits[ digits.size() - 2 ];

  switch (_typ)
  {
    case 22:
      return c22;
      break;
    case 23:
      return c23;
      break;
    default:
      return c32;
      break;
  }
}


double interactionType::getInvolvedInitial( const FLAVOR_TYPE _F, const int _id )
{
  if ( _F == gluon )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 221 :
      return 2;
      break;
    case 222 :
      return 2;
      break;
    case 223 :
      return 1;
      break;
    case 2230 :
      return 1;
      break;
      // charm sector
    case 2210 :
      return 2;
      break;
    case 2214 :
      return 1;
      break;
    case 22140 :
      return 1;
      break;
      // bottom sector
    case 2220 :
      return 2;
      break;
    case 2224 :
      return 1;
      break;
    case 22240 :
      return 1;
      break;
      // psi sector
    case 2240 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 231 :
      return 2;
      break;
    case 232 :
      return 2;
      break;
    case 233 :
      return 1;
      break;
    case 2330 :
      return 1;
      break;
      // charm sector
    case 2310 :
      return 2;
      break;
    case 2314 :
      return 1;
      break;
    case 23140 :
      return 1;
      break;
      // bottom sector
    case 2320 :
      return 2;
      break;
    case 2324 :
      return 1;
      break;
    case 23240 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 321 :
      return 3;
      break;
    case 322 :
      return 3;
      break;
    case 323 :
      return 2;
      break;
    case 3230 :
      return 2;
      break;
    case 324 :
      return 1;
      break;
    case 325 :
      return 1;
      break;
    case 326:
      return 1;
      break;
    case 327 :
      return 1;
      break;
    case 3270 :
      return 1;
      break;
    case 328 :
      return 1;
      break;
    case 3280 :
      return 1;
      break;
    case 32800 :
      return 1;
      break;
      // charm sector
    case 3210 :
      return 3;
      break;
    case 3211 :
      return 1;
      break;
    case 3212 :
      return 1;
      break;
    case 3213 :
      return 1;
      break;
    case 3214 :
      return 2;
      break;
    case 32140 :
      return 2;
      break;
    case 3215 :
      return 1;
      break;
    case 32150 :
      return 1;
      break;
    case 321500 :
      return 1;
      break;
    case 3215000 :
      return 1;
      break;
      // bottom sector
    case 3220 :
      return 3;
      break;
    case 3221 :
      return 1;
      break;
    case 3222 :
      return 1;
      break;
    case 3223 :
      return 1;
      break;
    case 3224 :
      return 2;
      break;
    case 32240 :
      return 2;
      break;
    case 3225 :
      return 1;
      break;
    case 32250 :
      return 1;
      break;
    case 322500 :
      return 1;
      break;
    case 3225000 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == light_quark || _F == up || _F == down || _F == strange )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 223 :
      return 1;
      break;
    case 224 :
      return 1;
      break;
    case 225 :
      return 1;
      break;
    case 226 :
      return 1;
      break;
    case 227 :
      return 2;
      break;
    case 228 :
      return 2;
      break;
    case 2280:
      return 1;
      break;
      // charm sector
    case 2211 :
      return 1;
      break;
    case 2215 :
      return 1;
      break;
    case 22150 :
      return 1;
      break;
      // bottom sector
    case 2221 :
      return 1;
      break;
    case 2225 :
      return 1;
      break;
    case 22250 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 233 :
      return 1;
      break;
    case 234 :
      return 1;
      break;
    case 235 :
      return 1;
      break;
    case 236 :
      return 1;
      break;
    case 237 :
      return 2;
      break;
    case 238 :
      return 2;
      break;
    case 2380:
      return 1;
      break;
      // charm sector
    case 2311 :
      return 1;
      break;
    case 2315 :
      return 1;
      break;
    case 23150 :
      return 1;
      break;
      // bottom sector
    case 2321 :
      return 1;
      break;
    case 2325 :
      return 1;
      break;
    case 23250 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 323 :
      return 1;
      break;
    case 324 :
      return 1;
      break;
    case 325 :
      return 1;
      break;
    case 326 :
      return 1;
      break;
    case 327 :
      return 2;
      break;
    case 328 :
      return 2;
      break;
    case 3280:
      return 1;
      break;
      // charm sector
    case 3211 :
      return 1;
      break;
    case 3215 :
      return 1;
      break;
    case 32150 :
      return 1;
      break;
      // bottom sector
    case 3221 :
      return 1;
      break;
    case 3225 :
      return 1;
      break;
    case 32250 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == anti_light_quark || _F == anti_up || _F == anti_down || _F == anti_strange )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 2230 :
      return 1;
      break;
    case 224 :
      return 1;
      break;
    case 225 :
      return 1;
      break;
    case 226 :
      return 1;
      break;
    case 2270 :
      return 2;
      break;
    case 2280:
      return 1;
      break;
    case 22800:
      return 2;
      break;
      // charm sector
    case 2211 :
      return 1;
      break;
    case 221500 :
      return 1;
      break;
    case 2215000 :
      return 1;
      break;
      // bottom sector
    case 2221 :
      return 1;
      break;
    case 222500 :
      return 1;
      break;
    case 2225000 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 2330 :
      return 1;
      break;
    case 234 :
      return 1;
      break;
    case 235 :
      return 1;
      break;
    case 236 :
      return 1;
      break;
    case 2370 :
      return 2;
      break;
    case 2380:
      return 1;
      break;
    case 23800:
      return 2;
      break;
      // charm sector
    case 2311 :
      return 1;
      break;
    case 231500 :
      return 1;
      break;
    case 2315000 :
      return 1;
      break;
      // bottom sector
    case 2321 :
      return 1;
      break;
    case 232500 :
      return 1;
      break;
    case 2325000 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 3230 :
      return 1;
      break;
    case 324 :
      return 1;
      break;
    case 325 :
      return 1;
      break;
    case 326 :
      return 1;
      break;
    case 3270 :
      return 2;
      break;
    case 3280:
      return 1;
      break;
    case 32800:
      return 2;
      break;
      // charm sector
    case 3211 :
      return 1;
      break;
    case 321500 :
      return 1;
      break;
    case 3215000 :
      return 1;
      break;
      // bottom sector
    case 3221 :
      return 1;
      break;
    case 322500 :
      return 1;
      break;
    case 3225000 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == charm )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 2212 :
      return 1;
      break;
    case 2213 :
      return 1;
      break;
    case 2214 :
      return 1;
      break;
    case 2215 :
      return 1;
      break;
    case 221500 :
      return 1;
      break;
    case 2241 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 2312 :
      return 1;
      break;
    case 2313 :
      return 1;
      break;
    case 2314 :
      return 1;
      break;
    case 2315 :
      return 1;
      break;
    case 231500 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 3212 :
      return 1;
      break;
    case 3213 :
      return 1;
      break;
    case 3214 :
      return 1;
      break;
    case 3215 :
      return 1;
      break;
    case 321500 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == anti_charm )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 2212 :
      return 1;
      break;
    case 2213 :
      return 1;
      break;
    case 22140 :
      return 1;
      break;
    case 22150 :
      return 1;
      break;
    case 2215000 :
      return 1;
      break;
    case 2241 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 2312 :
      return 1;
      break;
    case 2313 :
      return 1;
      break;
    case 23140 :
      return 1;
      break;
    case 23150 :
      return 1;
      break;
    case 2315000 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 3212 :
      return 1;
      break;
    case 3213 :
      return 1;
      break;
    case 32140 :
      return 1;
      break;
    case 32150 :
      return 1;
      break;
    case 3215000 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == bottom )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 2222 :
      return 1;
      break;
    case 2223 :
      return 1;
      break;
    case 2224 :
      return 1;
      break;
    case 2225 :
      return 1;
      break;
    case 222500 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 2322 :
      return 1;
      break;
    case 2323 :
      return 1;
      break;
    case 2324 :
      return 1;
      break;
    case 2325 :
      return 1;
      break;
    case 232500 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 3222 :
      return 1;
      break;
    case 3223 :
      return 1;
      break;
    case 3224 :
      return 1;
      break;
    case 3225 :
      return 1;
      break;
    case 322500 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == anti_bottom )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 2222 :
      return 1;
      break;
    case 2223 :
      return 1;
      break;
    case 22240 :
      return 1;
      break;
    case 22250 :
      return 1;
      break;
    case 2225000 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 2322 :
      return 1;
      break;
    case 2323 :
      return 1;
      break;
    case 23240 :
      return 1;
      break;
    case 23250 :
      return 1;
      break;
    case 2325000 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 3222 :
      return 1;
      break;
    case 3223 :
      return 1;
      break;
    case 32240 :
      return 1;
      break;
    case 32250 :
      return 1;
      break;
    case 3225000 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else
  {
    return -1;
  }

}



double interactionType::getInvolvedInitialInclusive( const FLAVOR_TYPE _F, const int _id )
{
  std::map<int, int>::const_iterator find_iter;
  
  FLAVOR_TYPE F_gen = ParticlePrototype::mapToGenericFlavorType( _F );

  switch ( F_gen )
  {
  case gluon:
    find_iter = involvedInitialGluons.find( _id );
    if ( find_iter != involvedInitialGluons.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of inclusive process index for gluon failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case light_quark:
    find_iter = involvedInitialLightQuarks.find( _id );
    if ( find_iter != involvedInitialLightQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of inclusive process index for light quarks failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case anti_light_quark:
    find_iter = involvedInitialAntiLightQuarks.find( _id );
    if ( find_iter != involvedInitialAntiLightQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of inclusive process index for light anti-quark failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case charm:
    find_iter = involvedInitialCharmQuarks.find( _id );
    if ( find_iter != involvedInitialCharmQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of inclusive process index for charm quark failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case anti_charm:
    find_iter = involvedInitialAntiCharmQuarks.find( _id );
    if ( find_iter != involvedInitialAntiCharmQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of inclusive process index for charm anti-quark failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case bottom:
    find_iter = involvedInitialBottomQuarks.find( _id );
    if ( find_iter != involvedInitialBottomQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of inclusive process index for bottom quark failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case anti_bottom:
    find_iter = involvedInitialAntiBottomQuarks.find( _id );
    if ( find_iter != involvedInitialAntiBottomQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of inclusive process index for bottom anti-quark failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  default:
    std::string errMsg = "Determination of inclusive process index failed. Wrong flavor. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}


std::vector<int> interactionType::getInteractionTypes( const FLAVOR_TYPE _F )
{
  std::vector<int> _types;
  _types.reserve( 20 );

  if ( _F == gluon )
  {
    int mytypes[] = { 221, 222, 223, 2230, 2210, 2214, 22140, 2220, 2224, 22240, 2240,
                      231, 232, 233, 2330, 2310, 2314, 23140, 2320, 2324, 23240,
                      321, 322, 323, 3230, 324, 325, 326, 327, 3270, 328, 3280, 32800, 3210, 3211, 3212, 3213, 3214, 32140, 3215, 32150, 321500, 3215000, 3220, 3221, 3222, 3223, 3224, 32240, 3225, 32250, 322500, 3225000 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == light_quark || _F == up || _F == down || _F == strange )
  {
    int mytypes[] = { 223, 224, 225, 226, 227, 228, 2280, 2211, 2215, 22150, 2221, 2225, 22250, 
                      233, 234, 235, 236, 237, 238, 2380, 2311, 2315, 23150, 2321, 2325, 23250, 
                      323, 324, 325, 326, 327, 328, 3280, 3211, 3215, 32150, 3221, 3225, 32250 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == anti_light_quark || _F == anti_up || _F == anti_down || _F == anti_strange )
  {
    int mytypes[] = { 2230, 224, 225, 226, 2270, 2280, 22800, 2211, 221500, 2215000, 2221, 222500, 2225000, 
                      2330, 234, 235, 236, 2370, 2380, 23800, 2311, 231500, 2315000, 2321, 232500, 2325000, 
                      3230, 324, 325, 326, 3270, 3280, 32800, 3211, 321500, 3215000, 3221, 322500, 3225000 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == charm )
  {
    int mytypes[] = { 2212, 2213, 2214, 2215, 221500, 2241,
                      2312, 2313, 2314, 2315, 231500, 
                      3212, 3213, 3214, 3215, 321500 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == anti_charm )
  {
    int mytypes[] = { 2212, 2213, 22140, 22150, 2215000, 2241,
                      2312, 2313, 23140, 23150, 2315000, 
                      3212, 3213, 32140, 32150, 3215000 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == bottom )
  {
    int mytypes[] = { 2222, 2223, 2224, 2225, 222500, 
                      2322, 2323, 2324, 2325, 232500, 
                      3222, 3223, 3224, 3225, 322500 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == anti_bottom )
  {
    int mytypes[] = { 2222, 2223, 22240, 22250, 2225000, 
                      2322, 2323, 23240, 23250, 2325000, 
                      3222, 3223, 32240, 32250, 3225000 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }

  return _types;
}




// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
