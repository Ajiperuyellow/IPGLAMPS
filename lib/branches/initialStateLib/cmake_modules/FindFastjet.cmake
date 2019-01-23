# Try to find the Fastjet library
#
# Florian Senzel, Kai Gallmeister ( based on FindCuba.cmake by Oliver Fochler)
#
# This package defines
#  Fastjet_FOUND - The fastjet library has been found
#  Fastjet_INCLUDE_DIRS - The directory in which the Fastjet headers reside
#  Fastjet_LIBRARIES - The Fastjet library
#  CGAL_FOUND - The CGAL library has been found
# in addition for the plugins:
#  FASTJET_SISCONE_FOUND - the siscone plugin has been found
#  FASTJET_SISCONE_SPHERICAL_FOUND - the siscone_sphercal plugin has been found
#  FASTJET_CMSITERATIVECONE_FOUND - the CMSIterativeCone plugin has been found

INCLUDE(FindPackageHandleStandardArgs)

# find the path of the include files, based on two headers
FIND_PATH( Fastjet_INCLUDE_DIR 
    NAMES PseudoJet.hh
    PATHS $ENV{HOME}/usr/include/    # suggest a user based include tree
    PATH_SUFFIXES fastjet FASTJET Fastjet FastJet fastJet  # suggest some path suffixes in which the headers could be located
)

# find the CGAL library
FIND_LIBRARY( CGAL_LIBRARY
    NAMES CGAL
    PATHS $ENV{HOME}/usr/lib /usr/local/lib         # suggest a user based include tree
    PATH_SUFFIXES CGAL cgal   # suggest some path suffixes in which the headers could be located
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( CGAL DEFAULT_MSG CGAL_LIBRARY)

# find the GMP library
FIND_LIBRARY( GMP_LIBRARY
    NAMES gmp libgmp-10
    PATHS $ENV{HOME}/usr/lib /usr/local/lib         # suggest a user based include tree
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GMP DEFAULT_MSG GMP_LIBRARY)
        
# find the GMPXX library
FIND_LIBRARY( GMPXX_LIBRARY
    NAMES gmpxx
    PATHS $ENV{HOME}/usr/lib /usr/local/lib         # suggest a user based include tree
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GMPXX DEFAULT_MSG GMPXX_LIBRARY)
                
# find the MPFR library
FIND_LIBRARY( MPFR_LIBRARY
    NAMES mpfr libmpfr-4 libmpfr-1
    PATHS $ENV{HOME}/usr/lib /usr/local/lib         # suggest a user based include tree
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( MPFR DEFAULT_MSG MPFR_LIBRARY)
                        
# find the fastjet library
FIND_LIBRARY( Fastjet_LIBRARY
    NAMES fastjet
    PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
    PATH_SUFFIXES fastjet FASTJET Fastjet FastJet fastJet   # suggest some path suffixes in which the headers could be located
)

# find the fastjet plugins library
FIND_LIBRARY( Fastjet_Plugins_LIBRARY
    NAMES fastjetplugins
    PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
    PATH_SUFFIXES fastjet FASTJET Fastjet FastJet fastJet   # suggest some path suffixes in which the headers could be located
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Fastjet_Plugins DEFAULT_MSG Fastjet_Plugins_LIBRARY)

# find some fastjet plugins:
FIND_LIBRARY( siscone_LIBRARY
   NAMES siscone
   PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Fastjet_siscone DEFAULT_MSG siscone_LIBRARY)

FIND_LIBRARY( siscone_spherical_LIBRARY
   NAMES siscone_spherical
   PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Fastjet_siscone_spherical DEFAULT_MSG siscone_spherical_LIBRARY)

FIND_PATH( CMSIterativeCone_LIBRARY
    NAMES CMSIterativeConePlugin.hh
    PATHS $ENV{HOME}/usr/include          # suggest a user based include tree
    PATH_SUFFIXES fastjet FASTJET Fastjet FastJet fastJet  # suggest some path suffixes in which the headers could be located
)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Fastjet_CMSIterativeCone DEFAULT_MSG CMSIterativeCone_LIBRARY)

# adhere to the standard nomenclature for find_package routines and set some variables
SET( Fastjet_LIBRARIES ${Fastjet_LIBRARY} )
IF( CGAL_FOUND AND GMP_FOUND AND GMPXX_FOUND AND MPFR_FOUND )
    SET( Fastjet_LIBRARIES ${Fastjet_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY} ${MPFR_LIBRARY}) 
ENDIF()
SET( Fastjet_LIBRARIES ${Fastjet_LIBRARIES} ${Fastjet_Plugins_LIBRARY} )
IF( FASTJET_SISCONE_FOUND )
    SET( Fastjet_LIBRARIES ${Fastjet_LIBRARIES} ${siscone_LIBRARY} )
ENDIF()
IF( FASTJET_SISCONE_SPHERICAL_FOUND )
    SET( Fastjet_LIBRARIES ${Fastjet_LIBRARIES} ${siscone_spherical_LIBRARY} )
ENDIF()

SET( Fastjet_INCLUDE_DIRS ${Fastjet_INCLUDE_DIR} )

# handle the QUIETLY and REQUIRED arguments and set Fastjet_FOUND to TRUE
# if all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Fastjet DEFAULT_MSG Fastjet_LIBRARY Fastjet_Plugins_LIBRARY Fastjet_INCLUDE_DIR)


## display some status information
#SET ( Fastjet_FOUND Fastjet_FOUND CACHE INTERNAL "Provide Fastjet_FOUND in addition to Fastjet_FOUND" FORCE )
#SET ( siscone_FOUND siscone_FOUND CACHE INTERNAL "Provide siscone_FOUND in addition to Fastjet_FOUND" FORCE )
#IF( Fastjet_FOUND )
# MESSAGE( STATUS "** Fastjet sub-libraries:" )
# MESSAGE( STATUS "**     Plugins:           ${Fastjet_Plugins_LIBRARY}" )
# MESSAGE( STATUS "**     Tools:             ${Fastjet_Tools_LIBRARY}" )
# MESSAGE( STATUS "**     SIScone:           ${siscone_LIBRARY}" )
# MESSAGE( STATUS "**     SIScone Spherical: ${siscone_spherical_LIBRARY}" )
# MESSAGE( STATUS "**     CMSIterativeCone:  ${CMSIterativeCone_LIBRARY}" )
# MESSAGE( STATUS "**     CGAL:              ${CGAL_LIBRARY}" )
# 
# MESSAGE( STATUS "** Fastjet library: ${Fastjet_LIBRARIES}" )
# MESSAGE( STATUS "** Fastjet include: ${Fastjet_INCLUDE_DIRS}" )
#ENDIF( Fastjet_FOUND )

# the variables will only show up in the GUI in the "advanced" view
MARK_AS_ADVANCED(Fastjet_INCLUDE_DIR Fastjet_LIBRARY CGAL_FOUND GMP_FOUND GMPXX_FOUND MPFR_FOUND)
