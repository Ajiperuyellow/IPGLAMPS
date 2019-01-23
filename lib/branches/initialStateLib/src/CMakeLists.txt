SET(bamps_lib_SRCS
binary_cross_sections.cpp
binning.cpp
binning2.cpp
binning3.cpp
centralityclass.cpp
configurationbase.cpp
coupling.cpp
hydroParticleType.cpp
hydroTensor.cpp
integrand23.cpp
integrand32.cpp
integrate.cpp
interactiontype.cpp
interpolationGSL.cpp
interpolation_n_dimensions.cpp
interpolation23.cpp
interpolation22.cpp
iterate_mfp.cpp
logger.cpp
lorentz.cpp
mfp_data.cpp
particleprototype.cpp
random.cpp
rangen_distr.cpp
ratesmanager.cpp
rotation.cpp
scattering22.cpp
scattering22_hydro.cpp
scattering23.cpp
scattering32.cpp 
signalhandler.cpp
tauDistro.cpp
thermal.cpp
vegas.cpp
)

# include the output directory, where the svnversion.h file is generated
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

######### create a config.h #########
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/configBAMPS.h.in ${CMAKE_CURRENT_BINARY_DIR}/configBAMPS.h)

ADD_LIBRARY(bamps STATIC ${bamps_lib_SRCS})

### activate interprocedural optimization if the compiler supports it
IF( CMAKE_BUILD_TYPE STREQUAL "Release" )
  SET_PROPERTY(TARGET bamps PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)    # -ipo option for intel compilers 
ENDIF()

#TARGET_LINK_LIBRARIES(bamps ${Boost_LIBRARIES} ${Cuba_LIBRARIES})

IF( GENERATE_BAMPS_LIB_INSTALL_TARGETS )
  FILE(GLOB headerFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  INSTALL(FILES ${headerFiles} DESTINATION include/BAMPS)
  INSTALL(TARGETS bamps ARCHIVE DESTINATION lib)
ENDIF()

######### create CPPCHECK targets #########
INCLUDE(CppcheckTargets)
add_cppcheck(bamps STYLE)
