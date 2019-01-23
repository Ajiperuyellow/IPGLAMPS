# This generates an up-to-date file 'revision.h' at runtime
#
# You may give the following defines as input:
# * MAIN_DIR --- The directory, where 'src' and 'tests' reside
# * SOURCE_DIR -- The source directory, probably '${MAIN_DIR}/src'
# * LIB_SOURCE_DIR -- The library's source dir, probably '${MAIN_DIR}/lib'


# the FindSubversion.cmake module is part of the standard distribution
INCLUDE(FindSubversion)

# get the physical path of the library dir (without symlink):
execute_process(
  COMMAND pwd -P
  WORKING_DIRECTORY ${LIB_SOURCE_DIR} 
  OUTPUT_VARIABLE LIB_SOURCE_DIR_ABS )

# remove trailing newlines:
STRING(REGEX REPLACE "(\r?\n)+$" "" LIB_SOURCE_DIR_ABS "${LIB_SOURCE_DIR_ABS}")

# extract working copy information for SOURCE_DIR into MY_XXX variables
IF ( NOT DISABLE_SVN_EXTRAS )
  Subversion_WC_INFO(${MAIN_DIR} MY)
  MESSAGE(STATUS "svn revision: ${MY_WC_REVISION}")
  MESSAGE(STATUS "running svn command: ${Subversion_SVN_EXECUTABLE} status -q ${MAIN_DIR}")
  execute_process(COMMAND ${Subversion_SVN_EXECUTABLE} status -q ${MAIN_DIR} OUTPUT_VARIABLE SVNSTATUS)

  SET( SVNSTATUSLIB "no info available" )
  IF( NOT LINK_TO_BAMPS_LIB )
    MESSAGE(STATUS "running svn command: ${Subversion_SVN_EXECUTABLE} status -q ${LIB_SOURCE_DIR_ABS}")
    execute_process(COMMAND ${Subversion_SVN_EXECUTABLE} status -q ${LIB_SOURCE_DIR_ABS} OUTPUT_VARIABLE SVNSTATUSLIB)
  ENDIF()

  # in the following line, the line-break is intended!
  MESSAGE(STATUS "svn status -q: (SRC)
${SVNSTATUS}")
  # in the following line, the line-break is intended!
  MESSAGE(STATUS "svn status -q: (LIB)
${SVNSTATUSLIB}")

ELSE ()
  SET( MY_WC_REVISION "\"no svn revision info available\"" )
  SET( SVNSTATUS "no info available" )
  SET( SVNSTATUSLIB "no info available" )

ENDIF()


# Prepare the 'svn status' output for C++:
STRING(REGEX REPLACE "\n" "\\\\n# " C_SVNSTATUS "${SVNSTATUS}")
STRING(REGEX REPLACE "\n" "\\\\n# " C_SVNSTATUSLIB "${SVNSTATUSLIB}")

# write a temporary header file with the SVNVERSION define
FILE(WRITE revision.h.temp "#ifndef REVISION_H
#define REVISION_H

#define SVN_REVISION ${MY_WC_REVISION}

#define SVN_STATUS_SRC \"# ${C_SVNSTATUS}\"

#define SVN_STATUS_LIB \"# ${C_SVNSTATUSLIB}\"

#endif")

# copy the file to the final header only if the version changes
# reduces needless rebuilds
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different revision.h.temp revision.h)