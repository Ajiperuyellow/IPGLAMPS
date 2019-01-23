
SET( BAMPS_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/src/" "${CMAKE_CURRENT_BINARY_DIR}/../src/" )
SET( BAMPS_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/../src/" )
SET( BAMPS_LIBRARIES -lbamps )

INCLUDE_DIRECTORIES( ${BAMPS_INCLUDE_DIRS} )

macro(bamps_add_prog _name)
   set(name ${_name})
   add_executable(${_name} ${_name}.cpp)
   target_link_libraries(${_name} bamps 
      ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
   IF( CMAKE_SYSTEM_NAME STREQUAL "Linux" ) 
      target_link_libraries( ${_name} rt )
   ENDIF()
   IF ( Cuba_FOUND )
     TARGET_LINK_LIBRARIES( ${_name} ${Cuba_LIBRARIES} )
   ENDIF ( Cuba_FOUND )
#   add_test(${_name} "${CMAKE_CURRENT_BINARY_DIR}/${_name}")
endmacro(bamps_add_prog)

bamps_add_prog(xx_VectorTiming)
bamps_add_prog(xx_OpenMP_TaskReduction)

######### Copy additional files #########

function(copy_files)
  foreach (filename ${ARGV})
#    MESSAGE(STATUS "copying: ${filename}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename})
  endforeach (filename)
endfunction(copy_files)

copy_files( testNf3_32configurations.dat )
