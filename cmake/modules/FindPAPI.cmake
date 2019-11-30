# Module that checks whether the PAPI performance measurement library
# is available or not.
#
# Sets the following variables:
# HAVE_PAPI
# PAPI_LIBRARIES
#
# perform tests
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
include(CheckCXXCompilerFlag)

find_path (PAPI_INCLUDE_DIR
  NAMES "papi.h"
  HINTS "${PAPI_DIR}"
  PATHS "${PROJECT_SOURCE_DIR}"
  PATH_SUFFIXES "include" "src"
  DOC "Path to papi library header files"
  ${_no_default_path} )

find_library (PAPI_LIBRARY
  NAMES "libpapi.a"
  HINTS "${PAPI_DIR}"
  PATH_SUFFIXES "src" "lib" "lib${_BITS}" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
  DOC "Path to papi library archive/shared object files"
  ${_no_default_path} )

set(PAPI_INCLUDE_DIRS ${PAPI_INCLUDE_DIR})
set(PAPI_LIBRARIES ${PAPI_LIBRARY} ${PAPI_PFM_LIBRARY})

# setup list of all required libraries to link with papi
if(PAPI_INCLUDE_DIRS AND PAPI_LIBRARIES)
  cmake_push_check_state(RESET)
  set (CMAKE_REQUIRED_INCLUDES ${PAPI_INCLUDE_DIRS})
  set (CMAKE_REQUIRED_LIBRARIES ${PAPI_LIBRARIES})
  CHECK_CXX_SOURCE_COMPILES("
#include <papi.h>

int main(void){
   PAPI_library_init(PAPI_VER_CURRENT);
   return 0;
}" PAPI_FOUND)
  cmake_pop_check_state()

  if (PAPI_FOUND)
    set(HAVE_PAPI "${PAPI_FOUND}")
    set (CMAKE_REQUIRED_INCLUDES ${PAPI_INCLUDE_DIRS})
    set (CMAKE_REQUIRED_LIBRARIES ${PAPI_LIBRARIES})
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Papi
  DEFAULT_MSG
  PAPI_INCLUDE_DIRS
  PAPI_LIBRARIES
  HAVE_PAPI)
