# .. cmake_module::
#
# This module's content is executed whenever a Dune module requires or
# suggests ewoms-numerics!
#

# PAPI this is unused in the vanilla code, but when performance
# measurement it is handy to specify PAPI probes without switching the
# branch of the build system repository.
find_package(PAPI)

if(PAPI_FOUND)
  set(HAVE_PAPI 1)
  dune_register_package_flags(
    INCLUDE_DIRS "${PAPI_INCLUDE_DIR}"
    LIBRARIES "${PAPI_LIBRARIES}")
endif()

# the HWLOC library is also unused by vanilla ewoms-numerics. Like PAPI, it is
# useful for debugging purposes.
find_package(Hwloc)

if(Hwloc_FOUND)
  set(HAVE_HWLOC 1)
  dune_register_package_flags(
    INCLUDE_DIRS "${Hwloc_INCLUDE_DIRS}"
    LIBRARIES "${Hwloc_LIBRARIES}")
endif()
