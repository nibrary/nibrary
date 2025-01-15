list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(NIBRARY_CMAKE_SOURCE_DIR        ${CMAKE_SOURCE_DIR}     CACHE INTERNAL "")
set(NIBRARY_CMAKE_INSTALL_PREFIX    ${CMAKE_INSTALL_PREFIX} CACHE INTERNAL "")

set(NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/nibrary_install_stash CACHE INTERNAL "")

include(ExternalEigen)
include(ExternalZlib)
include(ExternalGeogram)
include(ExternalLibigl)
include(ExternalSimde)
include(ExternalProxsuite)
include(ExternalDcm2niix)
