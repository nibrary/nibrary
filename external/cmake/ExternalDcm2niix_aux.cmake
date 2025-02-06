include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# CMAKE_INSTALL_PREFIX is NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX

# Create directories
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}")
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/dcm2niix")

# Rename headers
conditional_copy_file("${NIBRARY_CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix++.h" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/dcm2niix/dcm2niix++.h")

# Rename libraries
if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    conditional_copy_file("${CMAKE_BINARY_DIR}/dcm2niix/src/build_dcm2niix-build/lib/dcm2niix++.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/dcm2niix++.lib")
else()
    conditional_copy_file("${CMAKE_BINARY_DIR}/dcm2niix/src/build_dcm2niix-build/lib/libdcm2niix++.a" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libdcm2niix++.a")
endif()


