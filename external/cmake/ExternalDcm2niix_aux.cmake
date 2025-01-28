include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# CMAKE_INSTALL_PREFIX is NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX

# Create directories
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}")
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/dcm2niix")

# Rename headers
conditional_copy_file("${NIBRARY_CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix_fswrapper.h" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/dcm2niix/dcm2niix_fswrapper.h")

# Rename libraries
conditional_copy_file("${CMAKE_BINARY_DIR}/dcm2niix/src/build_dcm2niix-build/lib/libdcm2niixfs.a" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libdcm2niixfs.a")
conditional_copy_file("${CMAKE_BINARY_DIR}/dcm2niix/src/build_dcm2niix-build/lib/libdcm2niixfs.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libdcm2niixfs.lib")


