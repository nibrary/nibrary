include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# CMAKE_INSTALL_PREFIX passed here is NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX (temp install dir)
# NIBRARY_CMAKE_INSTALL_PREFIX is layout install dir

# Create directories
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}")
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/trx-cpp/trx")

# Rename/Copy headers
conditional_copy_directory("${CMAKE_INSTALL_PREFIX}/include/trx" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/trx-cpp/trx")
conditional_copy_directory("${CMAKE_INSTALL_PREFIX}/include/mio" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/trx-cpp/mio")
conditional_copy_file("${CMAKE_INSTALL_PREFIX}/include/json11.hpp" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/trx-cpp/json11.hpp")

# Copy libraries
if (NOT BUILD_SHARED_LIBS)
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libtrx.a" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libtrx.a")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/trx.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/trx.lib")
else()
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libtrx.so" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libtrx.so")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libtrx.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libtrx.dylib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/trx.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/trx.lib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/trx.dll" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/trx.dll")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/bin/trx.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/trx.lib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/bin/trx.dll" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/trx.dll")
endif()
