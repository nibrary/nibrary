include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# CMAKE_INSTALL_PREFIX is NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX

# Create directories
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}")
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/libzip")

# Rename headers
conditional_copy_file("${CMAKE_INSTALL_PREFIX}/include/zipconf.h" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/libzip/zipconf.h")
conditional_copy_file("${CMAKE_INSTALL_PREFIX}/include/zip.h" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/libzip/zip.h")

# Rename libraries
if (NOT BUILD_SHARED_LIBS)
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.a" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.a")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/zipstatic.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zipstatic.lib")

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.so")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.so.1")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.so.${LIBZIP_MIN_VERSION}")

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.dylib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.1.dylib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.${LIBZIP_MIN_VERSION}.dylib")

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/zip.lib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/zip1.dll")
else()
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.so" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.so")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.so.1" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.so.1")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.so.${LIBZIP_MIN_VERSION}" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.so.${LIBZIP_MIN_VERSION}")

    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.dylib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.1.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.1.dylib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libzip.${LIBZIP_MIN_VERSION}.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.${LIBZIP_MIN_VERSION}.dylib")

    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/zip.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zip.lib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/zip1.dll" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zip1.dll")
endif()

