include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# CMAKE_INSTALL_PREFIX is NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX

# Create directories
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}")
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/zlib")

# Rename headers
conditional_copy_file("${CMAKE_INSTALL_PREFIX}/include/zconf.h" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/zlib/zconf.h")
conditional_copy_file("${CMAKE_INSTALL_PREFIX}/include/zlib.h" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/zlib/zlib.h")

conditional_copy_file("${CMAKE_INSTALL_PREFIX}/include/zconf.h" "${CMAKE_INSTALL_PREFIX}/zconf.h")
conditional_copy_file("${CMAKE_INSTALL_PREFIX}/include/zlib.h" "${CMAKE_INSTALL_PREFIX}/zlib.h")

# Rename libraries
if (NOT BUILD_SHARED_LIBS)
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libz.a" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.a")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/zlibstatic.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlibstatic.lib")

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so.1")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so.${ZLIB_MIN_VERSION}")

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.dylib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.1.dylib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.${ZLIB_MIN_VERSION}.dylib")

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/zlib.lib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/zlib1.dll")
else()
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so.1" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so.1")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so.${ZLIB_MIN_VERSION}" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so.${ZLIB_MIN_VERSION}")

    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libz.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.dylib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libz.1.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.1.dylib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/libz.${ZLIB_MIN_VERSION}.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.${ZLIB_MIN_VERSION}.dylib")

    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/zlib.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlib.lib")
    conditional_copy_file("${CMAKE_INSTALL_PREFIX}/lib/zlib1.dll" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlib1.dll")
endif()

