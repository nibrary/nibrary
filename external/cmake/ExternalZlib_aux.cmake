include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# Create directories
conditional_make_directory("${CMAKE_INSTALL_PREFIX}/include/${nibrary}/zlib")
conditional_make_directory("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}")

# Rename headers
conditional_rename("${CMAKE_INSTALL_PREFIX}/include/zconf.h" "${CMAKE_INSTALL_PREFIX}/include/${nibrary}/zlib/zconf.h")
conditional_rename("${CMAKE_INSTALL_PREFIX}/include/zlib.h" "${CMAKE_INSTALL_PREFIX}/include/${nibrary}/zlib/zlib.h")

# Rename libraries
conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libz.a" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.a")

conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libz.so" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so")
conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libz.so.1" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so.1")
conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libz.so.${ZLIB_MIN_VERSION}" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so.${ZLIB_MIN_VERSION}")

conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libz.dylib" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.dylib")
conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libz.1.dylib" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.1.dylib")
conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libz.${ZLIB_MIN_VERSION}.dylib" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.${ZLIB_MIN_VERSION}.dylib")

conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/zlib.lib" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlib.lib")
conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/zlib1.dll" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlib1.dll")
conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/zlibstatic.lib" "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlibstatic.lib")

# Remove libraries
conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.a")

conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so")
conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so.1")
conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.so.${ZLIB_MIN_VERSION}")

conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.dylib")
conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.1.dylib")
conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/libz.${ZLIB_MIN_VERSION}.dylib")

conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/zlib.lib")
conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/zlib1.dll")
conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/zlibstatic.lib")


# Removed shared libraries if BUILD_SHARED_LIBS is OFF
if (NOT BUILD_SHARED_LIBS)

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so.1")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so.${ZLIB_MIN_VERSION}")

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.dylib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.1.dylib")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.${ZLIB_MIN_VERSION}.dylib")

else()

    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.a")
    conditional_remove_file("${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlibstatic.lib")

endif()
