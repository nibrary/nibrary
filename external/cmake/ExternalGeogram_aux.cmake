include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# CMAKE_INSTALL_PREFIX is NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX

# Create directories
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}")
conditional_make_directory("${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/geogram")

# Rename headers
conditional_make_directory("${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram")
conditional_copy_directory("${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram" "${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/geogram")


# Rename libraries
if (NOT BUILD_SHARED_LIBS)
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libgeogram.a" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.a")
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/geogram.lib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/geogram.lib")
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/geogram.dll" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/geogram.dll")
else()
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.so")
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.so.1")
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.${GEOGRAM_MIN_VERSION}" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.so.${GEOGRAM_MIN_VERSION}")

    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libgeogram.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.dylib")
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libgeogram.1.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.1.dylib")
    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/libgeogram.${GEOGRAM_MIN_VERSION}.dylib" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.${GEOGRAM_MIN_VERSION}.dylib")

    conditional_rename("${CMAKE_INSTALL_PREFIX}/lib/geogram.dll" "${NIBRARY_CMAKE_INSTALL_PREFIX}/lib/${nibrary}/geogram.dll")
endif()

