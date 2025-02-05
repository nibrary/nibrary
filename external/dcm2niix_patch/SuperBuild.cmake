# Check if git exists
find_package(Git)
if(NOT GIT_FOUND)
    message(FATAL_ERROR "Cannot find Git. Git is required for Superbuild")
endif()

include(${CMAKE_SOURCE_DIR}/cmake/dcm2niixInitializeBuildType.cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(USE_STATIC_RUNTIME "Use static runtime" ON)

if(USE_STATIC_RUNTIME)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        find_file(STATIC_LIBCXX "libstdc++.a" ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
        mark_as_advanced(STATIC_LIBCXX)
        if(NOT STATIC_LIBCXX)
            unset(STATIC_LIBCXX CACHE)
            # Only on some Centos/Redhat systems
            message(FATAL_ERROR
                "\"USE_STATIC_RUNTIME\" set to ON but \"libstdc++.a\" not found! Set it to OFF or \
                 \"yum install libstdc++-static\" to resolve the error.")
        endif()
    endif()
endif()

if(APPLE)
    set(OSX_ARCHITECTURES "-DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}")
endif()

option(USE_TURBOJPEG "Use TurboJPEG to decode classic JPEG" OFF)
option(USE_JASPER "Build with JPEG2000 support using Jasper" OFF)
set(USE_OPENJPEG "OFF" CACHE STRING "Build with JPEG2000 support using OpenJPEG.")
set_property(CACHE USE_OPENJPEG PROPERTY STRINGS  "OFF;GitHub;System;Custom")
option(USE_JPEGLS "Build with JPEG-LS support using CharLS" OFF)
option(USE_JNIFTI "Build with JNIFTI support" ON)

option(BATCH_VERSION "Build dcm2niibatch for multiple conversions" OFF)

option(BUILD_DCM2NIIXFSLIB "Build libdcm2niixfs.a" OFF)

if(BUILD_DCM2NIIXFSLIB)
    if(USE_OPENJPEG OR USE_TURBOJPEG OR USE_JASPER)
        message("-- Set BUILD_DCM2NIIXFSLIB to OFF since USE_TURBOJPEG/USE_JASPER/USE_OPENJPEG is ON.")
        set(BUILD_DCM2NIIXFSLIB OFF CACHE BOOL "Build libdcm2niixfs.a" FORCE)
    else()
        message("-- Build libdcm2niixfs.a: ${BUILD_DCM2NIIXFSLIB}")
    endif()
endif()

option(BUILD_DCM2NIIX_LIB "Build dcm2niix library" OFF)

if(BUILD_DCM2NIIX_LIB)
    if(USE_OPENJPEG OR USE_TURBOJPEG OR USE_JASPER)
        message("-- Set BUILD_DCM2NIIX_LIB to OFF since USE_TURBOJPEG/USE_JASPER/USE_OPENJPEG is ON.")
        set(BUILD_DCM2NIIX_LIB OFF CACHE BOOL "Build dcm2niix library" FORCE)
    else()
        message("-- Build dcm2niix library: ${BUILD_DCM2NIIX_LIB}")
    endif()
endif()


include(ExternalProject)

set(DEPENDENCIES)

option(INSTALL_DEPENDENCIES "Optionally install built dependent libraries (OpenJPEG and yaml-cpp) for future use." OFF)

if(INSTALL_DEPENDENCIES)
    set(DEP_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
else()
    set(DEP_INSTALL_DIR ${CMAKE_BINARY_DIR})
endif()

if(NOT ${USE_OPENJPEG} STREQUAL "OFF")
    message("-- Build with OpenJPEG: ${USE_OPENJPEG}")

    if(OpenJPEG_DIR)
        set(USE_OPENJPEG "Custom" CACHE STRING "Build with JPEG2000 support using OpenJPEG." FORCE)
        set(OpenJPEG_DIR "${OpenJPEG_DIR}" CACHE PATH "Path to OpenJPEG configuration file"  FORCE)
        message("--   Using OpenJPEG library from ${OpenJPEG_DIR}")
    elseif(${USE_OPENJPEG} STREQUAL "System")
        find_package(PkgConfig)
        if(PKG_CONFIG_FOUND)
            pkg_check_modules(OPENJPEG libopenjp2)
        endif()

        if(OPENJPEG_FOUND)
            if(NOT ${OPENJPEG_INCLUDE_DIRS} MATCHES "gdcmopenjpeg")
        	    set(OpenJPEG_DIR ${OPENJPEG_LIBDIR}/openjpeg-${OPENJPEG_VERSION} CACHE PATH "Path to OpenJPEG configuration file" FORCE)
                message("--   Using OpenJPEG library from ${OpenJPEG_DIR}")
            else()
                message("--   Unable to use GDCM's internal OpenJPEG")
            endif()
        endif()
    endif()

    if(${USE_OPENJPEG} STREQUAL "GitHub" OR NOT OpenJPEG_DIR)
        set(USE_OPENJPEG "GitHub" CACHE STRING "Build with JPEG2000 support using OpenJPEG." FORCE)
        include(${CMAKE_SOURCE_DIR}/SuperBuild/External-OPENJPEG.cmake)
        list(APPEND DEPENDENCIES openjpeg)
        set(BUILD_OPENJPEG TRUE)
        message("--   Will build OpenJPEG library from GitHub")
    endif()
endif()

if(BATCH_VERSION)
    message("-- Build dcm2niibatch: ${BATCH_VERSION}")

    set(YAML-CPP_IMPLEMENTATION "GitHub" CACHE STRING "Choose yaml-cpp implementation.")
    set_property(CACHE YAML-CPP_IMPLEMENTATION PROPERTY STRINGS  "GitHub;System;Custom")

    if(YAML-CPP_DIR)
        set(YAML-CPP_IMPLEMENTATION "Custom" CACHE STRING "Choose yaml-cpp implementation." FORCE)
        set(YAML-CPP_DIR ${YAML-CPP_DIR} CACHE PATH "Path to yaml-cpp configuration file"  FORCE)
        message("--   Using yaml-cpp library from ${YAML-CPP_DIR}")
    elseif(${YAML-CPP_IMPLEMENTATION} STREQUAL "System")
        find_package(PkgConfig)
        if(PKG_CONFIG_FOUND)
            pkg_check_modules(YAML-CPP yaml-cpp)
        endif()

        # Build from GitHub if not found or version < 0.5.3
        if(YAML-CPP_FOUND AND NOT (YAML-CPP_VERSION VERSION_LESS "0.5.3"))
            set(YAML-CPP_DIR ${YAML-CPP_LIBDIR}/cmake/yaml-cpp CACHE PATH "Path to yaml-cpp configuration file"  FORCE)
            message("--   Using yaml-cpp library from ${YAML-CPP_DIR}")
        endif()
    endif()

    if(${YAML-CPP_IMPLEMENTATION} STREQUAL "GitHub" OR NOT YAML-CPP_DIR)
        set(YAML-CPP_IMPLEMENTATION "GitHub" CACHE STRING "Choose yaml-cpp implementation." FORCE)
        include(${CMAKE_SOURCE_DIR}/SuperBuild/External-YAML-CPP.cmake)
        list(APPEND DEPENDENCIES yaml-cpp)
        set(BUILD_YAML-CPP TRUE)
        message("--   Will build yaml-cpp library from GitHub")
    endif()
endif()

set(ZLIB_IMPLEMENTATION "Miniz" CACHE STRING "Choose zlib implementation.")
set_property(CACHE ZLIB_IMPLEMENTATION PROPERTY STRINGS  "Miniz;System;Cloudflare;Custom")
if(${ZLIB_IMPLEMENTATION} STREQUAL "Cloudflare")
    message("-- Build with Cloudflare zlib: ON")
    include(${CMAKE_SOURCE_DIR}/SuperBuild/External-CLOUDFLARE-ZLIB.cmake)
    list(APPEND DEPENDENCIES zlib)
    set(BUILD_CLOUDFLARE-ZLIB TRUE)
    message("--   Will build Cloudflare zlib from GitHub")
elseif(${ZLIB_IMPLEMENTATION} STREQUAL "Custom")
    set(ZLIB_ROOT ${ZLIB_ROOT} CACHE PATH "Specify custom zlib root directory.")
    if(NOT ZLIB_ROOT)
        message(FATAL_ERROR "ZLIB_ROOT needs to be set to locate custom zlib!")
    endif()
endif()

set(MOD_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
set(MOD_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")

if (CMAKE_SYSTEM_NAME STREQUAL "Windows") 
    if(USE_STATIC_RUNTIME)
        set(MOD_C_FLAGS   "-fPIC $<$<CONFIG:Debug>:/MTd;/MT> ${CMAKE_C_FLAGS}")
        set(MOD_CXX_FLAGS "-fPIC $<$<CONFIG:Debug>:/MTd;/MT> ${CMAKE_CXX_FLAGS}")
    else()
        set(MOD_C_FLAGS   "-fPIC $<$<CONFIG:Debug>:/MDd;/MD> ${CMAKE_C_FLAGS}")
        set(MOD_CXX_FLAGS "-fPIC $<$<CONFIG:Debug>:/MDd;/MD> ${CMAKE_CXX_FLAGS}")
    endif()
endif()

ExternalProject_Add(console
    DEPENDS ${DEPENDENCIES}
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/console
    BINARY_DIR console-build
    CMAKE_ARGS
        -Wno-dev
        --no-warn-unused-cli
        ${EXTERNAL_PROJECT_BUILD_TYPE_CMAKE_ARGS}
        ${OSX_ARCHITECTURES}
        # Install directories
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}
        # Compiler settings
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_FLAGS=${MOD_C_FLAGS}
        -DCMAKE_CXX_FLAGS=${MOD_CXX_FLAGS}
        # Options
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=${CMAKE_VERBOSE_MAKEFILE}
        -DUSE_STATIC_RUNTIME:BOOL=${USE_STATIC_RUNTIME}
        -DUSE_TURBOJPEG:BOOL=${USE_TURBOJPEG}
        -DUSE_JASPER:BOOL=${USE_JASPER}
        -DUSE_JPEGLS:BOOL=${USE_JPEGLS}
        -DUSE_JNIFTI:BOOL=${USE_JNIFTI}
        # ZLIB
        -DZLIB_IMPLEMENTATION:STRING=${ZLIB_IMPLEMENTATION}
        -DZLIB_ROOT:PATH=${ZLIB_ROOT}
         # OpenJPEG
        -DUSE_OPENJPEG:BOOL=${USE_OPENJPEG}
        -DOpenJPEG_DIR:PATH=${OpenJPEG_DIR}
        # yaml-cpp
        -DBATCH_VERSION:BOOL=${BATCH_VERSION}
        -DYAML-CPP_DIR:PATH=${YAML-CPP_DIR}
        # Build libdcm2niixfs.a
        -DBUILD_DCM2NIIXFSLIB:BOOL=${BUILD_DCM2NIIXFSLIB}
        # Build dcm2niix library
        -DBUILD_DCM2NIIX_LIB:BOOL=${BUILD_DCM2NIIX_LIB}
)

if(SKBUILD)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/ DESTINATION dcm2niix USE_SOURCE_PERMISSIONS)
endif()
install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/ DESTINATION bin USE_SOURCE_PERMISSIONS)

option(BUILD_DOCS "Build documentation (manpages)" OFF)
if(BUILD_DOCS)
    add_subdirectory(docs)
endif()
