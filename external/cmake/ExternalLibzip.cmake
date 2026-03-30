# libzip

SET(LIBZIP_MIN_VERSION "1.11.4") 

include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# Try to use the system library if exists
if (USE_SYSTEM_LIBZIP)

    find_package(libzip ${LIBZIP_MIN_VERSION} QUIET)

    if (libzip_FOUND)
        set(BUILDING_LIBZIP_FROM_SOURCE FALSE CACHE INTERNAL "Using system libzip")
        set(LIBZIP_INCLUDE_DIR ${_libzip_includes} CACHE INTERNAL "")
        message(STATUS "Using system libzip from ${LIBZIP_INCLUDE_DIR}")
    endif()

endif()

# If library does not exist in the system then compile from source
if(NOT libzip_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/libzip/CMakeLists.txt")

        set(LIBZIP_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/libzip")
        set(BUILDING_LIBZIP_FROM_SOURCE TRUE CACHE INTERNAL "libzip will be built from local source")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/libzip-${LIBZIP_MIN_VERSION}/CMakeLists.txt")

        set(LIBZIP_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/libzip-${LIBZIP_MIN_VERSION}")
        set(BUILDING_LIBZIP_FROM_SOURCE TRUE CACHE INTERNAL "libzip will be built from local source")

    else()    

        set(BUILDING_LIBZIP_FROM_SOURCE TRUE CACHE INTERNAL "libzip will be downloaded and built from source")
        set(LIBZIP_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/libzip-${LIBZIP_MIN_VERSION}")

        set(DOWNLOAD_FNAME "libzip-${LIBZIP_MIN_VERSION}.tar.gz")
        set(DOWNLOAD_URL   "https://github.com/nih-at/libzip/releases/download/v${LIBZIP_MIN_VERSION}/${DOWNLOAD_FNAME}")
        set(DOWNLOAD_PATH  "${CMAKE_BINARY_DIR}/external/download/${DOWNLOAD_FNAME}")
        
        if (NOT EXISTS ${DOWNLOAD_PATH})
            file(DOWNLOAD ${DOWNLOAD_URL} ${DOWNLOAD_PATH}
                SHOW_PROGRESS
                STATUS download_status
                LOG download_log)
            list(GET download_status 0 status_code)
            if(NOT status_code EQUAL 0)
                message(FATAL_ERROR "Error downloading ${DOWNLOAD_URL}: ${download_log}")
            endif()

            execute_process(
                COMMAND ${CMAKE_COMMAND} -E tar xzf ${DOWNLOAD_PATH}
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external"
                RESULT_VARIABLE extract_result
            )
            if(NOT extract_result EQUAL 0)
                message(FATAL_ERROR "Error extracting ${DOWNLOAD_FNAME}: ${extract_result}")
            endif()
        endif()
        
    endif()

endif()

if(BUILDING_LIBZIP_FROM_SOURCE)

    include(ExternalProject)

    if(NOT ZLIB_FOUND)
        message(STATUS "Building libzip with nibrary's zlib")

        ExternalProject_Add(build_libzip

        SOURCE_DIR "${LIBZIP_SOURCE_DIR}"

        PREFIX ${CMAKE_BINARY_DIR}/external/libzip

        CMAKE_ARGS  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                    -DCMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
                    -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
                    -DZLIB_LIBRARY:PATH=${ZLIB_LIBRARIES}
                    -DBUILD_EXAMPLES=OFF
                    -DBUILD_DOC=OFF
                    -DBUILD_TOOLS=OFF
                    -DBUILD_REGRESS=OFF
                    -DENABLE_COMMONCRYPTO=OFF
                    -DENABLE_GNUTLS=OFF
                    -DENABLE_OPENSSL=OFF
                    -DENABLE_WINDOWS_CRYPTO=OFF
                    -DENABLE_BZIP2=OFF
                    -DENABLE_LZMA=OFF
                    -DENABLE_ZSTD=OFF
                    -DENABLE_FDOPEN=OFF
                    -DENABLE_COVERAGE=OFF
                    -DENABLE_MBEDTLS=OFF
        )


    else()
        message(STATUS "Building libzip with system zlib")

        ExternalProject_Add(build_libzip

        SOURCE_DIR "${LIBZIP_SOURCE_DIR}"

        PREFIX ${CMAKE_BINARY_DIR}/external/libzip

        CMAKE_ARGS  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                    -DCMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
                    -DBUILD_EXAMPLES=OFF
                    -DBUILD_DOC=OFF
                    -DBUILD_TOOLS=OFF
                    -DBUILD_REGRESS=OFF
                    -DENABLE_COMMONCRYPTO=OFF
                    -DENABLE_GNUTLS=OFF
                    -DENABLE_OPENSSL=OFF
                    -DENABLE_WINDOWS_CRYPTO=OFF
                    -DENABLE_BZIP2=OFF
                    -DENABLE_LZMA=OFF
                    -DENABLE_ZSTD=OFF
                    -DENABLE_FDOPEN=OFF
                    -DENABLE_COVERAGE=OFF
                    -DENABLE_MBEDTLS=OFF
        )
    endif()




    ExternalProject_Add_Step(build_libzip move_libzip_headers
        COMMENT "Moving libzip headers"
        DEPENDEES install
        COMMAND ${CMAKE_COMMAND} 
            -D nibrary=${nibrary} 
            -D NIBRARY_CMAKE_INSTALL_PREFIX=${NIBRARY_CMAKE_INSTALL_PREFIX} 
            -D CMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
            -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} 
            -D LIBZIP_MIN_VERSION=${LIBZIP_MIN_VERSION}
            -P "${CMAKE_CURRENT_LIST_DIR}/ExternalLibzip_aux.cmake"
        ALWAYS 0 # This ensures the step is only run when the DEPENDEES are updated, not every build
    )

    set(LIBZIP_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/${nibrary}/libzip CACHE INTERNAL "Libzip include directory")
    
    if(BUILD_SHARED_LIBS)

        if(UNIX)
            if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
                set(LIBZIP_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.so")
            elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
                set(LIBZIP_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.dylib")
            endif()
        else()
            set(LIBZIP_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zip.lib")
        endif()


    else()

        if(UNIX)
            set(LIBZIP_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libzip.a")
        else()
            set(LIBZIP_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zip.lib")
        endif()
    
    endif()

    set(LIBZIP_LIBRARY "${LIBZIP_LIBRARY_TO_USE}" CACHE INTERNAL "Libzip libraries to link against")

   
endif()