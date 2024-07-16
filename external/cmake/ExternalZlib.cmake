# ZLIB

SET(ZLIB_MIN_VERSION "1.2.13")

# Try to use the system library if exists
if (USE_SYSTEM_ZLIB)

    find_package(ZLIB ${ZLIB_MIN_VERSION} QUIET)

    if (ZLIB_FOUND)
        set(BUILDING_ZLIB_FROM_SOURCE FALSE CACHE INTERNAL "Using system zlib")
        set(ZLIB_INCLUDE_DIRS ${ZLIB_INCLUDE_DIRS} CACHE INTERNAL "")
        set(ZLIB_LIBRARIES ${ZLIB_LIBRARIES} CACHE INTERNAL "")
        message(STATUS "Using system zlib")
    endif()

endif()


# If library does not exist in the system then compile from source
if(NOT ZLIB_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/zlib/CMakeLists.txt")

        set(ZLIB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/zlib")

        set(BUILDING_ZLIB_FROM_SOURCE TRUE CACHE INTERNAL "zlib will be built from local source")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/zlib-${ZLIB_MIN_VERSION}/CMakeLists.txt")

        set(ZLIB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/zlib-${ZLIB_MIN_VERSION}")

        set(BUILDING_ZLIB_FROM_SOURCE TRUE CACHE INTERNAL "zlib will be built from local source")

    else()    

        set(BUILDING_ZLIB_FROM_SOURCE TRUE CACHE INTERNAL "zlib will be downloaded and built from source")

        set(ZLIB_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/zlib-${ZLIB_MIN_VERSION}")
        
        set(DOWNLOAD_FNAME "zlib-${ZLIB_MIN_VERSION}.tar.gz")
        set(DOWNLOAD_URL   "https://github.com/madler/zlib/releases/download/v${ZLIB_MIN_VERSION}/${DOWNLOAD_FNAME}")
        set(DOWNLOAD_PATH  "${CMAKE_BINARY_DIR}/external/download/zlib-${ZLIB_MIN_VERSION}.zip")
        
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


if(BUILDING_ZLIB_FROM_SOURCE)

    include(ExternalProject)
    ExternalProject_Add(build_zlib

        SOURCE_DIR "${ZLIB_SOURCE_DIR}"

        PREFIX ${CMAKE_BINARY_DIR}/external/zlib

        CMAKE_ARGS  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}

    )

    # Define variables based on conditions
    if(BUILD_SHARED_LIBS)
        if(UNIX)
            if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
                set(ZLIB_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.so")
            elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
                set(ZLIB_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.dylib")
            else()
                message(FATAL_ERROR "Unsupported processor operating system")
            endif()
        else()
            set(ZLIB_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlib.lib")
        endif()
    else()
        if(UNIX)
            if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
                set(ZLIB_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.a")
            elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
                set(ZLIB_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libz.a")
            else()
                message(FATAL_ERROR "Unsupported processor operating system")
            endif()
        else()
            set(ZLIB_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/zlibstatic.lib")
        endif()
    endif()

    ExternalProject_Add_Step(build_zlib POST_BUILD
        COMMENT "Moving Zlib headers and libraries"
        DEPENDEES install
        COMMAND ${CMAKE_COMMAND} -D nibrary=${nibrary} -D KEEP_SHARED=${BUILD_SHARED_LIBS} -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -D ZLIB_MIN_VERSION=${ZLIB_MIN_VERSION} -P "${CMAKE_CURRENT_LIST_DIR}/ExternalZlib_aux.cmake"
        ALWAYS 0
    )

    set(ZLIB_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include/${nibrary}/zlib CACHE INTERNAL "")
    set(ZLIB_LIBRARIES "${ZLIB_LIBRARY_TO_USE}" CACHE INTERNAL "ZLIB libraries to link against")


endif()


