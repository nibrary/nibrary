# TRX

SET(TRX_MIN_VERSION "0.1.0" CACHE STRING "Minimum trx version")

include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

if (USE_SYSTEM_TRX)
    find_package(trx ${TRX_MIN_VERSION} QUIET)
    if (trx_FOUND)
        set(BUILDING_TRX_CPP_FROM_SOURCE FALSE CACHE INTERNAL "Using system trx")
        set(TRX_INCLUDE_DIR ${TRX_INCLUDE_DIRS} CACHE INTERNAL "")
        message(STATUS "Using system trx")
    else()
        set(BUILDING_TRX_CPP_FROM_SOURCE TRUE CACHE INTERNAL "trx-cpp will be built from local source")
    endif()
endif()

if (NOT trx_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/trx-cpp/CMakeLists.txt")

        message(STATUS "trx-cpp will be built from local source")
        set(TRX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/trx-cpp")

        # Patching logic
        set(TRX_PATCH_DIR "${CMAKE_SOURCE_DIR}/external/trx_patch")
        if(EXISTS "${TRX_PATCH_DIR}/json11.cpp")
            set(JSON11_TARGET "${TRX_SOURCE_DIR}/third_party/json11/json11.cpp")
            set(JSON11_BACKUP "${TRX_SOURCE_DIR}/third_party/json11/json11_cpp_bak")
            conditional_move("${JSON11_TARGET}" "${JSON11_BACKUP}")
            conditional_copy_file("${TRX_PATCH_DIR}/json11.cpp" "${JSON11_TARGET}")
        endif()

        set(BUILDING_TRX_CPP_FROM_SOURCE TRUE CACHE INTERNAL "trx-cpp will be built from local source")

    else()
        
        ## Placeholder for downloading trx-cpp from source
        ## trx-cpp will be downloaded and built from source if not found
        # message(STATUS "trx-cpp will be downloaded and built from source")
        # set(TRX_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/trx-cpp-${TRX_MIN_VERSION}")

        # set(DOWNLOAD_FNAME "v${TRX_MIN_VERSION}.zip")
        # set(DOWNLOAD_URL "https://github.com/tee-ar-ex/trx-cpp/archive/refs/tags/${DOWNLOAD_FNAME}")
        # set(DOWNLOAD_PATH "${CMAKE_BINARY_DIR}/external/download/trx-cpp-${DOWNLOAD_FNAME}")

        # if (NOT EXISTS ${DOWNLOAD_PATH})
        #     file(DOWNLOAD ${DOWNLOAD_URL} ${DOWNLOAD_PATH}
        #         SHOW_PROGRESS
        #         STATUS download_status
        #         LOG download_log)
        #     list(GET download_status 0 status_code)
        #     if(NOT status_code EQUAL 0)
        #         message(FATAL_ERROR "Error downloading ${DOWNLOAD_URL}: ${download_log}")
        #     endif()

        #     execute_process(
        #         COMMAND ${CMAKE_COMMAND} -E tar xzf ${DOWNLOAD_PATH}
        #         WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external"
        #         RESULT_VARIABLE extract_result
        #     )
        #     if(NOT extract_result EQUAL 0)
        #         message(FATAL_ERROR "Error extracting ${DOWNLOAD_FNAME}: ${extract_result}")
        #     endif()
        # endif()

        # # Patching logic for downloaded source
        # set(TRX_PATCH_DIR "${CMAKE_SOURCE_DIR}/external/trx_patch")
        # if(EXISTS "${TRX_PATCH_DIR}/json11.cpp")
        #     set(JSON11_TARGET "${TRX_SOURCE_DIR}/third_party/json11/json11.cpp")
        #     set(JSON11_BACKUP "${TRX_SOURCE_DIR}/third_party/json11/json11_cpp_bak")
        #     conditional_move("${JSON11_TARGET}" "${JSON11_BACKUP}")
        #     conditional_copy_file("${TRX_PATCH_DIR}/json11.cpp" "${JSON11_TARGET}")
        # endif()

        # set(BUILDING_TRX_CPP_FROM_SOURCE TRUE CACHE INTERNAL "trx-cpp will be downloaded")

    endif()

endif()

if(BUILDING_TRX_CPP_FROM_SOURCE)

    include(ExternalProject)
    
    set(TRX_CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
        # Prefer the nibrary install layout so trx's find_package() prefers our internal libs
        -DCMAKE_PREFIX_PATH=${NIBRARY_CMAKE_INSTALL_PREFIX}
        # Help trx find our internal Eigen (trx checks Eigen3_ROOT first)
        -DEigen3_ROOT=${NIBRARY_CMAKE_INSTALL_PREFIX}
        -DEIGEN3_INCLUDE_DIR=${NIBRARY_CMAKE_INSTALL_PREFIX}/include
        # Provide libzip include path (trx will use TRX_LIBZIP_INCLUDE_DIR as a fallback)
        -DTRX_LIBZIP_INCLUDE_DIR=${NIBRARY_CMAKE_INSTALL_PREFIX}/include/${nibrary}/libzip
        -DTRX_BUILD_TESTS=OFF
        -DTRX_BUILD_EXAMPLES=OFF
        -DTRX_ENABLE_INSTALL=ON
    )

    ExternalProject_Add(build_trx_cpp
        SOURCE_DIR "${TRX_SOURCE_DIR}"
        PREFIX ${CMAKE_BINARY_DIR}/external/trx-cpp
        CMAKE_ARGS ${TRX_CMAKE_ARGS}
    )

    ExternalProject_Add_Step(build_trx_cpp POST_BUILD
        COMMENT "Moving trx headers and libraries"
        DEPENDEES install
        COMMAND ${CMAKE_COMMAND} 
            -D nibrary=${nibrary} 
            -D NIBRARY_CMAKE_INSTALL_PREFIX=${NIBRARY_CMAKE_INSTALL_PREFIX}
            -D CMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
            -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            -D TRX_MIN_VERSION=${TRX_MIN_VERSION}
            -P "${CMAKE_CURRENT_LIST_DIR}/ExternalTrx_aux.cmake"
        ALWAYS 0
    )

    # configure_file("${CMAKE_SOURCE_DIR}/external/trx_patch/trx.h" "${CMAKE_INSTALL_PREFIX}/include/${nibrary}/trx-cpp/trx/trx.h" COPYONLY)

    # Define variables based on conditions
    if(BUILD_SHARED_LIBS)
        if(UNIX)
            if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
                set(TRX_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libtrx.so")
            elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
                set(TRX_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libtrx.dylib")
            endif()
        else()
            set(TRX_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/trx.lib")
        endif()
    else()
        if(UNIX)
            set(TRX_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libtrx.a")
        else()
            set(TRX_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/trx.lib")
        endif()
    endif()

    set(TRX_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/${nibrary}/trx-cpp CACHE INTERNAL "")
    set(TRX_LIBRARY "${TRX_LIBRARY_TO_USE}" CACHE INTERNAL "TRX library to link against")

endif()
