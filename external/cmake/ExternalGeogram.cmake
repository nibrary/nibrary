# GEOGRAM

SET(GEOGRAM_MIN_VERSION "1.9.7" CACHE STRING "Minimum geogram version") 

include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

# Try to use the system library if exists
if (USE_SYSTEM_GEOGRAM)

    find_package(Geogram ${GEOGRAM_MIN_VERSION} QUIET)

    if (GEOGRAM_FOUND)
        set(BUILDING_GEOGRAM_FROM_SOURCE FALSE CACHE INTERNAL "Using system geogram")
        set(GEOGRAM_INCLUDE_DIR ${GEOGRAM_INCLUDE_DIR} CACHE INTERNAL "")
        set(GEOGRAM_LIBRARY ${GEOGRAM_LIBRARY} CACHE INTERNAL "")
        
        # mesh_fill_holes function actually uses a geogram .cpp file, which is not ideal but
        # for now, this was found to be good solution. Until we have a better solution,
        # we will download the source and compile geogram. So system package is not used for now.
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/mesh_fill_holes.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/mesh/mesh_fill_holes.cpp")
        message(STATUS "Using system geogram")
    else()
        message(STATUS "Geogram not found")
    endif()

endif()

# If library does not exist in the system then compile from source
if(NOT GEOGRAM_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/geogram/CMakeLists.txt")

        set(GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/geogram")

        set(BUILDING_GEOGRAM_FROM_SOURCE TRUE CACHE INTERNAL "Geogram will be built from local source")
        
        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/assert.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/assert_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/assert.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/assert.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/process_unix.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/process_unix_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/process_unix.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/process_unix.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/CMakeOptions.txt" "${CMAKE_SOURCE_DIR}/external/geogram/CMakeOptions_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/CMakeOptions.txt" "${CMAKE_SOURCE_DIR}/external/geogram/CMakeOptions.txt")

        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/src/lib/third_party/glfw/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/third_party/glfw/CMakeLists_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/third_party/glfw/CMakeLists.txt")

        message(STATUS "Geogram will be built from local source")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/geogram_${GEOGRAM_MIN_VERSION}/CMakeLists.txt")

        set(GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/geogram_${GEOGRAM_MIN_VERSION}")

        set(BUILDING_GEOGRAM_FROM_SOURCE TRUE CACHE INTERNAL "Geogram will be built from local source")

        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/assert.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/assert_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/assert.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/assert.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/process_unix.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/process_unix_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/process_unix.cpp" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/geogram/basic/process_unix.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/CMakeOptions.txt" "${CMAKE_SOURCE_DIR}/external/geogram/CMakeOptions_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/CMakeOptions.txt" "${CMAKE_SOURCE_DIR}/external/geogram/CMakeOptions.txt")

        conditional_move("${CMAKE_SOURCE_DIR}/external/geogram/src/lib/third_party/glfw/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/third_party/glfw/CMakeLists_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/geogram/src/lib/third_party/glfw/CMakeLists.txt")

        message(STATUS "Geogram will be built from local source")

    else()    

        set(BUILDING_GEOGRAM_FROM_SOURCE TRUE CACHE INTERNAL "Geogram will be downloaded and built from source")

        message(STATUS "Geogram will be downloaded and built from source")

        set(GEOGRAM_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/geogram_${GEOGRAM_MIN_VERSION}")

        set(DOWNLOAD_FNAME "geogram_${GEOGRAM_MIN_VERSION}.zip")
        set(DOWNLOAD_URL   "https://github.com/BrunoLevy/geogram/releases/download/v${GEOGRAM_MIN_VERSION}/${DOWNLOAD_FNAME}")
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

        conditional_move("${GEOGRAM_SOURCE_DIR}/src/lib/geogram/basic/assert.cpp" "${GEOGRAM_SOURCE_DIR}/src/lib/geogram/basic/assert_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/assert.cpp" "${GEOGRAM_SOURCE_DIR}/src/lib/geogram/basic/assert.cpp")

        conditional_move("${GEOGRAM_SOURCE_DIR}/src/lib/geogram/basic/process_unix.cpp" "${GEOGRAM_SOURCE_DIR}/src/lib/geogram/basic/process_unix_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/process_unix.cpp" "${GEOGRAM_SOURCE_DIR}/src/lib/geogram/basic/process_unix.cpp")

        conditional_move("${GEOGRAM_SOURCE_DIR}/CMakeOptions.txt" "${GEOGRAM_SOURCE_DIR}/CMakeOptions_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/CMakeOptions.txt" "${GEOGRAM_SOURCE_DIR}/CMakeOptions.txt")

        conditional_move("${CMAKE_SOURCE_DIR}/src/lib/third_party/glfw/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/src/lib/third_party/glfw/CMakeLists_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/geogram_patch/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/src/lib/third_party/glfw/CMakeLists.txt")


    endif()

endif()


if(BUILDING_GEOGRAM_FROM_SOURCE)

    include(ExternalProject)
    ExternalProject_Add(build_geogram

        SOURCE_DIR "${GEOGRAM_SOURCE_DIR}"

        PREFIX ${CMAKE_BINARY_DIR}/external/geogram

        CMAKE_ARGS  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                    -DCMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
                    -DGEOGRAM_LIB_ONLY=ON
                    -DGEOGRAM_WITH_TETGEN=OFF
                    -DGEOGRAM_WITH_EXPLORAGRAM=OFF
                    -DGEOGRAM_WITH_GRAPHICS=OFF
                    -DGEOGRAM_WITH_LEGACY_NUMERICS=OFF
                    -DGEOGRAM_WITH_HLBFGS=ON
                    -DGEOGRAM_WITH_TRIANGLE=OFF
                    -DGEOGRAM_WITH_LUA=OFF
                    -DGEOGRAM_WITH_FPG=OFF
                    -DGEOGRAM_WITH_GARGANTUA=OFF
                    -DGEOGRAM_WITH_GEOGRAMPLUS=OFF
                    -DGEOGRAM_WITH_VORPALINE=OFF
    )

    ExternalProject_Add_Step(build_geogram move_geogram_headers
        COMMENT "Moving Geogram headers from /include/geogram1/geogram to /include/${nibrary}/geogram"
        DEPENDEES install # Assuming 'install' is the last step of building geogram, adjust as needed
        COMMAND ${CMAKE_COMMAND} 
            -D nibrary=${nibrary} 
            -D NIBRARY_CMAKE_INSTALL_PREFIX=${NIBRARY_CMAKE_INSTALL_PREFIX} 
            -D CMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
            -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} 
            -D GEOGRAM_MIN_VERSION=${GEOGRAM_MIN_VERSION}
            -P "${CMAKE_CURRENT_LIST_DIR}/ExternalGeogram_aux.cmake"
        ALWAYS 0 # This ensures the step is only run when the DEPENDEES are updated, not every build
    )

    # We will for now set this as follows, not as ${CMAKE_INSTALL_PREFIX}/include/geogram
    # This is because mesh_fill_holes function actually uses a geogram .cpp file, which is not ideal but
    # for now, this was found to be working solution. Until we have a better solution, let's keep the line below.
    set(GEOGRAM_INCLUDE_DIR ${GEOGRAM_SOURCE_DIR}/src/lib CACHE INTERNAL "Geogram include directory")
    
    if(BUILD_SHARED_LIBS)

        if(UNIX)
            if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
                set(GEOGRAM_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.so")
            elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
                set(GEOGRAM_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.dylib")
            endif()
        else()
            set(GEOGRAM_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/geogram.lib")
        endif()


    else()

        if(UNIX)
            set(GEOGRAM_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libgeogram.a")
        else()
            set(GEOGRAM_LIBRARY_TO_USE "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/geogram.lib")
        endif()
    
    endif()

    set(GEOGRAM_LIBRARY "${GEOGRAM_LIBRARY_TO_USE}" CACHE INTERNAL "Geogram libraries to link against")

   
endif()

