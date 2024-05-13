# LIBIGL

SET(LIBIGL_MIN_VERSION "2.5.0")

# Option to use the system library if it exists
if (USE_SYSTEM_LIBIGL)

    find_package(LIBIGL ${LIBIGL_MIN_VERSION} QUIET)

    if (LIBIGL_FOUND)
        set(BUILDING_LIBIGL_FROM_SOURCE FALSE CACHE INTERNAL "Using system libigl")
        set(LIBIGL_INCLUDE_DIR ${LIBIGL_INCLUDE_DIRS} CACHE INTERNAL "")
        message(STATUS "Using system libigl")
    endif()

endif()


if(NOT LIBIGL_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/libigl/CMakeLists.txt")

        set(LIBIGL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/libigl")

        set(BUILDING_LIBIGL_FROM_SOURCE TRUE CACHE INTERNAL "Local libigl source will be used")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/libigl-${LIBIGL_MIN_VERSION}/CMakeLists.txt")

        set(LIBIGL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/libigl-${LIBIGL_MIN_VERSION}")

        set(BUILDING_LIBIGL_FROM_SOURCE TRUE CACHE INTERNAL "Local libigl source will be used")

    else()    

        set(BUILDING_LIBIGL_FROM_SOURCE TRUE CACHE INTERNAL "libigl will be downloaded")

        set(LIBIGL_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/libigl-${LIBIGL_MIN_VERSION}")
        
        set(DOWNLOAD_FNAME "libigl-v${LIBIGL_MIN_VERSION}.zip")
        set(DOWNLOAD_URL "https://github.com/libigl/libigl/archive/refs/tags/v${LIBIGL_MIN_VERSION}.zip")
        set(DOWNLOAD_PATH "${CMAKE_BINARY_DIR}/external/download/${DOWNLOAD_FNAME}")
        
        if (NOT EXISTS ${DOWNLOAD_PATH})
            file(DOWNLOAD ${DOWNLOAD_URL} ${DOWNLOAD_PATH}
                SHOW_PROGRESS
                STATUS download_status
                LOG download_log)
            list(GET download_status 0 status_code)
            if(NOT status_code EQUAL 0)
                message(FATAL_ERROR "Error downloading libigl from ${DOWNLOAD_URL}: ${download_log}")
            endif()

            execute_process(
                COMMAND ${CMAKE_COMMAND} -E tar xzf ${DOWNLOAD_PATH}
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external"
                RESULT_VARIABLE extract_result
            )
            if(NOT extract_result EQUAL 0)
                message(FATAL_ERROR "Error extracting libigl from ${DOWNLOAD_FNAME}: ${extract_result}")
            endif()
        endif()

    endif()

    set(LIBIGL_INCLUDE_DIR ${LIBIGL_SOURCE_DIR}/include CACHE INTERNAL "libigl include directory")

endif()

