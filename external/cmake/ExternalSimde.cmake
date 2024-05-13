# SIMDe

SET(SIMDE_MIN_VERSION "0.7.6")

# Option to use the system library if it exists
if (USE_SYSTEM_SIMDE)

    # Check for SIMDe header file
    find_path(SIMDE_INCLUDE_DIR NAMES simde/simde-features.h)

    if(SIMDE_INCLUDE_DIR)
        set(SIMDE_FOUND TRUE CACHE INTERNAL "System has SIMDe installed")
        message(STATUS "Found SIMDe: ${SIMDE_INCLUDE_DIR}")
    else()
        set(SIMDE_FOUND FALSE CACHE INTERNAL "System does not have SIMDe installed")
        message(STATUS "SIMDe not found in system include paths")
    endif()

    if (SIMDE_FOUND)
        set(BUILDING_SIMDE_FROM_SOURCE FALSE CACHE INTERNAL "Using system SIMDe")
        set(SIMDE_INCLUDE_DIR ${SIMDE_INCLUDE_DIR} CACHE INTERNAL "")
        message(STATUS "Using system SIMDe from ${SIMDE_INCLUDE_DIR}")
    endif()

endif()


if(NOT SIMDE_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/simde/simde/simde-features.h")

        set(SIMDE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/simde")

        set(BUILDING_SIMDE_FROM_SOURCE TRUE CACHE INTERNAL "Local SIMDe source will be used")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/simde-v${SIMDE_MIN_VERSION}/simde/simde-features.h")

        set(SIMDE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/simde-v${SIMDE_MIN_VERSION}")

        set(BUILDING_SIMDE_FROM_SOURCE TRUE CACHE INTERNAL "Local SIMDe source will be used")

    else()    

        set(BUILDING_SIMDE_FROM_SOURCE TRUE CACHE INTERNAL "SIMDe will be downloaded")

        set(SIMDE_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/simde-${SIMDE_MIN_VERSION}")
        
        set(DOWNLOAD_FNAME "simde-v${SIMDE_MIN_VERSION}.zip")
        set(DOWNLOAD_URL "https://github.com/simd-everywhere/simde/archive/refs/tags/v${SIMDE_MIN_VERSION}.zip")
        set(DOWNLOAD_PATH "${CMAKE_BINARY_DIR}/external/download/${DOWNLOAD_FNAME}")
        
        if (NOT EXISTS ${DOWNLOAD_PATH})
            file(DOWNLOAD ${DOWNLOAD_URL} ${DOWNLOAD_PATH}
                SHOW_PROGRESS
                STATUS download_status
                LOG download_log)
            list(GET download_status 0 status_code)
            if(NOT status_code EQUAL 0)
                message(FATAL_ERROR "Error downloading SIMDe from ${DOWNLOAD_URL}: ${download_log}")
            endif()

            execute_process(
                COMMAND ${CMAKE_COMMAND} -E tar xzf ${DOWNLOAD_PATH}
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external"
                RESULT_VARIABLE extract_result
            )
            if(NOT extract_result EQUAL 0)
                message(FATAL_ERROR "Error extracting SIMDe from ${DOWNLOAD_FNAME}: ${extract_result}")
            endif()
        endif()

    endif()

    set(SIMDE_INCLUDE_DIR ${SIMDE_SOURCE_DIR} CACHE INTERNAL "SIMDe include directory")

endif()

