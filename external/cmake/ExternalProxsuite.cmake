# PROXSUITE

SET(PROXSUITE_MIN_VERSION "0.6.3") 

# Try to use the system library if exists
if (USE_SYSTEM_PROXSUITE)

    find_package(PROXSUITE ${PROXSUITE_MIN_VERSION} QUIET)

    if (PROXSUITE_FOUND)
        set(BUILDING_PROXSUITE_FROM_SOURCE FALSE CACHE INTERNAL "Using system proxsuite")
        set(PROXSUITE_INCLUDE_DIR ${PROXSUITE_INCLUDE_DIR} CACHE INTERNAL "")
        message(STATUS "Using system proxsuite from ${PROXSUITE_INCLUDE_DIR}")
    endif()

endif()

# If library does not exist in the system then compile from source
if(NOT PROXSUITE_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/proxsuite/CMakeLists.txt")

        set(PROXSUITE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/proxsuite/include")
        file(COPY "${CMAKE_SOURCE_DIR}/external/proxsuite_config/config.hpp" DESTINATION "${PROXSUITE_SOURCE_DIR}/proxsuite")
        set(BUILDING_PROXSUITE_FROM_SOURCE TRUE CACHE INTERNAL "Proxsuite will be built from local source")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/proxsuite-${PROXSUITE_MIN_VERSION}/CMakeLists.txt")

        set(PROXSUITE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/proxsuite-${PROXSUITE_MIN_VERSION}/include")
        file(COPY "${CMAKE_SOURCE_DIR}/external/proxsuite_config/config.hpp" DESTINATION "${PROXSUITE_SOURCE_DIR}/proxsuite")
        set(BUILDING_PROXSUITE_FROM_SOURCE TRUE CACHE INTERNAL "Proxsuite will be built from local source")

    else()    

        set(BUILDING_PROXSUITE_FROM_SOURCE TRUE CACHE INTERNAL "Proxsuite will be downloaded and built from source")
        set(PROXSUITE_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/proxsuite-${PROXSUITE_MIN_VERSION}/include")

        set(DOWNLOAD_FNAME "proxsuite-${PROXSUITE_MIN_VERSION}.tar.gz")
        set(DOWNLOAD_URL   "https://github.com/Simple-Robotics/proxsuite/releases/download/v${PROXSUITE_MIN_VERSION}/${DOWNLOAD_FNAME}")
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

        file(COPY "${CMAKE_SOURCE_DIR}/external/proxsuite_config/config.hpp" DESTINATION "${PROXSUITE_SOURCE_DIR}/proxsuite")

    endif()

    set(PROXSUITE_INCLUDE_DIR ${PROXSUITE_SOURCE_DIR} CACHE INTERNAL "Proxsuite include directory")

endif()
