# EIGEN

SET(EIGEN_MIN_VERSION "3.4.0")

# Option to use the system library if it exists
if (USE_SYSTEM_EIGEN)

    find_package(Eigen3 ${EIGEN_MIN_VERSION} QUIET)

    if (Eigen3_FOUND)
        set(BUILDING_EIGEN_FROM_SOURCE FALSE CACHE INTERNAL "Using system Eigen")
        set(EIGEN_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR} CACHE INTERNAL "")
        message(STATUS "Using system Eigen")
    endif()

endif()

# If Eigen does not exist in the system then use the local version
if(NOT Eigen3_FOUND)

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/eigen/CMakeLists.txt")

        set(Eigen_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/eigen")

        set(BUILDING_EIGEN_FROM_SOURCE TRUE CACHE INTERNAL "Local Eigen source will be used")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/eigen-${EIGEN_MIN_VERSION}/CMakeLists.txt")

        set(Eigen_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/eigen-${EIGEN_MIN_VERSION}")

        set(BUILDING_EIGEN_FROM_SOURCE TRUE CACHE INTERNAL "Local Eigen source will be used")

    else()    

        set(BUILDING_EIGEN_FROM_SOURCE TRUE CACHE INTERNAL "Eigen will be downloaded")

        set(Eigen_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/eigen-${EIGEN_MIN_VERSION}")
        
        set(DOWNLOAD_FNAME "eigen-${EIGEN_MIN_VERSION}.zip")
        set(DOWNLOAD_URL "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_MIN_VERSION}/${DOWNLOAD_FNAME}")
        set(DOWNLOAD_PATH "${CMAKE_BINARY_DIR}/external/download/eigen-${EIGEN_MIN_VERSION}.zip")
        
        if (NOT EXISTS ${DOWNLOAD_PATH})
            file(DOWNLOAD ${DOWNLOAD_URL} ${DOWNLOAD_PATH}
                SHOW_PROGRESS
                STATUS download_status
                LOG download_log)
            list(GET download_status 0 status_code)
            if(NOT status_code EQUAL 0)
                message(FATAL_ERROR "Error downloading Eigen from ${DOWNLOAD_URL}: ${download_log}")
            endif()

            execute_process(
                COMMAND ${CMAKE_COMMAND} -E tar xzf ${DOWNLOAD_PATH}
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external"
                RESULT_VARIABLE extract_result
            )
            if(NOT extract_result EQUAL 0)
                message(FATAL_ERROR "Error extracting Eigen from ${DOWNLOAD_FNAME}: ${extract_result}")
            endif()
        endif()

    endif()

    set(EIGEN_INCLUDE_DIR ${Eigen_SOURCE_DIR} CACHE INTERNAL "Eigen include directory")

endif()

