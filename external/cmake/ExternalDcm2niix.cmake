# DCM2NIIX

SET(DCM2NIIX_VERSION "1.0.20241211" CACHE STRING "Minimum dcm2niix version") 

include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

set(DCM2NIIX_LIBRARY        ${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libdcm2niixfs.a  CACHE INTERNAL "")
set(DCM2NIIX_LIBRARY_DIR    ${CMAKE_INSTALL_PREFIX}/lib/${nibrary}                  CACHE INTERNAL "")
set(DCM2NIIX_INCLUDE_DIR    ${CMAKE_INSTALL_PREFIX}/include/${nibrary}/dcm2niix     CACHE INTERNAL "")

if (EXISTS "${CMAKE_SOURCE_DIR}/external/dcm2niix/CMakeLists.txt")
    
    message(STATUS "dcm2niix will be built from local source")

    set(DCM2NIIX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/dcm2niix")

    set(BUILDING_DCM2NIIX_FROM_SOURCE TRUE CACHE INTERNAL "dcm2niix will be built from local source")

elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/dcm2niix_v${DCM2NIIX_VERSION}/CMakeLists.txt")

    message(STATUS "dcm2niix will be built from local source")

    set(DCM2NIIX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/dcm2niix_v${DCM2NIIX_VERSION}")

    set(BUILDING_DCM2NIIX_FROM_SOURCE TRUE CACHE INTERNAL "dcm2niix will be built from local source")

else()    

    set(BUILDING_DCM2NIIX_FROM_SOURCE TRUE CACHE INTERNAL "dcm2niix will be downloaded and built from source")

    message(STATUS "dcm2niix will be downloaded and built from source")

    set(DCM2NIIX_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/dcm2niix_v${DCM2NIIX_VERSION}")

    set(DOWNLOAD_FNAME "dcm2niix_v${DCM2NIIX_VERSION}.zip")
    set(DOWNLOAD_URL   "https://github.com/rordenlab/dcm2niix/archive/refs/tags/${DOWNLOAD_FNAME}")
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

if(BUILDING_DCM2NIIX_FROM_SOURCE)

    include(ExternalProject)
    
    ExternalProject_Add(build_dcm2niix

        SOURCE_DIR "${DCM2NIIX_SOURCE_DIR}"

        PREFIX ${CMAKE_BINARY_DIR}/external/dcm2niix

        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            -DCMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DBUILD_DCM2NIIXFSLIB=ON
    )

    ExternalProject_Add_Step(build_dcm2niix POST_BUILD
        COMMENT "Moving dcm2niix headers and libraries"
        DEPENDEES install
        COMMAND ${CMAKE_COMMAND} 
            -D nibrary=${nibrary} 
            -D NIBRARY_CMAKE_INSTALL_PREFIX=${NIBRARY_CMAKE_INSTALL_PREFIX} 
            -D NIBRARY_CMAKE_SOURCE_DIR=${NIBRARY_CMAKE_SOURCE_DIR} 
            -D CMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX} 
            -D DCM2NIIX_VERSION=${DCM2NIIX_VERSION} 
            -P "${CMAKE_CURRENT_LIST_DIR}/ExternalDcm2niix_aux.cmake"
        ALWAYS 0
    )

endif()

