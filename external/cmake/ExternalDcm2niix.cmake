# DCM2NIIX

if (BUILD_DCM2NIIX)

    SET(DCM2NIIX_VERSION "1.0.20241211" CACHE STRING "Minimum dcm2niix version") 

    include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

    if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
        set(DCM2NIIX_LIBRARY        ${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/dcm2niix++.lib  CACHE INTERNAL "")
    else()
        set(DCM2NIIX_LIBRARY        ${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/libdcm2niix++.a  CACHE INTERNAL "")
    endif()

    set(DCM2NIIX_LIBRARY_DIR    ${CMAKE_INSTALL_PREFIX}/lib/${nibrary}                  CACHE INTERNAL "")
    set(DCM2NIIX_INCLUDE_DIR    ${CMAKE_INSTALL_PREFIX}/include/${nibrary}/dcm2niix     CACHE INTERNAL "")

    if (EXISTS "${CMAKE_SOURCE_DIR}/external/dcm2niix/CMakeLists.txt")
        
        message(STATUS "dcm2niix will be built from local source")

        set(DCM2NIIX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/dcm2niix")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nifti1_io_core.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nifti1_io_core_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nifti1_io_core.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nifti1_io_core.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nii_dicom_batch.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch_h_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nii_dicom_batch.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.h")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_fswrapper.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_fswrapper_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix_fswrapper.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_fswrapper.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_h_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix++.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.h")
        
        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix++.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/CMakeLists_txt_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/CMakeLists.txt")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/SuperBuild/SuperBuild.cmake" "${CMAKE_SOURCE_DIR}/external/dcm2niix/SuperBuild/SuperBuild_cmake_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/SuperBuild.cmake" "${CMAKE_SOURCE_DIR}/external/dcm2niix/SuperBuild/SuperBuild.cmake")

        set(BUILDING_DCM2NIIX_FROM_SOURCE TRUE CACHE INTERNAL "dcm2niix will be built from local source")

    elseif (EXISTS "${CMAKE_SOURCE_DIR}/external/dcm2niix_v${DCM2NIIX_VERSION}/CMakeLists.txt")

        message(STATUS "dcm2niix will be built from local source")

        set(DCM2NIIX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/dcm2niix_v${DCM2NIIX_VERSION}")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nifti1_io_core.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nifti1_io_core_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nifti1_io_core.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nifti1_io_core.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nii_dicom_batch.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch_h_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nii_dicom_batch.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/nii_dicom_batch.h")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_fswrapper.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_fswrapper_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix_fswrapper.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_fswrapper.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_h_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix++.h" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.h")
        
        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix_cpp_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix++.cpp" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/dcm2niix++.cpp")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/console/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/CMakeLists_txt_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/CMakeLists.txt" "${CMAKE_SOURCE_DIR}/external/dcm2niix/console/CMakeLists.txt")

        conditional_move("${CMAKE_SOURCE_DIR}/external/dcm2niix/SuperBuild/SuperBuild.cmake" "${CMAKE_SOURCE_DIR}/external/dcm2niix/SuperBuild/SuperBuild_cmake_nibr_bak")
        conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/SuperBuild.cmake" "${CMAKE_SOURCE_DIR}/external/dcm2niix/SuperBuild/SuperBuild.cmake")

        set(BUILDING_DCM2NIIX_FROM_SOURCE TRUE CACHE INTERNAL "dcm2niix will be built from local source")

    else()    

        set(BUILDING_DCM2NIIX_FROM_SOURCE TRUE CACHE INTERNAL "dcm2niix will be downloaded and built from source")

        message(STATUS "dcm2niix will be downloaded and built from source")

        set(DCM2NIIX_SOURCE_DIR "${CMAKE_BINARY_DIR}/external/dcm2niix-${DCM2NIIX_VERSION}")

        set(DOWNLOAD_FNAME "v${DCM2NIIX_VERSION}.zip")
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

            conditional_move("${DCM2NIIX_SOURCE_DIR}/console/nifti1_io_core.cpp" "${DCM2NIIX_SOURCE_DIR}/console/nifti1_io_core_cpp_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nifti1_io_core.cpp" "${DCM2NIIX_SOURCE_DIR}/console/nifti1_io_core.cpp")

            conditional_move("${DCM2NIIX_SOURCE_DIR}/console/nii_dicom_batch.cpp" "${DCM2NIIX_SOURCE_DIR}/console/nii_dicom_batch_cpp_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nii_dicom_batch.cpp" "${DCM2NIIX_SOURCE_DIR}/console/nii_dicom_batch.cpp")

            conditional_move("${DCM2NIIX_SOURCE_DIR}/console/nii_dicom_batch.h" "${DCM2NIIX_SOURCE_DIR}/console/nii_dicom_batch_h_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/nii_dicom_batch.h" "${DCM2NIIX_SOURCE_DIR}/console/nii_dicom_batch.h")

            conditional_move("${DCM2NIIX_SOURCE_DIR}/console/dcm2niix_fswrapper.cpp" "${DCM2NIIX_SOURCE_DIR}/console/dcm2niix_fswrapper_cpp_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix_fswrapper.cpp" "${DCM2NIIX_SOURCE_DIR}/console/dcm2niix_fswrapper.cpp")

            conditional_move("${DCM2NIIX_SOURCE_DIR}/console/dcm2niix++.h" "${DCM2NIIX_SOURCE_DIR}/console/dcm2niix_h_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix++.h" "${DCM2NIIX_SOURCE_DIR}/console/dcm2niix++.h")
            
            conditional_move("${DCM2NIIX_SOURCE_DIR}/console/dcm2niix++.cpp" "${DCM2NIIX_SOURCE_DIR}/console/dcm2niix_cpp_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/dcm2niix++.cpp" "${DCM2NIIX_SOURCE_DIR}/console/dcm2niix++.cpp")
            
            conditional_move("${DCM2NIIX_SOURCE_DIR}/console/CMakeLists.txt" "${DCM2NIIX_SOURCE_DIR}/console/CMakeLists_txt_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/CMakeLists.txt" "${DCM2NIIX_SOURCE_DIR}/console/CMakeLists.txt")

            conditional_move("${DCM2NIIX_SOURCE_DIR}/SuperBuild/SuperBuild.cmake" "${DCM2NIIX_SOURCE_DIR}/SuperBuild/SuperBuild_cmake_nibr_bak")
            conditional_copy_file("${CMAKE_SOURCE_DIR}/external/dcm2niix_patch/SuperBuild.cmake" "${DCM2NIIX_SOURCE_DIR}/SuperBuild/SuperBuild.cmake")

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
                -DUSE_STATIC_RUNTIME=$<$<BOOL:${BUILD_SHARED_LIBS}>:OFF;ON>
                $<$<BOOL:${BUILD_SHARED_LIBS}>:
                    -DCMAKE_C_FLAGS=-fPIC /MD
                    -DCMAKE_CXX_FLAGS=-fPIC /MD
                >
                -DUSE_JNIFTI=OFF
                -DBUILD_DCM2NIIX_LIB=ON
        )

        ExternalProject_Add_Step(build_dcm2niix POST_BUILD
            COMMENT "Moving dcm2niix headers and libraries"
            DEPENDEES install
            COMMAND ${CMAKE_COMMAND} 
                -D nibrary=${nibrary} 
                -D NIBRARY_CMAKE_INSTALL_PREFIX=${NIBRARY_CMAKE_INSTALL_PREFIX} 
                -D NIBRARY_CMAKE_SOURCE_DIR=${NIBRARY_CMAKE_SOURCE_DIR} 
                -D CMAKE_INSTALL_PREFIX=${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}
                -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                -D CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME} 
                -D DCM2NIIX_VERSION=${DCM2NIIX_VERSION} 
                -P "${CMAKE_CURRENT_LIST_DIR}/ExternalDcm2niix_aux.cmake"
            ALWAYS 0
        )

    endif()

endif()

