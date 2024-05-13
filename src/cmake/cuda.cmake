# CUDA

IF(ENABLE_CUDA)

    include(CheckLanguage)

    check_language(CUDA)
    
    FIND_PACKAGE(CUDA REQUIRED)
    
    INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
    
    SET(CUDA_INCLUDE_DIR ${CUDA_INCLUDE_DIRS} CACHE INTERNAL "")

    if(CMAKE_CUDA_COMPILER)

        enable_language(CUDA)

        add_definitions(-D HAVE_CUDA)

        if(NOT DEFINED CMAKE_CUDA_STANDARD)
            set(CMAKE_CUDA_STANDARD 11)
            set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        endif()

        message(STATUS "Enabled CUDA support")

    else()

        message(STATUS "Disabled CUDA support")

    endif()

ENDIF(ENABLE_CUDA)
