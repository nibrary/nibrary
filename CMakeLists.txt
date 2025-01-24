cmake_minimum_required(VERSION 3.15.0)

PROJECT(nibrary VERSION 0.3.0 LANGUAGES CXX C)

set(nibrary "nibrary_v${nibrary_VERSION}")
message(STATUS ">> Nibrary::Building ${nibrary}")

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")

option(BUILD_SHARED_LIBS        "Build shared libraries"                ON)
option(USE_SYSTEM_EIGEN         "Use system eigen if exists"            OFF)
option(USE_SYSTEM_ZLIB          "Use system zlib if exists"             OFF)
# option(USE_SYSTEM_GEOGRAM       "Use system geogram if exists"          OFF) - Remove support for this until the dependency of mesh_fill_holes.cpp is resolved
option(USE_SYSTEM_LIBIGL        "Use system libigl if exists"           OFF)
option(USE_SYSTEM_SIMDE         "Use system simde if exists"            OFF)
option(USE_SYSTEM_PROXSUITE     "Use system proxsuite if exists"        OFF)
option(ENABLE_CUDA              "Enable CUDA functions"                 OFF)
option(BUILD_DCM2NIIX           "Build dcm2niix (Unix only)"            ON)

# Using C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

# Cuda
include(cuda)

# Include external libraries/code
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/external)
include_directories(${CMAKE_SOURCE_DIR}/external/expat)
include_directories(${CMAKE_SOURCE_DIR}/external/marchingCubes)
include_directories(${CMAKE_SOURCE_DIR}/external/triangleVoxelIntersection)

# List source files
file(GLOB_RECURSE SRCS CONFIGURE_DEPENDS RELATIVE ${CMAKE_SOURCE_DIR}
    "src/*.c*"
    "external/niftilib/*.c*"
    "external/marchingCubes/*.c*"
    "external/gifticlib/*.c*"
    "external/expat/*.c*"
    "external/triangleVoxelIntersection/*.c*")

# Remove CUDA sources to be compiled by C++ compiler
file(GLOB_RECURSE CUDA_SRCS CONFIGURE_DEPENDS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cu*")
foreach(FILE ${CUDA_SRCS})
    list(REMOVE_ITEM SRCS ${FILE})
endforeach()

# Compiler definitions
add_compile_definitions(EIGEN_MPL2_ONLY)
add_compile_definitions(HAVE_ZLIB)
add_compile_definitions(PROXSUITE_VECTORIZE)

if(NOT BUILD_SHARED_LIBS)
    add_compile_definitions(XML_STATIC)
    add_compile_definitions(PROXSUITE_STATIC)
endif()

# Set the library
set(libNibrary "Nibrary")

# Compile the sources
add_library(OBJS OBJECT ${SRCS})

# Include external dependencies
add_subdirectory(external)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(${GEOGRAM_INCLUDE_DIR})
include_directories(${LIBIGL_INCLUDE_DIR})
include_directories(${SIMDE_INCLUDE_DIR})
include_directories(${PROXSUITE_INCLUDE_DIR})
include_directories(${CUDA_INCLUDE_DIR})

if (BUILD_DCM2NIIX AND UNIX AND NOT BUILD_SHARED_LIBS)  # For now let's experiment only in Unix
    message(STATUS ">> Nibrary::Building dcm2niix")
    add_compile_definitions(HAS_DCM2NIIX)
    include_directories(${DCM2NIIX_INCLUDE_DIR})
    add_dependencies(OBJS build_dcm2niix)
else()
    message(STATUS ">> Nibrary::dcm2niix will not be built")
endif()

if(BUILDING_ZLIB_FROM_SOURCE)
    add_dependencies(OBJS build_zlib)
endif()

if(BUILDING_GEOGRAM_FROM_SOURCE)
    add_dependencies(OBJS build_geogram)
endif()


IF(BUILD_SHARED_LIBS)
    message(STATUS ">> Nibrary::Building shared library")
    set_property(TARGET OBJS PROPERTY POSITION_INDEPENDENT_CODE ON)
    add_library(${libNibrary} SHARED $<TARGET_OBJECTS:OBJS>)
ELSE()
    message(STATUS ">> Nibrary::Building static library")
    add_library(${libNibrary} STATIC $<TARGET_OBJECTS:OBJS>)
ENDIF()

find_package(OpenMP REQUIRED)

if(OpenMP_FOUND)
    target_compile_options(${libNibrary} PUBLIC ${OpenMP_CXX_FLAGS})
    target_link_libraries(${libNibrary} PUBLIC OpenMP::OpenMP_CXX)
    target_include_directories(${libNibrary} PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})
endif()

target_link_libraries(${libNibrary} PUBLIC ${ZLIB_LIBRARIES} ${GEOGRAM_LIBRARY})

IF(ENABLE_CUDA)
    add_library(CUDA_OBJS OBJECT ${CUDA_SRCS})
    set(libNibrary "${libNibrary}CUDA_${CMAKE_CUDA_COMPILER_VERSION}")
    target_link_libraries(${libNibrary} CUDA_OBJS)
ENDIF(ENABLE_CUDA)



if ((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang") )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wno-long-long -Wno-stringop-truncation -Wno-dev -Wno-unknown-pragmas -Wno-unused-result -Wall -Wextra -pedantic -pedantic-errors" CACHE STRING "C++ Compiler Flags" FORCE)
endif()

if ((CMAKE_CXX_COMPILER_ID MATCHES "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9"))
    target_link_libraries(${libNibrary} PUBLIC stdc++fs)
endif()

if (UNIX)
    # Debug messages
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    endif()
endif()

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    message(STATUS ">> Nibrary::Building for Unix")

elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")

    message(STATUS ">> Nibrary::Building for MacOs")

elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")

    message(STATUS ">> Nibrary::Building for Windows")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /permissive- /bigobj")

    # For definitions
    add_compile_definitions(_USE_MATH_DEFINES)

    # Disable some compiler warnings
    set(DISABLED_WARNINGS "/wd4244 /wd4267 /wd4996 /wd4305 /wd4101 /wd4068 /wd4661 /wd4477 /wd4804 /wd4700")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DISABLED_WARNINGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DISABLED_WARNINGS}")

    set_target_properties(${libNibrary} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

else()
    message(FATAL_ERROR ">> Nibrary::This operating system is not supported")
endif()

install(TARGETS ${libNibrary} LIBRARY DESTINATION lib/${nibrary} ARCHIVE DESTINATION lib/${nibrary})

install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/ DESTINATION include/${nibrary} FILES_MATCHING PATTERN "*.h" PATTERN "*.cuh")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/expat DESTINATION include/${nibrary} FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/fast-marching-method DESTINATION include/${nibrary} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/gifticlib DESTINATION include/${nibrary} FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/marchingCubes DESTINATION include/${nibrary} FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/niftilib DESTINATION include/${nibrary} FILES_MATCHING PATTERN "*.h")

if(BUILDING_EIGEN_FROM_SOURCE)
    install(DIRECTORY ${EIGEN_INCLUDE_DIR}/Eigen DESTINATION include/${nibrary})
endif()

if(BUILDING_LIBIGL_FROM_SOURCE)
    install(DIRECTORY ${LIBIGL_INCLUDE_DIR}/igl DESTINATION include/${nibrary})
endif()

if(BUILDING_SIMDE_FROM_SOURCE)
    install(DIRECTORY ${SIMDE_INCLUDE_DIR}/simde DESTINATION include/${nibrary})
endif()

if(BUILDING_PROXSUITE_FROM_SOURCE)
    install(DIRECTORY ${PROXSUITE_INCLUDE_DIR}/proxsuite DESTINATION include/${nibrary})
endif()

install(CODE "file(REMOVE_RECURSE \"${NIBRARY_EXTERNAL_CMAKE_INSTALL_PREFIX}\")")

add_custom_target(uninstall 
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/include/${nibrary}"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}"
)


