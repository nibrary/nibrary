cmake_minimum_required(VERSION 3.15.0)

PROJECT(nibrary VERSION 0.1 LANGUAGES CXX C)

set(nibrary "nibrary_v${nibrary_VERSION}")
message(STATUS ">> Nibrary::Building ${nibrary}")

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")

option(BUILD_SHARED_LIBS        "Build shared libraries"                ON)
option(USE_SYSTEM_EIGEN         "Use system eigen if exists"            OFF)
option(USE_SYSTEM_ZLIB          "Use system zlib if exists"             OFF)
option(USE_SYSTEM_GEOGRAM       "Use system geogram if exists"          OFF)
option(USE_SYSTEM_LIBIGL        "Use system libigl if exists"           OFF)
option(USE_SYSTEM_SIMDE         "Use system simde if exists"            OFF)
option(USE_SYSTEM_PROXSUITE     "Use system proxsuite if exists"        OFF)
option(ENABLE_CUDA              "Enable CUDA functions"                 OFF)


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


# Include external dependencies
add_subdirectory(external)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(${GEOGRAM_INCLUDE_DIR})
include_directories(${LIBIGL_INCLUDE_DIR})
include_directories(${SIMDE_INCLUDE_DIR})
include_directories(${PROXSUITE_INCLUDE_DIR})
include_directories(${CUDA_INCLUDE_DIR})

find_package(OpenMP REQUIRED) # Needed for geogram
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Compile the sources
add_library(OBJS OBJECT ${SRCS})

if(BUILDING_ZLIB_FROM_SOURCE)
    add_dependencies(OBJS build_zlib)
endif()

if(BUILDING_GEOGRAM_FROM_SOURCE)
    add_dependencies(OBJS build_geogram)
endif()

# Set the library
set(libNibrary "Nibrary")

IF(BUILD_SHARED_LIBS)
    message(STATUS ">> Nibrary::Building shared library")
    set_property(TARGET OBJS PROPERTY POSITION_INDEPENDENT_CODE ON)
    add_library(${libNibrary} SHARED $<TARGET_OBJECTS:OBJS>)
ELSE()
    message(STATUS ">> Nibrary::Building static library")
    add_library(${libNibrary} STATIC $<TARGET_OBJECTS:OBJS>)
ENDIF()

target_link_libraries(${libNibrary} PUBLIC ${ZLIB_LIBRARIES} ${GEOGRAM_LIBRARY})


IF(ENABLE_CUDA)
    add_library(CUDA_OBJS OBJECT ${CUDA_SRCS})
    set(libNibrary "${libNibrary}CUDA_${CMAKE_CUDA_COMPILER_VERSION}")
    target_link_libraries(${libNibrary} CUDA_OBJS)
ENDIF(ENABLE_CUDA)


if (CMAKE_SYSTEM_NAME STREQUAL "Linux")

    message(STATUS ">> Nibrary::Building for Unix")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -march=native")

    if(CMAKE_BUILD_TYPE MATCHES Debug)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wno-long-long -Wno-stringop-truncation -Wno-dev -Wno-unknown-pragmas -Wno-unused-result -Wall -Wextra -pedantic -pedantic-errors" CACHE STRING "C++ Compiler Flags" FORCE)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9")
        target_link_libraries(${libNibrary} PUBLIC stdc++fs)
    endif()

elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")

    message(STATUS ">> Nibrary::Building for MacOs")
    add_compile_definitions(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -march=native")

elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")

    message(STATUS ">> Nibrary::Building for Windows")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /permissive- /bigobj")

    # For definitions
    add_compile_definitions(BUILD_FOR_WINDOWS _USE_MATH_DEFINES _WIN32 WIN32 _WIN64 WIN64)

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

install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/doc\")")
install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/share\")")
install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/cmake\")")
install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/lib/${nibrary}/pkgconfig\")")

add_custom_target(uninstall 
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/include/${nibrary}"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/lib/${nibrary}"
)


