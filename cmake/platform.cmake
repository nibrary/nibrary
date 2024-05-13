# Set platform specific options

# Using a separate .cmake as in the below trigger the recompilation of 
# the whole code base when cmake is run. So we are moving the below content under the
# main CMakeLists.txt for now

# if(UNIX)

#     message(STATUS "Building for Unix")


#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -march=native" CACHE STRING "C++ Compiler Flags" FORCE)

#     if(CMAKE_BUILD_TYPE MATCHES Debug)
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address" CACHE STRING "C Compiler Flags" FORCE)
#         set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address" CACHE STRING "Executable Linker Flags" FORCE)
#     endif()

#     if(CMAKE_CXX_COMPILER_ID MATCHES "gcc")

#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wno-long-long -Wno-stringop-truncation -Wno-dev -Wno-unknown-pragmas -Wno-unused-result -Wall -Wextra -pedantic -pedantic-errors" CACHE STRING "C++ Compiler Flags" FORCE)

#         if(NOT BUILD_SHARED_LIBS)
#             set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static" CACHE STRING "Executable Linker Flags" FORCE)
#         endif()

#     elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w" CACHE STRING "C++ Compiler Flags" FORCE)
#         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w" CACHE STRING "C Compiler Flags" FORCE)

#     endif()

# elseif(MSVC)

#     message(STATUS "Building for Windows")

#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /permissive- /bigobj" CACHE STRING "C++ Compiler Flags" FORCE)

#     # add_definitions(-D BUILD_FOR_WINDOWS)
#     # add_definitions(-D _USE_MATH_DEFINES)
#     # add_definitions(-D _WIN32)
#     # add_definitions(-D WIN32)
#     # add_definitions(-D _WIN64)
#     # add_definitions(-D WIN64)

#     # For definitions
#     set(ADDITIONAL_DEFINITIONS "-D BUILD_FOR_WINDOWS -D _USE_MATH_DEFINES -D _WIN32 -D WIN32 -D _WIN64 -D WIN64")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_DEFINITIONS}" CACHE STRING "C++ Compiler Flags" FORCE)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADDITIONAL_DEFINITIONS}" CACHE STRING "C Compiler Flags" FORCE)


#     # Disable some compiler warnings for cl
#     # add_compile_options(/wd4244) # conversion from 'double' to 'x', possible loss of data
#     # add_compile_options(/wd4267) # conversion from 'size_t' to 'x', possible loss of data
#     # add_compile_options(/wd4996) # 'sprintf': This function or variable may be unsafe. Consider using sprintf_s instead.
#     # add_compile_options(/wd4305) # truncation from 'double' to 'float'
#     # add_compile_options(/wd4101) # unreferenced local variable   
#     # add_compile_options(/wd4068) # unknown pragma
#     # add_compile_options(/wd4661) 
#     # add_compile_options(/wd4477) # 'sprintf' : format string '%lu' requires an argument of type 'unsigned long'
#     # add_compile_options(/wd4804) # unsafe use of type 'bool' in operation
#     # add_compile_options(/wd4700) # uninitialized local variable used

#     set(DISABLED_WARNINGS "/wd4244 /wd4267 /wd4996 /wd4305 /wd4101 /wd4068 /wd4661 /wd4477 /wd4804 /wd4700")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DISABLED_WARNINGS}" CACHE STRING "C++ Compiler Flags" FORCE)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DISABLED_WARNINGS}" CACHE STRING "C Compiler Flags" FORCE)

# else()
#     message(FATAL_ERROR "This operating system is not supported")
# endif()