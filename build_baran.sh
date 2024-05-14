#!/bin/bash

# conda deactivate

export PATH=/usr/bin:${PATH}

cmakeExe=cmake
buildType=Release #Release or Debug
buildShared=OFF
buildDir=build-static

c_compiler=/bin/gcc
cxx_compiler=/bin/g++

# c_compiler=clang
# cxx_compiler=clang++


rm -rf ${buildDir}
mkdir -p ${buildDir}
cd ${buildDir}

${cmakeExe} \
-DCMAKE_C_COMPILER=${c_compiler} \
-DCMAKE_CXX_COMPILER=${cxx_compiler} \
-DCMAKE_BUILD_TYPE=${buildType} \
-DBUILD_SHARED_LIBS=${buildShared} \
..

${cmakeExe} --build . --config ${buildType} --target install --parallel 16

cd ..
