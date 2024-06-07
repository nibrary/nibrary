#!/bin/bash

cmakeExe=cmake
buildType=Release
buildShared=OFF
buildDir=build-static

c_compiler=/bin/gcc-12
cxx_compiler=/bin/g++-12

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
