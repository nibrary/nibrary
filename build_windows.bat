@SET cmakeExe=cmake
@SET buildType=Release

@rmdir build_windows_shared /s /q
@mkdir build_windows_shared
@cd build_windows_shared

@%cmakeExe% -DCMAKE_C_COMPILER=%c_compiler% -DCMAKE_CXX_COMPILER=%cxx_compiler% -DCMAKE_BUILD_TYPE=%buildType% -DBUILD_SHARED_LIBS=ON ..
@%cmakeExe% --build . --config %buildType% --target install -- /m:10

@cd ..
