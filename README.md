nibrary: neuroimaging library
=============================

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)[![Shared build](https://github.com/nibrary/nibrary/actions/workflows/build_shared.yml/badge.svg)](https://github.com/nibrary/nibrary/actions/workflows/build_shared.yml) [![Static build](https://github.com/nibrary/nibrary/actions/workflows/build_static.yml/badge.svg)](https://github.com/nibrary/nibrary/actions/workflows/build_static.yml)


nibrary is an evolving collection of functions designed to support research in computational neuroimaging across various imaging modalities for developing new methods, tools, and applications. 


### Installation

nibrary is natively supported on Linux, Windows, and macOS. nibrary can be installed with a minimal set of standard development tools.

*   **CMake**: min v3.15.0
*   **OpenMP**
*   **C/C++ Compiler**: The following compilers have been tested and are known to work:
    *   **GCC** (min v9.0)
    *   **Clang**: v18.0.0 (v19 is known not to work)
    *   **Microsoft Visual C++ (MSVC)**: Visual Studio 2022

To install nibrary, first clone the repository with all the submodules as shown below. 

```bash
git clone --recurse-submodules https://github.com/nibrary/nibrary
```

Then follow the instructions below to build nibrary from source for different operating systems.


#
### Linux

#### 1. Install dependencies:

Debian/Ubuntu:
```bash
sudo apt install cmake libomp-dev build-essential
```

Fedora/CentOS/RHEL:
```bash
sudo dnf install cmake openmp-devel
sudo dnf group install "Development Tools"
```

Arch Linux/Manjaro:
```bash
sudo pacman -S cmake openmp
sudo pacman -S base-devel
```

#### 2. Run the build script:

The following will install a statically built library under the `build-static` folder using the default compiler in the system.
```bash
cd nibrary
sh build_linux.sh
```


Edit the `build_linux.sh` script to customize your installation.

#
### macOS

The provided installation script for macOS, `build_mac.sh`, will install the dependencies, set the environment variables, compile and install nibrary under the `build-static` folder.

```bash
cd nibrary
sh build_linux.sh
```

Note that nibrary requires `llvm` and `libomp` in macOS. We have successfully tested `llvm` v18 to build nibrary; however, compilation with the more recent version, v19, fails. The provided build script installs `llvm@18` and `libomp` in the system. However, it does not permanently set environment variables for future use. The `build_mac.sh` script contains information about how to permanently set these environment variables if needed.

Edit the `build_mac.sh` script to customize your installation.

#
### Windows

Install Visual Studio 2022 (other versions might work too but they have not been tested). Open command window and use the following to install a statically built library under the `build-static` folder:

```cmd
cd nibrary
call build_win.bat
```

Edit the `build_win.bat` script to customize your installation.

#
### Third-Party software

nibrary utilizes the following third-party libraries, which are included as Git submodules and are built *locally* as part of the nibrary build process:

*   [Eigen](https://eigen.tuxfamily.org)
*   [Geogram](https://github.com/BrunoLevy/geogram)
*   [libigl](https://libigl.github.io/)
*   [ProxSuite](https://github.com/Simple-Robotics/proxsuite)
*   [SIMDe](https://github.com/simd-everywhere/simde)
*   [zlib](http://zlib.net/)
*   [dcm2niix](https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage)

**Important:** These dependencies are built within the nibrary source tree and do not affect or modify your system's existing installations. They will co-exist peacefully with any previously installed versions of these libraries. You do not need to install these libraries separately.

#
### License

nibrary is licensed under the [BSD 3-Clause License](LICENSE.md). Third-party software included in the `external` directory have their own respective licenses, detailed in the [`external/LICENSE.md`](external/LICENSE.md) file.




&copy; 2025 Dogu Baran Aydogan, baran.aydogan@uef.fi



