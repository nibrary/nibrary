nibrary: neuroimaging library
=============================

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)[![Shared build](https://github.com/nibrary/nibrary/actions/workflows/build_shared.yml/badge.svg)](https://github.com/nibrary/nibrary/actions/workflows/build_shared.yml) [![Static build](https://github.com/nibrary/nibrary/actions/workflows/build_static.yml/badge.svg)](https://github.com/nibrary/nibrary/actions/workflows/build_static.yml)


nibrary is an evolving collection of functions designed to support research in computational neuroimaging across various imaging modalities for developing new methods, tools, and applications.


### Building nibrary from source

nibrary is natively supported on Linux, Windows, and macOS.

#### 1. Clone the Repository:

```bash
git clone --recurse-submodules https://github.com/nibrary/nibrary
```

#### 2. Install dependencies:

nibrary can be installed with a minimal set of standard development tools:

*   **CMake**: min v3.15.0
*   **OpenMP**
*   **C/C++ Compiler**: The following compilers have been tested and are known to work:
    *   **GCC** (min v9.0)
    *   **Clang**: v18.0.0 (v19 is known not to work)
    *   **Microsoft Visual C++ (MSVC)**: Visual Studio 2022

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

macOS:
```bash
brew install llvm@18 libomp
```

#### 3. Build

Edit and run the provided build scripts.

* **Linux:** `build_linux.sh`
* **macOS:** `build_mac.sh`
* **Windows:** `build_win.bat`

### Notes

1. nibrary bundles together the following external tools as submodules:

    *   [Eigen](https://eigen.tuxfamily.org)
    *   [Geogram](https://github.com/BrunoLevy/geogram)
    *   [libigl](https://libigl.github.io/)
    *   [ProxSuite](https://github.com/Simple-Robotics/proxsuite)
    *   [SIMDe](https://github.com/simd-everywhere/simde)
    *   [zlib](http://zlib.net/)
    *   [dcm2niix](https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage)

    To use some of the existing system libraries for installation please modify nibrary's cmake options. System libraries cannot be used for some of the dependencies.

2. For static compilations, make sure that your system libraries are static as well.


### License information and use of third-party software

nibrary uses a BSD 3-Clause License. However, we use a variety of third-party tools, which are included in the [external](./external) folder. Please read the [license](./external/LICENSE.md) file under the external folder for details.




&copy; 2025 Dogu Baran Aydogan, baran.aydogan@uef.fi



