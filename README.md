nibrary: neuroimaging library
=============================

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![Shared build](https://github.com/baranaydogan/nibrary/actions/workflows/build_shared.yml/badge.svg)](https://github.com/baranaydogan/nibrary/actions/workflows/build_shared.yml) [![Static build](https://github.com/baranaydogan/nibrary/actions/workflows/build_static.yml/badge.svg)](https://github.com/baranaydogan/nibrary/actions/workflows/build_static.yml)


nibrary is an evolving collection of functions designed to support research in computational neuroimaging across various imaging modalities for developing new methods, tools, and applications.


### Installation

nibrary natively supports Linux, Windows and Mac. nibrary depends on several external tools. You can either install everything together with nibrary, i.e. superbuild, or install the dependencies in your system separately:

1. **Superbuild (Recommended)**

    This method will install nibrary along with all required dependencies.

    1. **Clone the repository with submodules:**

        ```bash
        git clone --recurse-submodules https://github.com/nibrary/nibrary
        ```

    2. **Build the project**

        * **Linux and Mac:** Edit and run the `build.sh` file.
        * **Windows:** Edit and run the `build_windows.sh` file.

2. **Using system libraries**

    If you prefer to use pre-installed system libraries.

    1. **Clone the repository:**

        ```bash
        git clone https://github.com/nibrary/nibrary
        ```

    2. **Install dependencies manually:**

        *   [Eigen](https://eigen.tuxfamily.org)
        *   [Geogram](https://github.com/BrunoLevy/geogram)
        *   [libigl](https://libigl.github.io/)
        *   [ProxSuite](https://github.com/Simple-Robotics/proxsuite)
        *   [SIMDe](https://github.com/simd-everywhere/simde)
        *   [zlib](http://zlib.net/)

    3. **Build the project**
    
        * **Linux and Mac:** Edit and run the `build.sh` file.
        * **Windows:** Edit and run the `build_windows.sh` file.

        *Note:* For static compilations, make sure that your system libraries are linked statically. 


### License information and use of third-party software

nibrary uses a BSD 3-Clause License. However, we use a variety of third-party tools, which are included in the [external](./external/README.md) folder. Please read the [license](./external/LICENSE.md) file under the external folder for details.




&copy; 2024 Dogu Baran Aydogan, baran.aydogan@uef.fi



