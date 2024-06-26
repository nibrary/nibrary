[Efficient implementation of marching cubes cases with topological guarantees - Thomas Lewiner](http://thomas.lewiner.org/publication_page.php%EF%B9%96pubkey=marching_cubes_jgt.html)

1. Authors
2. Introduction
3. Installation with the Makefile
4. Usage
5. Implementation

Authors
-------

Thomas Lewiner,
Hélio Lopes,
Antônio Wílson Vieira and
Geovan Tavares.

Introduction
------------

  Those programs are simple illustrations of the Marching Cubes
algorithm and our extension of it. We present a small graphical
interface to look at each of the 733 cases of our extended lookup
table (luttest), a commandline and a graphical interface to
arching Cubes.

Requirements
------------

  The graphical interfaces require the openGL libraries with the GLU,
glut and glui extension, available at:
    http://www.cs.unc.edu/~rademach/glui/

  The glui library we used is v2.1.

Installation with the Makefile
------------------------------

The make file has 3x2 targets:

- mcd      :
- mc       : the commandline interface to the Marching Cubes algorithm.
  
             in debug mode (mcd) or optimized mode -O3 (mc).
- mcGLd
- mcGL     : the graphical interface to the Marching Cubes algorithm.
  
             in debug mode (mcGLd) or optimized mode -O3 (mcGL).
- luttestd :
- luttest  : the graphical interface to see the 733 cases of our
  
             extended lookup table, in debug mode (luttestd) or
             optimized mode -O3 (luttest).

Usage
-----

  The luttest program simply display in 3D the table entry given by
the corresponding case, subcase and configuration.

  The mc program builds the approximation of an implicit surface at
a fixed resolution. The formula and the resolution are hard coded,
but they are easily accessible through the 'main.cpp' file.

  The mcGL program offers different kinds of input, and visualizes or
save the output as a triangulated surface. The input can be an implicit
function given by its formula, preselected, or extracted from a CSG tree.
It can also be a trilinear function given by its coordinates at the cube
vertices, in order to simulate compare the ideal surface with its
triangulation in our tiling table. It can finally be an iso archive,
simply written in binary with the resolution over the X,Y and Z axis, and
the sequence of float values. The isosurface can be generated with the
classical Marching Cubes algorithm or with our extension.
usage ./mcGL [-i file.iso] [-c file.csg] [-f 'formula'] [-r[xyz] res] -R [-si file.iv] [-sp file.ply] [-q]

Implementation
--------------

  The Marching Cubes algorithm is implemented as described in the paper :

@article{lewiner-lopes-vieira-tavares-03,
  author  = {Thomas Lewiner and Hélio Lopes and Antônio Wilson Vieira and Geovan Tavares},
  year    = 2003,
  title   = {Efficient Implementation of Marching CubesŽ Cases with Topological Guarantees},
  journal = {Journal of Graphics Tools},publisher = {A.K.Peters},
  volume  = 8,
  number  = 2,
  pages   = {1--15}
  url     = {http://www-sop.inria.fr/prisme/personnel/Thomas.Lewiner/JGT.pdf},
}

  The function parser is the one of Warp :
http://www.students.tut.fi/~warp/FunctionParser/

  The code is not I/O optimized, but it can be used as is for other
applications with the following use:

- first, set the global variables 'size_x', 'size_y' and 'size_z' to
  the resolution of your grid.

- then, call the function 'init_all()', and fill the array 'data' with
  the accessor 'set_data(float val, int I, int J, int K )'.

- the Marching Cubes algorithm can then be run on the grid by calling
  'applyMC()'. The resulting triangulated surface can be read in the
  arrays 'vertices' and 'triangles', conataining respectively 'nverts'
  and 'ntrigs' valid elements. The mesh can be exported in Greg Turk's
  PLY format or in Open Inventor / VRML 1.0 format using respectively
  'writePLY(const char *fn, bool bin = false )', and
  'writeIV (const char *fn )'.

- The intermediate data can be cleaned using 'clean_temps()', while the
  resulting arrays are cleaned with 'clean_all()'.

Have fun!

  Tomlew.


