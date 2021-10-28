Delaunay Triangulation
======================

This repository features a header-only library for constructing 2D delaunay triangulations.
It uses [Eigen 3](https://eigen.tuxfamily.org/index.php) as backend and supports both 32-bit and 64-bit floating
point precision (and in fact whatever other types for which Eigen implements `determinant` and basic vector operations).

The library implements a triangle primitive based variant of divide-and-conquer algorithm proposed by
Guibas & Stolfi in *Primitives for the Manipulation of General Subdivisions and the Computation of Voronoi Diagrams*
(1983) (DOI:10.1145/282918.282923). More info on the triangle data structure can be found [here](
http://www.cs.cmu.edu/~quake/tripaper/triangle2.html).

Usage:
------

Linux (apt package manager):
```
sudo apt install libeigen3-dev
sudo apt install libopencv-dev
git clone git@github.com:Lehdari/Delaunay.git
cd Delaunay
mkdir build
cd build
cmake ..
make
./delaunay_demo
```
OpenCV is only required for building the demo. Demo build can be disabled from CMake.
As it's a header only library, copy it directly to your project or use as a submodule.

The interface is trivial: pass a vector of your favourite 2D vectors to `delaunayTriangulate` and it spits out
vector of index triplets indicating triangles. For the convex edges, third index of a triplet is -1.
