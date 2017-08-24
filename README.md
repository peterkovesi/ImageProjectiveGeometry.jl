ImageProjectiveGeometry
=======================

[![Build Status](https://travis-ci.org/peterkovesi/ImageProjectiveGeometry.jl.svg?branch=master)](https://travis-ci.org/peterkovesi/ImageProjectiveGeometry.jl)

----------------------------------------------

![banner image](doc/banner.png)

## Installation

Install via the package manager.  At the moment you should checkout
the current master rather than using Pkg.add().  While I believe the
package runs cleanly under Julia 0.5 and 0.6 there are still some
fixes I want to make before tagging a version.

```
julia> Pkg.checkout("ImageProjectiveGeometry")
```

```
help?> ImageProjectiveGeometry  # Lists a summary of the package functions 
```

## Summary

This Image Projective Geometry package is intended as a starting point
for the development of a library of projective geometry functions for
computer vision in Julia.

Currently the package consists of a number of components which
ultimately could/should be separated off into individual packages or
contributed to other existing packages.  Also, some of these
functions, no doubt, duplicate existing functions in other packages
and these should be eventually rationalised.  However at this stage,
given that Julia and its packages are still subject to some change, I
have chosen to keep all these components in this package to minimise
external dependencies and make it as self contained as possible.

## Function Reference

* [**projective**](doc/projective.md) Defines a camera structure,
implements image projection functions, functions for computing
homographies and fundamental matrices, stereo solution, etc.
* [**cornerfeatures**](doc/cornerfeatures.md) Implementations of a number
of corner detectors.  Ultimately could be contributed to the Images
package.
* [**ransac**](doc/ransac.md) A generic implementation of RANSAC along
with a collection of specific functions that use RANAC for robust
estimation of homographies and fundamental matrices, and for fitting
lines and planes etc.
* [**transforms**](doc/transforms.md) Functions for constructing,
applying, and decomposing homogeneous transforms, angle-axis
descriptors, and quaternions.
* [**utilities**](doc/utilities.md) Miscellaneous image processing functions
including nonmaximal suppression, image derivative computation and
efficient morphological dilation and erosion using rectangular and
octagonal structuring elements.  There is also a basic correlation
matcher.  Ultimately could be contributed to the Images package.
* [**geometry**](doc/geometry.md) Functions for some basic
geometric operations: minimum distance between 3D rays, intersection
of circles, convex polygons etc.
* [**ransacdemo**](doc/ransacdemo.md) Functions demonstrating the use
of ransac() to fit lines, planes, fundamental matrices and
homographies.

Also, within src, there are demo scripts for the different corner
detectors and the morphological functions.

## Contribute

* There is much that is missing.  For example there is no code for
camera calibration, computation of trifocal tensors, or bundle
adjustment.  While there is code for detecting corners there is
nothing for matching them other than a basic correlation matcher,
though feature matching probably belongs in its own package(s)

* These functions are mostly ported from MATLAB code at
 [http://www.peterkovesi.com/matlabfns](http://www.peterkovesi.com/matlabfns/index.html)
 Accordingly some of the code is still MATLABesque in nature.  There
 are, no doubt, many optimisations that could be made and type
 instabilities to be eliminated. Pull requests to make the code more Julian
 are welcome.


## Supplementary Material

* [The Fundamental Matrix Song](http://danielwedge.com/fmatrix/) by Daniel Wedge
* [The RANSAC Song](http://danielwedge.com/ransac/) by Daniel Wedge

