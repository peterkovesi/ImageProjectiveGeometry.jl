#=----------------------------------------------------------------------------

Image Projective Geometry

Functions supporting projective geometry for computer vision.

Copyright (c) 2016 Peter Kovesi
pk@peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK February  2016
September 2018 - Updates for v0.7/v1.0


----------------------------------------------------------------------------=#
"""
**ImageProjectiveGeometry**

Functions supporting projective geometry for computer vision.

Peter Kovesi

[peterkovesi.com](http://peterkovesi.com)


*Projective*

* Camera - Structure defining parameters of a camera.
* cameraproject - Projects 3D points into camera image.
* imagept2plane - Project image points to a plane and return their 3D locations.
* imagecorners2plane - Get the positions of image corners projected onto a plane.
* imagept2ray - Compute viewing ray corresponding to an image point.
* imagept2ray! - Compute viewing ray corresponding to an image point.
* camera2projmatrix - Generate a camera projection matrix from a Camera structure.
* decomposecamera -  Decomposition of a camera projection matrix.
* rq3 - RQ decomposition of 3x3 matrix.
* makehomogeneous - Appends a scale of 1 to an array inhomogeneous coordinates.
* makeinhomogeneous - Converts homogeneous coords to inhomogeneous coordinates.
* hnormalise - Normalises array of homogeneous coordinates to a scale of 1.
* hnormalise! - In-place normalisation of homogeneous coordinates to a scale of 1.
* homography1d - Computes 1D homography.
* homography2d - Computes 2D homography.
* solveaffine - Solve affine transformation between two sets of 2D points.
* normalise1dpts - Normalises 1D homogeneous points.
* normalise2dpts - Normalises 2D homogeneous points.
* skew - Constructs 3x3 skew-symmetric matrix from 3-vector.
* hcross - Homogeneous cross product, result normalised to s = 1.
* fundmatrix - Computes fundamental matrix from 8 or more points.
* affinefundmatrix - Computes affine fundamental matrix from 4 or more points.
* fundfromcameras - Fundamental matrix from camera matrices or structures.
* stereorectify - Rectify a stereo pair of images.
* stereorectifytransforms - Compute homographies that transform an image pair into a stereorectified pair.
* transformedimagebounds - Find where the corners of an image are transformed to by transform H and return the bounds.
* imgtrans - Homogeneous transformation of an image.
* imgtrans! - Homogeneous transformation of an image.
* idealimagepts - Ideal image points with no distortion.
* solvestereopt - Homogeneous linear solution of a stereo point.
* undistortimage - Removes lens distortion from an image.
* hline - Plot a 2D line defined in homogeneous coordinates.
* mapimage2plane! - Projects an image onto a plane in 3D.
* plotcamera - Plots graphical representation of camera(s) showing pose.

*Corner Features*

* shi_tomasi - Shi-Tomasi corner detector.
* harris - Harris corner detector.
* noble - Noble's variant of the Harris corner detector.
* coherence - Compute image coherence from structure tensor.
* structuretensor - Compute structure tensor values over an image.
* hessianfeatures  - Computes determiant of hessian features in an image.
* fastradial - Loy and Zelinski's fast radial feature detector.

*RANSAC*

* ransac - Robustly fits a model to data with the RANSAC algorithm.
* ransacfithomography - Fits 2D homography using RANSAC.
* ransacfitfundmatrix - Fits fundamental matrix using RANSAC.
* ransacfitaffinefundmatrix - Fits affine fundamental matrix using RANSAC.
* ransacfitplane - Fits plane to 3D array of points using RANSAC.
* ransacfitline - Fits line to 3D array of points using RANSAC.
* iscolinear - Are 3 points colinear.
* fitline2d - Least squares fit of a line to a set of 2D points.
* fitline3d - Fits a line to a set of 3D points.
* fitplane - Solves coefficients of plane fitted to 3 or more points.

*RANSAC demos*

* fitlinedemo - Demonstrates RANSAC line fitting.
* fitplanedemo - Demonstrates RANSAC plane fitting.
* fitfunddemo - Example of fundamental matrix computation.
* fithomogdemo - Example of finding a homography.

*Transforms*

* trans - Homogeneous transformation for a translation by x, y, z.
* rotx - Homogeneous transformation for a rotation about the x axis.
* roty - Homogeneous transformation for a rotation about the y axis.
* rotz - Homogeneous transformation for a rotation about the z axis.
* dhtrans - Computes Denavit Hartenberg matrix.
* homotrans - Homogeneous transformation of points/lines.
* invht - Inverse of a homogeneous transformation matrix.
* inveuler - Inverse of Euler transform.
* invrpy - Inverse of Roll Pitch Yaw transform.
* angleaxis - Constructs angle-axis descriptor.
* normaliseangleaxis - Normalises angle-axis descriptor.
* angleaxisrotate - Uses angle axis descriptor to rotate vectors.
* angleaxis2matrix - Converts angle-axis descriptor to 4x4 homogeneous transformation  matrix.
* matrix2angleandaxis - Decompose homogeneous matrix to angle and axis.
* matrix2angleaxis - Homogeneous matrix to angle-axis description.
* quaternion - Construct quaternion.
* quaternionconjugate - Conjugate of a quaternion.
* quaternionproduct - Computes product of two quaternions.
* quaternionrotate - Rotates a 3D vector by a quaternion.
* quaternion2matrix - Quaternion to a 4x4 homogeneous transformation matrix.
* matrix2quaternion - Homogeneous matrix to quaternion.
* vector2quaternion - Embeds 3-vector in a quaternion representation.

*Geometry*

* ray2raydist -  Minimum distance between two 3D rays.
* circleintersectionpts - Finds intersections of two circles.
* rectintersect - Test if rectangles intersect.
* pointinconvexpoly - Determine if a 2D point is within a convex polygon.
* convexpolyintersect - Test if convex polygons intersect.
* isconvexpoly - Test if 2D polygon is convex.

*Utilities*

* nonmaxsuppts - Non-maximal suppression for features/corners.
* derivative3 - 3-Tap discrete derivative filters.
* derivative5 - 5-Tap 1st and 2nd discrete derivatives.
* derivative7 - 7-Tap 1st and 2nd discrete derivatives.
* gaussfilt -  Small wrapper function for convenient Gaussian filtering.
* dilate1d - 1D morphological dilation of a signal.
* erode1d - 1D morphological erosion of a signal.
* imdilate - Image morpholgical dilation.
* imerode - Image morpholgical erosion.
* circularstruct - Generate circular structuring element for morphological operations.
* medfilt2 - Convenience wrapper for median filtering.
* stdfilt2 - Compute local standard deviation across an image.
* floatyx - Convert 2D AbstractImage to 2D float array with y x spatial order.
* histtruncate - Truncates ends of an image histogram.
* imgnormalise/imgnormalize - Normalises image values to 0-1, or to desired mean and variance.
* matchbycorrelation - Match image feature points by correlation.
* grey2census - Convert image grey scale values to census values.
* briefcoords - Compute BRIEF descriptor sampling coordinates within a patch.
* grey2lbp - Convert image grey scale values to local binary pattern.
* grey2lbp! - Convert image grey scale values to local binary pattern.
* keypause - Wait for user to hit return before continuing.

"""
module ImageProjectiveGeometry

include("projective.jl")
include("transforms3d.jl")
include("ransac.jl")
include("cornerfeatures.jl")
include("utilities.jl")
include("geometry.jl")
include("plotting.jl")

end  # module
