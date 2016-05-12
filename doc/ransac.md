RANSAC Function Reference
==========================

## Index

* [ransac](#ransac) - Robustly fits a model to data with the RANSAC algorithm.
* [ransacfithomography](#ransacfithomography) - Fits 2D homography using RANSAC.
* [ransacfitfundmatrix](#ransacfitfundmatrix) - Fits fundamental matrix using RANSAC.
* [ransacfitaffinefundmatrix](#ransacfitaffinefundmatrix) - Fits affine fundamental matrix using RANSAC.
* [ransacfitplane](#ransacfitplane) - Fits plane to 3D array of points using RANSAC.
* [ransacfitline](#ransacfitline) - Fits line to 3D array of points using RANSAC.
* [iscolinear](#iscolinear) - Are 3 points colinear.
* [fitline2d](#fitline2d) - Least squares fit of a line to a set of 2D points.
* [fitline3d](#fitline3d) - Fits a line to a set of 3D points.
* [fitplane](#fitplane) - Solves coefficients of plane fitted to 3 or more points.



_______________________


## ransac - Robustly fits a model to data with the RANSAC algorithm.

```
Usage:

(M, inliers) = ransac(x, fittingfn, distfn, degenfn s, t, feedback,
                      maxDataTrials, maxTrials)

Arguments:
    x         - Data sets to which we are seeking to fit a model M
                It is assumed that x is of size [d x Npts]
                where d is the dimensionality of the data and Npts is
                the number of data points.

    fittingfn - Reference to a function that fits a model to s
                data from x.  It is assumed that the function is of the
                form: 
                   M = fittingfn(x)
                Note it is possible that the fitting function can return
                multiple models (for example up to 3 fundamental matrices
                can be fitted to 7 pairs of matched points).  In this 
                case it is assumed that the fitting function returns an
                array of models.
                If this function cannot fit a model it should return M as
                an empty matrix.

    distfn    - Reference to a function that evaluates the
                distances from the model to data x.
                It is assumed that the function is of the form:
                   (inliers, M) = distfn(M, x, t)
                This function must evaluate the distances between points
                and the model returning the indices of elements in x that
                are inliers, that is, the points that are within distance
                't' of the model.  Additionally, if M is an array of
                possible models distfn() will return the model that has the
                most inliers.  If there is only one model this function
                must still copy the model to the output.  After this call M
                will be a single object representing only one model. 

    degenfn   - Reference to a function that determines whether a
                set of datapoints will produce a degenerate model.
                This is used to discard random samples that do not
                result in useful models.
                It is assumed that degenfn is a boolean function of
                the form: 
                   r = degenfn(x)
                It may be that you cannot devise a test for degeneracy in
                which case you should write a dummy function that always
                returns a value of true and rely on fittingfn() to return
                an empty model should the data set be degenerate.

    s         - The minimum number of samples from x required by
                fittingfn to fit a model.

    t         - The distance threshold between a data point and the model
                used to decide whether the point is an inlier or not.

    feedback  - An optional boolean flag. If set to true the trial count and the
                estimated total number of trials required is printed out at
                each step.  Defaults to false.

maxDataTrials - Maximum number of attempts to select a non-degenerate
                data set. This parameter is optional and defaults to 100.

    maxTrials - Maximum number of iterations. This parameter is optional and
                defaults to 1000.

    p         - Desired probability of choosing at least one sample
                free from outliers, defaults to 0.99

Returns:
    M         - The model having the greatest number of inliers.
    inliers   - An array of indices of the elements of x that were
                the inliers for the best model.
```

For an example of the use of this function see ransacfithomography() or
ransacfitplane() 


References:

   M.A. Fishler and  R.C. Boles. "Random sample concensus: A paradigm
   for model fitting with applications to image analysis and automated
   cartography". Comm. Assoc. Comp, Mach., Vol 24, No 6, pp 381-395, 1981

   Richard Hartley and Andrew Zisserman. "Multiple View Geometry in
   Computer Vision". pp 101-113. Cambridge University Press, 2001


## ransacfithomography - Fits 2D homography using RANSAC.

```
Usage:   (H, inliers) = ransacfithomography(x1, x2, t)

Arguments:
         x1  - 2xN or 3xN set of homogeneous points.  If the data is
               2xN it is assumed the homogeneous scale factor is 1.
         x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
         t   - The distance threshold between data point and the model
               used to decide whether a point is an inlier or not. 
               Note that point coordinates are normalised to that their
               mean distance from the origin is sqrt(2).  The value of
               t should be set relative to this, say in the range 
               0.001 - 0.01  
```

Note that it is assumed that the matching of x1 and x2 are putative and it
is expected that a percentage of matches will be wrong.

```
Returns:
         H       - The 3x3 homography such that x2 = H*x1.
         inliers - An array of indices of the elements of x1, x2 that were
                   the inliers for the best model.
```

See Also: ransac(), homography2d(), homography1d()


## ransacfitfundmatrix - Fits fundamental matrix using RANSAC.

```
Usage:   (F, inliers) = ransacfitfundmatrix(x1, x2, t)

Arguments:
         x1  - 2xN or 3xN set of homogeneous points.  If the data is
               2xN it is assumed the homogeneous scale factor is 1.
         x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
         t   - The distance threshold between data point and the model
               used to decide whether a point is an inlier or not. 
               Note that point coordinates are normalised to that their
               mean distance from the origin is sqrt(2).  The value of
               t should be set relative to this, say in the range 
               0.001 - 0.01  
```

Note that it is assumed that the matching of x1 and x2 are putative and it
is expected that a percentage of matches will be wrong.

```
Returns:
         F       - The 3x3 fundamental matrix such that x2'Fx1 = 0.
         inliers - An array of indices of the elements of x1, x2 that were
                   the inliers for the best model.
```

See also: ransac(), fundmatrix()


## ransacfitaffinefundmatrix - Fits affine fundamental matrix using RANSAC.

```
Usage:   (F, inliers) = ransacfitaffinefundmatrix(x1, x2, t)

Arguments:
         x1  - 2xN or 3xN set of homogeneous points.  If the data is
               2xN it is assumed the homogeneous scale factor is 1.
         x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
         t   - The distance threshold between data point and the model
               used to decide whether a point is an inlier or not. 
               Note that point coordinates are normalised to that their
               mean distance from the origin is sqrt(2).  The value of
               t should be set relative to this, say in the range 
               0.001 - 0.01  
```

Note that it is assumed that the matching of x1 and x2 are putative and it
is expected that a percentage of matches will be wrong.

```
Returns:
         F       - The 3x3 fundamental matrix such that x2'Fx1 = 0.
         inliers - An array of indices of the elements of x1, x2 that were
                   the inliers for the best model.
```

See also: ransac(), fundmatrix(), affinefundmatrix()


## ransacfitplane - Fits plane to 3D array of points using RANSAC.

```
Usage  (B, P, inliers) = ransacfitplane(XYZ, t, feedback)

This function uses the RANSAC algorithm to robustly fit a plane
to a set of 3D data points.

Arguments:
         XYZ - 3xNpts array of xyz coordinates to fit plane to.
         t   - The distance threshold between data point and the plane
               used to decide whether a point is an inlier or not.
         feedback - Optional flag 0 or 1 to turn on RANSAC feedback
                    information.

Returns:
          B - 4x1 array of plane coefficients in the form
              b[1]*X + b[2]*Y +b[3]*Z + b[4] = 0
              The magnitude of B is 1.
              This plane is obtained by a least squares fit to all the
              points that were considered to be inliers, hence this
              plane will be slightly different to that defined by P below.
          P - The three points in the data set that were found to
              define a plane having the most number of inliers.
              The three columns of P defining the three points
    inliers - The indices of the points that were considered
              inliers to the fitted plane.
```

See also:  ransac(), fitplane()


## ransacfitline - Fits line to 3D array of points using RANSAC.

```
Usage  (L, inliers) = ransacfitline(XYZ, t, feedback)

This function uses the RANSAC algorithm to robustly fit a line to a
set of 3D data points or to a set of 2D homogeneous points with scale
value of 1

Arguments:
         XYZ - 3xNpts array of xyz coordinates to fit line to.
         t   - The distance threshold between data point and the line
               used to decide whether a point is an inlier or not.
    feedback - Optional boolean flag to turn on RANSAC feedback
               information.

Returns:.
          V - Line obtained by a simple fitting on the points that
              are considered inliers.  The line goes through the
              calculated mean of the inlier points, and is parallel to
              the principal eigenvector.  The line is scaled by the
              square root of the largest eigenvalue.
              This line is returned as a nx2 matrix.  The first column 
              is the beginning point, the second column is the end point 
              of the line.
          L - The two points in the data set that were found to
              define a line having the most number of inliers.
              The two columns of L defining the two points.
    inliers - The indices of the points that were considered
              inliers to the fitted line.
```

See also:  ransac(), fitline3d(), ransacfitplane()


## iscolinear - Are 3 points colinear.

```
Usage:  r = iscolinear(p1, p2, p3, flag)

Arguments:
       p1, p2, p3 - Points in 2D or 3D.
       homog      - An optional boolean flag indicating that p1, p2, p3 
                    are homogneeous coordinates with arbitrary scale.  
                    The default is false, that is it is assumed that the
                    points are inhomogeneous, or that they are homogeneous 
                    ewith qual scale.

Returns:
       r  -  true if points are co-linear, false otherwise
```

## fitline2d - Least squares fit of a line to a set of 2D points.

```
Usage:   (C, dist) = fitline2d(XY)

Where:   XY  - 2xNpts array of xy coordinates to fit line to data of
               the form 
               [x1 x2 x3 ... xN
                y1 y2 y3 ... yN]
                
               XY can also be a 3xNpts array of homogeneous coordinates.

Returns: C    - 3x1 array of line coefficients in the form
                c[1]*X + c[2]*Y + c[3] = 0
         dist - Array of distances from the fitted line to the supplied
                data points.  Note that dist is only calculated if the
                function is called with two output arguments.
```

The magnitude of C is scaled so that line equation corresponds to

```
   sin(theta)*X + (-cos(theta))*Y + rho = 0
```

where theta is the angle between the line and the x axis and rho is the
perpendicular distance from the origin to the line.  Rescaling the
coefficients in this manner allows the perpendicular distance from any
point (x,y) to the line to be simply calculated as

```
   r = abs(c[1]*X + c[2]*Y + c[3])
```

If you want to convert this line representation to the classical form 

```
     Y = a*X + b
use
  a = -c[1]/c[2]
  b = -c[3]/c[2]
```

Note, however, that this assumes c[2] is not zero


## fitline3d - Fits a line to a set of 3D points.

```
Usage:   L = fitline3d(XYZ)

Where: XYZ - 3xNpts array of XYZ coordinates
              [x1 x2 x3 ... xN;
               y1 y2 y3 ... yN;
               z1 z2 z3 ... zN]

Returns: L - 3x2 matrix consisting of the two endpoints of the line
             that fits the points.  The line is centered about the
             mean of the points, and extends in the directions of the
             principal eigenvectors, with scale determined by the
             eigenvalues.
```

## fitplane - Solves coefficients of plane fitted to 3 or more points.

```
Usage:   B = fitplane(XYZ)

Where:   XYZ - 3xNpts array of xyz coordinates to fit plane to.   
               If Npts is greater than 3 a least squares solution 
               is generated.

Returns: B   - 4x1 array of plane coefficients in the form
               B[1]*X + B[2]*Y + B[3]*Z + B[4] = 0
               The magnitude of B is 1.
```

See also: ransacfitplane()