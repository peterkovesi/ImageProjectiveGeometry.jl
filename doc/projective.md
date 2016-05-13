Projective Function Reference
==============================

## Index

* [Camera](#Camera) - Structure defining parameters of a camera.
* [cameraproject](#cameraproject) - Projects 3D points into camera image.
* [imagept2plane](#imagept2plane) - Project image points to a plane and return their 3D locations.
* [camera2projmatrix](#camera2projmatrix) - Generate a camera projection matrix from a Camera structure.
* [decomposecamera](#decomposecamera) -  Decomposition of a camera projection matrix.
* [rq3](#rq3) - RQ decomposition of 3x3 matrix.
* [makehomogeneous](#makehomogeneous) - Appends a scale of 1 to an array inhomogeneous coordinates.
* [makeinhomogeneous](#makeinhomogeneous) - Converts homogeneous coords to inhomogeneous coordinates.
* [hnormalise](#hnormalise) - Normalises array of homogeneous coordinates to a scale of 1.
* [homography1d](#homography1d) - Computes 1D homography.
* [homography2d](#homography2d) - Computes 2D homography.
* [normalise1dpts](#normalise1dpts) - Normalises 1D homogeneous points.
* [normalise2dpts](#normalise2dpts) - Normalises 2D homogeneous points.
* [skew](#skew) - Constructs 3x3 skew-symmetric matrix from 3-vector.
* [hcross](#hcross) - Homogeneous cross product, result normalised to s = 1.
* [fundmatrix](#fundmatrix) - Computes fundamental matrix from 8 or more points.
* [affinefundmatrix](#affinefundmatrix) - Computes affine fundamental matrix from 4 or more points.
* [fundfromcameras](#fundfromcameras) - Fundamental matrix from camera matrices or structures.
* [idealimagepts](#idealimagepts) - Ideal image points with no distortion.
* [solvestereopt](#solvestereopt) - Homogeneous linear solution of a stereo point.
* [hline](#hline) - Plot a 2D line defined in homogeneous coordinates.


_________________


## Camera - Structure defining parameters of a camera.

```
      fx::Real      # Focal length.
      fy::Real
     ppx::Real      # Principal point.
     ppy::Real
      k1::Real      # Radial lens distortion parameters.
      k2::Real
      k3::Real
      p1::Real      # Tangential lens distortion parameters.
      p2::Real
    skew::Real
    rows::Integer   # Optional, providing rows and columns allows detection
    cols::Integer   # of points being projected out of image bounds.
       P::Array     # Camera position in world coordinates.
    Rc_w::Array     # Rotation matrix defining world orientation with 
                    # respect to the camera frame.
```

Note the ordering of the tangential distortion parameters p1 and p2 is
not always consistent in the literature.  Within cameraproject() they
are used in the following order, where x and y refer to column and
row coordinates repectively.

```
   dx = 2*p1*x*y          +  p2*(r^2 + 2*x^2)
   dy = p1*(r^2 + 2*y^2)  +  2*p2*x*y

```

Constructors:

```
      Camera(keyword = value, ....)

      Camera(P)     # Where P is a 3x4 projection matrix.  
                    # In this case only fx, fy, skew, P and Rc_w can be defined, 
                    # all other parameters are set to 0.
```

## cameraproject - Projects 3D points into camera image.

```
Usage 1:  xy = cameraproject(P, pt)

Arguments:
             P - 3x4 camera projection matrix.
             C - Camera structure.
            pt - 3xN matrix of 3D points to project into the image.

Returns:  
            xy - 2xN matrix of projected image positions.



Usage 2:  (xy, visible) = cameraproject(C, pt)

Arguments: 
             C - Camera structure.
            pt - 3xN matrix of 3D points to project into the image.
Returns: 
       xy      - 2xN matrix of projected image positions
       visible - Array of values 1/0 indicating whether the point is
                 within the field of view.  This is only evaluated if
                 the camera structure has non zero values for its
                 'rows' and 'cols' fields. Otherwise an empty matrix
                 is returned.
```

See also: Camera(), camstruct2projmatrix()


## imagept2plane - Project image points to a plane and return their 3D locations

```
Usage:  pt = image2plane(C, xy, planeP, planeN)

Arguments:  
         C - Camera structure, Alternatively
             C can be a 3x4 camera projection matrix.
        xy - Image points specified as 2 x N array (x,y) / (col,row)
    planeP - Some point on the plane.
    planeN - Plane normal vector.

Returns:
        pt - 3xN array of 3D points on the plane that the image points
             correspond to. 
```

Note that the plane is specified in terms of the world frame by defining a
3D point on the plane, planeP, and a surface normal, planeN.

Lens distortion is handled by using the standard lens distortion
parameters assuming locally that the distortion is constant, computing
the forward distortion and then subtracting the distortion.

See also Camera(), cameraproject()


## camera2projmatrix - Generate a camera projection matrix from a Camera structure.

```
Usage:     P = camera2projmatrix(C)

Argument:  C - Camera structure

Returns:   P - 3x4 camera projection matrix that maps homogeneous 3D world 
           coordinates to homogeneous image coordinates.
```

Function takes a camera structure and returns its equivalent projection matrix
ignoring lens distortion parameters etc.

See also: Camera(), projmatrix2camera(), cameraproject()


## decomposecamera -  Decomposition of a camera projection matrix.

```
Usage:  K, Rc_w, Pc, pp, pv = decomposecamera(P)

   P is decomposed into the form P = K*(R -R*Pc)

Argument:  P - 3 x 4 camera projection matrix
Returns:   
           K - Calibration matrix of the form
                 |  ax   s   ppx |
                 |   0   ay  ppy |
                 |   0   0    1  |

               Where: 
               ax = f/pixel_width and ay = f/pixel_height,
               ppx and ppy define the principal point in pixels,
               s is the camera skew.
        Rc_w - 3 x 3 rotation matrix defining the world coordinate frame
               in terms of the camera frame. Columns of R transposed define
               the directions of the camera X, Y and Z axes in world
               coordinates. 
          Pc - Camera centre position in world coordinates.
          pp - Image principal point.
          pv - Principal vector  from the camera centre C through pp
               pointing out from the camera.  This may not be the same as  
               R'(:,3) if the principal point is not at the centre of the
               image, but it should be similar. 
```

See also: rq3()


## rq3 - RQ decomposition of 3x3 matrix.

```
Usage:  R, Q = rq3(A)

Argument:  A - 3 x 3 matrix.
Returns:   R - Upper triangular 3 x 3 matrix
           Q - 3 x 3 orthonormal rotation matrix
               Such that  R*Q = A.
```

The signs of the rows and columns of R and Q are chosen so that the diagonal
elements of R are +ve.

See also: decomposecamera()


## makehomogeneous - Appends a scale of 1 to an array inhomogeneous coordinates.

```
Usage:  hx = makehomogeneous(x)

Argument:
        x  - an N x npts array of inhomogeneous coordinates.

Returns:
        hx - an (N+1) x npts array of homogeneous coordinates with the
             homogeneous scale set to 1
```

See also: makeinhomogeneous(), hnormalise()


## makeinhomogeneous - Converts homogeneous coords to inhomogeneous coordinates.

```
Usage:  x = makehomogeneous(hx)

Argument:
        hx  - an N x npts array of homogeneous coordinates.

Returns:
        x - an (N-1) x npts array of inhomogeneous coordinates
```

Warning:  If there are any points at infinity (scale = 0) the coordinates
of these points are simply returned minus their scale coordinate.

See also: makehomogeneous(), hnormalise()


## hnormalise - Normalises array of homogeneous coordinates to a scale of 1.

```
Usage:  nx = hnormalise(x)

Argument:
        x  - an Nxnpts array of homogeneous coordinates.

Returns:
        nx - an Nxnpts array of homogeneous coordinates rescaled so
             that the scale values nx(N,:) are all 1.
```

Note that any homogeneous coordinates at infinity (having a scale value of
0) are left unchanged.


## homography1d - Computes 1D homography.

```
Usage:   H = homography1d(x1, x2)

Arguments:
         x1  - 2xN set of homogeneous points
         x2  - 2xN set of homogeneous points such that x1<->x2
Returns:
          H - the 2x2 homography such that x2 = H*x1
```

## homography2d - Computes 2D homography.

```
Usage 1:   H = homography2d(x1, x2)

Arguments:
         x1  - 3xN set of homogeneous points
         x2  - 3xN set of homogeneous points such that x1<->x2


Usage 2:    H = homography2d(x)
Argument:
          x  - If a single argument is supplied it is assumed that it
               is in the form x = [x1; x2]  

Returns:
         H - the 3x3 homography such that x2 = H*x1
```

Usage 2 is intended for use with ransac()


## normalise1dpts - Normalises 1D homogeneous points.

Function translates and normalises a set of 1D homogeneous points 
so that their centroid is at the origin and their mean distance from 
the origin is 1.  

```
Usage:   (newpts, T) = normalise1dpts(pts)

Argument:
  pts -  2xN array of 2D homogeneous coordinates

Returns:
  newpts -  2xN array of transformed 1D homogeneous coordinates
  T      -  The 2x2 transformation matrix, newpts = T*pts
```           

Note that if one of the points is at infinity no normalisation
is possible.  In this case a warning is printed and pts is
returned as newpts and T is the identity matrix.

See also: normalise2dpts()


## normalise2dpts - Normalises 2D homogeneous points.

```
Usage:   (newpts, T) = normalise2dpts(pts)

Argument:
     pts -  3xN array of 2D homogeneous coordinates.

Returns:
  newpts -  3xN array of transformed 2D homogeneous coordinates.  The
            scaling parameter is normalised to 1 unless the point is at
            infinity. 
  T      -  The 3x3 transformation matrix, newpts = T*pts.
```          

Function translates and normalises a set of 2D homogeneous points 
so that their centroid is at the origin and their mean distance from 
the origin is sqrt(2).  This process typically improves the
conditioning of any equations used to solve homographies, fundamental
matrices etc.

If there are some points at infinity the normalisation transform
is calculated using just the finite points.  Being a scaling and
translating transform this will not affect the points at infinity.

## skew - Constructs 3x3 skew-symmetric matrix from 3-vector.

```
Usage:  s = skew(v)

Argument:  v - 3-vector
Returns:   s - 3x3 skew-symmetric matrix
```

The cross product between two vectors, a x b can be implemented as a matrix
product  skew(a)*b

## hcross - Homogeneous cross product, result normalised to s = 1.

Function to form cross product between two points, or lines,
in homogeneous coodinates.  The result is normalised to lie
in the scale = 1 plane.

``` 
Usage: c = hcross(a,b)

Arguments:  a, b  - 3x1 arrays or vectors
Returns:       c  - 3-vector
```

## fundmatrix - Computes fundamental matrix from 8 or more points.

Function computes the fundamental matrix from 8 or more matching points in
a stereo pair of images.  The normalised 8 point algorithm given by
Hartley and Zisserman p265 is used.  To achieve accurate results it is
recommended that 12 or more points are used

```
Usage:   (F, e1, e2) = fundmatrix(x1, x2)
         (F, e1, e2) = fundmatrix(x)

Arguments:
         x1, x2 - Two sets of corresponding 3xN set of homogeneous
         points.
        
         x      - If a single argument is supplied it is assumed that it
                  is in the form x = [x1; x2]
Returns:
         F      - The 3x3 fundamental matrix such that x2'*F*x1 = 0.
         e1     - The epipole in image 1 such that F*e1 = 0
         e2     - The epipole in image 2 such that F'*e2 = 0
```

Usage with a single argument is intended for the use of this function
with ransac() 

See also: affinefundmatrix()


## affinefundmatrix - Computes affine fundamental matrix from 4 or more points.

Function computes the affine fundamental matrix from 4 or more matching
points in a stereo pair of images.  

```
Usage:   (F, e1, e2) = affinefundmatrix(x1, x2)
         (F, e1, e2) = affinefundmatrix(x)

Arguments:
         x1, x2 - Two sets of corresponding points defined as
                  2xN arrays of inhomogeneous image coordinates.
         
         x      - If a single argument is supplied it is assumed that it
                  is in the form x = [x1; x2]
Returns:
         F      - The 3x3 fundamental matrix such that x2'*F*x1 = 0.
         e1     - The epipole in image 1 such that F*e1 = 0
         e2     - The epipole in image 2 such that F'*e2 = 0
```
Usage with a single argument is intended for the use of this function
with ransac() 

The Gold Standard algorithm given by Hartley and Zisserman p351 (2nd
Ed.) is used.

See also: fundmatrix()


## fundfromcameras - Fundamental matrix from camera matrices or structures.

```
Usage: F = fundfromcameras(P1, P2)
       F = fundfromcameras(C1, C2)

Arguments:  P1, P2 - Two 3x4 camera projection matrices or
            C1, C2 - Two Camera structres

Returns:    F      - Fundamental matrix relating the two camera views.
``

See also: fundmatrix(), affinefundmatrix(), Camera()


## idealimagepts - Ideal image points with no distortion.

```
Usage:  xyideal = idealimagepts(C, xy)

Arguments:  
         C - Camera structure
        xy - Image points specified as 2 x N array (x,y) / (col,row)

Returns:
   xyideal - Ideal image points.  These points correspond to the image
             locations that would be obtained if the camera had no lens
             distortion.  That is, if they had been projected using an ideal
             projection matrix computed from fx, fy, ppx, ppy, skew.
```

See also Camera(), cameraproject()


## solvestereopt - Homogeneous linear solution of a stereo point.

```
Usage:  (pt, xy_reproj) = solvestereopt(xy, P)
        (pt, xy_reproj) = solvestereopt(xy, C)

Multiview stereo: Solves 3D location of a point given image coordinates of
that point in two, or more, images.

Arguments:    xy - 2xN matrix of x, y image coordinates, one column for
                   each camera.
               P - N array of corresponding 3x4 image projection
                   matrices, or
               C - an N array of Camera structures

Returns:      pt - 3D location in space returned in normalised
                   homogeneous coordinates (a 4-vector with last element = 1)
       xy_reproj - 2xN matrix of reprojected image coordinates.
```

See also: idealimagepts(), camstruct2projmatrix(), Camera()


## hline - Plot a 2D line defined in homogeneous coordinates.

Function for ploting 2D homogeneous lines defined by 2 points
or a line defined by a single homogeneous vector

```
Usage 1:   hline(p1,p2, linestyle="b-")
Arguments: p1, p2 -  Two 3-vectors defining points in homogeneous 
                     coordinates

Usage 2:    hline(l, linestyle="b-")  
Argument:   l - A 3-vector defining a line in homogeneous coordinates
```

Note that in the case where a homogeneous line is supplied as the
argument the extent of the line drawn depends on the current axis
limits.  This will require you to set the desired limits with a call
to PyPlot.axis() prior to calling this function.

Side effect: PyPlot hold() state will be set to true