#=--------------------------------------------------------------------

projective - Functions supporting projective geometry for computer vision.

Part of the ImageProjectiveGeometry Module

Copyright (c) 2016 Peter Kovesi
pk@peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK February 2016

---------------------------------------------------------------------=#

export Camera
export cameraproject, camera2projmatrix, imagept2plane, decomposecamera, rq3
export makehomogeneous, makeinhomogeneous, hnormalise
export homography1d, homography2d, normalise1dpts, normalise2dpts 
export skew, hcross
export fundmatrix, affinefundmatrix, fundfromcameras
export idealimagepts, solvestereopt
export hline, plotcamera

import PyPlot

#=
imTrans.m 
imTransD.m
equalAngleConstraint.m 
knownAngleConstraint.m 
lengthRatioConstraint 
homoTrans 2D homogeneous transformation of points/lines.
undistortimage.m
=#

#------------------------------------------------------------------------
"""
Camera - Structure defining parameters of a camera
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

"""

type Camera{T<:AbstractFloat}
    fx::T           # Focal length.
    fy::T
    ppx::T          # Principal point.
    ppy::T
    k1::T           # Radial lens distortion.
    k2::T
    k3::T
    p1::T           # Tangential lens distortion.
    p2::T
    skew::T
    rows::Int       # Optional, providing rows and columns allows detection
    cols::Int       # of points being projected out of image bounds.
    P::Vector{T}    # Camera position in world coordinates.
    Rc_w::Matrix{T} # Rotation matrix defining world orientation with
                    # respect to the camera frame.
end

# Keyword-value constructor
function Camera(;fx=1.0, fy=1.0, ppx=0.0, ppy=0.0,
                k1=0.0, k2=0.0, k3=0.0, p1=0.0, p2=0.0, skew=0.0,
                rows=0, cols=0,
                P=[0.0, 0.0, 0.0], Rc_w=eye(3))
    if size(P) != (3,)
        error("Camera position must be a 3x1 array")
    end

    if size(Rc_w) != (3,3)
        error("Camera orientation must be a 3x3 rotation matrix")
    end

     # Get the promotion type
     args = promote(fx,fy,ppx,ppy,k1,k2,k3,p1,p2,skew,P[1],Rc_w[1,1])
     T = eltype(args)

     return Camera{T}(fx, fy, ppx, ppy, k1, k2, k3, p1, p2, skew, rows, cols, P, Rc_w)
end

# Contructor that takes a projection matrix
function Camera{T<:AbstractFloat}(P::Matrix{T})
    if size(P) != (3,4)
        error("Projection matrix must be 3x4")
    end

    (K, Rc_w, Pc, pp, pv) = decomposecamera(P)
    fx = K[1,1]
    fy = K[2,2]
    skew = K[1,2]
    k1 = 0.0
    k2 = 0.0
    k3 = 0.0
    p1 = 0.0
    p2 = 0.0
    rows = 0
    cols = 0
    Camera{T}(fx, fy, pp[1], pp[2], k1, k2, k3, p1, p2, skew, rows, cols, Pc, Rc_w)
end


#-------------------------------------------------------------------------
"""
cameraproject - Projects 3D points into camera image 
```
Usage:  xy = cameraproject(P, pt)

Arguments:
             P - 3x4 camera projection matrix.
            pt - 3xN matrix of 3D points to project into the image.

Returns:  
            xy - 2xN matrix of projected image positions.
```
See also: Camera(), camstruct2projmatrix()

"""
#=
PK September 2008
        June 2015    - Ported to Julia
=#
# Projection matrix version

function cameraproject(P::Array, pta::Array)

    pt = copy(pta) # Code should be altered so that this is not needed

    if size(P) != (3,4)
        error("Projection matrix must be 3x4")
    end

    rows = size(pt,1)
    if rows != 3
        error("Points must be in a 3xN array")
    end
    
    xy = P*makehomogeneous(pt)
    xy = makeinhomogeneous(xy)
    visible = []

    return xy, visible
end    

#--------------------------------------------------------------------------
"""
cameraproject - Projects 3D points into camera image 
```
Usage:  (xy, visible) = cameraproject(C, pt)

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
"""

# Camera structure version

function cameraproject(C::Camera, pta::Array)

    pt = copy(pta) # Make local copy
    rows = size(pt,1)
    
    if rows != 3
        error("Points must be in a 3xN array")
    end
    
    # Transformation of a ground point from world coordinates to camera coords
    # can be thought of as a translation, then rotation as follows
    #
    #   Gc =  Tc_p *  Tp_w * Gw    
    #
    #  | . . .   | | 1     -x | |Gx|
    #  | Rc_w    | |   1   -y | |Gy|
    #  | . . .   | |     1 -z | |Gz|
    #  |       1 | |        1 | | 1|
    #      
    # Subscripts:
    #   w - world frame
    #   p - world frame translated to camera origin
    #   c - camera frame

    # First translate world frame origin to match camera origin.  This is
    # the Tp_w bit.  
    for n = 1:3
        pt[n,:] = pt[n,:] - C.P[n] 
    end
    
    # Then rotate to camera frame using Rc_w
    pt = C.Rc_w*pt
    
    # Follow classical projection process    
    # Generate normalized coords
    x_n = pt[1:1,:]./pt[3:3,:]
    y_n = pt[2:2,:]./pt[3:3,:]
    rsqrd = x_n.^2 + y_n.^2  # radius squared from centre
    
    # Radial distortion factor
    r_d = 1.0 + C.k1*rsqrd + C.k2*rsqrd.^2 + C.k3*rsqrd.^3
    
    # Tangential distortion component
    dtx = 2*C.p1*x_n.*y_n         + C.p2*(rsqrd + 2*x_n.^2)
    dty = C.p1*(rsqrd + 2*y_n.^2) + 2*C.p2*x_n.*y_n
    
    # Obtain lens distorted coordinates
    x_d = r_d.*x_n + dtx
    y_d = r_d.*y_n + dty   
    
    # Finally project to pixel coordinates
    # Note skew represents the 2D shearing coefficient times fx
    x_p = C.ppx + x_d*C.fx + y_d*C.skew
    y_p = C.ppy + y_d*C.fy
    
    xy = [x_p
	  y_p]
    
    # If C.rows and C.cols not empty determine points that are within image bounds
    if C.rows !=0 && C.cols != 0
        visible = (x_p .>= 1.0) & (x_p .<= convert(Float64,C.cols)) &
                  (y_p .>= 1.0) & (y_p .<= convert(Float64,C.rows))
    else
        visible = []
    end

    return xy, visible
end

#--------------------------------------------------------------------------
"""
imagept2plane - Project image points to a plane and return their 3D locations
```
Usage:  pt = imagept2plane(C, xy, planeP, planeN)

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
"""
#=
PK May 2015
=#

function imagept2plane(C::Camera, xy::Array, planeP::Vector, planeN::Vector)

    dim = size(xy,1)
    if isa(xy, Vector)
        N = 1
    else    
        N = size(xy,2)
    end

    if dim != 2
        error("Image xy data must be a 2 x N array")
    end
    
    # Reverse the projection process as used by cameraproject()

    # Subtract principal point and divide by focal length to get normalised,
    # distorted image coordinates.  Note skew represents the 2D shearing
    # coefficient times fx
    y_d = (xy[2:2,:] - C.ppy)/C.fy
    x_d = (xy[1:1,:] - C.ppx - y_d*C.skew)/C.fx

    # Radial distortion factor.  Here the squared radius is computed from the
    # already distorted coordinates.  The approximation we are making here is to
    # assume that the distortion is locally constant.
    rsqrd = x_d.^2 + y_d.^2
    r_d = 1 + C.k1*rsqrd + C.k2*rsqrd.^2 + C.k3*rsqrd.^3
    
    # Tangential distortion component, again computed from the already distorted
    # coords.
    dtx = 2*C.p1*x_d.*y_d         + C.p2*(rsqrd + 2*x_d.^2)
    dty = C.p1*(rsqrd + 2*y_d.^2) + 2*C.p2*x_d.*y_d
    
    # Subtract the tangential distortion components and divide by the radial
    # distortion factor to get an approximation of the undistorted normalised
    # image coordinates (with no skew)
    x_n = (x_d - dtx)./r_d
    y_n = (y_d - dty)./r_d

    # Define a set of points at the normalised distance of z = 1 from the
    # principal point, these define the viewing rays in terms of the camera
    # frame.
    ray = [x_n; y_n; ones(1,N)]  
                          
    # Rotate to get the viewing rays in the world frame
    ray = C.Rc_w'*ray
    
    # The point of intersection of each ray with the plane will be 
    #   pt = C.P + k*ray    where k is to be determined
    #
    # Noting that the vector planeP -> pt will be perpendicular to planeN.
    # Hence the dot product between these two vectors will be 0.  
    #  dot((( C.P + k*ray ) - planeP) , planeN)  = 0
    # Rearranging this equation allows k to be solved
    pt = zeros(3,N)
    for n = 1:N
        k = (dot(planeP, planeN) - dot(C.P, planeN)) / dot(ray[:,n], planeN)
        pt[:,n] = C.P + k*ray[:,n]
    end

    return pt
end


# Case when the camera is represented by a projection matrix

function imagept2plane(P::Array, xy::Array, planeP::Vector, planeN::Vector)
    
    if size(P) != (3,4) 
        error("Projection matrix must be 3 x 4")
    end

    return pt = imagept2plane(Camera(P), xy, planeP, planeN)
end

#--------------------------------------------------------------------------
"""
camera2projmatrix - Generate a camera projection matrix from a Camera structure
```
Usage:     P = camera2projmatrix(C)

Argument:  C - Camera structure

Returns:   P - 3x4 camera projection matrix that maps homogeneous 3D world 
           coordinates to homogeneous image coordinates.
```
Function takes a camera structure and returns its equivalent projection matrix
ignoring lens distortion parameters etc.

See also: Camera(), projmatrix2camera(), cameraproject()
"""
function camera2projmatrix(C::Camera)
    
    K = [C.fx  C.skew   C.ppx  0
          0    C.fy     C.ppy  0
          0     0         1    0]
    
    T = [ C.Rc_w  -C.Rc_w*C.P
          0  0  0       1     ]
    
    return P = K*T
end

#--------------------------------------------------------------------------
"""
decomposecamera -  Decomposition of a camera projection matrix
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
"""
#=
Reference: Hartley and Zisserman 2nd Ed. pp 155-164

October  2010  Original version
November 2013  Description of rotation matrix R corrected (transposed)
=#

function decomposecamera{T1<:Real}(P::Array{T1,2})

    if size(P) != (3,4) 
        error("Projection matrix must be 3 x 4")
    end
    
    # Convenience variables for the columns of P
    p1 = P[:,1]
    p2 = P[:,2]
    p3 = P[:,3]
    p4 = P[:,4]

    M = [p1 p2 p3];
    m3 = M[3:3,:]'
    
    # Camera centre, analytic solution
    X =  det([p2 p3 p4])
    Y = -det([p1 p3 p4])
    Z =  det([p1 p2 p4])
    T = -det([p1 p2 p3])    
    
    Pc = [X;Y;Z;T]
    Pc = Pc/Pc[4]
    Pc = Pc[1:3]     # Make inhomogeneous
    
    # Pc = null(P,'r'); # numerical way of computing C
    
    # Principal point
    pp = M*m3
    pp = pp/pp[3]
    pp = pp[1:2]   # Make inhomogeneous
    
    # Principal ray pointing out of camera
    pv = det(M)*m3
    pv = pv/norm(pv)
    
    # Perform RQ decomposition of M matrix. Note that rq3 returns K with +ve
    # diagonal elements, as required for the calibration marix.
    (K, Rc_w) = rq3(M)
    
    # Check that R is right handed, if not give warning
    if dot(cross(Rc_w[:,1], Rc_w[:,2]), Rc_w[:,3]) < 0.0
        info("Warning: rotation matrix is left handed")
    end

    return K, Rc_w, Pc, pp, pv
end

#--------------------------------------------------------------------------
"""
rq3 - RQ decomposition of 3x3 matrix
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
"""
#=
Follows algorithm given by Hartley and Zisserman 2nd Ed. A4.1 p 579

October  2010
February 2014  Incorporated modifications suggested by Mathias Rothermel to
               avoid potential division by zero problems
=#

function rq3{T<:Real}(AA::Array{T,2})
    
    if size(AA) != (3,3)
        error("A must be 3x3")
    end
    
    A = copy(AA)  # Avoid side effects when we add eps(1e6) to A[3,3]

    # Find rotation Qx to set A[3,2] to 0
    A[3,3] = A[3,3] + eps(1e6)
    c = -A[3,3]/sqrt(A[3,3]^2+A[3,2]^2)
    s =  A[3,2]/sqrt(A[3,3]^2+A[3,2]^2)
    Qx = [1 0 0; 0 c -s; 0 s c]
    R = A*Qx
    
    # Find rotation Qy to set A[3,1] to 0
    R[3,3] = R[3,3] + eps(1e6)
    c = R[3,3]/sqrt(R[3,3]^2+R[3,1]^2)
    s = R[3,1]/sqrt(R[3,3]^2+R[3,1]^2)
    Qy = [c 0 s; 0 1 0;-s 0 c]
    R = R*Qy
    
    # Find rotation Qz to set A[2,1] to 0
    R[2,2] = R[2,2] + eps(1e6)
    c = -R[2,2]/sqrt(R[2,2]^2+R[2,1]^2)
    s =  R[2,1]/sqrt(R[2,2]^2+R[2,1]^2)    
    Qz = [c -s 0; s c 0; 0 0 1]
    R = R*Qz
    
    Q = Qz'*Qy'*Qx'
    
    # Adjust R and Q so that the diagonal elements of R are +ve
    for n = 1:3
        if R[n,n] < 0
            R[:,n] = -R[:,n]
            Q[n,:] = -Q[n,:]
        end
    end

    return R, Q
end    

#--------------------------------------------------------------------------
"""
makehomogeneous - Appends a scale of 1 to array inhomogeneous coordinates 
```
Usage:  hx = makehomogeneous(x)

Argument:
        x  - an N x npts array of inhomogeneous coordinates.

Returns:
        hx - an (N+1) x npts array of homogeneous coordinates with the
             homogeneous scale set to 1
```
See also: makeinhomogeneous(), hnormalise()
"""
# April 2010

function makehomogeneous{T<:Real}(x::Array{T,2})
    return [x; ones(1,size(x,2))]
end

function makehomogeneous{T<:Real}(x::Vector{T})
    return [x;1]
end

#--------------------------------------------------------------------------
"""
makeinhomogeneous - Converts homogeneous coords to inhomogeneous coordinates 
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
"""
function makeinhomogeneous{T<:Real}(hx::Array{T,2})
    
    x = hnormalise(hx)  # Normalise to scale of one
    x = x[1:end-1,:]    # Extract all but the last row
    return x
end

function makeinhomogeneous{T<:Real}(hx::Vector{T})
    
    x = hnormalise(hx)  # Normalise to scale of one
    x = x[1:end-1]      # Extract all but the last row
    return x
end

#---------------------------------------------------------------------------
"""
hnormalise - Normalises array of homogeneous coordinates to a scale of 1
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
"""
# February 2004

function hnormalise{T<:Real}(x::Array{T,2})
    nx = copy(x)    
    (rows,npts) = size(nx)

    # Find the indices of the points that are not at infinity
    finiteind = find(abs(x[rows,:]) .> eps())

    # Normalise points not at infinity
    for r = 1:rows-1
	nx[r,finiteind] = x[r,finiteind]./x[rows,finiteind]
    end

    nx[rows,finiteind] = 1.0
    return nx
end

function hnormalise{T<:Real}(x::Vector{T})
    nx = copy(x)

    # Only normalise if point is not at infinity
    if abs(nx[end]) > eps()
        nx /= nx[end]
    end

    return nx
end

#---------------------------------------------------------------------------
"""
homography1d - Computes 1D homography
```
Usage:   H = homography1d(x1, x2)

Arguments:
         x1  - 2xN set of homogeneous points
         x2  - 2xN set of homogeneous points such that x1<->x2
Returns:
          H - the 2x2 homography such that x2 = H*x1
```
This code is modelled after the normalised direct linear transformation
algorithm for the 2D homography given by Hartley and Zisserman p92.
"""

# May 2003
# February 2016 ported to Julia

function homography1d{T1<:Real,T2<:Real}(x1::Array{T1,2}, x2::Array{T2,2})
    
    if size(x1) != size(x2)
        error("x1 and x2 must have same dimensions")
    end
    
    (dim, Npts) = size(x1)

    # Attempt to normalise each set of points so that the origin 
    # is at centroid and mean distance from origin is 1.
    (x1n, Tx1) = normalise1dpts(x1)
    (x2n, Tx2) = normalise1dpts(x2)
    
    # Note that it may have not been possible to normalise
    # the points if one was at infinity so the following does not
    # assume that scale parameter w = 1.
    A = zeros(2*Npts,4)
    
    for n = 1:Npts
        X = x1n[:,n]'  
        x = x2n[1,n]
        w = x2n[2,n]
        A[n,:] = [-w*X x*X]
    end
    
    (U,S,V) = svd(A) 
    
    # Extract homography
    H = reshape(V[:,4],2,2)'
    
    # Denormalise
    return H = Tx2\H*Tx1
end

#---------------------------------------------------------------------------
"""
homography2d - Computes 2D homography
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

This code follows the normalised direct linear transformation 
algorithm given by Hartley and Zisserman "Multiple View Geometry in
Computer Vision" p92.
"""
#=
May 2003  - Original version.
Feb 2004  - Single argument allowed for to enable use with RANSAC.
June 2015 - Ported to Julia
=#

function homography2d{T1<:Real,T2<:Real}(x1::Array{T1,2}, x2::Array{T2,2})

    if size(x1) != size(x2)
        error("x1 and x2 must have same dimensions")
    end

    (dim, Npts) = size(x1)

    if dim != 3 
        error("pts must be 3xN")
    end

    if Npts < 4
        error("Number of input points must match and be >= 4")
    end

    # Attempt to normalise each set of points so that the origin 
    # is at centroid and mean distance from origin is sqrt(2).
    (x1n, Tx1) = normalise2dpts(x1)
    (x2n, Tx2) = normalise2dpts(x2)
    
    # Note that it may have not been possible to normalise
    # the points if one was at infinity so the following does not
    # assume that the scale parameter w = 1.
    A = zeros(3*Npts, 9)
    
    O = [0. 0. 0.]
    for n = 1:Npts
	X = x1n[:,n]'
	x = x2n[1,n]
        y = x2n[2,n]
        w = x2n[3,n]
	A[3*n-2,:] = [  O  -w*X  y*X]
	A[3*n-1,:] = [ w*X   O  -x*X]
	A[3*n  ,:] = [-y*X  x*X   O ]
    end
    
    (U,S,V) = svd(A)
    
    H = reshape(V[:,9],3,3)'  # Extract homography
    return H = Tx2\H*Tx1               # and denormalise
end

# Version with a single argument for use with ransac()
function homography2d{T<:Real}(x::Array{T,2})

    (dim, npts) = size(x)

    if dim != 6 
        error("pts must be a 6xN array")
    end

    if npts < 4
        error("Number of input points must match and be >= 4")
    end
    
    return homography2d(x[1:3,:], x[4:6,:])
end

#---------------------------------------------------------------------------
"""
normalise1dpts - Normalises 1D homogeneous points

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
"""

# ? Can this be unified with normalise2dpts ?

function normalise1dpts{T1<:Real}(ptsa::Array{T1,2})

    pts = copy(ptsa) # Copy because we alter ptsa (should be able to fix this)

    if size(pts,1) != 2
        error("pts must be 2xN")
    end

    if any(pts[2,:] .== 0)
        warning("Attempt to normalise a point at infinity")
        return pts, eye(2)
    end
    
    # Ensure homogeneous coords have scale of 1
    pts[1,:] = pts[1,:]./pts[2,:]
    
    c = mean(pts[1,:])        # Centroid.
    newp = pts[1,:] - c       # Shift origin to centroid.
    
    meandist = mean(abs(newp))
    scale = 1/meandist
    
    T = [scale    -scale*c
         0         1      ]
    
    newpts = T*pts

    return newpts, T
end

#---------------------------------------------------------------------------
"""
normalise2dpts - Normalises 2D homogeneous points
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
"""
#=
May 2003      - Original version
February 2004 - Modified to deal with points at infinity.
June 2015     - Ported to Julia
=#

function normalise2dpts{Ty<:Real}(ptsa::Array{Ty,2})

    pts = copy(ptsa) # Copy because we alter ptsa (should be able to fix this)
    newp = zeros(pts)

    if size(pts,1) != 3
        error("pts must be 3xN")
    end
    
    # Find the indices of the points that are not at infinity
    finiteind = find(abs(pts[3,:]) .> eps())
    
    # Should this worning be made?
    if length(finiteind) != size(pts,2)
        warning("Some points are at infinity")
    end
    
    # For the finite points ensure homogeneous coords have scale of 1
    pts[1,finiteind] = pts[1,finiteind]./pts[3,finiteind]
    pts[2,finiteind] = pts[2,finiteind]./pts[3,finiteind]
    pts[3,finiteind] = 1.0
    
    c = mean(pts[1:2,finiteind],2)            # Centroid of finite points
    newp[1,finiteind] = pts[1,finiteind]-c[1] # Shift origin to centroid.
    newp[2,finiteind] = pts[2,finiteind]-c[2]
    
    dist = sqrt(newp[1,finiteind].^2 + newp[2,finiteind].^2)
    meandist = mean(dist) 
    
    scale = sqrt(2.0)/meandist
    
    T = [scale   0.  -scale*c[1]
         0.    scale -scale*c[2]
         0.      0.     1.      ]
    
    newpts = T*pts
    
    return newpts, T 
end

#------------------------------------------------------------------------
"""
skew - Constructs 3x3 skew-symmetric matrix from 3-vector
```
Usage:  s = skew(v)

Argument:  v - 3-vector
Returns:   s - 3x3 skew-symmetric matrix
```
The cross product between two vectors, a x b can be implemented as a matrix
product  skew(a)*b
"""
function skew{T<:Real}(v::Array{T})
    
    @assert length(v) == 3  "Input must be a 3x1 array or vector"
    
    s = [ 0.  -v[3]  v[2]
         v[3]   0.  -v[1]
        -v[2]  v[1]   0. ]
end

#------------------------------------------------------------------------
"""
hcross - Homogeneous cross product, result normalised to s = 1.

Function to form cross product between two points, or lines,
in homogeneous coodinates.  The result is normalised to lie
in the scale = 1 plane.
``` 
Usage: c = hcross(a,b)

Arguments:  a, b  - 3x1 arrays or vectors
Returns:       c  - 3-vector
```
"""
function hcross{T1<:Real, T2<:Real}(a::Array{T1}, b::Array{T2})
    @assert length(a) == 3 && length(b) == 3  "Input must be a 3x1 arrays or vectors"
    c = cross(vec(a), vec(b))
    c = c/c[3]
end

#------------------------------------------------------------------------
"""
fundmatrix - Computes fundamental matrix from 8 or more points

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
"""
#=
Feb 2002  - Original version.
May 2003  - Tidied up and numerically improved.
Feb 2004  - Single argument allowed to enable use with RANSAC.
Mar 2005  - Epipole calculation added, 'economy' SVD used.
Feb 2016  - Ported to Julia
=#

function  fundmatrix{T1<:Real,T2<:Real}(x1i::Array{T1,2}, x2i::Array{T2,2}; epipoles=false)

    if size(x1i) != size(x2i)
        error("x1 and x2 must have same dimensions")
    end

    (dim, npts) = size(x1i)

    if dim != 3 
        error("pts must be a 3xN array of homogeneous coords")
    end

    if npts < 8
        error("Number of input points must match and be >= 8")
    end
    
    # Normalise each set of points so that the origin 
    # is at centroid and mean distance from origin is sqrt(2). 
    # normalise2dpts also ensures the scale parameter is 1.
    (x1, Tx1) = normalise2dpts(x1i)
    (x2, Tx2) = normalise2dpts(x2i)
    
    # Build the constraint matrix
    # The following hopefully works for both v0.4 and v0.5
    A = [x2[1:1,:]'.*x1[1:1,:]'   x2[1:1,:]'.*x1[2:2,:]'  x2[1:1,:]'  x2[2:2,:]'.*x1[1:1,:]'   x2[2:2,:]'.*x1[2:2,:]'  x2[2:2,:]'   x1[1:1,:]'   x1[2:2,:]'  ones(npts,1) ]

    # Build A transposed so that we can format it in code, and then transpose it
#=
    A = [x2[1:1,:].*x1[1:1,:]
         x2[1:1,:].*x1[2:2,:]
         x2[1:1,:] 
         x2[2:2,:].*x1[1:1,:]
         x2[2:2,:].*x1[2:2,:]
         x2[2:2,:]   
         x1[1:1,:]   
         x1[2:2,:]
         ones(1,npts) ]'

println(A)
=#
    (U,D,V) = svd(A,thin=false)

    # Extract fundamental matrix from the column of V corresponding to
    # smallest singular value.
    F = reshape(V[:,9], 3, 3)'
    
    # Enforce constraint that fundamental matrix has rank 2 by performing
    # a svd and then reconstructing with the two largest singular values.
    (U,D,V) = svd(F)
    F = U*diagm([D[1], D[2], 0])*V'
    
    # Denormalise
    F = Tx2'*F*Tx1
    
    if epipoles
        # Solve for epipoles
        (U,D,V) = svd(F)
        e1 = hnormalise(V[:,3])
        e2 = hnormalise(U[:,3])
        
        return F, e1, e2
    else
        return F
    end
end


# Version with a single argument for use with ransac()
function fundmatrix{T<:Real}(x::Array{T,2})

    (dim, npts) = size(x)

    if dim != 6 
        error("pts must be a 6xN array")
    end

    if npts < 8
        error("Number of input points must match and be >= 8")
    end
    
    return fundmatrix(x[1:3,:], x[4:6,:], epipoles=false)
end
    
#--------------------------------------------------------------------------
"""
affinefundmatrix - Computes affine fundamental matrix from 4 or more points

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
"""

# ** Inconsistent in that arguments are inhomogeneous but to use them
# ** with the fundamental matrix you have to make then homogeneous **

# Feb 2005  Original version
# Feb 2016  Ported to Julia

function affinefundmatrix{T1<:Real,T2<:Real}(x1::Array{T1,2}, x2::Array{T2,2})

    if size(x1) != size(x2)
        error("x1 and x2 must have same dimensions")
    end

    (dim, npts) = size(x1)

    if dim != 2
        error("pts must be a 2xN array of inhomogeneous image coords")
    end

    if npts < 4
        error("Number of input points must match and be >= 4")
    end
    
    X = [x2; x1]       # Form vectors of correspondences
    Xmean = mean(X,2)  # Mean 

    deltaX = zeros(size(X))
    for k = 1:4
	deltaX[k,:] = X[k,:] - Xmean[k]
    end
        
    (U,D,V) = svd(deltaX')
    
    # Form  fundamental matrix from the column of V corresponding to
    # smallest singular value.
    v = V[:,4]
    F = [ 0    0    v[1]
	  0    0    v[2]
	 v[3] v[4] -v'*Xmean]
    
    # Solve for epipoles
    (U,D,V) = svd(F)
    e1 = V[:,3]
    e2 = U[:,3]
    
    return F, e1, e2
end

# Version with a single argument for use with ransac()
function affinefundmatrix{T<:Real}(x::Array{T,2})

    (dim, npts) = size(x)

    if dim != 4
        error("pts must be a 4xN array")
    end

    if npts < 4
        error("Number of input points must match and be >= 4")
    end
    
    return affinefundmatrix(x[1:2,:], x[3:4,:])
end
    
#--------------------------------------------------------------------    
"""
fundfromcameras - Fundamental matrix from camera matrices or structures
```
Usage: F = fundfromcameras(P1, P2)
       F = fundfromcameras(C1, C2)

Arguments:  P1, P2 - Two 3x4 camera projection matrices or
            C1, C2 - Two Camera structres

Returns:    F      - Fundamental matrix relating the two camera views.
```
See also: fundmatrix(), affinefundmatrix(), Camera()
"""
# Reference: Hartley and Zisserman p244

# Version for projection matrices
function fundfromcameras{T1<:Real, T2<:Real}(P1::Array{T1,2}, P2::Array{T2,2})
    
    if (size(P1) != (3,4)) || (size(P2) != (3,4))
        error("Camera matrices must be 3x4")
    end
    
    C1 = nullspace(P1)  # Camera centre 1 is the null space of P1
    e2 = P2*C1          # epipole in camera 2
    
    e2x = [  0   -e2[3] e2[2]    # Skew symmetric matrix from e2
            e2[3]    0  -e2[1]
           -e2[2]  e2[1]   0  ]
    
    return F = e2x*P2*pinv(P1)
end

# Version using Camera structures
function fundfromcameras(C1::Camera, C2::Camera)
    return fundfromcameras(camera2projmatrix(C1), camera2projmatrix(C2))
end
    
#----------------------------------------------------------------------
"""
idealimagepts - Ideal image points with no distortion.
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
"""

# August   2015 Adapted from imagept2plane
# February 2016 Refined inversion of distortion process slightly
#               Ported to Julia

function idealimagepts{T<:Real}(C::Camera, xy::Union{Array{T,2}, Vector{T}})
    
    dim = size(xy,1)
    if dim != 2
        error("Image xy data must be a 2 x N array")
    end
    
    # Reverse the projection process as used by cameraproject()

    # Subtract principal point and divide by focal length to get normalised,
    # distorted image coordinates.  Note skew represents the 2D shearing
    # coefficient times fx
    y_d = (xy[2,:] - C.ppy)/C.fy
    x_d = (xy[1,:] - C.ppx - y_d*C.skew)/C.fx
    
    # Compute the inverse of the lens distortion effect by computing the
    # 'forward' direction lens distortion at this point and then subtracting
    # this from the current point.  

    # Radial distortion factor. Here the squared radius is computed from the
    # already distorted coordinates.  The approximation we are making here is to
    # assume that the distortion is locally constant.
    rsqrd = x_d.^2 + y_d.^2
    r_d = 1 + C.k1*rsqrd + C.k2*rsqrd.^2 + C.k3*rsqrd^3
    
    # Tangential distortion component, again computed from the already distorted
    # coords.
    dtx = 2*C.p1*x_d.*y_d         + C.p2*(rsqrd + 2*x_d.^2)
    dty = C.p1*(rsqrd + 2*y_d.^2) + 2*C.p2*x_d.*y_d

    # Subtract the tangential distortion components and divide by the radial
    # distortion factor to get an approximation of the undistorted normalised
    # image coordinates (with no skew)
    x_n = (x_d - dtx)./r_d
    y_n = (y_d - dty)./r_d
    
    # Finally project back to pixel coordinates.
    x_p = C.ppx + x_n*C.fx + y_n*C.skew
    y_p = C.ppy + y_n*C.fy
    
    return xyideal = [x_p
                      y_p]
end

#--------------------------------------------------------------------------
"""
solvestereopt - Homogeneous linear solution of a stereo point
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
"""

function solvestereopt{T1<:Real,T2<:Real}(xy::Array{T1,2}, P::Array{Array{T2,2}})

    (dim,N) = size(xy)
    @assert N == length(P) "Number of points and cameras must match"
    @assert dim == 2  "Image coordinates must be a 2xN array of inhomogeneous coords"
    @assert N >= 2  "Must have at least 2 camera views"
    
    # Build eqn of the form A*pt = 0
    A = zeros(2*N, 4)
    for n = 1:N
	A[2*n-1,:] = xy[1,n]*P[n][3,:] - P[n][1,:]
	A[2*n  ,:] = xy[2,n]*P[n][3,:] - P[n][2,:]
    end
	
    (~,~,v) = svd(A)
    pt = hnormalise(v[:,4])
    
    # Project the point back into the source images to determine the residual
    # error
    xy_reproj = zeros(size(xy))
    for n = 1:N
        xy_reproj[:,n] = cameraproject(P[n], pt[1:3])
    end    

    return pt, xy_reproj
end

# Version where the camera data is an array of Camera structures
function solvestereopt{T<:Real}(xy::Array{T,2}, C::Array{Camera})

    (dim,N) = size(xy)
    @assert N == length(C) "Number of points and cameras must match"
    @assert dim == 2  "Image coordinates must be a 2xN array of inhomogeneous coords"
    @assert N >= 2  "Must have at least 2 camera views"

    # Correct image points from each camera so that they correspond to image
    # points from an ideal camera with no lens distortion.  Also generate the
    # corresponding ideal projection matrices for each Camera structure
    xyideal = zeros(size(xy)) # copy xy so that xy is not changed
    for n = 1:N
        xyideal[:,n] = idealimagepts(C[n], xy[:,n])
        P[n] = camstruct2projmatrix(C[n])
    end

    return solvestereopt(xyideal, P)
end

#--------------------------------------------------------------------------
"""
hline - Plot a 2D line defined in homogeneous coordinates.

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
"""

# Using PyPlot

# Case when 2 homogeneous points are supplied
function hline(p1i::Vector, p2i::Vector, linestyle::ASCIIString="b-")

    p1 = p1i./p1i[3]    # make sure homogeneous points lie in z=1 plane
    p2 = p2i./p2i[3]

    hold(true)    
    plot([p1[1], p2[1]], [p1[2], p2[2]], linestyle);
end

# Case when homogeneous line is supplied
function hline(li::Vector, linestyle::ASCIIString="b-")

    l = li./li[3]   # ensure line in z = 1 plane (not needed?)
    
    if abs(l[1]) > abs(l[2])   # line is more vertical
        p1 = hcross(l, [0, -1, PyPlot.ylim()[1]])
        p2 = hcross(l, [0, -1, PyPlot.ylim()[2]])    
        
    else                       # line more horizontal
        p1 = hcross(l, [-1, 0, PyPlot.xlim()[1]])
        p2 = hcross(l, [-1, 0, PyPlot.xlim()[2]])
    end

    hold(true)
    plot([p1[1], p2[1]], [p1[2], p2[2]], linestyle);
end

#--------------------------------------------------------------------------
"""
plotcamera - Plots graphical representation of camera(s) showing pose
```
Usage: plotcamera(C, l; col=[0,0,1], plotCamPath=false, fig=1)

Arguments:
           C - Camera structure (or array of Camera structures)
           l - The length of the sides of the rectangular cone indicating
               the camera's field of view.
Keyword Arguments:
         col - Optional three element vector specifying the RGB colour to
               use. Defaults to blue.
 plotCamPath - Optional flag true/false to plot line joining camera centre
               positions. If omitted or empty defaults to false.
         fig - Optional figure number to be used. Defaults to 1.
```
The function plots into the current figure a graphical representation of one
or more cameras showing their pose.  This consists of a rectangular cone,
with its vertex at the camera centre, indicating the camera's field of view.
The camera's coordinate X and Y axes are also plotted at the camera centre.

See also: Camera
"""

function plotcamera(Ci, l; col=[0,0,1], plotCamPath=false, fig=1)
    
    if isa(Ci, Array)
        C = Ci
    else
        C = Vector(1)
        C[1] = Ci
    end

    figure(fig)
    hold(true)
    for i = 1:length(C)
        
        if C[i].rows == 0 || C[i].cols == 0
            warning("Camera rows and cols not specified")
            continue
        end
        
        f = C[i].fx  # Use fx as the focal length

        if i > 1 && plotCamPath
            plot3D([C[i-1].P[1], C[i].P[1]],
                   [C[i-1].P(2), C[i].P[2]],
                   [C[i-1].P(3), C[i].P[3]])
        end
        
        # Construct transform from camera coordinates to world coords
        Tw_c = [C[i].Rc_w'  C[i].P
                 0 0 0        1  ]
        
        # Generate the 4 viewing rays that emanate from the principal point and
        # pass through the corners of the image.
        ray = zeros(3,4)
        ray[:,1] = [-C[i].cols/2, -C[i].rows/2, f]
        ray[:,2] = [ C[i].cols/2, -C[i].rows/2, f]
        ray[:,3] = [ C[i].cols/2,  C[i].rows/2, f]
        ray[:,4] = [-C[i].cols/2,  C[i].rows/2, f]

        # Scale rays to distance l from the focal plane and make homogeneous        
        ray = makehomogeneous(ray*l/f)
        ray = Tw_c*ray                 # Transform to world coords
 
        for n = 1:4             # Draw the rays
            plot3D([C[i].P[1], ray[1,n]],
                   [C[i].P[2], ray[2,n]],
                   [C[i].P[3], ray[3,n]],
                   color=col)
        end
        
        # Draw rectangle joining ends of rays
        plot3D([ray[1,1], ray[1,2], ray[1,3], ray[1,4], ray[1,1]],
               [ray[2,1], ray[2,2], ray[2,3], ray[2,4], ray[2,1]],
               [ray[3,1], ray[3,2], ray[3,3], ray[3,4], ray[3,1]],
               color=col)
        
        # Draw and label axes
        X = Tw_c[1:3,1]*l + C[i].P
        Y = Tw_c[1:3,2]*l + C[i].P    
        Z = Tw_c[1:3,3]*l + C[i].P            
        
        plot3D([C[i].P[1], X[1,1]], [C[i].P[2], X[2,1]], [C[i].P[3], X[3,1]],
               color=col)
        plot3D([C[i].P[1], Y[1,1]], [C[i].P[2], Y[2,1]], [C[i].P[3], Y[3,1]],
               color=col)
        #    plot3D([C[i].P[1], Z(1,1)], [C[i].P[2], Z(2,1)], [C[i].P[3], Z(3,1)],...
        #           color=col)
        text3D(X[1], X[2], X[3], "X", color=col)
        text3D(Y[1], Y[2], Y[3], "Y", color=col)
        
    end

end
