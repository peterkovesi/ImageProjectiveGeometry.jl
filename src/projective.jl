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

PK February  2016
   April     2017 Extra functions added, code cleanup.
   September 2018 updated for v0.7/v1.0

---------------------------------------------------------------------=#

export Camera
export cameraproject, camera2projmatrix, imagept2plane, decomposecamera, rq3

export mapimage2plane!, imagecorners2plane
export imagept2ray, imagept2ray!

export makehomogeneous, makeinhomogeneous, hnormalise, hnormalise!
export solveaffine, imgtrans
export homography1d, homography2d, normalise1dpts, normalise2dpts
export skew, hcross
export fundmatrix, affinefundmatrix, fundfromcameras

export stereorectify, stereorectifytransforms
export transformedimagebounds

export idealimagepts, solvestereopt, undistortimage
export hline, plotcamera

using LinearAlgebra, Statistics, Printf, Interpolations, PyPlot

#=
imTrans.m
imTransD.m
equalAngleConstraint.m
knownAngleConstraint.m
lengthRatioConstraint

=#

#------------------------------------------------------------------------

# To consider in the future:
# Define an AbstractCamera type so that subtypes of Camera,
# FishEyeCamera, RollingShutterCamera etc etc. can be defined

"""
Camera - Structure defining parameters of a camera following the Brown
         distortion model
```
      fx::Float64      # Focal length.
      fy::Float64
     ppx::Float64      # Principal point.
     ppy::Float64
      k1::Float64      # Radial lens distortion parameters.
      k2::Float64
      k3::Float64
      p1::Float64      # Tangential lens distortion parameters.
      p2::Float64
    skew::Float64
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
                    # all other (distortion) parameters are set to 0.
```

"""
mutable struct Camera
    fx::Float64           # Focal length.
    fy::Float64
    ppx::Float64          # Principal point.
    ppy::Float64
    k1::Float64           # Radial lens distortion.
    k2::Float64
    k3::Float64
    p1::Float64           # Tangential lens distortion.
    p2::Float64
    skew::Float64
    rows::Int       # Optional, providing rows and columns allows detection
    cols::Int       # of points being projected out of image bounds.
    P::Vector{Float64}    # Camera position in world coordinates.
    Rc_w::Matrix{Float64}  # Rotation matrix defining world orientation with
                          # respect to the camera frame.
end

# Keyword-value constructor
function Camera(;fx=1.0, fy=1.0, ppx=0.0, ppy=0.0,
                k1=0.0, k2=0.0, k3=0.0, p1=0.0, p2=0.0, skew=0.0,
                rows=0, cols=0,
                P=[0.0, 0.0, 0.0], Rc_w=I(3))
    if size(P) != (3,)
        error("Camera position must be a 3x1 array")
    end

    if size(Rc_w) != (3,3)
        error("Camera orientation must be a 3x3 rotation matrix")
    end

    return Camera(float(fx), float(fy), float(ppx), float(ppy),
                  float(k1), float(k2), float(k3), float(p1), float(p2),
                  float(skew), rows, cols, float(P), float(Rc_w))
end

# Contructor that takes a projection matrix
function Camera(P::Array{T}) where T <: AbstractFloat
    if size(P) != (3,4)
        error("Projection matrix must be 3x4")
    end

    (K, Rc_w, Pc, pp, pv) = decomposecamera(P)
    Camera(fx=K[1,1], fy=K[2,2], ppx=pp[1], ppy=pp[2], skew=K[1,2], P=Pc, Rc_w=Rc_w)
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
function cameraproject(P::Array, pt::Array) # Projection matrix version

    if size(P) != (3,4)
        error("Projection matrix must be 3x4")
    end

    if size(pt, 1) != 3
        error("Points must be in a 3xN array")
    end

    nPts = size(pt,2)
    xy = zeros(2,nPts)

    for i in 1:nPts
        s = P[3,1]*pt[1,i] + P[3,2]*pt[2,i] + P[3,3]*pt[3,i] + P[3,4]
        xy[1,i] = (P[1,1]*pt[1,i] + P[1,2]*pt[2,i] + P[1,3]*pt[3,i] + P[1,4])/s
        xy[2,i] = (P[2,1]*pt[1,i] + P[2,2]*pt[2,i] + P[2,3]*pt[3,i] + P[2,4])/s
    end

    return xy
end

#--------------------------------------------------------------------------
"""
cameraproject - Projects 3D points into camera image
```
Usage:             xy = cameraproject(C, pt, computevisibility = false)
        (xy, visible) = cameraproject(C, pt, computevisibility = true)

Arguments:
             C - Camera structure.
            pt - 3xN matrix of 3D points to project into the image.

Keyword Argument:
 computevisibility - If set to true point visibility is returned
                     The default value is false.

Returns:
       xy      - 2xN matrix of projected image positions
       visible - Boolean array indicating whether the point is
                 within the field of view.  This is only evaluated if
                 'computevisibility is true and the camera structure has
                 non zero values for its 'rows' and 'cols' fields.

```
See also: Camera(), camstruct2projmatrix()
"""
function cameraproject(C::Camera, pta::Array; computevisibility = false)
    # Camera structure version

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
        pt[n,:] = pt[n,:] .- C.P[n]
    end

    # Then rotate to camera frame using Rc_w
    pt = C.Rc_w*pt

    # Follow classical projection process
    # Generate normalized coords
    x_n = pt[1:1,:]./pt[3:3,:]
    y_n = pt[2:2,:]./pt[3:3,:]
    rsqrd = x_n.^2 .+ y_n.^2  # radius squared from centre

    # Radial distortion factor
    r_d = 1.0 .+ C.k1*rsqrd .+ C.k2*rsqrd.^2 .+ C.k3*rsqrd.^3

    # Tangential distortion component
    dtx = 2*C.p1*x_n.*y_n          .+ C.p2*(rsqrd .+ 2*x_n.^2)
    dty = C.p1*(rsqrd .+ 2*y_n.^2) .+ 2*C.p2*x_n.*y_n

    # Obtain lens distorted coordinates
    x_d = r_d.*x_n .+ dtx
    y_d = r_d.*y_n .+ dty

    # Finally project to pixel coordinates
    # Note skew represents the 2D shearing coefficient times fx
    x_p = C.ppx .+ x_d*C.fx .+ y_d*C.skew
    y_p = C.ppy .+ y_d*C.fy

    xy = [x_p
	  y_p]

    # If C.rows and C.cols not empty determine points that are within image bounds
    if computevisibility
        if C.rows !=0 && C.cols != 0
           visible = (x_p .>= 1.0) .& (x_p .<= convert(Float64,C.cols)) .&
                     (y_p .>= 1.0) .& (y_p .<= convert(Float64,C.rows))
        else
            @warn("Point visibility requested but Camera structure has no image size data")
            visible = Array{Bool}(undef, 0)
        end

        return xy, visible
    else
        return xy
    end

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
function imagept2plane(C::Camera, xy::Array, planeP::Vector, planeN::Vector)

    (dim, N) = (size(xy,1), size(xy,2))  # xy might be a Vector
    if dim != 2
        error("Image xy data must be a 2 x N array")
    end

    # Reverse the projection process as used by cameraproject()

    # Subtract principal point and divide by focal length to get normalised,
    # distorted image coordinates.  Note skew represents the 2D shearing
    # coefficient times fx
    y_d = (xy[2:2,:] .- C.ppy)/C.fy
    x_d = (xy[1:1,:] .- C.ppx .- y_d*C.skew)/C.fx

    # Radial distortion factor.  Here the squared radius is computed from the
    # already distorted coordinates.  The approximation we are making here is to
    # assume that the distortion is locally constant.
    rsqrd = x_d.^2 .+ y_d.^2
    r_d = 1 .+ C.k1*rsqrd .+ C.k2*rsqrd.^2 .+ C.k3*rsqrd.^3

    # Tangential distortion component, again computed from the already distorted
    # coords.
    dtx = 2*C.p1*x_d.*y_d          .+ C.p2*(rsqrd .+ 2*x_d.^2)
    dty = C.p1*(rsqrd .+ 2*y_d.^2) .+ 2*C.p2*x_d.*y_d

    # Subtract the tangential distortion components and divide by the radial
    # distortion factor to get an approximation of the undistorted normalised
    # image coordinates (with no skew)

    # inverse of the radial distortion is calculated by iteration
    # (see e.g., Drap&Lefèvre 2015, doi: 10.3390/s16060807)
    # non-convergence of the mothod happens if an image point has no corresponding
    # point on the plane, due to distortion. This is not caught currently

    x_n = copy(x_d)
    y_n = copy(y_d)
    for i=1:N
        trd=0.0
        rsq=0.0
        for _=1:100
            rsq=x_n[i]^2+y_n[i]^2
            trd = 1 + C.k1*rsq + C.k2*rsq^2 + C.k3*rsq^3
            abs(x_n[i] - x_d[i]/trd) <1e-8 && abs(y_n[i] - y_d[i]/trd) <1e-8 && break #convergence check
            x_n[i] = x_d[i] / trd
            y_n[i] = y_d[i] / trd
        end
        x_n[i] = (x_d[i] - dtx[i]) / trd
        y_n[i] = (y_d[i] - dty[i]) / trd
    end

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
        pt[:,n] .= C.P .+ k*ray[:,n]
    end

    return pt
end


# Case when the camera is represented by a projection matrix

function imagept2plane(Proj::Array, cameraP::Vector, xy::Array, planeP::Vector, planeN::Vector)
    (dim, N) = (size(xy,1), size(xy,2))  # xy might be a Vector
    if dim != 2
        error("Image xy data must be a 2 x N array")
    end
    res = zeros(3, N)
    for i in 1:N
        res[:,N] = imagept2plane(Proj, cameraP, xy[:,i], planeP, planeN)
    end
    res
end

function imagept2plane(Proj::Array, cameraP::Vector, xy::Vector, planeP::Vector, planeN::Vector)

    if size(Proj) != (3,4)
        error("Projection matrix must be 3 x 4")
    end

    # Solve for point on the viewing ray closest to the origin via
    # pseudoinverse of Proj.
    pt = Proj\[xy[1], xy[2], 1]
    pt /= pt[4]
    ray = pt[1:3] .- cameraP   # Ray from camera centre to point

    # The point of intersection of the ray with the plane will be
    #   pt = cameraP + k*ray    where k is to be determined
    #
    # Noting that the vector planeP -> pt will be perpendicular to planeN.
    # Hence the dot product between these two vectors will be 0.
    #  dot((( cameraP + k*ray ) - planeP) , planeN)  = 0
    # Rearranging this equation allows k to be solved
    k = (dot(planeP, planeN) - dot(cameraP, planeN)) / dot(ray, planeN)

    return  cameraP .+ k*ray
end



#-------------------------------------------------------------------------------------
"""
mapimage2plane! - Projects an image onto a plane in 3D

```
Usage:  mapimage2plane!(mappedimg, rect3Dpts, img, P)

Arguments:
       mappedimg - Buffer for storing the resulting projected image. ::Array{Float64,2}
          rect3D - 3x4 array specifying the 3D coordinates of the four
                   corners of the image rectangle. Points are ordered
                   clockwise from the top-left corner.
             img - Image to be projected. ::Array{Float64,2}
               P - 3x4 projection matrix for img.

```

The four 3D points specifying the corners of the rectangle in the
plane are projected into the image to obtain their corresponding image
coordinates.  A homography is computed between these points and the
corner coordinates of the mappedimage buffer.  This is then applied to
img to project it into mappedimg.

"""
function mapimage2plane!(mappedimg::Array{Float64,2}, rect3Dpts::Array{Float64,2},
                                         img::Array{Float64,2}, P::Array{Float64,2})

    (rows,cols) = size(mappedimg)
    mappedxy = float( [1    cols  cols   1
                       1     1    rows  rows
                       1     1     1     1  ] )

    # Project the 3D corner locations of rectangle to the image
    if size(rect3Dpts,1) == 3
        imgxy = P*makehomogeneous(rect3Dpts)
    elseif size(rect3Dpts,1) == 4  # Assume homogeneous coords supplied
        imgxy = P*rect3Dpts
    else
        error("rect3Dpts must be 3x4 inhomogeneous, or 4x4 homogeneous points")
    end

    H = homography2d(imgxy, mappedxy)  # Generate homography
    imtrans!(mappedimg, img, H)        # and apply homography to img.

    return nothing
end


#-------------------------------------------------------------------------------------
"""
imagecorners2plane - Get the positions of image corners projected onto a plane.

```
Usage: pt = imagecorners2plane(C, planeP, planeN)

Arguments:
        C - Camera structure or 3x4 projection matrix.
   planeP - 3-vector defining a point on the plane.
   planeN - 3-vector defining the normal of the plane.

Returns:
        pt - 3x4 array of 3D points on the plane that the image corners
             project to. Points are orderer clockwise from the top-left
             corner of the image.

```

See also: imagept2plane()
"""
function imagecorners2plane(C::Camera, planeP::Vector, planeN::Vector)

    xy = [0  C.cols  C.cols   0
          0    0     C.rows  C.rows]

    return imagept2plane(C, xy, planeP, planeN)
end


#--------------------------------------------------------------------------
"""
imagept2ray! - Compute viewing ray corresponding to an image point.

```
Usage:  ray = imagept2ray!(ray, C, x, y)

Arguments:
       ray - 3-vector to store the result.
         C - Camera structure, Alternatively
             C can be a 3x4 camera projection matrix.
      x, y - Image point  (col,row)

Returns:
       ray - 3-vector corresponding to the image point
```

Lens distortion is handled by using the standard lens distortion
parameters assuming locally that the distortion is constant, computing
the forward distortion and then subtracting the distortion.

See also imagept2ray(), Camera(), cameraproject()
"""
function imagept2ray!(ray, C::Camera, x, y)

    # Reverse the projection process as used by cameraproject()

    # Subtract principal point and divide by focal length to get normalised,
    # distorted image coordinates.  Note skew represents the 2D shearing
    # coefficient times fx
    y_d = (y .- C.ppy)/C.fy
    x_d = (x .- C.ppx .- y_d*C.skew)/C.fx

    # Radial distortion factor.  Here the squared radius is computed from the
    # already distorted coordinates.  The approximation we are making here is to
    # assume that the distortion is locally constant.
    rsqrd = x_d.^2 .+ y_d.^2
    r_d = 1 .+ C.k1*rsqrd .+ C.k2*rsqrd.^2 .+ C.k3*rsqrd.^3

    # Tangential distortion component, again computed from the already distorted
    # coords.
    dtx = 2*C.p1*x_d.*y_d          .+ C.p2*(rsqrd .+ 2*x_d.^2)
    dty = C.p1*(rsqrd .+ 2*y_d.^2) .+ 2*C.p2*x_d.*y_d

    # Subtract the tangential distortion components and divide by the radial
    # distortion factor to get an approximation of the undistorted normalised
    # image coordinates (with no skew)
    x_n = (x_d .- dtx)./r_d
    y_n = (y_d .- dty)./r_d

    # [x_n, y_n, 1] define a set of points at the normalised distance
    # of z = 1 from the principal point, these define the viewing rays
    # in terms of the camera frame.  Rotate to get the viewing rays in
    # the world frame
    # ray = C.Rc_w'*ray
    mul!(ray, C.Rc_w', [x_n, y_n, 1])

    normalize!(ray)
    return ray
end

"""
imagept2ray - Compute viewing ray corresponding to an image point.

```
Usage:  ray = imagept2ray(C, x, y)

Arguments:
         C - Camera structure, Alternatively
             C can be a 3x4 camera projection matrix.
      x, y - Image point  (col,row)

Returns:
       ray - 3-vector corresponding to the image point
```

Lens distortion is handled by using the standard lens distortion
parameters assuming locally that the distortion is constant, computing
the forward distortion and then subtracting the distortion.

See also imagept2ray!(), Camera(), cameraproject()
"""
function imagept2ray(C, x, y)
    ray = zeros(3)
    imagept2ray!(ray, C, x, y)
    return ray
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
          0     0        1     0]

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
function decomposecamera(P::Array{T1,2}) where T1 <: AbstractFloat
    # Reference: Hartley and Zisserman 2nd Ed. pp 155-164

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

    Pc = [X,Y,Z,T]
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
        @warn("Note rotation matrix is left handed")
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
function rq3(A::Array{T,2}) where T <: Real
    #=
    Follows algorithm given by Hartley and Zisserman 2nd Ed. A4.1 p 579

    October  2010
    February 2014  Incorporated modifications suggested by Mathias Rothermel to
    avoid potential division by zero problems
    =#

    if size(A) != (3,3)
        error("A must be 3x3")
    end

    # Find rotation Qx required to set A[3,2] to 0
    A33 = A[3,3] + eps(1e6)  # Avoid side effects when we add eps(1e6)
    c = -A33/sqrt(A33^2+A[3,2]^2)
    s =  A[3,2]/sqrt(A33^2+A[3,2]^2)
    Qx = [1 0 0; 0 c -s; 0 s c]
    R = A*Qx

    # Find rotation Qy required to set A[3,1] to 0
    R[3,3] = R[3,3] + eps(1e6)
    c = R[3,3]/sqrt(R[3,3]^2+R[3,1]^2)
    s = R[3,1]/sqrt(R[3,3]^2+R[3,1]^2)
    Qy = [c 0 s; 0 1 0;-s 0 c]
    R = R*Qy

    # Find rotation Qz required to set A[2,1] to 0
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
function makehomogeneous(x::Array{T,2}) where T <: Real
    return [x; ones(1,size(x,2))]
end

function makehomogeneous(x::Vector{T}) where T <: Real
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
function makeinhomogeneous(hx::Array{T,2}) where T <: Real

    x = hnormalise(hx)  # Normalise to scale of one
    x = x[1:end-1,:]    # Extract all but the last row
    return x
end

function makeinhomogeneous(hx::Vector{T}) where T <: Real

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
             that the scale values nx[end,:] are all 1.
```
Note that any homogeneous coordinates at infinity (having a scale value of
0) are left unchanged.

See also: hnormalise!()
"""
function hnormalise(x::Union{Array{T,2},Vector{T}}) where T <: Real
    hnormalise!(copy(x))
end

#---------------------------------------------------------------------------
"""
hnormalise! - In-place normalisation of homogeneous coordinates to a scale of 1
```
Usage:   hnormalise!(x)

Argument:
        x - An Nxnpts array of homogeneous coordinates.

Return value and side effect:
        x - The homogeneous coordinates in x are rescaled so
            that the scale values x[end,:] are all 1.
```
Note that any homogeneous coordinates at infinity (having a scale value of
0) are left unchanged.

See also: hnormalise()

"""
function hnormalise!(x::Union{Array{T,2},Vector{T}}) where T <: Real

    (rows, npts) = (size(x,1), size(x,2))  # This handles 2D arrays and vectors

    for n = 1:npts
        if abs(x[rows, n]) > eps(T)   # point not at infinity
            for r = 1:rows-1
                x[r, n] /= x[rows, n]
            end
            x[rows, n] = 1
        end
    end

    return x
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
function homography1d(x1::Array{T1,2}, x2::Array{T2,2}) where {T1 <: Real, T2 <: Real}

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
function homography2d(x1::Array{T1,2}, x2::Array{T2,2}) where {T1 <: Real, T2 <: Real}
    #=
    May 2003  - Original version.
    Feb 2004  - Single argument allowed for to enable use with RANSAC.
    June 2015 - Ported to Julia
    =#

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
function homography2d(x::Array{T,2}) where T <: Real

    (dim, npts) = size(x)

    if dim != 6
        error("pts must be a 6xN array")
    end

    if npts < 4
        error("Number of input points must match and be >= 4")
    end

    return homography2d(x[1:3,:], x[4:6,:])
end

#-------------------------------------------------------------------------------
"""
solveaffine - Solve affine transformation between two sets of 2D points.

```
Usage: A = solveaffine(xy1, xy2)

Arguments:
   xy1, xy2 - 2xN arrays of corresponding 2D points

Returns:
     A - 3x3 affine transformation matrix such that xy2 = A*xy1
         (assuming xy1 and xy2 are in homogeneous coords)

    [ x2        [ a  b  c    [ x1
      y2    =     d  e  f      y1
       1 ]        0  0  1 ]     1 ]

```
"""
function solveaffine(xy1, xy2)

    if size(xy1) != size(xy2)
        error("inputs point sets must be equal in size")
    end

    (dim, npts) = size(xy1)
    X = zeros(2*npts,6)

    for n = 1:npts
        X[2*n-1,:] = [xy1[1,n] xy1[2,n]  1    0        0     0]
        X[2*n  ,:] = [   0        0      0 xy1[1,n] xy1[2,n] 1]
    end

    a = X \ xy2[:]

    A = [ a[1] a[2] a[3]
          a[4] a[5] a[6]
           0    0    1  ]

    return A
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
  pts -  2xN array of 1D homogeneous coordinates

Returns:
  newpts -  2xN array of transformed 1D homogeneous coordinates
  T      -  The 2x2 transformation matrix, newpts = T*pts
```
Note that if one of the points is at infinity no normalisation
is possible.  In this case a warning is printed and pts is
returned as newpts and T is the identity matrix.

See also: normalise2dpts()
"""
function normalise1dpts(ptsa::Array{T1,2}) where T1 <: Real
    # ? Can this be unified with normalise2dpts ?

    pts = copy(ptsa)

    if size(pts,1) != 2
        error("pts must be 2xN")
    end

    if any(pts[2,:] .== 0)
        @warn("Attempt to normalise a point at infinity")
        return pts, I(2)
    end

    # Ensure homogeneous coords have scale of 1
    hnormalise!(pts)

    # Get centroid and mean distance from centroid
    c = mean(view(pts, 1,:))
    meandist = mean(abs.(pts[1,:] .- c))
    scale = 1/meandist

    T = [scale    -scale*c
         0         1      ]

    pts = T*pts

    return pts, T
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
function normalise2dpts(ptsa::Array{T1,2}) where T1 <: Real

    pts = copy(ptsa) # Copy because we alter ptsa (should be able to fix this)
    newp = zeros(size(pts))

    if size(pts,1) != 3
        error("pts must be 3xN")
    end

    # Find the indices of the points that are not at infinity
    finiteind = findall(abs.(pts[3,:]) .> eps())

    # Should this warning be made?
    if length(finiteind) != size(pts,2)
        @warn("Some points are at infinity")
    end

    # For the finite points ensure homogeneous coords have scale of 1
    pts[1,finiteind] = pts[1,finiteind]./pts[3,finiteind]
    pts[2,finiteind] = pts[2,finiteind]./pts[3,finiteind]
    pts[3,finiteind] .= 1.0

    c = mean(pts[1:2,finiteind],dims=2)          # Centroid of finite points
    newp[1,finiteind] = pts[1,finiteind] .- c[1] # Shift origin to centroid.
    newp[2,finiteind] = pts[2,finiteind] .- c[2]

    dist = sqrt.(newp[1,finiteind].^2 .+ newp[2,finiteind].^2)
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
function skew(v::Array{T}) where T <: Real

    @assert length(v) == 3  "Input must be a 3x1 array or vector"

    s = [ zero(T)  -v[3]    v[2]
           v[3]   zero(T)  -v[1]
          -v[2]     v[1]   zero(T) ]
end

#------------------------------------------------------------------------
"""
hcross - Homogeneous cross product, result normalised to s = 1.

Function to form cross product between two points, or lines,
in homogeneous coodinates.  The result is normalised to lie
in the scale = 1 plane.
```
Usage: c = hcross(a,b)

Arguments:  a, b  - 3-vectors
Returns:       c  - 3-vector
```
"""
function hcross(a::Array{T1}, b::Array{T2}) where {T1 <: Real, T2 <: Real}
    @assert length(a) == 3 && length(b) == 3  "Input must be a 3-vectors"
    c = cross(a, b)
    c /= c[3]
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
function  fundmatrix(x1i::Array{T1,2}, x2i::Array{T2,2}; epipoles=false) where {T1<:Real,T2<:Real}

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
    # The following cumbersome construction hopefully works for both v0.4 and v0.5 ** check for 0.7 **
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
    (U,D,V) = svd(A, full=true)

    # Extract fundamental matrix from the column of V corresponding to
    # smallest singular value.
    F = reshape(V[:,9], 3, 3)'

    # Enforce constraint that fundamental matrix has rank 2 by performing
    # a svd and then reconstructing with the two largest singular values.
    (U,D,V) = svd(F)
    F = U*diagm(0 => [D[1], D[2], 0])*V'

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
function fundmatrix(x::Array{T,2}) where T <: Real

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
function affinefundmatrix(x1::Array{T1,2}, x2::Array{T2,2}) where {T1<:Real,T2<:Real}
    # ** Inconsistent in that arguments are inhomogeneous but to use them
    # ** with the fundamental matrix you have to make then homogeneous **

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
    Xmean = mean(X, dims=2)  # Mean

    deltaX = zeros(size(X))
    for k = 1:4
	deltaX[k,:] = X[k,:] .- Xmean[k]
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
function affinefundmatrix(x::Array{T,2}) where T <: Real

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
function fundfromcameras(P1::Array{T1,2}, P2::Array{T2,2}) where {T1<:Real, T2<:Real}
    # Reference: Hartley and Zisserman p244
    # Version for projection matrices

    if (size(P1) != (3,4)) || (size(P2) != (3,4))
        error("Camera matrices must be 3x4")
    end

    C1 = nullspace(P1)  # Camera centre 1 is the null space of P1
    e2 = P2*C1          # epipole in camera 2

    e2x = [   0   -e2[3] e2[2]    # Skew symmetric matrix from e2
            e2[3]    0  -e2[1]
           -e2[2]  e2[1]   0  ]

    return F = e2x*P2*pinv(P1)
end

# Version using Camera structures
function fundfromcameras(C1::Camera, C2::Camera)
    return fundfromcameras(camera2projmatrix(C1), camera2projmatrix(C2))
end


#------------------------------------------------------------------------------
"""
stereorectify - Rectify a stereo pair of images.

```
Usage: (P1r, img1r, H1, P2r, img2r, H2, dmin, dmax, dmed) = ...
                    stereorectify(P1, img1, P2, img2, xy1=zeros(0), xy2=zeros(0),
                     scale=1.0, disparitytruncation=0, diagnostics=false)

Arguments:
          P1, P2 - Projection matrices of the images to be rectified.
      img1, img2 - The two images, may be multi-channel (assumed same size).
         xy1,xy2 - Optional: Two sets of corresponding homogeneous image
                   coordinates in the images [3 x N] in size.  If supplied
                   these are used to construct a rectification that seeks
                   to minimise the relative displacement between the
                   rectified image points.
Keyword arguments:
                scale - The desired scaling of the image, defaults to 1.0
  disparitytruncation - The percentage to be clipped off the ends of the histogram
                        of image coordiante disparities for the determination of
                        the disparity range.
          diagnostics - If true diagnostic plots are generated and disparity range
                        printed on screen

Returns:
        P1r, P2r - The projection matrices of the rectified images.
    img1r, img2r - The rectified images.
          H1, H2 - Homographies that were applied to img1 and img2 to obtain
                   img1r, img2r.
      dmin, dmax - The range of disparities between the images derived from
                   the rectified image coordinates (if supplied).
            dmed - Median of disparities.
```

Note the value of supplying xy1 and xy2 is that the local spatial
arrangement of pixels at any matching point between images is likely
to be more similar and hence image descriptors are likely to work
better.  Also, the disparity ranges between the images will be kept as
uniform as possible.

The range of disparities is derived from a histogram of the
disparities between the supplied image coordinates with a specified
truncation applied to the ends of the histograms to exclude outliers.

See also: stereorectifytransforms()
"""
function stereorectify(P1, img1, P2, img2, xy1=zeros(0), xy2=zeros(0);
                       scale::Real=1.0, disparitytruncation::Real=0.0, diagnostics::Bool=false)
    # Reference:  Hartley and Zisserman 2nd Ed. Section 11.12, p302

    # Get the transformations required to perform stereorectification
    (P1r, H1, P2r, H2, dmin, dmax, dmed) = stereorectifytransforms(P1, img1, P2, img2,
                                   xy1, xy2, scale=scale, disparitytruncation=disparitytruncation)

    # Apply homographies to obtaoin rectified images.
    img1r = imtrans(img1, H1)
    img2r = imtrans(img2, H2)

    if diagnostics
        figure(21); clf()
        imshow(img1r, cmap=ColorMap("gray"))
        figure(22); clf()
        imshow(img2r, cmap=ColorMap("gray"))

        # Draw some horizontal lines on the rectified images to allow diagnostic
        # checks
        rows = max(size(img1r,1), size(img2r,1))
        cols = max(size(img1r,2), size(img2r,2))
        for fig = [21, 22]
            figure(fig)
#            hold(true)
            for y = 250:250:rows
                plot([0, cols], [y, y], "r-")
            end
#            hold(false)
        end

        if !isempty(xy1) && !isempty(xy2)  # We have affine correction data
            # Apply rectification transforms to the image points to determine the range
            # of disparities in the x direction and to check that there are no serious
            # shifts in the y direction.
            H1xy1 = hnormalise(H1*xy1)
            H2xy2 = hnormalise(H2*xy2)


            # Generate histogram of y disparities, should be almost zero
            ydisparity = vec(H2xy2[2,:] - H1xy1[2,:])
            plothistogram(ydisparity, linspace(minimum(ydisparity), maximum(ydisparity), 50),
                          fig = 23, title = "Y Disparities of rectified measurement points")
            @printf("Y disparity: mean %f   std dev %f\n", mean(ydisparity), std(ydisparity))

            @printf("Min, max and median disparity of 2nd image wrt to first is %d,  %d,  %d\n",
                    dmin, dmax, dmed)

            # Plot rectified tiepoint image positions on rectified images
            #            figure(21), hold on
            #            plot(H1xy1[1,:], H1xy1[2,:], "r+"); hold off
            #            figure(22), hold on
            #            plot(H2xy2[1,:], H2xy2[2,:], "r+"); hold off
        end
    end

    return (P1r, img1r, H1, P2r, img2r, H2, dmin, dmax, dmed)
end

#-----------------------------------------------------------------------------

"""
stereorectifytransforms

Compute homographies that transform an image pair into a stereorectified pair.

```
Usage:  (P1r, H1, P2r, H2, dmin, dmax, dmed) = stereorectifytransforms(P1, img1, P2, img2,
                                   xy1=zeros(0), xy2=zeros(0); scale=1.0, disparitytruncation=0)

Arguments:
          P1, P2 - Projection matrices of the images to be rectified.
      img1, img2 - The two images, may be multi-channel (assumed same size).
         xy1,xy2 - Optional: Two sets of corresponding homogeneous image
                   coordinates in the images [3 x N] in size.  If supplied
                   these are used to construct a rectification that seeks
                   to minimise the relative displacement between the
                   rectified image points.
Keyword argument:
                scale - The desired scaling of the image, defaults to 1.0
  disparitytruncation - The percentage to be clipped off the ends of the histogram
                        of image coordiante disparities for the determination of
                        the disparity range.

Returns:
        P1r, P2r - The projection matrices of the rectified images.
          H1, H2 - Homographies that should be applied to img1 and img2 to
                   obtain a stereorectifed pair
      dmin, dmax - The range of disparities between the images derived from
                   the rectified image coordinates (if supplied).
            dmed - Median of disparities.
```

Note the value of supplying xy1 and xy2 is that the local spatial
arrangement of pixels at any matching point between images is likely
to be more similar and hence image descriptors are likely to work
better.  Also, the disparity ranges between the images will be kept as
uniform as possible.

The range of disparities is derived from a histogram of the
disparities between the supplied image coordinates with a specified
truncation applied to the ends of the histograms to exclude outliers.

See also: stereorectify()
"""
function stereorectifytransforms(P1, img1, P2, img2, xy1=zeros(0), xy2=zeros(0);
                                 scale::Real=1.0, disparitytruncation::Real=0.0)
    # Reference:  Hartley and Zisserman 2nd Ed. Section 11.12, p302

    AFFINECORRECT = !isempty(xy1) && !isempty(xy2)

    (rows1, cols1, chan) = (size(img1, 1), size(img1,2), size(img1,3))

    # Get Fundamental matrix and epipoles
    # F *e1 = 0
    # F'*e2 = 0
    F = fundfromcameras(P1, P2)
    (U,D,V) = svd(F)
    e1 = hnormalise(V[:,3])
    e2 = hnormalise(U[:,3])

    # First build the transform that is applied to image1

    # Tranform that takes centre of image to origin
    T = [1.0   0.0 -cols1/2
         0.0   1.0 -rows1/2
         0.0   0.0   1.0   ]

    # Rotation that takes epipole to [f,0,1]
    Te1 = T*e1                    # Location of epipole after translation.
    Te1 = Te1[1:2]/norm(Te1[1:2]) # Unit direction vetor.

    R = [ Te1[1]  Te1[2]  0.0
         -Te1[2]  Te1[1]  0.0
           0.0     0.0    1.0]

    # Location of epipole after applying T and R.  Position will be of the
    # form [f,0,1]
    RTe1 = hnormalise(R*T*e1)
    f = RTe1[1]

    # Tranform that maps rotated epipole [f,0,1] to infinity at [f,0,0]
    G = [ 1.0   0.0   0.0
          0.0   1.0   0.0
         -1/f   0.0   1.0]

    H1 = G*R*T     # Transformation to be applied to image 1.

    # We now need to find a matching transformation H2 that maps epipolar lines
    # in image 2 so that they match up with the epipolar lines in image 1.
    # Additionally we want H2 such that the displacement of corresponding
    # image points in rectified images 1 and 2 is minimised.

    # H2 is of the form Ha*Ho .
    # Where: Ho = H1*M ,  M is the transfer plane such that F = skew(e1)*M
    # and Ha is an affine matrix of the form [a b c; 0 1 0; 0 0 1].  This
    # attempts to scale, shear and offset the image in the x direction so
    # that corresponding image points have minimal displacement.
    M = P1*pinv(P2)
    Ho = H1*M

    if AFFINECORRECT     # Solve for a, b, c in Ha such that corresponding
                         # point displacements are minimised.
        H1xy1 = hnormalise(H1*xy1)   # Rectified points in image 1
        Hoxy2 = hnormalise(Ho*xy2)   # Image 2 points rectified by Ho
        abc = Hoxy2'\H1xy1[1:1,:]'
        H2 = [abc'; 0.0 1.0 0.0; 0.0 0.0 1.0]*Ho
    else
        H2 = Ho
    end

    # Scale the transforms to required value
    S = [scale 0 0; 0 scale 0; 0 0 1]
    H1 = S * H1
    H2 = S * H2

    if AFFINECORRECT
        # Apply rectification transforms to the image points to
        # determine the range of disparities in the x direction and to
        # check that there are no serious shifts in the y direction.
        H1xy1 = hnormalise(H1*xy1)
        H2xy2 = hnormalise(H2*xy2)

        # Get min and max row and column coords of the transformed
        # image points and cut the images down to the minimum size
        # that covers the shared region
        minx1 = minimum(H1xy1[1,:])
        maxx1 = maximum(H1xy1[1,:])
        miny1 = minimum(H1xy1[2,:])
        maxy1 = maximum(H1xy1[2,:])

        minx2 = minimum(H2xy2[1,:])
        maxx2 = maximum(H2xy2[1,:])
        miny2 = minimum(H2xy2[2,:])
        maxy2 = maximum(H2xy2[2,:])

        # Note the range in y must be common to both images
        miny = max(miny1, miny2)
        maxy = min(maxy1, maxy2)

    else
        # Apply a final translation adjustment to both transforms so that the final
        # rectified images do not end up with -ve pixel coordinates.  We do this by
        # determining the bounds of both images with the existing rectification and
        # applying a translation so that we end up with non-negative pixel
        # coordinates.  Note that the y translation must be common to both
        # images.  Also, we want the minimum x and y coordinates to be 1, not zero.
        (minx1, maxx1, miny1, maxy1) = transformedimagebounds(img1, H1)
        (minx2, maxx2, miny2, maxy2) = transformedimagebounds(img2, H2)
        #    miny = min(miny1, miny2)

        # Set  miny to be the top of the common area
        miny = max(miny1, miny2)
    end

    T1 = [1.0  0.0  -minx1+1
          0.0  1.0  -miny+1
          0.0  0.0    1.0   ]

    T2 = [1.0  0.0  -minx2+1
          0.0  1.0  -miny+1
          0.0  0.0    1.0   ]

    H1 = T1*H1                  # Final rectification transforms.
    H2 = T2*H2

    P1r = H1*P1                 # Adjusted projection matrices.
    P2r = H2*P2


    if AFFINECORRECT
        # Apply rectification transforms to the image points to
        # determine the range of disparities in the x direction and to
        # check that there are no serious shifts in the y direction.
        H1xy1 = hnormalise(H1*xy1)
        H2xy2 = hnormalise(H2*xy2)

        # Truncate specfied percentage off the ends of the histogram
        # of x disparity values to eliminate outliers and establish a
        # range of disparities that need to be searched for 3D
        # reconstruction.
        disparities = H2xy2[1,:] - H1xy1[1,:]
        h = histtruncate(disparities, disparitytruncation)
        dmin = floor(Int, minimum(h))
        dmax = ceil(Int, maximum(h))
        dmed = round(Int, median(disparities))

    else
        dmin = zeros(Int,0); dmax= zeros(Int,0); dmed=zeros(Int,0);
    end

    return (P1r, H1, P2r, H2, dmin, dmax, dmed)
end


#---------------------------------------------------------------------
"""
transformedimagebounds

Find where the corners of an image are transformed to by transform H and
return the bounds.

```
Usage:  (minx, maxx, miny, maxy) = transformedimagebounds(img::Array, H::Array)
        (minx, maxx, miny, maxy) = transformedimagebounds(sze::Tuple, H::Array)

Arguments:   img - An array storing the image or
             sze - A tuple giving the size of the image.
               H - The transforming homography, a 3x3 matrix.

Returns:
      minx, maxx - The range of x, y coords (range of column and row coords) of
      miny, maxy   the transformed image.

```
"""
function  transformedimagebounds(img, H)
    transformedimagebounds(size(img), H)
end

function  transformedimagebounds(sze::Tuple, H::Array)

    rows, cols = sze[1], sze[2]

    # Image corners
    xy = [1  cols cols   1
          1   1   rows  rows
          1   1    1     1   ]

    # Transformed image corners
    Hxy = hnormalise(H*xy)

    minx = minimum(Hxy[1,:]); maxx = maximum(Hxy[1,:])
    miny = minimum(Hxy[2,:]); maxy = maximum(Hxy[2,:])

    return (minx, maxx, miny, maxy)
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
function idealimagepts(C::Camera, xy::Union{Array{T,2}, Vector{T}}) where T <: Real

    if size(xy, 1) != 2
        error("Image xy data must be a 2 x N array")
    end

    # Reverse the projection process as used by cameraproject()

    # Subtract principal point and divide by focal length to get normalised,
    # distorted image coordinates.  Note skew represents the 2D shearing
    # coefficient times fx
    y_d = (xy[2,:] .- C.ppy)/C.fy
    x_d = (xy[1,:] .- C.ppx - y_d*C.skew)/C.fx

    # Compute the inverse of the lens distortion effect by computing the
    # 'forward' direction lens distortion at this point and then subtracting
    # this from the current point.

    # Radial distortion factor. Here the squared radius is computed from the
    # already distorted coordinates.  The approximation we are making here is to
    # assume that the distortion is locally constant.
    rsqrd = x_d.^2 .+ y_d.^2
    r_d = 1 .+ C.k1*rsqrd .+ C.k2*rsqrd.^2 .+ C.k3*rsqrd.^3

    # Tangential distortion component, again computed from the already distorted
    # coords.
    dtx = 2*C.p1*x_d.*y_d          .+ C.p2*(rsqrd .+ 2*x_d.^2)
    dty = C.p1*(rsqrd .+ 2*y_d.^2) .+ 2*C.p2*x_d.*y_d

    # Subtract the tangential distortion components and divide by the radial
    # distortion factor to get an approximation of the undistorted normalised
    # image coordinates (with no skew)
    x_n = (x_d .- dtx)./r_d
    y_n = (y_d .- dty)./r_d

    # Finally project back to pixel coordinates.
    x_p = C.ppx .+ x_n*C.fx + y_n*C.skew
    y_p = C.ppy .+ y_n*C.fy

    return xyideal = [x_p
                      y_p]
end

#--------------------------------------------------------------------------
"""
solvestereopt - Homogeneous linear solution of a stereo point
```
Usage:
        pt = solvestereopt(xy, P)
        pt = solvestereopt(xy, C)
        (pt, xy_reproj) = solvestereopt(xy, P, reprojecterror=true)
        (pt, xy_reproj) = solvestereopt(xy, C, reprojecterror=true)

Multiview stereo: Solves 3D location of a point given image coordinates of
that point in two, or more, images.

Arguments:    xy - 2xN matrix of x, y image coordinates, one column for
                   each camera.
               P - N array of corresponding 3x4 image projection
                   matrices, or
               C - An N array of Camera structures
  reprojecterror - Boolean flag indicating whether reprojection errors should be retunred.

Returns:      pt - 3D location in space returned in normalised
                   homogeneous coordinates (a 4-vector with last element = 1)
       xy_reproj - 2xN matrix of reprojected image coordinates.
```

See also: idealimagepts(), camstruct2projmatrix(), Camera()
"""
function solvestereopt(xy::Array{T1,2}, P::Array{Array{T2,2}}; reprojecterror=false) where {T1<:Real,T2<:Real}

    (dim,N) = size(xy)
    @assert N == length(P) "Number of points and cameras must match"
    @assert dim == 2  "Image coordinates must be a 2xN array of inhomogeneous coords"
    @assert N >= 2  "Must have at least 2 camera views"

    # Build eqn of the form A*pt = 0
    A = zeros(2*N, 4)
    for n = 1:N
	A[2*n-1,:] = xy[1,n]*P[n][3,:] .- P[n][1,:]
	A[2*n  ,:] = xy[2,n]*P[n][3,:] .- P[n][2,:]
    end

    (~,~,v) = svd(A)
    pt = hnormalise(v[:,4])

    if reprojecterror
        # Project the point back into the source images to determine the residual
        # error
        xy_reproj = zeros(size(xy))
        for n = 1:N
            xy_reproj[:,n] = cameraproject(P[n], pt[1:3])
        end

        return pt, xy_reproj
    else
        return pt
    end

end

# Version where the camera data is an array of Camera structures
function solvestereopt(xy::Array{T,2}, C::Array{Camera}; reprojecterror=false) where T<:Real

    (dim,N) = size(xy)
    @assert N == length(C) "Number of points and cameras must match"
    @assert dim == 2  "Image coordinates must be a 2xN array of inhomogeneous coords"
    @assert N >= 2  "Must have at least 2 camera views"

    # Correct image points from each camera so that they correspond to image
    # points from an ideal camera with no lens distortion.  Also generate the
    # corresponding ideal projection matrices for each Camera structure
    xyideal = similar(xy) # copy xy so that xy is not changed
    for n = 1:N
        xyideal[:,n] = idealimagepts(C[n], xy[:,n])
    end
    P = [camera2projmatrix(c) for c in C]

    return solvestereopt(xyideal, P, reprojecterror=reprojecterror)
end

#---------------------------------------------------------------------
"""
undistortimage - Removes lens distortion from an image

```
Usage 1:  nimg = undistortimage(img, f, ppx, ppy, k1, k2, k3, p1, p2)

Arguments:
         img - Image to be corrected.
           f - Focal length in terms of pixel units
               (focal_length_mm/pixel_size_mm)
    ppx, ppy - Principal point location in pixels.
  k1, k2, k3 - Radial lens distortion parameters.
      p1, p2 - Tangential lens distortion parameters.

Usage 2.  nimg = undistortimage(img, camera)

Arguments:
         img - Image to be corrected.
         cam - Camera structure.

Returns:
         nimg - Corrected image.
```

It is assumed that radial and tangential distortion parameters are
computed/defined with respect to normalised image coordinates corresponding
to an image plane 1 unit from the projection centre.  This is why the
focal length is required.

"""
function undistortimage(img::Array{T,2}, f::Real,
                       ppx::Real, ppy::Real, k1::Real, k2::Real, k3::Real,
                       p1::Real, p2::Real) where T<:Real

    # ** Allow for separate valuees of fx and fy ??
    #
    # Version for image represented by 2D array

    # Strategy: Generate a grid of coordinate values corresponding to
    # an ideal undistorted image.  We then apply the imaging process
    # to these coordinates, including lens distortion, to obtain the
    # actual distorted image locations.  In this process these
    # distorted image coordinates end up being stored in a matrix that
    # is indexed via the original ideal, undistorted coords.  Thus for
    # every undistorted pixel location we can determine the location
    # in the distorted image that we should map the grey value from.

    (rows, cols) = size(img)
    nimg = zeros(size(img))

    # Set up interpolant
    itp = interpolate(img, BSpline(Linear()), OnCell())

    # Interpolate values from distorted image locations to their
    # ideal locations.
    for xu = 1:cols, yu = 1:rows
        (xd, yd) = distortedcoords(xu, yu, f, ppx, ppy, k1, k2, k3, p1, p2)
        nimg[yu,xu] = itp[yd, xd]
    end

    return nimg
end

function undistortimage(img, cam::Camera)
    undistortimage(img, cam.fx, cam.ppx, cam.ppy, cam.k1, cam.k2, cam.k3, cam.p1, cam.p2)
end

# Version for image represented by 3D array
function undistortimage(img::Array{T,3}, f::Real,
                       ppx::Real, ppy::Real, k1::Real, k2::Real, k3::Real,
                       p1::Real, p2::Real) where T<:Real

    (rows, cols, channels) = size(img)
    nimg = zeros(size(img))

    # Interpolate values from distorted image locations to their
    # ideal locations. ** should construct array of interpolants **
    for chan = 1:channels
        itp = interpolate(img[:,:,chan], BSpline(Linear()), OnCell())

        for xu = 1:cols, yu = 1:rows
            (xd, yd) = distortedcoords(xu, yu, f, ppx, ppy, k1, k2, k3, p1, p2)

            if xd >= 1 && xd <= cols && yd >= 1 && yd <= rows
                nimg[yu, xu, chan] = itp[yd, xd]
            end
        end
    end

    return nimg
end

# Unexported function used by undistortimage()

function distortedcoords(xu::Real, yu::Real, f::Real, ppx::Real, ppy::Real,
                         k1::Real, k2::Real, k3::Real, p1::Real, p2::Real)

    # Convert undistorted x, y values to normalised values with the
    # origin at the principal point.  Dividing pixel coordinates by
    # the focal length (defined in pixels) gives us normalised coords
    # corresponding to z = 1
    x = (xu-ppx)/f
    y = (yu-ppy)/f

    # Radial lens distortion component
    r2 = x^2 + y^2                  # Squared normalized radius.
    dr = k1*r2 + k2*r2^2 + k3*r2^3  # Distortion scaling factor.

    # Tangential distortion component (Beware of different p1,p2
    # orderings used in the literature)
    dtx =    2*p1*x*y      +  p2*(r2 + 2*x^2)
    dty = p1*(r2 + 2*y^2)  +    2*p2*x*y

    # Apply the radial and tangential distortion components to x and y
    x = x + dr*x + dtx
    y = y + dr*y + dty

    # Now rescale by f and add the principal point back to get
    # distorted x and y coordinates.
    xd = x*f + ppx
    yd = y*f + ppy

    return (xd, yd)

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

"""
function hline(p1i::Vector, p2i::Vector, linestyle::String="b-")
    # Case when 2 homogeneous points are supplied

    p1 = p1i./p1i[3]    # make sure homogeneous points lie in z=1 plane
    p2 = p2i./p2i[3]

#    hold(true)
    plot([p1[1], p2[1]], [p1[2], p2[2]], linestyle);
end

# Case when homogeneous line is supplied
function hline(li::Vector, linestyle::String="b-")

    l = li./li[3]   # ensure line in z = 1 plane (not needed?)

    if abs(l[1]) > abs(l[2])   # line is more vertical
        p1 = hcross(l, [0, -1, PyPlot.ylim()[1]])
        p2 = hcross(l, [0, -1, PyPlot.ylim()[2]])

    else                       # line more horizontal
        p1 = hcross(l, [-1, 0, PyPlot.xlim()[1]])
        p2 = hcross(l, [-1, 0, PyPlot.xlim()[2]])
    end

#    hold(true)
    plot([p1[1], p2[1]], [p1[2], p2[2]], linestyle);
end

#--------------------------------------------------------------------------
"""
plotcamera - Plots graphical representation of camera(s) showing pose.

```
Usage: plotcamera(C, l; col=[0,0,1], plotCamPath=false, fig=nothing)

Arguments:
           C - Camera structure (or array of Camera structures)
           l - The length of the sides of the rectangular cone indicating
               the camera's field of view.
Keyword Arguments:
         col - Optional three element vector specifying the RGB colour to
               use. Defaults to blue.
 plotCamPath - Optional flag true/false to plot line joining camera centre
               positions. If omitted or empty defaults to false.
         fig - Optional figure number to be used. If not specified a new
               figure is created.
```

The function plots into the current figure a graphical representation of one
or more cameras showing their pose.  This consists of a rectangular cone,
with its vertex at the camera centre, indicating the camera's field of view.
The camera's coordinate X and Y axes are also plotted at the camera centre.

See also: Camera
"""
function plotcamera(Ci, l; col=[0,0,1], plotCamPath=false, fig=nothing)

    if isa(Ci, Array)
        C = Ci
    else
        C = Vector(1)
        C[1] = Ci
    end

    figure(fig)
#    hold(true)
    for i = 1:length(C)

        if C[i].rows == 0 || C[i].cols == 0
            @warn("Camera rows and cols not specified")
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
        X = Tw_c[1:3,1]*l .+ C[i].P
        Y = Tw_c[1:3,2]*l .+ C[i].P
        Z = Tw_c[1:3,3]*l .+ C[i].P

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

#---------------------------------------------------------------------
"""
imgtrans - Homogeneous transformation of an image - no image scaling.

```
Applies a geometric transform to an image

Usage: newimg = imgtrans(img, T)

Arguments:
       img    - The image to be transformed.
       T      - The 3x3 homogeneous transformation matrix.

Returns:
      newimg  - The transformed image.
```

See also: imgtrans!()
"""
function imgtrans(img::Array{T,3}, H::Array{Float64,2}) where T <: Real
    # Version for image represented by 3D array

    # ** To do: Add options to allow/disallow rescaling translating etc
    # If rescaling and translating occur then the new homography should be returned.

    (rows, cols, channels) = size(img)
    invH = inv(H)

    # Determine where the corners of the image go
    (minx::Float64, maxx::Float64, miny::Float64, maxy::Float64) = transformedimagebounds(img, H)

#=
    if minx < 0.5 || miny < 0.5  # Allow 1/2 pixel 'slop'
        @printf("Warning: Part of the transformed image has coordinates < 1\n")
        @printf("minx %f maxx %f miny %f maxy %f\n", minx, maxx, miny, maxy)
    end
=#
    nrows = ceil(Int, maxy)
    ncols = ceil(Int, maxx)

    newimg = zeros(Float64, nrows, ncols, channels)
    Hxy = zeros(3)

    for chan = 1:channels
        # Set up interpolant Seems to be faster to repeat the loops
        # for each channel than it is to slice down through the
        # channels
        itp = interpolate(img[:,:,chan], BSpline(Linear()), OnCell())

        for c = 1:ncols, r = 1:nrows
            #  Hxy = invH*[c,r,1]  # The following is _much_ faster
            Hxy[1] = invH[1,1]*c + invH[1,2]*r + invH[1,3]
            Hxy[2] = invH[2,1]*c + invH[2,2]*r + invH[2,3]
            Hxy[3] = invH[3,1]*c + invH[3,2]*r + invH[3,3]
            x = Hxy[1]/Hxy[3]
            y = Hxy[2]/Hxy[3]
            if x >=1 && x <= cols && y >=1 && y <= rows
               newimg[r,c, chan] = itp[y,x]
            end
        end
    end

    return newimg
end

# Version for image represented by 2D array
function imgtrans(img::Array{T,2}, H::Array{Float64,2}) where T <: Real

    (rows,cols) = size(img)
    invH = inv(H)

    # Determine where the corners of the image go
    (minx::Float64, maxx::Float64, miny::Float64, maxy::Float64) = transformedimagebounds(img, H)

#=
    if minx < 0.5 || miny < 0.5   # Allow 1/2 pixel 'slop'
        @printf("Warning: Part of the image has transformed coordinates < 1\n")
        @printf("minx %f maxx %f miny %f maxy %f\n", minx, maxx, miny, maxy)
    end
=#
    nrows = ceil(Int, maxy)
    ncols = ceil(Int, maxx)

    newimg = zeros(Float64, nrows, ncols)

    # Set up interpolant
    itp = interpolate(img, BSpline(Linear()), OnCell())
    Hxy = zeros(3)

    for c = 1:ncols, r = 1:nrows
#        Hxy = invH*[c,r,1]  # The following is _much_ faster
        Hxy[1] = invH[1,1]*c + invH[1,2]*r + invH[1,3]
        Hxy[2] = invH[2,1]*c + invH[2,2]*r + invH[2,3]
        Hxy[3] = invH[3,1]*c + invH[3,2]*r + invH[3,3]

        x = Hxy[1]/Hxy[3]
        y = Hxy[2]/Hxy[3]

        if x >=1 && x <= cols && y >=1 && y <= rows
            newimg[r,c] = itp[y,x]
        end
    end

    return newimg
end

"""
imgtrans! - Homogeneous transformation of an image - no image scaling.

```
Applies a geometric transform to an image

Usage: imgtrans!(newimg, img, T)

Arguments:
      newimg - Buffer for storing the transformed image.
         img - The image to be transformed.
           T - The 3x3 homogeneous transformation matrix.

Returns: nothing

```
See also: imgtrans()
"""
function imgtrans!(newimg::Array{T,2}, img::Array{T,2}, H::Array{Float64,2}) where T <: Real
    # Version for image represented by 2D array inplace

    (rows,cols) = size(img)
    invH = inv(H)

    (nrows, ncols) = size(newimg)

    # Set up interpolant
    itp = interpolate(img, BSpline(Linear()), OnCell())
    Hxy = zeros(3)

    for c = 1:ncols, r = 1:nrows
#        Hxy = invH*[c,r,1]  # The following is _much_ faster
        Hxy[1] = invH[1,1]*c + invH[1,2]*r + invH[1,3]
        Hxy[2] = invH[2,1]*c + invH[2,2]*r + invH[2,3]
        Hxy[3] = invH[3,1]*c + invH[3,2]*r + invH[3,3]

        x = Hxy[1]/Hxy[3]
        y = Hxy[2]/Hxy[3]

        if x >=1 && x <= cols && y >=1 && y <= rows
            newimg[r,c] = itp[y,x]
        else
            newimg[r,c] = zero(T)
        end
    end

    return nothing
end

#=
# Images.Image version  *** Have to sort out rows cols order !!! ***
function imgtrans(img::Images.Image, H::Array)

    (rows,cols) = size(img)

    invH = inv(H)

    # Determine where the corners of the image go
    (minx, maxx, miny, maxy) = transformedimagebounds(img, H)

    if minx < 1 || miny < 1
        @printf("Warning: Part of the image has transformed coordinates < 1\n")
    end

    nrows = ceil(Int, maxy)
    ncols = ceil(Int, maxx)

    # Need to construct a zero image of same size and type as img
    newimg = zeros(Float64, nrows, ncols)

    # Set up interpolant
    itp = interpolate(img, BSpline(Linear()), OnCell())

    for c = 1:ncols, r = 1:nrows
        Hxy = invH*[c,r,1]
        newimg[r,c] = itp[Hxy[2]/Hxy[3], Hxy[1]/Hxy[3]]
    end

    return newimg
end
=#
