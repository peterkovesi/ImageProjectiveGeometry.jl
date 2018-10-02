#=--------------------------------------------------------------------

transforms3d - Functions for performing 3D transformations


Copyright (c) 2015-2017 Peter Kovesi
pk@peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK August    2015
   April     2017  Cleanups and optimisations
   September 2018  Updated for v0.7/v1.0

---------------------------------------------------------------------=#

export trans, rotx, roty, rotz, invht
export dhtrans, homotrans
export inveuler, invrpy
export angleaxis, angleaxis2matrix, angleaxisrotate, normaliseangleaxis
export matrix2angleaxis, matrix2quaternion, quaternion2matrix
export quaternion, quaternionconjugate, quaternionproduct, quaternionrotate
export vector2quaternion

using LinearAlgebra

#---------------------------------------------------------------------
"""
trans - Homogeneous transformation for a translation by x, y, z
```
Usage: T = trans(x, y, z)
       T = trans(v)

Arguments:  x,y,z - translations in x,y and z, or alternatively
            v     - 3x1 vector or array defining x, y and z.
Returns:    T     - 4x4 homogeneous transformation matrix.
```
See also: rotx(), roty(), rotz(), invht()
"""
function trans(x::Real, y::Real, z::Real)

  T = [ 1.0  0.0  0.0   float(x)
        0.0  1.0  0.0   float(y)
        0.0  0.0  1.0   float(z)
        0.0  0.0  0.0  1.0 ]
end  

function trans(t::Array{T1}) where T1<:Real
    
    if length(t) != 3
        error("Translation must be a 3x1 vector or array")
    end
    
    T = [ 1.0  0.0  0.0  float(t[1])
          0.0  1.0  0.0  float(t[2])
          0.0  0.0  1.0  float(t[3])
          0.0  0.0  0.0  1.0 ]
end  

#----------------------------------------------------------------------
"""
rotx - Homogeneous transformation for a rotation about the x axis
```
Usage: T = rotx(theta)

Argument:  theta  - rotation about x axis
Returns:    T     - 4x4 homogeneous transformation matrix
```
See also: trans(), roty(), rotz(), invht()
"""
function rotx(theta::Real)
    T = [ 1.0     0.0         0.0      0.0
          0.0  cos(theta) -sin(theta)  0.0
          0.0  sin(theta)  cos(theta)  0.0
          0.0     0.0         0.0      1.0 ]
end

#---------------------------------------------------------------------
"""
roty - Homogeneous transformation for a rotation about the y axis
```
Usage: T = roty(theta)

Argument:  theta  - rotation about y axis
Returns:    T     - 4x4 homogeneous transformation matrix
```
See also: trans(), rotx(), rotz(), invht()
"""
function roty(theta::Real)
    T = [ cos(theta)  0.0  sin(theta)  0.0
              0.0     1.0      0.0     0.0
         -sin(theta)  0.0  cos(theta)  0.0
              0.0     0.0      0.0     1.0 ]
end

#---------------------------------------------------------------------
"""
rotz - Homogeneous transformation for a rotation about the z axis
```
Usage: T = rotz(theta)

Argument:  theta  - rotation about z axis
Returns:    T     - 4x4 homogeneous transformation matrix
```
See also: trans(), rotx(), roty(), invht()
"""
function rotz(theta::Real)
    T = [ cos(theta) -sin(theta)  0.0   0.0
          sin(theta)  cos(theta)  0.0   0.0
              0.0        0.0      1.0   0.0
              0.0        0.0      0.0   1.0 ]
end

#---------------------------------------------------------------------
"""
invht - Inverse of a homogeneous transformation matrix
```
Usage:  Tinv = invht(T)

Argument:   T    - 4x4 homogeneous transformation matrix
Returns:    Tinv - inverse
```
See also: trans(), rotx(), roty(), rotz()
"""
function invht(T::Array{T1,2}) where T1 <: Real
    
    if size(T) != (4,4)
        error("T must be 4 x 4")
    end
    
    A = T[1:3,1:3]
    
    Tinv = [      A'        -A'*T[1:3,4]
            0.0  0.0  0.0      1.0      ]
end

#---------------------------------------------------------------------
"""
angleaxis2matrix - Converts angle-axis descriptor to 4x4 homogeneous
transformation  matrix
```
Usage:     T = angleaxis2matrix(t)

Argument:  t - 3x1 vector specifying the rotation axis and having 
               magnitude equal to the rotation angle in radians.
Returns:   T - 4x4 Homogeneous transformation matrix.
```
See also: matrix2angleaxis(), angleaxisrotate(), angleaxis(), normaliseangleaxis()
"""
function  angleaxis2matrix(t::Vector{T1}) where T1 <: Real

    if length(t) != 3
        error("Angle-axis must be a 3-vector")
    end

    theta = norm(t)
    if theta < eps()    # If the rotation is very small...
        T = [  1.0  -t[3]   t[2]   0.0
              t[3]   1.0   -t[1]   0.0
             -t[2]   t[1]    1.0   0.0
               0.0   0.0     0.0   1.0]
        
        return T
    end
    
    # Otherwise set up standard matrix, first setting up some convenience
    # variables
    tn = t./theta;  x = tn[1]; y = tn[2]; z = tn[3];
    
    c = cos(theta); s = sin(theta); C = 1.0-c
    xs = x*s;   ys = y*s;   zs = z*s
    xC = x*C;   yC = y*C;   zC = z*C
    xyC = x*yC; yzC = y*zC; zxC = z*xC

    T = [ x*xC+c   xyC-zs   zxC+ys  0.0
          xyC+zs   y*yC+c   yzC-xs  0.0
          zxC-ys   yzC+xs   z*zC+c  0.0
           0.0       0.0     0.0    1.0]
end
    
#---------------------------------------------------------------------
"""    
angleaxisrotate - Uses angle axis descriptor to rotate vectors
```
Usage: v2 = angleaxisrotate(t, v)

Arguments:  t  - 3-vector defining rotation axis and having magnitude 
                 equal to the rotation angle in radians.
            v  - 4xn matrix of homogeneous 4-vectors to be rotated or
                 3xn matrix of inhomogeneous 3-vectors to be rotated.
Returns:    v2 - The rotated vectors. 
```
See also: matrix2angleaxis(), angleaxis(), angleaxis2matrix(), normaliseangleaxis()
"""
function angleaxisrotate(t::Vector, v::Array)
    
    ndim = size(v,1)
    T = angleaxis2matrix(t)

    if ndim == 3
      v2 = T[1:3,1:3]*v

    elseif ndim == 4
      v2 = T*v

    else
      error("v must be 4xN or 3xN");
    end

    return v2
end

#---------------------------------------------------------------------
"""
dhtrans - Computes Denavit Hartenberg matrix

This function calculates the 4x4 homogeneous transformation matrix,
representing the Denavit Hartenberg matrix of a robot arm link, given
link parameters of joint angle, length, joint offset and twist.
```
Usage: T = dhtrans(theta, offset, length, twist)
 
Arguments:  theta - joint angle (rotation about local z)
           offset - offset (translation along z)
           length - translation along link x axis
            twist - rotation about link x axis

Returns:        T - 4x4 Homogeneous transformation matrix.
```
See also: trans(), rotx(), roty(), rotz(), invht()
"""
function dhtrans(theta::Real, offset::Real, length::Real, twist::Real)

T = [ cos(theta) -sin(theta)*cos(twist)  sin(theta)*sin(twist) length*cos(theta)
      sin(theta)  cos(theta)*cos(twist) -cos(theta)*sin(twist) length*sin(theta)
          0.0           sin(twist)             cos(twist)         offset
          0.0              0.0                    0.0               1.0          ]
end

#---------------------------------------------------------------------
"""
homotrans - Homogeneous transformation of points/lines

Function to perform a transformation on 2D or 3D homogeneous or
inhomogeneous coordinates. The resulting coordinates are normalised to
have a homogeneous scale of 1.

```
Usage:     t = homotrans(P, v; checkinf = false)

Arguments:
          P  - 3 x 3 or 4 x 4 homogeneous transformation matrix
          v  - 3 x n or 4 x n matrix of homogeneous coordinates or ...
               2 x n or 3 x n matrix of inhomogeneous coordinates
Keyword Argument:
    checkinf - If true code will check if the homogeneous coordinates
               of any transformed point is at infinity (scale parameter
               is zero). If so, no normalisation of the point is attempted.

Returns   t  - Transformed coordinates. Homogeneous or inhomogeneous 
               to match input data v.
```

The input v is assumed inhomogeneous if the number of rows is equal to
size(P,2) - 1

"""
function homotrans(P::Array, v::Array; checkinf = false)
    # ? Should normalisation of scale test for points at infinity ?
    
    (dim, nPts) = size(v,1,2)

    if size(P) == (dim+1, dim+1)   # Assume data is inhomogeneous
        t = P[1:dim, 1:dim]*v
        scale = P[dim+1,1:dim]*v
        
        # Add the translation component of the transform to the scale
        # then add the translation component to t and divide by the
        # scale to normalise the result.
        if checkinf   # Only normalise if the scale parameter is non-zero
            for j = 1:nPts
                scale[j] += P[dim+1,dim+1]
                if abs(scale[j]) > eps()
                    for i = 1:dim
                        t[i,j] = (t[i,j] + P[i,dim+1])/scale[j]
                    end
                end
            end

        else   # Do not check for points at infinity
            for j = 1:nPts
                scale[j] += P[dim+1,dim+1]
                
                for i = 1:dim
                    t[i,j] = (t[i,j] + P[i,dim+1])/scale[j]
                end
            end
        end
        
    elseif size(P) == (dim, dim)  # Homogeneous case
        t = P*v       
        
        # Normalise to scale of 1 
        if checkinf   # Only normalise if the scale parameter is non-zero
            for i = 1:nPts
                if abs(t[dim,i]) > eps()
                    for r = 1:dim-1  
                        t[r,i] /= t[dim,i]
                    end
                    t[dim,i] = 1.0
                end
            end

        else   # Do not check for points at infinity
            for i = 1:nPts
                for r = 1:dim-1  
                    t[r,i] /= t[dim,i]
                end
                t[dim,i] = 1.0
            end
        end

    else
        error("Transformation matrix and point data dimensions do not match")
    end

    return t
end
    
#---------------------------------------------------------------------
"""    
inveuler - Inverse of Euler transform
```
Usage:  (euler1, euler2) = inveuler(T)

Argument:  T - 4x4 Homogeneous transformation matrix or 3x3 rotation matrix

Returns: euler1 = [phi1, theta1, psi1] - the 1st solution and,
         euler2 = [phi2, theta2, psi2] - the 2nd solution

     rotz(phi1)   * roty(theta1)    * rotz(psi1)      = T
  rotz(euler1[1]) * roty(euler1[2]) * rotz(euler1[3]) = T
```
See also: invrpy(), invht(), rotx(), roty(), rotz()
"""
function inveuler(T::Array{T1,2}) where T1 <: Real
    #=
    Reference: Richard P. Paul  
    Robot Manipulators: Mathematics, Programming and Control.
    MIT Press 1981. Page 68
    =#

    phi1 = atan(T[2,3], T[1,3])
    phi2 = phi1 + pi
    
    theta1 = atan(cos(phi1)*T[1,3] + sin(phi1)*T[2,3], T[3,3])
    theta2 = atan(cos(phi2)*T[1,3] + sin(phi2)*T[2,3], T[3,3])
    
    psi1 = atan(-sin(phi1)*T[1,1] + cos(phi1)*T[2,1], 
                 -sin(phi1)*T[1,2] + cos(phi1)*T[2,2])
    psi2 = atan(-sin(phi2)*T[1,1] + cos(phi2)*T[2,1], 
                 -sin(phi2)*T[1,2] + cos(phi2)*T[2,2])
    
    euler1 = [phi1, theta1, psi1]
    euler2 = [phi2, theta2, psi2]
    return (euler1, euler2)
end
    
#---------------------------------------------------------------------
"""
invrpy - Inverse of Roll Pitch Yaw transform
```
Usage:  (rpy1, rpy2) = invrpy(RPY)
 
Argument:  RPY - 4x4 Homogeneous transformation matrix or 3x3 rotation matrix

Returns:  rpy1 = [phi1, theta1, psi1] - the 1st solution and
          rpy2 = [phi2, theta2, psi2] - the 2nd solution

   rotz(phi1)  * roty(theta1)  * rotx(psi1)    = RPY
 rotz(rpy1[1]) * roty(rpy1[2]) * rotx(rpy1[3]) = RPY
```
See also: inveuler(), invht(), rotx(), roty(), rotz()
"""
function invrpy(RPY::Array{T,2}) where T <: Real
    #=
    Reference: Richard P. Paul  
    Robot Manipulators: Mathematics, Programming and Control.
    MIT Press 1981. Page 70
    =#
    # Z rotation 
    phi1 = atan(RPY[2,1], RPY[1,1])
    phi2 = phi1 + pi
    
    # Y rotation
    theta1 = atan(-RPY[3,1], cos(phi1)*RPY[1,1] + sin(phi1)*RPY[2,1])
    theta2 = atan(-RPY[3,1], cos(phi2)*RPY[1,1] + sin(phi2)*RPY[2,1])
    
    # X rotation
    psi1 = atan(sin(phi1)*RPY[1,3] - cos(phi1)*RPY[2,3], 
                -sin(phi1)*RPY[1,2] + cos(phi1)*RPY[2,2]);
    psi2 = atan(sin(phi2)*RPY[1,3] - cos(phi2)*RPY[2,3], 
                -sin(phi2)*RPY[1,2] + cos(phi2)*RPY[2,2]);
    
    rpy1 = [phi1, theta1, psi1]
    rpy2 = [phi2, theta2, psi2]
    return (rpy1, rpy2)
end
    
#---------------------------------------------------------------------
"""
matrix2angleandaxis - Decompose homogeneous matrix to angle and axis 
```
Usage: (ang, axis) = matrix2angleandaxis(T)

Argument:   T - 4x4 Homogeneous transformation matrix, or 3x3 rotation matrix.
Returns:  ang - Rotation angle in radians.
         axis - Unit 3-vector defining the rotation axis.
```
Note that only the top left 3x3 rotation component of T is used, any
translation component in T is ignored.

See also: angleaxis2matrix(), angleaxisrotate(), angleaxis(), normaliseangleaxis()
"""
function matrix2angleandaxis(T::Array)

    # This code follows the implementation suggested by Hartley and Zisserman
    R = T[1:3, 1:3]   # Extract rotation part of T
    
    # Trap case where rotation is very small.  ( See angleaxis2matrix() )
    Reye = R - I
    if norm(Reye) < 1e-8
        t = [T[3,2], T[1,3], T[2,1]]
        return norm(Reye), t
    end

    # Otherwise find rotation axis as the eigenvector having unit eigenvalue
    # Solve (R-I)v = 0
#    (eigval, eigvec) = eig(Reye)
    F = eigen(Reye)

    # Find index of eigenvalue with smallest magnitude, the
    # corresponding eigenvector will be the axis.
    absval = abs.(F.values)
    ind = sortperm(absval)          # Find index of smallest one
    if absval[ind[1]] > 0.001       # Hopefully it is close to 0
        @warn("Rotation matrix is dubious")
    end
    
    axis = real(F.vectors[:,ind[1]])    # Extract appropriate eigenvector

    if abs(norm(axis) - 1.0) > .0001    # Debug
        @warn("Non unit rotation axis")
    end
    
    # Now determine the rotation angle
    twocostheta = tr(R)-1.0
    twosinthetav = [R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2]]
    twosintheta = dot(axis, twosinthetav)

    # Note use of twosinetheta to convert it from a 1x1 matrix to a scalar
    theta = atan(twosintheta, twocostheta)

    return theta, axis
end
    
#---------------------------------------------------------------------
"""
matrix2angleaxis - Homogeneous matrix to angle-axis description
```
Usage: t = matrix2angleaxis(T)

Argument:   T - 4x4 Homogeneous transformation matrix, or 3x3 rotation matrix.
Returns:    t - 3x1 column vector giving rotation axis with magnitude equal
                to the rotation angle in radians.
```
Note that only the top left 3x3 rotation component of T is used, any
translation component in T is ignored.

See also: angleaxis2matrix(),  angleaxisrotate(), angleaxis(), normaliseangleaxis()

"""
function matrix2angleaxis(T::Array)
    (theta, axis) = matrix2angleandaxis(T)
    return theta*axis
end
    
#---------------------------------------------------------------------
"""
matrix2quaternion - Homogeneous matrix to quaternion

Converts 4x4 homogeneous rotation matrix to quaternion
```
Usage: Q = matrix2quaternion(T)

Argument:   T - 4x4 Homogeneous transformation matrix
Returns:    Q - 4-vector quaternion in the form [w, xi, yj, zk]
```
See also: quaternion2matrix()
"""
function matrix2quaternion(T::Array)
    (theta, axis) = matrix2angleandaxis(T)
    return [cos(theta/2); axis*sin(theta/2)]
end
    
#---------------------------------------------------------------------
"""
angleaxis - Constructs angle-axis descriptor
```
Usage: t = angleaxis(theta, axis)

Arguments: theta - angle of rotation.
           axis  - 3x1 vector defining axis of rotation.
Returns:   t     - 3-vector giving rotation axis with magnitude equal to the
                   rotation angle in radians.
```
See also: matrix2angleaxis(), angleaxisrotate(), angleaxis2matrix(), normaliseangleaxis()
"""
function angleaxis(theta::Real, axis::Vector)
    
    if length(axis) != 3
        error("Axis must be a 3 vector or array")
    end
    
    ax = normalize(axis)  # Ensure unit magnitude
    
    # Normalise theta to lie in the range -pi to pi to ensure one-to-one mapping
    # between angle-axis descriptor and resulting rotation. 
    theta = rem(theta, 2*pi)  # Remove multiples of 2pi
    
    if theta > pi
        theta = theta - 2*pi
    elseif  theta < -pi
        theta = theta + 2*pi
    end
    
    return theta*ax
end

#---------------------------------------------------------------------
"""
quaternion  - Construct quaternion 

```
Usage:  Q = quaternion(theta, axis)

Arguments: theta - angle of rotation
           axis  - 3-vector defining axis of rotation.
Returns:   Q     - a 4-vector quaternion in the form [w, xi, yj, zk]
```
See also:  quaternion2matrix(), matrix2quaternion(), quaternionrotate()
"""
function quaternion(theta::Real, axis::Vector)
    
    if length(axis) != 3
        error("Axis must be a 3 vector or array")
    end

    Q = zeros(4)    
    Q[1] = cos(theta/2)
    Q[2:4] = sin(theta/2)*normalize(axis)
    return Q
end
    
#---------------------------------------------------------------------
"""    
normaliseangleaxis - Normalises angle-axis descriptor.

Function normalises theta so that it lies in the range -pi to pi to ensure
one-to-one mapping between angle-axis descriptor and resulting rotation.
```
Usage: t2 = normaliseangleaxis(t)

Argument:   t  - 3-vector giving rotation axis with magnitude equal to the
                 rotation angle in radians.
Returns:    t2 - Normalised angle-axis descriptor
```
See also: matrix2angleaxis(), angleaxis(), angleaxis2matrix(), angleaxisrotate()
"""
function normaliseangleaxis(t::Vector)
    
    if length(t) != 3
        error("Axis must be a 3-vector")
    end
    
    theta = norm(t)
    axis = t/theta
    
    theta = rem(theta, 2pi)  # Remove multiples of 2pi
    
    if theta > pi              # Note theta cannot be -ve
        theta = theta - 2pi; 
    end
    
    return theta*axis
end

#---------------------------------------------------------------------
"""
quaternion2matrix - Quaternion to a 4x4 homogeneous transformation matrix
```
Usage:  T = quaternion2matrix(Q)

Argument: Q - a quaternion in the form [w, xi, yj, zk].
Returns:  T - 4x4 Homogeneous rotation matrix.
``` 
See also: matrix2quaternion(), quaternion(), quaternionrotate()
"""
function quaternion2matrix(Qu::Vector)
    
    if length(Qu) != 4
        error("Quaternion must be a 4 vector")
    end

    Q = normalize(Qu)    # Ensure Q has unit norm
    
    # Set up convenience variables
    w = Q[1]; x = Q[2]; y = Q[3]; z = Q[4]
    w2 = w^2; x2 = x^2; y2 = y^2; z2 = z^2
    xy = x*y; xz = x*z; yz = y*z
    wx = w*x; wy = w*y; wz = w*z
    
    T = [w2+x2-y2-z2  2*(xy - wz)  2*(wy + xz)   0.0
         2*(wz + xy)  w2-x2+y2-z2  2*(yz - wx)   0.0
         2*(xz - wy)  2*(wx + yz)  w2-x2-y2+z2   0.0
             0.0           0.0          0.0      1.0]

    return T
end
    
#---------------------------------------------------------------------
"""
quaternionconjugate - Conjugate of a quaternion
```
Usage: Qconj = quaternionconjugate(Q)

Argument: Q     - Quaternions in the form  Q = [Qw, Qi, Qj, Qk].
Returns:  Qconj - Conjugate.
```
See also: quaternion(), quaternionrotate(), quaternionproduct()
"""
function quaternionconjugate(Q::Vector)
    
    if length(Q) != 4
        error("Quaternion must be a 4 vector")
    end

    return [Q[1], -Q[2], -Q[3], -Q[4]]
end

#---------------------------------------------------------------------
"""
quaternionproduct - Computes product of two quaternions
```
Usage: Q = quaternionproduct(A, B)

Arguments: A, B - Quaternions assumed to be 4-vectors in the
                  form  A = [Aw, Ai, Aj, Ak].
Returns:   Q    - Quaternion product.
```
See also: quaternion(), quaternionrotate(), quaternionconjugate()
"""
function  quaternionproduct(A::Vector, B::Vector)

    if length(A) != 4 || length(B) != 4
        error("Quaternions must be 4 vectors")
    end

    Q = zeros(4,1)
    Q[1]  =  A[1]*B[1]  -  A[2]*B[2]  -  A[3]*B[3]  -  A[4]*B[4]
    Q[2]  =  A[1]*B[2]  +  A[2]*B[1]  +  A[3]*B[4]  -  A[4]*B[3]
    Q[3]  =  A[1]*B[3]  -  A[2]*B[4]  +  A[3]*B[1]  +  A[4]*B[2]
    Q[4]  =  A[1]*B[4]  +  A[2]*B[3]  -  A[3]*B[2]  +  A[4]*B[1]

    return Q
end

#---------------------------------------------------------------------
"""
quaternionrotate - Rotates a 3D vector by a quaternion 
```
Usage:   vnew = quaternionrotate(Q, v)

Arguments: Q - a quaternion in the form [w, xi, yj, zk]
           v - a vector to rotate, either an inhomogeneous 3-vector or a
               homogeneous 4-vector.
Returns:   vnew - rotated vector.
```
See also: matrix2quaternion(), quaternion2matrix(), quaternion()
"""
function quaternionrotate(Q::Vector, v::Vector)
    #=
    % Code forms the equivalent 3x3 rotation matrix from the quaternion and
    % applies it to a vector 
    %
    % Note that Qw^2 + Qi^2 + Qj^2 + Qk^2 = 1
    % So the top-left entry of the rotation matrix of
    %   Qw^2 + Qi^2 - Qj^2 - Qk^2
    % can be rewritten as
    %   Qw^2 + Qi^2 + Qj^2 + Qk^2 - 2Qj^2 - 2Qk^2
    % = 1 - 2Qj^2 - 2Qk^2
    %
    % Similar optimization applies to the other diagonal elements
    =#
    
    if length(Q) != 4 
        error("Quaternion must be a 4 vector")
    end

    if length(v) != 4 && length(v) != 3
        error("Input vector must be of length 3 or 4")
    end

    # Copy v to vnew to allocate space.  If v is a 4 element
    # homogeneous vector this also sets the homogeneous scale factor
    # of vnew. Ensure vnew is float so that we avoid Inexact errors
    vnew = float(copy(v))

    Qw = Q[1];  Qi = Q[2];  Qj = Q[3];  Qk = Q[4]
    
    t2 =   Qw*Qi
    t3 =   Qw*Qj
    t4 =   Qw*Qk
    t5 =  -Qi*Qi
    t6 =   Qi*Qj
    t7 =   Qi*Qk
    t8 =  -Qj*Qj
    t9 =   Qj*Qk
    t10 = -Qk*Qk
    vnew[1] = 2*( (t8 + t10)*v[1] + (t6 -  t4)*v[2] + (t3 + t7)*v[3] ) + v[1]
    vnew[2] = 2*( (t4 +  t6)*v[1] + (t5 + t10)*v[2] + (t9 - t2)*v[3] ) + v[2]
    vnew[3] = 2*( (t7 -  t3)*v[1] + (t2 +  t9)*v[2] + (t5 + t8)*v[3] ) + v[3]

    return vnew
end
    
#---------------------------------------------------------------------
""" 
vector2quaternion - Embeds 3-vector in a quaternion representation
```
Usage: Q = vector2quaternion(v)

Argument:  v - 3-vector.
Returns:   Q - Quaternion given by [0; v]
```
See also: quaternion(), quaternionrotate(), quaternionproduct(), quaternionconjugate()
"""
function vector2quaternion(v::Vector)

  if length(v) != 3
      error("v must be a 3-vector")
  end

  return [0.0, v[1], v[2], v[3]]
end
