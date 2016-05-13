Transforms Function Reference
==============================

## Index

* [trans](#trans) - Homogeneous transformation for a translation by x, y, z.
* [rotx](#rotx) - Homogeneous transformation for a rotation about the x axis.
* [roty](#roty) - Homogeneous transformation for a rotation about the y axis.
* [rotz](#rotz) - Homogeneous transformation for a rotation about the z axis.
* [dhtrans](#dhtrans) - Computes Denavit Hartenberg matrix.
* [homotrans](#homotrans) - Homogeneous transformation of points/lines.
* [invht](#invht) - Inverse of a homogeneous transformation matrix.
* [inveuler](#inveuler) - Inverse of Euler transform.
* [invrpy](#invrpy) - Inverse of Roll Pitch Yaw transform.
* [angleaxis](#angleaxis) - Constructs angle-axis descriptor.
* [normaliseangleaxis](#normaliseangleaxis) - Normalises angle-axis descriptor.
* [angleaxisrotate](#angleaxisrotate) - Uses angle axis descriptor to rotate 
vectors.
* [angleaxis2matrix](#angleaxis2matrix) - Converts angle-axis descriptor to 4x4 homogeneous transformation  matrix.
* [matrix2angleandaxis](#matrix2angleandaxis) - Decompose homogeneous matrix 
to angle and axis.
* [matrix2angleaxis](#matrix2angleaxis) - Homogeneous matrix to angle-axis description.
* [quaternion](#quaternion) - Construct quaternion.
* [quaternionconjugate](#quaternionconjugate) - Conjugate of a quaternion.
* [quaternionproduct](#quaternionproduct) - Computes product of two quaternions.
* [quaternionrotate](#quaternionrotate) - Rotates a 3D vector by a quaternion.
* [quaternion2matrix](#quaternion2matrix) - Quaternion to a 4x4 homogeneous transformation matrix.
* [matrix2quaternion](#matrix2quaternion) - Homogeneous matrix to quaternion.
* [vector2quaternion](#vector2quaternion) - Embeds 3-vector in a quaternion representation.

____

## trans - Homogeneous transformation for a translation by x, y, z

```
Usage: T = trans(x, y, z)
       T = trans(v)

Arguments:  x,y,z - translations in x,y and z, or alternatively
            v     - 3x1 vector or array defining x, y and z.
Returns:    T     - 4x4 homogeneous transformation matrix.
```

## rotx - Homogeneous transformation for a rotation about the x axis

```
Usage: T = rotx(theta)

Argument:  theta  - rotation about x axis
Returns:    T     - 4x4 homogeneous transformation matrix
```

## roty - Homogeneous transformation for a rotation about the y axis

```
Usage: T = roty(theta)

Argument:  theta  - rotation about y axis
Returns:    T     - 4x4 homogeneous transformation matrix
```

## rotz - Homogeneous transformation for a rotation about the z axis
```
Usage: T = rotz(theta)

Argument:  theta  - rotation about z axis
Returns:    T     - 4x4 homogeneous transformation matrix
```

## dhtrans - Computes Denavit Hartenberg matrix.

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

## homotrans - Homogeneous transformation of points/lines.

Function to perform a transformation on 2D or 3D homogeneous coordinates
The resulting coordinates are normalised to have a homogeneous scale of 1
```
Usage:     t = homotrans(P, v)

Arguments:
          P  - 3 x 3 or 4 x 4 homogeneous transformation matrix.
          v  - 3 x n or 4 x n matrix of homogeneous coordinates.

Returns   t  - Transformed homogeneous coordinates
```

## invht - Inverse of a homogeneous transformation matrix

```
Usage:  Tinv = invht(T)

Argument:   T    - 4x4 homogeneous transformation matrix
Returns:    Tinv - inverse
```

## invrpy - Inverse of Roll Pitch Yaw transform.

```
Usage:  (rpy1, rpy2) = invrpy(RPY)
 
Argument:  RPY - 4x4 Homogeneous transformation matrix or 3x3 rotation matrix

Returns:  rpy1 = [phi1, theta1, psi1] - the 1st solution and
          rpy2 = [phi2, theta2, psi2] - the 2nd solution

   rotz(phi1)  * roty(theta1)  * rotx(psi1)    = RPY
 rotz(rpy1[1]) * roty(rpy1[2]) * rotx(rpy1[3]) = RPY
```

## inveuler - Inverse of Euler transform.
```
Usage:  (euler1, euler2) = inveuler(T)

Argument:  T - 4x4 Homogeneous transformation matrix or 3x3 rotation matrix

Returns: euler1 = [phi1, theta1, psi1] - the 1st solution and,
         euler2 = [phi2, theta2, psi2] - the 2nd solution

     rotz(phi1)   * roty(theta1)    * rotz(psi1)      = T
  rotz(euler1[1]) * roty(euler1[2]) * rotz(euler1[3]) = T
```


## angleaxis - Constructs angle-axis descriptor.

```
Usage: t = angleaxis(theta, axis)

Arguments: theta - angle of rotation.
           axis  - 3x1 vector or array defining axis of rotation.
Returns:   t     - 3-vector giving rotation axis with magnitude equal to the
                   rotation angle in radians.
```

## normaliseangleaxis - Normalises angle-axis descriptor.

Function normalises theta so that it lies in the range -pi to pi to ensure
one-to-one mapping between angle-axis descriptor and resulting rotation.

```
Usage: t2 = normaliseangleaxis(t)

Argument:   t  - 3-vector giving rotation axis with magnitude equal to the
                 rotation angle in radians.
Returns:    t2 - Normalised angle-axis descriptor
```

## angleaxisrotate - Uses angle axis descriptor to rotate vectors

```
Usage: v2 = angleaxisrotate(t, v)

Arguments:  t  - 3-vector defining rotation axis and having magnitude 
                 equal to the rotation angle in radians.
            v  - 4xn matrix of homogeneous 4-vectors to be rotated or
                 3xn matrix of inhomogeneous 3-vectors to be rotated.
Returns:    v2 - The rotated vectors. 
```

## angleaxis2matrix - Converts angle-axis descriptor to 4x4 homogeneous transformation matrix

```
Usage:     T = angleaxis2matrix(t)

Argument:  t - 3x1 vector or array specifying the rotation axis and having 
               magnitude equal to the rotation angle in radians.
Returns:   T - 4x4 Homogeneous transformation matrix.
```

## matrix2angleandaxis - Decompose homogeneous matrix to angle and axis.

```
Usage: (ang, axis) = matrix2angleandaxis(T)

Argument:   T - 4x4 Homogeneous transformation matrix, or 3x3 rotation matrix.
Returns:  ang - Rotation angle in radians.
         axis - Unit 3-vector defining the rotation axis.
```

Note that only the top left 3x3 rotation component of T is used, any
translation component in T is ignored.


## matrix2angleaxis - Homogeneous matrix to angle-axis description.

```
Usage: t = matrix2angleaxis(T)

Argument:   T - 4x4 Homogeneous transformation matrix, or 3x3 rotation matrix.
Returns:    t - 3x1 column vector giving rotation axis with magnitude equal
                to the rotation angle in radians.
```

Note that only the top left 3x3 rotation component of T is used, any
translation component in T is ignored.


## quaternion  - Construct quaternion.

```
Usage:  Q = quaternion(theta, axis)

Arguments: theta - angle of rotation
           axis  - 3x1 vector or array defining axis of rotation.
Returns:   Q     - a 4-vector quaternion in the form [w, xi, yj, zk]
```

## quaternionconjugate - Conjugate of a quaternion.

```
Usage: Qconj = quaternionconjugate(Q)

Argument: Q     - Quaternions in the form  Q = [Qw, Qi, Qj, Qk].
Returns:  Qconj - Conjugate.
```

## quaternionproduct - Computes product of two quaternions.

```
Usage: Q = quaternionproduct(A, B)

Arguments: A, B - Quaternions assumed to be 4-vectors in the
                  form  A = [Aw, Ai, Aj, Ak].
Returns:   Q    - Quaternion product.
```

## quaternionrotate - Rotates a 3D vector by a quaternion.

```
Usage:   vnew = quaternionrotate(Q, v)

Arguments: Q - a quaternion in the form [w, xi, yj, zk]
           v - a vector to rotate, either an inhomogeneous 3-vector or a
               homogeneous 4-vector.
Returns:   vnew - rotated vector.
```

## quaternion2matrix - Quaternion to a 4x4 homogeneous transformation matrix.

```
Usage:  T = quaternion2matrix(Q)

Argument: Q - a quaternion in the form [w, xi, yj, zk].
Returns:  T - 4x4 Homogeneous rotation matrix.
``` 

## matrix2quaternion - Homogeneous matrix to quaternion.

Converts 4x4 homogeneous rotation matrix to quaternion

```
Usage: Q = matrix2quaternion(T)

Argument:   T - 4x4 Homogeneous transformation matrix
Returns:    Q - 4-vector quaternion in the form [w, xi, yj, zk]
```

## vector2quaternion - Embeds 3-vector in a quaternion representation.

```
Usage: Q = vector2quaternion(v)

Argument:  v - 3-vector.
Returns:   Q - Quaternion given by [0; v]
```