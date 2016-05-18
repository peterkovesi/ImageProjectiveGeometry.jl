# Test transforms.jl

println("testing transforms")

tol = 1e-10

# rotx, roty, rotz, invrpy
roll = 0.3
pitch = -1.0
yaw = 0.8
T = rotz(yaw)*roty(pitch)*rotx(roll)

(rpy1, rpy2) = invrpy(T)

# Check 1st solution matches roll pitch and yaw
@test abs(rpy1[1]-yaw) < tol
@test abs(rpy1[2]-pitch) < tol
@test abs(rpy1[3]-roll) < tol

# Check 2nd solution produces the same transform
T2 = rotz(rpy2[1])*roty(rpy2[2])*rotx(rpy2[3])
@test maximum(abs(T2-T)) < tol

# inveuler
phi = 0.5
theta = -0.2
psi = 0.1
T = rotz(phi) * roty(theta) * rotz(psi)     
(euler1, euler2) = inveuler(T)

# Check both solutions produce the same transform
T1 = rotz(euler1[1])*roty(euler1[2])*rotz(euler1[3])
@test maximum(abs(T1-T)) < tol

T2 = rotz(euler2[1])*roty(euler2[2])*rotz(euler2[3])
@test maximum(abs(T2-T)) < tol

# invht
v = [1,-2,3]
T = rotz(yaw)*roty(pitch)*rotx(roll) * trans(v)

@test maximum(abs(invht(T)*T - eye(4))) < tol
@test maximum(abs(trans(v) - trans(v[1],v[2],v[3]))) < tol

# angleaxis, angleaxis2matrix, matrix2angleaxis
# matrix2quaternion, quaternion2matrix
T1 = rotz(yaw)*roty(pitch)*rotx(roll) 
ax = matrix2angleaxis(T)
T2 = angleaxis2matrix(ax)
@test maximum(abs(T2-T1)) < tol

theta = 0.35
axis = T1[1:3,1]
ax = angleaxis(theta, axis)
Q = quaternion(theta,axis)

T1 = quaternion2matrix(Q)
T2 = angleaxis2matrix(ax)
@test maximum(abs(T2-T1)) < tol

Q2 = matrix2quaternion(T1)
@test maximum(abs(Q-Q2)) < tol

# angleaxisrotate, quaternionrotate
ax = angleaxis(theta, [1,0,0])
Q = quaternion(theta, [1,0,0])
v = [3,-2,7,1]
vnew1 = quaternionrotate(Q,v)
vnew2 = angleaxisrotate(ax,v)
vnew3 = rotx(theta)*v

@test maximum(abs(vnew3-vnew1)) < tol
@test maximum(abs(vnew3-vnew2)) < tol

# quaternionconjugate, quaternionproduct, 
# vector2quaternion
