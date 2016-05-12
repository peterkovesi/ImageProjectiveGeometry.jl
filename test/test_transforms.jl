# Test transforms.jl


println("testing transforms")

tol = 1e-10

# rotx, roty, rotz, invrpy
roll = .3
pitch = -1
yaw = .8
T = rotz(yaw)*roty(pitch)*rotx(roll)

(rpy1, rpy2) = invrpy(T)

# check 1st solution matches roll pitch and yaw
@test abs(rpy1(1)-yaw) < tol
@test abs(rpy1(2)-pitch) < tol
@test abs(rpy1(3)-roll) < tol

# check end solution produces the same transform
T2 = rotz(rpy2(1))*roty(rpy2(2))*rotx(rpy2(3))
@test maximum(abs(T2-T)) < tol

