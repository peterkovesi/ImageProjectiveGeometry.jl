# Test projective.jl
using Test, Statistics, LinearAlgebra

println("testing projective")

tol = 1e-8

# Camera
# cameraproject, camera2projmatrix, imagept2plane, decomposecamera, rq3
# makehomogeneous, makeinhomogeneous, hnormalise
# homography1d, homography2d, normalise1dpts, normalise2dpts 
# skew, hcross
# fundmatrix, affinefundmatrix, fundfromcameras
# idealimagepts, solvestereopt
# hline, plotcamera

# Generate Camera with some odd paraeters
f = 4000          
rows = 3000;  cols = 4000
ppx = cols*0.3; ppy = rows*0.6;
skew = 0.3
k1 = 0.0
k3 = 0.0
X = 0; Y = 100; Z = 1000
Rc_w = rotx(pi-0.1)[1:3,1:3]     # camera points down(ish)

Cam = Camera(P=[X, Y, Z], Rc_w=Rc_w, fx=f, fy=f, ppx=ppx, ppy=ppy, 
          skew=skew, k1=k1, k3=k3, rows=rows, cols=cols)

# ground points
gpt1 = [0.0, 0.0, 0.0]
gpt2 = [100.0, -200.0, 0.0]

(xy, visible) = cameraproject(Cam, [gpt1 gpt2], computevisibility=true)  # project to image

# Then project image points to ground plane to see if we reconstruct them
# ** Note this does not seem to work with lens distortion
planeP = [0.0,0.0,0.0]
planeN = [0.0,0.0,1.0]
planept = imagept2plane(Cam, xy, planeP, planeN)
@test maximum(abs.(planept - [gpt1 gpt2])) < f*tol

# Convert Cam structure to projection matrix and decompose to see if we get the same parameters back
P = camera2projmatrix(Cam)

(K, Rc_w_d, Pc_d, pp_d, pv) = decomposecamera(P)

@test abs(K[1,1] - f) < f*tol
@test abs(K[2,2] - f) < f*tol
@test abs(K[1,2] - skew) < f*tol
@test abs(K[1,3] - ppx) < rows*tol
@test abs(K[2,3] - ppy) < rows*tol

# makehomogeneous makeinhomoheneous hnormalise
x = rand(2,4)
hx = makehomogeneous(x)
x2 = makeinhomogeneous(3*hx)
@test maximum(abs.(x-x2)) < tol

# normalise1dpts
hx = hnormalise(100*rand(2,10))
(xn, T) = normalise1dpts(hx)
@test abs(mean(xn[1,:])) < tol  # mean 0
@test abs(mean(abs.(xn[1,:])) - 1) < tol  # mean deviation 1

#  normalise2dpts 
hx = hnormalise(100*rand(3,10))
(xn, T) = normalise2dpts(hx)
@test abs(mean(xn[1,:])) < tol  # mean 0
@test abs(mean(xn[2,:])) < tol

# mean distance from origin is sqrt(2)
@test abs(mean(sqrt.(sum(xn[1:2,:].^2, dims=1))) - sqrt(2)) < tol

# check transform works
@test maximum(abs.(xn - T*hx)) < tol

# skew, hcross
a = rand(3)
b = rand(3)
aXb = ImageProjectiveGeometry.skew(a)*b
@test maximum(abs.(aXb .- cross(a,b))) < tol

@test maximum(abs.(hnormalise(aXb) .- hcross(a,b))) < tol

# homography1d, homography2d

# Define a 1D homography
a = 1; b = 2; c = 3
H = [a b
     c 1.]

x1 = [1 2 3 4 5
      1 1 1 1 1.]
x2 = hnormalise(H*x1)

Hfit = homography1d(x1,x2)
Hfit = Hfit/Hfit[2,2]
@test maximum(abs.(Hfit - H)) < tol

# Define a 2D homography
a = 1; b = 2; c = 3; d = -1; e = -2; f = -3; g = 2; h = -1
H = [a b c
     d e f
     g h 1.]

x1 = makehomogeneous(rand(2, 6))
x2 = hnormalise(H*x1)
Hfit = homography2d(x1,x2)
Hfit = Hfit/Hfit[3,3]
@test maximum(abs.(Hfit - H)) < tol


# fundmatrix, affinefundmatrix, fundfromcameras

# Generate a set of 3D points
pts = rand(3,12)

# Define 2 Cameras
f = 4000
rows = 3000;  cols = 4000
ppx = cols/2; ppy = rows/2;
X = 0; Y = 0; Z = 10
Rc_w = (rotx(pi)*roty(0.1))[1:3,1:3]    

Cam1 = Camera(P=[X, Y, Z], Rc_w=Rc_w, fx=f, fy=f, ppx=ppx, ppy=ppy, 
                         rows=rows, cols=cols)
Rc_w = (rotx(pi-.1)*roty(-0.1))[1:3,1:3]    
Cam2 = Camera(P=[X+2, Y, Z], Rc_w=Rc_w, fx=f, fy=f, ppx=ppx, ppy=ppy, 
                         rows=rows, cols=cols)

xy1 = cameraproject(Cam1, pts)  # Project points into images
xy2 = cameraproject(Cam2, pts)

# Form fundamental matrix from corresponding image points
F1 = fundmatrix(makehomogeneous(xy1), makehomogeneous(xy2))

# Form fundamental matrix from the camera projection matrices 
F2 = fundfromcameras(Cam1, Cam2)

F1 = F1/F1[3,3]  # Adjust matrices to the same scale
F2 = F2/F2[3,3]

@test maximum(abs.(F1 - F2)) < tol

# solvestereopt
res = []
for i in 1:12
          push!(res, maximum(makeinhomogeneous(solvestereopt([xy1[:,i] xy2[:,i]], [Cam1, Cam2])) .- pts[:,i] ) < tol)
end
@test all(res)

# imagept2ray
# construct rays for projected image coordinates and compare with original points
res = []
for i in 1:12
     dist = norm(pts[:,i] - Cam1.P)
     ray = normalize(imagept2ray(Cam1, xy1[1,i], xy1[2,i]))
     target = dist * ray + Cam1.P
     push!(res, norm(target - pts[:,i]) < tol)
end
@test all(res)

