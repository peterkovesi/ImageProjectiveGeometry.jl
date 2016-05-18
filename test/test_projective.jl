# Test projective.jl

println("testing projective")

tol = 1e-10

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
skew = 0.01
X = 0; Y = 100; Z = 1000
Rc_w = rotx(pi-0.1)[1:3,1:3]     # camera points down(ish)

Cam = Camera(P=[X, Y, Z], Rc_w=Rc_w, fx=f, fy=f, ppx=ppx, ppy=ppy, 
                         rows=rows, cols=cols)

# ground points
gpt1 = [0, 0, 0]
gpt2 = [100, -200, 0]

(xy, visible) = cameraproject(Cam, [gpt1 gpt2])  # project to image

# Then project image points to ground plane to see if we reconstruct them
planeP = [0,0,0]
planeN = [0,0,1]
planept = imagept2plane(Cam, xy, planeP, planeN)

@test maximum(abs(planept - [gpt1 gpt2])) < tol
