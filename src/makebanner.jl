# Generate some kind of plot to use for the banner image

using ImageProjectiveGeometry
using PyPlot
  
figure(1)
clf()

# Generate 100 points
nGpts = 100
dx = 700; dy = 700; dz = 100

gx = (2*rand(nGpts)-1) * dx
gy = (2*rand(nGpts)-1) * dy
gz = (2*rand(nGpts)-1) * dz

plot3D(gx,gy,gz,"k*")
hold(true)

# Generate Cameras 
f = 4000          
rows = 3000;  cols = 4000
ppx = cols/2; ppy = rows/2;
Z = 1000

Ca = []
for X = -500:500:500
    for Y = -800:400:800
        Rc_w = rotx(pi)[1:3,1:3]     # camera points down
        push!(Ca, Camera(P=[X, Y, Z], Rc_w=Rc_w, fx=f, fy=f, ppx=ppx, ppy=ppy, 
                         rows=rows, cols=cols))
    end
end

plotcamera(Ca,150)
axis("equal")
