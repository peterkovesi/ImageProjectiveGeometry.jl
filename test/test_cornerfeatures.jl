# Test cornerfeatures.jl

println("testing cornerfeatures")

using ImageProjectiveGeometry
using Test

tol = 1;

# structuretensor, harris, noble, shi_tomasi, coherence
# hessianfeatures
# fastradial
# nonmaxsuppts

# Create an image with an isolated point
img = zeros(100,100)
r0 = 20;
c0 = 40;
img[r0,c0] = 10

(cimg, r, c) = harris(img, N=1,radius=5,subpixel=false)
@test abs(r[1]-r0) < tol && abs(c[1]-c0) < tol

(cimg, r, c) = noble(img, N=1,radius=5,subpixel=false)
@test abs(r[1]-r0) < tol && abs(c[1]-c0) < tol

(cimg, r, c) = shi_tomasi(img, N=1,radius=5,subpixel=false)
@test abs(r[1]-r0) < tol && abs(c[1]-c0) < tol

hdet = hessianfeatures(img, 1)
(r,c) = nonmaxsuppts(hdet, radius=5, N=1)
@test abs(r[1]-r0) < tol && abs(c[1]-c0) < tol

(S,So) = fastradial(img, [1, 3, 5])
(r,c) = nonmaxsuppts(S, radius=5, thresh=maximum(S)/2);
@test abs(r[1]-r0) < tol && abs(c[1]-c0) < tol

# Check nonmaxsuppts does subpixel precision. Create a 2x2 block of
# pixels in the image and locate the centre of that
img[r0:r0+1, c0:c0+1] .= 10
(cimg, r, c) = shi_tomasi(img, N=1,radius=5,subpixel=true)
tol = 0.1
@test abs(r[1]-(r0+0.5)) < tol && abs(c[1]-(c0+0.5)) < tol
