# Test utililties.jl

println("testing utilites")

tol = 1e-2   # Fairly sloppy, needed for the derivative tests

# nonmaxsuppts - tested in test_cornerfeatures.jl

# derivative 3 5 7
# Generate an image of constant gradient
img = [c for r=1:50, c=1:50]
(gx,gy) = derivative3(img, ("x", "y"))

# Test interior of derivatives matrix matches expected values
@test all(abs.(gx[5:end-5,5:end-5] - 1) .< tol)
@test all(abs.(gy[5:end-5,5:end-5]) .< tol)

(gx,gy, gxx, gyy, gxy) = derivative5(img, ("x", "y", "xx", "yy", "xy"))
# Test interior of derivatives matrix matches expected values
@test all(abs.(gx[5:end-5,5:end-5] - 1) .< tol)
@test all(abs.(gy[5:end-5,5:end-5]) .< tol)
@test all(abs.(gxx[5:end-5,5:end-5]) .< tol)
@test all(abs.(gyy[5:end-5,5:end-5]) .< tol)
@test all(abs.(gxy[5:end-5,5:end-5]) .< tol)

(gx,gy, gxx, gyy, gxy) = derivative7(img, ("x", "y", "xx", "yy", "xy"))
# Test interior of derivatives matrix matches expected values
@test all(abs.(gx[5:end-5,5:end-5] - 1) .< tol)
@test all(abs.(gy[5:end-5,5:end-5]) .< tol)
@test all(abs.(gxx[5:end-5,5:end-5]) .< tol)
@test all(abs.(gyy[5:end-5,5:end-5]) .< tol)
@test all(abs.(gxy[5:end-5,5:end-5]) .< tol)


# gaussfilt
img = zeros(100,100)
img[50,50] = 1
simg = gaussfilt(img,2)
@test abs(sum(simg) - 1) < tol

# dilate1d erode1d
# Generate a 1D signal with a single spike
v = zeros(20)
v[10] = 1          
k = 5
vd = dilate1d(v,k)
@test all(vd[8:12] .== 1) && all(vd[1:7] .== 0) && all(vd[13:end] .== 0)

ve = erode1d(vd,k)   # Erode back and check we reconstruct v
@test all(ve .== v)

# imdilate imerode
img = zeros(100,100)
img[50,50] = 1
dimg = imdilate(img, "rect", [5, 3])
@test all(dimg[47:47,48:52] .== [0 0 0 0 0])
@test all(dimg[50:50,48:52] .== [0 1 1 1 0])
@test all(dimg[53:53,48:52] .== [0 0 0 0 0])
@test all(dimg[47:53,50] .== [0, 1, 1, 1, 1, 1, 0])

eimg = imerode(dimg, "rect", [5, 3])
@test all(eimg[49:51,49:51] .== [0 0 0; 0 1 0; 0 0 0])

# circularstruct

# floatxy

# matchbycorrelation

# grey2census

# imgnormalise
a = rand(5,3)
an = imgnormalise(a)
@test isapprox(minimum(an), 0)
@test isapprox(maximum(an), 1)

reqmean = 5
reqvar = 7
an = imgnormalise(a, reqmean, reqvar)
@test isapprox(mean(an), reqmean)
@test isapprox(var(an), reqvar)

# histtruncate
rows = 100;
cols = 200;
img = rand(rows, cols)
lHistCut = 2
uHistCut = 4
htimg = histtruncate(img, lHistCut, uHistCut)
minval = minimum(htimg)
maxval = maximum(htimg)
# Check that the number of saturated pixels at each extreme is
# approximately lHistCut and uHistCut. Note some integer rounding has
# to occur.
v = find(abs.(htimg-minval).<eps())
@test abs(length(v)/(rows*cols) * 100 - lHistCut) < 1

v = find(abs.(htimg-maxval).<eps())
@test abs(length(v)/(rows*cols) * 100 - uHistCut) < 1

