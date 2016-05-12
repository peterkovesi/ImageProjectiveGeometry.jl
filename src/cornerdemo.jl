
# Demonstrate/test functions in cornerfeatures

using ImageProjectiveGeometry
using PyPlot
using PerceptualColourMaps


img = PyPlot.imread("lena.tif");
figure(10), PyPlot.imshow(img);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))
PyPlot.title("Input image")

println("Harris corner detector")
(cimg, r, c) = ImageProjectiveGeometry.harris(img, 1,k= 0.04, radius=2, N=100, subpixel=false, img=img)
cimg[cimg .< 0] = 0;   # eliminate -ve values
figure(2), PyPlot.imshow(cimg.^0.25);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))
ImageProjectiveGeometry.keypause()

println("Noble corner detector")
(cimg, r, c) = ImageProjectiveGeometry.noble(img, 1, radius=2, N=100, subpixel=false, img=img);
figure(2), PyPlot.imshow(cimg.^0.25);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))
ImageProjectiveGeometry.keypause()

println("Shi-Tomasi corner detector")
(cimg, r, c) = ImageProjectiveGeometry.shi_tomasi(img, 1, radius=2, N=100, subpixel=false, img=img);
figure(2), PyPlot.imshow(cimg.^0.25);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))
ImageProjectiveGeometry.keypause()

println("-ve Hessian features")
hdet = ImageProjectiveGeometry.hessianfeatures(img, 1)
(r, c) = ImageProjectiveGeometry.nonmaxsuppts(-hdet, N=100, img=img)
figure(2), PyPlot.imshow(-hdet);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))
ImageProjectiveGeometry.keypause()

println("+ve Hessian features")
hdet = ImageProjectiveGeometry.hessianfeatures(img, 1)
(r, c) = ImageProjectiveGeometry.nonmaxsuppts(+hdet, N=100, img=img)
figure(2), PyPlot.imshow(hdet);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))
ImageProjectiveGeometry.keypause()


println("Coherence")
cimg = ImageProjectiveGeometry.coherence(img, 2)
(r, c) = ImageProjectiveGeometry.nonmaxsuppts(cimg, N=10, img=img)
figure(2), PyPlot.imshow(cimg);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))
ImageProjectiveGeometry.keypause()

println("Fast radial")
(S, So) = ImageProjectiveGeometry.fastradial(img, [1, 3, 5, 7])
(r, c) = ImageProjectiveGeometry.nonmaxsuppts(-S, N=10, img=img)
figure(2), PyPlot.imshow(-S);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))

figure(3), PyPlot.imshow(-So);
PyPlot.set_cmap(PyPlot.ColorMap("gray"))


