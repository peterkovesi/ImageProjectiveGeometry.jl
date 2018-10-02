# Demonstrate/test morphological functions 

using ImageProjectiveGeometry
using FileIO, PyPlot

#img = imread("lena.tif");
img = float(load("../sampleimages/lena.tif"))
figure(10), imshow(img);
set_cmap(ColorMap("gray"))
title("Input image")

img2 = rand(512,512) .> 0.999;
figure(11), imshow(img2);
title("Dots input image")

imgd = imerode(img, "RECT", [9, 15]);
figure(12), imshow(imgd);
title("Rectangular dilate")
keypause()

imgd = imerode(img, "RECT", 21);
figure(12), imshow(imgd);
title("Rectangular erode")
keypause()

imgd = imdilate(img, "OCT", 9);
figure(12), imshow(imgd);
title("Octagonal dilate")
keypause()

imgd = imerode(img, "OCT", 9);
figure(12), imshow(imgd);
title("Octagonal erode")
keypause()

img2d = imdilate(img2, "RECT", [19, 15]);
figure(12), imshow(img2d);
title("Dots: Rectangular dilate")
keypause()

img2d = imdilate(img2, "OCT", 9);
figure(12), imshow(img2d);
title("Dots: Octagonal dilate")
keypause()

se = circularstruct(7);
img2d = imdilate(img2, se);
figure(12), imshow(img2d);
title("Dots: Circular dilate")

nothing
