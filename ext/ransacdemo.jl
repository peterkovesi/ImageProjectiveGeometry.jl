#=---------------------------------------------------------------------

ransacdemo - Functions demonstrating the use of RANSAC to fit lines,
planes, fundamental matrices and homographies.

Copyright (c) 2016 Peter Kovesi
pk@peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK March     2016
   September 2018 Updates for v0.7/v1.0

---------------------------------------------------------------------=#

export fitlinedemo, fitplanedemo
export fitfunddemo, fithomogdemo

using ImageProjectiveGeometry, PyPlot, Printf, FileIO

#-----------------------------------------------------------------------
"""
fitlinedemo - Demonstrates RANSAC line fitting.

Function generates a noisy set of data points with outliers and uses
RANSAC to fit a line.

```
Usage: fitlinedemo(outliers, sigma, t, feedback=false)

Arguments:
              outliers - Fraction specifying how many points are to be
                         outliers.
              sigma    - Standard deviation of inlying points from the
                         true line.
              t        - Distance threshold to be used by the RANSAC
                         algorithm for deciding whether a point is an
                         inlier. 
              feedback - Optional Boolean flag turn on RANSAC feedback
                         information.

Try using:  fitlinedemo(0.3, 0.05, 0.05)
```
See also: ransacfitplane(), fitplane()
"""
function fitlinedemo(outliers, sigma, t::Real, feedback::Bool = false)
    # August 2006  testfitline created from testfitplane
    #              author: Felix Duvallet
    
    # Hard wire some constants - vary these as you wish
    npts = 100  # Number of 3D data points	
    
    # Define a line:
    #    Y = m*X
    #    Z = n*X + Y + b
    # This definition needs fixing, but it works for now
    m = 6
    n = -3
    b = -4
    
    outsigma = 30*sigma  # outlying points have a distribution that is
                         # 30 times as spread as the inlying points
    
    vpts = round(Int, (1-outliers)*npts)  # No of valid points
    opts = npts - vpts                    # No of outlying points
    
    # Generate npts points in the line
    X = rand(Float64, 1,npts)
    
    Y = m*X
    Z = n*X .+ Y .+ b
    Z = zeros(size(Y))
    
    XYZ =  [X
    	    Y
    	    Z]

    # Add uniform noise of +/-sigma
    XYZ = XYZ .+ (2*rand(Float64, size(XYZ)) .- 1)*sigma
    
    # Generate opts random outliers
    n = size(XYZ,2)
    ind = randperm(n)[1:opts]  # get a random set of point indices of length opts 
      
    # Add uniform noise of outsigma to the points chosen to be outliers.  
    XYZ[:,ind] = XYZ[:,ind] .+ sign.(rand(Float64, 3,opts).-0.5).*(rand(Float64, 3,opts).+1)*outsigma    
    
    # Perform RANSAC fitting of the line
    (V, P, inliers) = ransacfitline(XYZ, t, feedback)
    
    if feedback
        @printf("Number of Inliers: %d\n", length(inliers))
    end

    # Plot the inlier points blue, with the outlier points in red.
    # Display the cloud of outlier points
    figure(1); clf()
    outliers = setdiff(1:npts, inliers)
    plot3D(vec(XYZ[1,outliers]), vec(XYZ[2, outliers]), vec(XYZ[3, outliers]), "r*")
#    hold(true)

    # Plot the inliers as blue points
    plot3D(vec(XYZ[1,inliers]), vec(XYZ[2, inliers]), vec(XYZ[3, inliers]), "b*")

    # Display the line formed by the 2 points that gave the
    # line of maximum consensus as a green line
    plot3D(vec(P[1,:]), vec(P[2,:]), vec(P[3,:]), "g-", linewidth=4)
    
    #Display the line formed by the covariance fitting in magenta
    plot3D(vec(V[1,:]), vec(V[2, :]), vec(V[3,:]), "m-", linewidth=4)
    box("on"), grid("on")
#    hold(false)
    
end    
#-----------------------------------------------------------------------
"""
fitplanedemo - Demonstrates RANSAC plane fitting

Function generates a noisy set of data points with outliers and uses
RANSAC to fit a plane.

```
Usage: fitplanedemo(outliers, sigma, t, feedback)

Arguments:
              outliers - Fraction specifying how many points are to be
                         outliers.
              sigma    - Standard deviation of inlying points from the
                         true plane.
              t        - Distance threshold to be used by the RANSAC
                         algorithm for deciding whether a point is an
                         inlier. 
              feedback - Optional flag 0 or 1 to turn on RANSAC feedback
                         information.

Try using:  fitplanedemo(0.3, 0.02, 0.05)
```
See also: ransacfitplane(), fitplane()
"""
function fitplanedemo(outliers, sigma, t, feedback::Bool = false)

    # Hard wire some constants - vary these as you wish
    npts = 100  # Number of 3D data points	
    
    # Define a plane  ax + by + cz + d = 0
    a = 10; b = -3; c = 5; d = 1
    
    B = [a, b, c, d]
    
    outsigma = 30*sigma  # outlying points have a distribution that is
                         # 30 times as spread as the inlying points
    
    vpts = round(Int, (1-outliers)*npts)  # No of valid points
    opts = npts - vpts                    # No of outlying points
    
    # Generate npts points in the plane
    X = rand(Float64, npts)'
    Y = rand(Float64, npts)'
    Z = (-a*X .- b*Y .-d)/c
    
    XYZ =  [X
	    Y
	    Z]
    
    # Add uniform noise of +/-sigma
    XYZ = XYZ .+ (2*rand(Float64, size(XYZ)) .- 1)*sigma
    
    # Generate opts random outliers
    
    n = size(XYZ,2)
    ind = randperm(n)[1:opts]  # get a random set of point indices of length opts
    
    # Add uniform noise of outsigma to the points chosen to be outliers.
    #  XYZ(:,ind) = XYZ(:,ind)  + (2*rand(3,opts)-1)*outsigma
    XYZ[:,ind] = XYZ[:,ind] .+ sign.(rand(Float64, 3,opts).-0.5).*(rand(Float64, 3,opts).+1)*outsigma    
    
    # Display the cloud of points (vec() needed for 0.4)
    figure(1); clf(); plot3D(vec(XYZ[1,:]),vec(XYZ[2,:]),vec(XYZ[3,:]), "r*")
    box("on")
    
    # Perform RANSAC fitting of the plane
    (Bfitted, P, inliers) = ransacfitplane(XYZ, t, feedback)

    Bfitted = Bfitted/Bfitted[4]
    @printf("Original plane coefficients:  ")
    @printf("%8.3f %8.3f %8.3f %8.3f \n",B[1], B[2], B[3], B[4])
    @printf("Fitted plane coefficients:    ")
    @printf("%8.3f %8.3f %8.3f %8.3f \n",Bfitted[1], Bfitted[2], Bfitted[3], Bfitted[4])

    # Display the triangular patch formed by the 3 points that gave the
    # plane of maximum consensus
#    hold(true)
    pts = [P P[:,1]]'
    plot3D(pts[:,1], pts[:,2], pts[:,3], "k-")
#    hold(false)
    
    @printf("\nRotate image so that the triangular patch is seen edge on\n")
    @printf("These are the points that form the plane of max consensus.\n\n")

end

#-----------------------------------------------------------------------
"""
fitfunddemo - Example of fundamental matrix computation

Demonstration of feature matching via simple correlation, and then using
RANSAC to estimate the fundamental matrix and at the same time identify
(mostly) inlying matches
```
Usage:  fitfunddemo           - Demonstrates fundamental matrix calculation
                                on two default images.
        funddemo(img1,img2)   - Computes fundamental matrix on two supplied
                                images.
```
Edit code as necessary to tweak parameters

See also: ransacfitfundmatrix(), fundmatrix()
"""
function fitfunddemo(img1=[], img2=[])
    
    if isempty(img1)
	#img1 = PyPlot.imread("img02.jpg")
	#img2 = PyPlot.imread("img03.jpg")
        img1 = float(load("img02.jpg"))
	img2 = float(load("img03.jpg"))
    end

    (rows1,cols1) = size(img1)
    (rows2,cols2) = size(img2)

    nonmaxrad = 3  # Non-maximal suppression radius
    dmax = 100     # Maximum search distance for matching
    w = 11         # Window size for correlation matching
    
    # Find Shi Tomasi corners in image1 and image2
    (cim1, r1, c1) = shi_tomasi(img1, 1, radius=nonmaxrad, N=200, img=img1, fig=1)
    (cim2, r2, c2) = shi_tomasi(img2, 1, radius=nonmaxrad, N=200, img=img2, fig=2)
    figure(1); axis("off"); 
    figure(2); axis("off"); 
    keypause()

    @printf("Matching features...\n")
    (m1,m2) = matchbycorrelation(copy(img1), [r1';c1'], copy(img2), [r2';c2'], w, dmax)
    
    # Display putative matches
    figure(3); clf(); imshow(img1); axis("off")  # hold(true)
    nMatch = size(m1,2)
    for n = 1:nMatch
	plot([m1[2,n], m2[2,n]], [m1[1,n], m2[1,n]],"r-")
    end
    title("Putative matches")
#    hold(false)

    # Assemble homogeneous feature coordinates for fitting of the
    # fundamental matrix, note that [x,y] corresponds to [col, row]
    x1 = [m1[2:2,:]; m1[1:1,:]; ones(1,nMatch)]
    x2 = [m2[2:2,:]; m2[1:1,:]; ones(1,nMatch)]    

    t = .002  # Distance threshold for deciding outliers
    
    # Change the commenting on the lines below to switch between the 8
    # point fundamental matrix solution or the affine fundamental
    # matrix solution.
    (F, inliers) = ransacfitfundmatrix(x1, x2, t, true)
#   (F, inliers) = ransacfitaffinefundmatrix(x1, x2, t, true)    

    @printf("Number of inliers was %d (%d%%) \n", 
	    length(inliers),round(Int, 100*length(inliers)/nMatch))
    @printf("Number of putative matches was %d \n", nMatch)
    
    # Display both images overlayed with inlying matched feature points
    figure(4); clf();  imshow(img1);  axis("off") # hold(true)
    plot(m1[2,inliers], m1[1,inliers],"r+")
    plot(m2[2,inliers], m2[1,inliers],"g+")    
    for n = inliers
	plot([m1[2,n], m2[2,n]], [m1[1,n], m2[1,n]], "b-")
    end
    title("Inlying matches")
#    hold(false)

    @printf("Step through each epipolar line [y/n]?\n")
    response = readline()
    if response[1] == 'n'
	return
    end    

    # Step through each matched pair of points and display the
    # corresponding epipolar lines on the two images.
    l2 = F*x1    # Epipolar lines in image2
    l1 = F'*x2   # Epipolar lines in image1
    
    # Solve for epipoles
    (U,D,V) = svd(F)
    e1 = hnormalise(V[:,3])
    e2 = hnormalise(U[:,3])

    for n = inliers
	figure(1); clf(); imshow(img1)
        axis([1,cols1,rows1,1])
#        hold(true)
        plot(x1[1,n],x1[2,n],"r+",markersize=8)
	hline(l1[:,n]);
        plot(e1[1], e1[2], "g*"); 
#        hold(false)

	figure(2); clf(); imshow(img2); 
        axis([1,cols2,rows2,1])
 #       hold(true); 
        plot(x2[1,n],x2[2,n],"r+",markersize=8)
	hline(l2[:,n]); 
        plot(e2[1], e2[2], "g*")
#        hold(false)
	keypause()
    end

    @printf("                                         \n")

end    

#-----------------------------------------------------------------------
"""
fithomogdemo - Example of finding a homography

Demonstration of feature matching via simple correlation, and then using
RANSAC to estimate the homography between two images and at the same time
identify (mostly) inlying matches
```
Usage:  fithomogdemo            - Demonstrates homography calculation on two 
                                  default images
        fithomogdemo(img1,img2) - Computes homography on two supplied images

Edit code as necessary to tweak parameters
```
See also: ransacfithomography(), homography2d()
"""
function fithomogdemo(img1=[], img2=[])

    if isempty(img1)
	#img1 = PyPlot.imread("boats.jpg")
	#img2 = PyPlot.imread("boatsrot.jpg")
        img1 = float(load("boats.jpg"))
	img2 = float(load("boatsrot.jpg"))
    end

    nonmaxrad = 3  # Non-maximal suppression radius
    dmax = 100     # Maximum search distance for matching
    w = 11         # Window size for correlation matching

    # Find Shi Tomasi corners in image1 and image2
    (cim1, r1, c1) = shi_tomasi(img1, 1, radius=nonmaxrad, N=100, img=img1, fig=1)
    (cim2, r2, c2) = shi_tomasi(img2, 1, radius=nonmaxrad, N=100, img=img2, fig=2)
    keypause()
    @printf("Matching features...\n")
    (m1,m2) = matchbycorrelation(img1, [r1';c1'], img2, [r2';c2'], w, dmax)

    # Display putative matches
    figure(3); clf();  imshow(img1);  # hold(true)
    nMatch = size(m1,2)
    for n = 1:nMatch
	plot([m1[2,n], m2[2,n]], [m1[1,n], m2[1,n]],"r-")
    end
    title("Putative matches")
#    hold(false)

    # Assemble homogeneous feature coordinates for fitting of the
    # fundamental matrix, note that [x,y] corresponds to [col, row]
    x1 = [m1[2:2,:]; m1[1:1,:]; ones(1,nMatch)]
    x2 = [m2[2:2,:]; m2[1:1,:]; ones(1,nMatch)]    

    t = .001  # Distance threshold for deciding outliers
    (H, inliers) = ransacfithomography(x1, x2, t)

    @printf("Number of inliers was %d (%d%%) \n", 
	    length(inliers),round(Int, 100*length(inliers)/nMatch))
    @printf("Number of putative matches was %d \n", nMatch) 
    
    # Display both images overlayed with inlying matched feature points
    figure(4); clf();  imshow(img1);  # hold(true)
    plot(m1[2,inliers], m1[1,inliers],"r+")
    plot(m2[2,inliers], m2[1,inliers],"g+")    
    for n = inliers
	plot([m1[2,n], m2[2,n]], [m1[1,n], m2[1,n]], "b-")
    end

    title("Inlying matches")
#    hold(false)
end

