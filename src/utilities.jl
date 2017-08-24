#=--------------------------------------------------------------------

utilities - Various utility functions for computer viison


Copyright (c) 2015-2017 Peter Kovesi
peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK April 2016 Initial version
   April 2017 Updates

---------------------------------------------------------------------=#

export derivative3, derivative5, derivative7
export nonmaxsuppts
export dilate1d, erode1d, imdilate, imerode, circularstruct
export gaussfilt
export imgnormalise, imgnormalize, histtruncate
export floatyx
export matchbycorrelation
export medfilt2, stdfilt2
export grey2census, grey2lbp, grey2lbp!
export keypause

import Images, ImageFiltering, PyPlot

#----------------------------------------------------------------------
"""
derivative3 - 3-Tap discrete derivative filters

This function computes 1st derivatives of an image using the 3-tap
coefficients given by Farid and Simoncelli.  The results are significantly
more accurate than simple gradient functions on edges that are at angles
other than vertical or horizontal. This, in turn, improves gradient orientation
estimation enormously.  If you are after extreme accuracy try using derivative7().
```
Usage:  (gx, gy) = derivative3(img, derivative specifiers)

Arguments:
                     img - Image to compute derivatives from. ::Array{T,2}
   derivative specifiers - A tuple of character strings
                           that can be any of "x" or "y"
                           These can be in any order, the order of the
                           computed output arguments will match the order
                           of the derivative specifier strings.
Returns:
 Function returns requested derivatives which can be:
    gx, gy   - 1st derivative in x and y

 Example:
   Compute 1st derivatives in x and y
   (gx, gy) = derivative3(img, ("x", "y"))
```                                           
See also: derivative5(), derivative7()
"""
# Reference: Hany Farid and Eero Simoncelli "Differentiation of Discrete
# Multi-Dimensional Signals" IEEE Trans. Image Processing. 13(4): 496-508 (2004)

# April 2010

function derivative3{T<:Real}(img::Array{T,2}, spec)

    # Convert to float if needed
    if isa(img[1,1], Integer)
        fimg = float(img) 
    else
        fimg = img
    end

    # 3-Tap filters
    p = [0.229879,  0.540242,  0.229879]
    
    # ?? Coefficients seem too small added the division...??
    d1 =[0.425287,  0.000000, -0.425287]/0.8506
    
    # Compute derivatives.  Note that in the 1st call below conv2
    # function performs a 1D convolution down the columns using p then a 1D
    # convolution along the rows using d1. etc etc.
    G = Array{Array}(length(spec))
    
    for n = 1:length(spec)
      if spec[n] == "x"
          G[n] = conv2(p, d1, fimg)[2:end-1, 2:end-1]
      elseif spec[n] == "y"
          G[n] = conv2(d1, p, fimg)[2:end-1, 2:end-1]
      else
          badspec = spec[n]
          error(" $badspec is an unrecognized derivative option")
      end
    end
    
    return tuple(G...)
end

#----------------------------------------------------------------------
"""
derivative5 - 5-Tap 1st and 2nd discrete derivatives

This function computes 1st and 2nd derivatives of an image using the 5-tap
coefficients given by Farid and Simoncelli.  The results are significantly
more accurate than simple gradient functions on edges that are at angles
other than vertical or horizontal. This, in turn, improves gradient orientation
estimation enormously.  If you are after extreme accuracy try using derivative7().
```
Usage:  (gx, gy, gxx, gyy, gxy) = derivative5(img, derivative_specifiers)

Arguments:
                     img - Image to compute derivatives from. ::Array{T,2}
   derivative specifiers - A tuple of character strings
                           that can be any of "x", "y", "xx", "yy" or "xy"
                           These can be in any order, the order of the
                           computed output arguments will match the order
                           of the derivative specifier strings.
Returns:
 Function returns requested derivatives which can be:
    gx, gy   - 1st derivative in x and y
    gxx, gyy - 2nd derivative in x and y
    gxy      - 1st derivative in y of 1st derivative in x

 Examples:
   Just compute 1st derivatives in x and y
   (gx, gy) = derivative5(img, ("x", "y"))
                                           
   Compute 2nd derivative in x, 1st derivative in y and 2nd derivative in y
   (gxx, gy, gyy) = derivative5(img, ("xx", "y", "yy"))
```
See also: derivative7()
"""
#=
Reference: Hany Farid and Eero Simoncelli "Differentiation of Discrete
Multi-Dimensional Signals" IEEE Trans. Image Processing. 13(4): 496-508 (2004)

April 2010
September 2015 Ported to Julia
=#

function derivative5{T<:Real}(img::Array{T,2}, spec)

    # Convert to float if needed
    if isa(img[1,1], Integer)
        fimg = float(img) 
    else
        fimg = img
    end

    # Check if we are just computing 1st derivatives.  If so use the
    # interpolant and derivative filters optimized for 1st derivatives, else
    # use 2nd derivative filters and interpolant coefficients.
    # Detection is done by seeing if any of the derivative specifier
    # arguments is longer than 1 char, this implies 2nd derivative needed.
    secondDeriv = false    
    for n = 1:length(spec)
        if length(spec[n]) > 1
            secondDeriv = true
            break
        end
    end
    
    if !secondDeriv
        # 5 tap 1st derivative cofficients.  These are optimal if you are just
        # seeking the 1st deriavtives
        p = [0.037659,  0.249153,  0.426375,  0.249153,  0.037659]
        d1 =[0.109604,  0.276691,  0.000000, -0.276691, -0.109604]
    else         
        # 5-tap 2nd derivative coefficients. The associated 1st derivative
        # coefficients are not quite as optimal as the ones above but are
        # consistent with the 2nd derivative interpolator p and thus are
        # appropriate to use if you are after both 1st and 2nd derivatives.
        p  = [0.030320,  0.249724,  0.439911,  0.249724,  0.030320]
        d1 = [0.104550,  0.292315,  0.000000, -0.292315, -0.104550]
        d2 = [0.232905,  0.002668, -0.471147,  0.002668,  0.232905]
    end

    # Compute derivatives.  Note that in the 1st call below conv2
    # function performs a 1D convolution down the columns using p then
    # a 1D convolution along the rows using d1. etc etc.
    # Note that the result has to be trimmed becasue conv2 expands the 
    # result due to size of the filter
    gx = false
    gxn = 0
    G = Array{Array}(length(spec))

    for n = 1:length(spec)
      if spec[n] == "x"
          G[n] = conv2(p, d1, fimg)[3:end-2, 3:end-2]
          gx = true   # Record that gx is available for gxy if needed
          gxn = n
      elseif spec[n] == "y"
          G[n] = conv2(d1, p, fimg)[3:end-2, 3:end-2]
      elseif spec[n] == "xx"
          G[n] = conv2(p, d2, fimg)[3:end-2, 3:end-2]
      elseif spec[n] == "yy"
          G[n] = conv2(d2, p, fimg)[3:end-2, 3:end-2]
      elseif spec[n] == "xy" || spec[n] == "yx"
          if gx
              G[n] = conv2(d1, p, G[gxn])[3:end-2, 3:end-2]
          else
              gx = conv2(p, d1, fimg)[3:end-2, 3:end-2]
              G[n] = conv2(d1, p, gx)[3:end-2, 3:end-2]
          end
      else
          badspec = spec[n]
          error(" $badspec is an unrecognized derivative option")
      end
    end
   
    return tuple(G...)
end

#----------------------------------------------------------------------
"""
derivative7 - 7-Tap 1st and 2nd discrete derivatives

This function computes 1st and 2nd derivatives of an image using the 7-tap
coefficients given by Farid and Simoncelli.  The results are significantly
more accurate than simple gradient functions on edges that are at angles
other than vertical or horizontal. This, in turn, improves gradient orientation
estimation enormously.
```
Usage:  (gx, gy, gxx, gyy, gxy) = derivative7(img, derivative specifiers)

Arguments:
                      im - Image to compute derivatives from. ::Array{T,2}
   derivative specifiers - A tuple of character strings
                           that can be any of "x", "y", "xx", "yy" or "xy"
                           These can be in any order, the order of the
                           computed output arguments will match the order
                           of the derivative specifier strings.
Returns:
 Function returns requested derivatives which can be:
    gx, gy   - 1st derivative in x and y
    gxx, gyy - 2nd derivative in x and y
    gxy      - 1st derivative in y of 1st derivative in x

 Examples:
   Just compute 1st derivatives in x and y
   (gx, gy) = derivative7(img, ("x", "y"))  
                                           
   Compute 2nd derivative in x, 1st derivative in y and 2nd derivative in y
   (gxx, gy, gyy) = derivative7(img, ("xx", "y", "yy"))
```
See also: derivative5()
"""
# Reference: Hany Farid and Eero Simoncelli "Differentiation of Discrete
# Multi-Dimensional Signals" IEEE Trans. Image Processing. 13(4): 496-508 (2004)

# April 2010

function derivative7{T<:Real}(img::Array{T,2}, spec)

    # Convert to float if needed
    if isa(img[1,1], Integer)
        fimg = float(img) 
    else
        fimg = img
    end

    # 7-tap interpolant and 1st and 2nd derivative coefficients
    p  = [ 0.004711,  0.069321,  0.245410,  0.361117,  0.245410,  0.069321,  0.004711]
    d1 = [ 0.018708,  0.125376,  0.193091,  0.000000, -0.193091, -0.125376, -0.018708]
    d2 = [ 0.055336,  0.137778, -0.056554, -0.273118, -0.056554,  0.137778,  0.055336]
    
    # Compute derivatives.  Note that in the 1st call below conv2
    # function performs a 1D convolution down the columns using p then a 1D
    # convolution along the rows using d1. etc etc.
    # Note that the result has to be trimmed becasue conv2 expands the 
    # result due to size of the filter
    gx = false
    gxn = 0
    G = Array{Array}(length(spec))

    for n = 1:length(spec)
      if spec[n] == "x"
          G[n] = conv2(p, d1, fimg)[4:end-3, 4:end-3]
          gx = true    # Record that gx is available for gxy if needed
          gxn = n
      elseif spec[n] == "y"
          G[n] = conv2(d1, p, fimg)[4:end-3, 4:end-3]
      elseif spec[n] == "xx"
          G[n] = conv2(p, d2, fimg)[4:end-3, 4:end-3]
      elseif spec[n] == "yy"
          G[n] = conv2(d2, p, fimg)[4:end-3, 4:end-3]
      elseif spec[n] == "xy" || spec[n] == "yx"
          if gx
              G[n] = conv2(d1, p, G[gxn])[4:end-3, 4:end-3]
          else
              gx = conv2(p, d1, fimg)[4:end-3, 4:end-3]
              G[n] = conv2(d1, p, gx)[4:end-3, 4:end-3]
          end
      else
          badspec = spec[n]
          error(" $badspec is an unrecognized derivative option")
      end
    end
    
    return tuple(G...)
end

#-------------------------------------------------------------------
"""
nonmaxsuppts - Non-maximal suppression for features/corners

Non maxima suppression and thresholding for points generated by a feature
or corner detector.
```
Usage:   (r,c) = nonmaxsuppts(cimg; radius=1, thresh=0, N=Inf, subpix=false, img=[], fig=nothing)
                                                            
Argument:
           cimg   - corner strength image. ::Array{T,2}

Keyword Arguments:
           radius - radius of region considered in non-maximal
                    suppression. Typical values to use might
                    be 1-3 pixels. Default (and minimum value) is 1.
           thresh - Threshold, only features with value greater than
                    threshold are returned. Default is 0.
                N - Maximum number of features to return.  In this case the
                    N strongest features with value above 'thresh' are
                    returned. Default is Inf.
         subpixel - If set to true features are localised to subpixel
                    precision. Default is false.
             img  - Optional image data.  If an image is supplied the
                    detected corners are overlayed on this image. This can
                    be useful for parameter tuning. Default is [].
              fig - Optional figure number to display image and corner points. 
                    If not specified a new figure is generated.

Returns:
           r,c    - row and column coordinates of corner points.
```

Example of use:
```
 > hes = hessianfeatures(img, 1)   # Compute Hessian feature image in image img
```
Find the 1000 strongest features to subpixel precision using a non-maximal
suppression radius of 2 and overlay the detected corners on the origin image.
```
 > (r,c) = nonmaxsuppts(abs(hes), radius=2, N=1000, img=img, subpixel=true)
```
Note: An issue with integer valued images is that if there are multiple pixels
all with the same value within distance 2*radius of each other then they will
all be marked as local maxima. 

See also: harris(), noble(), shi_tomasi(), hessianfeatures()
"""
#=
September 2003  Original MATLAB version
August    2015  Ported to Julia
January   2016  Reworked to allow the N strongest features to be returned.
=#

function nonmaxsuppts{T<:Real}(cimg::Array{T,2}; radius::Real=1, thresh::Real=0, 
                               N::Int=typemax(Int), subpixel::Bool=false, img=[], fig = nothing)

    if radius < 1
        error("Radius must be 1 or greater")
    end

    (rows,cols) = size(cimg,1,2)
    
    # Extract local maxima by performing a grey scale morphological
    # dilation and then finding points in the corner strength image that
    # match the dilated image and are also greater than the threshold.
#    mx = imdilate(cimg, circularstruct(radius)) 
    mx = imdilate(cimg, "OCT", radius) # Much faster

    # Make mask to exclude points within radius of the image boundary. 
    bordermask = zeros(Bool, size(cimg))
    bordermask[radius+1:end-radius, radius+1:end-radius] = true
    
    # Find maxima, threshold, and apply bordermask
    cimgmx = (cimg.==mx) .& (cimg.>thresh) .& bordermask
    
    # Get row, col coords of points
    (r, c) = ind2sub(size(cimgmx),find(cimgmx)) 

    # Check if we want to limit the number of maxima
    if isfinite(N)   
        mxval = cimg[cimgmx]       # Get values of maxima and sort them.
        ind = sortperm(mxval,rev=true) # Permutation order of sorted values
        
        r = r[ind]                 # Reorder r and c arrays
        c = c[ind]                 # to match.
        
        if length(r) > N           # Grab the N strongest features.
            r = r[1:N]
            c = c[1:N]
        end
    end
    
    if isempty(r)     
        warn("No maxima above threshold found")
        r = zeros(Int,0)
        c = zeros(Int,0)
        return r, c
    end
    
    # If requested compute local maxima to sub pixel accuracy.  We fit
    # a 1D quadratic vertically and horizontally to the image data
    # using the pixels above and below and left and right of the
    # detected maxima and adjust the location of the maxima to match
    # maxima of the fitted quadartic.  Note that since the nonmaxima
    # radius is >= 1 we do not need to check whether we step outside
    # the image bounds when fitting the quadratics.
    if subpixel
        rsubpix = float(r)
        csubpix = float(c)

        for n = 1:length(r)
            # Solve for quadratic vertically
            cy = cimg[r[n],c[n]]
            ay = (cimg[r[n]-1,c[n]] + cimg[r[n]+1,c[n]])/2 - cy
            by = ay + cy - cimg[r[n]-1,c[n]]
            if abs(ay) > eps()
                rsubpix[n] += -by/(2*ay)  # Maxima of quadradic
            end

            # Solve for quadratic across columns
            cx = cimg[r[n],c[n]]
            ax = (cimg[r[n],c[n]-1] + cimg[r[n],c[n]+1])/2 - cx
            bx = ax + cx - cimg[r[n],c[n]-1]    
            if abs(ax) > eps()
                csubpix[n] += -bx/(2*ax)  # Maxima of quadradic
            end

        end      
    end

    # If an image has been supplied display it and overlay corners.
    if !isempty(img)  
        PyPlot.figure(fig)
        PyPlot.clf()
        PyPlot.imshow(img)
        PyPlot.set_cmap(PyPlot.ColorMap("gray"))
#        PyPlot.hold(true)  # No longer needed
        if subpixel
            PyPlot.plot(csubpix,rsubpix,"r+")
        else
            PyPlot.plot(c,r,"r+")
        end
        PyPlot.axis([1,cols,rows,1])
        PyPlot.title("Corners detected")
#        PyPlot.hold(false)  # No longer needed
    end
    
    if subpixel
        return rsubpix, csubpix
    else
        return r, c
    end
end


#---------------------------------------------------------------------
"""
erode1d - 1D morphological erosion of a signal.
```
Usage:  ef = erode1d(f, k)

Arguments:  f - 1D array of values to be eroded
            k - Size of erosion kernel 

Returns:   ef - Array of eroded values
```
Note that if the size of the structuring element is even then the centre pixel
is taken to be the integer pixel location to the right of the ideal.

This function uses Marcel van Herk's algorithm to run in linear time with
respect to the length of the signal, irrespective of the structing element
size.

See also: dilate1d(), imerode(), imdilate()
"""

#=
Reference: Marcel van Herk, "A fast algorithm for local minimum and maximum
filters on rectangular and octagonal kernels".  Pattern Recognition Letters 13
(1992) pp 517-521

PK April 2016
=#

# Two separate functions here because I am having trouble forming the
# function template using a parameterised Union as follows:
# function erode1d{T<:Real}(f::Union{Array{T,1}, BitArray{1}}, ki::Real)

function erode1d{T<:Real}(f::Array{T,1}, ki::Real)
    return erode_dilate1d(f, ki, "ERODE")
end

function erode1d(f::BitArray{1}, ki::Real)
    return erode_dilate1d(f, ki, "ERODE")
end

#---------------------------------------------------------------------
"""
dilate1d - 1D morphological dilation of a signal.
```
Usage:  df = dilate1d(f, k)

Arguments:  f - 1D array of values to be dilated
            k - Size of dilation kernel 

Returns:   df - Array of dilated values
```
Note that if the size of the structuring element is even then the centre pixel
is taken to be the integer pixel location to the right of the ideal.

This function uses Marcel van Herk's algorithm to run in linear time with
respect to the length of the signal, irrespective of the structing element
size.

See also: erode1d(), imdilate(), imerode()
"""

#=
Reference: Marcel van Herk, "A fast algorithm for local minimum and maximum
filters on rectangular and octagonal kernels".  Pattern Recognition Letters 13
(1992) pp 517-521

PK April 2016
=#

function dilate1d{T<:Real}(f::Array{T,1}, ki::Real)
    return erode_dilate1d(f, ki, "DILATE")
end

function dilate1d(f::BitArray{1}, ki::Real)
    return erode_dilate1d(f, ki, "DILATE")
end

#---------------------------------------------------------------------

# Unexported function that performs 1D erosion or dilation 

# ** This function is slow because it repeatedly allocates and frees
# ** space.  Need to develop a version erode_dilate1d!()

#function erode_dilate1d{T<:Real}(f::Array{T,1}, ki::Real, 
#                           erode_dilate::String)
function erode_dilate1d(f, ki::Real, erode_dilate::String)

 
    if uppercase(erode_dilate) == "ERODE"
        min_max = min
        typemin_max = typemax
    elseif  uppercase(erode_dilate) == "DILATE"
        min_max = max
        typemin_max = typemin
    else
        error("Option must be 'erode' or 'dilate' ")
    end
   
    df = zeros(size(f))
    Nf = length(f)
    k = round(Int, ki)
    
    # Determine the 'radius' of structuring element.  If the size is odd then
    # there is a centre pixel and the 'radius' is symmetric.  If the size is
    # even then the centre pixel is taken to be the integer pixel location to
    # the right of the ideal. In this case we have a 'left radius', rad1 and
    # a 'right radius' rad2.
    if isodd(k)
        rad1 = round(Int, (k-1)/2)
        rad2 = rad1
    else
        rad1 = round(Int, k/2)
        rad2 = rad1-1
    end
    
    # Pad ends of f with rad values of -Inf and ensure overall length of padded
    # array is a multiple of k
    extrapad = round(Int, mod(Nf + rad1+rad2, k))
    padval = typemin_max(typeof(f[1]));
    fpad = [repmat([padval], rad1); f[:]; repmat([padval], rad2+extrapad)]

    N = length(fpad)
    g=copy(fpad)
    h=copy(fpad)
    
    # Generate the intermediate arrays
    for o = 0:k:(N-k)
        for n = 2:k
            g[o+n] = min_max(g[o+n], g[o+n-1])
        end
        
        for n = (k-1):-1:1
            h[o+n] = min_max(h[o+n], h[o+n+1])
        end            
    end
    
    # Combine the intermediate arrays g and h to obtain the dilation.  Note
    # that for even sized structuring elements there is a wasted min() or max()
    # operation every k steps.  However in the interests of keeping the code
    # clean and simple this is accepted.
    for n = rad1+1:rad1+Nf
        df[n-rad1] = min_max(g[n+rad2], h[n-rad1])
    end

    return df
end


#---------------------------------------------------------------------
"""
imdilate - 2D morphological dilation with rectangular or octagonal structing element
```
Usage: dimg = imdilate(img, seType, seSize)

Arguments: img - Image to be dilated
        seType - String either "rectangle" or "octagon" specifying the
                 structuring element type. Can be shortened to "rect" or
                 "oct".
        seSize - Structuring element size.  
                 If the seType is 'rect' seSize can be a 2 element vector
                 indicating the size of the rectangular structuring element
                 [rows, cols], or if it is a single value the structuring
                 element is square [seSize x seSize]
                 If seType is 'oct' then seSize is the nominal radius of
                 the octagonal structuring element. 

Returns:  dimg - The dilated image.
```
Note that the radius of the octagonal structuring element is somewhat nominal
due to discrete approximations.  If anything the structuring element may end
up with a slightly larger radius than specified rather than being smaller.

This function uses Marcel van Herk's algorithm to run in linear time with
respect to the size of the image, irrespective of the structing element size.

See also imerode(), erode1d(), dilate1d()
"""

function imdilate{T<:Real, T2<:Integer}(img::Array{T,2}, seType::String, 
                                        seSize::Union{Array{T2,1},T2})
    imerode_dilate(img, seType, seSize, "DILATE")
end

function imdilate{T2<:Integer}(img::BitArray{2}, seType::String, 
                               seSize::Union{Array{T2,1},T2})
    imerode_dilate(img, seType, seSize, "DILATE")
end

#---------------------------------------------------------------------
"""
imerode - 2D morphological erosion with rectangular or octagonal structing element
```
Usage: eimg = imerode(img, seType, seSize)

Arguments: img - Image to be eroded
        seType - String either "rectangle" or "octagon" specifying the
                 structuring element type. Can be shortened to "rect" or
                 "oct".
        seSize - Structuring element size.  
                 If the seType is 'rect' seSize can be a 2 element vector
                 indicating the size of the rectangular structuring element
                 [rows, cols], or if it is a single value the structuring
                 element is square [seSize x seSize]
                 If seType is 'oct' then seSize is the nominal radius of
                 the octagonal structuring element. 

Returns:  eimg - The eroded image.
```
Note that the radius of the octagonal structuring element is somewhat nominal
due to discrete approximations.  If anything the structuring element may end
up with a slightly larger radius than specified rather than being smaller.

This function uses Marcel van Herk's algorithm to run in linear time with
respect to the size of the image, irrespective of the structing element size.

See also imdilate(), erode1d(), dilate1d()
"""

function imerode{T<:Real, T2<:Integer}(img::Array{T,2}, seType::String, 
                                       seSize::Union{Array{T2,1},T2})
    imerode_dilate(img, seType, seSize, "ERODE")
end

function imerode{T2<:Integer}(img::BitArray{2}, seType::String, 
                 seSize::Union{Array{T2,1},T2})
    imerode_dilate(img, seType, seSize, "ERODE")
end

#------------------------------------------------------------------------

# Unexported function that performs 1D erosion or dilation 

#=
Reference: Marcel van Herk, "A fast algorithm for local minimum and maximum
filters on rectangular and octagonal kernels".  Pattern Recognition Letters 13
(1992) pp 517-521

PK April 2016

Code is ugly and should eventually be restructured to accept a structuring
element defined via decomposition into a set of linear dilations.
=#

#function imdilate{T<:Real}(img::Union{Array{T,2}, BitArray{2}}, seType::String, seSize)
function imerode_dilate(img, seType::String, seSize, erode_dilate::String)

    if ndims(img) == 3
        error("Image must be binary or greyscale")
    end

    if !(uppercase(erode_dilate) == "ERODE"  || uppercase(erode_dilate) == "DILATE")
        error("Option must be 'erode' or 'dilate' ")
    end

    (rows,cols) = size(img,1,2)
    
    # Rectangular structuring element    
    if uppercase(seType[1:3]) == "REC"
        if length(seSize) == 1
            k = round.([Int], [seSize, seSize])  # [Int] for v0.5
        else
            k = round.([Int], seSize)
        end
        
        dimg = copy(img)
        
        for c = 1:cols
            dimg[:,c] = erode_dilate1d(img[:,c], k[1], erode_dilate)
        end
        
        for r = 1:rows
            dimg[r,:] = erode_dilate1d(vec(dimg[r,:]), k[2], erode_dilate)
            # vec() needed for 0.4 in line above
        end    
        
    # Octagonal structuring element        
    elseif uppercase(seType[1:3]) == "OCT"
     
        # The size of the linear structuring element for the
        # vertical and horizontal dilations matches the desired radius.
        k = round(Int, seSize)
        # The size of the linear structuring element for the dilations in the
        # diagonal direction is 1/sqrt(2) of k
        dk = round(Int, k/sqrt(2))
        # Make both structuring element sizes odd so that there is a centre
        # pixel and the resulting octagon is symmetric (not needed)
#        if iseven(dk), dk = dk+1; end
#        if iseven(round(k)), k = round(k)+1; end
        
        # First do 2D square dilation
        dimg = imerode_dilate(img, "rect", k, erode_dilate)
        
        # Extract diagonal lines of pixels from dimg and perform 1d dilation on
        # them
        
        # NE lines, emanating from left edge of image
        for r = 2:rows
            cmax = min(r,cols)
            l = zeros(cmax)
            for c = 1:cmax
                l[c] = dimg[r-c+1, c]
            end
            
            dl = erode_dilate1d(l, dk, erode_dilate)
            for c = 1:cmax
                dimg[r-c+1,c] = dl[c]
            end        
        end


        # NE lines, emanating from bottom of image    
        for c = 2:cols-1
            rmin = max(1, rows-(cols-c))
            l = zeros(rows-rmin+1)
            for r = rows:-1:rmin
                l[r-rmin+1] = dimg[r,c+rows-r]
            end
            
            dl = erode_dilate1d(l, dk, erode_dilate) 
            for r = rows:-1:rmin
                dimg[r,c+rows-r] = dl[r-rmin+1]
            end   
        end
        

        # SE lines, emanating from left edge of image
        for r = 1:rows-1
            cmax = min(rows-r+1,cols)
            l = zeros(cmax)
            for c = 1:cmax
                l[c] = dimg[r+c-1, c]
            end
            
            dl = erode_dilate1d(l, dk, erode_dilate)
            for c = 1:cmax
                dimg[r+c-1,c] = dl[c]
            end        
        end
        
        # SE lines, emanating from top of image    
        for c = 2:cols
            rmax = min(rows, cols-c+1)
            l = zeros(rmax)
            for r = 1:rmax
                l[r] = dimg[r,c+r-1]
            end
            
            dl = erode_dilate1d(l, dk, erode_dilate)
            for r = 1:rmax
                dimg[r,c+r-1] = dl[r]
            end        
        end
        
    else
        error("Structure element type must be 'rect' or 'oct' ")
    end
    
    return dimg
end


#------------------------------------------------------------------------
"""
imerode - 2D morpholgical erosion with arbitrary structuring element
```
Usage:   eimg = imerode(img, se)

Arguments:
          img - Image to be eroded.  May be greyscale or binary
                ::Array{T,2} or ::BitArray{2}
           se - Structuring element defined by a 2D Bool array

Returns:
         eimg - eroded image
```
See also: imdilate(), circularstruct()
"""

function imerode{T<:Real}(img::Array{T,2}, se)
    return imerode_dilate(img, se, "ERODE")
end

function imerode(img::BitArray{2}, se)
    return imerode_dilate(img, se, "ERODE")
end

#------------------------------------------------------------------------
"""
imdilate - 2D morpholgical dilation with arbitrary structuring element
```
Usage:   dimg = imdilate(img, se)

Arguments:
          img - Image to be dilated.  May be greyscale or binary
                ::Array{T,2} or ::BitArray{2}
           se - Structuring element defined by a 2D Bool array

Returns:
         dimg - dilated image
```
See also: imerode(), circularstruct()
"""

function imdilate{T<:Real}(img::Array{T,2}, se)
    return imerode_dilate(img, se, "DILATE")
end

function imdilate(img::BitArray{2}, se)
    return imerode_dilate(img, se, "DILATE")
end

#------------------------------------------------------------------
# Unexported function to perform image erosion or dilation depending
# on the whether minimum() or maximum() functions are passed for min_max

# Inefficient code!
# April 2016

function imerode_dilate(img, se, erode_dilate::String)

 
    if uppercase(erode_dilate) == "ERODE"
        min_max = minimum
    elseif  uppercase(erode_dilate) == "DILATE"
        min_max = maximum
    else
        error("Option must be 'erode' or 'dilate' ")
    end

    (rows,cols) = size(img,1,2)
    (sr,sc) = size(se)
    roff = round(Int, (sr-1)/2)
    coff = round(Int, (sc-1)/2)
    dimg = zeros(img)

    for r = 1:rows-sr, c = 1:cols-sc
        dimg[r+roff, c+coff] = min_max(img[r:r+sr-1,  c:c+sc-1].*se)
    end
    
    return dimg
end

#------------------------------------------------------------------------
"""
circularstruct - Generate circular structuring element for morphological operations
```
Usage:  se = circularstruct(radius)

Argument:  radius - Desired radius of structuring element
Returns:       se - Bool structuring element of size (2*radius+1)^2

```
See also: imdilate(), imerode()
"""
# April 2016

function circularstruct(radius::Real)
    se = [ sqrt(x^2+y^2) .<= radius  for x = -radius:radius, y = -radius:radius ]
end


#----------------------------------------------------------------------
"""
gaussfilt -  Small wrapper function for convenient Gaussian filtering
```
Usage:  smimg = gaussfilt(img, sigma)

Arguments:  img - Image to be smoothed, can be multi-channel.
                  ::Array{T,2}  or ::Array{T,3}
          sigma - Standard deviation of Gaussian filter.

Returns:  smimg - Smoothed image.
```
If called with sigma = 0 the function immediately returns with img assigned
to smimg

"""
# March 2010
# June  2013  - Provision for multi-channel images

function gaussfilt(img::Array, sigma::Real)
 
    if sigma < eps()
        return copy(img)   # ?? Return img or copy(img) ?
    end
    
    h = ImageFiltering.Kernel.gaussian(sigma)
#    h = Images.KernelFactors.IIRGaussian(sigma)

    if ndims(img) == 2      # Greyscale image
        return Images.imfilter(img,h)

    elseif ndims(img) == 3  # Multi-channel image
        nchan = size(img,3)
        smimg = zeros(size(img))

        for n = 1:nchan
            smimg[:,:,n] = Images.imfilter(img[:,:,n],h)
        end
   
        return smimg
    end 
end


#----------------------------------------------------------------------
"""
floatyx - Convert 2D ImageMeta to 2D float array with y x spatial order
```
Usage:  (fimg, prop) = floatyx(img)

Argument:  img - ::ImageMeta{T,2}

Returns:  fimg - ::Array{Float64,2} in "y" "x" spatial order.
          prop - A copy of the properties dictionary of the input image 
                 with "spatialorder" adjusted (if this was needed).
```
Most image processing functions expect 2D arrays in (row, column)
format and most feature localisation functions return corner and edge
coordinates in terms of (row, column) coordinates.

This convenience function takes a 2D ImageMeta, extracts the data,
converts it to a 2D Float64 array and, if necessary, transposes the data
so that it is in "y" "x" (row, column) spatial order.  The
ImageMeta spatial order property is set to ["y","x"] and the
properties returned.
"""

function floatyx{T}(img::Images.ImageMeta{T,2})

    fimg = float(Images.data(img))
    prop = copy(img.properties)

    if prop["spatialorder"] == ["x","y"]
        fimg = fimg'
        prop["spatialorder"] = ["y","x"]
    elseif img.properties["spatialorder"] != ["y","x"]
        error("Unable to handle the spatial order of this image")
    end

    return fimg, prop
end

#----------------------------------------------------------------------
"""
imgnormalise/imgnormalize - Normalises image values to 0-1, or to desired mean and variance

```
Usage 1:      nimg = imgnormalise(img)
```
Offsets and rescales image so that the minimum value is 0
and the maximum value is 1.  

```
Usage 2:      nimg = imgnormalise(img, reqmean, reqvar)

Arguments:  img     - A grey-level input image.
            reqmean - The required mean value of the image.
            reqvar  - The required variance of the image.
```
Offsets and rescales image so that nimg has mean reqmean and variance
reqvar.  
"""

# Normalise 0 - 1
function imgnormalise(img::Array) 
    n = img - minimum(img)
    return n = n/maximum(n)
end

# Normalise to desired mean and variance
function imgnormalise(img::Array, reqmean::Real, reqvar::Real)
    n = img - mean(img)
    n = n/std(img)      # Zero mean, unit std dev
    return n = reqmean + n*sqrt(reqvar)
end

# For those who spell normalise with a 'z'
"""
imgnormalize - Normalizes image values to 0-1, or to desired mean and variance
```
Usage 1:      nimg = imgnormalize(img)
```
Offsets and rescales image so that the minimum value is 0
and the maximum value is 1.  
```
Usage 2:      nimg = imgnormalize(img, reqmean, reqvar)

Arguments:  img     - A grey-level input image.
            reqmean - The required mean value of the image.
            reqvar  - The required variance of the image.
```
Offsets and rescales image so that nimg has mean reqmean and variance
reqvar.  
"""

function imgnormalize(img::Array) 
    return imgnormalise(img)
end

function imgnormalize(img::Array, reqmean::Real, reqvar::Real)
    return imgnormalise(img, reqmean, reqvar)
end

#----------------------------------------------------------------------
"""
histtruncate - Truncates ends of an image histogram.

Function truncates a specified percentage of the lower and
upper ends of an image histogram.

This operation allows grey levels to be distributed across
the primary part of the histogram.  This solves the problem
when one has, say, a few very bright values in the image which
have the overall effect of darkening the rest of the image after
rescaling.

```
Usage: 
1)   newimg = histtruncate(img, lHistCut, uHistCut)
2)   newimg = histtruncate(img, HistCut)

Arguments:
 Usage 1)
   img         -  Image to be processed.
   lHistCut    -  Percentage of the lower end of the histogram
                  to saturate.
   uHistCut    -  Percentage of the upper end of the histogram
                  to saturate.  If omitted or empty defaults to the value
                  for lHistCut.
 Usage 2)
   HistCut     -  Percentage of upper and lower ends of the histogram to cut.

Returns:
   newimg      -  Image with values clipped at the specified histogram
                  fraction values.  If the input image was colour the
                  lightness values are clipped and stretched to the range
                  0-1.  If the input image is greyscale no stretching is
                  applied. You may want to use imgnormalise() to achieve this.
```
See also: imgnormalise()
"""

# July      2001 - Original version
# February  2012 - Added handling of NaN values in image
# February  2014 - Code cleanup
# September 2014 - Default for uHistCut + cleanup

function  histtruncate(img::Array, lHistCut::Real, uHistCut::Real)
    
    if lHistCut < 0 || lHistCut > 100 || uHistCut < 0 || uHistCut > 100
	error("Histogram truncation values must be between 0 and 100")
    end
    
    if ndims(img) > 2
	error("histtruncate only defined for grey scale images")
    end

    newimg = copy(img)    
    sortv = sort(newimg[:])   # Generate a sorted array of pixel values.

    # Any NaN values will end up at the end of the sorted list. We
    # need to ignore these.
    N = sum(.!isnan.(sortv))  # Number of non NaN values.
    
    # Compute indicies corresponding to specified upper and lower fractions
    # of the histogram.
    lind = floor(Int, 1 + N*lHistCut/100)
    hind =  ceil(Int, N - N*uHistCut/100)

    low_val  = sortv[lind]
    high_val = sortv[hind]

    # Adjust image
    newimg[newimg .< low_val] = low_val
    newimg[newimg .> high_val] = high_val
    
    return newimg
end


function  histtruncate(img::Array, HistCut::Real)
    return histtruncate(img, HistCut, HistCut)
end


#----------------------------------------------------------------------
"""
matchbycorrelation - Match image feature points by correlation

Function generates putative matches between previously detected
feature points in two images by looking for points that are maximally
correlated with each other within windows surrounding each point.
Only points that correlate most strongly with each other in both
directions are returned.

This implements normalised cross correlation in its most basic form.
There is no attempt to deal with scale or orientation differences
between the images and thus it will only work if the images are not
too dissimilar.

This is a simple-minded N^2 comparison.

```
Usage: (m1, m2, p1ind, p2ind, cormat) = 
                matchbycorrelation(img1, p1, img2, p2, w, dmax)

Arguments:
      img1, img2 - Images containing points that we wish to match.
        p1, p2   - Coordinates of feature pointed detected in img1 and
                   img2 respectively using a corner detector (say Harris
                   or phasecong3).  p1 and p2 are [2 x npts] arrays though
                   p1 and p2 are not expected to have the same number
                   of points.  The first row of p1 and p2 gives the row
                   coordinate of each feature point, the second row
                   gives the column of each point.
        w        - Window size (in pixels) over which the correlation
                   around each feature point is performed.  This should
                   be an odd number.
        dmax     - (Optional) Maximum search radius for matching
                   points.  Used to improve speed when there is little
                   disparity between images. Even setting it to a generous
                   value of 1/4 of the image size gives a useful
                   speedup. If this parameter is omitted it defaults to Inf. 

Returns:
        m1, m2   - Coordinates of points selected from p1 and p2
                   respectively such that (putatively) m1[:,i] matches
                   m2[:,i]. m1 and m2 are [2 x npts] arrays defining the
                   points in each of the images in the form [row;col].
  p1ind, p2ind   - Indices of points in p1 and p2 that form a match.  Thus,
                   m1 = p1[:,p1ind] and m2 = p2[:,p2ind]
        cormat   - Correlation matrix; rows correspond to points in p1,
                   columns correspond to points in p2
```
This function is slow as mud! Needs some attention.
"""
# Probably needs devectorisation to improve speed 

#=
February 2004    - Original version
May      2004    - Speed improvements + constraint on search radius for
                   additional speed
August   2004    - Vectorized distance calculation for more speed
                   (thanks to Daniel Wedge)
December 2009    - Added return of indices of matching points from original
                   point arrays
May 2016         - Ported to Julia
=#

function matchbycorrelation(img1i, p1, img2i, p2, w, dmax=Inf)

    #-------------------------------------------------------------------------
    # Function that does the work.  This function builds a correlation matrix
    # that holds the correlation strength of every point relative to every
    # other point.  While this seems a bit wasteful we need all this data if
    # we want to find pairs of points that correlate maximally in both
    # directions.
    #
    # This code assumes img1 and img2 have zero mean.  This speeds the
    # calculation of the normalised correlation measure.
    
    function correlationmatrix(img1, p1, img2, p2, w, dmax)

        if iseven(w)
            error("Window size should be odd")
        end
        
        (rows1, npts1) = size(p1,1,2)
        (rows2, npts2) = size(p2,1,2)    
        if rows1 != 2 || rows2 != 2
            error("Feature points must be specified in 2xN arrays")
        end
        
        # Initialize correlation matrix values to -infinty
        cormat = -Inf*ones(npts1,npts2)
        
        (im1rows, im1cols) = size(img1,1,2)
        (im2rows, im2cols) = size(img2,1,2)    

        r = (w-1)÷2   # 'radius' of correlation window

        # For every feature point in the first image extract a window of data
        # and correlate with a window corresponding to every feature point in
        # the other image.  Any feature point less than distance 'r' from the
        # boundary of an image is not considered.
        
        # Find indices of points that are distance 'r' or greater from
        # boundary on image1 and image2;
        n1ind = find((p1[1,:].>r) .& (p1[1,:].<im1rows+1-r) .& (p1[2,:].>r) .& (p1[2,:].<im1cols+1-r))
        n2ind = find((p2[1,:].>r) .& (p2[1,:].<im2rows+1-r) .& (p2[2,:].>r) .& (p2[2,:].<im2cols+1-r))
        
        dists2 = zeros(length(n2ind))

        for n1 = n1ind            
            # Generate window in 1st image   	
            w1 = vec(img1[p1[1,n1]-r:p1[1,n1]+r, p1[2,n1]-r:p1[2,n1]+r])
            # Pre-normalise w1 to a unit vector.
            w1 = w1/sqrt(w1'*w1)
            
            # Identify the indices of points in p2 that we need to consider.
            if dmax == Inf
	        n2indmod = n2ind # We have to consider all of n2ind
                
            else     # Compute distances from p1[:,n1] to all available p2.
                ind = 0
                for i = n2ind
                    ind += 1
                    dists2[ind] = sum((p1[:,n1]-p2[:,i]).^2)
                end
	        n2indmod = n2ind[find(dists2 .< dmax^2)]
            end
            
            # Calculate normalised correlation measure.  Note this gives
            # significantly better matches than the unnormalised one.
            for n2 = n2indmod 
                # Generate window in 2nd image
                w2 = vec(img2[p2[1,n2]-r:p2[1,n2]+r, p2[2,n2]-r:p2[2,n2]+r])
                cormat[n1,n2] = (w1'*w2/sqrt(w2'*w2))[1]  # [1] to create scalar
            end
        end

        return cormat
    end
    #-----------------------------------------------------------------
    # The main function
    
    img1 = float(img1i)
    img2 = float(img2i)
    
    # Subtract image smoothed with an Gaussian filter with sigma = w/3 from
    # each of the images.  This attempts to compensate for brightness
    # differences in each image.  Doing it now allows faster correlation
    # calculation.
    sigma = w/3   # A guess as to what might be appropriate...
    img1 .-= gaussfilt(img1, sigma)
    img2 .-= gaussfilt(img2, sigma)
    
    # Generate correlation matrix
    cormat = correlationmatrix(img1, p1, img2, p2, w, dmax)
    (corrows,corcols) = size(cormat,1,2)

    # Find max along rows give strongest match in p2 for each p1
    mp2forp1 = zeros(corrows)
    colp2forp1 = zeros(Int,corrows)
    for r = 1:corrows
        (mp2forp1[r], colp2forp1[r]) = findmax(cormat[r,:])  
    end
    # Find max down cols give strongest match in p1 for each p2 
    mp1forp2 = zeros(corcols)
    rowp1forp2 = zeros(Int,corcols)
    for c = 1:corcols
        (mp1forp2[c], rowp1forp2[c]) = findmax(cormat[:,c])
    end
    
    # Now find matches that were consistent in both directions
    p1ind = zeros(Int, size(p1,2))  # Arrays for storing matched indices
    p2ind = zeros(Int, size(p2,2))
    indcount = 0
    for n = 1:corrows
        if rowp1forp2[colp2forp1[n]] == n  # consistent both ways
            indcount += 1
            p1ind[indcount] = n
            p2ind[indcount] = colp2forp1[n]
        end
    end
    
    # Trim arrays of indices of matched points
    p1ind = p1ind[1:indcount]
    p2ind = p2ind[1:indcount]
    
    # Extract matched points from original arrays
    m1 = p1[:,p1ind]
    m2 = p2[:,p2ind]
    
    return (m1, m2, p1ind, p2ind, cormat)
    
end # of matchbycorrelation

#----------------------------------------------------------------------
"""
grey2census - Convert image grey scale values to census values
```
Usage:  cimg = grey2census(img, window)

Arguments:
          img - greyscale image to be processed
       window - Optional 2-vector [rows,cols] specifying the size of the
                window to be considered in computing the census
                transform. Defaults to [7, 9].  The values must be odd and
                their product must not be greater than 64.
Returns:
         cimg - Census encoded UInt64 image.
```
Each pixel is encoded with a bit pattern formed by comparing the pixel with
the pixels in the window around it. If a window pixel is less than the centre
pixel its corresponding bit in the census encoding is set.  This provides an
encoding that describes a pixel in terms of its surrounding pixels in a way
that is invariant to lighting variations.  Note, however, that the encoding is
dependent on the image orientation.

Use the Hamming distance to compare encoded pixel values when matching.
```
  hammingDist = count_ones(img1[r,c] \$ img2[r,c])
```
"""

# Reference: Ramin Zabih and John Woodfill. Non-parametric Local
# Transforms for Computing Visual Correspondence. European Conference
# on Computer Vision, Stockholm, Sweden, May 1994, pages 151-158

# ** There is a type instability that makes this slower than it should
# ** be and the array slicing must be costly.


function grey2census{T<:Real,T2<:Integer}(img::Array{T,2}; window::Vector{T2}=[7,9], medianref = true)

    (rows,cols) = size(img,1,2)    

    nBits = window[1]*window[2]
    if nBits > 64
        error("Window must not have more than 64 elements")
    end
    
    if iseven(window[1]) || iseven(window[2])
        error("Window sizes must be odd")
    end
    
    rrad = round(Int,(window[1]-1)/2)
    crad = round(Int,(window[2]-1)/2)

    
    # Define the reference values.  If we use the median of the window we
    # apply a median filter to obtain the reference values, else use the
    # original image
    if medianref
        imref = medfilt2(img, window)
    else
        imref = copy(img)
    end
    
    # Build table of values corresponding to each bit
    bitval = UInt64[2^(b-1) for b=1:nBits]
    
    cimg = zeros(UInt64, rows, cols)

    # Build meshgrid of coordinates corresponding to all the offsets/bits
    coff = [c for r = -rrad:rrad, c = -crad:crad]
    roff = [r for r = -rrad:rrad, c = -crad:crad]

    for n = 1:nBits
        cimg[rrad+1:rows-rrad, crad+1:cols-crad] = 
        cimg[rrad+1:rows-rrad, crad+1:cols-crad] + 
#=
        bitval[n]*(
                img[rrad+1+roff[n]:rows-rrad+roff[n], crad+1+coff[n]:cols-crad+coff[n]] .< 
                img[rrad+1:rows-rrad, crad+1:cols-crad])
=#
        bitval[n]*(
                view(img,rrad+1+roff[n]:rows-rrad+roff[n], crad+1+coff[n]:cols-crad+coff[n]) .< 
                view(img,rrad+1:rows-rrad, crad+1:cols-crad))

    end
    
    return cimg
end


#-----------------------------------------------------------------
"""
briefcoords - Compute BRIEF descriptor sampling coordinates within a patch

```
Usage:  rc = briefcoords(S, nPairs, UorG; disp=false)

Arguments:  S - An odd integer specifying the patch size (S x S).  
       nPairs - Number of point pairs required.  Typically this is a power
                of 2, 128, 256, 512 for packing the descriptor values into
                a bit string.
         UorG - Character 'U' or 'G' indicating whether a uniform or
                gaussian distribution of point pairings is formed within
                the patch.  
         disp - Optional boolean flag indicating whether a plot of the point
                pairings should be displayed.

Returns:   rc - [2 x 2*nPairs] array of integer (row; col) coordinates
                with all values in the range -(S-1)/2..(S-1)/2.  Each
                successive pair of columns is intended to provide a pair of
                points for comparing image grey values against each other
                for forming a BRIEF descriptor of a patch about some
                central feature location. 
```

Use of the Gaussian distribution of point pairings corresponds to the approach
suggested by Calonder et al in their original paper.  Note, however that for
small patches and a large number of pairings one will get repeated pairings
which reduces the degrees of freedom in the final descriptor.  The uniformly
distributed option does not result in repeated pairings making it the
preferred option for small patch sizes.

See also: grey2lbp(), grey2lbp!()
"""

# Reference:
# Michael Calonder, Vincent Lepetit, Christoph Strecha, and Pascal Fua.
# BRIEF: Binary Robust Independent Elementary Features

# To do: Constrain point selection to be within a circular region
# around the centre pixel, rather than a square one.  This will make
# the size and localisation invariant to orientation of the underlying
# pixel data.

function briefcoords(S, nPairs, UorG; disp=false)
    
    if iseven(S)
        error("Patch size must be an odd integer")
    end

    # Fix the initialisation of the random number generator so that we get
    # repeatable sampling point locations.
    rng = MersenneTwister(1)

    R = round(Int, (S-1)/2)    # 'Radius' of patch

    # Option that generates a uniformly distributed pairings of
    # coordinates within the patch.    
    if uppercase(UorG[1]) == 'U'
        rc = rand(rng, -R:R, (2, 2*nPairs))
        
    # Option that generates a set of normally distributed points within the
    # patch and forms pairings between them as suggested by the orginal paper.
    elseif uppercase(UorG[1]) == 'G'
        sigma = S/5  # Standard deviation suggested by Calonder et al.
        rc = round(Int, sigma * randn(rng, 2, 2*nPairs))
        
        # Clamp coordinates to patch bounds
        rc[rc .> R] = R
        rc[rc .< -R] = -R
        
    else
        error("Option must be 'U'niform or 'G'aussian")
    end

    # Diagnostic display
    if disp
        figure(200); clf
        R = (S-1)/2 + 1
        axis([-R, R, -R, R])

        for n = 1:nPairs
            plot(rc[2, 2*n-1:2*n], rc[1, 2*n-1:2*n], color = rand(3), 
                 linewidth = 3)
            hold(true)
        end

        axis("equal")
        hold(false)

        # Determine distances between pairs and display histogram
        sdist = zeros(nPairs)
        for n = 1:nPairs
            sdist[n] = norm(rc[:,2*n-1] - rc[:,2*n])
        end
        
        (e, counts) = hist(sdist, 0:.5:S)
        figure(30); clf()
        PyPlot.bar(left=e[1:end-1], height = counts, width=e[2]-e[1])
        title("Histogram of distances between pairs")
    end

    return rc
end

#-----------------------------------------------------------------
"""
grey2lbp - Convert image grey scale values to local binary pattern

```
Usage:  lbpim = grey2lbp(img, rc, window)

Arguments:
          img - greyscale image to be processed
           rc - [2 x 2*nPairs] array of integer (row;col) coordinates.
                Each successive pair of columns provides a pair of
                points for comparing image grey values against each other
                for forming a binary descriptor of a patch about some
                central feature location. Array rc must not have more than
                128 columns corresponding to an encoding of 64 bits,
       window - 2-vector specifying the size of the window to be considered
                in computing the local standard deviation for determining
                regions to be encoded with a 0.  The window size values
                should match the range of values in rc and be odd.

Returns:
        lbpim - Local Binary Pattern encoded UInt64 image.
```

Each pixel is encoded with a bit pattern formed by comparing pairs of pixels.
If the difference between pixels is -ve they are are encoded with 1, else 0.
This provides an encoding that describes a pixel in terms of its surrounding
pixels in a way that is invariant to lighting variations.  Note, however, that
the encoding is dependent on the image orientation.   The encoding is also
sensitive to noise on near constant regions.

Use the Hamming distance to compare encoded pixel values when matching.

See also: grey2lbp!(), grey2census(), briefcoords()
"""

function grey2lbp(img::Array, rc::Array, window::Vector)

    if ndims(img) == 3
        error("Image must be grey scale")
    end

    (rows,cols) = size(img,1,2)
    lbpim = zeros(UInt64, rows, cols)
    
    nBits = div(size(rc,2),2)       # No of coordiante pairs / No of bits
    
    if nBits > 64
        error("Array 'rc' must not have more than 128 columns")
    end
    
    if iseven(window[1]) || iseven(window[2])
        error("Window sizes must be odd")
    end
    
    rrad = div((window[1]-1),2)
    crad = div((window[2]-1),2)

    # Build table of values corresponding to each bit. 
    bitval = [UInt64(1) << n for n = 0:63]

    # Perform encoding
    # Jumping around in img is not cache friendly
    for b = 1:nBits
        twob = 2*b
        twobm1 = 2*b-1
        roff1 = rc[1,twob]
        roff2 = rc[1,twobm1]
        for c = crad+1:cols-crad
            c1 = c+rc[2,twob]
            c2 = c+rc[2,twobm1]
            
            for r = rrad+1:rows-rrad
                if img[r+roff1, c1] < img[r+roff2, c2] 
                    lbpim[r,c] |= bitval[b]  
                end
            end
        end
    end

    return lbpim
end

#----------------------------
# Inplace version

function grey2lbp!(lbpim::Array{UInt64,2}, img::Array, rc::Array, window::Vector)

    if ndims(img) == 3
        error("Image must be grey scale")
    end

    (rows,cols) = size(img,1,2)
    fill!(lbpim, zero(UInt64))
    
    nBits = div(size(rc,2),2)       # No of coordiante pairs / No of bits
    
    if nBits > 64
        error("Array 'rc' must not have more than 128 columns")
    end
    
    if iseven(window[1]) || iseven(window[2])
        error("Window sizes must be odd")
    end
    
    rrad = div((window[1]-1),2)
    crad = div((window[2]-1),2)

    # Build table of values corresponding to each bit. 
    bitval = [UInt64(1) << n for n = 0:63]

    rrange = rrad+1:rows-rrad
    crange = crad+1:cols-crad

    nrows = length(rrange)

    # Perform encoding
    # Jumping around in img is not cache friendly
    @inbounds begin
        for b = 1:nBits
            twob = 2*b
            twobm1 = twob-1
            roff1 = rc[1,twob]
            roff2 = rc[1,twobm1]

            for c = crange
                c1 = c+rc[2,twob]
                c2 = c+rc[2,twobm1]
                
                # This loop is the killer for time, what can one do?
                for r = rrange
                    if img[r+roff1, c1] < img[r+roff2, c2] 
                        lbpim[r,c] |= bitval[b]
                    end
                end
            end
        end
    end # @inbounds

    return nothing
end

#------------------------------------------------------------------    
"""
medfilt2 - Convenience wrapper for median filtering

```
Usage:  medimg = medfilt2(img, h::Int, w::Int)
        medimg = medfilt2(img, h::Int)   
        medimg = medfilt2(img, hw::Tuple)

Arguments:   img - Image to be processed, Array{T,2}
            h, w - Height and width of rectangular window over which the 
                   median is to be computed. Values must be odd.
                   h and w may be specified as a size tuple, or as a
                   single value in which case w and h are made equal.

Returns:  medimg - Median filtered image.

```
"""

# ** To be rewritten using Huang's algorithm or Perreault and Hebert's
# ** algorithm

function medfilt2{T<:Real}(img::Array{T,2}, h::Int, w::Int)

    if iseven(h) || iseven(w)
        error("Median filter region size must be odd")
    end

    return mapwindow(median!, img, (h,w))
end

function medfilt2(img::Array, sze::Tuple{Int, Int})
    return medfilt(img, sze[1], sze[2])
end

function medfilt2{T<:Real}(img::Array{T,2}, sze::Int)
    medfilt2(img, sze, sze)
end

#------------------------------------------------------------------------------
"""
stdfilt2 - Compute local standard deviation across an image

```
Usage:  stdimg = stdfilt2(img, h::Int, w::Int)
        stdimg = stdfilt2(img, h::Int)
        stdimg = stdfilt2(img, hw::Tuple)

Arguments:   img - Image to be processed, Array{T,2}
            h, w - Height and width of rectangular window over which the 
                   standard deviation is to be computed. Values must be odd.
                   h and w may be specified as a size tuple, or as a
                   single value in which case w and h are made equal.

Returns:  stdimg - Standard deviation image.

```
"""

function stdfilt2(img::Array, h::Int, w::Int)

    if ndims(img) != 2
        error("Image must be 2D")
    end
    
    (rows,cols) = size(img,1,2)
    
    # Get mean image
    kern = Images.centered(Images.imaverage((h,w)))
    meanimg::Array{Float64,2} = Images.imfilter(img, kern)

    # Subtract the original image from the mean image and square
    meanimg .= (meanimg .- img).^2

    # Perform averaging to get local variance
    stdimg::Array{Float64,2} = Images.imfilter(meanimg, kern)

    # Take sqrt to get standard deviation
    stdimg .= sqrt.(stdimg)

    return stdimg
end

function stdfilt2(img::Array, sze::Tuple{Int, Int})
    return stdfilt(img, sze[1], sze[2])
end

function stdfilt2(img::Array, sze::Int)
    return stdfilt(img, sze, sze)
end


#----------------------------------------------------------------------
"""
keypause - Wait for user to hit return before continuing
```
Usage:  keypause()
```
"""

function keypause()
    println("Hit return to continue, or 'x' to exit")
    a = readline()
    if a[1] == 'x'
        error("Exiting")  # Should be a nice way to do this
    else
        return 
    end
end

#----------------------------------------------------------------------

function testinput{T<:Real}(img::Union{T, Array{T,2}, BitArray{2},Vector{T}})
println(typeof(img))
end
