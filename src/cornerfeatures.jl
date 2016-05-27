#=--------------------------------------------------------------------

cornerfeatures - Functions for the detection of 'corner' features 


Copyright (c) 2015-2016 Peter Kovesi
peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK August 2015

---------------------------------------------------------------------=#

export structuretensor, harris, noble, shi_tomasi, coherence
export hessianfeatures
export fastradial

import Images

#----------------------------------------------------------------------
"""
harris - Harris corner detector
```
Usage:                cimg = harris(img, sigma=1; k=0.04)
              (cimg, r, c) = harris(img, sigma=1; k=0.04, args...)

Arguments:   
           img    - Image to be processed.
                    ::Array{T,2} or ::AbstractImage{T,2} 
           sigma  - Standard deviation of Gaussian summation window. 
                    Typical values to use might be 1-3. Default is 1.

Keyword arguments:
                k - The scaling factor to apply to the trace of the 
                    structure tensor. Defaults to the traditional value 
                    of 0.04 

          args... - A sequence of keyword arguments for performing non-maximal suppression.
                    These are: (see nonmaxsuppts() for more details)
           radius - Radius of region considered in non-maximal
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

Returns:
           cimg   - Corner strength image.  
           r      - Row coordinates of corner points.
           c      - Column coordinates of corner points.
```
If the input is an AbstractImage cimg will be an AbstractImage with
property spatialorder set to ["y","x"] so that it is consistent with
the coordinates returned in r, c

With none of the optional nonmaximal keyword arguments specified only
`cimg` is returned as a raw corner strength image.  You may then want
to look at the values within 'cimg' to determine the appropriate
threshold value to use. Note that the Harris corner strength varies
with the intensity gradient raised to the 4th power.  Small changes in
input image contrast result in huge changes in the appropriate
threshold.

The Harris measure is det(M) - k*trace^2(M), where k is a parameter
you have to set (traditionally k = 0.04) and M is the structure
tensor.  Use Noble's measure if you wish to avoid the need to set a
parameter k.  However the Shi-Tomasi measure is probably what you
really want.

See also: noble(), shi_tomasi(), hessianfeatures(), nonmaxsuppts(), derivative5()
"""
#=
References: 
C.G. Harris and M.J. Stephens. "A combined corner and edge detector", 
Proceedings Fourth Alvey Vision Conference, Manchester.
pp 147-151, 1988.

March     2002 - Original version
August    2005 - Changed so that code calls nonmaxsuppts
August    2010 - Changed to use Farid and Simoncelli's derivative filters
September 2015 - Ported to Julia
January   2016 - noble() made distict from harris() and argument handling 
                 changed 
=#

function  harris{T<:Real}(img::Array{T,2}, sigma::Real=1; k::Real=0.04, args...)

    (Ix2, Iy2, Ixy) = structuretensor(img, sigma)
    cimg = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2 

    if isempty(args)
        return cimg
    else
        (r,c) = nonmaxsuppts(cimg; args...)
        return cimg, r, c
    end

end

# Case for img::AbstractImage 
function harris{T}(img::Images.AbstractImage{T,2}, sigma::Real=1; k::Real=0.04, args...)

    # Extract 2D float array in y x spatial order
    (dimg, prop) = floatyx(img)

    if isempty(args)
        cimg = harris(dimg, sigma, k=k)
        return Images.Image(cimg, prop)
    else
        (cimg, r, c) = harris(dimg, sigma, k=k; args...)
        return Images.Image(cimg, prop), r, c
    end
end


#----------------------------------------------------------------------
"""
noble - Noble's variant of the Harris corner detector
```

Usage:                cimg = noble(img, sigma=1)
              (cimg, r, c) = noble(img, sigma=1; args...)

Arguments:   
           img    - Image to be processed.
                    ::Array{T,2} or ::AbstractImage{T,2} 
           sigma  - Standard deviation of Gaussian summation window. 
                    Typical values to use might be 1-3. Default is 1.
Keyword arguments:
          args... - A sequence of keyword arguments for performing non-maximal suppression.
                    These are: (see nonmaxsuppts() for more details)
           radius - Radius of region considered in non-maximal
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

Returns:
           cimg   - Corner strength image.
           r      - Row coordinates of corner points.
           c      - Column coordinates of corner points.
```
If the input is an AbstractImage cimg will be an AbstractImage with
property spatialorder set to ["y","x"] so that it is consistent with
the coordinates returned in r, c

With none of the optional keyword arguments specified only `cimg` is
returned as a raw corner strength image.  You may then want to look at
the values within 'cimg' to determine the appropriate threshold value
to use.

Noble's variation of the Harris operator avoids the need to set the
parameter `k`.  However the Shi-Tomasi measure is probably what you
really want.

See also: harris(), shi_tomasi(), hessianfeatures(), nonmaxsuppts(),
structuretensor(), derivative5()
"""
#=
References: 
C.G. Harris and M.J. Stephens. "A combined corner and edge detector", 
Proceedings Fourth Alvey Vision Conference, Manchester.
pp 147-151, 1988.

Alison Noble, "Descriptions of Image Surfaces", PhD thesis, Department
of Engineering Science, Oxford University 1989, p45.

March     2002 - Original version
August    2005 - Changed so that code calls nonmaxsuppts
August    2010 - Changed to use Farid and Simoncelli's derivative filters
September 2015 - Ported to Julia
January   2016 - noble() made distict from harris() and argument handling 
                 changed 
=#

function  noble{T<:Real}(img::Array{T,2}, sigma::Real=1; args...)

    (Ix2, Iy2, Ixy) = structuretensor(img, sigma)
    cimg = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps()) 

    if !isempty(args)
        (r,c) = nonmaxsuppts(cimg; args...)
        return cimg, r, c
    else
        return cimg
    end

end

# Case for img::AbstractImage 
function noble{T}(img::Images.AbstractImage{T,2}, sigma::Real=1; args...)

    # Extract 2D float array in y x spatial order
    (dimg, prop) = floatyx(img)

    if isempty(args)
        cimg = noble(dimg, sigma)
        return Images.Image(cimg, prop)
    else
        (cimg, r, c) = noble(dimg, sigma; args...)
        return Images.Image(cimg, prop), r, c
    end
end


#----------------------------------------------------------------------
"""
shi_tomasi - Shi - Tomasi corner detector
```
Usage:               cimg = shi_tomasi(img, sigma=1)
             (cimg, r, c) = shi_tomasi(img, sigma=1; args...)

Arguments:   
           img    - Image to be processed.
                    ::Array{T,2} or ::AbstractImage{T,2} 
           sigma  - Standard deviation of Gaussian summation window. 
                    Typical values to use might be 1-3. Default is 1.
Keyword arguments:
          args... - A sequence of keyword arguments for performing non-maximal suppression.
                    These are: (see nonmaxsuppts() for more details)
           radius - Radius of region considered in non-maximal
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

Returns:
           cimg   - Corner strength image.
           r      - Row coordinates of corner points.
           c      - Column coordinates of corner points.
```
If the input is an AbstractImage cimg will be an AbstractImage with
property spatialorder set to ["y","x"] so that it is consistent with
the coordinates returned in r, c

With no keyword arguments supplied only 'cimg' is returned as a raw
corner strength image.  You may then want to look at the values within
'cimg' to determine the appropriate threshold value to use.

The Shi - Tomasi measure returns the minimum eigenvalue of the
structure tensor.  This represents the ideal that the Harris and Noble
detectors attempt to approximate.  Back in 1988 Harris wanted to avoid
the computational cost of taking a square root when computing
features.  This is no longer relevant today!

See also: harris(), noble(), hessianfeatures(), structuretensor(),
nonmaxsuppts(), derivative5()
"""
#=
Reference: 
J. Shi and C. Tomasi. "Good Features to Track,". 9th IEEE Conference on
Computer Vision and Pattern Recognition. 1994.

January 2016 - Original version
=#

function shi_tomasi{T1<:Real}(img::Array{T1,2}, sigma::Real=1; args...)

    (Ix2, Iy2, Ixy) = structuretensor(img, sigma)    

    T = Ix2 + Iy2                 # trace
    D = Ix2.*Iy2 - Ixy.^2         # determinant
    
    # The two eigenvalues of the 2x2 structure tensor are:
    # L1 = T/2 + sqrt(T.^2/4 - D)
    # L2 = T/2 - sqrt(T.^2/4 - D)

    # We just want the minimum eigenvalue
    cimg = T/2 - sqrt(T.^2/4 - D + eps())
    
    if !isempty(args)
        (r,c) = nonmaxsuppts(cimg; args...)
        return cimg, r, c
    else
        return cimg
    end

end    

# Case for img::AbstractImage 
function shi_tomasi{T}(img::Images.AbstractImage{T,2}, sigma::Real=1; args...)

    # Extract 2D float array in y x spatial order
    (dimg, prop) = floatyx(img)

    if isempty(args)
        cimg = shi_tomasi(dimg, sigma)
        return Images.Image(cimg, prop)
    else
        (cimg, r, c) = shi_tomasi(dimg, sigma; args...)
        return Images.Image(cimg, prop), r, c
    end
end


#----------------------------------------------------------------------
"""
coherence - Compute image coherence from structure tensor
```
Usage:      cimg = coherence(img, sigma=1)

Arguments:   
           img    - Image to be processed.
                    ::Array{T,2} or ::AbstractImage{T,2} 
           sigma  - Standard deviation of Gaussian summation window. 
                    Typical values to use might be 1-3. Default is 1.
Returns:
           cimg   - Coherence image.
```

Coherence is defined as the difference between the eigenvalues of the
structure tensor divided by the sum of the eigenvalues, all squared.  
```
     ((L1-L2)./(L1+L2)).^2
```
A small value indicates that the eigenvalues are similar in magnitude
and hence the local structure has a strong 2D component.

If the input is an AbstractImage cimg will be an AbstractImage with
property spatialorder set to ["y","x"].

See also: structuretensor()
"""

function coherence{T1<:Real}(img::Array{T1,2}, sigma::Real=1)

    (Ix2, Iy2, Ixy) = structuretensor(img, sigma)    

    T = Ix2 + Iy2                 # trace
    D = Ix2.*Iy2 - Ixy.^2         # determinant
    
    # The two eigenvalues of the 2x2 structure tensor are:
    L1 = T/2 + sqrt(T.^2/4 - D)
    L2 = T/2 - sqrt(T.^2/4 - D)

    # Coherence is defined as ((L1-L2)./(L1+L2)).^2
    # this can be reduced to the following
    return 4*(T.^2/4 - D) ./ T.^2

end

# Case for img::AbstractImage 
function coherence{T}(img::Images.AbstractImage{T,2}, sigma::Real=1)

    # Extract 2D float array in y x spatial order
    (dimg, prop) = floatyx(img)
    
    cimg = coherence(dimg, sigma)
    return Images.Image(cimg, prop)
end


#----------------------------------------------------------------------
"""
structuretensor - Compute structure tensor values over an image
```
Usage:  (Ix2, Iy2, Ixy) = structuretensor(img, sigma=1)

Arguments: 
           img   - Image to be processed  ::Array{T,2}
           sigma - Standard deviation of Gaussian summation window. Typical
                   values to use might be 1-3. Default is 1.
Returns:
           Ix2   - 2nd moment of gradients in x
           Iy2   - 2nd moment of gradients in y
           Ixy   - 2nd moment of gradients in xy
```
See also: shi_tomasi(), harris(), noble(), coherence(), derivative5()
"""
# April 2016

function structuretensor{T<:Real}(img::Array{T,2}, sigma::Real=1)

    # Convert to float if needed
    if isa(img[1,1], Integer)
        fimg = float(img) 
    else
        fimg = img
    end

    # Compute derivatives and elements of the structure tensor.
    (Ix, Iy) = derivative5(fimg, ("x", "y"))
    Ix2 = Images.imfilter_gaussian(Ix.^2,  [sigma, sigma])
    Iy2 = Images.imfilter_gaussian(Iy.^2,  [sigma, sigma])    
    Ixy = Images.imfilter_gaussian(Ix.*Iy, [sigma, sigma])    

    return Ix2, Iy2, Ixy
end


#----------------------------------------------------------------------
"""
hessianfeatures  - Computes determiant of hessian features in an image.
```
Usage: hdet = hessianfeatures(img, sigma)

Arguments:
            img      - Greyscale image to be processed.
                       ::Array{T,2} or ::AbstractImage{T,2} 
            sigma    - Defines smoothing scale.

Returns:    hdet     - Matrix of determinants of Hessian
```
The local maxima of hdet tend to mark the centres of dark or light blobs.
However, the point that gets localised can be dependent on scale.  If the
blobs you wish to detect are large you will need to use a value of sigma that
is comparable in magnitude.

The local minima of hdet is useful for marking the intersection points of a
camera calibration checkerbaord pattern.  These saddle features seem to be
more stable under scale variations.

For example to get the 100 strongest saddle features in image `img` use:
```
 > hdet = hessianfeatures(img, 1)      # sigma = 1
 > (r, c) = nonmaxsuppts(-hdet, N=100)
```
If the input is an AbstractImage cimg will be an AbstractImage with
property spatialorder set to ["y","x"].

See also: harris(), noble(), shi_tomasi(), derivative5(), nonmaxsuppts()
"""
function  hessianfeatures{T<:Real}(img::Array{T,2}, sigma::Real=1)

    if sigma > 0    # Convolve with Gaussian at desired sigma
        Gimg = Images.imfilter_gaussian(img,  [sigma, sigma])
    else           # No smoothing
        Gimg = img
        sigma = 1  # Needed for normalisation later
    end
    
    # Take 2nd derivatives in x and y
    (Lxx, Lxy, Lyy) = derivative5(Gimg, ("xx", "xy", "yy"))
    
    # Apply normalizing scaling factor of sigma^2 to 2nd derivatives        
    Lxx = Lxx*sigma^2
    Lyy = Lyy*sigma^2
    Lxy = Lxy*sigma^2
    
    # Determinant
    return hdet = Lxx.*Lyy - Lxy.^2
end

# Case for img::AbstractImage 
function hessianfeatures{T}(img::Images.AbstractImage{T,2}, sigma::Real=1)

    # Extract 2D float array in y x spatial order
    (dimg, prop) = floatyx(img)
    
    hdet = hessianfeatures(dimg, sigma)
    return Images.Image(cimg, prop)
end


#----------------------------------------------------------------------
"""
fastradial - Loy and Zelinski's fast radial feature detector
```
Usage: (S, So) = fastradial(img, radii, alpha=2, beta=0)

Arguments:
           img   - Image to be analysed
                   ::Array{T,2} or ::AbstractImage{T,2} 
           radii - Vector of integer radius values to be processed
                   suggested radii might be [1, 3, 5]
           alpha - Radial strictness parameter.
                   1 - slack, accepts features with bilateral symmetry.
                   2 - a reasonable compromise.
                   3 - strict, only accepts radial symmetry.
                       ... and you can go higher
           beta  - Gradient threshold.  Gradients below this threshold do
                   not contribute to symmetry measure, defaults to 0.

Returns    S     - Symmetry map.  Bright points with high symmetry are
                   marked with large positive values. Dark points of
                   high symmetry marked with large -ve values.
           So    - Symmetry map based on orientation only.
```
To localize points use nonmaxsuppts() on S, -S or abs(S) depending on
what you are seeking to find.

A good algorithm for detecting eyes in face images.

If the input is an AbstractImage cimg will be an AbstractImage with
property spatialorder set to ["y","x"].

Reference:  
Loy, G.  Zelinsky, A.  Fast radial symmetry for detecting points of
interest.  IEEE PAMI, Vol. 25, No. 8, August 2003. pp 959-973.
"""

# November 2004  - original version
# December 2009  - Gradient threshold added + minor code cleanup
# July     2010  - Gradients computed via Farid and Simoncelli's 5 tap
#                  derivative filters
# August   2015  - Ported to Julia

function fastradial{T<:Real}(img::Array{T,2}, radii::Vector, alpha::Real=2, beta::Real=0)
    
    if !isa(radii, Vector{Int}) && !isa(radii, Int) || minimum(radii) < 1
        error("radii must be a vector of integers and > 1")
    end
    
    (rows,cols)=size(img)
    
    # Compute derivatives in x and y via Farid and Simoncelli's 5 tap
    # derivative filters
    (imgx, imgy) = derivative5(img, ("x", "y"))
    mag = sqrt(imgx.^2 + imgy.^2)+eps() # (+eps to avoid division by 0)
    
    # Normalise gradient values so that [imgx imgy] form unit 
    # direction vectors.
    imgx = imgx./mag   
    imgy = imgy./mag
    
    S = zeros(rows,cols)  # Symmetry matrix
    So = zeros(rows,cols) # Orientation only symmetry matrix    
    
    for n in radii
        M = zeros(rows,cols)  # Magnitude projection image
        O = zeros(rows,cols)  # Orientation projection image

        # Form the orientation and magnitude projection matrices
        for r = 1:rows
            for c = 1:cols
                if mag[r,c] > beta

                    # Coordinates of 'positively' and 'negatively' affected pixels
                    posx = c + round(Int, n*imgx[r,c])
                    posy = r + round(Int, n*imgy[r,c])

                    negx = c - round(Int, n*imgx[r,c])
                    negy = r - round(Int, n*imgy[r,c])

                    # clamp to image limits
                    posx = max(min(posx,cols),1)
                    posy = max(min(posy,rows),1)
                    negx = max(min(negx,cols),1)
                    negy = max(min(negy,rows),1)

                    O[posy, posx] += 1
                    O[negy, negx] -= 1
                    
                    M[posy, posx] += mag[r,c]
                    M[negy, negx] -= mag[r,c]
                end
            end
        end
        
        # Clamp Orientation projection matrix values to a maximum of 
        # +/-kappa,  but first set the normalization parameter kappa to the
        # values suggested by Loy and Zelinski
        if n == 1
            kappa = 8
        else 
            kappa = 9.9
        end
        
        O[O .>  kappa] =  kappa  
        O[O .< -kappa] = -kappa  
        
        # Unsmoothed symmetry measure at this radius value
        F = M./kappa .* (abs(O)/kappa).^alpha
        Fo = sign(O) .* (abs(O)/kappa).^alpha   # Orientation only based measure
        
        # Smooth and spread the symmetry measure with a Gaussian
        # proportional to n.  Also scale the smoothed result by n so
        # that large scales do not lose their relative weighting.
        # Note we use the function gaussfilt because it seems that
        # small sigma values (< 1) cause problems for
        # Images.imfilt_gaussian()
        # ** But if sigma < 1 do we need to smooth anyway? **
        S  +=  gaussfilt(F,  0.25*n) * n
        So +=  gaussfilt(Fo, 0.25*n) * n        
        
    end  # for each radius
    
    S  = S /length(radii)  # Average
    So = So/length(radii) 

    return S, So
end

# Case for img::AbstractImage 
function fastradial{T}(img::Images.AbstractImage{T,2}, radii::Vector, alpha::Real=2, beta::Real=0)

    # Extract 2D float array in y x spatial order
    (dimg, prop) = floatyx(img)
    
    (S, So) = fastradial(dimg, radii, alpha, beta)
    return Images.Image(S, prop), Images.Image(So, prop)
end



