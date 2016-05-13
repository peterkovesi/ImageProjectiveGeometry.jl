Corner Features Function Reference
==================================

## Index

* [shi_tomasi](#shi_tomasi) - Shi - Tomasi corner detector.
* [harris](#harris) - Harris corner detector.
* [noble](#noble) - Noble's variant of the Harris corner detector.
* [coherence](#coherence) - Compute image coherence from structure tensor
* [structuretensor](#structuretensor) - Compute structure tensor values over an image
* [hessianfeatures](#hessianfeatures)  - Computes determiant of hessian features in an image.
* [fastradial](#fastradial) - Loy and Zelinski's fast radial feature detector.

___________


## shi_tomasi 

Shi - Tomasi corner detector

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


## harris 

Harris corner detector.

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
to look at the values within `cimg` to determine the appropriate
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


## noble 

Noble's variant of the Harris corner detector.

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


## coherence 

Compute image coherence from structure tensor

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

If the input is an AbstractImage cimg will be an AbstractImage with
property spatialorder set to ["y","x"].

A small value indicates that the eigenvalues are similar in magnitude
and hence the local structure has a strong 2D component.

See also: structuretensor()


## structuretensor 

Compute structure tensor values over an image

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


## hessianfeatures 

Computes determiant of hessian features in an image.

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


## fastradial 

Loy and Zelinski's fast radial feature detector.

```
Usage: (S, So) = fastradial(img, radii, alpha, beta)

Arguments:
           img   - Image to be analysed
                   ::Array{T,2} or ::AbstractImage{T,2} 
           radii - Array of integer radius values to be processed
                   suggested radii might be [1 3 5]
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

