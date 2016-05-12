function dilate1d_old{T<:Real}(f::Array{T,1}, ki::Real)
    
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
    padval = typemin(typeof(f[1]));
    fpad = [repmat([padval], rad1); f[:]; repmat([padval], rad2+extrapad)]

    N = length(fpad)
    g=copy(fpad)
    h=copy(fpad)
    
    # Generate the intermediate arrays
    for o = 0:k:(N-k)
        for n = 2:k
            g[o+n] = max(g[o+n], g[o+n-1])
        end
        
        for n = (k-1):-1:1
            h[o+n] = max(h[o+n], h[o+n+1])
        end            
    end
    
    # Combine the intermediate arrays g and h to obtain the dilation.  Note
    # that for even sized structuring elements there is a wasted max()
    # operation every k steps.  However in the interests of keeping the code
    # clean and simple this is accepted.
    for n = rad1+1:rad1+Nf
        df[n-rad1] = max(g[n+rad2], h[n-rad1])
    end

    return df
end




function imdilate{T<:Real}(img::Array{T,2}, se)
   
    (rows,cols) = size(img)
    (sr,sc) = size(se)
    roff = round(Int, (sr-1)/2)
    coff = round(Int, (sc-1)/2)
    dimg = zeros(size(img))

    for r = 1:rows-sr, c = 1:cols-sc
        dimg[r+roff, c+coff] = maximum(img[r:r+sr-1,  c:c+sc-1].*se)
    end
    
    return dimg
end



#------------------------------------------------------------------------
"""
imerode - Image morpholgical erosion
```
Usage:   dimg = imerode(img, se)

Arguments:
          img - Image to be eroded.  May be greyscale or binary
                ::Array{T,2}
           se - Structuring element defined by a 2D array

Returns:
         dimg - eroded image
```
See also: imdilate(), circularstruct()
"""

# Inefficient code!  Only handles square structuring elements
# April 2016

function imerode{T<:Real}(img::Array{T,2}, se)
   
    (rows,cols) = size(img)
    (sr,sc) = size(se)
    roff = round(Int32, (sr-1)/2)
    coff = round(Int32, (sc-1)/2)
    dimg = zeros(size(img))

    for r = 1:rows-sr, c = 1:cols-sc
        dimg[r+roff, c+coff] = minimum(img[r:r+sr-1,  c:c+sc-1])
    end
    
    return dimg
end



#function imdilate{T<:Real}(img::Union{Array{T,2}, BitArray{2}}, seType::ASCIIString, seSize)
function imdilate(img, seType::ASCIIString, seSize)

    if ndims(img) == 3
        error("Image must be binary or greyscale")
    end

    (rows,cols) = size(img)
    
    # Rectangular structuring element    
    if uppercase(seType[1:3]) == "REC"
        if length(seSize) == 1
            k = round(Int, [seSize, seSize])
        else
            k = round(Int, seSize)
        end
        
        dimg = copy(img)
        
        for c = 1:cols
            dimg[:,c] = dilate1d(img[:,c], k[1])
        end
        
        for r = 1:rows
            dimg[r,:] = dilate1d(vec(dimg[r,:]), k[2]) # 0.4
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
        dimg = imdilate(img, "rect", k)
        
        # Extract diagonal lines of pixels from dimg and perform 1d dilation on
        # them
        
        # NE lines, emanating from left edge of image
        for r = 2:rows
            cmax = min(r,cols)
            l = zeros(cmax)
            for c = 1:cmax
                l[c] = dimg[r-c+1, c]
            end
            
            dl = dilate1d(l, dk)
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
            
            dl = dilate1d(l, dk) 
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
            
            dl = dilate1d(l, dk)
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
            
            dl = dilate1d(l, dk)
            for r = 1:rmax
                dimg[r,c+r-1] = dl[r]
            end        
        end
        
    else
        error("Structure element type must be 'rect' or 'oct' ")
    end
    
    return dimg
end

