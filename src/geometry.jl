#=--------------------------------------------------------------------

geometry - Functions for some basic geometric operations

Part of the ImageProjectiveGeometry Module

Copyright (c) 2016 Peter Kovesi
pk@peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK June 2016

---------------------------------------------------------------------=#

export ray2raydist, circleintersectionpts
export rectintersect, pointinconvexpoly, isconvexpoly
export convexpolyintersect

#--------------------------------------------------------------------------
"""
ray2raydist -  Minimum distance between two 3D rays

```
Usage: d = ray2raydist(p1, v1, p2, v2)

Arguments:
      p1, p2 - 3D points that lie on rays 1 and 2.
      v1, v2 - 3D vectors defining the direction of each ray.

Returns:
           d - The minimum distance between the rays.
```

Each ray is defined by a point on the ray and a vector giving the direction
of the ray.  Thus a point on ray 1 will be given by  p1 + alpha*v1  where
alpha is some scalar.
"""

function ray2raydist(p1::Vector, v1::Vector, p2::Vector, v2::Vector)
    
    # Get vector perpendicular to both rays
    n = cross(v1, v2)
    nnorm = sqrt(vecdot(n,n))

    # Check if lines are parallel. If so, form a vector perpendicular to v1
    # that is within the plane formed by the parallel rays.
    if nnorm < eps()
        n = cross(v1, p1-p2)  # Vector perpendicular to plane formed by pair
                              # of rays.
        n = cross(v1, n)      # Vector perpendicular to v1 within the plane
                              # formed by the 2 rays.
        nnorm = sqrt(vecdot(n,n))

        if nnorm < eps()      # p1 and p2 must be coincident
            return  0.0
        end
    end                     
    
    n /= nnorm         # Unit vector
    
    # The vector joining the two closest points on the rays is:
    #    d*n = p2 + beta*v2 - (p1 + alpha*v1) 
    # for some unknown values of alpha and beta.   
    #
    # Take dot product of n on both sides
    #    d*n.n = p2.n + beta*v2.n - p1.n - alpha*v1.n
    #
    # as n is prependicular to v1 and v2 this reduces to
    #   d = (p2 - p1).n
    return abs(vecdot((p2 - p1), n))
end

#--------------------------------------------------------------------------
"""
circleintersectionpts - Finds intersections of two circles.

Function computes the intersection points between two circles given
their centres and radii.

```
Usage: (i1, i2) = circleintersectionpts(c1, r1, c2, r2, lr)
             i1 = circleintersectionpts(c1, r1, c2, r2, lr)

Arguments:
        c1, c2 - 2-vectors specifying the centres of the two circles.
        r1, r2 - Radii of the two circles.
            lr - Optional string specifying what solutions are wanted:
                "l" for the solution to the left of the line from c1 to c2
                "r" for the solution to the right of the line from c1 to c2
                "lr" if both solutions are wanted (default).

Returns:
     (i1, i2) - Tuple of 2D intersection points if both solutions requested.
          i1  - Single 2D intersection point if a single solution requested.
```

If no solution exists empty vectors are returned for i1, i2.  If the
two circles are identical there are an infinite number of solutions.
In this case we simply choose to subtract/add r1 in the x direction to
c1 to obtain the 'left' and 'right' solutions.
"""

function circleintersectionpts(c1i::Vector, r1::Real, c2i::Vector, r2::Real, lr="lr")

    c1 = float(c1i) # Ensure c1 and c2 are floats for type stability
    c2 = float(c2i)

    maxmag = maximum(abs.([c1; c2; r1; r2]))  # Maximum value in data input
    tol = (maxmag+1)*eps()              # Tolerance used for deciding whether values are equal
    
    bv = c2 - c1         # Baseline vector from c1 to c2
    b = norm(bv)         # The distance between the centres
    
    # Trap case of baseline of zero length.  If r1 == r2 we have a
    # valid geometric situation, but there are an infinite number of
    # solutions.  Here we simply choose to subtract/add r1 in the x
    # direction to c1 to obtain the 'left' and 'right' solutions.
    if b < tol && abs(r1-r2) < tol  
        i1 = c1 - [r1, 0]
        i2 = c1 + [r1, 0]
    else  
        bv = bv/b             # Normalised baseline vector.
        bvp = [-bv[2], bv[1]] # Vector perpendicular to baseline
        
        # Trap the degenerate cases where one of the radii are zero, or nearly zero
        if r1 < tol && abs(b-r2) < tol
            i1 = c1
            i2 = c1
        elseif r2 < tol && abs(b-r1) < tol
            i1 = c2
            i2 = c2
            
        # Check triangle inequality
        elseif b > (r1+r2) || r1 > (b+r2) || r2 > (b+r1)
            i1 = zeros(eltype(c1), 0)
            i2 = zeros(eltype(c1), 0)
            
        else   # Normal solution
            cosR2 = (b^2 + r1^2 - r2^2)/(2*b*r1)
            sinR2 = sqrt(1-cosR2^2)
            
            i1 = c1 + r1*cosR2*bv + r1*sinR2*bvp   # 'left' solution
            i2 = c1 + r1*cosR2*bv - r1*sinR2*bvp   # and 'right solution
        end    
    end

    # Return requested solution(s)
    if lr == "lr"
        return (i1, i2)

    elseif lr == "l"
        return i1

    elseif lr == "r"
        return i2

    else
        error("illegal left/right solution request")
    end
    
end

#--------------------------------------------------------------------------
"""
rectintersect - Test if rectangles intersect

```
Usage:  intersect = rectintersect(r1, r2) 

Arguments:  r1, r2 - The two rectangles to be tested where the rectangles
                     are defined using 4-vectors with the elements
                     defining  [left, bottom, width, height]

Returns:             Boolean result
```

See also: convexpolyintersect()
"""

function rectintersect(r1::Vector, r2::Vector) 
    
    intersect = (r1[1] <= r2[1]+r2[3] && 
                 r2[1] <= r1[1]+r1[3] && 
                 r1[2] <= r2[2]+r2[4] && 
                 r2[2] <= r1[2]+r1[4])
end

#--------------------------------------------------------------------------
"""
pointinconvexpoly - Determine if a 2D point is within a convex polygon

```
Usage:  v = pointinconvexpoly(p, poly)

Arguments:   p - 2-vector specifying the point.
          poly - Convex polygon defined as a series of vertices in
                 sequence, clockwise or anticlockwise around the polygon as
                 a 2 x N array.

Returns:    v - +1 if within the polygon
                -1 if outside
                 0 if on the boundary
```

Note:  There is no check to see if the polygon is indeed convex, use
isconvexpoly() if needed

See also: isconvexpoly()
"""

# Strategy: Determine whether, in traveling from p to a vertex and
# then to the following vertex, we turn clockwise or anticlockwise.
# If for every vertex we turn consistently clockwise, or consistently
# anticlockwise we are inside the polygon.  If for one of the vertices
# we did not turn clockwise or anticlockwise then we must be on the
# boundary.

function pointinconvexpoly(p::Vector, poly::Array)
    
    (dim, N) = size(poly)
    
    if length(p) != 2 || dim != 2
        error("Data must be 2D")
    end
    
    # Determine whether, in traveling from p to a vertex and then to the
    # following vertex, we turn clockwise or anticlockwise
    c = zeros(Int, N)
    for n = 1:N-1
        c[n] = clockwise(p, view(poly,:,n), view(poly,:,n+1))
    end
    c[N] = clockwise(p, view(poly,:,N), view(poly,:,1))
    
    # If for every vertex we turn consistently clockwise, or
    # consistently anticlockwise we are inside the polygon.  If for
    # one of the vertices we did not turn clockwise or anticlockwise
    # then we must be on the boundary.
    if all(c .>= 0) || all(c .<= 0)
        if any(c .== 0)  # We are on the boundary
            v = 0
        else             # We are inside
            v = 1
        end
    else                 # We are outside
        v = -1
    end

    return v
end
    
#----------------------------------------------------------------------
# Determine whether, in traveling from p1 to p2 to p3 we turn clockwise or
# anticlockwise.  Returns +1 for clockwise, -1 for anticlockwise, and 0 for
# p1, p2, p3 in a straight line.
    
function clockwise(p1::AbstractArray, p2::AbstractArray, p3::AbstractArray)
    
    # Form vectors p1->p2 and p2->p3 with z component = 0, form cross product
    # if the resulting z value is -ve the vectors turn clockwise, if +ve
    # anticlockwise, and if 0 the points are in a line.
        
    return -sign((p2[1]-p1[1])*(p3[2]-p2[2]) - (p2[2]-p1[2])*(p3[1]-p2[1]))

end

#----------------------------------------------------------------------
"""
convexpolyintersect - Returns true if convex polygons intersect.

```
Usage:  intersect = convexpolyintersect(p1, p2) 

Arguments:  p1, p2 - The two convex polygons to be tested where the
                     polygons are defined as 2xN arrays of vertex
                     coordinates in sequence around the polygon.

Returns:             Boolean result.
```

Note: There is no check to see if the polygons are indeed convex, use
isconvexpoly() if needed

See also: isconvexpoly(), rectintersect()
"""

function convexpolyintersect(p1::Array, p2::Array)

    # Function to test if two 1D segments overlap.
    # Each segment is defined by a 2 element vector with 1st element
    # being the minimum x value and 2nd element the maximum x value.
    function  overlap(s1::Vector, s2::Vector)
        return  s2[2] >= s1[1] && s1[2] >= s2[1]
    end

    
    if size(p1,1) != 2 || size(p2,1) != 2
        error("Polygon data must be in 2xN arrays")
    end
    
    # Build a sequence of vertices that trace the edges of both
    # polygons.  We prepend the last vertex to each polygon list to
    # form a complete cycle.  We also concatenate the two polygon
    # arrays together, this adds an extra 'edge' that connects the
    # last vertex of the 1st polygon with the last vertex of the 2nd
    # polygon.  This allows us to run a single loop of tests below.
    p = [p1[:,end] p1  p2[:,end] p2]

    for n = 2:size(p,2)
        # Generate a 2x2 rotation matrix describing a frame with the y axis
        # aligned with the polygon side.
        y = p[:,n] - p[:,n-1]
        normy = sqrt(y[1]^2 + y[2]^2)
        if normy < eps()  # repeated vertex
            continue
        end
        y = y./normy
        x = [y[2], -y[1]]
        R = [x y]'
        
        # Rotate all vertices to frame defined by this side
        p1r = R*p1
        p2r = R*p2
        
        # Check whether x coordinates of rotated p1 vertices overlap
        # with x coordinates of rotated p2 vertices.  If they do not
        # overlap we have found a separating axis indicating that the
        # polygons do not intersect.
        if !overlap([minimum(p1r[1,:]), maximum(p1r[1,:])], 
                    [minimum(p2r[1,:]), maximum(p2r[1,:])])
            return false
        end
    end

    return true        # No separating axis found
end

#--------------------------------------------------------------------------
"""
isconvexpoly - Test if 2D polygon is convex.

```
Usage: v = isconvexpoly(p)

Argument:  p - 2D polygon are defined as 2xN array of vertices.

Returns:   true if polygon is convex.

```

See also: convexpolyintersect(), pointinconvexpoly()
"""

function isconvexpoly(poly::Array)

    (dim, N) = size(poly)

    if dim != 2
        error("Polygon data must be a 2xN array")
    end

    # Strategy: Determine turning direction as one travels through
    # each set of three successive vertices.  If, for all vertices,
    # the turning direction is always clockwise or straight, or always
    # anti-clockwise or straight then the polygon is convex.

    c = zeros(Int, N)
    for n = 1:N-2
        c[n] = clockwise(view(poly,:,n), view(poly,:,n+1), view(poly,:,n+2))
    end
    c[N-1] = clockwise(view(poly,:,N-1), view(poly,:,N), view(poly,:,1))
    c[N] = clockwise(view(poly,:,N), view(poly,:,1), view(poly,:,2))

    if all(c .>= 0) || all(c .<= 0)
        return true
    else
        return false
    end
end
