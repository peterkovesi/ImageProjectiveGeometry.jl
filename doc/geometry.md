Geometry Function Reference
==================================

## Index

* [ray2raydist](#ray2raydist) -  Minimum distance between two 3D rays.
* [circleintersectionpts](#circleintersectionpts) - Finds intersections of two circles.
* [rectintersect](#rectintersect) - Test if rectangles intersect.
* [pointinconvexpoly](#pointinconvexpoly) - Determine if a 2D point is within a convex polygon.
* [convexpolyintersect](#convexpolyintersect) - Test if convex polygons intersect.
* [isconvexpoly](#isconvexpoly) - Test if 2D polygon is convex.

_____________________

## ray2raydist 

 Minimum distance between two 3D rays.

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


## circleintersectionpts 

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


## rectintersect 

Test if rectangles intersect

```
Usage:  intersect = rectintersect(r1, r2) 

Arguments:  r1, r2 - The two rectangles to be tested where the rectangles
                     are defined using 4-vectors with the elements
                     defining  [left, bottom, width, height]

Returns:             Boolean result.
```

See also: convexpolyintersect()


## pointinconvexpoly 

Determine if a 2D point is within a convex polygon.

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


## convexpolyintersect 

Test if convex polygons intersect.

```
Usage:  intersect = convexpolyintersect(p1, p2) 

Arguments:  p1, p2 - The two convex polygons to be tested where the
                     polygons are defined as 2xN arrays of vertex
                     coordinates in sequence around the polygon.

Returns:             Boolean result
```

Note: There is no check to see if the polygons are indeed convex, use
isconvexpoly() if needed

See also: isconvexpoly(), rectintersect()


## isconvexpoly 

Test if 2D polygon is convex.

```
Usage: v = isconvexpoly(p)

Argument:  p - 2D polygon are defined as 2xN array of vertices.

Returns:   true if polygon is convex.

```

See also: convexpolyintersect(), pointinconvexpoly()