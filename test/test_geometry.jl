# Test geometry.jl

println("testing geometry")

tol = 1e-8

# ray2raydist
p1 = [0,0,0]
p2 = [0,0,1]
v1 = [1,0,0]
v2 = [0,1,0]
d = ray2raydist(p1,v1,p2,v2) # basic test
@test abs(d-1) < tol

d = ray2raydist(p1,v1,p2,v1) # Test parallel rays
@test abs(d-1) < tol

d = ray2raydist(p1,v1,p1,v2) # Test intersecting rays
@test abs(d) < tol
d = ray2raydist(p1,v1,p1,v1) # Test identical rays
@test abs(d) < tol

# circleintersectionpts
c1 = [-1,0]
r1 = 2
c2 = [1,0]
r2 = 2
(i1, i2) = circleintersectionpts(c1, r1, c2, r2, "lr")
@test abs(i1[1] < eps())
@test abs(i1[2] - sqrt(3)) < eps()
@test abs(i2[1] < eps())
@test abs(i2[2] + sqrt(3)) < eps()

# Check left and right solutions
ileft = circleintersectionpts(c1, r1, c2, r2, "l")
@test all(abs(ileft - i1) .< eps())
iright = circleintersectionpts(c1, r1, c2, r2, "r")
@test all(abs(iright - i2) .< eps())

# no intersection
(i1, i2) = circleintersectionpts(c1, r1, c2+10, r2, "lr")
@test isempty(i1)
@test isempty(i2)

# circles just touching
r1 = 1
r2 = 1
il = circleintersectionpts(c1, r1, c2, r2, "l")
@test all(abs(il) .< eps())

# circle of zero radius
c2 = [0,0]
r2 = 0
il = circleintersectionpts(c1, r1, c2, r2, "l")
@test all(abs(il) .< eps())



# rectintersect
r1 = [0, 0, 2, 4]
r2 = [0, 0, 1, 1]
@test rectintersect(r1, r2) # r2 inside r1
r2 = [2, 4, 1, 1]
@test rectintersect(r1, r2) # r2 touching r1
r2 = [10, 20, 1, 1]
@test !rectintersect(r1, r2) # r2 outside r1
@test rectintersect(r1, r1)  # r2 same as r1

# pointinconvexpoly
poly = [0  1  2  2  0
        0  0  0  1  1]

@test pointinconvexpoly([2,1], poly) == 0      # on boundary
@test pointinconvexpoly([0.5,0.5], poly) == 1  # inside
@test pointinconvexpoly([10,10], poly) == -1   # outside

# isconvexpoly
@test isconvexpoly(poly)  # convex
poly2 = [0    1  2  3  2  0
         0  -10  0  0  1  1]
@test !isconvexpoly(poly2)  # non-convex


# convexpolyintersect
poly1 = [0  2  3  0 
         0  0  1  1]

poly2 = [0  2  1
         0  0  1]
@test convexpolyintersect(poly1, poly2)      # intersect
@test !convexpolyintersect(poly1, poly2+10)  # not intersect

poly3 = [-2  0  -1
          0  0   1]
@test convexpolyintersect(poly2, poly3)      # touching
@test convexpolyintersect(poly1, poly1)      # identical

