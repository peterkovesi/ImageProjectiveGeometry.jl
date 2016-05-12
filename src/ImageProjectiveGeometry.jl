#=----------------------------------------------------------------------------

Image Projective Geometry

Functions supporting projective geometry for computer vision.

Copyright (c) 2016 Peter Kovesi
pk@peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK February 2016


----------------------------------------------------------------------------=#
"""
**ImageProjectiveGeometry**

Functions supporting projective geometry for computer vision.

Peter Kovesi  

[peterkovesi.com](http://peterkovesi.com)

"""
module ImageProjectiveGeometry

include("projective.jl")
include("transforms3d.jl")
include("ransac.jl")
include("cornerfeatures.jl")
include("utilities.jl")
include("ransacdemo.jl")

end  # module

