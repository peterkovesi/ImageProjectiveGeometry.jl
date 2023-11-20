module ImageProjectiveGeometryPyPlotExt
using LinearAlgebra
using PyPlot

isdefined(Base, :get_extension) ? (using ImageProjectiveGeometry) : (using ..ImageProjectiveGeometry)

include("ransacdemo.jl")

function plot_briefcoords(S, nPairs, rc)
    figure(200); clf
        R = (S.-1)/2 .+ 1
        axis([-R, R, -R, R])

        for n = 1:nPairs
            plot(rc[2, 2*n-1:2*n], rc[1, 2*n-1:2*n], color = rand(3), 
                 linewidth = 3)
        end

        axis("equal")

        # Determine distances between pairs and display histogram
        sdist = zeros(nPairs)
        for n = 1:nPairs
            sdist[n] = norm(rc[:,2*n-1] - rc[:,2*n])
        end
        
        (e, counts) = hist(sdist, 0:.5:S)
        e = e .+ 0.5
        figure(30); clf()
        PyPlot.bar(x=e, height = counts[2:end], width=e[2]-e[1])
        title("Histogram of distances between pairs")
end

function plot_nonmaxsuppts(fig, img, c, r, cols, rows)
    # If an image has been supplied display it and overlay corners.
    if !isempty(img)
        print("Plotting corners")
        PyPlot.figure(fig)
        PyPlot.clf()
        PyPlot.imshow(img)
        PyPlot.set_cmap(PyPlot.ColorMap("gray"))
        PyPlot.plot(c,r,"r+")
        PyPlot.axis([1,cols,rows,1])
        PyPlot.title("Corners detected")
    end
end


"""
hline - Plot a 2D line defined in homogeneous coordinates.

Function for ploting 2D homogeneous lines defined by 2 points
or a line defined by a single homogeneous vector
```
Usage 1:   hline(p1,p2, linestyle="b-")
Arguments: p1, p2 -  Two 3-vectors defining points in homogeneous
                     coordinates

Usage 2:    hline(l, linestyle="b-")
Argument:   l - A 3-vector defining a line in homogeneous coordinates

```
Note that in the case where a homogeneous line is supplied as the
argument the extent of the line drawn depends on the current axis
limits.  This will require you to set the desired limits with a call
to PyPlot.axis() prior to calling this function.

"""
function hline(p1i::Vector, p2i::Vector, linestyle::String="b-")
    # Case when 2 homogeneous points are supplied

    p1 = p1i./p1i[3]    # make sure homogeneous points lie in z=1 plane
    p2 = p2i./p2i[3]

#    hold(true)
    plot([p1[1], p2[1]], [p1[2], p2[2]], linestyle);
end

# Case when homogeneous line is supplied
function hline(li::Vector, linestyle::String="b-")

    l = li./li[3]   # ensure line in z = 1 plane (not needed?)

    if abs(l[1]) > abs(l[2])   # line is more vertical
        p1 = hcross(l, [0, -1, PyPlot.ylim()[1]])
        p2 = hcross(l, [0, -1, PyPlot.ylim()[2]])

    else                       # line more horizontal
        p1 = hcross(l, [-1, 0, PyPlot.xlim()[1]])
        p2 = hcross(l, [-1, 0, PyPlot.xlim()[2]])
    end

#    hold(true)
    plot([p1[1], p2[1]], [p1[2], p2[2]], linestyle);
end


"""
plotcamera - Plots graphical representation of camera(s) showing pose.

```
Usage: plotcamera(C, l; col=[0,0,1], plotCamPath=false, fig=nothing)

Arguments:
           C - Camera structure (or array of Camera structures)
           l - The length of the sides of the rectangular cone indicating
               the camera's field of view.
Keyword Arguments:
         col - Optional three element vector specifying the RGB colour to
               use. Defaults to blue.
 plotCamPath - Optional flag true/false to plot line joining camera centre
               positions. If omitted or empty defaults to false.
         fig - Optional figure number to be used. If not specified a new
               figure is created.
```

The function plots into the current figure a graphical representation of one
or more cameras showing their pose.  This consists of a rectangular cone,
with its vertex at the camera centre, indicating the camera's field of view.
The camera's coordinate X and Y axes are also plotted at the camera centre.

See also: Camera
"""
function plotcamera(Ci, l, col=[0,0,1], plotCamPath=false, fig=nothing)

    if isa(Ci, Array)
        C = Ci
    else
        C = [Ci]
    end

    figure(fig)
    for i = 1:eachindex(C)

        if C[i].rows == 0 || C[i].cols == 0
            @warn("Camera rows and cols not specified")
            continue
        end

        f = C[i].fx  # Use fx as the focal length

        if i > 1 && plotCamPath
            plot3D([C[i-1].P[1], C[i].P[1]],
                   [C[i-1].P(2), C[i].P[2]],
                   [C[i-1].P(3), C[i].P[3]])
        end

        # Construct transform from camera coordinates to world coords
        Tw_c = [C[i].Rc_w'  C[i].P
                 0 0 0        1  ]

        # Generate the 4 viewing rays that emanate from the principal point and
        # pass through the corners of the image.
        ray = zeros(3,4)
        ray[:,1] = [-C[i].cols/2, -C[i].rows/2, f]
        ray[:,2] = [ C[i].cols/2, -C[i].rows/2, f]
        ray[:,3] = [ C[i].cols/2,  C[i].rows/2, f]
        ray[:,4] = [-C[i].cols/2,  C[i].rows/2, f]

        # Scale rays to distance l from the focal plane and make homogeneous
        ray = makehomogeneous(ray*l/f)
        ray = Tw_c*ray                 # Transform to world coords

        for n = 1:4             # Draw the rays
            plot3D([C[i].P[1], ray[1,n]],
                   [C[i].P[2], ray[2,n]],
                   [C[i].P[3], ray[3,n]],
                   color=col)
        end

        # Draw rectangle joining ends of rays
        plot3D([ray[1,1], ray[1,2], ray[1,3], ray[1,4], ray[1,1]],
               [ray[2,1], ray[2,2], ray[2,3], ray[2,4], ray[2,1]],
               [ray[3,1], ray[3,2], ray[3,3], ray[3,4], ray[3,1]],
               color=col)

        # Draw and label axes
        X = Tw_c[1:3,1]*l .+ C[i].P
        Y = Tw_c[1:3,2]*l .+ C[i].P
        Z = Tw_c[1:3,3]*l .+ C[i].P

        plot3D([C[i].P[1], X[1,1]], [C[i].P[2], X[2,1]], [C[i].P[3], X[3,1]],
               color=col)
        plot3D([C[i].P[1], Y[1,1]], [C[i].P[2], Y[2,1]], [C[i].P[3], Y[3,1]],
               color=col)
        #    plot3D([C[i].P[1], Z(1,1)], [C[i].P[2], Z(2,1)], [C[i].P[3], Z(3,1)],...
        #           color=col)
        text3D(X[1], X[2], X[3], "X", color=col)
        text3D(Y[1], Y[2], Y[3], "Y", color=col)

    end

end

end # module

