export fitlinedemo, fitplanedemo
export fitfunddemo, fithomogdemo
export hline, plotcamera
export plot_briefcoords, plot_nonmaxsuppts


# This function is used to not have the extensions have to overwrite the dfferent
# plotting fucntions, which would raise a precompilation warning.
# This also allows implementing a smarter way to select the correct plutting
# backend in case multiple backends are loaded (I would however doubt this 
# would ever happen).
# If another extension module is implemented (Makie for example), implement
# checking for it here. 
function get_plot_backend()
    ext = Base.get_extension(@__MODULE__, :ImageProjectiveGeometryPyPlotExt)
    if ext === nothing
        error("PyPlot extension not loaded...")
    else
        ext
    end
end

function plot_briefcoords(S, nPairs, rc)
    ext = get_plot_backend()
    return ext.plot_briefcoords(S, nPairs, rc)
end

function plot_nonmaxsuppts(fig, img, c, r, cols, rows)
    ext = get_plot_backend()
    return ext.plot_nonmaxsuppts(fig, img, c, r, cols, rows)
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
    ext = get_plot_backend()
    return ext.hline(p1i, p2i, linestyle)
end

# Case when homogeneous line is supplied
function hline(li::Vector, linestyle::String="b-")
    ext = get_plot_backend()
    return ext.hline(li, linestyle)
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
function plotcamera(Ci, l; col=[0,0,1], plotCamPath=false, fig=nothing)
    ext = get_plot_backend()
    return ext.plotcamera(Ci, l, col, plotCamPath, fig)
end

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
    ext = get_plot_backend()
    return ext.fitlinedemo(outliers, sigma, t, feedback)
end

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
    ext = get_plot_backend()
    return ext.fitplanedemo(outliers, sigma, t, feedback)
end

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
    ext = get_plot_backend()
    return ext.fitfunddemo(img1, img2)
end

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
    ext = get_plot_backend()
    return ext.fithomogdemo(img1, img2)
end

