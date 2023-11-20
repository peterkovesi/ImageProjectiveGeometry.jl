
function get_plot_backend()
    if !isdefined(Base, :get_extension)
        error("PyPlot extension not loaded...")
    end

    ext = Base.get_extension(@__MODULE__, :ImageProjectiveGeometryPyPlotExt)
    # TODO: add different plotting backends
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

function fitlinedemo(outliers, sigma, t::Real, feedback::Bool = false)
    ext = get_plot_backend()
    return ext.fitlinedemo(outliers, sigma, t, feedback)
end

function fitplanedemo(outliers, sigma, t, feedback::Bool = false)
    ext = get_plot_backend()
    return ext.fitplanedemo(outliers, sigma, t, feedback)
end

function fitfunddemo(img1=[], img2=[])
    ext = get_plot_backend()
    return ext.fitfunddemo(img1, img2)
end

function fithomogdemo(img1=[], img2=[])
    ext = get_plot_backend()
    return ext.fithomogdemo(img1, img2)
end
