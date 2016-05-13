Ransacdemo Function Reference
=============================

## Index

* [fitlinedemo](#fitlinedemo) - Demonstrates RANSAC line fitting.
* [fitplanedemo](#fitplanedemo) - Demonstrates RANSAC plane fitting.
* [fitfunddemo](#fitfunddemo) - Example of fundamental matrix computation.
* [fithomogdemo](#fithomogdemo) - Example of finding a homography.

____________________________________________

## fitlinedemo 

Demonstrates RANSAC line fitting.

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


## fitplanedemo 

Demonstrates RANSAC plane fitting.

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


## fitfunddemo 

Example of fundamental matrix computation.

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


## fithomogdemo 

Example of finding a homography.

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