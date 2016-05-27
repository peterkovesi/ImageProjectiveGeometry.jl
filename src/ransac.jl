#=--------------------------------------------------------------------

ransac - Functions for robust estimation using RANSAC


Copyright (c) 2016 Peter Kovesi
pk@peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK February 2016

---------------------------------------------------------------------=#

export ransac
export ransacfithomography, ransacfitfundmatrix, ransacfitaffinefundmatrix
export ransacfitplane, ransacfitline
export iscolinear, fitline2d, fitline3d, fitplane

#---------------------------------------------------------------------
"""
ransac - Robustly fits a model to data with the RANSAC algorithm
```
Usage:

(M, inliers) = ransac(x, fittingfn, distfn, degenfn s, t, feedback,
                      maxDataTrials, maxTrials)

Arguments:
    x         - Data sets to which we are seeking to fit a model M
                It is assumed that x is of size [d x Npts]
                where d is the dimensionality of the data and Npts is
                the number of data points.

    fittingfn - Reference to a function that fits a model to s
                data from x.  It is assumed that the function is of the
                form: 
                   M = fittingfn(x)
                Note it is possible that the fitting function can return
                multiple models (for example up to 3 fundamental matrices
                can be fitted to 7 pairs of matched points).  In this 
                case it is assumed that the fitting function returns an
                array of models.
                If this function cannot fit a model it should return M as
                an empty matrix.

    distfn    - Reference to a function that evaluates the
                distances from the model to data x.
                It is assumed that the function is of the form:
                   (inliers, M) = distfn(M, x, t)
                This function must evaluate the distances between points
                and the model returning the indices of elements in x that
                are inliers, that is, the points that are within distance
                't' of the model.  Additionally, if M is an array of
                possible models distfn() will return the model that has the
                most inliers.  If there is only one model this function
                must still copy the model to the output.  After this call M
                will be a single object representing only one model. 

    degenfn   - Reference to a function that determines whether a
                set of datapoints will produce a degenerate model.
                This is used to discard random samples that do not
                result in useful models.
                It is assumed that degenfn is a boolean function of
                the form: 
                   r = degenfn(x)
                It may be that you cannot devise a test for degeneracy in
                which case you should write a dummy function that always
                returns a value of true and rely on fittingfn() to return
                an empty model should the data set be degenerate.

    s         - The minimum number of samples from x required by
                fittingfn to fit a model.

    t         - The distance threshold between a data point and the model
                used to decide whether the point is an inlier or not.

    feedback  - An optional boolean flag. If set to true the trial count and the
                estimated total number of trials required is printed out at
                each step.  Defaults to false.

maxDataTrials - Maximum number of attempts to select a non-degenerate
                data set. This parameter is optional and defaults to 100.

    maxTrials - Maximum number of iterations. This parameter is optional and
                defaults to 1000.

    p         - Desired probability of choosing at least one sample
                free from outliers, defaults to 0.99

Returns:
    M         - The model having the greatest number of inliers.
    inliers   - An array of indices of the elements of x that were
                the inliers for the best model.
```
For an example of the use of this function see ransacfithomography() or
ransacfitplane() 
"""
#=
References:
   M.A. Fishler and  R.C. Boles. "Random sample concensus: A paradigm
   for model fitting with applications to image analysis and automated
   cartography". Comm. Assoc. Comp, Mach., Vol 24, No 6, pp 381-395, 1981

   Richard Hartley and Andrew Zisserman. "Multiple View Geometry in
   Computer Vision". pp 101-113. Cambridge University Press, 2001

May      2003 - Original version in MATLAB
February 2004 - Tidied up.
August   2005 - Specification of distfn changed to allow model fitter to
                return multiple models from which the best must be selected
Sept     2006 - Random selection of data points changed to ensure duplicate
                points are not selected.
December 2008 - Octave compatibility mods
June     2009 - Argument 'MaxTrials' corrected to 'maxTrials'!
February 2016 - Ported to Julia
=#

function ransac(x, fittingfn, distfn, degenfn, s, t, feedback = false, 
                maxDataTrials = 1000, maxTrials = 1000, p = 0.99)

    (rows, npts) = size(x)

    if npts < s
        error("Not enough data points to form a model")
    end
    
    bestM = NaN      # Sentinel value allowing detection of solution failure.
    trialcount = 0
    bestscore =  0
    bestinliers = []  # Ensure global to function
    bestM = []
    N = 1            # Dummy initialisation for number of trials.
    
    while N > trialcount
        
        # Select at random s datapoints to form a trial model, M.
        # In selecting these points we have to check that they are not in
        # a degenerate configuration.
        degenerate = true
        count = 1
        while degenerate
            # Generate s random indicies in the range 1..npts
            ind = randperm(npts)[1:s]

            # Test that these points are not a degenerate configuration.
            degenerate = degenfn(x[:,ind])

            if !degenerate
                # Fit model to this random selection of data points.
                # Note that M may represent a set of models that fit the data in
                # this case M will be a cell array of models
                M = fittingfn(x[:,ind])

                # Depending on your problem it might be that the only way you
                # can determine whether a data set is degenerate or not is to
                # try to fit a model and see if it succeeds.  If it fails we
                # reset degenerate to true.
                if isempty(M)
                    degenerate = true
                end
            end
            
            # Safeguard against being stuck in this loop forever
            count += 1
            if count > maxDataTrials
                warn("Unable to select a nondegenerate data set")
                break
            end
        end
        
        # Once we are out here we should have some kind of model...
        # Evaluate distances between points and model returning the indices
        # of elements in x that are inliers.  Additionally, if M is a cell
        # array of possible models 'distfn' will return the model that has
        # the most inliers.  After this call M will be a non-cell object
        # representing only one model.
        (inliers, M) = distfn(M, x, t)

        # Find the number of inliers to this model.
        ninliers = length(inliers)
        
        if ninliers > bestscore    # Largest set of inliers so far...
            bestscore = ninliers   # Record data for this model
            bestinliers = copy(inliers)
            bestM = copy(M)
            
            # Update estimate of N, the number of trials to ensure we pick,
            # with probability p, a data set with no outliers.
            fracinliers =  ninliers/npts
            pNoOutliers = 1 -  fracinliers^s
            pNoOutliers = max(eps(), pNoOutliers)  # Avoid division by -Inf
            pNoOutliers = min(1-eps(), pNoOutliers)# Avoid division by 0.
            N = log(1-p)/log(pNoOutliers)
        end
        
        trialcount += 1
        if feedback
            @printf("trial %d out of %d         \r",trialcount, ceil(N))
        end

        # Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
            @printf("RANSAC reached the maximum number of %d trials", maxTrials)
            break
        end
    end
    
    if feedback; @printf("\n"); end

    if !any(isnan(bestM))   # We got a solution
        M = bestM
        inliers = bestinliers
    else
        M = []
        inliers = []
        error("RANSAC was unable to find a useful solution")
    end

    return M, inliers
end    

#-----------------------------------------------------------------------
"""
ransacfithomography - Fits 2D homography using RANSAC
```
Usage:   (H, inliers) = ransacfithomography(x1, x2, t)

Arguments:
         x1  - 2xN or 3xN set of homogeneous points.  If the data is
               2xN it is assumed the homogeneous scale factor is 1.
         x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
         t   - The distance threshold between data point and the model
               used to decide whether a point is an inlier or not. 
               Note that point coordinates are normalised to that their
               mean distance from the origin is sqrt(2).  The value of
               t should be set relative to this, say in the range 
               0.001 - 0.01  
```
Note that it is assumed that the matching of x1 and x2 are putative and it
is expected that a percentage of matches will be wrong.
```
Returns:
         H       - The 3x3 homography such that x2 = H*x1.
         inliers - An array of indices of the elements of x1, x2 that were
                   the inliers for the best model.
```
See Also: ransac(), homography2d(), homography1d()
"""

#=
February 2004 - original version
July     2004 - error in denormalising corrected (thanks to Andrew Stein)
August   2005 - homogdist2d modified to fit new ransac specification.
=#

function ransacfithomography(x1i, x2i, t::Real)
    
    # Define some subfunctions

    #----------------------------------------------------------------------
    # Function to evaluate the symmetric transfer error of a homography with
    # respect to a set of matched points as needed by RANSAC.
    
    function homogdist2d(H, x, t)
        
        x1 = x[1:3,:]   # Extract x1 and x2 from x
        x2 = x[4:6,:]    
        
        # Calculate, in both directions, the transfered points    
        Hx1    = H*x1
        invHx2 = H\x2
        
        # Normalise so that the homogeneous scale parameter for all coordinates
        # is 1.
        x1     = hnormalise(x1)
        x2     = hnormalise(x2)     
        Hx1    = hnormalise(Hx1)
        invHx2 = hnormalise(invHx2) 
        
        d2 = sum((x1-invHx2).^2, 1)  + sum((x2-Hx1).^2, 1)
        inliers = find(abs(d2) .< t)    
        return inliers, H
    end    
    
    #----------------------------------------------------------------------
    # Function to determine if a set of 4 pairs of matched  points give rise
    # to a degeneracy in the calculation of a homography as needed by RANSAC.
    # This involves testing whether any 3 of the 4 points in each set is
    # colinear. 
    
    function isdegenerate(x)
        
        (rows, npts) = size(x)
        assert(npts == 4)

        x1 = x[1:3,:]    # Extract x1 and x2 from x
        x2 = x[4:6,:]    
        
        r = 
        iscolinear(x1[:,1],x1[:,2],x1[:,3]) ||
        iscolinear(x1[:,1],x1[:,2],x1[:,4]) ||
        iscolinear(x1[:,1],x1[:,3],x1[:,4]) ||
        iscolinear(x1[:,2],x1[:,3],x1[:,4]) ||
        iscolinear(x2[:,1],x2[:,2],x2[:,3]) ||
        iscolinear(x2[:,1],x2[:,2],x2[:,4]) ||
        iscolinear(x2[:,1],x2[:,3],x2[:,4]) ||
        iscolinear(x2[:,2],x2[:,3],x2[:,4])
        
        return r        
    end
    #----------------------------------------------------------

    # The main function

    if size(x1i) != size(x2i)
        error("Data sets x1 and x2 must have the same dimensions")
    end
    
    (rows,npts) = size(x1i)
    if !(rows == 2 || rows == 3)
        error("x1 and x2 must have 2 or 3 rows")
    end
    
    if npts < 4
        error("Must have at least 4 points to fit homography")
    end
    
    if rows == 2    # Pad data with homogeneous scale factor of 1
        x1 = [x1i; ones(1,npts)]
        x2 = [x2i; ones(1,npts)]        
    else
        x1 = copy(x1i)
        x2 = copy(x2i)
    end
        
    # Normalise each set of points so that the origin is at centroid and
    # mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    # scale parameter is 1.  Note that 'homography2d' will also call
    # 'normalise2dpts' but the code in 'ransac' that calls the distance
    # function will not - so it is best that we normalise beforehand.
    (x1, T1) = normalise2dpts(x1)
    (x2, T2) = normalise2dpts(x2)
    
    s = 4  # Minimum No of points needed to fit a homography.
    
    fittingfn = homography2d
    distfn    = homogdist2d
    degenfn   = isdegenerate
    # x1 and x2 are 'stacked' to create a 6xN array for ransac
    (H, inliers) = ransac([x1; x2], fittingfn, distfn, degenfn, s, t)
    
    # Now do a final least squares fit on the data points considered to
    # be inliers.
    H = homography2d(x1[:,inliers], x2[:,inliers])
    
    # Denormalise
    H = T2\H*T1    

    return H, inliers

end # of ransacfithomography


#-----------------------------------------------------------------------
"""
ransacfitfundmatrix - Fits fundamental matrix using RANSAC
```
Usage:   (F, inliers) = ransacfitfundmatrix(x1, x2, t)

Arguments:
         x1  - 2xN or 3xN set of homogeneous points.  If the data is
               2xN it is assumed the homogeneous scale factor is 1.
         x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
         t   - The distance threshold between data point and the model
               used to decide whether a point is an inlier or not. 
               Note that point coordinates are normalised to that their
               mean distance from the origin is sqrt(2).  The value of
               t should be set relative to this, say in the range 
               0.001 - 0.01  
```
Note that it is assumed that the matching of x1 and x2 are putative and it
is expected that a percentage of matches will be wrong.
```
Returns:
         F       - The 3x3 fundamental matrix such that x2'Fx1 = 0.
         inliers - An array of indices of the elements of x1, x2 that were
                   the inliers for the best model.
```
See also: ransac(), fundmatrix()
"""

# February 2004  Original version
# August   2005  Distance error function changed to match changes in RANSAC

function ransacfitfundmatrix(x1i, x2i, t::Real, feedback::Bool = false)

    # Define some subfunctions

    #--------------------------------------------------------------------------
    # Function to evaluate the first order approximation of the geometric error
    # (Sampson distance) of the fit of a fundamental matrix with respect to a
    # set of matched points as needed by RANSAC.  See: Hartley and Zisserman,
    # 'Multiple View Geometry in Computer Vision', page 270.
    #
    # Note that this code allows for F being an array of fundamental matrices of
    # which we have to pick the best one. (A 7 point solution can return up to 3
    # solutions)
    
    function funddist(F, x, t)
        
        npts = size(x, 2)
        x1 = x[1:3,:]    # Extract x1 and x2 from x
        x2 = x[4:6,:]
        
        if isa(F, Array{Any})  # We have several solutions each of which must be tested
	    nF = length(F)   # Number of solutions to test
	    bestF = F[1]     # Initial allocation of best solution
	    ninliers = 0     # Number of inliers
	
	    for k = 1:nF
	        x2tFx1 = zeros(npts)
	        for n = 1:npts
		    x2tFx1[n] = (x2[:,n:n]'*F[k]*x1[:,n:n])[1]
	        end
	        
	        Fx1 = F[k]*x1
	        Ftx2 = F[k]'*x2     
                
	        # Evaluate distances
	        d =  x2tFx1.^2 ./ 
		vec(Fx1[1,:].^2 + Fx1[2,:].^2 + Ftx2[1,:].^2 + Ftx2[2,:].^2)
	    
	        inliers = find(abs(d) .< t)     # Indices of inlying points
	        
	        if length(inliers) > ninliers   # Record best solution
		    ninliers = length(inliers)
		    bestF = copy(F[k])
		    bestInliers = copy(inliers)
	        end
	    end
            
        else     # We just have one solution
	    x2tFx1 = zeros(npts)
	    for n = 1:npts
	        x2tFx1[n] = (x2[:,n:n]'*F*x1[:,n:n])[1] # Ugly thanks to v0.5
	    end
	    
	    Fx1 = F*x1
	    Ftx2 = F'*x2     

	    # Evaluate distances (The call to vec() needed for 0.4)
	    d =  x2tFx1.^2 ./ 
	    vec(Fx1[1,:].^2 + Fx1[2,:].^2 + Ftx2[1,:].^2 + Ftx2[2,:].^2)
	    
	    bestInliers = find(abs(d) .< t)     # Indices of inlying points
	    bestF = F                          # Copy F directly to bestF
        end

        return bestInliers, bestF
    end

    #----------------------------------------------------------------------
    # (Degenerate!) function to determine if a set of matched points will result
    # in a degeneracy in the calculation of a fundamental matrix as needed by
    # RANSAC.  This function assumes this cannot happen...
    
    function isdegenerate(x)
        return false
    end    
    
    #----------------------------------------------------------------------
    # The main function

    if size(x1i) != size(x2i)
        error("Data sets x1 and x2 must have the same dimension")
    end
    
    (rows, npts) = size(x1i)
    if !(rows==2 || rows==3)
        error("x1 and x2 must have 2 or 3 rows")
    end
    
    if rows == 2    # Pad data with homogeneous scale factor of 1
        x1 = makehomogeneous(x1i)
        x2 = makehomogeneous(x2i)
    else
        x1 = copy(x1i)
        x2 = copy(x2i)
    end
    
    # Normalise each set of points so that the origin is at centroid and
    # mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    # scale parameter is 1.  Note that fundmatrix() will also call
    # normalise2dpts() but the code in ransac() that calls the distance
    # function will not - so it is best that we normalise beforehand.
    (x1, T1) = normalise2dpts(x1)
    (x2, T2) = normalise2dpts(x2)

    s = 8  # Number of points needed to fit a fundamental matrix. Note that
           # only 7 are needed but the function fundmatrix() only
           # implements the 8-point solution.
    
    fittingfn = fundmatrix
    distfn    = funddist
    degenfn   = isdegenerate
    # x1 and x2 are 'stacked' to create a 6xN array for ransac
    (F, inliers) = ransac([x1; x2], fittingfn, distfn, degenfn, s, t, feedback)

    # Now do a final least squares fit on the data points considered to
    # be inliers.
    F = fundmatrix(x1[:,inliers], x2[:,inliers])
    
    # Denormalise
    F = T2'*F*T1

    return F, inliers

end # of ransacfitfundmatrix()


#-----------------------------------------------------------------------
"""
ransacfitaffinefundmatrix - Fits affine fundamental matrix using RANSAC
```
Usage:   (F, inliers) = ransacfitaffinefundmatrix(x1, x2, t)

Arguments:
         x1  - 2xN or 3xN set of homogeneous points.  If the data is
               2xN it is assumed the homogeneous scale factor is 1.
         x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
         t   - The distance threshold between data point and the model
               used to decide whether a point is an inlier or not. 
               Note that point coordinates are normalised to that their
               mean distance from the origin is sqrt(2).  The value of
               t should be set relative to this, say in the range 
               0.001 - 0.01  
```
Note that it is assumed that the matching of x1 and x2 are putative and it
is expected that a percentage of matches will be wrong.
```
Returns:
         F       - The 3x3 fundamental matrix such that x2'Fx1 = 0.
         inliers - An array of indices of the elements of x1, x2 that were
                   the inliers for the best model.
```
See also: ransac(), fundmatrix(), affinefundmatrix()
"""

# February 2004  Original version
# August   2005  Distance error function changed to match changes in RANSAC

function ransacfitaffinefundmatrix(x1i, x2i, t::Real, feedback::Bool = false)

    # Define some subfunctions
    #--------------------------------------------------------------------------
    # Function to evaluate the first order approximation of the geometric error
    # (Sampson distance) of the fit of a fundamental matrix with respect to a
    # set of matched points as needed by RANSAC.  See: Hartley and Zisserman,
    # 'Multiple View Geometry in Computer Vision', page 270.
    
    function funddist(F, x, t)
    
        ( _ , npts) = size(x)
        x1 = x[1:3,:]    # Extract x1 and x2 from x
        x2 = x[4:6,:]
        
        x2tFx1 = zeros(npts)
        for n = 1:npts
	    x2tFx1[n] = (x2[:,n]'*F*x1[:,n])[1]
        end
        
        Fx1 = F*x1
        Ftx2 = F'*x2     
        
        # Evaluate distances
        d =  x2tFx1.^2 ./ 
	vec(Fx1[1,:].^2 + Fx1[2,:].^2 + Ftx2[1,:].^2 + Ftx2[2,:].^2)
        
        inliers = find(abs(d) .< t)     # Indices of inlying points
        return inliers, F
    end

    #----------------------------------------------------------------------
    # (Degenerate!) function to determine if a set of matched points will result
    # in a degeneracy in the calculation of a fundamental matrix as needed by
    # RANSAC.  This function assumes this cannot happen...
    
    function isdegenerate(x)
        return false
    end    
    #----------------------------------------------------------------------
    # The main function

    if size(x1i) != size(x2i)
        error("Data sets x1 and x2 must have the same dimension")
    end
    
    (rows, npts) = size(x1i)
    if !(rows ==2 || rows ==3)
        error("x1 and x2 must have 2 or 3 rows")
    end
    
    if rows == 2    # Pad data with homogeneous scale factor of 1
        x1 = makehomogeneous(x1i)
        x2 = makehomogeneous(x2i)
    else
        x1 = copy(x1i)
        x2 = copy(x2i)
    end
    
    # Normalise each set of points so that the origin is at centroid and
    # mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    # scale parameter is 1.  Note that 'fundmatrix' will also call
    # 'normalise2dpts' but the code in 'ransac' that calls the distance
    # function will not - so it is best that we normalise beforehand.
    (x1, T1) = normalise2dpts(x1)
    (x2, T2) = normalise2dpts(x2)

    s = 4  # Number of points needed to fit an affine fundamental matrix. 
    
    fittingfn = affinefundmatrix
    distfn    = funddist
    degenfn   = isdegenerate
    # x1 and x2 are 'stacked' to create a 6xN array for ransac
    (F, inliers) = ransac([x1; x2], fittingfn, distfn, degenfn, s, t, feedback)

    # Now do a final least squares fit on the data points considered to
    # be inliers.
    F = affinefundmatrix(x1[:,inliers], x2[:,inliers])
    
    # Denormalise
    F = T2'*F*T1

    return F, inliers

end # of ransacfitaffinefundmatrix


#-----------------------------------------------------------------------
"""
ransacfitplane - Fits plane to 3D array of points using RANSAC
```
Usage  (B, P, inliers) = ransacfitplane(XYZ, t, feedback)

This function uses the RANSAC algorithm to robustly fit a plane
to a set of 3D data points.

Arguments:
         XYZ - 3xNpts array of xyz coordinates to fit plane to.
         t   - The distance threshold between data point and the plane
               used to decide whether a point is an inlier or not.
         feedback - Optional flag 0 or 1 to turn on RANSAC feedback
                    information.

Returns:
          B - 4x1 array of plane coefficients in the form
              b[1]*X + b[2]*Y +b[3]*Z + b[4] = 0
              The magnitude of B is 1.
              This plane is obtained by a least squares fit to all the
              points that were considered to be inliers, hence this
              plane will be slightly different to that defined by P below.
          P - The three points in the data set that were found to
              define a plane having the most number of inliers.
              The three columns of P defining the three points
    inliers - The indices of the points that were considered
              inliers to the fitted plane.
```
See also:  ransac(), fitplane()
"""
#=
June 2003 - Original version.
Feb  2004 - Modified to use separate ransac function
Aug  2005 - planeptdist modified to fit new ransac specification
Dec  2008 - Much faster distance calculation in planeptdist (thanks to
            Alastair Harrison) 
=#

function ransacfitplane(XYZ, t::Real, feedback::Bool = false)

    # Define some sub functions
    #------------------------------------------------------------------------
    # Function to define a plane given 3 data points as required by
    # RANSAC. In our case we use the 3 points directly to define the plane.
    function defineplane(X)
        return X
    end 
   
    #------------------------------------------------------------------------
    # Function to calculate distances between a plane and a an array of points.
    # The plane is defined by a 3x3 matrix, P.  The three columns of P defining
    # three points that are within the plane.
    function planeptdist(P, X, t)
        
        n = cross(P[:,2]-P[:,1], P[:,3]-P[:,1]) # Plane normal.
        n = n/norm(n)                           # Make it a unit vector.
        
        npts = size(X,2)
        d = zeros(npts)   # d will be an array of distance values.
        
        # The following loop builds up the dot product between a vector from P[:,1]
        # to every X[:,i] with the unit plane normal.  This will be the
        # perpendicular distance from the plane for each point
        for i=1:3
	    d = d + vec(X[i,:]-P[i,1])*n[i] 
        end
        
        inliers = find(abs(d) .< t)
        
        return inliers, P
    end    
    #------------------------------------------------------------------------
    # Function to determine whether a set of 3 points are in a degenerate
    # configuration for fitting a plane as required by RANSAC.  In this case
    # they are degenerate if they are colinear.
    
    function isdegenerate(X)
        
        # The three columns of X are the coords of the 3 points.
        r = iscolinear(X[:,1],X[:,2],X[:,3])
    end
    
    #------------------------------------------------------------------------
    # The main function
    
    (rows, npts) = size(XYZ)
    
    if rows !=3
        error("Data is not 3D")
    end
    
    if npts < 3
        error("Too few points to fit plane")
    end
    
    s = 3  # Minimum No of points needed to fit a plane.
        
    fittingfn = defineplane
    distfn    = planeptdist
    degenfn   = isdegenerate

    (P, inliers) = ransac(XYZ, fittingfn, distfn, degenfn, s, t, feedback)
    
    # Perform least squares fit to the inlying points
    B = fitplane(XYZ[:,inliers])

    return B, P, inliers
end # of ransacfitplane


#-----------------------------------------------------------------------
"""
ransacfitline - Fits line to 3D array of points using RANSAC
```
Usage  (V, L, inliers) = ransacfitline(XYZ, t, feedback)

This function uses the RANSAC algorithm to robustly fit a line to a
set of 3D data points or to a set of 2D homogeneous points with scale
value of 1

Arguments:
         XYZ - 3xNpts array of xyz coordinates to fit line to.
         t   - The distance threshold between data point and the line
               used to decide whether a point is an inlier or not.
    feedback - Optional boolean flag to turn on RANSAC feedback
               information.

Returns:.
          V - Line obtained by a simple fitting on the points that
              are considered inliers.  The line goes through the
              calculated mean of the inlier points, and is parallel to
              the principal eigenvector.  The line is scaled by the
              square root of the largest eigenvalue.
              This line is returned as a nx2 matrix.  The first column 
              is the beginning point, the second column is the end point 
              of the line.
          L - The two points in the data set that were found to
              define a line having the most number of inliers.
              The two columns of L defining the two points.
    inliers - The indices of the points that were considered
              inliers to the fitted line.
```
See also:  ransac(), fitplane(), ransacfitplane()
"""
# Aug  2006 - created ransacfitline from ransacfitplane
#             author: Felix Duvallet

function ransacfitline(XYZ, t::Real, feedback::Bool = false)

    # First define some internal functions
    # ------------------------------------------------------------------------
    # Function to define a line given 2 data points as required by
    # RANSAC.  X is a 2 column matrix where each column defines a
    # point. Here we use the 2 points directly to define the line
    
    function defineline(X)
        return X
    end    

    #------------------------------------------------------------------------
    # Function to calculate distances between a line and an array of points.
    # The line is defined by a 3x2 matrix, L.  The two columns of L defining
    # two points that are the endpoints of the line.
    #
    # A line can be defined with two points as:
    #        lambda*p1 + (1-lambda)*p2
    # Then, the distance between the line and another point (p3) is:
    #        norm( lambda*p1 + (1-lambda)*p2 - p3 )
    # where
    #                  (p2-p1).(p2-p3)
    #        lambda =  ---------------
    #                  (p1-p2).(p1-p2)
    #
    # lambda can be found by taking the derivative of:
    #      (lambda*p1 + (1-lambda)*p2 - p3)*(lambda*p1 + (1-lambda)*p2 - p3)
    # with respect to lambda and setting it equal to zero
    
    function lineptdist(L, X, t)
        
        p1 = L[:,1]
        p2 = L[:,2]
        
        npts = size(X,2)
        d = zeros(npts)
        
        for i = 1:npts
            p3 = X[:,i]
            
            lambda = dot((p2 - p1), (p2-p3)) / dot( (p1-p2), (p1-p2) )
            
            d[i] = norm(lambda*p1 + (1-lambda)*p2 - p3)
        end
        
        inliers = find(abs(d) .< t)
        return inliers, L
    end    

    #------------------------------------------------------------------------
    # Function to determine whether a set of 2 points are in a degenerate
    # configuration for fitting a line as required by RANSAC.
    # In this case two points are degenerate if they are the same point
    # or if they are exceedingly close together.
    
    function isdegenerate(X)
        # Find the norm of the difference of the two points. 
        # If this is zero we are degenerate
        return norm(X[:,1] - X[:,2]) < eps()
    end

    #------------------------------------------------------------------------
    # Now the main function ...
    
    (rows, npts) = size(XYZ)
    
    if rows !=3
        error("Data is not 3D")
    end
    
    if npts < 2
        error("Too few points to fit line")
    end
    
    s = 2  # Minimum No of points needed to fit a line.
        
    fittingfn = defineline
    distfn    = lineptdist
    degenfn   = isdegenerate

    (L, inliers) = ransac(XYZ, fittingfn, distfn, degenfn, s, t, feedback)
    
    # Find the line going through the mean, parallel to the major
    # eigenvector
    V = fitline3d(XYZ[:, inliers])

    return V, L, inliers

    
end # of ransacfitline

#-----------------------------------------------------------------------
"""
iscolinear - Are 3 points colinear
```
Usage:  r = iscolinear(p1, p2, p3, flag)

Arguments:
       p1, p2, p3 - Points in 2D or 3D.
       homog      - An optional boolean flag indicating that p1, p2, p3 
                    are homogneeous coordinates with arbitrary scale.  
                    The default is false, that is it is assumed that the
                    points are inhomogeneous, or that they are homogeneous 
                    ewith qual scale.

Returns:
       r  -  true if points are co-linear, false otherwise
```
"""
# February 2004
# January  2005 - modified to allow for homogeneous points of arbitrary
#                 scale (thanks to Michael Kirchhof)

function iscolinear(p1i::Vector, p2i::Vector, p3i::Vector, homog::Bool = false)

    if (length(p1i) != length(p2i))  || (length(p1i) != length(p3i)) || !(length(p1i) == 2 || length(p1i) == 3)
        error("Points must have the same dimension of 2 or 3")
    end
    
    # If data is 2D, assume they are 2D inhomogeneous coords. Make them
    # homogeneous with scale 1.
    if length(p1i) == 2    
        p1 = [p1i; 1]; p2 =[p2i; 1]; p3 = [p3i; 1]
    else
        p1 = p1i; p2 = p2i; p3 = p3i
    end

    if homog
	# Apply test that allows for homogeneous coords with arbitrary
        # scale.  p1 X p2 generates a normal vector to plane defined by
        # origin, p1 and p2.  If the dot product of this normal with p3
        # is zero then p3 also lies in the plane, hence co-linear.
	r = abs(dot(cross(p1, p2),p3)) < eps()
    else
	# Assume inhomogeneous coords, or homogeneous coords with equal
        # scale.
	r = norm(cross(p2-p1, p3-p1)) < eps()
    end

    return r    
end

#-----------------------------------------------------------------------
"""
fitline2d - Least squares fit of a line to a set of 2D points
```
Usage:   (C, dist) = fitline2d(XY)

Where:   XY  - 2xNpts array of xy coordinates to fit line to data of
               the form 
               [x1 x2 x3 ... xN
                y1 y2 y3 ... yN]
                
               XY can also be a 3xNpts array of homogeneous coordinates.

Returns: C    - 3x1 array of line coefficients in the form
                c[1]*X + c[2]*Y + c[3] = 0
         dist - Array of distances from the fitted line to the supplied
                data points.  Note that dist is only calculated if the
                function is called with two output arguments.
```
The magnitude of C is scaled so that line equation corresponds to
```
   sin(theta)*X + (-cos(theta))*Y + rho = 0
```
where theta is the angle between the line and the x axis and rho is the
perpendicular distance from the origin to the line.  Rescaling the
coefficients in this manner allows the perpendicular distance from any
point (x,y) to the line to be simply calculated as
```
   r = abs(c[1]*X + c[2]*Y + c[3])
```
If you want to convert this line representation to the classical form 
```
     Y = a*X + b
use
  a = -c[1]/c[2]
  b = -c[3]/c[2]
```
Note, however, that this assumes c[2] is not zero
"""
# June      2003 - Original version
# September 2004 - Rescaling to allow simple distance calculation.
# November  2008 - Normalising of coordinates added to condition the solution.

function fitline2d(XYi)
    
    (rows,npts) = size(XYi)
    
    if npts < 2
        error("Too few points to fit line")
    end  
    
    if rows == 2    # Add homogeneous scale coordinate of 1 
        XY = makehomogneous(XYi)
    else
        XY = copy(XYi)
    end
    
    if npts == 2    # Pad XY with a third column of zeros
        XY = [XY zeros(3,1)] 
    end
    
    # Normalise points so that centroid is at origin and mean distance from
    # origin is sqrt(2).  This conditions the equations
    (XYn, T) = normalise2dpts(XY)
    
    # Set up constraint equations of the form  XYn'*C = 0,
    # where C is a column vector of the line coefficients
    # in the form   c(1)*X + c(2)*Y + c(3) = 0.
    
    (u, d, v) = svd(XYn')   # Singular value decomposition.
    C = v[:,3]              # Solution is last column of v.
    
    # Denormalise the solution
    C = T'*C
    
    # Rescale coefficients so that line equation corresponds to
    #   sin(theta)*X + (-cos(theta))*Y + rho = 0
    # so that the perpendicular distance from any point [x,y] to the line
    # to be simply calculated as 
    #   r = abs(c[1]*X + c[2]*Y + c[3])
    theta = atan2(C[1], -C[2])
    
    # Find the scaling (but avoid dividing by zero)
    if abs(sin(theta)) > abs(cos(theta))
        k = C[1]/sin(theta)
    else
        k = -C[2]/cos(theta)
    end
    
    C = C/k
    
    # Calculate the distances from the fitted line to the supplied
    # data points
    dist = abs(C[1]*XY[1,:] + C[2]*XY[2,:] + C[3])
    
    return C, dist
end

#-----------------------------------------------------------------------
"""
fitline3d - Fits a line to a set of 3D points
```
Usage:   L = fitline3d(XYZ)

Where: XYZ - 3xNpts array of XYZ coordinates
              [x1 x2 x3 ... xN;
               y1 y2 y3 ... yN;
               z1 z2 z3 ... zN]

Returns: L - 3x2 matrix consisting of the two endpoints of the line
             that fits the points.  The line is centered about the
             mean of the points, and extends in the directions of the
             principal eigenvectors, with scale determined by the
             eigenvalues.
```
"""
# Author: Felix Duvallet (CMU)
# August 2006

function fitline3d(XYZ)

    # Since the covariance matrix should be 3x3 (not NxN), need
    # to take the transpose of the points.
    
    # find mean of the points
    mu = mean(XYZ', 1)
    
    # covariance matrix
    C = cov(XYZ')
    
    # get the eigenvalues and eigenvectors
    (D, V) = eig(C)
    
    # largest eigenvector is in the last column
    col = size(V, 2)  # get the number of columns
    
    # get the last eigenvector column and the last eigenvalue
    eVec = V[:, col]
    eVal = D[col]

    L = zeros(3,2)    
    # start point - center about mean and scale eVector by eValue
    L[:, 1] = mu' - sqrt(eVal)*eVec
    L[:, 2] = mu' + sqrt(eVal)*eVec      # end point
    
    return L
end
#-----------------------------------------------------------------------
"""
fitplane - Solves coefficients of plane fitted to 3 or more points
```
Usage:   B = fitplane(XYZ)

Where:   XYZ - 3xNpts array of xyz coordinates to fit plane to.   
               If Npts is greater than 3 a least squares solution 
               is generated.

Returns: B   - 4x1 array of plane coefficients in the form
               B[1]*X + B[2]*Y + B[3]*Z + B[4] = 0
               The magnitude of B is 1.
```
See also: ransacfitplane()
"""
# June 2003

function fitplane{T<:Real}(XYZ::Array{T,2})
    
    (rows, npts) = size(XYZ)    
    
    if rows != 3
        error("Data is not 3D")
    end
    
    if npts < 3
        error("Too few points to fit plane")
    end
    
    # Set up constraint equations of the form  AB = 0,
    # where B is a column vector of the plane coefficients
    # in the form   B[1]*X + B[2]*Y + B[3]*Z + B[4] = 0
    A = [XYZ' ones(npts,1)]  # Build constraint matrix
    
    if npts == 3             # Pad A with zeros
        A = [A; zeros(1,4)] 
    end
    
    (u, d, v) = svd(A)       # Singular value decomposition.

    return B = v[:,4]        # Solution is last column of v.
end

