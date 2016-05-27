
println("testing ransac")

# Just a simple line fitting test for the moment
# Adapted from fitlinedemo()
outliers = 0.2
sigma = 0.02
outsigma = 30*sigma  # outlying points have a distribution that is
                     # 30 times as spread as the inlying points
npts = 100           # Number of 3D data points	
vpts = round(Int, (1-outliers)*npts)  # No of valid points
opts = npts - vpts                    # No of outlying points

# Define a line:
#    Y = m*X
#    Z = 0
m = 6    
# Generate npts points in the line
X = rand(1,npts)
Y = m*X
Z = zeros(size(Y))

XYZ =  [X
    	Y
    	Z]

# Add uniform noise of +/-sigma
XYZ = XYZ + (2*rand(size(XYZ))-1)*sigma

# Generate opts random outliers
n = size(XYZ,2)
ind = randperm(n)[1:opts]  # get a random set of point indices of length opts 


# Perform RANSAC fitting of the line
t = 0.02  # inlier tolerance
(V, P, inliers) = ransacfitline(XYZ, t)

# get slope of line and test values with a generous tolerance given
# the stochasic nature of things
slope = (V[2,2]-V[2,1])/(V[1,2]-V[1,1])
@test abs(slope-m) < 0.2
@test maximum(V[3,:]) < 0.1  # Check Z is close to 0
