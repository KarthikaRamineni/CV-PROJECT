## Colour Homography Colour Correction

Colors across a change in viewing condition (changing light color, shading and camera) are related by a homography.
Colour correction : Given pixel - pixel correspondaces, colour homography is applied for shading invariant colour mapping.


### Methods
+ Least Squares
+ Alternate Least Squares
+ RANSAC

To test least squares and alternate least squares : `test.m`
To test RANSAC : `test_ransac.m`


### Evaluation
The above methods are compared using Mean, Median, 95th percentile, Max in the following colour spaces.
+ CIE Lab
+ CIE Luv
+ RGB

To evaluate the methods : `evaluate.m`


### Dataset
[HG Colour Checker](www2.cmp.uea.ac.uk/~ybb15eau/db/HG_ColourChecker.zip)


### Requirements
MATLAB


### Reference
Graham D. Finlayson, Han Gong, and Robert B. Fisher, 2016. Color Homography Color Correction


