# GEM-detector-Alignment

# linalg.h
Structs in place for:
  * Point
      * Access doubles with p.x, p.y, p.z, p.dx, p.dy, p.dz (for uncertainties)
  * Track 
      * Access Points with t[0], t[1], t[2]
  * Vector 
      * Access doubles with v[0], v[1], v[2]
  * Matrix
      * Access Column Vectors with m[0], m[1], m[2]
      * Operations assume 3x3 matrices

Useful selected methods:

  * getAngle 
    * returns the positive angle between two vectors
    * getAngle(Point, Point, Point)
    * getAngle(Vector, Vector)
  * getRotation 
    * returns a rotation matrix for the input radians about each axis
    * getRotation(double xAxis=0, double yAxis=0, double zAxis=0)
# convertRaw.c
Converts histograms (x: adc channel y: signal) into tracks. option skipOffset skips uncertainty 
  * For every 12 histograms attepmts to reconstruct atleast one track 
    * if there are mutliple signal peaks per histogram (assumed gaussian), it splits them
    * if there is an inconsistent amount of peaks or the peaks are too close (in any xy plane) to resolve ambiguity, skip the data
    * used a signal weighted integral to determine the centroid of the peak and uncertainty associated with this centroid
    * if too many or not enough track combinations are possible within uncertainty, skip the data.
    * writes these tracks to file and visualizes them using checkLine.c
# checkLine.c
Used for 3d track visualization
  * initVisualize
    * sets up the graph and draws the gem planes
    * initVisualize()
      * Uses default planes and range z = 0, 100, 200, range = (-5, 5)
    * initVisualize(** double)
      * takes array of 3 planes where [][0] is z, [][1] is length, and [][2] is width
    * initVisualize(Point, ** double)
      * takes origin point. all planes offset by this value
      * takes array as above
  * visualize
    * draws a track to the 3d graph with uncertainty bars
    * visualize(Point, Point, Point)
    * visualize(Track)
# regress.c
Used to determine track all offsets 
WIP

# Z_Rotation_Regression.C
Calculates Z rotation offset and associated translation constants

# trackgenerator.C
Creates ROOT file to be read by Z_Rotation_Regression.C
Uniform muon angular distribution
