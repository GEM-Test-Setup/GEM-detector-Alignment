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

# makeTrack.c
Generates tracks using cos^2 sin distribution and stores them as histograms. 

# convertRaw.c
Converts histograms (x: adc channel y: signal) into tracks. option skipOffset skips uncertainty 
  * Initial conditions must be well fit to the data to get physically possible tracks in a reasonable amount of time 
    * signalHalfWidth (bins) is used to debounce peaks and determine integration range for the centroid
    * nGems number of gems)
    * nGemReadouts number of gem readout planes (likely nGems x 2, 1 readout for x and y)
    * signalMinIntegral used to determine multiplicity, assuming the distribtuion allows
    * sigLevel the maximum allowed difference in charge integral for two correlated planes 
    * Tracks are resolved by generating all charge-correlated permutations then resolving geometricaly (after applying offsets)
    * This has a poor time complexity, but due to the charge cut, n remains small.
    * In a denser tracking environment with a narrower charge range a less strict and more approximative method may be preferred
    * However, this code can succeed in a denser tracking environment given a strict and accurate charge correlation cut (see sigLevel)

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

# Z_Rotation_Regression.C
Calculates Z rotation offset and associated translation constants

# trackgenerator.C
Creates ROOT file to be read by Z_Rotation_Regression.C
Uniform muon angular distribution
