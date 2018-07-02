# GEM-detector-Alignment


# linalg.h
##Structs in place for:
###Point    
...Access doubles with p.x, p.y, p.z, p.dx, p.dy, p.dz (for uncertainties)
###Track
...Access Points with t[0], t[1], t[2]
###Vector
...Access doubles with v[0], v[1], v[2]
###Matrix
...Access Column Vectors with m[0], m[1], m[2]
...Operations assume 3x3 matrices
##Useful nontrivial methods:

###getAngle returns the positive angle between two vectors
...getAngle(Point, Point, Point)
...getAngle(Vector, Vector)

##getRotation returns a rotation matrix for the input radians about each axis
...getRotation(double xAxis=0, double yAxis=0, double zAxis=0)
