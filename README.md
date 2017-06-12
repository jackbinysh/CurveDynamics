# CurveDynamics
Some C code to run curve dynamics

The update of the curve is done over 4 functions:

UpdateState(): this just calls UpdateGeometryFromPosition and UpdateVelocityFromGeometry

UpdateGeometryFromPosition(): Given 3D spatial points, compute curvature torsion, tangent normal binormal.

UpdateVelocityFromGeometry(): ALL THE MODEL DYNAMICS GOES HERE. compute curve velocity from the above geometric quanitites

UpdatePositionFromVelocity(): Update the 3D spatial points, given the velocities computed above.

To alter the curve dynamics being considered, only the function UpdateVelocityFromGeometry need be modified. Its currently set to the LIA.
