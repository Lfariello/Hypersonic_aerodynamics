# Hypersonic_aerodynamics
Comparison between the flow field around a cone (3D) and a wedge (2D). Analysis of the Taylor-Maccoll flow field within the supersonic boundary layer of a cone at zero angle of attack in a hypersonic flow.  
Results are compared with those obtained for a wedge under the same asymptotic conditions.  
The associated function (AngoloUrto.m) is able to calculate the shock angle using the Newton-Raphson method (iterative method).  
The input and output angles are expressed in degrees.  
The other associated function (TaylorMaccoll.m) is able to solve the hypersonic conical flow field by integrating  the Taylor-Maccoll equations using the direct method. The function takes as input the asymptotic conditions and the flow deflection angle (in radians)  
and returns the thermo-fluid dynamic conditions within the shock layer.
