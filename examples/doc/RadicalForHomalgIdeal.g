#! @System RadicalOfIdeal
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
Q := HomalgFieldOfRationalsInSingular( );
#! Q
A := Q * "x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z", A );
#! <A torsion-free ideal given by 6 generators>
RI := RadicalOfIdeal( I );
#! <A torsion-free ideal given by 6 generators>
Display( last );
#! z^3-z,
#! y^3-y,
#! x^3-x 
#!
#! An ideal generated by the 3 entries of the above matrix
RI = I;
#! false
#! @EndExample
