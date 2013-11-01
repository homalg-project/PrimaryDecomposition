#! @System BasisOverCoefficientsRing

#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
Q := HomalgFieldOfRationalsInSingular( );
#! Q
A := Q * "x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "x^3-x, y*x^2-y,y^2-x^2,z-x*y", A );
#! <A torsion-free ideal given by 4 generators>
R := A / I;
#! Q[x,y,z]/( y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z )
bas := BasisOverCoefficientsRing( R );
#! <A 5 x 1 matrix over a residue class ring>
Display( bas );
#! 1,
#! x,
#! y,
#! z,
#! z^2
#! 
#! modulo [ y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z ]
#! @EndExample
