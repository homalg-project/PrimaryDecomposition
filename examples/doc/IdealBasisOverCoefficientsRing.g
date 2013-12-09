#! @System IdealBasisOverCoefficientsRing
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
Q := HomalgFieldOfRationalsInSingular( );
#! Q
A := Q * "x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z", A );
#! <A torsion-free ideal given by 6 generators>
R := A / I;
#! Q[x,y,z]/( y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z )
G := HomalgMatrix( "[ z + 1, y - 1 ]", 1, 1, R );
#! <A 1 x 1 matrix over a residue class ring>
B := IdealBasisOverCoefficientRing( G );         
#! <An unevaluated 3 x 1 matrix over a residue class ring>
Display( last );
#! z+1, 
#! x+y, 
#! z^2+z
#!
#! modulo [ y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z ]
#! @EndExample
