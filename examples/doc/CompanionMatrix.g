#! @System CompanionMatrix
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular() * "x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z", A );
#! <A torsion-free ideal given by 6 generators>
R := A / I;
#! Q[x,y,z]/( y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z )
r := "x^2+x*y+z" / R;
#! |[ x^2+x*y+z ]|
mu := MinimalPolynomial( r );
#! t^3-2*t^2-3*t
M := CompanionMatrix( mu );
#! <An unevaluated non-zero 3 x 3 matrix over an external ring>
Display( last );
#! 2,3,0,
#! 1,0,0,
#! 0,1,0 
#! @EndExample
