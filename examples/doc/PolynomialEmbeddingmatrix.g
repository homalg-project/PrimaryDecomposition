#! @System PolynomialEmbeddingmatrix
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular() * "x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "x^3-x, y*x^2-y,y^2-x^2,z-x*y", A );
#! <A torsion-free ideal given by 4 generators>
R := A / I;
#! Q[x,y,z]/( x^3-x, x^2*y-y, -x^2+y^2, -x*y+z )
r := "x^2+x*y+z" / R;
#! |[ x^2+x*y+z ]|
mu := MinimalPolynomial( r );
#! t^3-2*t^2-3*t
M := PolynomialEmbeddingmatrix( mu );
#! <An unevaluated non-zero 3 x 3 matrix over an external ring>
Display( last );
#! 2,3,0,
#! 1,0,0,
#! 0,1,0 
#! @EndExample
