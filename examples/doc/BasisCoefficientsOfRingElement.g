#! @System BasisCoefficientsOfRingElement
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular( ) * "x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z", A );
#! <A torsion-free ideal given by 6 generators>
R := A / I;
#! Q[x,y,z]/( y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z )
r := "x^4 * y^2 + z*x" / R;
#! |[ x^4*y^2+x*z ]|
coeffs := BasisCoefficientsOfRingElement( r );
#! <A 1 x 5 matrix over an external ring>
Display(coeffs);
#! 0,0,1,0,1
bas := BasisOverCoefficientsRing( R );
#! <A 5 x 1 matrix over a residue class ring>
Display( bas );
#! 1, 
#! z, 
#! y, 
#! x, 
#! z^2
#!
#! modulo [ y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z ]
R * coeffs * bas ;
#! <An unevaluated 1 x 1 matrix over a residue class ring>
Display(last);
#! z^2+y
#!
#! modulo [ y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z ]
#! @EndExample

