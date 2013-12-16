#! @System RepresentationOverCoefficientsRing
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
r := HomalgRingElement( "x^2 +xy +z", R );
#! |[ x^2+x*y+z ]|
RepresentationOverCoefficientsRing( r );
#! <A 5 x 5 matrix over an external ring>
Display(last);
#! 0,2,0,0,1,
#! 0,1,0,0,2,
#! 0,0,1,2,0,
#! 0,0,2,1,0,
#! 0,2,0,0,1 
#! @EndExample
