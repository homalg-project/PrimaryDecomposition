#! @System RepresentationOverCoefficientsRing
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
#! Q[x,y,z]/( x^3-x, x^2*y-y, -x^2+y^2, -x*y+z )
r := HomalgRingElement( "x^2 +xy +z", R );
#! |[ x^2+x*y+z ]|
RepresentationOverCoefficientsRing( r );
#! <A 5 x 5 matrix over an external ring>
Display(last);
#! 0,0,0,2,1,
#! 0,1,2,0,0,
#! 0,2,1,0,0,
#! 0,0,0,1,2,
#! 0,0,0,2,1 
#! @EndExample
