LoadPackage( "PrimaryDecomposition" );

Q := HomalgFieldOfRationalsInSingular( );

A := Q * "x,y,z";

I := LeftSubmodule( "y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z", A );

R := A / I;

r := "x^2+x*y+z" / R;

p := MinimalPolynomial( r );
