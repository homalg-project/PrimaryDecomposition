LoadPackage( "PrimaryDecomposition" );

Q := HomalgFieldOfRationalsInSingular( );

A := Q * "x,y,z";

I := LeftSubmodule( "x^3-x, y*x^2-y,y^2-x^2,z-x*y", A );

R := A / I;
