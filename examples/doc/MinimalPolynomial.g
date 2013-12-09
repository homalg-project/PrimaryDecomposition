#! @System MinimalPolynomial
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
r := "x^2+x*y+z" / R;
#! |[ x^2+x*y+z ]|
s := "x*y+z"/ R;
#! |[ x*y+z ]|
mu_r := MinimalPolynomial( r );
#! t^3-2*t^2-3*t
mu_s := MinimalPolynomial( s );
#! t^3-4*t
mu_r + mu_s;
#! 2*t^3-2*t^2-7*t
kt := HomalgRing( mu_r );
#! Q[t]
t := Indeterminates( kt )[1];
#! t
mu := MinimalPolynomial( r, t );
#! t^3-2*t^2-3*t
mu = mu_r;
#! true
mu := MinimalPolynomial( r, "t" );
#! t^3-2*t^2-3*t
mu = mu_r;
#! true
mu := MinimalPolynomial( r, "a" );
#! a^3-2*a^2-3*a
#! @EndExample
