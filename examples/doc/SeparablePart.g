#! @System SeparablePart
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular() * "t";
#! Q[t]
r := "t^4-t^3-5*t^2-3*t" / A;
#! t^4-t^3-5*t^2-3*t
SquareFreeFactors( r );
#! [ t, t+1, t-3 ]
SeparablePart( r );
#! t^3-2*t^2-3*t
A := HomalgRingOfIntegersInSingular( 3, "t,s")* "x";
#! GF(3)(t,s)[x]
r := "( x^3 - s ) * ( x^3 + t )" / A;
#! x^6+(t-s)*x^3+(-t*s)
SeparablePart( r );
#! x^2+(t-s)*x+(-t*s)
#! @EndExample
