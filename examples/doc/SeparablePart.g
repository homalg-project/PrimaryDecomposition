#! @System SeparablePart
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular() * "t";
#! Q[t]
r := "t^4-t^3-5*t^2-3*t" / A;
#! t^4-t^3-5*t^2-3*t
SquareFreeFactors( r );
#! [ t+1, t, t-3 ]
SeparablePart( r );
#! t^3-2*t^2-3*t
#! @EndExample
