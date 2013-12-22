#! @System Derivative
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
R := HomalgFieldOfRationalsInSingular() * "t";
#! Q[t]
r := "t^4-t^3-5*t^2-3*t" / R;
#! t^4-t^3-5*t^2-3*t
Derivative( r );
#! 4*t^3-3*t^2-10*t-3
S := HomalgRingOfIntegersInSingular( 3 ) * "t";
#! GF(3)[t]
r := "t^5 + 3*t^4 - t^3 + 2*t" / S;  
#! t^5-t^3-t
Derivative( r );
#! -t^4-1
#! @EndExample
