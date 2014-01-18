#! @System CoefficientsTransformation
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
C := HomalgRingOfIntegersInSingular( 3, "t,s");
#! GF(3)(t,s)
deg := 2;
#! 2
M := HomalgMatrix( "[t,s*t, s^2, 1]", 4, 1, C );
#! <A 4 x 1 matrix over an external ring>
Display( M );
#! (t),  
#! (t*s),
#! (s^2),
#! 1     
CoefficientsTransformation( M, deg);
#! <A 4 x 1 matrix over an external ring>
Display( last );
#! (t^9),
#! (t^9*s^9),
#! (s^18),
#! 1         
#! @EndExample
