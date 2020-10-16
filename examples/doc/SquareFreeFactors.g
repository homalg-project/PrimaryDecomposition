#! @Chunk SquareFreeFactors
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular() * "t";
#! Q[t]
p := "t^3-2*t^2-3*t" / A;
#! t^3-2*t^2-3*t
s := SquareFreeFactors( p );
#! [ t, t+1, t-3 ]
Product( s );
#! t^3-2*t^2-3*t
#! @EndExample                        
