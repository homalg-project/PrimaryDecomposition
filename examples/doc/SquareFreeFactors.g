#! @System SquareFreeFactors
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular() * "t";
#! Q[t]
p := "t^3-2*t^2-3*t" / A;
#! t^3-2*t^2-3*t
SquareFreeFactors( p );
#! [ t, t+1, t-3 ]
Product( last );
#! t^3-2*t^2-3*t
#! @EndExample                        
