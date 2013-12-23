#! @System IsPrimaryZeroDim
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular( ) * "x";
#! Q[x]
I := LeftSubmodule( "x^2 - 2*x + 1", A);
#! <A principal torsion-free ideal given by a cyclic generator>
IsPrimaryZeroDim( I );
#! true
IsPrimeZeroDim( I );
#! false
A := HomalgFieldOfRationalsInSingular( ) * " x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "x^3-x, y*x^2-y,y^2-x^2,z-x*y", A );
#! <A torsion-free ideal given by 4 generators>
IsPrimaryZeroDim( I );
#! false
IsBound( I!.AZeroDivisor );
#! true
I!.AZeroDivisor;
#! |[ x ]|
A := HomalgRingOfIntegersInSingular( 3, "t,s" ) * "x,y";
#! GF(3)(t,s)[x,y]
I := LeftSubmodule( "x- s*t, y-s", A );
#! <A torsion-free ideal given by 2 generators>
IsPrimaryZeroDim( I^3 );
#! true
#! @EndExample
