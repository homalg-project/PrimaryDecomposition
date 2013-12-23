#! @System IsPrimeZeroDim
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular( ) * "x";
#! Q[x]
I1 := LeftSubmodule( "x^2 + 1", A);             
#! <A principal torsion-free ideal given by a cyclic generator>
IsPrimeZeroDim( I1 );
#! true
I2 := LeftSubmodule( "(x-1)^2", A);
#! <A principal torsion-free ideal given by a cyclic generator>
IsPrimeZeroDim( I2 );
#! false
I2!.AZeroDivisor;
#! |[ x ]|
I2!.ANilpotentElement;
#! |[ x ]|
#! @EndExample
