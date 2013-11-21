#! @System PrimaryDecompositionZeroDim
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular( ) * " x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "x^3-x, y*x^2-y,y^2-x^2,z-x*y", A );
#! <A torsion-free ideal given by 4 generators>
Decomp := PrimaryDecompositionZeroDim( I );
#! [ <A torsion-free ideal given by 3 generators>, 
#!   <A torsion-free ideal given by 3 generators>, 
#!   <A torsion-free ideal given by 3 generators>, 
#!   <A torsion-free ideal given by 3 generators>, 
#!   <A torsion-free ideal given by 3 generators> ]
Perform( Decomp, Display );
#! z-1,
#! y-1,
#! x-1 
#!
#! An ideal generated by the 3 entries of the above matrix
#! z+1,
#! y+1,
#! x-1 
#!
#! An ideal generated by the 3 entries of the above matrix
#! z,
#! y,
#! x 
#!
#! An ideal generated by the 3 entries of the above matrix
#! z+1,
#! y-1,
#! x+1 
#!
#! An ideal generated by the 3 entries of the above matrix
#! z-1,
#! y+1,
#! x+1 
#!
#! An ideal generated by the 3 entries of the above matrix
GF3 := HomalgRingOfIntegersInSingular( 3 );
#! GF(3)
K := GF3 * "x";
#! GF(3)[x]
J := LeftSubmodule( "( x^2 - 1 ) * ( x + 2 )", K);
#! <A principal torsion-free ideal given by a cyclic generator>
Decomp := PrimaryDecompositionZeroDim( J );
#! [ <A principal ideal given by a cyclic generator>, 
#!   <A principal ideal given by a cyclic generator> ]
Perform( Decomp, Display );
#! x^2+x+1
#!
#! An ideal generated by the entry of the above matrix
#! x+1
#!
#! An ideal generated by the entry of the above matrix
A := HomalgFieldOfRationalsInSingular( ) * " x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "(x^2 +1)*(x^2 -2), y,z",A );
#! <A torsion-free ideal given by 3 generators>
Decomp := PrimaryDecompositionZeroDim( I );
#! [ <A torsion-free ideal given by 3 generators>, 
#!   <A torsion-free ideal given by 3 generators> ]
Perform( Decomp, Display );
#! z,   
#! y,   
#! x^2-2
#!
#! An ideal generated by the 3 entries of the above matrix
#! z,   
#! y,   
#! x^2+1
#!
#! An ideal generated by the 3 entries of the above matrix
IsBound(I!.ListOfNonLinearFactors);
#! true
Perform( I!.ListOfNonLinearFactors, Display);
#! t^2-2
#!
#! An ideal generated by the entry of the above matrix
#! t^2+1
#!
#! An ideal generated by the entry of the above matrix
#! @EndExample
