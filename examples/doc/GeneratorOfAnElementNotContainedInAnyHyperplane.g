#! @System GeneratorOfAnElementNotContainedInAnyHyperplane
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
Q := HomalgFieldOfRationalsInSingular( );
#! Q
e := CertainRows( HomalgIdentityMatrix( 5, Q), [ 1 .. 3 ] );
#! <An unevaluated diagonal right invertible sub-identity 3 x 5 matrix over an external ring>
lambda := GeneratorOfAnElementNotContainedInAnyHyperplane( e );
#! <A non-zero left regular 1 x 5 matrix over an external ring>
l := UnionOfRows( e, lambda );
#! <An unevaluated non-zero 4 x 5 matrix over an external ring>
lambda := GeneratorOfAnElementNotContainedInAnyHyperplane( l );
#! <A non-zero left regular 1 x 5 matrix over an external ring>
A := HomalgRingOfIntegersInSingular( 3, "t" );
#! GF(3)(t)
S := HomalgRingOfIntegersInSingular( 3, "t" );
#! GF(3)(t)
S!.RootOfBaseField:= 1;                       
#! 1
L := HomalgMatrix( "[ 2, t, 0, t, 0, 1 ]", 2, 3, A );
#! <A 2 x 3 matrix over an external ring>
GeneratorOfAnElementNotContainedInAnyHyperplane( L, S );
#! <A 1 x 3 matrix over an external ring>
#! @EndExample
