#! @System GeneratorOfAnElementNotContainedInAnyHyperplane
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
Q := HomalgFieldOfRationalsInSingular( );
#! Q
e := HomalgIdentityMatrix( 5, Q);
#! <An unevaluated 5 x 5 identity matrix over an external ring>
lambda := GeneratorOfAnElementNotContainedInAnyHyperplane( e );
#! <A non-zero left regular 1 x 5 matrix over an external ring>
l := UnionOfRows( e, lambda );
#! <An unevaluated non-zero left invertible 6 x 5 matrix over an external ring>
lambda := GeneratorOfAnElementNotContainedInAnyHyperplane( l );
#! <A non-zero left regular 1 x 5 matrix over an external ring>
#! @EndExample
