#! @System IsNotContainedInAnyHyperplane
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
Q := HomalgFieldOfRationalsInSingular( );
#! Q
L := HomalgIdentityMatrix( 5, Q );
#! <An unevaluated 5 x 5 identity matrix over an external ring>
Display( L );
#! 1,0,0,0,0,
#! 0,1,0,0,0,
#! 0,0,1,0,0,
#! 0,0,0,1,0,
#! 0,0,0,0,1 
lambda := HomalgMatrix( "[ 2,0,2,1,1]", 1, 5, Q );
#! <A 1 x 5 matrix over an external ring>
Display( lambda );
#! 2,0,2,1,1
IsNotContainedInAnyHyperplane( lambda, L );
#! false
lambda2 := GeneratorOfAnElementNotContainedInAnyHyperplane( L );;
IsNotContainedInAnyHyperplane( lambda2, L );
#! true
#! @EndExample                     
