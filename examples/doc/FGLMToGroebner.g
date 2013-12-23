#! @System FGLMToGroebner
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
A := HomalgFieldOfRationalsInSingular( ) * "x1,x2,x3";
#! Q[x1,x2,x3]
I := LeftSubmodule( "x2*x3-x1, x1*x3-x2, x2^2-x3^2, x1*x2-x3, x1^2-x3^2, x3^3-x3", A );
#! <A torsion-free ideal given by 6 generators>
R := A / I;
#! Q[x1,x2,x3]/( x2*x3-x1, x1*x3-x2, x2^2-x3^2, x1*x2-x3, x1^2-x3^2, x3^3-x3 )
M := FGLMdata( R );
#! [ <A 5 x 5 matrix over an external ring>, 
#!   <A 5 x 5 matrix over an external ring>, 
#!   <A 5 x 5 matrix over an external ring> ]
e := CertainRows( HomalgIdentityMatrix( NrRows( M[1] ), HomalgRing( M[1] ) ), [1] );
#! <An unevaluated diagonal right invertible sub-identity 1 x 5 matrix over an external ring>
FGLMToGroebner( M, e );
#! [ [ 1, x3, x2, x1, x3^2 ], 
#!   [ x2*x3-x1, x2^2-x3^2, x1*x3-x2, x1*x2-x3, x1^2-x3^2, x3^3-x3 ] ]
bas := BasisOverCoefficientsRing( R );
#! <A 5 x 1 matrix over a residue class ring>
Display( bas );
#! 1,  
#! x3, 
#! x2, 
#! x1, 
#! x3^2
#! 
#! modulo [ x2*x3-x1, x1*x3-x2, x2^2-x3^2, x1*x2-x3, x1^2-x3^2, x3^3-x3 ]
#! @EndExample
