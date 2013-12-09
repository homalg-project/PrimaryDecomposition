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
e := One( R );
#! |[ 1 ]|
BasisCoefficientsOfRingElement( e );
#! <An unevaluated 1 x 5 matrix over an external ring>
Display(last);
#! 1,0,0,0,0
FGLMToGroebner( M, e );
#! [ [ x2*x3-x1, x2^2-x3^2, x1*x3-x2, x1*x2-x3, x1^2-x3^2, x3^3-x3 ], 
#!   [ 1, x3, x2, x1, x3^2 ] ]
bas := BasisOverCoefficientsRing( R );
#! <A 5 x 1 matrix over a residue class ring>
Display( bas );
#! 1,  
#! x1, 
#! x2, 
#! x3, 
#! x3^2
#!
#! modulo [ x2*x3-x1, x1*x3-x2, x2^2-x3^2, x1*x2-x3, x1^2-x3^2, x3^3-x3 ]
#! @EndExample
