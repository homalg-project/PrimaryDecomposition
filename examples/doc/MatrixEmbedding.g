#! @System MatrixEmbedding
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
R := HomalgRingOfIntegersInSingular( 3, "t,s" ) * "x";
#! GF(3)(t,s)[x]
R!.RootOfBaseField:=1;
#! 1
f := "x^3 - t" / R;
#! x^3+(-t)
p := "x^3 - s" / R;
#! x^3+(-s)
A := CoefficientsRing( R );
#! GF(3)(t,s)
M := HomalgMatrix( "[t*s - s, 1, 0, s]", 2, 2, A );
#! <A 2 x 2 matrix over an external ring>
Display( M );
#! (t*s-s),1, 
#! 0,      (s)
N := MatrixEmbedding( M, f );
#! <An unevaluated non-zero 6 x 6 matrix over an external ring>
Display( N );
#! (-s),0,   (t*s),1,  0,  0, 
#! (s), (-s),0,    0,  1,  0, 
#! 0,   (s), (-s), 0,  0,  1, 
#! 0,   0,   0,    (s),0,  0, 
#! 0,   0,   0,    0,  (s),0, 
#! 0,   0,   0,    0,  0,  (s)
S := HomalgRingOfIntegersInSingular( 3, "t" ) * "x";
#! GF(3)(t)[x]
S!.RootOfBaseRing := 3;
#! 3
I := LeftSubmodule( "x^2 - t^3", S );
#! <A principal torsion-free ideal given by a cyclic generator>
L := S / I;
#! GF(3)(t)[x]/( x^2+(-t^3) )
FGLM := FGLMdata( L );
#! [ <A 2 x 2 matrix over an external ring> ]
Display( FGLM[1] );
#! 0,    1,
#! (t^3),0
f := "x^3 - t" / S;
#! x^3+(-t)
M := MatrixEmbedding( FGLM[1] , f );
#! <An unevaluated non-zero 6 x 6 matrix over an external ring>
Display( M );
#! 0,  0,  0,  1,0,0,
#! 0,  0,  0,  0,1,0,
#! 0,  0,  0,  0,0,1,
#! (t),0,  0,  0,0,0,
#! 0,  (t),0,  0,0,0,
#! 0,  0,  (t),0,0,0 
T := HomalgRing( M );
#! GF(3)(t)
e := CertainRows( HomalgIdentityMatrix( NrRows( M ), HomalgRing( M ) ), [1] );
#! <An unevaluated diagonal right invertible sub-identity 1 x 6 matrix over an external ring>
FGLMToGroebner( [ M ], e );
#! [ [ 1, x ], [ x^2+(-t) ] ]
#! @EndExample
