#! @System PolysOverTheSameRing
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
R1 := HomalgRingOfIntegersInSingular( 3, "t,s" ) * "x";
#! GF(3)(t,s)[x]
R1!.RootOfBaseField := 3;
#! 3
R2 := HomalgRingOfIntegersInSingular( 3, "t,s" ) * "x";
#! GF(3)(t,s)[x]
R2!.RootOfBaseField := 2;
#! 2
r1 := "(t*s)*x^3+(s^2)*x+(t)" / R1;
#! (t*s)*x^3+(s^2)*x+(t)
r2 := "(t)*x^2+(s)" / R2;
#! (t)*x^2+(s)
PolysOverTheSameRing( [ r1, r2 ] );
#! [ (t^27*s^27)*x^3+(s^54)*x+(t^27), (t^81)*x^2+(s^81) ]
R3 := HomalgRingOfIntegersInSingular( 3, "t,s" ) * "x";
#! GF(3)(t,s)[x]
R3!.RootOfBaseField := 0;
#! 0
r3 := "(t)*x^2+(-s)*x"/ R3;
#! (t)*x^2+(-s)*x
PolysOverTheSameRing( [ r1, r2, r3 ] );
#! [ (t^27*s^27)*x^3+(s^54)*x+(t^27), (t^81)*x^2+(s^81), (t^729)*x^2+(-s^729)*x ]
#! @EndExample
