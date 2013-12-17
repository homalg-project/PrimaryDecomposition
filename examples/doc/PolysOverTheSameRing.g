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
#! [ (t^6*s^6)*x^3+(s^12)*x+(t^6), (t^6)*x^2+(s^6) ]
R3 := HomalgRingOfIntegersInSingular( 3, "t,s" ) * "x";
#! GF(3)(t,s)[x]
R3!.RootOfBaseField := 0;
#! 0
r3 := "(t)*x^2+(-s)*x"/ R3;
#! (t)*x^2+(-s)*x
PolysOverTheSameRing( [ r1, r2, r3 ] );
#! [ (t^6*s^6)*x^3+(s^12)*x+(t^6), (t^6)*x^2+(s^6), (t^6)*x^2+(-s^6)*x ]
#! @EndExample
