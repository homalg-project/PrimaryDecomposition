#! @System FGLMdata
#! @Example
LoadPackage( "PrimaryDecomposition" );
#! true
Q := HomalgFieldOfRationalsInSingular( );
#! Q
A := Q * "x,y,z";
#! Q[x,y,z]
I := LeftSubmodule( "y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z", A );
#! <A torsion-free ideal given by 6 generators>
R := A / I;
#! Q[x,y,z]/( y*z-x, x*z-y, y^2-z^2, x*y-z, x^2-z^2, z^3-z )
F:= FGLMdata(R);
#! [ <A 5 x 5 matrix over an external ring>,
#!   <A 5 x 5 matrix over an external ring>,
#!   <A 5 x 5 matrix over an external ring> ]
Display(F[1]);
#! 0,0,0,1,0,
#! 0,0,1,0,0,
#! 0,1,0,0,0,
#! 0,0,0,0,1,
#! 0,0,0,1,0 
Display(F[2]);
#! 0,0,1,0,0,
#! 0,0,0,1,0,
#! 0,0,0,0,1,
#! 0,1,0,0,0,
#! 0,0,1,0,0 
Display(F[3]);
#! 0,1,0,0,0,
#! 0,0,0,0,1,
#! 0,0,0,1,0,
#! 0,0,1,0,0,
#! 0,1,0,0,0 
#! @EndExample
