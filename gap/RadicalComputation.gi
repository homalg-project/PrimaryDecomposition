#############################################################################
##
##  RadicalComputation.gi                       PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Implementation of some functions for radical computation.
##
#############################################################################

####################################
#
# methods for attributes:
#
####################################

##
InstallMethod( PreparationForRadicalOfIdeal,
	"for an ideal",
	[ IsFinitelyPresentedSubmoduleRep and ConstructedAsAnIdeal ],
	
  function( I )
    local A, R, indets, Sep, deg, J, M, i, d;
    
    A := HomalgRing( I );
    
    ## step 1 and 2:
    R := A / I ;
    
    indets := Indeterminates( R );
    
    Sep := List( indets, a -> SeparablePart( MinimalPolynomial( a ) ) );
    
    if not IsPerfect( CoefficientsRing( A ) ) then
        
        Sep := PolysOverTheSameRing( Sep );
        
        deg := HomalgRing( Sep[1] )!.RootOfBaseField;
        
        A := CoefficientsRing( HomalgRing( Sep[1] ) ) * Indeterminates( R );
        
        I := LeftSubmodule( EntriesOfHomalgMatrix( A * MatrixOfSubobjectGenerators( I ) ), A );
        
        R := A / I;
        
    fi;
    
    indets := List( indets, a -> a / A );
    
    Sep := List( [ 1 .. Length( Sep ) ], i -> Value( Sep[i], indets[i] ) );
    
    J := HomalgMatrix( Sep, Length( Sep ), 1, A );
    
    ## step 3 and 4:
    J := IdealBasisToGroebner( IdealBasisOverCoefficientRing( R * J ) );
    
    if IsPerfect( CoefficientsRing( A ) ) then
        return LeftSubmodule( EntriesOfHomalgMatrix( J ), A );
    fi;
    
    ## step 5:
    return [ FGLMdata( HomalgRing( J ) / J ), I ];
    
end );

##
InstallMethod( RadicalOfIdeal,
	"for a an ideal",
	[  IsHomalgModule ],

  function( I )
    local A, M, R, K, p, i, f, d;  
    
    A := HomalgRing( I );
    
    if IsPerfect( CoefficientsRing( A ) ) then
        return PreparationForRadicalOfIdeal( I );
    fi;
    
    M := PreparationForRadicalOfIdeal( I )[1];
    
    I := PreparationForRadicalOfIdeal( I )[2];
    
    A := AmbientRing( I );
    
    R := A / I;
    
    K := CoefficientsRing( A ) * "var";
    
    p := Characteristic( CoefficientsRing( CoefficientsRing( A ) ) );
    
    for i in [ 1 .. Length( RationalParameters( A ) ) ] do
    
        f := ("( var^p )^deg" / K ) - RationalParameters[i] / K;
         
        M := List( M, i -> MatrixEmbedding( i ) );
    
    od;
        
    M := List( M, n -> n );
    
    d := NrRows( M[1] );
    
    return FGLMToGroebner( M, CertainRows( HomalgIdentityMatrix( d, CoefficientsRing( R ) ), [1] ) );
    
end );

##
InstallMethod( CompanionMatrix,
	"for a univariate polynomial",
        [ IsHomalgRingElement ],
        
  function( mu )
    local R, c, d, m, L, e;
    
    R := HomalgRing( mu );
    
    mu := CoefficientsOfUnivariatePolynomial( mu );
    
    m := NrColumns( mu );
    
    mu := CertainColumns( -mu, Reversed( [ 1 .. m - 1 ] ) );
    
    L := HomalgIdentityMatrix( m - 2, R );
    
    e := HomalgZeroMatrix( m - 2, 1, R );
    
    L := UnionOfColumns( L, e );
    
    return UnionOfRows( mu, L );
    
end );
