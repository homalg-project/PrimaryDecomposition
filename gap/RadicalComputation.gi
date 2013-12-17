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
    local A, R, indets, Sep, deg, list, J, M, i, d;
    
    A := HomalgRing( I );
    
    ## step 1 and 2:
    R := A / I ;
    
    indets := Indeterminates( R );
    
    if IsPerfect( CoefficientsRing( A ) ) then
        Sep := List( indets, a -> SeparablePart( MinimalPolynomial( a ) ) );
    else
        Sep := List( indets, a -> SepUnvollkommen( MinimalPolynomial( a ) ) );
    fi;
    
    if not IsPerfect( CoefficientsRing( A ) ) then
        
        Sep := PolysOverTheSameRing( Sep );
        
        deg := HomalgRing( Sep[1] )!.RootOfBaseField;
        
        A := CoefficientsRing( HomalgRing( Sep[1] ) ) * Indeterminates( A );
        
        list := EntriesOfHomalgMatrix( MatrixOfSubobjectGenerators( I ) );
        list := Add( list, Sep[1] );
        list := PolysOverTheSameRing( list );
        list := List( list , i -> i / A );
        
        I := LeftSubmodule( list, A );
        
        R := A / I;
        
    fi;
    
    indets := List( indets, a -> a / R );
    
    Sep := List( [ 1 .. Length( Sep ) ], i -> Value( Sep[i], indets[i] ) );
    
    J := HomalgMatrix( Sep, Length( Sep ), 1, A );
    
    ## step 3 and 4:
    J := IdealBasisToGroebner( IdealBasisOverCoefficientRing( R * J ) );
    J := LeftSubmodule( EntriesOfHomalgMatrix( J ), A );
    
    if IsPerfect( CoefficientsRing( A ) ) then
        return J;
    fi;
    
    ## step 5:
    return [ FGLMdata( HomalgRing( J ) / J ), I, deg ];
    
end );

##
InstallMethod( RadicalOfIdeal,
	"for a an ideal",
	[  IsHomalgModule ],

  function( I )
    local A, p, M, deg, R, K, i, f, d, e;  
    
    A := HomalgRing( I );
    
    if IsPerfect( CoefficientsRing( A ) ) then
        return PreparationForRadicalOfIdeal( I );
    fi;
    
    p := PreparationForRadicalOfIdeal( I );
    
    M := p[1];
    I := p[2];
    deg := p[3];
    
    A := HomalgRing( I );
    
    R := A / I;
    
    K := CoefficientsRing( A ) * "x";
    
    p := Characteristic( CoefficientsRing( CoefficientsRing( A ) ) );
    
    for i in [ 1 .. Length( RationalParameters( A ) ) ] do
    
        f := (("x"/K)^p)^deg - RationalParameters( A )[i] / K;
         
        M := List( M, i -> MatrixEmbedding( i , f ) );
    
    od;
        
    M := List( M, n -> n );
    
    d := NrRows( M[1] );
    
    e :=CertainRows( HomalgIdentityMatrix( d, CoefficientsRing( R ) ) , [1] );
    
    return FGLMToGroebner( M, e );
    
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
