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
    
    Sep := List( indets, a -> SeparablePart( MinimalPolynomial( a ) ) );
    
    if not IsPerfect( CoefficientsRing( A ) ) then
        
        Sep := PolysOverTheSameRing( Sep );
                
        A := CoefficientsRing( HomalgRing( Sep[1] ) ) * Indeterminates( A );
        
        A!.RootOfBaseField := HomalgRing( Sep[1] )!.RootOfBaseField;
        
        list := EntriesOfHomalgMatrix( MatrixOfSubobjectGenerators( I ) );
        
        Add( list, Sep[1] );
        list := PolysOverTheSameRing( list );
        Unbind( list[ Length( list ) ] );
        
        list := List( list , i -> i / A );
        
        I := LeftSubmodule( list, A );
        
        R := A / I;
        
        R!.RootOfBaseField := A!.RootOfBaseField;
        
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
    return [ FGLMdata( HomalgRing( J ) / J ), I ];
    
end );

##
InstallMethod( RadicalOfIdeal,
	"for a an ideal",
	[  IsHomalgModule ],

  function( I )
    local A, p, M, J, deg, R, x, K, i, f, d, e, FGLM;  
    
    A := HomalgRing( I );
    
    if IsPerfect( CoefficientsRing( A ) ) then
        return PreparationForRadicalOfIdeal( I );
    fi;
    
    p := PreparationForRadicalOfIdeal( I );
    
    M := p[1];
    J := p[2];
    
    deg := HomalgRing( J )!.RootOfBaseField;
    
    A := HomalgRing( J );
    
    R := A / J;
    
    x := UnusedVariableName( CoefficientsRing( A ), "x" );
    
    K := CoefficientsRing( A ) * x;
    
    x := Indeterminates( K )[1];
    
    p := Characteristic( CoefficientsRing( CoefficientsRing( A ) ) );
    
    for i in [ 1 .. Length( RationalParameters( A ) ) ] do
    
        f := ( x^p)^deg - RationalParameters( A )[i] / K;
         
        M := List( M, i -> MatrixEmbedding( i , f ) );
    
    od;
        
    M := List( M, n -> n );
    
    d := NrRows( M[1] );
    
    e :=CertainRows( HomalgIdentityMatrix( d, CoefficientsRing( R ) ) , [1] );
    
    FGLM := FGLMToGroebner( M, e, List( Indeterminates( A ), Name ) )[2];
    
    return LeftSubmodule( FGLM, HomalgRing( I ) );
    
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
