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
InstallMethod( RadicalOfIdeal,
	"for an ideal",
	[ IsFinitelyPresentedSubmoduleRep and ConstructedAsAnIdeal ],
	
  function( I )
    local A, R, indets, Sep, L, RI;
    
    A := HomalgRing( I );
    
    if not IsPerfect( CoefficientsRing( A ) ) then
        TryNextMethod( );
    fi;
    
    R := A / I ;
    
    indets := Indeterminates( R );
    
    Sep := List( indets, a -> SeparablePart( MinimalPolynomial( a ) ) );
    
    #for the non perfect fields: L := CoefficientsRing( HomalgRing( Sep[ 1 ] ) ) * Indeterminates( A );
    
    indets := List( indets, a -> a / A );
    
    Sep := List( [ 1 .. Length( Sep ) ], i -> Value( Sep[i], indets[i] ) );
    
    RI := HomalgMatrix( Sep, Length( Sep ), 1, A );
    
    ## missing two steps for not perfect fields -> matrixmatrixembedding, fglmToGroebner
    ## then L <> A
    
    return LeftSubmodule( EntriesOfHomalgMatrix( IdealBasisToGroebner( IdealBasisOverCoefficientRing( R * RI ) ) ), A );
    ##return LeftSubmodule( BasisOfRows( RI * A ) );
    
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
