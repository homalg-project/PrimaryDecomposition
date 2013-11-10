#############################################################################
##
##  RadicalComputation.gi                                    PrimaryDecomposition package
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
InstallMethod( RadicalForHomalgIdeal, 
	" for an homalg ideal",
	[ IsHomalgObject ],
	
  function( I )
     local A, ind, R, Sep, L, RI;
  
     A := HomalgRing( I );

     R := A / I ;
     
     ind := Indeterminates( R );
     
     Sep := List( ind, a -> SeparablePart( MinimalPolynomial( a ) ) );

     #L := CoefficientsRing( HomalgRing( Sep[ 1 ] ) ) * Indeterminates( A );
     
     ind := List( ind, a -> a / A );
     
     Sep := List( [1 .. Length( Sep )], i -> Value( Sep[i], ind[i] ));

     RI := HomalgMatrix( Sep, Length( Sep ), 1, A );
     
     ## missing two steps for not perfect fields -> matrixmatrixembedding, fglmToGroebner           
     ## then L <> A
     
     return LeftSubmodule( BasisOfRows( RI * A) );

end );

##
InstallMethod( PolynomialEmbeddingmatrix,
        [ IsHomalgRingElement ],
        
  function( mu )
    local R, c, d, m, L, e;
    
    R := HomalgRing( mu );
    
    mu := CoefficientsOfUnivariatePolynomial( mu );
    
    m := NrColumns( mu );
    
    mu := CertainColumns( -mu, Reversed( [ 1 .. m-1 ] ) );
    
    L := HomalgIdentityMatrix( m -2, R );
    
    e := HomalgZeroMatrix( m-2, 1, R );
    
    L := UnionOfColumns( L, e );
    
    return ( UnionOfRows( mu, L ) );
    
end );

