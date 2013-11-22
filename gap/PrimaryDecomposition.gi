#############################################################################
##
##  PrimaryDecomposition.gi                     PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Implementation of some functions for primary decomposition.
##
#############################################################################

####################################
#
# methods for attributes:
#
####################################

##
InstallMethod( IsPrimeZeroDim,
	"for a zero dimensional ideal",
	[ IsFinitelyPresentedSubmoduleRep and ConstructedAsAnIdeal ],
	
  function( I )
    local A, C, R, indets, mu, sf_mu, degI, i, RadI, n, RmodRadI, degRadI, e, L, iter, 
          lambda, z, w, comb, bool, l, W, Wext;
    
    A := HomalgRing( I );
    
    C := CoefficientsRing( A );
    
    ## preliminary test for perfectness, should be replaced by IsPerfect
    if not IsPerfect( C ) then
        TryNextMethod( );
    fi;
    
    ## The first part of the algorithm computes the minimal polynomials of the 
    ## indeterminates of R and determines if at least one of them is irreducible 
    ## of degree dim_C( R/I ) or reducible. In the first case this element proves
    ## that the ideal I is a prime ideal, in the second case the ideal I cannot be
    ## a prime ideal.
    
    R := A / I;
    
    indets := Indeterminates( R );
    
    mu := List( indets, MinimalPolynomial );
    
    sf_mu := List( mu, SquareFreeFactors );
    
    degI := NrRows( BasisOverCoefficientsRing( R ) );
    
    for i in [ 1 .. Length( mu ) ] do
        if Length( sf_mu[i] ) > 1 then
            I!.AZeroDivisor := indets[i];
            return false;
        elif Degree( sf_mu[i][1] ) < Degree( mu[i] ) then
            I!.AZeroDivisor := indets[i];
            I!.ANilpotentElement := indets[i];
            return false;
        elif Degree( sf_mu[i][1] ) = degI then
            return true;
        fi;
    od;
    
    ## Now the algorithms asks if I is a radical ideal.
    ## If yes, it is a prime ideal. If not, then I cannot be a prime ideal.
    
    RadI := RadicalOfIdeal( I );
    
    if not IsSubset( I, RadI ) then
        return false;
    fi;
        
    ## The last part of the algorithm iterates over ring elements until
    ## the minimal polynomial of the element is either irreducible of degree 
    ## degRadI or reducible.
    ## It is splitted into two cases: Finite and infinite coefficient fields.
    
    ## Some values needed in both cases:
    
    RmodRadI := A / RadI;
    
    degRadI := NrRows( BasisOverCoefficientsRing( RmodRadI ) );
    
    n := Length( indets );
    
    e := HomalgIdentityMatrix( n, C );
    
    indets := HomalgMatrix( indets, Length( indets ), 1, R );
    
    ## First case: Coefficients ring is finite.
    ## The iteration goes over all ring elements.
    
    if IsFinite( C ) then
        
        L := [];
        L := List( [ 1 .. n ], i -> CertainRows( e, [i] ) );
        
        iter := Iterator( e );
        
        lambda := NextIterator( iter );  ## The zero element will be left out.
        
        while true do
            
            lambda := NextIterator( iter );

            if Position( L, lambda )= fail then
                
                w := ( R * lambda ) * indets;
                
                mu := MinimalPolynomial( w );
            
                if IsIrreducible( mu ) and Degree( mu ) = degRadI then
                
                    return true;
                
                elif not IsIrreducible( mu ) then
                
                    I!.AZeroDivisor := w;
                    return false;
                
                fi;
                
                Add( L , lambda);
                
            fi;
                        
        od;
        
    fi;
    
    ## Second case: Coefficients field is not finite.
    ## The iteration goes over elements, whose coefficients of the basis do 
    ## not lie in any hyperplane of the vectorspace over the coeffient field. 
    L := e;
     
    while true do
        
        lambda := GeneratorOfAnElementNotContainedInAnyHyperplane( L );
        
        w := ( R * lambda ) * indets;
            
        mu := MinimalPolynomial( w );
        
        if IsIrreducible( mu ) and Degree( mu ) = degRadI then
            
            return true;
            
        elif not IsIrreducible( mu ) then
        
            I!.AZeroDivisor := w;
            return false;
        
        fi;
        
        Add( L, lambda );
        
    od;
    
end );

##
InstallMethod ( IsPrimaryZeroDim,
	"for a zero dimensional ideal",
	[ IsFinitelyPresentedSubmoduleRep and ConstructedAsAnIdeal ],
        
  function( I )
    local bool, Rad;
    
    Rad := RadicalOfIdeal( I );
    
    bool := IsPrimeZeroDim( Rad );
        
    if IsBound( Rad!.AZeroDivisor ) then
        I!.AZeroDivisor := Rad!.AZeroDivisor;
    fi;
    
    if IsBound( Rad!.ANilpotentElement ) then
        I!.ANilpotentElement := Rad!.ANilpotentElement;
    fi;
    
    return bool;
    
end );

##
InstallMethod( PrimaryDecompositionZeroDim,
        [ IsHomalgObject ],

  function( I )
    local Decomp, A, R, a, fac, ListOfNonLinearFactors, factor, N, W, bas, i, J, M, j;
    
    Decomp := [ ]; 
    
    ## If I is already primary, then the algorithmus returns I itself.
    
    if IsPrimaryZeroDim( I ) then
        Decomp[1] := I;
        return Decomp;
    fi;
    
    ## If the ideal I is not primary, then the algorithmus should have found a 
    ## zerodivisor. If not, something went terribly wrong.
    
    if not IsBound( I!.AZeroDivisor ) then
        Error( "IsPrimaryZeroDim couldn't find a zero divisor" );
    fi;
    
    A := HomalgRing( I );
    
    R := A / I;
        
    a := I!.AZeroDivisor / R;
    
    ## Now the algorithm computes the pairwise coprime factors of the minimal
    ## polynomial of a.
    
    fac := PrimaryDecomposition( LeftSubmodule( MinimalPolynomial( a ) ) );
    
    ## If there are factors of the minimal polynomial of degree higher than one,
    ## they got recorded in a list, to simplify the complete primary decomposition
    ## afterwards.
    ListOfNonLinearFactors := [];
    
    for i in [ 1 .. Length( fac ) ] do
        factor := MatElm( MatrixOfSubobjectGenerators( fac[i][2] ), 1, 1 );
        if Degree( factor ) > 1 then
            
            ListOfNonLinearFactors[i] := fac[i][2];
            
            
        fi;
    od;
    
    if Length( ListOfNonLinearFactors ) > 0 then
        I!.ListOfNonLinearFactors:= ListOfNonLinearFactors;
    fi;
    
    fac := List( [ 1 .. Length( fac ) ], i -> MatrixOfSubobjectGenerators( fac[i][1] ) );
    fac := List( [ 1 .. Length( fac ) ], i -> MatElm( fac[i], 1 ,1 ) );
    
    ## Computation of the ring elements obtained by evaluating the factors of
    ## the minimal polynomial for the element a, their representation matrices 
    ## and their kernels.
    
    N := [ ];
    W := [ ];
        
    bas := BasisOverCoefficientsRing( R );
    
    a := a / R;
    
    for i in [ 1 .. Length( fac ) ] do
        N[i] := RepresentationOverCoefficientsRing( Value( fac[i] , a ) );
        W[i] := R * SyzygiesOfRows( N[i] ) * bas;
    od;
    
    ## Computing the Ideals M[i]. The Generators for M[ i ]: 
    ## all kernels W[j] without the ith and the generators of the ideal I.
    ## At least compute inductively the primary decompositions of the M[i] and
    ## return the union of them as the primary decomposition of the ideal I.
        
    J := Iterated( W, UnionOfRows );
    
    M := [ ];
    
    j := [ ];
    j[1] := 0;
    
    for i in [ 1 .. Length( fac ) ] do
        
        j[i + 1] := j[i] + NrRows( W[i] );
        
        M[i] := UnionOfRows( CertainRows( J , [ 1 .. j[i] ] ), CertainRows( J, [ j[i + 1] + 1 .. NrRows( J ) ] ) );
        
        M[i] := UnionOfRows( A * MatrixOfGenerators( I ), A * M[i] );
        
        M[i] := LeftSubmodule( BasisOfRows( M[i] ) );
        
        Append( Decomp, PrimaryDecompositionZeroDim( M[i] ) );
    od;
    
    return Decomp;

end );
