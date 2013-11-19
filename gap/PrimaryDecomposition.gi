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
    local R, C, indets, mu, sf_mu, degI, i, RadI, n, RmodRadI, degRadI, e, L, iter, 
          lambda, z, w, comb, bool, l, W, Wext;
    
    R := HomalgRing( I );
    
    C := CoefficientsRing( R );
    
    ## preliminary test for perfectness, should be replaced by IsPerfect
    if not IsPerfect( C ) then
        TryNextMethod( );
    fi;
    
    ## The first part of the algorithm computes the minimal polynomials of the 
    ## indeterminates of R and determines if at least one of them is irreducible 
    ## of degree dim_C( R/I ) or reducible. In the first case this element proves
    ## that the ideal I is a prime ideal, in the second case the ideal I cannot be
    ## a prime ideal.
    
    indets := Indeterminates( R / I );
    
    mu := List( indets, MinimalPolynomial );
    
    sf_mu := List( mu, SquareFreeFactors );
    
    degI := NrRows( BasisOverCoefficientsRing( R / I ) );
    
    for i in [ 1 .. Length( mu ) ] do
        if Length( sf_mu[i] ) > 1 then
            I!.WitnessOfExistenceOfZeroDivisor := indets[i];
            return false;
        elif Degree( sf_mu[i][1] ) < Degree( mu[i] ) then
            I!.WitnessOfExistenceOfZeroDivisor := indets[i];
            I!.WitnessForExistenceOfNilpotentElement := indets[i];
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
        
    ## The last part of the algorithm computes witness elements. 
    ## Splitted in two cases: The coefficients ring is finite or not.
    ## Some values needed in both cases:
    
    n := Length( indets );
    
    RmodRadI := R / RadI;
    
    degRadI := NrRows( BasisOverCoefficientsRing( RmodRadI ) );
    
    e := HomalgIdentityMatrix( n, C );
        
    L := [];
        
    L := List( [ 1 .. n ], i -> CertainRows( e, [i] ) );
    
    indets := HomalgMatrix( indets, Length( indets ), 1, R/I );
    
    ## First case: Coefficients ring is finite.
    
    if IsFinite( C ) then
        
        ## After initializing an iterator the algorithm repeats the loop until
        ## a suitable element has fulfilled the asked properties.
        
        iter := Iterator( e );
        
        lambda := NextIterator( iter );  ## The zero element will be left out.
        
        while true do
            
            lambda := NextIterator( iter );

            if Position( L, lambda )= fail then
                
                w := ( ( R / I ) * lambda ) * indets;
                
                mu := MinimalPolynomial( w );
            
                if IsIrreducible( mu ) and Degree( mu ) = degRadI then
                
                    return true;
                
                elif not IsIrreducible( mu ) then
                
                    I!.WitnessOfExistenceOfZeroDivisor := w;
                    return false;
                
                fi;
                
                Add( L , lambda);
                
            fi;
                        
        od;
        
    fi;
    
    ## Second case: Coefficients ring is not finite.
        
    lambda := Iterated( L, \+ );
    
    l := Set( [ 1 .. Length( L ) ] );
        
    comb := Combinations( l, n - 1 );
    
    while true do
        
        ## the first part of this loop ensures that lambda is not contained in any
        ## n - 1 dimensional C-subspace spanned by elements of L. In this case
        ## bool remains to be 0.
                
        bool := 0; 
        
        for i in comb do
            
            z := List( i, L[ i ] );
            W := Iterated( z, UnionOfRows );
            Wext := UnionOfRows( W, lambda);
            
            if RowRankOfMatrix( W ) = n - 1 and RowRankOfMatrix( Wext ) < n then
                bool := 1;
            fi;
        od;
        
        ## after finding a suitable lambda the algorithm computes its minimal
        ## polynomial and checks like in the first part of the algorithm, if it is
        ## irreducible of degree dim C ( R / RadI) or reducible.
        
        if bool = 0 then        
        
            w := ( ( R / I ) * lambda ) * indets;
            
            mu := MinimalPolynomial( w );
        
            if IsIrreducible( mu ) and Degree( mu ) = degRadI then
                    
                return true;
                
            elif not IsIrreducible( mu ) then
                
                I!.WitnessOfExistenceOfZeroDivisor := w;
                return false;
                
            fi;
            
            Add( L, lambda );
        
        fi;
        
        ## the last part of the loop creats a new lambda. It is not sure so far,
        ## if lambda is suitable. That will be tested in the beginning of the loop.
        
        lambda := HomalgZeroMatrix( 1, n, C );
        
        l := Set( [ 1 .. Length( L ) ] );
        
        comb := Combinations( l, n );
        
        z := Random( [ 1 .. Length( comb ) ] );
        
        l := Iterated( List( comb[ z ], i -> L[i] ), UnionOfRows );
        
        comb := Combinations( l, n - 1 );

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
        
    if IsBound( Rad!.WitnessOfExistenceOfZeroDivisor ) then
        I!.WitnessOfExistenceOfZeroDivisor := Rad!.WitnessOfExistenceOfZeroDivisor;
    fi;
    
    if IsBound( Rad!.WitnessForExistenceOfNilpotentElement ) then
        I!.WitnessForExistenceOfNilpotentElement := Rad!.WitnessForExistenceOfNilpotentElement;
    fi;
    
    return bool;
    
end );

##
InstallMethod( PrimaryDecompositionZeroDim,
        [ IsHomalgObject ],

  function( I )
    local Decomp, Rad, a, fac, N, W, i, R, RR, bas, J, j, M;
    
    Decomp := []; 
    
    ## If I is already primary, then the algorithmus returns I itself.
    
    if IsPrimaryZeroDim( I ) then
        Decomp[ 1 ] := I;
        return Decomp;
    fi;
    
    ## If the ideal I is not primary, then the algorithmus asks if IsPrimaryZeroDim has
    ## found a zerodivisor. If yes, then it decomposes the ideal.
    
    if IsBound( I!.WitnessOfExistenceOfZeroDivisor ) then
        a := I!.WitnessOfExistenceOfZeroDivisor;
    fi;
    
    fac := PrimaryDecomposition( LeftSubmodule( MinimalPolynomial( a ) ) );
    fac := List( [ 1 .. Length( fac ) ], i -> MatrixOfSubobjectGenerators( fac[ i ][ 1 ]) );
    fac := List( [ 1 .. Length( fac ) ], i -> MatElm( fac[ i ], 1 ,1 ) );
    
    N := [ ];
    W := [ ];
    
    R := HomalgRing( I );
    RR := R / I;
    
    bas := BasisOverCoefficientsRing( RR );
    
    a := a / RR;
    
    for i in [ 1 .. Length( fac ) ] do
        N[ i ] := RepresentationOverCoefficientsRing( Value( fac[ i ] , a ) );
        W[ i ] := SyzygiesOfRows( N[i] ) * RR;
    od;
    
    J := Iterated( W, UnionOfRows );
    
    M := [ ];
    
    j := [ ];
    j[ 1 ] := 0;
    
    for i in [ 1 .. Length( fac ) ] do
        
        j[ i + 1 ] := j[ i ] + NrRows( W[ i ] );
        M[ i ] := UnionOfRows( CertainRows( J , [ 1 .. j[i] ] ), CertainRows( J, [ j[ i + 1 ] + 1 .. NrRows( J ) ] ) );
        # M[ i ] := IdealBasisToGroebner( ( M[ i ] ) );
        
        M[ i ] := UnionOfRows( MatrixOfGenerators( I ) * R, ( M[ i ] * bas ) * R );
        M[ i ] := LeftSubmodule( BasisOfRows( M[ i ] ) );
        
        if not IsOne( M[ i ] ) then
            Append( Decomp, PrimaryDecompositionZeroDim( M[ i ] ) );
        fi;
    od;
    
    return Decomp;

end );
