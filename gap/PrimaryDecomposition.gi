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
          z, lambda, w;
    
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
    ## a prime ideal
    
    indets := Indeterminates( R / I );
    
    mu := List( indets, MinimalPolynomial );
    
    sf_mu := List( mu, SquareFreeFactors );
    
    degI := NrRows( BasisOverCoefficientsRing( R / I ) );
    
    for i in [ 1 .. Length( mu ) ] do
        if Length( sf_mu[i] ) > 1 then
            I!.WitnessOfExistenceOfZeroDivisor := indets[i];
            return false;
        elif Degree( sf_mu[i][1] ) < Degree( mu[i][1] ) then
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
        
    ## The last part of the algorithm computes witness elements
    
    n := Length( indets );
    
    RmodRadI := R / RadI;
    
    degRadI := NrRows( BasisOverCoefficientsRing( RmodRadI ) );
    
    if IsFinite( C ) then
        
        e := HomalgIdentityMatrix( n, C );
        
        L := [];
        
        L := List( [ 1 .. n ], i -> CertainRows( e, [i] ) );
        
        iter := Iterator( e );
        
        lambda := NextIterator( iter );  ## The zero element will be left out.
        
        while true do
        
            lambda := NextIterator( iter );

            if Position( L, lambda )= fail then
                
                lambda := ( ( R / I ) * lambda ) * indets;
                
                mu := MinimalPolynomial( lambda );
            
                if IsIrreducible( mu ) and Degree( mu ) = degRadI then
                
                    return true;
                
                elif not IsIrreducible( mu ) then
                
                    I!.WitnessOfExistenceOfZeroDivisor := lambda;
                    return false;
                
                fi;
                
                Add( L , lambda);
                
            fi;
                        
        od;
        
    fi;
    
    n := Length( indets );
    
    z := Zero( RmodRadI );
    
    while true do
        
        ## lambda has to be changed each time
        lambda := CertainRows( HomalgIdentityMatrix( n, C ), [ 1 ] );
        lambda := EntriesOfHomalgMatrix( RmodRadI * lambda );
        
        w := z;
        
        for i in [ 1 .. n ] do
             w := w + lambda[ i ] * ( indets[ i ] / RmodRadI );
        od;
        
        mu := MinimalPolynomial( w );
        
        if IsIrreducible( mu ) and Degree( mu ) = degRadI then
                
            return true;
            
        elif not IsIrreducible( mu ) then
            
            I!.WitnessOfNonPrimeness := w / R;
            return false;
            
        fi;
        
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
    bas := EntriesOfHomalgMatrix( bas );
    bas := HomalgMatrix( bas, 1, Length( bas ), HomalgRing( bas[ 1 ] ) );
    
    a := a / RR;
    
    for i in [ 1 .. Length( fac ) ] do
        N[ i ] := RepresentationOverCoefficientsRing( Value( fac[ i ] , a ) );
        W[ i ] := bas * (SyzygiesOfColumns( N[i] ) * RR);
    od;
    
    J := Iterated( W, UnionOfColumns );
    
    M := [ ];
    
    j := [ ];
    j[ 1 ] := 0;
    
    for i in [ 1 .. Length( fac ) ] do
        j[ i + 1 ] := j[ i ] + NrColumns( W[ i ] );
        M[ i ] := UnionOfColumns( CertainColumns( J , [ 1 .. j[i] ] ), CertainColumns( J, [ j[ i + 1 ] + 1 .. NrColumns( J ) ] ) );
        ## M[ i ] := PrimaryDecompositionZeroDim( LeftSubmodule( BasisOfColumns( M[ i ] ) ) );
    od;
    
    return Decomp;

end );
