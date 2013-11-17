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
    local R, C, indets, mu, sf_mu, degI, i, RadI, RmodRadI, degRadI, n,
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
    
    RmodRadI := R / RadI;
    
    ## The last part of the algorithm computes witness elements
    
    degRadI := NrRows( BasisOverCoefficientsRing( RmodRadI ) );  
    
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
    
    return IsPrimeZeroDim( RadicalOfIdeal( I ) );
    
end );
