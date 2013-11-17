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
    local R, ind, dI, mu, fac, i, RI, dRI, M;
    
    R := HomalgRing( I );
    
    ind := Indeterminates( R / I );

    dI := NrRows( BasisOverCoefficientsRing( R / I ) );

    mu := List( ind, MinimalPolynomial );

    fac := List( mu, SquareFreeFactors );

    for i in [1 .. Length(mu)] do
      if Length( fac[i] ) > 1 then
         return( [ false, [ ind[i] ] ]);
   
      elif ( Length( fac[i] ) = 1 and Degree(fac[1][1]) = dI ) then
               return ( [ true , [ ] ] );
      
      fi;
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
    
    RI := RadicalOfIdeal( I );
    ## Now the algorithms asks if I is a radical ideal.
    ## If yes, it is a prime ideal. If not, then I cannot be a prime ideal.
    
    if not I=RI then
       return ( [ false, [ ] ] );
    fi;

    dRI := NrRows( BasisOverCoefficientsRing( R / RI ) );
    
    M := HomalgIdentityMatrix( Length( ind ), CoefficientsRing( R ) );

    return(RI);

    ## The last part of the algorithm computes witness elements
    
end );

##
InstallMethod ( IsPrimaryZeroDim,
	"for a zerodimensional homalg ideal",
	[ IsFinitelyPresentedSubmoduleRep and ConstructedAsAnIdeal ],
        
  function( I );
    
    return IsPrimeZeroDim( RadicalOfIdeal( I ) );
    
end );
