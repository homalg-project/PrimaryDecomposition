#############################################################################
##
##  PrimaryDecomposition.gi                                    PrimaryDecomposition package
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
	" for a zerodimensional homalg ideal",
	[ IsHomalgObject ],
	
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
     od;
     
     RI := RadicalForHomalgIdeal( I );
     
     if not I=RI then
         return ( [ false, [ ] ] );
     fi;

     dRI := NrRows( BasisOverCoefficientsRing( R / RI ) );  
     
     M := HomalgIdentityMatrix( Length( ind ), CoefficientsRing( R ) );

     return(RI);

end );

##
InstallMethod ( IsPrimaryZeroDim,
	"for a zerodimensional homalg ideal",
	[ IsHomalgObject ],
     function( I );
     
        return IsPrime( RadicalForHomalgIdeal( I ) );

end );

