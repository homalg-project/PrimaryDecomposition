#############################################################################
##
##  ToolsForNonPerfectRings.gi                   PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Implementation of some tools for non perfect rings.
##
#############################################################################
    
####################################
#
# methods for attributes:
#
####################################

##
InstallMethod( PolysOverTheSameRing,
	"for a list",
	[ IsList ],

  function( Sep )
    local R, Rings, S, lcm, K, indets, i, coeffs, r, cocoeffs, b, j, a, monoms, k, deg;
    
    R := HomalgRing( Sep[1] );
    
    Rings := List( [ 1 .. Length( Sep ) ], i -> HomalgRing( Sep[i] )!.RootOfBaseField );
    S := [ ];
    
    ## if any polynomial is the one or if the related ring is not a root of the
    ## base field we left
    for i in [ 1 .. Length( Rings ) ] do
        if not IsZero( Rings[i] ) and not IsOne( Sep[i] ) then
            Add( S, Rings[i] );
        fi;
    od;
    
    lcm := Lcm( S ); 
    
    ## S is the common ring.
    S := R;
    S!.RootOfBaseField := lcm;
    
    ## K is necessary to compute the representation of the polynomials over the 
    ## new ring S. 
    K := CoefficientsRing( CoefficientsRing( R ) ) * RationalParameters( R );
    
    indets := Indeterminates( K );
    
    for i in [ 1 .. Length( Sep ) ] do
    
        coeffs := Coefficients( Sep[i] );
        
        Sep[i] := Zero( S );
        
        for r in [ 1 .. NrRows( coeffs ) ] do
        
            cocoeffs := Coefficients( MatElm( K * coeffs, r, 1 ) );
            
            b := Zero( S );
                        
            for j in [ 1 .. NrRows( cocoeffs ) ] do
                
                a := One( S );
                
                monoms := cocoeffs!.monomials[j];
                
                for k in [ 1 .. Length( indets ) ] do
                
                    deg := 0;
                
                    while not monoms / indets[k] = fail do
                    
                        monoms := monoms/indets[k];
                        
                        deg := deg + 1;
                    
                    od;
                    
                    if not IsZero( deg ) then
                        a := a * RationalParameters( S )[k]^( lcm * deg );
                    fi;
                                    
                od;
            
                a := a / S;
                b := b + MatElm( S * cocoeffs, j, 1 ) * a;
                        
            od;
            
            Sep[i] := Sep[i] + ( b / S ) * ( coeffs!.monomials[r] / S );
            
        od;
        
    od;
    
    return Sep;
    
end );

