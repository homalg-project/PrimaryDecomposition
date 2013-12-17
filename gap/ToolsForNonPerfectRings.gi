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
    local R, Rings, S, lcm, K, indets, i, coeffs, r, cocoeffs, b, j, a, monoms,
          k, deg;
    
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

##
InstallMethod( SepUnvollkommen,
	"for a ring element",
	[ IsHomalgRingElement ],

  function( f )
    local R, h, g1, h1, CoeffsRing, p, S, T, K, x, coeffs, monoms, i, a,
          coeffs2, monoms2, j, b, g2, g3;
    ## Kemper's algorithm:
    
    ## make sure that the polynomial is univariate
    R := HomalgRing( f );
        
    if Length( Indeterminates( R ) ) <> 1 then
        TryNextMethod( );
    fi;
        
    ## step 1:
    h := Gcd_UsingCayleyDeterminant( f, DerivativeSep( f ) );
    
    g1 := f / h;
        
    ## step 2:
    h1 := Zero( R );
    
    ## step 3:
    while h <> h1 do
        
        h1 := h;
        
        h := Gcd_UsingCayleyDeterminant( h, DerivativeSep( h ) );
    
    od;
    
    ## step 4:
    if Degree( h ) = 1 then
        
        if not IsBound( R!.RootOfBaseField ) then
            R!.RootOfBaseField := 0;
        fi;
        
        return g1;
    
    fi;
    
    ## step 5:
    ## write h( x ) = h( x^p ), replace the ti and the ai by their p-th root 
    
    ## Characteristic of the base field:
    CoeffsRing := CoefficientsRing( CoefficientsRing( R ) );
    
    p := Characteristic( CoeffsRing );
    
    S := CoeffsRing * RationalParameters( R );
    
    ## T is the new ring. The RootOfBaseField determines the power of the 
    ## characteristic p the ti of the base ring are send to the ti in the ring T.
    T := R;
    
    if IsBound( T!.RootOfBaseField ) then
        T!.RootOfBaseField := 1 + T!.RootOfBaseField;
    else
        T!.RootOfBaseField := 1;
    fi;
    
    ## For Computation of the p-th root of elements of the perfect base field.
    K := CoefficientsRing( CoefficientsRing( R ) ) * Indeterminates( R );
    
    x := Indeterminates( K )[1];
    
    ## coefficients of h are still polynomials in the rational parameters
    coeffs := Coefficients( h );
        
    monoms := coeffs!.monomials;
    
    h := Zero( T );
    
    for i in [ 1 .. NrRows( coeffs ) ] do
    
        a := MatElm( coeffs, i, 1 ) / S;
        
        ## coefficients of a depend only on the base field.
        coeffs2 := Coefficients( a );
        
        monoms2 := coeffs2!.monomials;
        
        a := Zero( T );
        
        for j in [ 1 .. NrRows( coeffs2 ) ] do
            
            ## Computing the p-th root of the coefficients
            b := x^p - MatElm( coeffs2, j, 1 ) / K;
            
            b := SquareFreeFactors( b )[1];
            
            b := - MatElm( Coefficients( b ), 2, 1 );
            
            b := b / T;
            monoms2[j] := monoms2[j] / T;
            
            ## Putting the polynomial together:
            a := a + b / T * monoms2[j] / T;
            
        od;
        
        ## writing h(x) = h(x^p).
        h := h + a * Indeterminates( T )[1]^( Degree( monoms[i] ) / p ) / T;
        
    od;
    
    ## step 6:
    g2 := SepUnvollkommen( h );
    
    ## step 7:
    h := Product( PolysOverTheSameRing( g1, g2) );
    g3 := SepUnvollkommen( g1 * g2 );
    
    return g3;

end );

####################################
#
# methods for operations:
#
####################################

##
InstallMethod( MatrixEmbedding,
        [ IsHomalgMatrix, IsHomalgRingElement ],
        
  function( M, f )
    local R, N, L, K, param, i, j, elm, S, coeffs, monoms, k, deg;
        
    R := HomalgRing( M );
    
    N := CompanionMatrix( f );
    
    L := [];
    
    ## K is necessary to get change the ti.
    K := CoefficientsRing( R ) * RationalParameters( R );
    
    param := - MatElm( Coefficients( f ), NrRows( Coefficients( f ) ), 1 ) / K;
    
    ## i and j are the variable for iterating over the matrix entries.
    for i in [ 1 .. NrRows( M ) ] do
        
        L[i] := HomalgZeroMatrix( NrRows( N ) , 0, R );
        
        for j in [ 1 .. NrColumns( M ) ] do
            
            elm := MatElm( M, i, j ) / K;
            
            ## S is the Matrix which replaces the matrix entry elm.
            S := HomalgZeroMatrix( NrRows( N ), NrRows( N ), R );
            
            if IsZero( elm ) then
                L[i] := UnionOfColumns( L[i], S );
            
            else
                coeffs := Coefficients( elm );
                monoms := coeffs!.monomials;
                
                for k in [ 1 .. NrRows( coeffs ) ] do
                    
                    deg := 0;
                    
                    while not monoms[k]/param = fail do
                        
                        monoms[k] := monoms[k] / param;
                        
                        deg := deg + 1;
                    
                    od;
                    
                    elm := ( monoms[k] / R ) * ( MatElm( coeffs, k, 1 ) / R )* ( R * N^deg );
                    
                    S := S + elm;
                
                od;
                
                L[i] := UnionOfColumns( L[i], S );
            
            fi;
            
        od;
        
    od;
    
    return Iterated( L, UnionOfRows );
    
end );


