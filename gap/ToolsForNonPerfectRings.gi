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

#! @Chunk PolysOverTheSameRing_info
#!  For each polynomial ring <M>R_j</M> containing <A>L_j</A> the number
#!  <M>e_j = R_j</M>!.RootOfBaseField defines the embedding map 
#!  <M>R \hookrightarrow R_j: a \mapsto a; t_i \mapsto t_i^{p^{e_j}}</M>
#!  where <M>p</M> is the characteristic of the base field, <M>a</M> a field
#!  element and <M>t_i</M> are the rational parameters and <M>R</M> the common
#!  base function field. The rings <M>R_j</M> and <M>R</M> have to be
#!  technically identical, such that even the variable names coincide.
#!  The common ring <M>S</M> the method defines is in turn a copy of <M>R</M>
#!  with <M>e_S = S</M>!.RootOfBaseField is the least common multiply of the 
#!  <M>e_j</M>, excluding zero. The polynomials <A>L_j</A> will be represented
#!  over <M>S</M> by taking the $t_i$ to the power of <M>p^{e_S-e_j}</M>.
#! @EndChunk

InstallMethod( PolysOverTheSameRing,
	"for a list",
	[ IsList ],

  function( L )
    local R, i, j, lcm, S, K, indets, coeffs, monoms, deg, C;
    
    if IsOne( Length( L ) ) then
        return L;
    fi;
    
    R := HomalgRing( L[1] );
    
    i := List( [ 1 .. Length( L ) ], i -> HomalgRing( L[i] )!.RootOfBaseField );
    
    j := [ ];
    
    Perform( [ 1 .. Length( L ) ], function( k ) if not IsZero( i[k] ) then Add( j, i[k] ); fi; end );
    
    if IsEmpty( j ) then
        return( L );
    fi;
    
    lcm := Lcm( j );
    
    ## S is the common ring.
    S := CoefficientsRing( R ) * List( Indeterminates( R ), Name );
    S!.RootOfBaseField := lcm;
    
    ## K is necessary to compute the representation of the polynomials over the 
    ## new ring S.
    K := CoefficientsRing( CoefficientsRing( R ) ) * RationalParameters( R );
    
    indets := Indeterminates( K );
    
    for i in [ 1 .. Length( L ) ] do
    
        coeffs := Coefficients( L[i] );
        
        monoms := coeffs!.monomials;
        
        deg := lcm - HomalgRing( L[i])!.RootOfBaseField;
        
        C := HomalgRing( coeffs );
        
        coeffs := CoefficientsTransformation( CoefficientsRing( C ) * coeffs, deg );
        
        coeffs := C * coeffs;
        
        L[ i ] := Involution( HomalgMatrix( monoms, Length( monoms ), 1, S ) ) * ( S * coeffs ); 
        
        L[ i ] := MatElm( L[i], 1, 1 );
        
    od;
    
    return L;
    
end );

#! @Chunk SeparablePart_info
#!  For polynomials over nonperfect rings the method uses an algorithm
#!  of Kemper (see <Cite Key="Kemper" />).
#! @EndChunk

InstallMethod( SeparablePart,
	"for a ring element",
	[ IsHomalgRingElement ],

  function( f )
    local R, h, g1, h1, CoeffsRing, p, S, T, K, x, coeffs, monoms, i, a,
          coeffs2, monoms2, j, b, g2, g3;
    
    ## make sure that the polynomial is univariate and the coefficients ring
    ## is not perfect.
    R := HomalgRing( f );
    
    if not IsBound( R!.RootOfBaseField ) then
            R!.RootOfBaseField := 0;
    fi;
    
    if IsPerfect( CoefficientsRing( R ) ) then
        TryNextMethod( );
    fi;
    
    if Length( Indeterminates( R ) ) <> 1 then
        TryNextMethod( );
    fi;
    
    ## Kemper's algorithm:
    
    ## step 1:
    h := Gcd_UsingCayleyDeterminant( f, Derivative( f ) );
    
    g1 := f / h;
    
    g1 := g1 / R;
    
    ## step 2:
    h1 := Zero( R );
    
    ## step 3:
    while h <> h1 do
        
        h1 := h;
        
        h := Gcd_UsingCayleyDeterminant( h, Derivative( h ) );
    
    od;
    
    ## step 4:
    if IsOne( h ) then
    
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
    T := CoefficientsRing( R ) * List( Indeterminates( R ), Name );
    
    T!.RootOfBaseField := 1 + R!.RootOfBaseField;
    
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
    g2 := SeparablePart( h );
    
    ## step 7:
    h := Product( PolysOverTheSameRing( [ g1, g2 ] ) );
    g3 := SeparablePart( h );
    
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
    
    ## K is necessary to change the ti.
    K := CoefficientsRing( R ) * RationalParameters( R );
    
#! @Chunk MatrixEmbedding_info
#!  The method provides that the generator of the field extension <M>\alpha</M>
#!  coincides withe the negativ of the constant coefficients of the minimal
#!  polynomial.
#! @EndChunk
    
    param := - MatElm( Coefficients( f ), NrRows( Coefficients( f ) ), 1 ) / K;
    
    ## i and j are the variable for iterating over the matrix entries.
    for i in [ 1 .. NrRows( M ) ] do
        
        ## L[i] becomes the i-th row of the new matrix.
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

##
InstallMethod( CoefficientsTransformation,
	"for a matrix and a integer",
	[ IsHomalgMatrix, IsInt ],

  function( M, deg )
    local R, p, K, indets, N, r, c, b, coeffs, monoms, l, a, k, d;
    
    if IsZero( deg ) then
        return M;
    fi;
    
    R := HomalgRing( M );
    
    p := Characteristic( CoefficientsRing( R ) );
    
    K := CoefficientsRing( R ) * RationalParameters( R );
    
    indets := Indeterminates( K );
    
    M := K * M;
    
    N := HomalgInitialMatrix( NrRows( M ), NrColumns( M ), R );
    
    for r in [ 1 .. NrRows( M ) ] do
    
        for c in [ 1 .. NrColumns( M ) ] do
            
            b := Zero( R );
            
            coeffs := Coefficients( MatElm( M, r, c ) );
            
            monoms := coeffs!.monomials;
            
            if not monoms = [ ] then
                
                for l in [ 1 .. NrRows( coeffs ) ] do
                    
                    a := One( R );
                    
                    for k in [ 1 .. Length( indets ) ] do
                    
                        d := 0;
                        
                        while not monoms[l] / indets[k] = fail do
                        
                            monoms[l] := monoms[l] / indets[k];
                            
                            d := d + 1;
                            
                        od;
                        
                        a := a * ( RationalParameters( R )[k] )^( p^deg * d );
                        
                    od;
                    
                    a := a / R;
                    
                    b := b + MatElm( R * coeffs, l, 1 ) * a;
                    
                od;
                
                SetMatElm( N, r, c, b );
            
            else
            
                SetMatElm( N, r, c, Zero( R ) );
            
            fi;
            
        od;
    
    od;
    
    MakeImmutable( N );
    
    return( N );

end );

