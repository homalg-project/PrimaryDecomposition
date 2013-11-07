#############################################################################
##
##  Tools.gi                                    PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Implementation of some tools.
##
#############################################################################

####################################
#
# methods for attributes:
#
####################################

##
InstallMethod( BasisOverCoefficientsRing,
        "for a residue class ring",
        [ IsHomalgResidueClassRingRep ],
        
  function( R )
    local A, I, M, S, bas, i, m;
    
    A := AmbientRing( R );
    
    I := DefiningIdeal( R );
    
    M := FactorObject( I );
    
    S := GradedRing( A );
    
    M := GradedModule( M, S );
    
    bas := HomalgZeroMatrix( 0, 1, A );
    
    i := 0;
    
    repeat
        
        m := HomogeneousPartOverCoefficientsRing( i, M );
        m := MatrixOfGenerators( UnderlyingModule( m ) );
        
        bas := UnionOfRows( bas, m );
        
        i := i + 1;
        
    until IsZero( m );
    
    return R * bas;
    
end );

##
InstallMethod( RepresentationOverCoefficientsRing,
        "for a ring element",
        [ IsHomalgRingElement ],
        
  function( r )
    local R, A, k, bas, d, mat, m, i, coeffs, monoms, pos;
    
    R := HomalgRing( r );
    
    A := AmbientRing( R );
    
    k := CoefficientsRing( A );
    
    bas := BasisOverCoefficientsRing( R );
    
    m := DecideZero( r * bas );
    m := EntriesOfHomalgMatrix( m );
    
    bas := EntriesOfHomalgMatrix( A * bas );
    
    d := Length( bas );
    
    mat := HomalgInitialMatrix( d, d, k );
    
    for i in [ 1 .. d ] do
        coeffs := Coefficients( m[i] / A );
        monoms := coeffs!.monomials;
        coeffs := k * coeffs;
        pos := List( monoms, a -> Position( bas, a ) );
        Perform( [ 1 .. Length( pos ) ], function( j ) SetMatElm( mat, i, pos[j], MatElm( coeffs, j, 1 ) / k ); end );
    od;
    
    return mat;
    
end );

##
InstallMethod( FGLMdata,
        "for a residue class ring",
        [ IsHomalgResidueClassRingRep ],
        
  function( R )
    
    return List( Indeterminates( R ), RepresentationOverCoefficientsRing );
    
end );

##
InstallMethod( MinimalPolynomial,
        "for a ring element",
	[ IsHomalgRingElement ],
        
  function( r )
    
    return MinimalPolynomial( r, "t" );
    
end );
    
####################################
#
# methods for operations:
#
####################################

##
InstallMethod( MinimalPolynomial,
        "for two ring elements",
	[ IsHomalgRingElement, IsHomalgRingElement ],
        
  function( r, t )
    local M, k, m, n, c, e, KT, T, p, i, f;
    
    M :=RepresentationOverCoefficientsRing(r);
    
    k :=HomalgRing(M);
    
    c := NrColumns( M );
    
    n := HomalgMatrix( [1], 1, c, k );
    
    m := n * M;
    
    while IsZero(DecideZeroRows( m, BasisOfRows( n ) ) ) = false do
        n := UnionOfRows( n, m );
        m := m * M;
    od;
    
    n := UnionOfRows(n, m);
    
    e := SyzygiesGeneratorsOfRows( n );
    e := HomalgRing( t ) * e;
    e := EntriesOfHomalgMatrix( e );
    
    return Sum( Reversed( [ 1 .. Length( e ) ] ), i -> e[i] * t^(i - 1) );
    
end );

##
InstallMethod( MinimalPolynomial,
        "for a ring element and a string",
	[ IsHomalgRingElement, IsString ],
        
  function( r, t )
    local R, u, kt;
    
    R := HomalgRing( r );
    
    if not IsBound( R!.UnivariatePolynomialRingOverCoefficientRing ) then
        R!.UnivariatePolynomialRingOverCoefficientRing := rec( );
    fi;
    
    u := R!.UnivariatePolynomialRingOverCoefficientRing;
    
    if IsBound( u.(t) ) then
        kt := u.(t);
    else
        kt := CoefficientsRing( R ) * t;
        u.(t) := kt;
    fi;
    
    t := Indeterminates( kt )[1];
    
    return MinimalPolynomial( r, t );
    
end );
