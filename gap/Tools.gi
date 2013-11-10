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
# methods for properties:
#
####################################

##
InstallMethod( IsIrreducibleHomalgRingElement,
        "for a ring element",
	[ IsHomalgRingElement ],
        
  function( r )
    
    return IsPrimeModule( LeftSubmodule( r ) );
    
end );
    
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
    
##
InstallMethod( SquareFreeFactors,
        "for a ring element",
	[ IsHomalgRingElement ],
        
  function( r )
    local rad;
    
    rad := RadicalDecomposition( LeftSubmodule( r ) );
    rad := List( rad, MatrixOfSubobjectGenerators );
    
    rad := List( rad, m -> MatElm( m, 1, 1 ) );
    
    Assert( 8, Set( List( rad, IsIrreducibleHomalgRingElement ) ) = [ true ] );
    Perform( rad, function( r ) SetIsIrreducibleHomalgRingElement( r, true ); end );
    
    return rad;
    
end );

##
InstallMethod( SeparablePart,
        "for a ring element",
	[ IsHomalgRingElement ],

  function( p )
    local R;
    
    R := CoefficientsRing( HomalgRing( p ) );
    
    if Characteristic( R ) = 0 or IsFinite( R ) then
       return Product( SquareFreeFactors( p ) );
    else
       TryNextMethod( );
    fi;
    
end );

####################################
#
# methods for operations:
#
####################################

#! @Chunk MinimalPolynomial_ring_element
#!  The installed method computes the representation matrix <M>M</M>
#!  of the ring element <A>r</A> and finds the first linear dependency
#!  among the vectors <M> e, eM, eM^2, ... </M>, where <M>e</M>
#!  is the identity of the ring.
#! @EndChunk

InstallMethod( MinimalPolynomial,
        "for two ring elements",
	[ IsHomalgRingElement, IsHomalgRingElement ],
        
  function( r, t )
    local R, bas, M, k, m, n, c, e;
    
    R := HomalgRing( r );
    
    bas := BasisOverCoefficientsRing( R );
    
    if not IsOne( MatElm( bas, 1, 1 ) ) then
        Error( "the one of the algebra is not the first element of the basis\n" );
    fi;
    
    M := RepresentationOverCoefficientsRing( r );
    
    k := HomalgRing( M );
    
    c := NrColumns( M );
    
    ## this relies on the convention that the identity of
    ## the algebra is the first basis vector
    n := CertainRows( HomalgIdentityMatrix( c, k ), [ 1 ] );
    
    #n := HomalgInitialMatrix( 1, c, k );
    #SetMatElm( n, 1, 1, One( k ) );
    
    m := n * M;
    
    while not IsZero( DecideZeroRows( m, n ) ) do
        n := UnionOfRows( n, m );
        m := m * M;
    od;
    
    n := UnionOfRows( n, m );
    
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

##
InstallMethod( Value,
        "for two homalg ring elements",
	[ IsHomalgRingElement, IsHomalgRingElement ],
        
  function( p, v )
    local R, t, coeffs, A, z;
    
    R := HomalgRing( p );
    
    if HasRelativeIndeterminatesOfPolynomialRing( R ) then
        t := RelativeIndeterminatesOfPolynomialRing( R );
    else
        t := Indeterminates( R );
    fi;
    
    if not Length( t ) = 1 then
        TryNextMethod( );
    fi;
    
    coeffs := CoefficientsOfUnivariatePolynomial( p );

    A := HomalgRing( v );
    
    coeffs := A * coeffs;
    coeffs := EntriesOfHomalgMatrix( coeffs );
    
    z := Zero( A );
    
    return Sum( [ 0 .. Length( coeffs ) - 1 ],
                function( i )
                  local a_i;
                  a_i := coeffs[i + 1];
                  if IsZero( a_i ) then
                      return z;
                  fi;
                  return a_i * v^i;
                end );
    
end );

##
InstallMethod( IsIrreducible,
        "for a ring element",
	[ IsHomalgRingElement ],
        
  IsIrreducibleHomalgRingElement );
