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
        
        bas := UnionOfRows( bas, CertainRows( m, Reversed( [ 1 .. NrRows( m ) ] ) ) );
        
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
    
    if HasIsInitialMatrix( mat ) and IsInitialMatrix( mat ) then
        mat := HomalgZeroMatrix( d, d, k );
    fi;
    
    MakeImmutable( mat );
    
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
    
    return SortedList( rad );
    
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

##
InstallMethod( IsPerfect,
        "for a ring",
        [ IsHomalgRing ],
  
  function( R )

    if Characteristic( R ) = 0 or IsFinite( R ) then
        return true;
    else
        return false;
    fi;

end );

##
InstallMethod( BasisCoefficientsOfRingElement, 
	[ IsHomalgRingElement],
	
  function( r )
    local R, M, K, e;
    
    R := HomalgRing( r );
    
    M := RepresentationOverCoefficientsRing( r );
    
    K := CoefficientsRing( AmbientRing ( R ) );
    
    e := CertainRows( HomalgIdentityMatrix( NrRows( M ) , K ) , [ 1 ] );
    
    return ( e * M );
    
end );

##
InstallMethod( IdealBasisOverCoefficientRing,
        [ IsHomalgMatrix ],
        
  function( G )
    local R, M, W, n, T, k, bool, N, j, bas;  
    
    R := HomalgRing( G );
    
    M := FGLMdata( R );
        
    W := List( [ 1 .. NrRows( G ) ], i -> BasisCoefficientsOfRingElement( MatElm( G, i, 1 ) / R ) );
    W := Iterated( W, UnionOfRows );
    
    if IsZero( W ) then
        return HomalgZeroMatrix( 1, 1, R );
    fi;
    
    n := Length( M );
    
    T := List( [ 1 .. NrRows( W ) ], i -> CertainRows( W, [ i ] ) );
    
    for k in [ 1 .. n ] do
        
        bool := 1;
        
        while bool = 1 do
            
            bool := 0;
            
            N := List( [ 1 .. Length( T ) ], i ->  T[i] * M[k] );
            
            T := [ ];
            
            for j in N do
            
                if  not IsZero( DecideZeroRows( j, W ) ) then
                    
                    Add( T, j );
                    
                    bool := 1;
                    
                    W := UnionOfRows( W, j );
                
                fi;
                
            od;
            
        od;
        
        T := List( [ 1 .. NrRows( W ) ], i -> CertainRows( W, [ i ] ) );
        
    od;
    
    bas := BasisOverCoefficientsRing( R );
    
    return Iterated( List( [ 1 .. NrRows( W ) ], i ->  R * CertainRows( W, [ i ] ) * bas ), UnionOfRows );
    
end );

##
InstallMethod( IdealBasisToGroebner, 
	[ IsHomalgMatrix ],

  function( M )
    local R, K, C, Ech, p, d, S, pos, bas, A, el, GJ, j, I, GI, i, lambda;
    
    R := HomalgRing( M );
    
    K := CoefficientsRing( AmbientRing( R ) );
    
    if IsZero( M ) then
        return MatrixOfGenerators( DefiningIdeal( R ) );
    fi;
    
    ## Computes the coefficients with respect to the basis of R.
    C := List( [ 1 .. NrRows( M ) ], i -> BasisCoefficientsOfRingElement( MatElm( M, i, 1 ) ) );
    C := Iterated( C, UnionOfRows );
    
    ## Reverse order of columns of the matrix instead of reverse order of basis
    Ech := CertainColumns( C, Reversed( [ 1 .. NrColumns( C ) ] ) );
    
    ## Regard the matrix over an internal GAP ring to compute the reduced
    ## echelon form
    if IsPerfect( K ) then
        
        p := Characteristic( K );
        
        d := DegreeOverPrimeField( K );
        
        if IsZero( p ) then
            S := HomalgFieldOfRationals( );
        else
            S := HomalgRingOfIntegers( p, d );
        fi;
        
        Ech := BasisOfRows( S * Ech );
        
        Ech := K * Ech;
        
    fi;
    
    ## The reduced echelon form is a reduction of the generators of the ideal.
    ## Now we iterate over the generators starting with the one which has the
    ## smallest leading monomial and add them to the GJ if they
    ## are not a multiply of the leading monomials of GJ.
    ## To proof this we defined the ideal I and the residue class ring S.
    
    C := CertainColumns( Ech, Reversed( [ 1 .. NrColumns( Ech ) ] ) );
    
    pos := PositionOfFirstNonZeroEntryPerRow( Ech );
       
    bas := BasisOverCoefficientsRing( R );
    
    A := AmbientRing( R );
    
    ## el saves the leading monomials of the GJ.
    el := [ Zero( A ) ];
    
    ## GJ is supposed to become the Groebner basis of the ideal.
    GJ := HomalgZeroMatrix( 0, 1, A );
    
    for j in Reversed( [ 1 .. NrRows( Ech ) ] ) do
    
        I := LeftSubmodule( el, A );
        
        S := A / I;
        
        if not IsZero( MatElm( bas, NrRows( bas ) + 1 - pos[ j ] , 1 ) / S ) then
        
            if IsZero( el[1] ) then
                el[1] := MatElm( bas, NrRows( bas ) + 1 - pos[ j ] , 1 ) / A;
            else
                Add( el, MatElm( bas, NrRows( bas ) + 1 - pos[ j ], 1 ) / A );
            fi;
            
            GJ := UnionOfRows( GJ, ( A * CertainRows( C, [j]) ) * ( A * bas ) );
        
        fi;
        
    od;
    
    I := LeftSubmodule( el, A );
    
    S := A / I;
    
    ## At least we have to add the generators of the defining ideal in which
    ## residue class ring the ideal was living in. We only add the generators
    ## whose leading monomial is not a multiply of the leading monomials of GJ.
    
    GI := MatrixOfSubobjectGenerators( DefiningIdeal( R ) );
    
    ## C contains the leading monomials of the generators of GI.
    C := List( [ 1 .. NrRows( GI ) ], i -> Coefficients( MatElm( GI, i, 1 ) )!.monomials[1] );
    
    for i in [ 1 .. NrRows( GI ) ] do
        
        if not IsZero( C[i] / S ) then
            
            ## if the leading monomial of a generator of the defining ideal is 
            ## not a multiply of the leading monomials of GJ, we add the
            ## we add the generator to the GJ.
            
            lambda := BasisCoefficientsOfRingElement( ( MatElm( GI, i, 1 ) - C[i] ) / R  );
            d := NrColumns( lambda );
            
            lambda := CertainColumns( lambda, [d, d - 1 .. 1 ] );
            
            lambda := DecideZeroRows( lambda, Ech );
            
            lambda := CertainColumns( lambda, [d, d - 1 .. 1 ] );
            
            lambda := C[i] + ( MatElm( R * lambda * bas, 1, 1 ) /A );
            
            lambda := HomalgMatrix( [lambda], 1, 1, A );
            
            GJ := UnionOfRows( GJ, lambda );
            
        fi;
        
    od;
    
    return GJ;
    
end );

##
InstallMethod( DerivativeSep,
	"for a ring element",
	[ IsHomalgRingElement ],

  function( p )
    local R, indets, coeffs, monoms, i;
    
    R := HomalgRing( p );
    
    indets := Indeterminates( R );
    
    if Length( indets ) <> 1 then
        TryNextMethod( );
    fi;
    
    indets := indets[1];
    
    coeffs := Coefficients( p );
    
    monoms := coeffs!.monomials;
    
    p := Zero( HomalgRing( p ) );
    
    for i in [ 1 .. NrRows( coeffs ) ] do
        
        if not IsOne( monoms[i] ) then
        
            p := p + MatElm( coeffs, i, 1) * Degree( monoms[i] ) * indets^( Degree( monoms[i] ) - 1 );
        
        fi;
    
    od;
    
    return p;

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

##
InstallMethod( IsNotContainedInAnyHyperplane,
	"for two matrices",
	[ IsHomalgMatrix, IsHomalgMatrix ],

  function( lambda, L )
    local n, l, comb, i, W;
    
    n := NrColumns( L );
    
    l := Set( [ 1 .. NrRows( L ) ] );
    
    comb := Combinations( l, n - 1 );
    
    for i in comb do
        
        W := CertainRows( L, i );
        
        if RowRankOfMatrix( W ) = n - 1 and IsZero( DecideZeroRows( lambda, W ) ) then
            return false;
        fi;
    
    od;
    
    return true;

end );

##
InstallMethod( GeneratorOfAnElementNotContainedInAnyHyperplane,
	"for a matrix",
	[ IsHomalgMatrix ],

  function( L )
    local lambda;
    
    while true do
    
        lambda := Iterated( List( [ 1 .. NrRows( L ) ], i -> Random( [ -10 .. 10 ] ) * CertainRows( L, [i]) ), \+ );
        
        if IsNotContainedInAnyHyperplane( lambda, L ) then
            return lambda;
        fi;
        
    od;
    
end );

##
InstallMethod( FGLMToGroebner,
        "for a list and a matrix",
        [ IsList, IsHomalgMatrix ],
        
  function( M, e )
    local n, A, K, indets, GK, G, B, BB, L, J, S, deg, bool, monoms, j, a, b, k, c, syz, i, l;
    
    ## Creating a Polynomialring with the right number of variables
    n := Length( M );
    
    A := HomalgRing( M[1] );
    K := A * Concatenation( "x1..", String( n ) );
    
    indets := Indeterminates( K );
    
    ## G is supposed to become the GroebnerBasis, B the monomial basis of the 
    ## residue class ring and BB contains all elements in B represented by the
    ## coefficients in the rows of BB.
    G := [ ];
    
    B := [ One( K ) ];
    BB := e;
    
    ## To decide wether a monomial is a multiple of a monomial
    ## in L, we create the residue class K modulo all elements of L.
    ## Since the ideal is generated by monomials, it is very easy to compute a
    ## Groebner basis and therefore easy to decide if an element of the residue
    ## class field is zero or not.
    L := [ Zero( K ) ];
    
    J := LeftSubmodule( L, K );
    S := K / J;
        
    ## To iterate over all monomials (monoms) of K, we create the graded ring
    ## of K and start with the smallest one except one.
    ## Gradually increasing the degree of the monomials we proof if they 
    ## are a multiple of an element of L. If all monomials of a 
    ## certain degree are mulitples of elements of L (that is the case if bool
    ## remains to be 1 in the for loop for), all monomials of higher 
    ## degree will be also multiples of elements of L and we can stop the 
    ## iteration.
    GK := GradedRing( K );
        
    deg := 0;
        
    bool := 0;
    
    while IsZero( bool ) do
        
        deg := deg + 1;
        
        monoms := MonomialMatrix( deg, GK );
        
        bool := 1;
        
        for j in Reversed( [ 1 .. NrRows( monoms ) ] ) do
        
            a := MatElm( monoms, j, 1 ) / K;
            b := e;
            
            ## this if statement is redundant, since we regard every monomial
            ## only one time.
            # if Position( B, a ) = fail then
            
                ## Determines if a is a multiple of an element of L:
                if not IsZero( a / S ) then
                
                    bool := 0;
                    
                    ## Computes the coefficient matrix b of the element a:
                    for k in [ 1 .. n ] do                    
                        
                        c := a;
                        
                        while c / indets[k] <> fail do
                        
                            b := b * M[k];
                            c := c / indets[k];
                        
                        od;
                        
                    od;
                    
                    ## Determines if a represented by b is contained in the 
                    ## k-span of the elements in B represented by the rows of 
                    ## BB. If not, a resp. b is a basis element of the residue
                    ## class fiel and will be added to B resp. BB.
                    ## If yes, the linear combination is a element of the
                    ## Groebner basis.
                    
                    if not IsZero( DecideZeroRows( b, BB ) ) then
                        
                        Add( B, a );
                        BB := UnionOfRows( BB, b );
                    
                    else
                        
                        syz := SyzygiesOfRows( UnionOfRows( BB, b ) );
                        
                        i := 0;
                        l := 1;
                        
                        while IsZero( i ) and ( (l - 1) < NrRows( syz ) ) do
                            if not IsZero( MatElm( syz, l, NrColumns( syz ) ) ) then
                                syz := CertainRows( syz, [l] );
                                i := 1;
                                l := l + 1;
                            fi;
                        od;
                        
                        
                        if IsZero( MatElm( syz, 1, NrColumns( syz ) ) ) then
                            
                            Error();
                        
                        elif not IsOne( MatElm( syz, 1, NrColumns( syz ) ) ) then
                        
                            Error();
                        
                        else
                            ## Every time we change L, J and S have to get
                            ## changed, too.
                            Add( L, a );
                            
                            J := LeftSubmodule( L, K );
                            
                            S := K / J;
                            
                            ## Add the linear combination to the Groebner basis:
                            c :=  MatElm( syz, 1, NrColumns( syz ) ) / K * a;
                            
                            for l in [ 1 .. NrColumns( syz ) - 1 ] do
                            
                                c :=c + ( MatElm( syz, 1, l ) / K ) * B[l];
                            
                            od;
                            
                            Add( G, c );
                            
                            
                        fi;
                    
                    fi;
                    
                fi;
             
            #fi;
        
        od;
    
    od;
    
    return[ B, G ];
    
end );
