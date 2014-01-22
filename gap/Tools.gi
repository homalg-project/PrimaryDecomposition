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

InstallValue( PRIMARY_DECOMPOSITION,
        rec(
              RandomSource := GlobalMersenneTwister,
            )
);

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
        
        ## Use reversed list to get a sorted output.
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
    local R, t;
    
    R := HomalgRing( r );
    
    if HasAmbientRing( R ) then
        R := AmbientRing( R );
    fi;
    
    t := UnusedVariableName( R, "t" );
    
    return MinimalPolynomial( r, t );
    
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
    
    ## We assume that the polynomials generating the radical components are 
    ## irreducible.
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
    
    if IsPerfect( R ) then
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
    local R, A, k, bas, m, d, mat, coeffs, monoms, pos;
    
    R := HomalgRing( r );
    
    A := AmbientRing( R );
    
    k := CoefficientsRing( A );
    
    m := DecideZero( r * One(R) );
    
    bas := EntriesOfHomalgMatrix( A * BasisOverCoefficientsRing( R ) );
    
    d := Length( bas );
    
    mat := HomalgInitialMatrix( 1, d, k );
    
    coeffs := Coefficients( m / A );
    
    monoms := coeffs!.monomials;
    
    coeffs := k * coeffs;
    
    pos := List( monoms, a -> Position( bas, a ) );
    
    Perform( [ 1 .. Length( pos ) ], function( j ) SetMatElm( mat, 1, pos[j], MatElm( coeffs, j, 1 ) / k ); end );
    
    if HasIsInitialMatrix( mat ) and IsInitialMatrix( mat ) then
        mat := HomalgZeroMatrix( 1, d, k );
    fi;
    
    MakeImmutable( mat );
    
    return mat;
    
    
end );

##
InstallMethod( IdealBasisOverCoefficientRing,
        [ IsHomalgMatrix ],
        
  function( G )
    local R, W, M, n, T, k, bool, N, j, bas;  
    
    R := HomalgRing( G );
        
    W := List( [ 1 .. NrRows( G ) ], i -> BasisCoefficientsOfRingElement( MatElm( G, i, 1 ) / R ) );
    W := Iterated( W, UnionOfRows );
    
    ## If W is the zero matrix, then the ideal generated by the
    ## elements of the matrix G is the zero ideal over the residue
    ## class ring.
    if IsZero( W ) then
        return HomalgZeroMatrix( 1, 1, R );
    fi;
    
    M := FGLMdata( R );
    
    n := Length( M );
    
    T := List( [ 1 .. NrRows( W ) ], i -> CertainRows( W, [ i ] ) );
    
    for k in [ 1 .. n ] do
        
        bool := 1;
        
        while bool = 1 do
            
            bool := 0;
            
            N := List( T, i ->  i * M[k] );
            
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
InstallMethod( GapInternalIsomorphicField,
	"for a ring",
	[ IsHomalgRing ],

  function( K )
    local p, d, S;
    
    p := Characteristic( K );
    
    d := DegreeOverPrimeField( K );
    
    if IsZero( p ) then
        S := HomalgFieldOfRationals( );
    else
        S := HomalgRingOfIntegers( p, d );
    fi;
    
    return S;

end );

##
InstallMethod( IdealBasisToGroebner, 
	[ IsHomalgMatrix ],

  function( M )
    local R, K, C, Ech, p, d, S, pos, bas, A, leadingmonoms, GJ, j, I, GI, i, lambda;
    
    R := HomalgRing( M );
    
    if not IsPerfect( CoefficientsRing( R ) ) then
        TryNextMethod( );
    fi;
    
    K := CoefficientsRing( AmbientRing( R ) );
    
    if IsZero( M ) then
        return MatrixOfSubobjectGenerators( DefiningIdeal( R ) );
    fi;
    
    ## Computes the coefficients with respect to the basis of the ring R.
    C := List( [ 1 .. NrRows( M ) ], i -> BasisCoefficientsOfRingElement( MatElm( M, i, 1 ) ) );
    C := Iterated( C, UnionOfRows );
    
    ## Reverse order of the columns of the matrix instead of reversed order
    ## of basis.
    Ech := CertainColumns( C, Reversed( [ 1 .. NrColumns( C ) ] ) ); 

    ## Regard the matrix over an GAP internal field to compute the reduced
    ## echelon form.
    S := GapInternalIsomorphicField( K );
    
    Ech := K * BasisOfRows( S * Ech );
    
    C := CertainColumns( Ech, Reversed( [ 1 .. NrColumns( Ech ) ] ) );
    

    ## A reduction of the generators of the ideal is computed by the
    ## reduced echelon form.
    ## The next part of the algorithms iterates over all generators 
    ## starting with the one which has the smallest leading monomial. The
    ## generator will be added to the set GJ if its leading monomial is not 
    ## a multiply of the leading monomials of GJ.
    ## A monomial is a multiply of the leading monomials of GJ if it is zero
    ## over the residue class ring S. 
    
    pos := PositionOfFirstNonZeroEntryPerRow( Ech );
       
    bas := BasisOverCoefficientsRing( R );
    
    A := AmbientRing( R );
    
    ## The list leadingmonoms stores the leading monomials of GJ, which define
    ## the ideal I of the residue class ring S.
    
    leadingmonoms := [ Zero( A ) ];
    
    ## GJ is supposed to become the Groebner Basis of the ideal.
    
    GJ := HomalgZeroMatrix( 0, 1, A );
    
    for j in Reversed( [ 1 .. NrRows( Ech ) ] ) do
    
        I := LeftSubmodule( leadingmonoms, A );
        
        S := A / I;
        
        if not IsZero( MatElm( bas, NrRows( bas ) + 1 - pos[ j ] , 1 ) / S ) then
        
            if IsZero( leadingmonoms[1] ) then
                leadingmonoms[1] := MatElm( bas, NrRows( bas ) + 1 - pos[ j ] , 1 ) / A;
            else
                Add( leadingmonoms, MatElm( bas, NrRows( bas ) + 1 - pos[ j ], 1 ) / A );
            fi;
            
            GJ := UnionOfRows( GJ, ( A * CertainRows( C, [j]) ) * ( A * bas ) );
        
        fi;
        
    od;
    
    I := LeftSubmodule( leadingmonoms, A );
    
    S := A / I;
    
    ## In the last part the algorithm adds the generators of the 
    ## defining ideal of the residue class ring R to the set GJ, whose 
    ## leading monomial is not a multiply of the leading monomials of GJ. 
    ## In the same step the generators get reduced with respect to GJ. 
        
    GI := MatrixOfSubobjectGenerators( DefiningIdeal( R ) );
    
    ## C contains the leading monomials of GI.
    C := List( [ 1 .. NrRows( GI ) ], i -> Coefficients( MatElm( GI, i, 1 ) )!.monomials[1] );
    
    for i in [ 1 .. NrRows( GI ) ] do
        
        if not IsZero( C[i] / S ) then
            
            lambda := BasisCoefficientsOfRingElement( ( MatElm( GI, i, 1 ) - C[i] ) / R  );
            
            d := NrColumns( lambda );
            
            lambda := CertainColumns( lambda, [ d, d - 1 .. 1 ] );
            lambda := DecideZeroRows( lambda, Ech );
            lambda := CertainColumns( lambda, [ d, d - 1 .. 1 ] );
            
            lambda := C[i] + ( MatElm( R * lambda * bas, 1, 1 ) / A );
            
            lambda := HomalgMatrix( [ lambda ], 1, 1, A );
            
            GJ := UnionOfRows( GJ, lambda );
            
        fi;
        
    od;
    
    ## Here we assume that the Gršbner basis oracle
    ## does not rely on its own sorting of a Gršbner basis.
    Assert( 6, BasisOfRows( GJ ) = GJ );
    
    SetIsBasisOfRowsMatrix( GJ, true );
    
    return GJ;
    
end );

#! @Chunk AppendToGroebner_info
#!  If the coefficients ring <M>K</M> is perfect, the method uses the 
#!  IdealBasisOverCoefficientsRing and IdealBasisToGroebner.
#!  If not, the method uses the usual Groebner basis algorithms, since 
#!  IdealBasisToGroebner only works for perfect coefficient rings.
#! @EndChunk

InstallMethod( AppendToGroebnerBasisOfZeroDimensionalIdeal,
	"for a matrix",
	[ IsHomalgMatrix ],

  function( G )
    local R, I, M;
    
    R := HomalgRing( G );
    
    if IsPerfect( CoefficientsRing( R ) ) then
        return IdealBasisToGroebner( IdealBasisOverCoefficientRing( G ) );
    fi;
    
    ## If the ring R has a non perfect coefficients ring the method uses the
    ## actual Groebner Basis algorithm since we cannot define a non perfect
    ## GAP internal ring, which is needed in the alternative method.
    
    I := MatrixOfSubobjectGenerators( DefiningIdeal( R ) );
        
    M := UnionOfRows( AmbientRing( R ) * G, I );
    
    return BasisOfRows( M );
    
end );

##
InstallMethod( Derivative,
	"for a ring element",
	[ IsHomalgRingElement ],

  function( f )
    local R, indets, coeffs, monoms, i;
    
    R := HomalgRing( f );
    
    indets := Indeterminates( R );
    
    if Length( indets ) <> 1 then
        TryNextMethod( );
    fi;
    
    coeffs := Coefficients( f );
    
    monoms := coeffs!.monomials;
    
    f := Zero( HomalgRing( f ) );
    
    for i in [ 1 .. NrRows( coeffs ) ] do
        
        if not IsOne( monoms[i] ) then
        
            f := f + MatElm( coeffs, i, 1) * Degree( monoms[i] ) * indets[1]^( Degree( monoms[i] ) - 1 );
        
        fi;
    
    od;
    
    return f;

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

##
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
    
    ## This relies on the convention that the identity of
    ## the algebra is the first basis vector.
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
    
    n := NrRows( L );
    
    l := Set( [ 1 .. n ] );
    
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
    local R, lambda, rs;
    
    while true do
        
        ## Generating a random lambda.
        
        R := HomalgRing( L );
        
        lambda := HomalgInitialMatrix( 1, NrColumns( L ), R );
        
        rs := PRIMARY_DECOMPOSITION.RandomSource;
        
        Perform( [ 1 .. NrColumns( L ) ], function(j) SetMatElm( lambda, 1, j , Random( rs, [ - 50 .. 50 ] ) / R ); end );
        
        MakeImmutable( lambda );
        
        if IsNotContainedInAnyHyperplane( lambda, L ) then
            return lambda;
        fi;
        
    od;
    
end );

##
InstallMethod( GeneratorOfAnElementNotContainedInAnyHyperplane,
	"for a matrix",
	[ IsHomalgMatrix , IsHomalgRing ],

  function( L , S )
    local R, lambda, rs, r, s, lambda2;
    
    while true do
        
        ## Generating a random lambda.
        
        R := HomalgRing( L );
        
        ## Todo: involve the rational parameters of R to generate lambda
        lambda := HomalgInitialMatrix( 1, NrColumns( L ), R );
        
        rs := PRIMARY_DECOMPOSITION.RandomSource;
        
        Perform( [ 1 .. NrColumns( L ) ], function(j) SetMatElm( lambda, 1, j , Random( rs, [ - 50 .. 50 ] ) / R ); end );
        
        MakeImmutable( lambda );
        
        if IsBound( R!.RootOfBaseField ) then
            r := R!.RootOfBaseField;
        else
            r := 0;
        fi;
        
        if IsBound( S!.RootOfBaseField ) then
            s := S!.RootOfBaseField;
        else
            s := 0;
        fi;
        
        ## Computing the representation over the new ring S.
        L := CoefficientsTransformation( L, s - r );
        
        lambda2 := CoefficientsTransformation( lambda, s - r );
        
        if IsNotContainedInAnyHyperplane( lambda2, L ) then
            return lambda;
        fi;
        
    od;
    
end );

##
InstallMethod( FGLMToGroebner,
	"for a list and a matrix",
	[ IsList, IsHomalgMatrix ],

  function( M, e )
    local n, R, x, l;
    
    n := Length( M );
    
    R := HomalgRing( M[1] );
    
    if HasAmbientRing( R ) then
        R := AmbientRing( R );
    fi;
    
    x := UnusedVariableName( R, "x" );
    
    if IsOne( n ) then
        return FGLMToGroebner( M, e, [x] );
    else
        ##TODO: list should depend on the UnusedVariable
        return FGLMToGroebner( M, e, Concatenation( "x1..", String( n ) ) );
    fi;
    
end );

##
InstallMethod( FGLMToGroebner,
        "for a list, a matrix and a list",
        [ IsList, IsHomalgMatrix, IsList ],
        
  function( M, e, l )
    local K, indets, G, B, BB, L, J, S, GK, deg, bool, n, monoms, j, a, b,
          c, k, syz, i;
    
    ## The list l should contain as many variables as the list M contains elements.
    K := HomalgRing( M[1] ) * l;
    
    indets := Indeterminates( K );
    
    ## G is supposed to become the Groebner Basis, B the monomial basis of the 
    ## residue class ring and BB contains all elements of B represented by the
    ## coefficients in the rows of BB.
    G := [ ];
    
    B := [ One( K ) ];
    BB := e;
    
    ## To decide whether a monomial is a multiple of a monomial
    ## in L, we create the residue class K modulo all elements of L.
    ## Since the ideal is generated by monomials, it is very easy to compute a
    ## Groebner Basis and therefore easy to decide whether an element of the 
    ## residue class field is zero or not.
    L := [ Zero( K ) ];
    
    J := LeftSubmodule( L, K );
    S := K / J;
        
    ## To iterate over all monomials (monoms) of K, we create the graded ring GK
    ## of K and start with the smallest one except one.
    ## Gradually increasing the degree of the monomials we proof whether they 
    ## are a multiple of an element of L. If all monomials of a 
    ## certain degree are mulitples of elements of L (that is the case if bool
    ## remains to be 1), all monomials of higher degree will be also multiples
    ## of elements of L and we can stop the iteration.
    GK := GradedRing( K );
        
    deg := 0;
        
    bool := 0;
    
    n := Length( M );
    
    while IsZero( bool ) do
        
        deg := deg + 1;
        
        monoms := MonomialMatrix( deg, GK );
        
        bool := 1;
        
        for j in Reversed( [ 1 .. NrRows( monoms ) ] ) do
        
            a := MatElm( monoms, j, 1 ) / K;
            b := e;
            
            ## Determines whether a is a multiple of an element of L.
             if not IsZero( a / S ) then
            
                bool := 0;
                
                c := a;
                
                ## Computes the coefficient matrix b of the element a.
                for k in [ 1 .. n ] do
                
                    while c / indets[k] <> fail do
                   
                        b := b * M[k];
                        c := c / indets[k];
                    
                    od;
                    
                od;
                
                ## Determines whether a represented by b is contained in the 
                ## k-span of the elements in B represented by the rows of 
                ## BB. If not, a resp. b is a basis element of the residue
                ## class fiel and will be added to B resp. BB.
                ## If yes, the linear combination is a element of the
                ## Groebner Basis.
                
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
                   
                    #elif not IsOne( MatElm( syz, 1, NrColumns( syz ) ) ) then
                    ## if not IsOne( MatElm( syz, 1, NrColumns( syz ) ) ) then the
                    ## leading coefficients of the linear combination which 
                    ## will be added to the Groebner Basis is not one. If one
                    ## wants to fix this there will probably appear some
                    ## coefficients depending on rational parameters with
                    ## negativ exponents.
                    
                    else
                        ## Every time we change L, J and S have to get
                        ## changed, too.
                        Add( L, a );
                        
                        J := LeftSubmodule( L, K );
                        
                        S := K / J;
                        
                        ## Add the linear combination to the Groebner basis.
                        c :=  MatElm( syz, 1, NrColumns( syz ) ) / K * a;
                        
                        for l in [ 1 .. NrColumns( syz ) - 1 ] do
                        
                            c :=c + ( MatElm( syz, 1, l ) / K ) * B[l];
                        
                        od;
                        
                        Add( G, c );
                                                 
                    fi;
                
                fi;
                
            fi;
        
        od;
    
    od;
    
    return[ B, G ];
    
end );
