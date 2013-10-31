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
        [ IsHomalgResidueClassRingRep ],
        
  function( R )
    
    return List( Indeterminates( R ), RepresentationOverCoefficientsRing );
    
end );
