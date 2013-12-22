#############################################################################
##
##  ToolsForNonPerfectRings.gd                   PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Declaration of some tools for non perfect rings.
##
#############################################################################

#! @Chapter Tools

####################################
#
#! @Section Attributes
#
####################################

#! @Arguments L
#! @Returns a list
#! @Description
#!  Defines a new polynomial ring <M>S</M>, whose coefficients ring is a root
#!  of the same field than the coefficients ring, in which the polynomials of
#!  the list <A>L</A> are defined. Furthermore a representation of the 
#!  polynomials will be computed. The degree <M>d</M> of the root for a ring <M>R</M>
#!  is stored in <M>R!.RootOfBaseField</M>, which is given as the exponent of the
#!  characeteristic <M>p</M> of the base field <M>d = p^e</M>.
DeclareAttribute( "PolysOverTheSameRing",
        IsList );
#! @InsertSystem PolysOverTheSameRing

####################################
#
#! @Section Operations
#
####################################

#! @Arguments M, f
#! @Returns a matrix
#! @Description
#!  Computes the embedding of the matrix over a field extension.
DeclareOperation( "MatrixEmbedding",
        [ IsHomalgMatrix, IsHomalgRingElement ] );
#! @InsertSystem MatrixEmbedding

#! @Arguments M, deg
#! @Returns a matrix
#! @Description
#!  Raises the rational parameters to the power <M>p^<A>deg</A></M> in the
#!  entries of the matrix <A>M</A>, where p is the characteristic of the
#!  base field.
DeclareOperation( "CoefficientsTransformation",
        [ IsHomalgMatrix, IsInt ] );
#! @InsertSystem CoefficientsTransformation
