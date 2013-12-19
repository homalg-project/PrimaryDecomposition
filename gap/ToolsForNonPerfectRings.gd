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

#! @Arguments Sep
#! @Returns a list
#! @Description
#!  Defines a common polynomial ring of polynomials who have the same base field
#!  and computes the representation of the polynomials in the list.
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
#!  Computes the embedding of an matrix over a field extension.
DeclareOperation( "MatrixEmbedding",
        [ IsHomalgMatrix, IsHomalgRingElement ] );
#! @InsertSystem MatrixEmbedding

#! @Arguments M, deg
#! @Returns a matrix
#! @Description
#!  Defines a common polynomial ring of polynomials who have the same base field
#!  and computes the representation of the polynomials in the list.
DeclareOperation( "CoefficientsTransformation",
        [ IsHomalgMatrix, IsInt ] );
#! @InsertSystem CoefficientsTransformation
