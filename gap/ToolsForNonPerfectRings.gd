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
#!  Defines a common polynomial ring and computes the representation of the 
#!  polynomials in the list.
DeclareAttribute( "PolysOverTheSameRing",
        IsList );

#! @Arguments f
#! @Returns a ring element
#! @Description
#!  Computes the separable part of a univariate polynomial over a non perfect 
#!  ring.
DeclareAttribute( "SepUnvollkommen",
        IsHomalgRingElement );

####################################
#
#! @Section Operations
#
####################################

#! @Arguments M, f
#! @Returns a matrix
#! @Description
#!  Computes the embedding of an matrix in another matrix
DeclareOperation( "MatrixEmbedding",
        [ IsHomalgMatrix, IsHomalgRingElement ] );

