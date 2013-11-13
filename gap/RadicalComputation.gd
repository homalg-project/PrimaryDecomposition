#############################################################################
##
##  RadicalComputation.gd                                    PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Declaration of some functions for radical computation.
##
#############################################################################

#! @Chapter RadicalComputation

####################################
#
#! @Section Attributes
#
####################################

#! @Description
#!  Computes the radical of the zerodimensional ideal <A>I</A> of a homalg ring. 
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareAttribute( "RadicalForHomalgIdeal",
	IsHomalgObject ); 
#! @InsertSystem RadicalForHomalgIdeal

#! @Arguments mu
#! @Returns a matrix
#! @Description
#!  Computes the embeddingmatrix for an monic univariate polynomial
#!  over a homalg ring.
DeclareAttribute( "PolynomialEmbeddingmatrix",
        IsHomalgRingElement );
#! @InsertSystem PolynomialEmbeddingmatrix

