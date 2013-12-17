#############################################################################
##
##  RadicalComputation.gd                       PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Declaration of some functions for radical computation.
##
#############################################################################

#! @Chapter Radical Computation

####################################
#
#! @Section Attributes
#
####################################

#! @Description
#!  Computes some stuff needed in RadicalOfIdeal and IsPrimeZeroDim.
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareAttribute( "PreparationForRadicalOfIdeal",
	IsHomalgModule );

#! @Description
#!  Computes the radical of a zero dimensional ideal <A>I</A>. 
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareAttribute( "RadicalOfIdeal",
	IsHomalgModule );
#! @InsertSystem RadicalOfIdeal

#! @Arguments mu
#! @Returns a matrix
#! @Description
#!  Computes the companion matrix of a monic univariate polynomial
#!  over a univariate ring.
DeclareAttribute( "CompanionMatrix",
        IsHomalgRingElement );
#! @InsertSystem CompanionMatrix
