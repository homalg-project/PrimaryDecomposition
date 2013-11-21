#############################################################################
##
##  PrimaryDecomposition.gd                     PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Declaration of some functions for primary decomposition.
##
#############################################################################

#! @Chapter Primary Decomposition

####################################
#
#! @Section Properties
#
####################################

#! @Description
#!  Determines if the zerodimensional ideal <A>I</A> is a primeideal 
#!  and eventually saves an element, which proofs that <A>I</A> is not prime.
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareProperty( "IsPrimeZeroDim",
	IsHomalgModule );
#! @InsertSystem IsPrimeZeroDim

#! @Description
#!  Determines if the zerodimensional ideal <A>I</A> is a primaryideal.
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareProperty( "IsPrimaryZeroDim",
	IsHomalgModule );
#! @InsertSystem IsPrimaryZeroDim

####################################
#
#! @Section Attributes
#
####################################

#! @Description
#!  Computes the primary decomposition of a zerodimensional ideal.
#! @Arguments I
#! @Returns a list
DeclareAttribute( "PrimaryDecompositionZeroDim",
	IsHomalgObject );
#! @InsertSystem PrimaryDecompositionZeroDim

