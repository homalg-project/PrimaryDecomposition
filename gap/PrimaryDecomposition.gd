#############################################################################
##
##  PrimaryDecomposition.gd                                    PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Declaration of some functions for primary decomposition.
##
#############################################################################

#! @Chapter RadicalComputation

####################################
#
#! @Section Attributes
#
####################################

#! @Description
#!  Determines if the zerodimensional ideal <A>I</A> is a primeideal 
#!  and returns eventually an element, which proofs that <A>I</A> is not prime.
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareAttribute( "IsPrimeZeroDim",
	IsHomalgObject );

#! @Description
#!  Determines if the zerodimensional ideal <A>I</A> is a primaryideal.
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareAttribute( "IsPrimaryZeroDim",
	IsHomalgObject );

