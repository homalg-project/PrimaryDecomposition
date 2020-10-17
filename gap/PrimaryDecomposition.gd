# SPDX-License-Identifier: GPL-2.0-or-later
# PrimaryDecomposition: Tools for primary decomposition
#
# Declarations
#

#! @Chapter Primary Decomposition

####################################
#
#! @Section Properties
#
####################################

#! @Description
#!  Determines if the zerodimensional ideal <A>I</A> is a prime ideal 
#!  and eventually saves an element, which proofs that <A>I</A> is not prime.
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareProperty( "IsPrimeZeroDim",
        IsHomalgModule );
#! @InsertChunk IsPrimeZeroDim

#! @Description
#!  Determines if the zerodimensional ideal <A>I</A> is a primary ideal.
#! @Arguments I
#! @Returns a LeftSubmodule
DeclareProperty( "IsPrimaryZeroDim",
        IsHomalgModule );
#! @InsertChunk IsPrimaryZeroDim

####################################
#
#! @Section Attributes
#
####################################

#! @Description
#!  Computes the primary decomposition of a zerodimensional ideal <A>I</A>.
#! @Arguments I
#! @Returns a list
DeclareAttribute( "PrimaryDecompositionZeroDim",
        IsHomalgObject );
#! @InsertChunk PrimaryDecompositionZeroDim

