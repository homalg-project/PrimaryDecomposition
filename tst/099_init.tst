# SPDX-License-Identifier: GPL-2.0-or-later
# PrimaryDecomposition: Tools for primary decomposition
gap> debug_info_level := InfoLevel( InfoDebug );;
gap> SetInfoLevel( InfoDebug, 0 );;
gap> LoadPackage( "PrimaryDecomposition", false );
true
gap> SetInfoLevel( InfoDebug, debug_info_level );;
gap> HOMALG.SuppressParityInViewObjForCommutativeStructureObjects := true;;
gap> PRIMARY_DECOMPOSITION.RandomSource := RandomSource( IsMersenneTwister );;
