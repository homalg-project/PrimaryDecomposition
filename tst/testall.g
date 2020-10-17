# SPDX-License-Identifier: GPL-2.0-or-later
# PrimaryDecomposition: Tools for primary decomposition
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
options := rec(
    exitGAP := true,
    testOptions := rec(
        compareFunction := "uptowhitespace",
    ),
);

LoadPackage( "PrimaryDecomposition" );
HOMALG.SuppressParityInViewObjForCommutativeStructureObjects := true;
PRIMARY_DECOMPOSITION.RandomSource := RandomSource( IsMersenneTwister );

TestDirectory( DirectoriesPackageLibrary( "PrimaryDecomposition", "tst" ), options );

FORCE_QUIT_GAP( 1 ); # if we ever get here, there was an error
