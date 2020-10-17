# SPDX-License-Identifier: GPL-2.0-or-later
# PrimaryDecomposition: Tools for primary decomposition
#
# This file tests if the package can be loaded without errors or warnings.
#
gap> package_loading_info_level := InfoLevel( InfoPackageLoading );;
gap> SetInfoLevel( InfoPackageLoading, PACKAGE_ERROR );;
gap> LoadPackage( "PrimaryDecomposition", false );
true
gap> LoadPackage( "IO_ForHomalg", false );
true
gap> SetInfoLevel( InfoPackageLoading, PACKAGE_INFO );;
gap> LoadPackage( "PrimaryDecomposition" );
true
gap> LoadPackage( "IO_ForHomalg" );
true
gap> SetInfoLevel( InfoPackageLoading, package_loading_info_level );;
gap> HOMALG_IO.show_banners := false;;
