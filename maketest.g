## This file is automatically generated
## Standard maketest.g for the homalg project

LoadPackage( "PrimaryDecomposition" );
LoadPackage( "IO_ForHomalg" );
HOMALG_IO.show_banners := false;
HOMALG_IO.use_common_stream := true;
example_tree := ExtractExamples( DirectoriesPackageLibrary( "PrimaryDecomposition", "doc" )[1]![1], "PrimaryDecomposition.xml", [ ], 500 );
RunExamples( example_tree, rec( compareFunction := "uptowhitespace" ) );
GAPDocManualLab( "PrimaryDecomposition" );
QUIT;
