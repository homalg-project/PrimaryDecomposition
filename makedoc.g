LoadPackage( "AutoDoc" );

AutoDoc( "PrimaryDecomposition" :
        
        scaffold := rec( entities := [ "homalg", "GAP4" ],
                         ),
        
        autodoc := true,
        
        maketest := rec( folder := ".",
                         commands :=
                         [ "LoadPackage( \"PrimaryDecomposition\" );",
                           "LoadPackage( \"IO_ForHomalg\" );",
                           "HOMALG_IO.show_banners := false;",
                           "HOMALG_IO.suppress_PID := true;",
                           "HOMALG_IO.use_common_stream := true;",
                           "HOMALG.SuppressParityInViewObjForCommutativeStructureObjects := true;",
                           "PRIMARY_DECOMPOSITION.RandomSource := RandomSource( IsMersenneTwister );",
                           ],
                         ),
        
        Bibliography := "PrimaryDecomposition.bib"
        
);

# Create VERSION file for "make towww"
PrintTo( "VERSION", PackageInfo( "PrimaryDecomposition" )[1].Version );

QUIT;
