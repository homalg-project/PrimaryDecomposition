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
                           ],
                         ),
        
        Bibliography := "PrimaryDecomposition.bib"
        
);

QUIT;
