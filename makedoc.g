LoadPackage( "AutoDoc" );

AutoDoc( "PrimaryDecomposition" :
        
        scaffold := true,
                    
        autodoc := true,

        EntityList :=
                    [ "<!ENTITY see '<Alt Only=\"LaTeX\">$\\to$</Alt><Alt Not=\"LaTeX\">--&gt;</Alt>'>",
                      "<!ENTITY GAP4 '<Package>GAP4</Package>'>",
                      "<!ENTITY Singular '<Package>Singular</Package>'>",
                      "<!ENTITY homalg '<Package>homalg</Package>'>",
                      ],
       
        Bibliography := "PrimaryDecomposition.bib"
);

QUIT;
