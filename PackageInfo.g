SetPackageInfo( rec(

PackageName := "PrimaryDecomposition",

Subtitle := "Tools for primary decomposition",

Version := Maximum( [
                   "2014.09.15", ## Mohamed's version
                   ## this line prevents merge conflicts
                   "2013.11.02", ## Eva Maria's version
                   ] ),

# this avoids git-merge conflicts
Date := ~.Version{[ 1 .. 10 ]},
Date := Concatenation( ~.Date{[ 9, 10 ]}, "/", ~.Date{[ 6, 7 ]}, "/", ~.Date{[ 1 .. 4 ]} ),

ArchiveURL := Concatenation( "http://homalg.math.rwth-aachen.de/~barakat/homalg-project/PrimaryDecomposition-", ~.Version ),

ArchiveFormats := ".tar.gz",

Persons := [
  rec( 
    LastName      := "Barakat",
    FirstNames    := "Mohamed",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "mohamed.barakat@rwth-aachen.de",
    WWWHome       := "http://www.mathematik.uni-kl.de/~barakat/",
    PostalAddress := Concatenation( [
                       "Lehrstuhl B fuer Mathematik, RWTH Aachen\n",
                       "Templergraben 64\n",
                       "52062 Aachen\n",
                       "Germany" ] ),
    Place         := "Aachen",
    Institution   := "RWTH Aachen University"
  ),
  rec(
    LastName      := "Hemmerling",
    FirstNames    := "Eva Maria",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "evamaria@hemmerling-nw.de",
    WWWHome       := "",
    PostalAddress := Concatenation( [
                       "Department of Mathematics\n",
                       "University of Kaiserslautern\n",
                       "67653 Kaiserslautern\n",
                       "Germany" ] ),
    Place         := "Kaiserslautern",
    Institution   := "University of Kaiserslautern"
  ),
  
],

Status := "dev",

README_URL := 
  "http://homalg.math.rwth-aachen.de/~barakat/homalg-project/PrimaryDecomposition/README.PrimaryDecomposition",
PackageInfoURL := 
  "http://homalg.math.rwth-aachen.de/~barakat/homalg-project/PrimaryDecomposition/PackageInfo.g",

PackageDoc := rec(
  BookName  := ~.PackageName,
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := ~.Subtitle,
  Autoload  := false
),

Dependencies := rec(
  GAP := ">=4.4",
  NeededOtherPackages := [
                   [ "AutoDoc", ">= 2013.11.06" ],
                   [ "RingsForHomalg", ">= 2013.08.22" ],
                   [ "GaussForHomalg", ">= 2013.06.26" ],
                   [ "MatricesForHomalg", ">= 2013.12.18" ],
                   [ "Modules", ">= 2013.09.14" ],
                   [ "GradedModules", "2013.09.30" ],
                   [ "GAPDoc", ">= 1.1" ]
                   ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ]
),

AvailabilityTest := function( )
    return true;
  end,

Autoload := false,

Keywords := [
             "primary decomposition",
             "radical ideal",
             "prime decomposition",
             "associated primes",
             "embedded primes" ]

));
