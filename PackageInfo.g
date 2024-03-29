# SPDX-License-Identifier: GPL-2.0-or-later
# PrimaryDecomposition: Tools for primary decomposition
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "PrimaryDecomposition",
Subtitle := "Tools for primary decomposition",
Version := "2022.04-01",

Date := ~.Version{[ 1 .. 10 ]},
Date := Concatenation( "01/", ~.Version{[ 6, 7 ]}, "/", ~.Version{[ 1 .. 4 ]} ),
License := "GPL-2.0-or-later",

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Mohamed",
    LastName := "Barakat",
    WWWHome := "https://mohamed-barakat.github.io",
    Email := "mohamed.barakat@uni-siegen.de",
    PostalAddress := Concatenation(
               "Walter-Flex-Str. 3\n",
               "57068 Siegen\n",
               "Germany" ),
    Place := "Siegen",
    Institution := "University of Siegen",
  ),
  rec(
    LastName      := "Hemmerling",
    FirstNames    := "Eva Maria",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "evamaria@hemmerling-nw.de",
    WWWHome       := "https://github.com/HemEM/",
    PostalAddress := Concatenation( [
                       "Department of Mathematics\n",
                       "University of Kaiserslautern\n",
                       "67653 Kaiserslautern\n",
                       "Germany" ] ),
    Place         := "Kaiserslautern",
    Institution   := "University of Kaiserslautern"
  ),
  
],

# BEGIN URLS
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/homalg-project/PrimaryDecomposition",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := "https://homalg-project.github.io/pkg/PrimaryDecomposition",
PackageInfoURL  := "https://homalg-project.github.io/PrimaryDecomposition/PackageInfo.g",
README_URL      := "https://homalg-project.github.io/PrimaryDecomposition/README.md",
ArchiveURL      := Concatenation( "https://github.com/homalg-project/PrimaryDecomposition/releases/download/v", ~.Version, "/PrimaryDecomposition-", ~.Version ),
# END URLS

ArchiveFormats := ".tar.gz .zip",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "PrimaryDecomposition",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Tools for primary decomposition",
),

Dependencies := rec(
  GAP := ">= 4.11.1",
  NeededOtherPackages := [
                   [ "AutoDoc", ">= 2013.11.06" ],
                   [ "RingsForHomalg", ">= 2020.04.17" ],
                   [ "GaussForHomalg", ">= 2013.06.26" ],
                   [ "MatricesForHomalg", ">= 2019.09.01" ],
                   [ "Modules", ">= 2019.09.01" ],
                   [ "GradedModules", "2013.09.30" ],
                   [ "GAPDoc", ">= 1.1" ]
                   ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ]
),

AvailabilityTest := function()
        return true;
    end,

TestFile := "tst/testall.g",

Keywords := [
             "primary decomposition",
             "radical ideal",
             "prime decomposition",
             "associated primes",
             "embedded primes",
             ],

));
