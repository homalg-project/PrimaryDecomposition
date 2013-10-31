#############################################################################
##
##  Tools.gd                                    PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Declaration of some tools.
##
#############################################################################

#! @Chapter Tools

#! @Section Attributes

####################################
#
# attributes:
#
####################################

#! @Arguments R
#! @Returns a &homalg; matrix
#! @Description
#!  Computes a basis of the ring <A>R</A>,
#!  provided it is free of finite rank over its coefficients ring.
DeclareAttribute( "BasisOverCoefficientsRing",
        IsHomalgRing );
#! @InsertSystem BasisOverCoefficientsRing

#! @Arguments r
#! @Returns a &homalg; matrix
#! @Description
#!  Computes a representation matrix of the ring element <A>r</A>,
#!  provided it is free of finite rank over its coefficients ring.
DeclareAttribute( "RepresentationOverCoefficientsRing",
        IsHomalgRingElement );

#! @Arguments R
#! @Returns a list
#! @Description
#!  Computes the FGLM data of the ring <A>R</A> (see <Cite Key="SJ"/>),
#!  provided it is free of finite rank over its coefficients ring.
DeclareAttribute( "FGLMdata",
        IsHomalgRing );

####################################
#
# global functions and operations:
#
####################################

# basic operations:

