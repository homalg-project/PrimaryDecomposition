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

####################################
#
#! @Section Attributes
#
####################################

#! @Description
#!  Computes a basis of the ring <A>R</A>,
#!  provided it is free of finite rank over its coefficients ring.
#! @Arguments R
#! @Returns a &homalg; matrix
DeclareAttribute( "BasisOverCoefficientsRing",
        IsHomalgRing );
#! @InsertSystem BasisOverCoefficientsRing

#! @Description
#!  Computes a representation matrix of the ring element <A>r</A>,
#!  provided it is free of finite rank over its coefficients ring.
#! @Arguments r
#! @Returns a &homalg; matrix
DeclareAttribute( "RepresentationOverCoefficientsRing",
        IsHomalgRingElement );
#! @InsertSystem RepresentationOverCoefficientsRing

#! @Description
#!  Computes the FGLM data of the ring <A>R</A> (see <Cite Key="SJ"/>),
#!  provided it is free of finite rank over its coefficients ring.
#!  The FGLM data of such a ring consists of the representationmatrices
#!  of all basiselements of <A>R</A> computed by BasisOverCoefficientsRing.
#! @Arguments R
#! @Returns a list
DeclareAttribute( "FGLMdata",
        IsHomalgRing );
#! @InsertSystem FGLMdata

####################################
#
# global functions and operations:
#
####################################

# basic operations:

