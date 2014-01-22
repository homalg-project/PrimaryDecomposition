#############################################################################
##
##  ToolsForNonPerfectRings.gd                   PrimaryDecomposition package
##
##  Copyright 2013,      Mohamed Barakat, University of Kaiserslautern
##                  Eva Maria Hemmerling, University of Kaiserslautern
##
##  Declaration of some tools for non perfect rings.
##
#############################################################################

#! @Chapter Tools

####################################
#
#! @Section Attributes
#
####################################

#! @Description
#!  Accepts a list <A>L</A> of polynomials which are defined over function
#!  fields, whose coefficients rings are a root of the same function field.
#!  The method computes a polynomial ring over a function field and the 
#!  representations of the polynomials in <A>L</A> over this ring.
#!  @InsertChunk PolysOverTheSameRing_info
#! @Arguments L
#! @Returns a list
DeclareAttribute( "PolysOverTheSameRing",
        IsList );
#! @InsertSystem PolysOverTheSameRing

####################################
#
#! @Section Operations
#
####################################

#! @Description
#!  Computes the image of the matrix <M><A>M</A>\in K(\alpha)^{r \times c})</M>
#!  where <M>K(\alpha)</M> with minimal polynomial <A>f</A>, under the  
#!  extension field <M>K(\alpha)</M> with minimal polynomial <A>f</A>
#!  embedding <M>K(\alpha)^{r \times c } \hookrightarrow
#!  K^{r deg(<A>f</A>) \times c deg(<A>f</A>)}</M> which is a the natural
#!  extension of <M>K(\alpha )\hookrightarrow 
#!  K^{deg(<A>f</A>)\times deg(<A>f</A>)}: \alpha \mapsto \texttt{CompanionMatrix}(\alpha)</M>.
#! @InsertChunk MatrixEmbedding_info
#! @Arguments M, f
#! @Returns a matrix
DeclareOperation( "MatrixEmbedding",
        [ IsHomalgMatrix, IsHomalgRingElement ] );
#! @InsertSystem MatrixEmbedding

#! @Arguments M, deg
#! @Returns a matrix
#! @Description
#!  Raises the rational parameters to the power <M>p^<A>deg</A></M> in the
#!  entries of the matrix <A>M</A>, where p is the characteristic of the
#!  base field.
DeclareOperation( "CoefficientsTransformation",
        [ IsHomalgMatrix, IsInt ] );
#! @InsertSystem CoefficientsTransformation

