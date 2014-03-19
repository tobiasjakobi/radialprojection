(********************************************************************)
(*                                                                  *)
(* :Title:       SphinxTiling                                       *)
(*                                                                  *)
(* :Authors:     Uwe Grimm                                          *)
(*                                                                  *)
(* :Context:     AperiodicTilings`SphinxTiling`                     *)
(*                                                                  *)
(* :Version:     1.02                                               *)
(*                                                                  *)
(* :Date:        January 12, 2005                                   *)
(*                                                                  *)
(* :Summary:     Implementation of the inflation procedure          *)
(*               for the sphinx tiling                              *)
(*                                                                  *)
(*               This package is part of a collection of            *)
(*               Mathematica programs that were originally          *)
(*               developed for a HERAEUS summer school on           *)
(*               quasicrystals held at Chemnitz University          *)
(*               of Technology in September 1997. For more          *)
(*               information, we refer to the proceedings:          *)
(*                                                                  *)
(*               J.-B. Suck, M. Schreiber, P. Haeussler,            *)
(*               Quasicrystals: An Introduction to structure,       *)
(*               physical properties, and applications,             *)
(*               Springer, Berlin (2002)                            *)
(*                                                                  *)
(* :References:  B. Gruenbaum and G.C. Shephard,                    *)
(*               Tilings and Patterns,                              *)
(*               W.H. Freeman, New York (1987),                     *)
(*               Figure 10.1.6 on page 523                          *)
(*                                                                  *)
(* :Description: As the chair tiling, the sphinx tiling also        *)
(*               consists of a single shape only, which, however,   *)
(*               shows up in six rotational orientations and        *)
(*               their mirror images. Thus there are altogether     *)
(*               twelve possible orientations of a single tile.     *)
(*               As in ChairTiling.m, a tile is encoded as a        *)
(*               list {tiletype,refpoint}, where now there are      *)
(*               twelve values of tiletype, 1,2,3,4,5,6 and         *)
(*               -1,-2,-3,-4,-5,-6 for the rotated versions of      *)
(*               the basic tile, and their mirror images,           *)
(*               respectively. Two tiles whose tiletypes differ     *)
(*               by a sign are mirror images with respect to the    *)
(*               vertical (y-direction) axis. The second argumant,  *)
(*               refpoint={x,y} contains the coordinates of a       *)
(*               reference point that describes the position of     *)
(*               the tile in space. That is, the first argument     *)
(*               determines the rotational and reflection degrees   *)
(*               of freedom, the second the translational degree.   *)
(*               The length of the egdes is fixed to be 1.          *)
(*                                                                  *)
(*               The inflation procedure is performed together      *)
(*               with a rescaling of the length such that the       *)
(*               edge lengths of the tiles are always the same.     *)
(*               In this way, all calculations can be done using    *)
(*               integer coordinates.                               *)
(*                                                                  *)
(* :Changes:     January 12, 2005 [UG] Version 1.01 -> 1.02:        *)
(*               -> just minor changes in the documentation         *)
(*               July 18, 2000 [UG] Version 1.0 -> 1.01:            *)
(*               -> corrected package name in BeginPackage          *)
(*                                                                  *)
(* :Bugs:        No bugs known so far.                              *)
(*               Please send bug reports to:                        *)
(*               u.grimm@physics.org                                *)
(*                                                                  *)
(* :Copyright:   This package may be copied and redistributed       *) 
(*               freely by anyone, but it may not be sold           *)
(*               commercially. Any changes to the original code     *)
(*               should be documented, or preferably suggested      *)
(*               to the authors of this package to make useful      *)
(*               changes or additions available to the community.   *)
(*                                                                  *)
(* :Disclaimer:  No guarantee for correctness of the program and    *)
(*               results obtained with the program is given. Use    *)
(*               at your own risk!                                  *)
(*                                                                  *)
(********************************************************************)

BeginPackage["AperiodicTilings`SphinxTiling`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

Inflation::usage =
"Inflation[tiling,n] performs an n-fold inflation of tiling which has
to consist of a list of single tiles in the form {tiletype,refpoint}."

PlotTiling::usage = 
"PlotTiling[tiling,tilecol,linecol,linewidth] produces a plot of
tiling which has to consist of a list of single tiles in the form
{tiletype,refpoint}. The six different rotational orientations of
tiles and their mirror images may be given different colours by a list
of length twelve for argument tilecol (white by default). The colour
of edges can be changed with linecol (black by default), and their
width by linewidth (1/200 by default)."

TileCoordinates::usage =
"TileCoordinates[{tiletype,refpoint}] gives the coordinates of the
vertices of the tiletype (1,2,3,4,5,6 and -1,-2,-3,-4,-5,-6 for the 6
rotated versions of the basic tile and its mirror image, respectively)
located at the position refpoint. The coordinates are given in terms
of integer linear combinations of the two basis vectors {1,0} and
{1/2,Sqrt[3]/2} which span the triangular lattice. Use TwoCoordinates
to convert to Cartesian coordinates."

TileInflation::usage =
"TileInflation[{tiletype,refpoint}] produces a list of the four tiles
obtained by inflating a tile of type tiletype (1,2,3,4,5,6 and
-1,-2,-3,-4,-5,-6 for the 6 rotated versions of the basic tile and its
mirror image, respectively) located at the position refpoint."

TwoCoordinates::usage =
"TwoCoordinates[{x,y}] yields Cartesian coordinates for points of the
triangular lattice which are given by their integer components x and y
in direction of the two basis vectors {1,0} and {1/2,Sqrt[3]/2},
respectively."

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[Inflation,
          PlotTiling,
          TileCoordinates,
          TileInflation,
          TwoCoordinates]

(********************************************************************)
(*                                                                  *)
(*                     Start of Private context                     *)
(*                                                                  *)
(********************************************************************)

Begin["`Private`"]

(*********************************)
(* Clearing previous definitions *)
(*********************************)

Clear[Inflation,
      PlotTiling,
      ScaleFactor,
      TileCoordinates,
      TileInflation,
      TileMod,
      TwoCoordinates,
      TwoVector];

(********************************************)
(* Scaling factor used in an inflation step *)
(********************************************)

ScaleFactor = 2

(*******************************************************************)
(* Some two-dimensional vectors used to determine vertex positions *)
(*******************************************************************)

TwoVector[0] = { 1, 0}
TwoVector[1] = { 0, 1}
TwoVector[2] = {-1, 1}
TwoVector[3] = -TwoVector[0]
TwoVector[4] = -TwoVector[1]
TwoVector[5] = -TwoVector[2]
TwoVector[num_Integer] := 
  TwoVector[num] = 
  TwoVector[Mod[num,6]]

TwoCoordinates[{xnum_Integer,ynum_Integer}] :=
  {xnum+ynum/2,Sqrt[3]*ynum/2}

(*******************************************)
(* Modified Mod-function to identify tiles *)
(*******************************************)

TileMod[num_Integer] := 
  TileMod[num] = Sign[num]*(#+6*Floor[(6-#)/6])&[Mod[Abs[num],6]]

(*****************************************)
(* Coordinates of the vertices of a tile *)
(*****************************************)

TileCoordinates[{tile_Integer,refpoint_List}] :=
  If[tile>0,
     Map[(refpoint+#)&,
         {0,#2,#1+#2,#1+2*#2,3*#1}]&[TwoVector[tile-1],
                                     TwoVector[tile]],
     Map[(refpoint+#)&,
         {0,#2,#1+#2,#1+2*#2,3*#1}]&[TwoVector[tile+4],
                                     TwoVector[tile+3]],
     Print["Error: Argument in function TileCoordinates is ",
           "neither positive nor negative"]]

(******************************)
(* Inflation of a single tile *)
(******************************)

TileInflation[{tile_Integer,refpoint_List}] :=
  If[tile>0,
     {{#6,ScaleFactor*(#7+3*TwoVector[#1]/2)},
      {#6,ScaleFactor*(#7+3*TwoVector[#1])},
      {#4,ScaleFactor*(#7+3*TwoVector[#1]+TwoVector[#3]/2)},
      {#5,ScaleFactor*(#7+TwoVector[#2])}}&[TileMod[tile-1],
                                            TileMod[tile],
                                            TileMod[tile+1],
                                            TileMod[tile+2],
                                            TileMod[tile-11],
                                            TileMod[tile-8],
                                            refpoint],
     {{#5,ScaleFactor*(#6+3*TwoVector[#1]/2)},
      {#5,ScaleFactor*(#6+3*TwoVector[#1])},
      {#1,ScaleFactor*(#6+3*TwoVector[#1]+TwoVector[#3]/2)},
      {#4,ScaleFactor*(#6+TwoVector[#2])}}&[TileMod[tile-2],
                                            TileMod[tile-3],
                                            TileMod[tile-4],
                                            TileMod[tile+11],
                                            TileMod[tile+8],
                                            refpoint],
     Print["Error: Argument in function TileInflation is ",
           "neither positive nor negative"]]

(*******************************)
(* Inflation of a set of tiles *)
(*******************************)

Inflation[tiling_List,
          num_Integer:1] /; num>=0 :=
  Nest[Flatten[Map[TileInflation,#],1]&,tiling,num]

(**************************************)
(* Graphical presentation of a tiling *)
(**************************************)

PlotTiling[tiling_List,
           tilecol_List:Table[GrayLevel[1],{12}],
           linecol_:GrayLevel[0],
           linewidth_:1/200] :=
  Graphics[Join[Table[{tilecol[[i]],
                       Map[Polygon,#],
                       linecol,Thickness[linewidth],
                       Map[Line[Join[#,Take[#,2]]]&,#]}&[
                       Map[TwoCoordinates,
                           Map[TileCoordinates,
                               Select[tiling,First[#]==i&]],{2}]],
                      {i,6}],
                Table[{tilecol[[6-i]],
                       Map[Polygon,#],
                       linecol,Thickness[linewidth],
                       Map[Line[Join[#,Take[#,2]]]&,#]}&[
                       Map[TwoCoordinates,
                           Map[TileCoordinates,
                               Select[tiling,First[#]==i&]],{2}]],
                      {i,-1,-6,-1}]],
           AspectRatio->Automatic]

(********************************************************************)
(*                                                                  *)
(*                      End of Private context                      *)
(*                                                                  *)
(********************************************************************)

End[]

(********************************************************************)
(*                                                                  *)
(*                    Protect exported functions                    *)
(*                                                                  *)
(********************************************************************)

Protect[Inflation,
        PlotTiling,
        TileCoordinates,
        TileInflation,
        TwoCoordinates]

(********************************************************************)
(*                                                                  *)
(*          End of package "AperiodicTilings`SphinxTiling`"         *)
(*                                                                  *)
(********************************************************************)

EndPackage[]

