(********************************************************************)
(*                                                                  *)
(* :Title:       ChairTiling                                        *)
(*                                                                  *)
(* :Authors:     Uwe Grimm                                          *)
(*                                                                  *)
(* :Context:     AperiodicTilings`ChairTiling`                      *)
(*                                                                  *)
(* :Version:     1.02                                               *)
(*                                                                  *)
(* :Date:        January 12, 2005                                   *)
(*                                                                  *)
(* :Summary:     Implementation of the inflation procedure          *)
(*               for the chair tiling                               *)
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
(*               W.H. Freeman, New York (1987)                      *)
(*               Figure 10.1.5 on page 523                          *)
(*                                                                  *)
(* :Description: The chair tiling consists only of a single tile,   *)
(*               which however occurs in four different             *)
(*               orientations. In this program, a single tile is    *)
(*               given as a list {tiletype,refpoint}, where         *)
(*               tiletype=0,1,2,3 for the 4 rotated versions of     *)
(*               the basic tile, and refpoint={x,y} contains the    *)
(*               coordinates of a reference point that describes    *)
(*               the position of the tile in space. That is, the    *)
(*               first argument determines the rotational degree    *)
(*               of freedom, the second the translational degree.   *)
(*               The length of the egdes is fixed to be 1 or 2,     *)
(*               respectively.                                      *)
(*                                                                  *)
(*               The inflation procedure is performed together      *)
(*               with a rescaling of the length such that the       *)
(*               edge lengths of the tiles are always the same.     *)
(*               In this way, all calculations can be done using    *)
(*               integer coordinates.                               *)
(*                                                                  *)
(* :Changes:     January 12, 2005 [UG] Version 1.01 -> 1.02:        *)
(*               -> just minor changes in the documentation         *)
(*               May 22, 2000 [UG] Version 1.0 -> 1.01:             *)
(*               -> just minor changes in the documentation         *)
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
(* :Disclaimer:  No guarantuee for correctness of the program and   *)
(*               results obtained with the program is given. Use    *)
(*               at your own risk!                                  *)
(*                                                                  *)
(********************************************************************)

BeginPackage["AperiodicTilings`ChairTiling`"]

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
{tiletype,refpoint}. The four different orientations of tiles may be
given different colours by a list of length four for argument tilecol
(white by default). The colour of edges can be changed with linecol
(black by default), and their width by linewidth (1/200 by default)."

TileCoordinates::usage =
"TileCoordinates[{tiletype,refpoint}] gives the coordinates of the
vertices of the tiletype (0,1,2,3, for the 4 rotated versions of the
basic tile) located at the position refpoint."

TileInflation::usage =
"TileInflation[{tiletype,refpoint}] produces a list of the four tiles
obtained by inflating a tile of type tiletype (0,1,2,3, for the 4
rotated versions of the basic tile) located at the position refpoint."

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[Inflation,
          PlotTiling,
          TileCoordinates,
          TileInflation]

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
      TwoVector]

(********************************************)
(* Scaling factor used in an inflation step *)
(********************************************)

ScaleFactor = 2

(*******************************************************************)
(* Some two-dimensional vectors used to determine vertex positions *)
(*******************************************************************)

TwoVector[num_Integer] :=
  TwoVector[num] =
  Dot[MatrixPower[{{0,-1},
                   {1, 0}},num],{1,1}]

(*****************************************)
(* Coordinates of the vertices of a tile *)
(*****************************************)

TileCoordinates[{tile_Integer,refpoint_List}] :=
  Map[(refpoint+#)&,
      {0,(#3+#4)/2,#4,#1,#2,(#2+#3)/2}]&[TwoVector[tile],
                                         TwoVector[tile+1],
                                         TwoVector[tile+2],
                                         TwoVector[tile+3]]

(******************************)
(* Inflation of a single tile *)
(******************************)

TileInflation[{tile_Integer,refpoint_List}] :=
  {{#1,#4},
   {#2,#4+TwoVector[#2]},
   {#1,#4+TwoVector[#1]},
   {#3,#4+TwoVector[#3]}}&[Mod[tile,4],
                           Mod[tile-1,4],
                           Mod[tile+1,4],
                           ScaleFactor*refpoint]

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
           tilecol_List:Table[GrayLevel[1],{4}],
           linecol_:GrayLevel[0],
           linewidth_:1/200] :=
  Graphics[Table[{tilecol[[i+1]],
                  Map[Polygon,#],
                  linecol,Thickness[linewidth],
                  Map[Line[Join[#,Take[#,2]]]&,#]}&[
                    Map[TileCoordinates,
                        Select[tiling,Mod[First[#],4]==i&]]],
                 {i,0,3}],
           AspectRatio->1]

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
        TileInflation]

(********************************************************************)
(*                                                                  *)
(*          End of package "AperiodicTilings`ChairTiling`"          *)
(*                                                                  *)
(********************************************************************)

EndPackage[]

