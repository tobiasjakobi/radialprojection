(********************************************************************)
(*                                                                  *)
(* :Title:       PinwheelTiling                                     *)
(*                                                                  *)
(* :Authors:     Uwe Grimm                                          *)
(*                                                                  *)
(* :Context:     AperiodicTilings`PinwheelTiling`                   *)
(*                                                                  *)
(* :Version:     1.0                                                *)
(*                                                                  *)
(* :Date:        January 12, 2005                                   *)
(*                                                                  *)
(* :Summary:     Implementation of the inflation procedure          *)
(*               for the pinwheel tiling                            *)
(*                                                                  *)
(*               This package is part of a collection of            *)
(*               Mathematica programs that were originally          *)
(*               developed for a HERAEUS summer school on           *)
(*               quasicrystals held at Chemnitz University          *)
(*               of Technology in September 1997. For more          *)
(*               information, we refer to the proceedings:          *)
(*                                                                  *)
(*               J.-B. Suck, M. Schreiber, P. Haeussler,            *)
(*               Quasicrystals - An Introduction to structure,      *)
(*               physical properties and application of             *)
(*               quasicrystalline alloys                            *)
(*                                                                  *)
(* :References:  C. Radin,                                          *)
(*               Miles of Tiles,                                    *)
(*               AMS, Providence, Rhode Island (1999)               *)
(*               Figure 30                                          *)
(*                                                                  *)
(* :Description: The pinwheel tiling consists only of a single      *)
(*               tile, which however occurs in infinitely many      *)
(*               different orientations.                            *)
(*                                                                  *)
(*               The inflation procedure is performed together      *)
(*               with a rescaling of the length such that the       *)
(*               edge lengths of the tiles are always the same.     *)
(*                                                                  *)
(* :Bugs:        No bugs known so far.                              *)
(*               Please send bug reports to:                        *)
(*               u.grimm@physik.tu-chemnitz.de                      *)
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

BeginPackage["AperiodicTilings`PinwheelTiling`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

Inflation::usage =
"Inflation[tiling,n] performs an n-fold inflation of tiling which has
to consist of a list of single tiles given by the vertex coordinates."

PlotColorTiling::usage =
"PlotColorTiling[tiling,lincol,linewidth] produces a plot of tiling which 
has to consist of a list of single tiles given by their vertices. The
tiles are colored according to their orientation. The colour of edges 
can be changed with linecol (black by default), and their width by 
linewidth (1/200 by default)."

PlotTiling::usage =
"PlotTiling[tiling,lincol,linewidth] produces a plot of tiling which 
has to consist of a list of single tiles given by their vertices. The 
colour of edges can be changed with linecol (black by default), and 
their width by linewidth (1/200 by default)."

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
          PlotColorTiling,
          PlotTiling,
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

Clear[ColorTile,
      Inflation,
      PlotTiling,
      ScaleFactor,
      TileInflation,
      TileOrientation]

(********************************************)
(* Scaling factor used in an inflation step *)
(********************************************)

ScaleFactor = Sqrt[5]

(******************************)
(* Inflation of a single tile *)
(******************************)

TileInflation[{tilevertex1_,tilevertex2_,tilevertex3_}] :=
 ScaleFactor*{{#1,#5,#4},
              {#4,#7,#2},
              {#4,#7,#6},
              {#6,#5,#4},
              {#2,#6,#3}}&[tilevertex1,
                           tilevertex2,
                           tilevertex3,
                           tilevertex1/2+tilevertex2/2,
                           3*tilevertex1/5+2*tilevertex3/5,
	                   tilevertex1/5+4*tilevertex3/5,
	                   tilevertex1/10+tilevertex2/2+2*tilevertex3/5]

(*******************************)
(* Inflation of a set of tiles *)
(*******************************)

Inflation[tiling_List,
          num_Integer:1] /; num>=0 :=
  Nest[Flatten[Map[TileInflation,#],1]&,tiling,num]


(*************************)
(* Orientation of a tile *)
(*************************)

TileOrientation[{tilevertex1_, tilevertex2_, tilevertex3_}] := 
  N[Arg[-2*tilevertex1/3 + tilevertex2/3 + tilevertex3/3 /. 
    List -> Complex]/(2*Pi) + 1/2]; 

(**************************************)
(* Graphical presentation of a tiling *)
(**************************************)

PlotTiling[tiling_List,
	   linecol_:GrayLevel[0],
           linewidth_:1/200] :=
  Graphics[{linecol,Thickness[linewidth],
            Map[Line[Join[#,Take[#,2]]]&,tiling]},
	   AspectRatio->Automatic]

PlotColorTiling[tiling_List,
	        linecol_:GrayLevel[0],
                linewidth_:1/200] :=
  Graphics[{Map[{Hue[TileOrientation[#]],Polygon[#]}&,tiling],
            linecol,Thickness[linewidth],
            Map[Line[Join[#,Take[#,2]]]&,tiling]},
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
        PlotColorTiling,
        PlotTiling,
        TileInflation]

(********************************************************************)
(*                                                                  *)
(*          End of package "AperiodicTilings`ChairTiling`"          *)
(*                                                                  *)
(********************************************************************)

EndPackage[]

