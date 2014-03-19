(* ::Package:: *)

(********************************************************************)
(*                                                                  *)
(* :Title:       OctagonalTiling                                    *)
(*                                                                  *)
(* :Authors:     Uwe Grimm                                          *)
(*                                                                  *)
(* :Context:     AperiodicTilings`OctagonalTiling`                  *)
(*                                                                  *)
(* :Version:     1.01                                               *)
(*                                                                  *)
(* :Date:        January 12, 2005                                   *)
(*                                                                  *)
(* :Summary:     This package is part of a collection of            *)
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
(* :References:  R. Ammann, B. Gruenbaum, and G.C. Shephard,        *)
(*               Aperiodic tiles,                                   *)
(*               Discrete Comput. Geom. 8 (1992) 1-25               *)
(*                                                                  *)
(*               M. Duneau, R. Mosseri, and C. Oguey,               *)
(*               Approximants of quasiperiodic structures           *)
(*               generated by the inflation mapping,                *)
(*               J. Phys. A: Math. Gen. 22 (1989) 4549-4564         *)
(*                                                                  *)
(*               B. Gruenbaum and G.C. Shephard,                    *)
(*               Tilings and Patterns,                              *)
(*               W.H. Freeman, New York (1987),                     *)
(*               Figure 10.4.14 on page 556                         *)
(*                                                                  *)
(* :Description: This package deals with two different ways of      *)
(*               constructing the octagonal tiling, also known      *)
(*               as the Ammann-Beenker tiling. This tiling          *)
(*               possesses an eightfold rotational symmetry and     *)
(*               contains two different kinds of tiles, squares     *)
(*               and 45-degree rhombs.                              *) 
(*                                                                  *)
(*               The first method employed here is the inflation    *)
(*               procedure. The main routine OctagonalInflation     *)
(*               inflates the tiling, preserving the lengths of     *)
(*               edges in the tiling. In this part of this          *)
(*               the squares of the tiling are treated as two       *)
(*               triangles in order to simplify the inflation       *)
(*               rules. In the routine PlotOctagonalTiling          *)
(*               that produces a graphical version of the tiling    *)
(*               obtained from OctagonalInflation, one can choose   *)
(*               whether the triangles are shown or not. Other      *)
(*               optional choices include colors of edges and       *) 
(*               tiles, and edge and vertex decorations that        *)
(*               encode the matching rules for the octagonal        *)
(*               tiling. The inflation procedure can nicely be      *)
(*               seen in a plot made with OctagonalInflationPlot,   *)
(*               which displays a tiling and its inflation next     *)
(*               to each other. In addition, some initial patches   *)
(*               of the tiling are provided, such as SquarePatch    *)
(*               or OctagonalPatch, and CutOctagonalPatch can be    *)
(*               used to cut out a smaller patch from a tiling.     *)
(*               The connectivity of the tilings as a graph can     *)
(*               also be obtained in matrix form with the command   *)
(*               OctagonalAdjacencyMatrix.                          *) 
(*                                                                  *)
(*               The second method used in this package is the      *)
(*               standard cut-and-project method, in this case      *) 
(*               from the four-dimensional hypercubic lattice.      *)      
(*               The acceptance domain is the projection of the     *)
(*               four-dimensional hypercube, the projections to     *)
(*               parallel and perpendicular space can be viewed     *)
(*               with the commands PlotCubeParallelProjection and   *)
(*               PlotCubeOrthogonalProjection. Actually, in this    *)
(*               particular case the projections turn out to be     *)
(*               the same, and the acceptance domain is a regular   *)
(*               octagon. A patch of the tiling is produced by      *)
(*               the command OctagonalProjectionTiling, starting    *)
(*               from an initial point of the hypercubic lattice    *)
(*               and considering all those lattice points that can  *)  
(*               be reached in a given number of steps along the    *)
(*               basic lattice vectors. The resulting tiling in     *)
(*               parallel space and the distribution of points in   *)
(*               the acceptance domain in perpendicular space can   *)
(*               be viewed with PlotParallelProjection and with     *)
(*               PlotOrthogonalProjection, respectively.            *)
(*                                                                  *)
(* :Notes:       While most calculations involve only integer       *)
(*               numbers, some of the checks are done using         *)
(*               real numbers. In this case, in order to decide     *)
(*               whether a real number equals zero or not, we       *)
(*               compare with the internal variable                 *)
(*               NumericalAccuracy which by default we chose to     *)
(*               to be two orders of magnitude larger than the      *)
(*               $MachinePrecision. If any problems occurs due      *)
(*               numerical inaccuracies, the working precision      *)
(*               NumericalPrecision should be adjusted accordingly  *)
(*               using the command SetNumericalPrecision.           *)
(*                                                                  *)
(*               Due to the different construction principles, the  *)
(*               tilings produced by inflation and by projection    *)
(*               are stored in very different ways. Therefore, it   *)
(*               is not easily possible to transfer from one        *)
(*               description to the other. This is also the reason  *)
(*               why two separate routines for the graphical        *)
(*               representation of the tilings have been included   *)
(*               in this package.                                   *)
(*                                                                  *)
(* :Changes:     January 12, 2005 [UG] Version 1.0 -> 1.01:         *)
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
(* :Disclaimer:  No guarantee for correctness of the program and    *)
(*               results obtained with the program is given. Use    *)
(*               at your own risk!                                  *)
(*                                                                  *)
(********************************************************************)

BeginPackage["AperiodicTilings`OctagonalTiling`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

OctagonalInflation::usage = 
"OctagonalInflation[tiling,n] performs an n-fold inflation of tiling
which has be to a list of triangles and rhombs, which in turn are
specified by three and by four coordinates, respectively. Here, each
coordinate is given as a list of two integers {a,b}, corresponding to
the real number a+b*Sqrt[2], compare TrianglePatch, RhombPatch or
SquarePatch for some examples. In this way, all calculations can be
done with integers. The default value for the number of inflations n
is one."

CutOctagonalPatch::usage =
"CutOctagonalPatch[tiling,centerpoint,radius] cuts a spherical patch
of given center point and radius from an inflation tiling. The tiling
has to be specified as in OctagonalInflation."

OctagonalAdjacencyMatrix::usage =
"OctagonalAdjacencyMatrix[tiling] returns the adjacency matrix of
tilings encoding the neighboring vertices in the tiling. The matrix
has elements zero if two vertices are not neighbors, and one if they
are neighbors, and is symmetric. The tiling has to be specified as in
OctagonalInflation."
          
PlotOctagonalTiling::usage =
"PlotOctagonalTiling[tiling,linethickness,psize,showtriangles,
showtrianglelines,showedgedeco,showvertexdeco,pointcolor,linecolor,
rhombcolor,trianglecolor,edgedecocolor,vertexdecocolor] produces
graphics output for an inflation tiling that has to be specified as in
OctagonalInflation, use Show to display the graphics. The argument
linethickness determines the thicknesses of lines, the default values
is 1/100. The argument psize determines the sizes of points
representing the vertices, its default value is 0 (no points). The two
arguments showtriangles and showtrianglelines determine whether
triangles and their edges are shown in the picture, their default
values are False.  Similarly, showedgedeco and showvertexdeco concern
the edge and vertex decorations corresponding to the matching rules,
these are also False by default. The remaining optional arguments
pointcolor, linecolor, rhombcolor, trianglecolor, edgedecocolor, and
vertexdecocolor determine the colors of the various elements, they can
be specified either by a single number num, resulting in
GrayLevel[num], or by a list of three of four numbers, in which cases
the corresponding RGB or CMYK color is used, respectively."

OctagonalInflationPlot::usage =
"OctagonalInflationPlot[tiling,linethickness,psize,showtriangles,
showtrianglelines,showedgedeco,showvertexdeco,pointcolor,linecolor,
rhombcolor,trianglecolor,edgedecocolor,vertexdecocolor] uses
PlotOctagonalTiling to produce graphics for tiling and its inflated
version made by OctagonalInflation. The meaning and default values of
the arguments are the same as in PlotOctagonalTiling."

TrianglePatch::usage =
"TrianglePatch provides an initial patch for the inflation procedure
consisting of a single triangle (i.e., of half of a square)."

RhombPatch::usage =
"RhombPatch provides an initial patch for the inflation procedure
consisting of a single rhomb."

SquarePatch::usage =
"SquarePatch provides an initial patch for the inflation procedure
consisting of a single square (i.e., of two triangles)."

RotatedSquarePatch::usage = 
"RotatedSquarePatch provides an initial patch for the inflation
procedure consisting of a single square (i.e., of two triangles),
rotated by 45 degrees with respect to the square produced by
SquarePatch."

OctagonalPatch::usage = 
"OctagonalPatch provides an initial patch for the inflation procedure
consisting of eight squares (i.e., sixteen triangles) and sixteen
rhombs. It has the shape of a regular octagon. The patches obtained by
inflation of this particular patch will have perfect eightfold
rotational symmetry with respect to the center vertex of the patch."

SectorPatch::usage =
"SectorPatch provides an initial patch for the inflation procedure
consisting of two triangles and two rhombs, corresponding to one
sector of the eightfold patch OctagonalPatch."

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[OctagonalInflation,
          CutOctagonalPatch,
          OctagonalAdjacencyMatrix,
          PlotOctagonalTiling,
          OctagonalInflationPlot,
          TrianglePatch,
          RhombPatch,
          SquarePatch,
          RotatedSquarePatch,
          OctagonalPatch,
          SectorPatch]

(********************************************************************)
(*                                                                  *)
(*                     Start of Private context                     *)
(*                                                                  *)
(********************************************************************)

Begin["`Private`"]


Clear[SilverCoordinate,
      NumericalSilverCoordinate,
      SilverMultiply,
      SilverFactor1,
      SilverFactor2,
      OctagonalInflation,
      OctagonalTileInflation,
      CutOctagonalPatch,
      TrianglePatch,
      RhombPatch,
      SquarePatch,
      RotatedSquarePatch,
      OctagonalPatch,
      SectorPatch,
      OctagonalVertices,
      VertexToNumber,
      NumberToVertex,
      OctagonalLines,
      OctagonalLineLengths,
      OctagonalShortLines,
      OctagonalLongLines,
      OctagonalAdjacencyMatrix,
      OctagonalRhombs,
      OctagonalTriangles,
      OctagonalSquares,
      OctagonalEdgeDecoration,
      OctagonalVertexDecoration,
      PlotOctagonalTiling,
      OctagonalInflationPlot]

(**********************)
(**********************)
(** INFLATION METHOD **)
(**********************)
(**********************)

(******************************************)
(* Coordinates and their numerical values *)
(******************************************)

SilverCoordinate[{a_,b_}] := 
   a+b*Sqrt[2]

NumericalSilverCoordinate[{a_,b_}] := 
   N[a+b*Sqrt[2]]

NumericalSilverCoordinate[{{a_,b_},{c_,d_}}] := 
   N[{a+b*Sqrt[2],c+d*Sqrt[2]}]

(*********************************************)
(* Multiplication of point by scalar factors *)
(*********************************************)

SilverMultiply[{a_,b_},{c_,d_}] := 
   {a*c+2*b*d,a*d+b*c}

SilverMultiply[{{ax_,bx_},{ay_,by_}},{c_,d_}] := 
   {{ax*c+2*bx*d,ax*d+bx*c},
    {ay*c+2*by*d,ay*d+by*c}}

(*****************************************)
(* Scalar factors used in inflation rule *)
(*****************************************)

SilverFactor1 =
   {1,1}

SilverFactor2 = 
   {0,1/2}

(*******************************************)
(* Inflation of tiling (length-preserving) *)
(*******************************************)

OctagonalInflation[tiling_,n_Integer:1] :=
   Nest[Flatten[Map[OctagonalTileInflation,#],1]&,tiling,n]

(****************************************************)
(* Inflation rules for triangles (halfs of squares) *)
(****************************************************)

OctagonalTileInflation[{a_,b_,c_}] := 
   {{#3+#6,#1,#1+#4},
    {#2-#4,#1+#4,#3+#5},
    {#2,#2-#7,#3+#5,#2-#4},
    {#2-#7,#3,#3+#5},
    {#3,#3+#6,#1+#4,#3+#5}}&[SilverMultiply[a,SilverFactor1],
                             SilverMultiply[b,SilverFactor1],
                             SilverMultiply[c,SilverFactor1],
                             SilverMultiply[b-a,SilverFactor2],
                             SilverMultiply[b+a-2*c,SilverFactor2],
                             a-c,
                             b-c]

(******************************)
(* Inflation rules for rhombs *)
(******************************)

OctagonalTileInflation[{a_,b_,c_,d_}] := 
   {{#3,#3-#7,#3+#5,#3+#6},
    {#3+#6,#2,#3+#5},
    {#1+#7,#2,#1-#5},   
    {#1,#1-#6,#1-#5,#1+#7},
    {#1-#6,#4,#1-#5},
    {#3-#7,#4,#3+#5},
    {#4,#3+#5,#2,#1-#5}}&[SilverMultiply[a,SilverFactor1],
                          SilverMultiply[b,SilverFactor1],
                          SilverMultiply[c,SilverFactor1],
                          SilverMultiply[d,SilverFactor1],
                          a-c,
                          b-c,
                          c-d]

(***************************)
(* Cut out part of a patch *)
(***************************)

CutOctagonalPatch[tiling_,
                  centerpoint_,
                  radius_] := 
   Select[tiling,
          Max[Map[(#.#)&,
              Map[NumericalSilverCoordinate,#]-
                  Table[centerpoint,{Length[#]}]]]<radius&];


(************************)
(* Some initial patches *)
(************************)

TrianglePatch =
   {{{{0,0},{0,1}},{{0,2},{0,1}},{{0,1},{0,2}}}}

RhombPatch =
   {{{{0,0},{0,0}},{{0,1},{0,1}},{{2,1},{0,1}},{{2,0},{0,0}}}}

SquarePatch =
   {{{{-1,0},{ 1,0}},{{1,0},{-1,0}},{{ 1,0},{ 1,0}}},
    {{{-1,0},{ 1,0}},{{1,0},{-1,0}},{{-1,0},{-1,0}}}}

RotatedSquarePatch =
   {{{{0,-1},{0,0}},{{0,1},{0,0}},{{0,0},{0, 1}}},
    {{{0,-1},{0,0}},{{0,1},{0,0}},{{0,0},{0,-1}}}}

OctagonalPatch =
   {{{{ 0, 1},{ 0, 1}},{{ 2, 1},{ 2, 1}},{{ 2, 1},{ 0, 1}}},
    {{{ 0, 0},{ 0, 0}},{{ 2, 0},{ 0, 0}},{{ 2, 1},{ 0, 1}},{{ 0, 1},{ 0, 1}}},
    {{{ 2, 0},{ 0, 0}},{{ 2, 2},{ 0, 0}},{{ 2, 1},{ 0, 1}}},
    {{{ 2, 2},{ 0, 0}},{{ 2, 2},{ 2, 0}},{{ 2, 1},{ 2, 1}},{{ 2, 1},{ 0, 1}}}, 
    {{{ 0,-1},{ 0, 1}},{{-2,-1},{ 2, 1}},{{ 0,-1},{ 2, 1}}},
    {{{ 0, 0},{ 0, 0}},{{ 0, 0},{ 2, 0}},{{ 0,-1},{ 2, 1}},{{ 0,-1},{ 0, 1}}},
    {{{ 0, 0},{ 2, 0}},{{ 0, 0},{ 2, 2}},{{ 0,-1},{ 2, 1}}},
    {{{ 0, 0},{ 2, 2}},{{-2, 0},{ 2, 2}},{{-2,-1},{ 2, 1}},{{ 0,-1},{ 2, 1}}},
    {{{ 0, 1},{ 0, 1}},{{ 2, 1},{ 2, 1}},{{ 0, 1},{ 2, 1}}},
    {{{ 0, 0},{ 0, 0}},{{ 0, 0},{ 2, 0}},{{ 0, 1},{ 2, 1}},{{ 0, 1},{ 0, 1}}},
    {{{ 0, 0},{ 2, 0}},{{ 0, 0},{ 2, 2}},{{ 0, 1},{ 2, 1}}},
    {{{ 0, 0},{ 2, 2}},{{ 2, 0},{ 2, 2}},{{ 2, 1},{ 2, 1}},{{ 0, 1},{ 2, 1}}},
    {{{ 0, 0},{-2, 0}},{{ 0, 0},{-2,-2}},{{ 0, 1},{-2,-1}}},
    {{{ 0, 0},{ 0, 0}},{{ 0, 1},{ 0,-1}},{{ 0, 1},{-2,-1}},{{ 0, 0},{-2, 0}}},
    {{{ 0, 1},{ 0,-1}},{{ 2, 1},{-2,-1}},{{ 0, 1},{-2,-1}}},
    {{{ 2, 1},{-2,-1}},{{ 2, 0},{-2,-2}},{{ 0, 0},{-2,-2}},{{ 0, 1},{-2,-1}}},
    {{{ 2, 0},{ 0, 0}},{{ 2, 2},{ 0, 0}},{{ 2, 1},{ 0,-1}}},
    {{{ 0, 0},{ 0, 0}},{{ 0, 1},{ 0,-1}},{{ 2, 1},{ 0,-1}},{{ 2, 0},{ 0, 0}}},
    {{{ 0, 1},{ 0,-1}},{{ 2, 1},{-2,-1}},{{ 2, 1},{ 0,-1}}},
    {{{ 2, 1},{-2,-1}},{{ 2, 2},{-2, 0}},{{ 2, 2},{ 0, 0}},{{ 2, 1},{ 0,-1}}},
    {{{ 0,-1},{ 0, 1}},{{-2,-1},{ 2, 1}},{{-2,-1},{ 0, 1}}},
    {{{ 0, 0},{ 0, 0}},{{-2, 0},{ 0, 0}},{{-2,-1},{ 0, 1}},{{ 0,-1},{ 0, 1}}},
    {{{-2, 0},{ 0, 0}},{{-2,-2},{ 0, 0}},{{-2,-1},{ 0, 1}}},
    {{{-2,-2},{ 0, 0}},{{-2,-2},{ 2, 0}},{{-2,-1},{ 2, 1}},{{-2,-1},{ 0, 1}}},
    {{{ 0, 0},{-2, 0}},{{ 0, 0},{-2,-2}},{{ 0,-1},{-2,-1}}},
    {{{ 0, 0},{ 0, 0}},{{ 0,-1},{ 0,-1}},{{ 0,-1},{-2,-1}},{{ 0, 0},{-2, 0}}},
    {{{ 0,-1},{ 0,-1}},{{-2,-1},{-2,-1}},{{ 0,-1},{-2,-1}}},
    {{{-2,-1},{-2,-1}},{{-2, 0},{-2,-2}},{{ 0, 0},{-2,-2}},{{ 0,-1},{-2,-1}}},
    {{{ 0, 0},{ 0, 0}},{{ 0,-1},{ 0,-1}},{{-2,-1},{ 0,-1}},{{-2, 0},{ 0, 0}}},
    {{{-2, 0},{ 0, 0}},{{-2,-2},{ 0, 0}},{{-2,-1},{ 0,-1}}},
    {{{ 0,-1},{ 0,-1}},{{-2,-1},{-2,-1}},{{-2,-1},{ 0,-1}}},
    {{{-2,-1},{-2,-1}},{{-2,-1},{ 0,-1}},{{-2,-2},{ 0, 0}},{{-2,-2},{-2, 0}}}}

SectorPatch = 
   {{{{0,1},{0,1}},{{2,1},{2,1}},{{2,1},{0,1}}},
    {{{0,0},{0,0}},{{2,0},{0,0}},{{2,1},{0,1}},{{0,1},{0,1}}},
    {{{2,0},{0,0}},{{2,2},{0,0}},{{2,1},{0,1}}},
    {{{2,2},{0,0}},{{2,2},{2,0}},{{2,1},{2,1}},{{2,1},{0,1}}}}

(**************************************)
(* Ordered list of vertices of tiling *)
(**************************************)

OctagonalVertices[tiling_] := 
   OctagonalVertices[tiling] =
   Module[{lv},
          lv = Union[Flatten[tiling,1]];
          lv = Sort[lv,OrderedQ[{Reverse[NumericalSilverCoordinate[#1]],
                                 Reverse[NumericalSilverCoordinate[#2]]}]&];
          lv]

(*******************************************************************)
(* Rules encoding relations between vertices and their coordinates *)
(*******************************************************************)

VertexToNumber[tiling_] :=
   VertexToNumber[tiling] =
   Dispatch[Table[Rule[OctagonalVertices[tiling][[i]],i],
                  {i,Length[OctagonalVertices[tiling]]}]]

NumberToVertex[tiling_] :=
   NumberToVertex[tiling] =
   Dispatch[Table[Rule[i,OctagonalVertices[tiling][[i]]],
                  {i,Length[OctagonalVertices[tiling]]}]]

(*************************************************)
(* Lines of tiling, given by their two endpoints *)
(*************************************************)

OctagonalLines[tiling_] :=
   OctagonalLines[tiling] = 
   Module[{lp},
          lp = tiling /. VertexToNumber[tiling];
          lp = Map[Append[#,First[#]]&,lp];
          lp = Join[Map[Partition[#,2]&,lp],
                    Map[Partition[#,2]&,Map[RotateLeft,lp]]];
          lp = Flatten[lp,1];
          lp = Union[lp,Map[Reverse,lp]];
          lp = Select[lp,#[[2]]>#[[1]]&];
          lp]

(****************************)
(* Squared lengths of lines *)
(****************************)

OctagonalLineLengths[tiling_] :=
   OctagonalLineLengths[tiling] =
   Module[{nvr,ll,len,p1,p2,x1,y1,x2,y2},
          nvr = NumberToVertex[tiling];
          ll[{p1_,p2_}]   := len[(p1/.nvr)-(p2/.nvr)];
          len[{{x1_,y1_},
               {x2_,y2_}}]:= {x1^2+2*y1^2+x2^2+2*y2^2,
                              2*(x1*y1+x2*y2)};
          Map[ll,OctagonalLines[tiling]]];


(********************************************************)
(* Short lines in tiling (by default those of length 2) *)
(********************************************************)

OctagonalShortLines[tiling_,shortlength_:4] := 
   OctagonalShortLines[tiling,shortlength] =
   Module[{nvr,lencheck,p1,p2,len,x1,x2,y1,y2},
          nvr = NumberToVertex[tiling];
          lencheck[{p1_,p2_}] := (len[(p1/.nvr)-(p2/.nvr)]==={shortlength,0});
          len[{{x1_,y1_},
               {x2_,y2_}}]:= {x1^2+2*y1^2+x2^2+2*y2^2,
                              2*(x1*y1+x2*y2)};
          Select[OctagonalLines[tiling],lencheck]]

(***************************************************************)
(* Long lines in tiling (by default those of length 2*Sqrt[2]) *)
(***************************************************************)

OctagonalLongLines[tiling_,longlength_:8] := 
   OctagonalLongLines[tiling,longlength] =
   Module[{nvr,lencheck,p1,p2,len,x1,x2,y1,y2},
          nvr = NumberToVertex[tiling];
          lencheck[{p1_,p2_}] := (len[(p1/.nvr)-(p2/.nvr)]==={longlength,0});
          len[{{x1_,y1_},
               {x2_,y2_}}]:= {x1^2+2*y1^2+x2^2+2*y2^2,
                              2*(x1*y1+x2*y2)};
          Select[OctagonalLines[tiling],lencheck]]

(*******************************)
(* Adjacency matrix for tiling *)
(*******************************)

OctagonalAdjacencyMatrix[tiling_] :=
   Module[{sl,adj},
          adj = Table[{},{Length[OctagonalVertices[tiling]]}];
          sl  = OctagonalShortLines[tiling];
          Do[AppendTo[adj[[sl[[i,1]]]],sl[[i,2]]];
             AppendTo[adj[[sl[[i,2]]]],sl[[i,1]]],
             {i,Length[sl]}];
          adj]

(******************************************************************)
(* Tiles of tiling: rhombs, squares, triangles (halfs of squares) *)
(******************************************************************)

OctagonalRhombs[tiling_] := 
   Select[tiling/.VertexToNumber[tiling],
          Length[#]==4&]

OctagonalTriangles[tiling_] := 
   Select[tiling/.VertexToNumber[tiling],
          Length[#]==3&]

OctagonalSquares[tiling_] := 
   Module[{tt=OctagonalTriangles[tiling],
           ll=OctagonalLongLines[tiling],
           t,sq={},newsq,p1,p2,p3,p4},
          newsq[{p1_,p2_},{p3_,p4_}] := {p1,p3,p2,p4};
          Do[t=Select[tt,
                      Complement[ll[[i]],#]=={}&];
             If[Length[t]==2,
                AppendTo[sq,
                         newsq[Complement[Flatten[t],
                                          ll[[i]]],
                               ll[[i]]]]],
             {i,Length[ll]}];
          sq]

(*********************************************)
(* Decoration of tiling by arrows and houses *)
(*********************************************)

Clear[OctagonalEdgeDecoration];
OctagonalEdgeDecoration[{p1_List, 
                         p2_List, 
                         p3_List}] := 
   Map[Polygon,
       {{#3-#1,#3+#1-#2,#3+#1},
        {#4-#2,#4+#1+#2,#4+#2}}&[(p2-p3)/10,
                                 (p3-p1)/10,
				 (p2+p3)/2,
                                 (p3+p1)/2]]

OctagonalEdgeDecoration[{p1_List, 
                         p2_List, 
                         p3_List, 
                         p4_List}] := 
   Map[Polygon,
       {{#3+#1,#3+2*#1-#7*#2,#3-#1},
        {#4-#2,#4+#7*#1-2*#2,#4+#2},
        {#5-#1,#5-2*#1+#7*#2,#5+#1},
        {#6+#2,#6-#7*#1+2*#2,#6-#2}}&[(p1-p2)/10,
                                      (p2-p3)/10,
				      (p1+p2)/2,
                                      (p2+p3)/2,
                                      (p3+p4)/2,
                                      (p4+p1)/2,
			              N[Sqrt[2]]]]

OctagonalVertexDecoration[{p1_List,
                           p2_List, 
                           p3_List}] := 
   Map[Polygon,
       {{#1,#1-#4,#1+#6},
        {#2,#2+#4,#2-#5},
        {#3,#3-#6,#3+#5-#6,#3+#5}}&[p1, 
                                    p2,
                                    p3,
				    (p1-p2)/5,
                                    (p2-p3)/5,
                                    (p3-p1)/5]]

OctagonalVertexDecoration[{p1_List,
                           p2_List,
                           p3_List,
                           p4_List}] := 
   Map[Polygon,
       {{#1,#1-#6,#1-#7*#5},
        {#2,#2+#5,#2+2*#5-#7*#6,#2+#5-#7*#6,
         #2+#7*#5-#7^2*#6,#2-#6},
	{#3,#3+#5,#3+#7*#6},
        {#4,#4+#6,#4+2*#6-#7*#5,#4+#6-#7*#5,
         #4+#7*#6-#7^2*#5,#4-#5}}&[p1, 
                                   p2, 
                                   p3, 
                                   p4,
				   (p1-p2)/5,
                                   (p2-p3)/5,
			           N[Sqrt[2]]]]

(****************************************)
(* Plot of tiling obtained by inflation *)
(****************************************)

PlotOctagonalTiling[tiling_,
                    linethickness_:1/100,
                    psize_:0,
                    showtriangles_:False,
                    showtrianglelines_:False,
                    showedgedeco_:False,
                    showvertexdeco_:False,
                    pointcolor_:{0,0,1},
                    linecolor_:0,
                    rhombcolor_:{0,1,0},
                    trianglecolor_:{0,0,1,0},
                    edgedecocolor_:0,
                    vertexdecocolor_:0] :=
Module[{ntv=NumberToVertex[tiling],
        plotpoints,
        plotlines,
        plotsquares,
        plottriangle,
        plotdeco,
        color,
        cc},
       color[cc_Integer] := GrayLevel[cc];       
       color[cc_Real] := GrayLevel[cc];
       color[cc_List] := Which[Length[cc]==3,
                               RGBColor[Apply[Sequence,cc]],
                               Length[cc]==4,
                               CMYKColor[Apply[Sequence,cc]],
                               True,
                               Print["wrong color specification"];
                               Return[]];
       plotrhombs  = Join[{color[rhombcolor]},
                          Map[Polygon[Map[NumericalSilverCoordinate,#]]&,
                              OctagonalRhombs[tiling]/.ntv]];
       plottriangle = Join[{color[trianglecolor]},
                          Map[Polygon[Map[NumericalSilverCoordinate,#]]&,
                                          If[showtriangles,
                                             OctagonalTriangles[tiling],
                                             OctagonalSquares[tiling]]/.ntv]];
       plotpoints  = Join[{color[pointcolor],
                           PointSize[psize]},
                          Map[Point[NumericalSilverCoordinate[#]]&,
                              OctagonalVertices[tiling]]];
       plotlines  = If[linethickness==0,
                       {},
                       Join[{color[linecolor],
                             Thickness[linethickness]},
                            Map[Line[Map[NumericalSilverCoordinate,#]]&,
                                         If[showtrianglelines,
                                            OctagonalLines[tiling],
                                            OctagonalShortLines[tiling]]/.
                                            ntv]]];
       plotdeco   = Join[{color[edgedecocolor]},
                         If[showedgedeco,
                            Map[OctagonalEdgeDecoration,#],
                            {}],
                         {color[vertexdecocolor]},
		         If[showvertexdeco,
                            Map[OctagonalVertexDecoration,#],
                            {}]]&[Map[NumericalSilverCoordinate,
                                      tiling,{2}]];
       Graphics[Join[plottriangle,
                     plotrhombs,
                     plotpoints,
                     plotlines, 
                     plotdeco],
                AspectRatio -> Automatic,
                PlotRange -> All]]

(**********************************************)
(* Plot of tiling together with its inflation *)
(**********************************************)

Clear[OctagonalInflationPlot];
OctagonalInflationPlot[tiling_,
                       linethickness_:1/100,
                       psize_:0,
                       showtriangles_:False,
                       showtrianglelines_:False,
                       showedgedeco_:False,
                       showvertexdeco_:False,
                       pointcolor_:{0,0,1},
                       linecolor_:0,
                       rhombcolor_:{0,1,0},
                       trianglecolor_:{0,0,1,0},
                       edgedecocolor_:0,
                       vertexdecocolor_:0] :=
   GraphicsArray[{PlotOctagonalTiling[tiling,
                                      linethickness,
                                      psize,
                                      showtriangles,
                                      showtrianglelines,
                                      showedgedeco,
                                      showvertexdeco,
                                      pointcolor,
                                      linecolor,
                                      rhombcolor,
                                      trianglecolor,
                                      edgedecocolor,
                                      vertexdecocolor],
                  Graphics[{GrayLevel[0],
                            Line[{{0,1/2},{4/5,1/2}}],
                            Polygon[{{1,1/2},{4/5,2/5},{4/5,3/5}}]},
                            AspectRatio -> Automatic],
                  PlotOctagonalTiling[OctagonalInflation[tiling],
                                      linethickness,
                                      psize,
                                      showtriangles,
                                      showtrianglelines,
                                      showedgedeco,
                                      showvertexdeco,
                                      pointcolor,
                                      linecolor,
                                      rhombcolor,
                                      trianglecolor,
                                      edgedecocolor,
                                      vertexdecocolor]}]


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

Protect[OctagonalInflation,
        CutOctagonalPatch,
        OctagonalAdjacencyMatrix,
        PlotOctagonalTiling,
        OctagonalInflationPlot,
        TrianglePatch,
        RhombPatch,
        SquarePatch,
        RotatedSquarePatch,
        OctagonalPatch,
        SectorPatch]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`OctagonalTiling`"        *)
(*                                                                  *)
(********************************************************************)

EndPackage[]

