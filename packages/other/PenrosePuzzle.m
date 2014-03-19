(********************************************************************)
(*                                                                  *)
(* :Title:       PenrosePuzzle                                      *)
(*                                                                  *)
(* :Authors:     Uwe Grimm                                          *)
(*                                                                  *)
(* :Context:     AperiodicTilings`PenrosePuzzle`                    *)
(*                                                                  *)
(* :Version:     1.03                                               *)
(*                                                                  *)
(* :Date:        June 14, 2004                                      *)
(*                                                                  *)
(* :Summary:     Implementation of the construction of Penrose's    *)
(*               tiling of the plane by matching rules              *)
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
(*               physical properties and applications,              *)
(*               Springer, Berlin (2002)                            *)
(*                                                                  *)
(* :References:  R. Penrose,                                        *)
(*               The Role of Aesthetics in Pure and Applied         *)
(*               Mathematical Research,                             *)
(*               Bull. Inst. Math. Applications 10 (1974) 266;      *)
(*               Math. Intell. 2 (1979) 32                          *)
(*                                                                  *)
(*               N. G. de Bruijn,                                   *)
(*               Algebraic theory of Penrose's non-periodic         *)
(*               tilings of the plane. I & II,                      *)
(*               Indagationes mathematicae (Proc. Kon. Ned. Akad.   *)
(*               Wet. Ser. A) 84 (1991) 39-52, 53-66                *)
(*                                                                  *)
(* :Description: This package provides a procedure PenrosePuzzle    *)
(*               which allows you to assemble a rhombic Penrose     *)
(*               tiling, adding or removing a single tile at        *)
(*               each step. The matching rules of the Penrose       *)
(*               tiling are encoded be single and double arrows     *)
(*               on the edges of the two rhombic tiles, and after   *)
(*               each addition of a tile to a surface edge of       *)
(*               the patch several checks are performed to ensure   *)
(*               that the matching rules are not violated and       *)
(*               that tiles do not overlap. In case that any        *)
(*               inconsistency is deteced, the user may choose      *)
(*               whether to add the tile nonetheless or whether     *)
(*               to try another choice.                             *)
(*                                                                  *)
(* :Notes:       While most calculations involve only integer       *)
(*               numbers, some of the checks are done using         *)
(*               real numbers. In this case, in order to decide     *)
(*               whether a real number equals zero or not, we       *)
(*               compare with the internal variable                 *)
(*               NumericalAccuracy which by default we chose to     *)
(*               to be two orders of magnitude larger than the      *)
(*               $MachinePrecision. If any problems occurs due      *)
(*               numerical inaccuracies, NumericalAccuracy and      *)
(*               the working precision (NumericalPrecision) should  *)
(*               be adjusted accordingly using the commands         *)
(*               SetNumericalAccuracy and SetNumericalPrecision.    *)
(*                                                                  *)
(* :Changes:     June 14, 2004 [UG] Version 1.02 -> 1.03:           *)
(*               -> Fixed problem with "CollCheck" in               *)
(*                  Mathematica 5.0                                 *)
(*               July 23, 1999 [UG] Version 1.01 -> 1.02:           *)
(*               -> Still fixing the bug in "RemoveTile"            *)
(*               -> Changing input interpretation on error          *)
(*                  messages                                        *) 
(*               July 22, 1999 [UG]  Version 1.0 -> 1.01:           *)
(*               -> "RemoveTile" did not work properly for larger   *)
(*                  patches. This has (hopefully) been fixed in     *)
(*                  version 1.01.                                   *)
(*               -> "Thickness" in plot routines was replaced by    *)
(*                  "AbsoluteThickness" to keep line thickness      *)
(*                  constant.                                       *)
(*                                                                  *)
(* :Bugs:        No further bugs known so far.                      *)
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

BeginPackage["AperiodicTilings`PenrosePuzzle`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

AddTile::usage =
"AddTile[{tiling,surface},edge,tiletype] add a single tile of type
tiletype at the specified edge of the tiling, updating the set of
surface edges accordingly."

NumericalAccuracy::usage =
"NumericalAccuracy is the numerical accuracy used in numerical checks
on the correctness of the tiling. This means that all values that are
smaller (in modulus) than NumericalAccuracy are assumed to be exactly
zero. The value of NumericalAccuray can be changed by SetNumericalAccuray,
its default value is 10^(2-$MachinePrecision)."

NumericalPrecision::usage =
"NumericalPrecision is the numerical precision used to compute numerical
values for the coordinates of vertices. Its value can be changed by
SetNumericalPresicion. The default value of NumericalPrecision is
$MachinePrecision."

PenrosePuzzle::usage =
"PenrosePuzzle[initpatch,colfat,colthin,coledge,colsurfedge,colarrow,
colsurfarrow,collabel,linewidth,arrowsize,labelsize,plrange] starts
the puzzle program from the initial patch initpatch (a single fat
rhomb by default). The initial patch has to be given as a list
containing  two parts. The first part is a list of tiles, where each
tile {tiletype,refpoint,orient} is specified by its type ('F' for the
fat and 'T' for the thin rhomb), a reference point given as a
four-dimensional vector, and its orientation which is a number between
1 and 4. The second part gives a list of surface edges of the patch, each
given as a set of two integers, the first referring to the number of the
tile in the list tiling, the second (between 1 and 4) specifying the
particular edge of that tile. The optional parameters allow to change
the default settings for the colours of the fat (colfat) and thin
(colthin) rhombs, the edges (coledge) and the edges on the surface of
the patch (colsurfedge), the arrows (colarrows) and the arrows on the
surface edges (colsurfarrows), and finally the colour of the labels of
surface edges (collabel). In addition, the width of the edges (linewidth),
the size of the arrows (arrowsize) and the labels (labelsize), and
the range of the plot (plrange) can be modified from their default
values. Upon executing PenrosePuzzle[], a plot of the two basic
tiles and of the initial patch appears, and at each step the programs
inquires for the action to be taken (adding or removing tiles or
exiting the program). As output, the final patch of the puzzle
process is returned."

PlotTiling::usage =
"PlotTiling[{tiling,surface},colfat,colthin,coledge,colsurfedge,
colarrow,colsurfarrow,collabel,linewidth,arrowsize,labelsize,plrange]
produces a plot of the list of tiles contained in tiling, where each
tile {tiletype,refpoint,orient} is specified by its type ('F' for the
fat and 'T' for the thin rhomb), a reference point given as a
four-dimensional vector, and its orientation which is a number between
1 and 4. The second part of the first argument gives a list of surface
edges of the patch, each given as a set of two integers, the first
referring to the number of the tile in the list tiling, the second
(between 1 and 4) specifying the particular edge of that tile. The
optional parameters allow to change the default settings for the colours
of the fat (colfat) and thin (colthin) rhombs, the edges (coledge) and
the edges on the surface of the patch (colsurfedge), the arrows
(colarrows) and the arrows on the surface edges (colsurfarrows),
and finally the colour of the labels of surface edges (collabel).
In addition, the width of the edges (linewidth), the size of the
arrows (arrowsize) and the labels (labelsize), and the range of the
plot (plrange) can be modified from their default values."

RemoveTile::usage =
"RemoveTile[[{tiling,surface},edge] removes a tile from the tiling,
specified by its surface edge, and updates the set of surface edges
accordingly."

SetNumericalAccuracy::usage =
"SetNumericalAccuracy[numacc] replaces the value of NumericalAccuray,
used in numerical checks on the correctness of the tiling, by numacc.
This means that all values that are smaller (in modulus) than numacc
are assumed to be exactly zero. The default value of NumericalAccuray
is 10^(2-$MachinePrecision)."

SetNumericalPrecision::usage =
"SetNumericalPrecision[numprec] replaces the value of NumericalPrecision,
used to compute numerical values for the coordinates of vertices, by
the integer numacc. The default value of NumericalPrecision is
$MachinePrecision."

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[AddTile,
          NumericalAccuracy,
          NumericalPrecision,
          PenrosePuzzle,
          PlotTiling,
          RemoveTile,
          SetNumericalAccuracy,
          SetNumericalPrecision]

(********************************************************************)
(*                                                                  *)
(*                     Start of Private context                     *)
(*                                                                  *)
(********************************************************************)

Begin["`Private`"]

(*********************************)
(* Clearing previous definitions *)
(*********************************)

Clear[AddTile,
      ArrowCheck,
      BaseVector,
      CollCheck,
      CommonEdge,
      EdgeOrientation,
      EdgeType,
      FatRhomb,
      FatRhombSurface,
      IntersectCheck,
      NumericalAccuracy,
      NumericalPrecision,
      OrthogonalTwoVector,
      OverlapCheck,
      ParallelProjection,
      PenrosePuzzle,
      PlotEdge,
      PlotLine,
      PlotTile,
      PlotTiling,
      RemoveTile,
      RhombPlot,
      RotateDirection,
      RotationAngle,
      SetNumericalAccuracy,
      SetNumericalPrecision,
      ThinRhomb,
      ThinRhombSurface,
      TileAngle,
      TileVertex,
      TwoCoordinates,
      TwoVector]

(*************************************************)
(* Numerical Accuracy used for checking purposes *)
(*************************************************)

NumericalPrecision = $MachinePrecision

NumericalAccuracy = 10.^(2 - $MachinePrecision)

SetNumericalPrecision[numprec_Integer] :=
  Module[{},
         Unprotect[NumericalPrecision];
         NumericalPrecision = numprec;
         Protect[NumericalPrecision];]

SetNumericalAccuracy[numacc_] :=
  Module[{},
         Unprotect[NumericalAccuracy];
         NumericalAccuracy = numprec;
         Protect[NumericalAccuracy];]

(***********************************************)
(* Definitions of base vectors and coordinates *)
(***********************************************)

BaseVector[1] = {1,0,0,0}

BaseVector[2] = {0,1,0,0}

BaseVector[3] = {0,0,1,0}

BaseVector[4] = {0,0,0,1}

BaseVector[5] = -Apply[Plus,Array[BaseVector,4]]

BaseVector[j_]/;(-5<=j<=-1) :=
  BaseVector[j] =
  -BaseVector[-j]

TwoVector[j_] :=
  TwoVector[j] =
  N[{Cos[2*Pi*(j-1)/5],Sin[2*Pi*(j-1)/5]},
    NumericalPrecision]

ParallelProjection :=
  Array[TwoVector,4]

TwoCoordinates[pt_List] :=
  Dot[pt,ParallelProjection]

OrthogonalTwoVector[{x_,y_}] :=
  {-y,x}

(****************************************)
(* Rotations by multiples of 36 degrees *)
(****************************************)

RotateDirection[j_] /; (1<=j<=5) :=
  RotateDirection[j] =
  RotateRight[{-1,-2,-3,-4,-5},3-j][[1]]

RotateDirection[j_] /; (-5<=j<=-1) :=
  RotateDirection[j] =
  -RotateDirection[-j]

RotateDirection[j_,n_] :=
  RotateDirection[j,n] =
  Nest[RotateDirection,j,n]

(**********************************************)
(* Rotation angles as multiples of 36 degrees *)
(**********************************************)

RotationAngle[j_] /;  1<=j<=5  :=
  RotationAngle[j] =
  2*(j-1)

RotationAngle[j_] /; -5<=j<=-1 :=
  RotationAngle[j] =
  Mod[RotationAngle[-j]+5,10]

(**********************************)
(* Single tiles and their surface *)
(**********************************)

FatRhomb         = {{"F",{0,0,0,0},1}}

FatRhombSurface  = {{1,1},{1,2},{1,3},{1,4}}

ThinRhomb        = {{"T",{0,0,0,0},1}}

ThinRhombSurface = {{1,1},{1,2},{1,3},{1,4}}

(*******************************************)
(* Angle of tiles (in units of 36 degrees) *)
(*******************************************)

TileAngle["T"] = 1

TileAngle["F"] = 2

(*************************************)
(* Directions and type of tile edges *)
(*************************************)

EdgeType[{"F",pt_,dir_},1] :=
  {pt,dir,1}

EdgeType[{"F",pt_,dir_},2] :=
  {pt+BaseVector[dir],RotateDirection[dir,2],2}

EdgeType[{"F",pt_,dir_},3] :=
  {pt+BaseVector[RotateDirection[dir,2]],dir,2}

EdgeType[{"F",pt_,dir_},4] :=
  {pt,RotateDirection[dir,2],1}

EdgeType[{"T",pt_,dir_},1] :=
  {pt,dir,1}

EdgeType[{"T",pt_,dir_},2] :=
  {pt+BaseVector[dir]+BaseVector[RotateDirection[dir,1]],
   -RotateDirection[dir,1],1}

EdgeType[{"T",pt_,dir_},3] :=
  {pt+BaseVector[dir]+BaseVector[RotateDirection[dir,1]],-dir,2}

EdgeType[{"T",pt_,dir_},4] :=
  {pt,RotateDirection[dir,1],2}

(*****************************)
(* Orientation of tile edges *)
(*****************************)

EdgeOrientation[{"F",_,_},1] = -1

EdgeOrientation[{"F",_,_},2] = -1

EdgeOrientation[{"F",_,_},3] =  1

EdgeOrientation[{"F",_,_},4] =  1

EdgeOrientation[{"T",_,_},1] = -1

EdgeOrientation[{"T",_,_},2] =  1

EdgeOrientation[{"T",_,_},3] = -1

EdgeOrientation[{"T",_,_},4] =  1

(*************************************)
(* Position of the vertices of tiles *)
(*************************************)

TileVertex[{tiletype_,pt_,dir_},1] :=
  pt

TileVertex[{tiletype_,pt_,dir_},2] :=
  pt+BaseVector[dir]

TileVertex[{tiletype_,pt_,dir_},3] :=
  pt+BaseVector[dir]+BaseVector[RotateDirection[dir,TileAngle[tiletype]]]

TileVertex[{tiletype_,pt_,dir_},4] :=
  pt+BaseVector[RotateDirection[dir,TileAngle[tiletype]]]

(**************************************************)
(* Possible common edges shared between two tiles *)
(**************************************************)

CommonEdge["F","F"] = {4,3,2,1}

CommonEdge["T","T"] = {2,1,4,3}

CommonEdge["F","T"] = {2,4,3,1}

CommonEdge["T","F"] = {4,1,3,2}

(**************************************)
(* Lines with single or double arrows *)
(**************************************)

PlotLine[{pt_,dir_,numarrows_Integer},
         numlabel_Integer:0,
         coledge_:GrayLevel[0],
         colarrow_:GrayLevel[0],
         collabel_:GrayLevel[0],
         linewidth_:3/2,
         arrowsize_:1/10,
         orientation_:1] :=
  Module[{cp1=TwoCoordinates[pt],
          cp2=TwoCoordinates[pt+BaseVector[dir]]},
         Join[If[numlabel!=0,
                 {collabel,
                  Text[ToString[numlabel],
                       #1-#2/(10*arrowsize)+
                       OrthogonalTwoVector[2*orientation*#2]]},
                 {}],
              {AbsoluteThickness[linewidth],
               coledge,
               Line[{cp1,cp2}],
               colarrow},
              If[numarrows==1,
                 {Line[{#1-#2+OrthogonalTwoVector[#2],
                        #1+#2,
                        #1-#2-OrthogonalTwoVector[#2]}]},
                 {Line[{#1-3*#2/2+OrthogonalTwoVector[#2],
                        #1+#2/2,
                        #1-3*#2/2-OrthogonalTwoVector[#2]}],
                  Line[{#1-#2/2+OrthogonalTwoVector[#2],
                        #1+3*#2/2,
                        #1-#2/2-OrthogonalTwoVector[#2]}]}]
             ]&[(cp2+cp1)/2,arrowsize*(cp2-cp1)]
        ]

(***********************************************)
(* Edges of a given tile with arrow decoration *)
(***********************************************)

PlotEdge[tile_,
         linenum_,
         numlabel_Integer:0,
         coledge_:GrayLevel[0],
         colarrow_:GrayLevel[0],
         collabel_:GrayLevel[0],
         linewidth_:3/2,
         arrowsize_:1/10] :=
  PlotLine[EdgeType[tile,linenum],
           numlabel,
           coledge,
           colarrow,
           collabel,
           linewidth,
           arrowsize,
           EdgeOrientation[tile,linenum]]

(********************************************************)
(* Plot of a single tile with arrow decoration of edges *)
(********************************************************)

PlotTile[tile_,
         colfat_:GrayLevel[1],
         colthin_:GrayLevel[1],
         coledge_:GrayLevel[0],
         colarrow_:GrayLevel[0],
         collabel_:GrayLevel[0],
         linewidth_:3/2,
         arrowsize_:1/10] :=
  Join[{If[First[tile]==="F",
           colfat,
           colthin],
        Polygon[
          Table[
            TwoCoordinates[
              TileVertex[tile,i]],
               {i,4}]
               ],
        coledge},
        Table[
          PlotEdge[tile,
                   i,
                   0,
                   coledge,
                   colarrow,
                   collabel,
                   linewidth,
                   arrowsize],
              {i,4}]]

(********************)
(* Plot of a tiling *)
(********************)

PlotTiling[{tiling_List,surface_List},
           colfat_:GrayLevel[1],
           colthin_:GrayLevel[1],
           coledge_:GrayLevel[0],
           colsurfedge_:GrayLevel[0],
           colarrow_:GrayLevel[0],
           colsurfarrow_:GrayLevel[0],
           collabel_:GrayLevel[0],
           linewidth_:3/2,
           arrowsize_:1/10,
           labelsize_:18,
           plrange_:All] :=
  Graphics[
    Join[
      Map[
        PlotTile[#,
                 colfat,
                 colthin,
                 coledge,
                 colarrow,
                 collabel,
                 linewidth,
                 arrowsize]&,
        tiling],
      Map[
        PlotEdge[tiling[[#[[1]]]],
                 #[[2]],
                 #[[3]],
                 colsurfedge,
                 colsurfarrow,
                 collabel,
                 linewidth,
                 arrowsize]&,
        MapIndexed[
          Join[#1,#2]&,
          surface]]
        ],
    AspectRatio -> Automatic,
    DefaultFont -> {"Times-Roman",labelsize},
    PlotRange   -> plrange]

(*********************************)
(* Plot of the two rhombic tiles *)
(*********************************)

RhombPlot[colfat_:GrayLevel[1],
          colthin_:GrayLevel[1],
          coledge_:GrayLevel[0],
          colarrow_:GrayLevel[0],
          linewidth_:3/2,
          arrowsize_:1/10] :=
  GraphicsArray[
    MapThread[
      Graphics[
        Join[
          PlotTile[{#1,{0,0,0,0},1},
                   colfat,
                   colthin,
                   coledge,
                   colarrow,
                   GrayLevel[0],
                   linewidth,
                   arrowsize],
          {GrayLevel[0],
           Text[FontForm[#1,{"Helvetica",24}],
                TwoCoordinates[#2]]}
            ],
      AspectRatio -> Automatic,
      PlotRange   -> #3]&,
    {{"F","T"},
     {{0.5,0.5,0,0},
      {0.5,0,0,-0.5}},
     {{{-0.5,1.7},
       {-0.2,1.2}},
      {{-0.2,2.},
       {-0.2,1.2}}}}]]

(**********************************************)
(* Check that arrows on identical edges match *)
(**********************************************)

ArrowCheck[tiling_] :=
  Module[{alledges,
          dbledges},
         alledges = Map[Table[EdgeType[#,i],{i,4}]&,
                        tiling];
         alledges = Map[Drop[#,-1]&,
                        alledges,
                        {2}];
         dbledges = Union[Flatten[alledges,1]];
         dbledges = Join[dbledges,
                         Map[{#[[1]]+BaseVector[#[[2]]],-#[[2]]}&,
                             dbledges]];
         dbledges = Union[Select[dbledges,
                                 Count[dbledges,#]>1&]];
         Flatten[Map[Position[alledges,#]&,
                     dbledges],
                 1]]

(*******************************************)
(* Check whether two vectors are collinear *)
(*******************************************)

CollCheck[vec1_,vec2_] :=
(NumericalAccuracy-Abs[Dot[vec1,OrthogonalTwoVector[vec2]]]>0)

(*********************************************)
(* Check whether two line segments intersect *)
(*********************************************)

IntersectCheck[{{pt1_,pt2_},{pt3_,pt4_}}] :=
  If[CollCheck[pt2-pt1,pt4-pt3],
     False,
     Or[And[-NumericalAccuracy<#[[1]],
            NumericalAccuracy<#[[2]]<1-NumericalAccuracy],
        And[NumericalAccuracy<#[[1]]<1-NumericalAccuracy,
            #[[2]]<1+NumericalAccuracy]
       ]&[Sort[{Dot[#1,pt3-pt1]/Dot[#,pt2-pt1]&[
                  OrthogonalTwoVector[pt4-pt3]],
                -Dot[#1,pt3-pt1]/Dot[#,pt4-pt3]&[
                  OrthogonalTwoVector[pt2-pt1]]}]]]

(***********************************)
(* Check whether two tiles overlap *)
(***********************************)

OverlapCheck[tile1_,
             tile2_] :=
  If[tile1 == tile2,
     True,
     MemberQ[
       Map[IntersectCheck,
           Flatten[
             Outer[List,
                   Partition[Append[#,First[#]],
                             2,1]&[Array[
                                     TwoCoordinates[
                                       TileVertex[tile1,#]]&,4]],
                   Partition[Append[#,First[#]],
                             2,1]&[Array[
                                     TwoCoordinates[
                                       TileVertex[tile2,#]]&,4]],
                   1],
                   1]],
             True]]

(*********************************************)
(* Remove one tile from a patch [which must  *)
(* have an edge on the surface of the patch] *)
(*********************************************)

RemoveTile[{til_List,surf_List},
           edge_] :=
Module[{alledges,
        checknew,
        i,
        j,
        misfit,
        newsurf,
        nnew,
        remtile,
        surface=surf,
        tilenum,
        tiling=til},
       tilenum  = surface[[edge,1]];
       misfit   = ArrowCheck[tiling];
       If[misfit=!={},
          surface = Union[surface,misfit]];
       remtile  = tiling[[tilenum]];
       tiling   = Delete[tiling,tilenum];
       newsurf  = Map[#[[2]]&,Select[surface,#[[1]]==tilenum&]];
       newsurf  = Complement[{1,2,3,4},newsurf];
       newsurf  = Map[EdgeType[remtile,#]&,newsurf];
       nnew     = Length[newsurf];
       alledges = Map[Table[EdgeType[#,i],{i,4}]&,tiling];
       newsurf  = Flatten[Map[Position[alledges,#]&,newsurf],1];
       checknew = Map[EdgeType[tiling[[#[[1]]]],#[[2]]]&,newsurf];
       checknew = Map[Count[checknew,#]&,checknew];
       If[Union[checknew]=!={1},
          checknew = Complement[Partition[Range[Length[checknew]],1],
                                Position[checknew,1]];
          newsurf  = Delete[newsurf,checknew]];
       surface  = Select[surface,#[[1]]!=tilenum&];
       surface  = surface /. Table[{i+1,j_}->{i,j},
                                   {i,tilenum,Length[tiling]}];
       surface  = Sort[Join[surface,newsurf]];
       If[misfit=!={},
          misfit = ArrowCheck[tiling]];
       If[misfit=!={},
          Print["Error warning: The tiling still contains arrow misfits!"];
          surface = Complement[surface,misfit]];
       {tiling,surface}]

(***************************)
(* Add one tile to a patch *)
(***************************)

AddTile[{til_List,surf_List},
        edge_,
        tiletype_] :=
Module[{angle,
        checkdouble,
        checkoverlap,
        checksurf,
        checktile,
        i,
        misfit,
        newsurface,
        newtile,
        newtilenum,
        surface=surf,
        surfedge,
        surfedgetype,
        surftile,
        tileedge,
        tileedgetype,
        tiling=til},
       surfedge     = surface[[edge]];
       surftile     = tiling[[surfedge[[1]]]];
       surfedgetype = EdgeType[surftile,surfedge[[2]]];
       tileedge     = CommonEdge[surftile[[1]],
                                 tiletype][[surfedge[[2]]]];
       tileedgetype = EdgeType[{tiletype,{0,0,0,0},1},tileedge];
       angle        = Mod[RotationAngle[surfedgetype[[2]]]-
                          RotationAngle[tileedgetype[[2]]]+10,10];
       newtile      = {tiletype,{0,0,0,0},RotateDirection[1,angle]};
       tileedgetype = EdgeType[newtile,tileedge];
       newtile[[2]] = surfedgetype[[1]]-tileedgetype[[1]];
       newtilenum   = Length[tiling]+1;
       surface      = Delete[surface,edge];
       newsurface   = Delete[Table[{newtilenum,i},{i,4}],tileedge];
       checksurf    = Map[EdgeType[tiling[[#[[1]]]],#[[2]]]&,
                          surface];
       checktile    = Map[EdgeType[newtile,#[[2]]]&,newsurface];
       checkdouble  = Intersection[checksurf,checktile];
       If[checkdouble=!={},
          surface    = Delete[surface,
                              Partition[Flatten[
                              Map[Position[checksurf,#]&,
                                           checkdouble]],1]];
          newsurface = Delete[newsurface,
                              Partition[Flatten[
                              Map[Position[checktile,#]&,
                                           checkdouble]],1]]];
       checksurf    = Flatten[Map[Table[EdgeType[#,i],{i,4}]&,
                                  tiling],1];
       checktile    = Map[EdgeType[newtile,#[[2]]]&,newsurface];
       checkdouble  = Intersection[checksurf,checktile];
       If[checkdouble=!={},
          If[ToString[
          Input["Error detected: "<>
                "New tile shares an edge with two other tiles!\n"<>
                "Do you want to continue adding the tile "<>
                "(answer YES or NO)?"]] =!=
                "YES",
                Return["NO"]];
          newsurface = Delete[newsurface,
                              Partition[Flatten[
                              Map[Position[checktile,#]&,
                                           checkdouble]],1]]];
       If[Union[Map[OverlapCheck[newtile,#]&,tiling]]=!={False},
          If[ToString[
             Input["Error detected: "<>
                   "New tile overlaps with another tiles or contains\n"<>
                   "                a vertex that lies on an edge "<>
                   "of another tile!\n"<>
                   "Do you want to continue adding the tile "<>
                   "(answer YES or NO)?"]] =!=
                   "YES",
             Return["NO"]]]; 
       tiling  = Append[tiling,newtile];
       surface = Join[surface,newsurface];
       misfit  = ArrowCheck[tiling];
       If[misfit=!={},
          If[ToString[
             Input["Error detected: "<>
                   "Arrows do not fit properly!\n"<>
                   "Do you want to continue adding the tile "<>
                   "(answer YES or NO)?"]] =!=
                   "YES",
             Return["NO"]];
          surface = Complement[surface,misfit]];
       {tiling,surface}]

(******************)
(* Puzzle Program *)
(******************)

PenrosePuzzle[initpatch_:"F",
              colfat_:RGBColor[1,1,0],
              colthin_:RGBColor[0,1,0],
              coledge_:RGBColor[0,0,1],
              colsurfedge_:RGBColor[1,0,0],
              colarrow_:RGBColor[0,0,1],
              colsurfarrow_:RGBColor[1,0,0],
              collabel_:RGBColor[1,0,0],
              linewidth_:3/2,
              arrowsize_:1/10,
              labelsize_:18,
              plrange_:All] :=
Module[{edge,
        inp,
        inplen,
        misfit={},
        newplot,
        newtiling,
        presenttiling,
        tiletype},
       Off[General::spell1];
       Clear[F,T,BYE,YES,NO];
       Show[RhombPlot[colfat,
                      colthin,
                      coledge,
                      colarrow,
                      linewidth,
                      arrowsize]];
       Which[initpatch==="F",
             presenttiling = {FatRhomb,FatRhombSurface},
             initpatch==="T",
             presenttiling = {ThinRhomb,ThinRhombSurface},
             True,
             presenttiling = initpatch;
             misfit = ArrowCheck[presenttiling[[1]]]];
       If[misfit=!={},
          If[ToString[Input["Error detected: "<>
                            "Arrows do not fit properly in initial patch!\n"<>
                   "Do you want to continue (answer YES or NO)?"]]=!="YES",
             Return[]]];
       Show[PlotTiling[presenttiling,
                       colfat,
                       colthin,
                       coledge,
                       colsurfedge,
                       colarrow,
                       colsurfarrow,
                       collabel,
                       linewidth,
                       arrowsize,
                       labelsize,
                       plrange]];
       While[inp=!="BYE",
             newplot = True;
             inp = ToString[Input["Enter action:\n"<>
                         "AF2 adds a fat tile to surface edge 2\n"<>
                         "R1  removes the tile with surface edge 1\n"<>
                         "BYE terminates the program"]];
             inplen = StringLength[inp];
             Which[StringTake[inp,1]==="A",
                   tiletype = StringTake[inp,{2}];
                   edge = ToExpression[StringTake[inp,2-inplen]];
                   If[And[MemberQ[{"F","T"},tiletype],
                          MemberQ[Range[Length[presenttiling[[2]]]],edge]],
                      newtiling = AddTile[presenttiling,edge,tiletype];
                      If[newtiling==="NO",
                         newplot=False,
                         presenttiling = newtiling],
                      Print["Wrong input: Wrong tile type or non-existing"<>
                            " edge specified"];
                      newplot = False],
                   StringTake[inp,1]==="R",
                   edge = ToExpression[StringTake[inp,1-inplen]];
                   If[MemberQ[Range[Length[presenttiling[[2]]]],edge],
                      presenttiling = RemoveTile[presenttiling,edge],
                      Print["Wrong input: No such edge"],
                      newplot = False],
                   inp==="BYE",
                   Print["Good bye!"];
                   On[General::spell1];
                   Return[presenttiling],
                   True,
                   Print["Wrong input!"];
                   newplot = False];
           If[presenttiling==={{},{}},
              Return["Empty Tiling!"]];
           If[newplot,
              Show[PlotTiling[presenttiling,
                              colfat,
                              colthin,
                              coledge,
                              colsurfedge,
                              colarrow,
                              colsurfarrow,
                              collabel,
                              linewidth,
                              arrowsize,
                              labelsize,
                              plrange]]]]]

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

Protect[AddTile,
        NumericalAccuracy,
        NumericalPrecision,
        PenrosePuzzle,
        PlotTiling,
        RemoveTile,
        SetNumericalAccuracy,
        SetNumericalPrecision]

(********************************************************************)
(*                                                                  *)
(*         End of package "AperiodicTilings`PenrosePuzzle`"         *)
(*                                                                  *)
(********************************************************************)

EndPackage[]

