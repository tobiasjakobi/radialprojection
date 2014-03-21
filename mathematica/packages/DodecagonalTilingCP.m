(* ::Package:: *)

BeginPackage["AperiodicTilings`DodecagonalTilingCP`",{"EuclideanAlgorithmZSqrt3`"}]

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[DodecagonalWindowTest,
          DodecagonalProjectionTilingVis,
          DodecagonalSelectVertices,
          PlotRadialProjectionDistribution,
          PlotParallelProjection,
          ParallelProjectionRaw,
          ParallelProjectionAnnotate,
          ParallelProjectionDynamic,
          PlotOrthogonalProjection,
          OrthogonalProjectionRaw,
          NumericalAccuracy,
          NumericalPrecision,
          CheckProjectionInWindow,
          CheckScaledProjInWindow,
          WindowCenterShift,
          ShiftWindowCenter,
          SetNumericalPrecision]

(********************************************************************)
(*                                                                  *)
(*                     Start of Private context                     *)
(*                                                                  *)
(********************************************************************)

Begin["`Private`"]


Clear[SetNumericalPrecision,
      NumericalPrecision,
      NumericalAccuracy,
      DodecagonalLambda,
      InnerDodecagonRadius,
      OuterDodecagonRadius,
      DodecagonVertices,
      ProjectionVector,
      LatticeToPhysicalMatrix,
      WindowCenterShift,
      ShiftWindowCenter,
      VectorParallelProjection,
      VectorOrthogonalProjection,
      VectorOrthogonalProjectionShift,
      PrintNumericalPrecisionWarning,
      CheckPosition,
      CheckProjectionInWindow,
      CheckScaledProjInWindow,
      CheckProjectionInSector,
      DodecagonalWindowTest,
      DodecagonalProjectionTilingVis,
      DodecagonalSelectVertices,
      SectorSelectorDodecagon,
      PlotRadialProjectionDistribution,
      PlotParallelProjection,
      ParallelProjectionRaw,
      CreateVertexAnnotation,
      ParallelProjectionAnnotate,
      ParallelProjectionDynamic,
      PlotOrthogonalProjection,
      OrthogonalProjectionRaw]


(*************************************************)
(* Numerical Accuracy used for checking purposes *)
(*************************************************)

NumericalPrecision = $MachinePrecision

NumericalAccuracy = 10.^(2 - $MachinePrecision)

SetNumericalPrecision[numprec_Integer] :=
  Module[{},
         Unprotect[NumericalPrecision,NumericalAccuracy];
         NumericalPrecision = numprec;
         NumericalAccuracy = 10.^(2-numprec);
         Protect[NumericalPrecision,NumericalAccuracy];]

(*******************************************)
(* The fundamental unit in the n=12 case.  *)
(*******************************************)

DodecagonalLambda = 
   2 + Sqrt[3]

(***********************************************)
(* Inner and outer radius of acceptance domain *)
(***********************************************)

InnerDodecagonRadius = DodecagonalLambda/2

OuterDodecagonRadius = Sqrt[DodecagonalLambda]

(* Regular dodecagon with edge length = 1. *)
DodecagonVerticesAlt =
   {{OuterDodecagonRadius, 0},
    {Sqrt[3]/2*OuterDodecagonRadius, OuterDodecagonRadius/2},
    {OuterDodecagonRadius/2, Sqrt[3]/2*OuterDodecagonRadius},
    {0, OuterDodecagonRadius},
    {-OuterDodecagonRadius/2, Sqrt[3]/2*OuterDodecagonRadius},
    {-Sqrt[3]/2*OuterDodecagonRadius, OuterDodecagonRadius/2},
    {-OuterDodecagonRadius, 0},
    {-Sqrt[3]/2*OuterDodecagonRadius, -OuterDodecagonRadius/2},
    {-OuterDodecagonRadius/2, -Sqrt[3]/2*OuterDodecagonRadius},
    {0, -OuterDodecagonRadius},
    {OuterDodecagonRadius/2, -Sqrt[3]/2*OuterDodecagonRadius},
    {Sqrt[3]/2*OuterDodecagonRadius, -OuterDodecagonRadius/2}
}

normTwoElement = 1 + Sqrt[3]

(* Dodecagon with orientation from the book *)
DodecagonVerticesBook =
   {{DodecagonalLambda/2, 1/2},
    {normTwoElement/2, normTwoElement/2},
    {1/2, DodecagonalLambda/2},
    {-1/2, DodecagonalLambda/2},
    {-normTwoElement/2, normTwoElement/2},
    {-DodecagonalLambda/2, 1/2},
    {-DodecagonalLambda/2, -1/2},
    {-normTwoElement/2, -normTwoElement/2},
    {-1/2, -DodecagonalLambda/2},
    {1/2, -DodecagonalLambda/2},
    {normTwoElement/2, -normTwoElement/2},
    {DodecagonalLambda/2, -1/2},
    {DodecagonalLambda/2, 1/2}
}

(* Current selection: book orientation *)
DodecagonVertices := DodecagonVerticesBook

(********************************)
(* Projections of basis vectors *)
(********************************)

ProjectionVector[1] = 
   {1, Sqrt[3]/2, 1/2, 0}

ProjectionVector[2] = 
   {0, 1/2, Sqrt[3]/2, 1}

ProjectionVector[3] = 
   {1, -Sqrt[3]/2, 1/2, 0}

ProjectionVector[4] = 
   {0, 1/2, -Sqrt[3]/2, 1}

(*******************************************************************)
(* Parallel and orthogonal projections of four-dimensional vectors *)
(*******************************************************************)

VectorParallelProjection[v_] :=
   {v.ProjectionVector[1], v.ProjectionVector[2]}

VectorOrthogonalProjection[v_] :=
   {v.ProjectionVector[3], v.ProjectionVector[4]}

VectorOrthogonalProjectionShift[v_] :=
   {v.ProjectionVector[3], v.ProjectionVector[4]} - WindowCenterShift

VectorOrthogonalProjectionShift[v_, scale_, invert_] :=
  If[invert,
   scale*{v.ProjectionVector[3], v.ProjectionVector[4]} + WindowCenterShift,
   scale*{v.ProjectionVector[3], v.ProjectionVector[4]} - WindowCenterShift]

(****************************************************************************)
(* Transformation matrix from lattice (L12) into physical space:            *)
(* Matrix applied to a vector from (4-dimensional) lattice space gives      *)
(* the corresponding projection in physical space in the form               *)
(*         (alpha1, alpha2, beta1, beta2)                                   *)
(* where                                                                    *)
(*        alpha = alpha1 + Sqrt[3] * alpha2,                                *)
(*        beta  = beta1  + Sqrt[3] * beta2,                                 *)
(* and                                                                      *)
(*        alpha * 1 + beta * xi  (xi = Exp[2*Pi*I/12])                      *)
(* is the projected vector in Z[xi] (physical space).                       *)
(* Transformation takes into account that                                   *)
(*        VectorParallelProjection(L12) = 1.0 * Z[xi]                       *)
(* (in contrast to L8 there is no scaling). This representation can be      *)
(* used to determine visibility of the vector in physical space (by         *)
(* checking Z[Sqrt[3]]-coprimality of alpha and beta.                       *)
(****************************************************************************)

LatticeToPhysicalMatrix =
   {{1, 0, -1, 0},
    {0, 0, 0, -1},
    {0, 1, 0, 2},
    {0, 0, 1, 0}}

(*********************************)
(* Location of acceptance domain *)
(*********************************)

WindowCenterShift = 
   {0,0}

ShiftWindowCenter[vec_] :=
  Module[{},
         Unprotect[WindowCenterShift];
         WindowCenterShift = vec;
         Protect[WindowCenterShift];]


(***********************************************************************)
(* Determine whether lattice points project into the acceptance domain *)
(***********************************************************************)

PrintNumericalPrecisionWarning[functionName_] :=
   Module[{},
          Print["Warning: Insufficient accuracy in module ",
                functionName];
          Print["         Result may be wrong!"];
          Print["         To avoid this problem increase ",
                "NumericalPrecision from its current value ",
                NumericalPrecision];
         ]

(* Check on which side of the line (going through a and b) the point c is. *)
CheckPosition[a_, b_, c_] :=
   (b[[1]] - a[[1]]) * (c[[2]] - a[[2]]) - (b[[2]] - a[[2]]) * (c[[1]] - a[[1]])

CheckProjectionInWindow[point_] := 
   Module[{pt = VectorOrthogonalProjectionShift[point],
          pt1},
          pt1 = pt.pt;
          Which[N[InnerDodecagonRadius^2 - pt1,
                  NumericalPrecision] > NumericalAccuracy,
                Return[True],
                N[OuterDodecagonRadius^2 - pt1,
                  NumericalPrecision] < -NumericalAccuracy,
                Return[False],
                True,
                Return[CheckProjectionInSector[pt]]]]

CheckScaledProjInWindow[point_, normTwo_] := 
   Module[{pt, pt1, scaler},
          If[normTwo, scaler = Sqrt[DodecagonalLambda/2],
            scaler = Sqrt[DodecagonalLambda*2]];
          pt = VectorOrthogonalProjectionShift[point, scaler, normTwo];
          pt1 = pt.pt;
          Which[N[InnerDodecagonRadius^2 - pt1,
                  NumericalPrecision] > NumericalAccuracy,
                Return[True],
                N[OuterDodecagonRadius^2 - pt1,
                  NumericalPrecision] < -NumericalAccuracy,
                Return[False],
                True,
                Return[CheckProjectionInSector[pt]]]]

CheckProjectionInSector[orthpoint_] :=
   Module[{pt2, test},
          pt2 = Through[{Min, Max}[Abs[orthpoint]]];

          test = N[CheckPosition[DodecagonVertices[[2]], DodecagonVertices[[3]], pt2],
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSector"]];

          test = N[CheckPosition[DodecagonVertices[[3]], DodecagonVertices[[4]], pt2],
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSector"]];

          Return[True]];

(* Graphical check tool for the dodecagonal window test. *)
DodecagonalWindowTest[initpoint_] :=
   DynamicModule[{position = initpoint},
       {Show[Graphics[Locator[Dynamic[position]], PlotRange -> 3],
             Graphics[Line[Append[DodecagonVertices, DodecagonVertices[[1]]]]],
             ListPlot[DodecagonVertices, AspectRatio -> 1]],
       {Dynamic[position], Dynamic[CheckProjectionInSector[position]]}}];

(*************************************)
(* Construction of projection tiling *)
(*************************************)

CreateHypersteps[onlySector_] :=
   Module[{posSteps, negSteps},
          posSteps = Join[NestList[RotateRight, #, 3]&[{1,0,0,0}], {{-1,0,1,0},{0,-1,0,1}}];
          negSteps = Join[NestList[RotateRight, #, 3]&[{-1,0,0,0}], {{1,0,-1,0},{0,1,0,-1}}];

          If[onlySector,
             Return[posSteps],
             Return[Join[posSteps, negSteps]]]
         ];

(*******************************************************************)
(* Construction of visible vertices from the projection tiling:    *)
(* This construct only vertices but no adjacency information.      *)
(* Both a full vertex list and a visible vertex list is returned.  *)
(*******************************************************************)

DodecagonalProjectionTilingVis[initpoint_, maxstep_,
                               silent_:False] := 
   Module[{p,
           nop0=0,
           nop0i=1,
           nop0f=0,
           nop1=1,
           nop1i=1,
           nop1f=1,
           nop2=0,
           nop2i=2,
           nop2f=1,
           hyperstep,
           tilingpoints = {initpoint},
           visiblepoints = {},
           i,
           j,
           pp,
           pos, gcd},
          (* The restrict-to-positive-direction trick (used in the octogonal *)
          (* case) doesn't work here. Always use the full set of hypersteps. *)
          hyperstep = CreateHypersteps[False];

          If[Not[CheckProjectionInWindow[initpoint]],
             Return[Print["Initial point not in projection window"]]];
          Do[
             Do[
                p = tilingpoints[[i]];
                Do[
                   pp = p+hyperstep[[j]];
                   If[CheckProjectionInWindow[pp],
                      (* This basically searches for a copy of pp in the tilingpoints list *)
                      (* while restricting the search to the interval [nop0i, nop0f].      *)
                      pos = Position[Table[pp,{nop0}]-
                                     Take[tilingpoints,{nop0i,nop0f}],
                                     {0,0,0,0}];
                      If[pos==={},
                         pos = Position[Table[pp,{nop1}]-
                                        Take[tilingpoints,{nop1i,nop1f}],
                                        {0,0,0,0}];
                         If[pos==={},
                            pos = Position[Table[pp,{nop2}]-
                                           Take[tilingpoints,{nop2i,nop2f}],
                                           {0,0,0,0}];
                            If[pos==={},
                               nop2f++;
                               nop2++;
                               AppendTo[tilingpoints, pp];
                               gcd = GCDZ3Alt[LatticeToPhysicalMatrix.pp];
                               If[NormZ3[gcd] <= 2, If[Not[CheckScaledProjInWindow[pp, NormZ3[gcd] == 2]],
                                 AppendTo[visiblepoints, pp]]];
                               ]]]],
                   {j,Length[hyperstep]}],
                {i,nop1i,nop1f}];
             nop0  = nop1;
             nop0i = nop1i;
             nop0f = nop1f;
             nop1  = nop2;
             nop1i = nop2i;
             nop1f = nop2f;
             nop2  = 0;
             nop2i = nop1f+1;
             nop2f = nop1f,
             {maxstep}];
             If[Not[silent], Print["constructed patch of dodecagonal tiling with ",
                                   Length[tilingpoints], " vertices and ",
                                   Length[visiblepoints], " visible ones."]];
             {tilingpoints, visiblepoints}];


(* Select vertices where the GCD of the coordinates (in physical space)         *)
(* has a specific norm. We use this to analyze the validity of                  *)
(* the window-scale-plus-coprime-test, which seems to be a bit harder           *)
(* for the dodecagonal case (compared to Ammann-Beenker and the decagonal one). *)
DodecagonalSelectVertices[til_, norm_] := Module[{selfct},
    selfct = Function[u, NormZ3[GCDZ3[u[[1;;2]], u[[3;;4]]]] == norm];
    Select[til, selfct[LatticeToPhysicalMatrix.#]&]
  ];


(******************************************************)
(* Plot the radial projection distribution of the     *)
(* visible vertices obtained through the              *)
(* OctagonalProjectionTilingVis function.             *)
(******************************************************)

SectorSelectorDodecagon[{x_, y_}] := Module[{},
   If[x <= 0 || y <= 0, Return[False], Null];
   Return[(y/x) <= Sqrt[3]]]

PlotRadialProjectionDistribution[tilingpoints_List,
                                 {CutOff_,Step_}] :=
   Module[{vertices,
           radius,
           radprojections,
           mdist,
           angledists},
          (* Project to physical space and extract single (eighth) sector. *)
          vertices = Select[Map[VectorParallelProjection, tilingpoints],
                            SectorSelectorDodecagon];
          (* Compute angles of the points, then sort and compute *)
          (* differences between neighbouring angles. *)
          radprojections = Sort[Map[ArcTan[#[[1]], #[[2]]]&, vertices], Greater];
          mdist = (radprojections[[1]] - radprojections[[-1]])/(Length[radprojections] - 1);
          Print["After radial projection and sector extraction ", Length[vertices],
                " datapoints remain."];
          Clear[vertices];
          angledists = Drop[radprojections - RotateLeft[radprojections], -1];
          Clear[radprojections];
          Histogram[angledists/mdist, {0, CutOff, Step}]
]


(******************************************)
(* Plot parallel projection of tiling:    *)
(******************************************)



PlotParallelProjection[tilingpoints_List,
                       psize_:1/200,
                       onlySector_:False,
                       pointcolor_:{0,0,1}] :=
   Module[{plotpoints,
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
          If[psize <= 0, Return[Null]];
          plotpoints = N[Map[VectorParallelProjection,
                             tilingpoints]];
          If[onlySector, plotpoints = Select[plotpoints, SectorSelectorDodecagon]];

          plotpoints = Join[{color[pointcolor],
                             PointSize[psize]},
                            Map[Point, plotpoints]];
          Graphics[plotpoints, AspectRatio -> 1]]


ParallelProjectionRaw[tilingpoints_List] :=
  Map[VectorParallelProjection, tilingpoints];

CreateVertexAnnotation[vtx_] := Module[{phy, point, gcd},
  point = Point[VectorParallelProjection[vtx]];
  phy = LatticeToPhysicalMatrix.vtx;
  gcd = GCDZ3[phy[[1;;2]], phy[[3;;4]]];
  Return[Tooltip[point, vtx]];
  ];

ParallelProjectionAnnotate[tilingpoints_List] :=
  Map[CreateVertexAnnotation, tilingpoints]

ParallelProjectionDynamic[tilingpoints_List, psize_] := DynamicModule[
  {selectedVertex = {0, 0, 0, 0},
   psize2 = psize * 2},
  GraphicsRow[{
  Graphics[{Black, PointSize[psize],
    Point[Map[VectorOrthogonalProjection, tilingpoints]],
    Red, PointSize[psize2], Point[{0, 0}],
    Green, PointSize[psize2], Dynamic[Point[VectorOrthogonalProjection[selectedVertex]]]}],
  Graphics[{Black, PointSize[psize],
    Map[EventHandler[
      Point[VectorParallelProjection[#]],
      {"MouseMoved" :> {selectedVertex = #}}]&,
    tilingpoints],
    Red, PointSize[psize2], Point[{0, 0}]}]
  }]];

(****************************************)
(* Plot orthogonal projection of tiling *)
(****************************************)

PlotOrthogonalProjection[tilingpoints_List,
                         psize_:1/25,
                         pointcolor_:{0,0,1}] :=
   Module[{plotpoints,
           plotlines,
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
          If[psize <= 0, Return[Null]];
          plotpoints = N[Map[VectorOrthogonalProjectionShift,
                             tilingpoints]];
          plotpoints = Join[{color[pointcolor],
                             PointSize[psize]},
                            Map[Point,plotpoints]];
          Graphics[plotpoints, AspectRatio -> 1, PlotRange -> All]]

(* Projects the tiling vertices into internal space. *)
(* Shift of the window is ignored (mainly used for debugging) *)
OrthogonalProjectionRaw[tilingpoints_List] :=
  Map[VectorOrthogonalProjection, tilingpoints];


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

Protect[DodecagonalWindowTest,
        DodecagonalProjectionTilingVis,
        DodecagonalSelectVertices,
        PlotRadialProjectionDistribution,
        PlotParallelProjection,
        ParallelProjectionRaw,
        ParallelProjectionAnnotate,
        ParallelProjectionDynamic,
        PlotOrthogonalProjection,
        OrthogonalProjectionRaw,
        NumericalAccuracy,
        NumericalPrecision,
        CheckProjectionInWindow,
        ShiftWindowCenter,
        WindowCenterShift,
        SetNumericalPrecision]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`DodecagonalTilingCP`"    *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
