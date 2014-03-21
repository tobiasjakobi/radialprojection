(* ::Package:: *)

BeginPackage["AperiodicTilings`DecagonalTilingCP`",{"EuclideanAlgorithmZTau`"}]

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[DecagonalProjectionTilingVis,
          PlotRadialProjectionDistribution,
          PlotParallelProjection,
          PlotOrthogonalProjection,
          NumericalAccuracy,
          NumericalPrecision,
          CheckProjectionInWindow,
          CheckScaledProjInWindow,
          DecagonWindowTest,
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
      InnerOctagonRadius,
      OuterOctagonRadius,
      DecagonVertices,
      DecagonVerticesAlt,
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
      CheckProjectionInSectorAlt,
      DoDecagonWindowTest,
      CreateHypersteps,
      DecagonalProjectionTilingVis,
      PlotRadialProjectionDistribution,
      SectorSelectorDecagon,
      PlotParallelProjection,
      PlotOrthogonalProjection]


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

(***********************************************)
(* Inner and outer radius of acceptance domain *)
(***********************************************)

InnerDecagonRadius =
   (GoldenRatio + 1)/2

OuterDecagonRadius =
   GoldenRatio^2/Sqrt[GoldenRatio^2 + 1]

DecagonVertices =
   {{OuterDecagonRadius, 0}, {Sqrt[(11*GoldenRatio + 7)/5]/2, GoldenRatio/2},
    {Sqrt[(GoldenRatio + 2)/5]/2, InnerDecagonRadius}, {-Sqrt[(GoldenRatio + 2)/5]/2, InnerDecagonRadius},
    {-Sqrt[(11*GoldenRatio + 7)/5]/2, GoldenRatio/2}, {-OuterDecagonRadius, 0},
    {-Sqrt[(11*GoldenRatio + 7)/5]/2, -GoldenRatio/2}, {-Sqrt[(GoldenRatio + 2)/5]/2, -InnerDecagonRadius},
    {Sqrt[(GoldenRatio + 2)/5]/2, -InnerDecagonRadius}, {Sqrt[(11*GoldenRatio + 7)/5]/2, -GoldenRatio/2}}

DecagonVerticesAlt =
   {{InnerDecagonRadius, Sqrt[(GoldenRatio+1)/(GoldenRatio+2)]/2},
    {GoldenRatio/2, Sqrt[(8*GoldenRatio+5)/(GoldenRatio+2)]/2},
    {0, OuterDecagonRadius},
    {-GoldenRatio/2, Sqrt[(8*GoldenRatio+5)/(GoldenRatio+2)]/2},
    {-InnerDecagonRadius, Sqrt[(GoldenRatio+1)/(GoldenRatio+2)]/2},
    {-InnerDecagonRadius, -Sqrt[(GoldenRatio+1)/(GoldenRatio+2)]/2},
    {-GoldenRatio/2, -Sqrt[(8*GoldenRatio+5)/(GoldenRatio+2)]/2},
    {0, -OuterDecagonRadius},
    {GoldenRatio/2, -Sqrt[(8*GoldenRatio+5)/(GoldenRatio+2)]/2},
    {InnerDecagonRadius, -Sqrt[(GoldenRatio+1)/(GoldenRatio+2)]/2}}

(********************************)
(* Projections of basis vectors *)
(********************************)

ProjectionVector[1] =
   {1, (-1 + Sqrt[5])/4, (-1 - Sqrt[5])/4, (-1 - Sqrt[5])/4}

ProjectionVector[2] =
   {0, Sqrt[(5 + Sqrt[5])/8], Sqrt[(5 - Sqrt[5])/8], -Sqrt[(5 - Sqrt[5])/8]}

ProjectionVector[3] =
   {1, (-1 - Sqrt[5])/4, (-1 + Sqrt[5])/4, (-1 + Sqrt[5])/4}

ProjectionVector[4] =
   {0, Sqrt[(5 - Sqrt[5])/8], -Sqrt[(5 + Sqrt[5])/8], Sqrt[(5 + Sqrt[5])/8]}

(*******************************************************************)
(* Parallel and orthogonal projections of four-dimensional vectors *)
(*******************************************************************)

VectorParallelProjection[v_] :=
   {v.ProjectionVector[1],v.ProjectionVector[2]}

VectorOrthogonalProjection[v_] :=
   {v.ProjectionVector[3],v.ProjectionVector[4]}

VectorOrthogonalProjectionShift[v_] :=
   {v.ProjectionVector[3],v.ProjectionVector[4]} - WindowCenterShift

(****************************************************************************)
(* Transformation matrix from lattice (L5) into physical space:             *)
(* Matrix applied to a vector from (4-dimensional) lattice space gives      *)
(* the corresponding projection in physical space in the form               *)
(*         (alpha1, alpha2, beta1, beta2)                                   *)
(* where                                                                    *)
(*        alpha = alpha1 + Tau * alpha2,                                    *)
(*        beta  = beta1  + Tau * beta2,                                     *)
(* and                                                                      *)
(*        alpha * 1 + beta * xi  (xi = Exp[2*Pi*I/5])                       *)
(* is the projected vector in Z[xi] (physical space).                       *)
(* Transformation takes into account that                                   *)
(*        VectorParallelProjection(L5) = 1.0 * Z[xi]                        *)
(* (in contrast to L8 there is no scaling). This representation can be      *)
(* used to determine visibility of the vector in physical space (by         *)
(* checking Z[Tau]-coprimality of alpha and beta).                          *)
(****************************************************************************)

LatticeToPhysicalMatrix =
   {{1, 0, -1, 1},
    {0, 0, 0, -1},
    {0, 1, -1, 1},
    {0, 0, 1, -1}}

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
          Which[N[InnerDecagonRadius^2 - pt1,
                  NumericalPrecision] > NumericalAccuracy,
                Return[True],
                N[OuterDecagonRadius^2 - pt1,
                  NumericalPrecision] < -NumericalAccuracy,
                Return[False],
                True,
                Return[CheckProjectionInSectorAlt[pt]]]]

CheckScaledProjInWindow[point_] := 
   Module[{pt = GoldenRatio*VectorOrthogonalProjectionShift[point],
           pt1},
          pt1 = pt.pt;
          Which[N[InnerDecagonRadius^2 - pt1,
                  NumericalPrecision] > NumericalAccuracy,
                Return[True],
                N[OuterDecagonRadius^2 - pt1,
                  NumericalPrecision] < -NumericalAccuracy,
                Return[False],
                True,
                Return[CheckProjectionInSectorAlt[pt]]]]

CheckProjectionInSector[orthpoint_] :=
   Module[{pt2 = Abs[orthpoint], test},

          test = N[CheckPosition[DecagonVertices[[1]], DecagonVertices[[2]], pt2],
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSector"]];

          test = N[CheckPosition[DecagonVertices[[2]], DecagonVertices[[3]], pt2],
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSector"]];

          test = N[CheckPosition[DecagonVertices[[3]], DecagonVertices[[4]], pt2], 
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSector"]];

          Return[True]];

CheckProjectionInSectorAlt[orthpoint_] :=
   Module[{pt2 = Abs[orthpoint], test},

          test = N[Simplify[CheckPosition[DecagonVerticesAlt[[10]], DecagonVerticesAlt[[1]], pt2]],
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSectorAlt"]];

          test = N[Simplify[CheckPosition[DecagonVerticesAlt[[1]], DecagonVerticesAlt[[2]], pt2]],
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSectorAlt"]];

          test = N[Simplify[CheckPosition[DecagonVerticesAlt[[2]], DecagonVerticesAlt[[3]], pt2]], 
                   NumericalPrecision];
          Which[test < -NumericalAccuracy, Return[False],
                test > NumericalAccuracy, Null,
                True, PrintNumericalPrecisionWarning["CheckProjectionInSectorAlt"]];

          Return[True]];

(* Graphical check tool for the decagonal window test. *)
DecagonWindowTest[initpoint_, alt_] :=
   DynamicModule[{position = initpoint},
     If[alt, {{Show[Graphics[Locator[Dynamic[position]], PlotRange -> 2],
               Graphics[Line[Append[DecagonVerticesAlt, DecagonVerticesAlt[[1]]]]],
               ListPlot[DecagonVerticesAlt, AspectRatio -> 1]]},
              {Dynamic[position], Dynamic[CheckProjectionInSectorAlt[position]]}},
             {{Show[Graphics[Locator[Dynamic[position]], PlotRange -> 2],
                Graphics[Line[Append[DecagonVertices, DecagonVertices[[1]]]]],
                ListPlot[DecagonVertices, AspectRatio -> 1]]},
              {Dynamic[position], Dynamic[CheckProjectionInSector[position]]}}]];

(*************************************)
(* Construction of projection tiling *)
(*************************************)

CreateHypersteps[onlySector_] :=
   Module[{posSteps, negSteps},
          posSteps = Append[NestList[RotateRight, #, 3]&[{1,0,0,0}], {1,1,1,1}];
          negSteps = Append[NestList[RotateRight, #, 3]&[{-1,0,0,0}], {-1,-1,-1,-1}];

          If[onlySector,
             Return[posSteps],
             Return[Join[posSteps, negSteps]]]
         ];

(*******************************************************************)
(* Construction of visible vertices from the projection tiling:    *)
(* This construct only vertices but no adjacency information.      *)
(* Both a full vertex list and a visible vertex list is returned.  *)
(*******************************************************************)

DecagonalProjectionTilingVis[initpoint_, maxstep_,
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
           pos},
          (* The restrict-to-positive-direction trick (used in the octogonal *)
          (* case) doesn't work here. Always use the full set of hypersteps. *)
          hyperstep = CreateHypersteps[False];

          If[Not[CheckProjectionInWindow[initpoint]],
             Return[Print["Initial point not in projection window"]]];
          Do[
             Do[
                p = tilingpoints[[i]];
                Do[
                   pp = p + hyperstep[[j]];
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
                               If[Not[CheckScaledProjInWindow[pp]] &&
                                  CoprimeZTauAlt[LatticeToPhysicalMatrix.pp],
                                  AppendTo[visiblepoints, pp]]]]]],
                   {j, Length[hyperstep]}],
                {i, nop1i, nop1f}];
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
             If[Not[silent], Print["constructed patch of octagonal tiling with ",
                                   Length[tilingpoints], " vertices and ",
                                   Length[visiblepoints], " visible ones."]];
             {tilingpoints, visiblepoints}];

(******************************************************)
(* Plot the radial projection distribution of the     *)
(* visible vertices obtained through the              *)
(* DecagonalProjectionTilingVis function.             *)
(******************************************************)

PlotRadialProjectionDistribution[tilingpoints_List,
                                 {CutOff_, Step_}] :=
   Module[{vertices,
           radius,
           radprojections,
           mdist,
           angledists},
          (* TODO: Extract sector after projection to physical space. *)
          vertices = Map[VectorParallelProjection, tilingpoints];
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

SectorSelectorDecagon[{x_, y_}] := Module[{},
   If[x <= 0 || y <= 0, Return[False], Null];
   Return[(y/x) <= Sqrt[5 + 2*Sqrt[5]]]]

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
          If[onlySector, plotpoints = Select[plotpoints, SectorSelectorDecagon]];

          plotpoints = Join[{color[pointcolor],
                             PointSize[psize]},
                            Map[Point, plotpoints]];
          Graphics[plotpoints, AspectRatio -> 1]]

(****************************************)
(* Plot orthogonal projection of tiling *)
(****************************************)

PlotOrthogonalProjection[tilingpoints_List,
                         psize_:1/25,
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
          plotpoints = N[Map[VectorOrthogonalProjectionShift,
                             tilingpoints]];
          plotpoints = Join[{color[pointcolor],
                             PointSize[psize]},
                            Map[Point,plotpoints]];
          Graphics[plotpoints, AspectRatio -> 1, PlotRange -> All]]


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

Protect[DecagonalProjectionTilingVis,
        PlotRadialProjectionDistribution,
        PlotParallelProjection,
        PlotOrthogonalProjection,
        NumericalAccuracy,
        NumericalPrecision,
        CheckProjectionInWindow,
        CheckScaledProjInWindow,
        DecagonWindowTest,
        ShiftWindowCenter,
        WindowCenterShift,
        SetNumericalPrecision]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`DecagonalTilingCP`"      *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
