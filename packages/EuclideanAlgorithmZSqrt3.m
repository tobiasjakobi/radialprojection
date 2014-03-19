(* ::Package:: *)

BeginPackage["EuclideanAlgorithmZSqrt3`"]

Unprotect[GCDZ3,GCDZ3Alt,ModuloZ3,DivZ3,MultZ3,Norm2Type,NormZ3,CoprimeZ3,CoprimeZ3Alt]

Begin["`Private`"]

Clear[ValZ3,PreNormZ3,GCDZ3,GCDZ3Alt,ModuloZ3,DivZ3,MultZ3,Norm2Type,NormZ3,CoprimeZ3,CoprimeZ3Alt]

(******************************************************************)
(* Package implements some operations in Z[Sqrt[3]].              *)
(* Elements are always represented as a two-dimensional           *)
(* vector (a, b) corresponding to the element a + Sqrt[3] * b.    *)
(******************************************************************)

ValZ3[{a_Integer, b_Integer}] := a + b*Sqrt[3];

PreNormZ3[{a_Integer, b_Integer}] :=
   a^2 - 3*b^2

NormZ3[{a_Integer, b_Integer}] :=
   Abs[a^2 - 3*b^2]

(* Modulo operation in Z[Sqrt[3]]: *)
(* No error checking is done for cases like a % b with b equals 0. *)
ModuloZ3[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[
     {A, r, s},
     A = PreNormZ3[{a2, b2}];
     r = Round[(a1*a2 - 3*b1*b2)/A];
     s = Round[(b1*a2 - a1*b2)/A];
     Return[{a1 - (r*a2 + 3*s*b2),
             b1 - (r*b2 + s*a2)}]
  ];

DivZ3[{a_Integer, b_Integer}, {c_Integer, d_Integer}] :=
  Module[{A},
    A = PreNormZ3[{c, d}];
    Return[{(a*c - 3*b*d)/A,
            (b*c - a*d)/A}];
  ];

MultZ3[{a_Integer, b_Integer}, {c_Integer, d_Integer}] :=
  {a*c + 3*b*d, a*d + b*c};

(* Element x with absolute norm 2 can be written as *)
(*      x = plusminus lambda^k phi                   *)
(* where lambda = 2 + sqrt[3] is the fundamental    *)
(* unit and phi = 1 + sqrt[3] is the fundamental    *)
(* element with absnorm 2. The function Norm2Type   *)
(* computes the exponent k of this representation.  *)
Norm2Type[{a_Integer, b_Integer}] := Module[{x, k=0, mult},
  If[NormZ3[{a, b}] != 2, Return[Null]];
  If[ValZ3[{a, b}] < 0, x = {-a, -b}, x = {a, b}];
  If[Abs[ValZ3[x]] >= ValZ3[{5,3}], mult = {2, -1}, mult = {2, 1}];
  While[x != {1, 1}, x = MultZ3[x, mult]; k++];
  If[mult == {2, 1}, Return[-k], Return[k]];
  ];

(* Simple implementation of the Euclidean algorithm in Z[Sqrt[3]]. *)
GCDZ3[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[
     {x, y, z},
     x = {a1, b1};
     y = {a2, b2};
     While[y != {0, 0},
       z = y;
       y = ModuloZ3[x, y];
       x = z;
     ];
     Return[x]
  ];

GCDZ3Alt[{a1_Integer, b1_Integer, a2_Integer, b2_Integer}] :=
  GCDZ3[{a1, b1}, {a2, b2}];

(* Check if the two elements are coprime in Z[Sqrt[3]]. *)
CoprimeZ3[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[{x},
    x = GCDZ3[{a1, b1}, {a2, b2}];
    Return[NormZ3[x]==1]
  ];

CoprimeZ3Alt[{a1_Integer, b1_Integer, a2_Integer, b2_Integer}] :=
  (NormZ3[GCDZ3[{a1, b1}, {a2, b2}]]==1);

End[]

Protect[GCDZ3,GCDZ3Alt,ModuloZ3,DivZ3,MultZ3,Norm2Type,NormZ3,CoprimeZ3,CoprimeZ3Alt]

EndPackage[]

