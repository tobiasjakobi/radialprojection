(* ::Package:: *)

BeginPackage["EuclideanAlgorithmZSqrt2`"]

Unprotect[GCDZ2,ModuloZ2,NormZ2,CoprimeZ2,CoprimeZ2Alt]

Begin["`Private`"]

Clear[PreNormZ2,GCDZ2,ModuloZ2,NormZ2,CoprimeZ2,CoprimeZ2Alt]

(******************************************************************)
(* Package implements some operations in Z[Sqrt[2]].              *)
(* Elements are always represented as a two-dimensional           *)
(* vector (a, b) corresponding to the element a + Sqrt[2] * b.    *)
(******************************************************************)

PreNormZ2[{a_Integer, b_Integer}] :=
   a^2 - 2*b^2

NormZ2[{a_Integer, b_Integer}] :=
   Abs[a^2 - 2*b^2]

(* Modulo operation in Z[Sqrt[2]]: *)
(* No error checking is done for cases like a % b with b equals 0. *)
ModuloZ2[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[
     {A, x, y, r, s},
     A := PreNormZ2[{a2, b2}];
     x := (a1*a2 - 2*b1*b2)/A;
     y := (b1*a2 - a1*b2)/A;
     r = x - Round[x];
     s = y - Round[y];
     Return[{a2*r + 2*b2*s, a2*s + b2*r}]
  ];

(* Simple implementation of the Euclidean algorithm in Z[Sqrt[2]]. *)
GCDZ2[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[
     {x, y, z},
     x = {a1, b1};
     y = {a2, b2};
     While[y != {0, 0},
       z = y;
       y = ModuloZ2[x, y];
       x = z;
     ];
     Return[x]
  ];

(* Check if the two elements are coprime in Z[Sqrt[2]]. *)
CoprimeZ2[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[{x},
    x = GCDZ2[{a1, b1}, {a2, b2}];
    Return[NormZ2[x]==1]
  ];

CoprimeZ2Alt[{a1_Integer, b1_Integer, a2_Integer, b2_Integer}] :=
  Module[{},
    Return[NormZ2[GCDZ2[{a1, b1}, {a2, b2}]]==1]
  ];

End[]

Protect[GCDZ2,ModuloZ2,NormZ2,CoprimeZ2,CoprimeZ2Alt]

EndPackage[]
