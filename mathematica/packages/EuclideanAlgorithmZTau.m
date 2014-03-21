(* ::Package:: *)

BeginPackage["EuclideanAlgorithmZTau`"]

Unprotect[GCDZTau,ModuloZTau,NormZTau,CoprimeZTau,CoprimeZTauAlt]

Begin["`Private`"]

Clear[PreNormZTau,GCDZTau,ModuloZTau,NormZTau,CoprimeZTau,CoprimeZTauAlt]

(**********************************************************************)
(* Package implements some operations in Z[Tau]. (Tau = golden ratio) *)
(* Elements are always represented as a two-dimensional               *)
(* vector (a, b) corresponding to the element a + Tau * b.            *)
(**********************************************************************)

PreNormZTau[{a_Integer, b_Integer}] :=
   a^2 + a*b - b^2

NormZTau[{a_Integer, b_Integer}] :=
   Abs[a^2 + a*b - b^2]

(* Modulo operation in Z[Tau]: *)
(* No error checking is done for cases like a % b with b equals 0. *)
ModuloZTau[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[
     {A, r, s},
     A = PreNormZTau[{a2, b2}];
     r = Round[(a1*a2 - b1*b2 + a1*b2)/A];
     s = Round[(b1*a2 - a1*b2)/A];
     Return[{a1 - (r*a2 + s*b2),
             b1 - (r*b2 + s*(a2 + b2))}]
  ];

(* Simple implementation of the Euclidean algorithm in Z[Tau]. *)
GCDZTau[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[
     {x, y, z},
     x = {a1, b1};
     y = {a2, b2};
     While[y != {0, 0},
       z = y;
       y = ModuloZTau[x, y];
       x = z;
     ];
     Return[x]
  ];

(* Check if the two elements are coprime in Z[Tau]. *)
CoprimeZTau[{a1_Integer, b1_Integer}, {a2_Integer, b2_Integer}] :=
  Module[{x},
    x = GCDZTau[{a1, b1}, {a2, b2}];
    Return[NormZTau[x]==1]
  ];

CoprimeZTauAlt[{a1_Integer, b1_Integer, a2_Integer, b2_Integer}] :=
  Module[{},
    Return[NormZTau[GCDZTau[{a1, b1}, {a2, b2}]]==1]
  ];

End[]

Protect[GCDZTau,ModuloZTau,NormZTau,CoprimeZTau,CoprimeZTauAlt]

EndPackage[]
