(********************************************************************)
(*                                                                  *)
(* :Title:       FibonacciChain                                     *)
(*                                                                  *)
(* :Authors:     Uwe Grimm                                          *)
(*                                                                  *)
(* :Context:     AperiodicTilings`FibonacciChain`                   *)
(*                                                                  *)
(* :Version:     1.01                                               *)
(*                                                                  *)
(* :Date:        January 12,2005                                    *)
(*                                                                  *)
(* :Summary:     Implementation of several methods to construct     *)
(*               one-dimensional aperiodic sequences, in            *) 
(*               particular the celebrated Fibonacci sequence       *)
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
(* :References:  Any textbook on quasicrystals introduces the       *)
(*               Fibonacci sequence as a paradigm of a one-         *)
(*               dimensional quasicrystal.                          *)
(*                                                                  *)
(* :Description: This package provides a number of definitions      *)
(*               that can be used to produce one-dimensional        *)   
(*               aperiodic sequences by different means,            *) 
(*               frequently refering to the example of the          *)
(*               Fibonacci sequence which is the default setting    *)  
(*               of the functions defined below. The program can    *)
(*               treat arbitrary n-letter substitution sequences,   *)   
(*               construct the Fibonacci sequence by recursion,     *) 
(*               give some statistical analysis (for instance       *)  
(*               by counting the frequences of subwords in the      *)
(*               sequence), and contains some basic definitions     *)
(*               on Fibonacci numbers. In addition, three           *)
(*               variants of the projection method are also         *)
(*               implemented, the results are shown graphically     *)
(*               and are presented in a way to make it easy         *)
(*               to understand the mechanism.                       *)
(*                                                                  *)
(* :Notes:       The graphical output can partly be changed by      *)
(*               optional arguments to the functions, if other      *)
(*               alterations are wished for these have to be        *)
(*               implemented by the user.                           *)
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

BeginPackage["AperiodicTilings`FibonacciChain`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

AtomicSurfaceProjection::usage =
"AtomicSurfaceProjection[nl,sl,sh,phas,cl] projects from the an
(nl x nl)-patch of the two-dimensional square lattice onto a one-dimensional
line with slope sl using `atomic surfaces' which are obtained as projections
of the unit square. The default values are nl=20 and sl=1/GoldenRatio,
yielding the Fibonacci sequence of short and long intervals with length
ratio GoldenRatio. The argument sh represents the shift of the cut space
with respect to the line through the origin, the default value is \"random\"
which draws a random value for the shift in the interval [-2,2]. The
cut space can be deformed by a function phas, the default is Function[0],
that is no deformation. Finally, the colors in the plots can be given
explicitly as the last argument cl, which must be a list of (at least)
seven elements. The default colors are {RGBColor[0,0,0],RGBColor[1,0,0],
RGBColor[0,1,0],RGBColor[0,0,1],RGBColor[0,1,1],RGBColor[1,1,0],
RGBColor[1,0,1]}."

FibonacciNumber::usage =
"FibonacciNumber[n] gives the n-th Fibonacci number, defined recursively
as the sum of the two previous Fibonacci numbers with initial values
FibonacciNumber[0] = 0 and FibonacciNumber[1] = 1. The definition is
extended to negative values of n."

FibonacciRecursionSequence::usage =
"FibonacciRecursionSequence[n] gives the word obtained by the concatenation
of FibonacciRecursionSequence[n-1] and FibonacciRecursionSequence[n-2],
as a list of single letters. The initial data are
FibonacciRecursionSequence[0] = {\"S\"} and
FibonacciRecursionSequence[1] = {\"L\"}."

FibonacciRecursionWord::usage =
"FibonacciRecursionWord[n] gives the word FibonacciRecursionSequence[n], but
as a string."

GoldenNumber::usage =
"GoldenNumber is the golden ratio (1+Sqrt[5])/2."

GoldenNumberApproximant::usage =
"GoldenNumberApproximant[n] gives the n-th rational approximant of
GoldenNumber defined as the ratio FibonacciNumber[n+1]/FibonacciNumber[n]."

KlotzConstruction::usage =
"KlotzConstructionAtomic[nl,sl,sh,cl] projects from the an (nl x nl)-patch
of the two-dimensional square lattice onto a one-dimensional line with slope
sl using the `Klotz contruction' or `dualization' method. Here, one first
constructs a new fundamental domain of the square lattice whose one-dimensional
boundaries are either parallel or perpendicular to the cut space. The way
this construction is performed is also shown in detail in a GraphicsArray
of four plots. The default values are nl=20 and sl=1/GoldenRatio, yielding
the Fibonacci sequence of short and long intervals with length ratio
GoldenRatio. The argument sh represents the shift of the cut space
with respect to the line through the origin, the default value is \"random\"
which draws a random value for the shift in the interval [-2,2]. Finally,
the colors in the plots can be given explicitly as the last argument cl,
which must be a list of (at least) seven elements. The default colors are
{RGBColor[0,0,0],RGBColor[1,0,0],RGBColor[0,1,0],RGBColor[0,0,1],
RGBColor[0,1,1],RGBColor[1,1,0],RGBColor[1,0,1]}."

LetterCount::usage =
"LetterCount[w] counts the number of letters in the word w, which can be
represented as a list of single letters or as a string. LetterCount[w,l]
counts the number of occurances of the single letter l in w."

LetterFrequency::usage =
"LetterFrequency[w,l] computes the frequency of the single letter l in the
word w. Here, w can be represented as a list of single letters or as a
string."

ListCount::usage =
"ListCount[w1,w2] counts the number of times the list w2 occurs as a 'sublist'
of the list w1, that is, how often w1 contains the same sequence of elements
as w2."

PlotLinearChain::usage =
"PlotLinearChain[wl,scl,ivl,lic,lec,les] plots a linear chain of intervals
of lengths given by the argument ivl, the default being {GoldenRatio,1}.
Here, wl is assumed to be a list of words, where the words are represented
either as lists of single letters, or as strings. The argument scl determines
the relative scaling of the chains for different words - if scl is False
(which is the default), the true lengths are used, if scl is True, the chain
for all words are scaled to the same total length. The colors of the line
segments and the letters can be altered by the arguments lic and lec, the
defaults are {RGBColor[1,0,0],RGBColor[0,1,0]} (red and green) for both.
Finally, the size of the letters is controlled by the argument les, the
corresponding default value is 8."

StringCount::usage =
"StringCount[w1,w2] counts the number of times the string w2 occurs as a
'substring' of the list w2, that is, how often w1 contains the same sequence
of letters as w2."

StripProjection::usage =
"StripProjection[nl,sl,sh,uc,cl] projects from the an (nl x nl)-patch
of the two-dimensional square lattice onto a one-dimensional line with
slope sl using the `strip' method. The strip is obtained by sliding
the unit square along the cut space. The default values are nl=20 and
sl=1/GoldenRatio, yielding the Fibonacci sequence of short and long
intervals with length ratio GoldenRatio. The argument sh represents the
shift of the cut space with respect to the line through the origin, the
default value is \"random\" which draws a random value for the shift in
the interval [-2,2]. The argument uc determines the position of the
unit square highlighted in the figure, its default value is 1/2. Finally,
the colors in the plots can be given explicitly as the last argument cl,
which must be a list of (at least) seven elements. The default colors are
{RGBColor[0,0,0],RGBColor[1,0,0],RGBColor[0,1,0], RGBColor[0,0,1],
RGBColor[0,1,1],RGBColor[1,1,0],RGBColor[1,0,1]}."

Substitution::usage =
"Substitution[w], defined by SubstitutionSystem, applies SubstitionRule
to the word w."

SubstitutionMatrix::usage =
"SubstitutionMatrix, defined by SubstitutionSystem, gives the substitution
matrix corresponding to the SubstitutionRule."

SubstitutionRule::usage =
"SubstitutionRule, defined by SubstitutionSystem, gives the substitution
rule of the sequence entered as the first argument to SubstitutionSystem."

SubstitutionSequence::usage =
"SubstitutionSequence[n], defined by SubstitutionSystem, gives the word
obtained by n-fold application of SubstitutionRule to the initial word
(which is defined to be SubstitutionSequence[0]), as a list of single
letters."

SubstitutionSystem::usage =
"SubstitutionSystem[subst,initword] defines SubstitutionRule,
SubstitutionMatrix, Substitution, SubstitutionSequence, and
SubstitutionWord for a substitution system. Here, subst specifies
the substitution rule, it must be a list of the replacements rules
for all letters, and initword is used as the initial word in the
substitution. The defaults are {\"S\"->{\"L\"},\"L\"->{\"L\",\"S\"}}
and {\"S\"} which yield the Fibonacci sequence. Note that the elements
of the alphabet and the initial word have to be strings."

SubstitutionWord::usage =
"SubstitutionWord[n], defined by SubstitutionSystem, gives the word
obtained by n-fold application of SubstitutionRule to the initial word
(which is defined to be SubstitutionSequence[0]), as a string."

SubWordCount::usage =
"SubWordCount[w1,w2] counts how often the word w2 occurs as a subword of
the word w1. Here, the words can be represented either as lists of single
letters, or as strings."

SubWordFrequency::usage =
"SubWordFrequency[w1,w2] determined the frequency of the word w2 as a subword
of the word w1. Here, the words can be represented either as lists of single
letters, or as strings."


(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[AtomicSurfaceProjection,
          GoldenNumber,
          GoldenNumberApproximant,
          KlotzConstruction,
          LetterCount,
          LetterFrequency,
          ListCount,
          PlotLinearChain,
          StringCount,
          StripProjection,
          SubstitutionSystem,
          SubWordCount,
          SubWordFrequency]

(********************************************************************)
(*                                                                  *)
(*                     Start of Private context                     *)
(*                                                                  *)
(********************************************************************)

Begin["`Private`"]

(*********************************)
(* Clearing previous definitions *)
(*********************************)

Clear[AtomicSurfaceProjection,
      GoldenNumber,
      GoldenNumberApproximant,
      KlotzConstruction,
      LetterCount,
      LetterFrequency,
      ListCount,
      PlotLinearChain,
      StringCount,
      StripProjection,
      SubstitutionSystem,
      SubWordCount,
      SubWordFrequency]

(**********************************)
(* Defining a substitution system *)
(**********************************)

SubstitutionSystem[subst_List:{"S"->{"L"},
                               "L"->{"L","S"}},
                   initword_:{"S"}] :=
  Module[{i,
          j,
          l=Length[subst],
          letters=Map[First,subst],
          substitutes=Map[Last,subst]},
          Clear[SubstitutionRule,
                SubstitutionMatrix,
                Substitution,
                SubstitutionSequence,
                SubstitutionWord];
          SubstitutionRule = subst;
          SubstitutionMatrix =
            Table[Count[substitutes[[j]],letters[[i]]],{i,l},{j,l}];
          Substitution[w_List] :=
            Flatten[w /. SubstitutionRule];
          SubstitutionSequence[0] =
            initword;
          SubstitutionSequence[n_Integer] :=
            SubstitutionSequence[n] =
            Substitution[SubstitutionSequence[n-1]];
          SubstitutionWord[n_Integer] :=
            Apply[StringJoin,SubstitutionSequence[n]];]

(**************************************************)
(* Recursive generation of the Fibonacci sequence *)
(**************************************************)

Clear[FibonacciRecursionSequence,
      FibonacciRecursionWord]

FibonacciRecursionSequence[0] = {"S"}

FibonacciRecursionSequence[1] = {"L"}

FibonacciRecursionSequence[n_Integer] :=
  FibonacciRecursionSequence[n] =
  Join[FibonacciRecursionSequence[n-1],
       FibonacciRecursionSequence[n-2]]

FibonacciRecursionWord[n_Integer] :=
  Apply[StringJoin,FibonacciRecursionSequence[n]]

(*********************)
(* Fibonacci numbers *)
(*********************)

Clear[FibonacciNumber]

FibonacciNumber[0] = 0

FibonacciNumber[1] = 1

FibonacciNumber[n_Integer /; Positive[n]] :=
  FibonacciNumber[n] =
  FibonacciNumber[n-1] + FibonacciNumber[n-2]

FibonacciNumber[n_Integer /; Negative[n]] :=
  FibonacciNumber[n] =
  FibonacciNumber[n+2] - FibonacciNumber[n+1]

GoldenNumber = (1+Sqrt[5])/2

GoldenNumberApproximant[n_Integer /; n>0] :=
  FibonacciNumber[n+1]/FibonacciNumber[n]

(**********************************)
(* Counting subwords of sequences *)
(**********************************)

StringCount[w1_String,
            w2_String] :=
  Module[{i,
          lw2=StringLength[w2]},
         Count[Table[StringTake[w1,{i,i+lw2-1}],
                     {i,StringLength[w1]-lw2+1}],w2]]

ListCount[w1_List,
          w2_List] :=
  Count[Partition[w1,Length[w2],1],w2]

LetterCount[w_List] :=
  Length[w]

LetterCount[w_String] :=
  StringLength[w]

LetterCount[w_List,
            letter_String] :=
  Count[w,letter]

LetterCount[w_String,
            letter_String] :=
  StringCount[w,letter]

SubWordCount[w1_String,
             w2_String] :=
  StringCount[w1,w2]

SubWordCount[w1_List,
             w2_List] :=
  ListCount[w1,w2]

LetterFrequency[w_,
               letter_String] :=
  LetterCount[w,letter]/LetterCount[w]

SubWordFrequency[w1_,
                w2_] :=
  SubWordCount[w1,w2]/LetterCount[w1]

(********************************************************************)
(* Plot of geometric Fibonacci chains with intervals of two lengths *)
(********************************************************************)

PlotLinearChain[words:{___List},
                scaletotallength_Symbol:False,
                intervallengths_List:{N[GoldenRatio],1},
                linecolorlist_List:{RGBColor[1,0,0],RGBColor[0,1,0]},
                lettercolorlist_List:{RGBColor[1,0,0],RGBColor[0,1,0]},
                lettersize_:8] :=
  Module[{i,
          interv,
          j,
          lettercolors,
          letters,
          linecolors,
          n,
          numletters,
          numwords,
          plotletters,
          plotlines,
          plotpart,
          pos,
          scalefac,
          wordlen},
         numwords     = Length[words];
         letters      = Union[Flatten[words]];
         numletters   = Length[letters];
         If[Length[intervallengths]<numletters,
            Return["Error in PlotLinearChain: Not enough lengths specified"]];
         If[Or[Length[linecolorlist]<numletters,
               Length[lettercolorlist]<numletters],
            Return["Error in PlotLinearChain: Not enough colors specified"]];
         interv       = words /.
                        MapIndexed[Rule[#1,intervallengths[[First[#2]]]]&,
                                   letters];
         linecolors   = words /.
                        MapIndexed[Rule[#1,linecolorlist[[First[#2]]]]&,
                                   letters];
         lettercolors = words /.
                        MapIndexed[Rule[#1,lettercolorlist[[First[#2]]]]&,
                                   letters];
         wordlen      = Map[Length,words];
         pos          = Map[Map[Apply[Plus,#]&,
                                Table[Take[#,i],{i,0,Length[#]}]]&,
                            interv];
         If[scaletotallength,
            scalefac  = Max[wordlen]/Map[Last,pos];
            pos       = scalefac*pos];
         plotlines = Table[{linecolors[[n,i]],
                            Line[{{pos[[n,i]],-n},
                                  {pos[[n,i+1]],-n}}]},
                           {n,numwords},{i,wordlen[[n]]}];
         plotpart  = Join[{GrayLevel[0]},
                           Table[Line[{{pos[[n,i]],-n-0.15},
                                       {pos[[n,i]],-n+0.15}}],
                                 {n,numwords},{i,wordlen[[n]]+1}]];
         plotletters   = Table[{lettercolors[[n,i]],
                                Text[FontForm[words[[n,i]],
                                              {"Times-Italic",lettersize}],
                                     {(pos[[n,i]]+pos[[n,i+1]])/2,-n-0.3}]},
                               {n,numwords},{i,wordlen[[n]]}];
         Show[Graphics[Join[plotlines,plotpart,plotletters]],
              AspectRatio -> Automatic,
              PlotRange   -> All];]

PlotLinearChain[words:{___String},
                scaletotallength_Symbol:False,
                intervallengths_List:{N[GoldenRatio],1},
                linecolorlist_List:{RGBColor[1,0,0],RGBColor[0,1,0]},
                lettercolorlist_List:{RGBColor[1,0,0],RGBColor[0,1,0]},
                lettersize_:8] :=
  PlotLinearChain[Map[Table[StringTake[#,{i}],{i,StringLength[#]}]&,words],
                  scaletotallength,
                  intervallengths,
                  linecolorlist,
                  lettercolorlist,
                  lettersize]

(********************************************************)
(* Projection from the square lattice I: `strip method' *)
(********************************************************)

StripProjection[nx_:20,
                slope_:1/GoldenRatio,
                shift_:"random",
                unitcellposition_:1/2,
                colorlist_:{RGBColor[0,0,0],
                            RGBColor[1,0,0],
                            RGBColor[0,1,0],
                            RGBColor[0,0,1],
                            RGBColor[0,1,1],
                            RGBColor[1,1,0],
                            RGBColor[1,0,1]}] :=
  Module[{clip,
          cutline,
          i,
          j,
          lattice,
          numshift,
          numslope,
          nymax,
          nymin,
          ofs=3/4,
          path,
          proclip,
          project,
          proline,
          propoint,
          prostep,
          step,
          stepcolor,
          strip,
          strippoint,
          unitcell},
         If[Length[colorlist]<7,
            Print["not enough colors specified"];
            Return[]];
         numslope  = N[slope];
         If[numslope<0,
            Print["only positive values of the slope are allowed"];
            Return[]];
         numshift = N[shift];
         If[shift==="random",
            numshift=Random[Real,{-2,2}];
            Print["note: variable shift chosen randomly: shift=",numshift],
            numshift = N[shift]];
         nymin = Floor[Min[0,-numslope+numshift]];
         nymax = Ceiling[Max[nx*numslope+numshift+1,nx*numslope]];
         clip[{x_,y_}] :=
           And[nx+1/2>=x>=-1/2,nymax>=y>=nymin];
         strippoint[{x_,y_}] :=
           And[y>=(x-1)*numslope+numshift,y<x*numslope+numshift+1];
         stepcolor[p1_,p2_] :=
           colorlist[[If[p2-p1=={1,0},2,3]]];
         step[p1_,p2_] :=
           {stepcolor[p1,p2],
            Thickness[8/(100*nx)],
            Line[{p1,p2}]};
         project[{x_,y_}] :=
           {x+numslope*y,
            numslope*(x+numslope*y)}/(1+numslope^2);
         prostep[p1_,p2_] :=
           {stepcolor[p1,p2],
            Thickness[8/(100*nx)],
            Line[{project[p1],project[p2]}],
            colorlist[[5]],
            Thickness[4/(100*nx)],
            Line[{p1,project[p1]}],
            Line[{p2,project[p2]}]};
         lattice  = Flatten[Table[{i,j},{i,0,nx},{j,nymin,nymax}],1];
         cutline  = {colorlist[[4]],
                     Thickness[8/(100*nx)],
                     Line[{{-ofs,-ofs*numslope},
                           {nx+ofs,(nx+ofs)*numslope}}]};
         strip    = {colorlist[[6]],
                     Polygon[{{-ofs,-(ofs+1)*numslope+numshift},
                              {nx+ofs,(nx+ofs-1)*numslope+numshift},
                              {nx+ofs,(nx+ofs)*numslope+numshift+1},
                              {-ofs,-ofs*numslope+numshift+1}}]};
         propoint = Sort[Select[lattice,strippoint[#]&]];
         unitcell = {colorlist[[7]],
                     Polygon[{{#,#*numslope+numshift},
                              {#+1,#*numslope+numshift},
                              {#+1,#*numslope+numshift+1},
                              {#,#*numslope+numshift+1}}]}&[unitcellposition];
         path     = Table[step[propoint[[i]],propoint[[i+1]]],
                          {i,Length[propoint]-1}];
         proclip  = Apply[And,
                          Partition[Map[clip[project[#]]&,propoint],2,1],1];
         proline  = Complement[Table[If[proclip[[i]],
                                        prostep[propoint[[i]],
                                                propoint[[i+1]]]],
                                     {i,Length[propoint]-1}],{Null}];
         lattice  = Join[{colorlist[[1]],PointSize[2/(10*nx)]},
                         Map[Point,lattice]];
         Show[GraphicsArray[{{Graphics[Join[strip,unitcell,lattice,cutline],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]},
                             {Graphics[Join[strip,path,lattice,cutline],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]},
                             {Graphics[Join[strip,path,lattice,proline],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]}}
                           ]];]

(************************************************************)
(* Projection from the square lattice II: `atomic surfaces' *)
(************************************************************)

AtomicSurfaceProjection[nx_:20,
                        slope_:1/GoldenRatio,
                        shift_:"random",
                        phasonstrain_:Function[0],
                        colorlist_:{RGBColor[0,0,0],
                                    RGBColor[1,0,0],
                                    RGBColor[0,1,0],
                                    RGBColor[0,0,1],
                                    RGBColor[0,1,1],
                                    RGBColor[1,1,0],
                                    RGBColor[1,0,1]}] :=
  Module[{colors,
          cut,
          cutline,
          cutph,
          ddot,
          hyper,
          i,
          j,
          lattice,
          lsh,
          nymax,
          nymin,
          numshift,
          numslope,
          ofs=3/4,
          path,
          prohyper,
          proline,
          prolineph,
          propoint,
          result,
          sh},
         If[Length[colorlist]<7,
            Print["not enough colors specified"];
            Return[]];
         numslope  = N[slope];
         If[numslope<0,
            Print["only positive values of the slope are allowed"];
            Return[]];
         If[shift==="random",
            numshift=Random[Real,{-2,2}];
            Print["note: variable shift chosen randomly: shift=",numshift],
            numshift  = N[shift]];
         nymin = Floor[Min[0,-numslope+numshift]];
         nymax = Ceiling[Max[nx*numslope+numshift+1,nx*numslope]];
         sh    = {numslope*(1+numslope),-(1+numslope)}/(2*(1+numslope^2));
         lsh   = Dot[sh,sh];
         ddot[x_] :=
           Dot[x,x];
         cut[{x_,y_}] :=
           {(x+numslope*(y-numshift))/(1+numslope^2),
            numslope*(x+numslope*(y-numshift))/(1+numslope^2)+numshift};
         cutph[{x_,y_}] :=
           (#+phasonstrain[First[#]]*sh)&[cut[{x,y}]];
         lattice   = Flatten[Table[{i,j},{i,0,nx},{j,nymin,nymax}],1];
         cutline   = {colorlist[[1]],
                      Thickness[8/(100*nx)],
                      Line[{{-ofs,-ofs*numslope},
                            {nx+ofs,(nx+ofs)*numslope}}]};
         proline   = {colorlist[[4]],
                      Thickness[8/(100*nx)],
                      Line[Table[{i,i*numslope+numshift},
                                 {i,-ofs,nx+ofs,1/100}]]};
         prolineph = {colorlist[[4]],
                      Thickness[8/(100*nx)],
                      Line[Table[{i,i*numslope+numshift}+
                                 phasonstrain[i]*sh,
                                 {i,-ofs,nx+ofs,1/100}]]};
         hyper     = {colorlist[[7]],
                      Thickness[4/(100*nx)],
                      Table[Line[{lattice[[i]]-sh,lattice[[i]]+sh}],
                            {i,Length[lattice]}]};
         propoint  = Select[lattice,ddot[#-cutph[#]]<=lsh&];
         prohyper  = {colorlist[[5]],
                      Thickness[4/(100*nx)],
                      Table[Line[{propoint[[i]]-sh,propoint[[i]]+sh}],
                            {i,Length[propoint]}]};
         colors    = Table[colorlist[[Which[Abs[propoint[[i+1]]-
                                                propoint[[i]]]=={1,0},
                                            2,
                                            Abs[propoint[[i+1]]-
                                                propoint[[i]]]=={0,1},
                                            3,
                                            True,
                                            6]]],
                           {i,Length[propoint]-1}];
         path      = Join[{Thickness[8/(100*nx)]},
                          Table[{colors[[i]],
                                Line[{propoint[[i]],propoint[[i+1]]}]},
                                {i,Length[propoint]-1}]];
         result    = Join[{Thickness[10/(100*nx)]},
                          Table[{colors[[i]],
                                Line[{cut[propoint[[i]]],
                                      cut[propoint[[i+1]]]}]},
                               {i,Length[propoint]-1}],
                          {Thickness[5/(100*nx)],GrayLevel[0]},
                          Table[Line[{cut[propoint[[i]]]-sh/4,
                                      cut[propoint[[i]]]+sh/4}],
                                {i,Length[propoint]}]];
         lattice   = Join[{colorlist[[1]],PointSize[2/(10*nx)]},
                          Map[Point,lattice]];
         propoint  = Join[{colorlist[[5]],PointSize[2/(10*nx)]},
                          Map[Point,propoint]];
         Show[GraphicsArray[{{Graphics[Join[hyper,lattice,cutline,prolineph,
                                            proline],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]},
                             {Graphics[Join[hyper,lattice,cutline,prolineph,
                                            proline,prohyper,propoint],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]},
                             {Graphics[Join[hyper,lattice,cutline,prolineph,
                                            path,prohyper,propoint],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]},
                             {Graphics[Join[hyper,lattice,prohyper,propoint,
                                            cutline,prolineph,result],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]},
                             {Graphics[Join[hyper,lattice,prohyper,propoint,
                                            result],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax-1+ofs}}]}}
                           ]];]


(****************************************************************)
(* Projection from the square lattice III: `klotz construction' *)
(****************************************************************)

KlotzConstruction[nx_:20,
                  slope_:1/GoldenRatio,
                  shift_:"random",
                  colorlist_:{RGBColor[0,0,0],
                              RGBColor[1,0,0],
                              RGBColor[0,1,0],
                              RGBColor[0,0,1],
                              RGBColor[0,1,1],
                              RGBColor[1,1,0],
                              RGBColor[1,0,1]}] :=
  Module[{cutline,
          fundomain,
          i,
          j,
          kloetze,
          klotz1,
          klotz2,
          lattice,
          newdomain,
          numc1,
          numc2,
          nymax,
          nymin,
          numshift,
          numslope,
          ofs=3/4,
          orthlines,
          plotgrid,
          plotline,
          profundom,
          proklotz,
          proline,
          result,
          strip1,
          strip2,
          tiling,
          x,
          y},
         If[Length[colorlist]<7,
            Print["not enough colors specified"];
            Return[]];
         numslope  = N[slope];
         If[Or[numslope>1,numslope<0],
            Print["the value of the slope must lie between zero and one"];
            Return[]];
         If[shift==="random",
            numshift=Random[Real,{-2,2}];
            Print["note: variable shift chosen randomly: shift=",numshift],
            numshift  = N[shift]];
         numc1     = 1/(1+numslope^2);
         numc2     = numslope/(1+numslope^2);
         nymin     = Floor[Min[0,-numslope+numshift]];
         nymax     = Ceiling[Max[nx*numslope+numshift+1,nx*numslope]];
         plotgrid  = Join[{colorlist[[1]],Thickness[7/1000]},
                          Map[Line,Flatten[Table[{{{0,i},{3,i}},{{i,0},{i,3}}},
                                                 {i,0,3}],1]]];
         plotline  = {colorlist[[4]],Thickness[1/100],
                      Line[{{0,1-numslope},{3,1+2*numslope}}]};
         fundomain = {colorlist[[6]],
                      Polygon[{{1,1},{2,1},{2,2},{1,2}}]};
         profundom = {colorlist[[6]],
                      Polygon[{{1,1},{1+numc1,1+numc2},{2-numslope,2},{1,2}}],
                      Polygon[{{2-numslope,1},{2,1},{1+numc1,1+numc2},
                               {1+numc1-numc2,numc1+numc2}}],
                      colorlist[[5]],
                      Polygon[{{1,2},{2-numslope,2},
                               {1+numc1-numc2,1+numc1+numc2}}],
                      Polygon[{{1,1},{2-numslope,1},
                               {1+numc1-numc2,numc1+numc2}}],
                      colorlist[[3]],
                      Polygon[{{1,1},{1,2},{1-numc2,1+numc1}}],
                      Polygon[{{2,1},{2,2},{2-numc2,1+numc1}}],
                      colorlist[[2]],
                      Polygon[{{2-numslope,2},{2-numc2,1+numc1},{2,2}}],
                      Polygon[{{2-numslope,1},{2-numc2,numc1},{2,1}}]};
         orthlines = {colorlist[[7]],
                      Line[{{0,1-2*numslope},{3,1+1*numslope}}],
                      Line[{{0,2-1*numslope},{3,2+2*numslope}}],
                      Line[{{0,2-2*numslope},{3,2+1*numslope}}],
                      Line[{{0,1+1/numslope},{3,1-2/numslope}}],
                      Line[{{0,1+2/numslope},{3,1-1/numslope}}],
                      Line[{{0,2+1/numslope},{3,2-2/numslope}}],
                      Line[{{0,0+2/numslope},{3,0-1/numslope}}]};
         newdomain[x_,y_] :=
           {colorlist[[7]],
            Line[{{x,y},{x+numc1,y+numc2},
                  {x+numc1-numc2,y+numc1+numc2},
                  {x-numc2,y+numc1},{x,y}}],
            Line[{{x+1,y},{x+numc1,y+numc2},
                  {x+numc1-numc2,y-1+numc1+numc2},
                  {x+1-numc2,y-1+numc1},{x+1,y}}]};
         klotz1[x_,y_] :=
           {colorlist[[6]],
            Polygon[{{x,y},{x+numc1,y+numc2},
                     {x+numc1-numc2,y+numc1+numc2},
                     {x-numc2,y+numc1}}]};
         klotz2[x_,y_] :=
           {colorlist[[6]],
            Polygon[{{x+1,y},{x+numc1,y+numc2},
                     {x+numc1-numc2,y-1+numc1+numc2},
                     {x+1-numc2,y-1+numc1}}]};
         kloetze[x_,y_] :=
           Join[klotz1[x,y],klotz2[x,y]];
         Show[GraphicsArray[{{Graphics[Join[fundomain,plotgrid,plotline],
                              PlotRange   -> {{0,3},{0,3}},
                              AspectRatio -> Automatic]},
                             {Graphics[Join[fundomain,plotgrid,orthlines,
                                            plotline],
                              PlotRange   -> {{0,3},{0,3}},
                              AspectRatio -> Automatic]},
                             {Graphics[Join[profundom,plotgrid,plotline,
                                            newdomain[1,1]],
                              PlotRange   -> {{0,3},{0,3}},
                              AspectRatio -> Automatic]},
                             {Graphics[Join[kloetze[1,1],plotgrid,plotline,
                                            Table[newdomain[i,j],
                                                  {i,0,3},{j,0,3}]],
                              PlotRange   -> {{0,3},{0,3}},
                              AspectRatio -> Automatic]}}]];
         lattice   = Flatten[Table[{i,j},{i,0,nx},{j,nymin,nymax}],1];
         cutline   = {colorlist[[1]],Thickness[8/(100*nx)],
                      Line[{{-ofs,-ofs*numslope},
                            {nx+ofs,(nx+ofs)*numslope}}]};
         proline   = {colorlist[[4]],Thickness[8/(100*nx)],
                      Line[{{-ofs,-ofs*numslope+numshift},
                            {nx+ofs,(nx+ofs)*numslope+numshift}}]};
         tiling    = Join[{Thickness[4/(100*nx)]},
                          Flatten[Table[newdomain[i,j],
                                        {i,-1,nx},{j,nymin-1,nymax}]]];
         strip1    = Select[lattice,
                            -1<=(#[[2]]-numslope*#[[1]]-numshift)<=0&];
         strip2    = Select[lattice,
                            0<=(#[[2]]-numslope*#[[1]]-numshift)<=numslope&];
         proklotz  = Join[Map[klotz1[#[[1]],#[[2]]]&,strip1],
                          Map[klotz2[#[[1]],#[[2]]]&,strip2]];
         result    = Join[{colorlist[[2]]},
                          Map[Line[{{(#[[2]]-numshift)*numc2+#[[1]]*numc1,
                                     (numslope*#[[2]]+#[[1]])*numc2+
                                     numshift*numc1},
                                    {(#[[2]]-numshift)*numc2+(#[[1]]+1)*numc1,
                                     (numslope*#[[2]]+#[[1]]+1)*numc2+
                                     numshift*numc1}}]&,
                              strip1],
                          {colorlist[[3]]},
                          Map[Line[{{(#[[2]]-numshift-1)*numc2+
                                     (#[[1]]+1)*numc1,
                                     (numslope*(#[[2]]-1)+#[[1]]+1)*numc2+
                                     numshift*numc1},
                                    {(#[[2]]-numshift)*numc2+(#[[1]]+1)*numc1,
                                     (numslope*#[[2]]+#[[1]]+1)*numc2+
                                     numshift*numc1}}]&,
                              strip2]];
         lattice   = Join[{colorlist[[1]],PointSize[2/(10*nx)]},
                          Map[Point,lattice]];
         Show[GraphicsArray[{{Graphics[Join[tiling,lattice,cutline,proline],
                                       Axes -> True,
                                       AspectRatio -> Automatic,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax+ofs}}]},
                             {Graphics[Join[proklotz,tiling,lattice,cutline,
                                            proline],
                                       AspectRatio -> Automatic,
                                       Axes -> True,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax+ofs}}]},
                             {Graphics[Join[proklotz,tiling,lattice,cutline,
                                            result],
                                       AspectRatio -> Automatic,
                                       Axes -> True,
                                       PlotRange -> {{-ofs,nx+ofs},
                                                     {nymin-ofs,nymax+ofs}}]}}
                           ]];]


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

Protect[AtomicSurfaceProjection,
        GoldenNumber,
        GoldenNumberApproximant,
        KlotzConstruction,
        LetterCount,
        LetterFrequency,
        ListCount,
        PlotLinearChain,
        StringCount,
        StripProjection,
        SubstitutionSystem,
        SubWordCount,
        SubWordFrequency]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`FibonacciChain`"         *)
(*                                                                  *)
(********************************************************************)

EndPackage[]

