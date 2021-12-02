(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "cosine"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*T*)


(* ::Text:: *)
(*See (r5.10) (Page r5-2).*)


tKnown[b_][x_, y_] := 1 - b Cos[x] Cosh[y];


(* ::Subsubsection:: *)
(*Location of straight contour*)


xStraight = Pi/2;


tKnown[\[FormalCapitalB]][xStraight, \[FormalY]] == 1


(* ::Subsubsection:: *)
(*P = \[PartialD]T/\[PartialD]x, Q = \[PartialD]T/\[PartialD]y*)


(* ::Text:: *)
(*See (r5.13) and (r5.14) (Page r5-2).*)


p[b_][x_, y_] := D[tKnown[b][x, y], x] // Evaluate;
q[b_][x_, y_] := D[tKnown[b][x, y], y] // Evaluate;


With[{x = \[FormalX], y = \[FormalY], b = \[FormalCapitalB]},
  {
    p[b][x, y] == b Sin[x] Cosh[y],
    q[b][x, y] == -b Cos[x] Sinh[y]
  }
]


(* ::Subsection:: *)
(*Flux function F*)


(* ::Text:: *)
(*See (r5.15) (Page r5-2).*)


f[a_, b_][x_, y_] := -tKnown[b][x, y]^4 / a // Evaluate;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


(* ::Text:: *)
(*See (r5.17) (Page r5-2).*)


vi[a_, b_][x_, y_] := p[b][x, y]^2 + q[b][x, y]^2 - f[a, b][x, y]^2 // Evaluate;


(* ::Subsection:: *)
(*Simple case (B = 1)*)


(* ::Subsubsection:: *)
(*Critical terminal points along y = 0*)


(* ::Text:: *)
(*Observe that \[CapitalPhi](y = 0) = (1 - c^2) - (1 - c)^8 / A^2, where c = cos(x).*)
(*The critical terminal point at x = 0 is trivial.*)


x0Simp[a_] :=
  With[{c = \[FormalC]},
    ArcCos[c]
      /. Last @ Solve[
        {
          (1 - c^2) - (1 - c)^8 / a^2 == 0,
          c < 1
        }, c, Reals
      ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Representative values of A*)


aValuesSimp = {1/10, 1, 2};


aNamesSimp = AssociationThread[
  aValuesSimp,
  {"small", "med", "large"}
];


aValuesSimpConvex = {8/10, 9/10};


aNamesSimpConvex = AssociationThread[
  aValuesSimpConvex,
  {"convex_8", "convex_9"}
];


(* ::Subsection:: *)
(*General case (B arbitrary)*)


(* ::Subsubsection:: *)
(*Critical terminal points along x = 0*)


(* ::Text:: *)
(*Along the y-axis, there are critical terminal points at y = \[PlusMinus]y_0 when B < 1.*)
(*When B = 1, \[PlusMinus]y_0 merges with x_\[Flat] at the origin.*)
(*When B > 1, we don't really care because \[PlusMinus]y_0 lie in the unphysical region*)
(*(note that things are a bit iffy because the non-viable region containing the origin*)
(*might join with other non-viable blobs for B sufficiently large).*)


(* ::Subsubsubsection:: *)
(*Critical terminal point y_0*)


(* ::Text:: *)
(*Observe that \[CapitalPhi](x = 0) = -(1 - B C)^8 / A^2 + B^2 (C^2 - 1), where C = cosh(y).*)
(*Note that for B < 1, the critical terminal point y_0 lies in the strictly physical region,*)
(*i.e. at C = cosh(y) < 1/B.*)


With[{c = \[FormalC], a = \[FormalCapitalA], b = \[FormalCapitalB]},
  Solve[
    {
      -(1 - b c)^8 / a^2 + b^2 (c^2 - 1) == 0,
      1 < c < 1/b,
      0 < a,
      0 < b < 1
    },
    c,
    Reals
  ]
]


y0CriticalA[b_] := Root[
  -65536
  + 458752 b^2
  - 1376256 b^4
  + 2293760 b^6
  - 2293760 b^8
  + 1376256 b^10
  - 458752 b^12
  + 65536 b^14
  + (729 - 77200 b^2 - 449856 b^4 - 283392 b^6 - 13824 b^8) #^2
  + 729 b^2 #^4 &,
  2
]


y0PolyC[a_, b_] := Function[
  1 + a^2 b^2 - 8 b #
  + (28 b^2 - a^2 b^2) #^2
  - 56 b^3 #^3
  + 70 b^4 #^4
  - 56 b^5 #^5
  + 28 b^6 #^6
  - 8 b^7 #^7
  + b^8 #^8
]


y0[a_, b_] /; 0 < b < 1 && 0 < a < y0CriticalA[b] := ArcCosh @ Root[y0PolyC[a, b], 1];
y0[a_, b_] /; 0 < b < 1 && a >= y0CriticalA[b] := ArcCosh @ Root[y0PolyC[a, b], 3];


(* ::Subsubsection:: *)
(*Critical terminal points along y = 0*)


(* ::Text:: *)
(*Observe that \[CapitalPhi](y = 0) = B^2 (1 - c^2) - (1 - B c)^8 / A^2, where c = cos(x).*)
(*For a given A, there exists a constant B_\[Natural](A) such that*)
(*only for B > B_\[Natural](A) will there be two critical terminal x-values,*)
(*x_\[Flat] < x_\[Sharp], otherwise there will be none.*)


(* ::Subsubsubsection:: *)
(*Exploratory analysis: solutions of the polynomial*)


With[{c = \[FormalC], a = \[FormalCapitalA], b = \[FormalCapitalB]},
  Solve[
    {
      b^2 (1 - c^2) - (1 - b c)^8 / a^2 == 0,
      a > 0,
      b > 0
    }, c, Reals
  ]
]


(* ::Subsubsubsection:: *)
(*Critical value B_\[Natural](A)*)


(* ::Text:: *)
(*From the above output,*)
(*we see that the critical constant B_\[Natural](A) is given by*)
(*the edge of the inequality A > Root[-65536 + ... + 729 B^2 #1^4 &, 2].*)
(*Inverting that expression for B in terms of A:*)


With[{a = \[FormalCapitalA], b = \[FormalCapitalB]},
  Solve[
    {
      a == Root[
        -65536 + 458752 b^2 - 1376256 b^4 + 2293760 b^6 - 2293760 b^8
        + 1376256 b^10 - 458752 b^12 + 65536 b^14
        + (-729 + 77200 b^2 + 449856 b^4 + 283392 b^6 + 13824 b^8) #1^2
        + 729 b^2 #1^4 &,
        2
      ],
      b > 0
    }, b, Reals
  ]
]


(* ::Text:: *)
(*Therefore:*)


bNat[a_] := Root[
  -65536 - 729 a^2
  + (458752 + 77200 a^2 + 729 a^4) #1^2
  + (-1376256 + 449856 a^2) #1^4
  + (2293760 + 283392 a^2) #1^6
  + (-2293760 + 13824 a^2) #1^8
  + 1376256 #1^10
  - 458752 #1^12
  + 65536 #1^14 &,
  2
]


(* ::Subsubsubsection:: *)
(*Critical terminal points x_\[Flat] and x_\[Sharp]*)


(* ::Text:: *)
(*From the above output (see Exploratory analysis: solutions of the polynomial),*)
(*we see that the two values of critical terminal cos(x)*)
(*are given by the Root[1 - A^2 B^2 - ... + B^8 #1^8 &, ...] expressions:*)


polyC[a_, b_] := Function[{c},
  1 - a^2 b^2
  - 8 b c
  + (28 b^2 + a^2 b^2) c^2
  - 56 b^3 c^3
  + 70 b^4 c^4
  - 56 b^5 c^5
  + 28 b^6 c^6
  - 8 b^7 c^7
  + b^8 c^8
];


xFlat[a_, b_] := ArcCos @ Root[polyC[a, b], 2] // Evaluate;


xSharp[a_, b_] := ArcCos @ Root[polyC[a, b], 1] // Evaluate;


(* ::Subsubsubsection:: *)
(*Critical x_\[Natural] at B = B_\[Natural](A)*)


xNat[a_] := xSharp[a, bNat[a]] // Evaluate;


(* ::Subsubsection:: *)
(*Representative values of A*)


aValuesGen = {3, 10};


aNamesGen = AssociationThread[
  aValuesGen,
  {"three", "ten"}
];


(* ::Subsubsection:: *)
(*Representative values of B for a given A*)


(* ::Text:: *)
(*The names of the regimes are:*)
(*(1) gentle, B < B_\[Natural](A)*)
(*(2) gentle-to-fair, B = B_\[Natural](A)*)
(*(3) fair, B_\[Natural](A) < B < 1*)
(*(4) fair-to-steep (or simple), B = 1*)
(*(5) steep, B > 1*)
(*See Page r5-8.*)


(* Gentle regime: choose B == 9/10 B_\[Natural](A) *)
bValueGen["gentle"][a_] := 9/10 bNat[a];


(* Gentle-to-fair transition: B == B_\[Natural](A) *)
bValueGen["gentle_fair"][a_] := bNat[a];


(* Fair regime: B_\[Natural] < B < 1 such that x_\[Flat] == 3/4 x_\[Natural] *)
bValueGen["fair"][a_] :=
  SeekRootBisection[
    xFlat[a, #] - 3/4 xNat[a] &,
    {bNat[a], 1}
  ];


(* Fair-to-steep transition: B == 1 *)
bValueGen["fair_steep"][a_] := 1;


(* Steep regime: choose B == 2 *)
bValueGen["steep"][a_] := 2;


(* ::Subsubsection:: *)
(*Representative pair (A, B) for an asymmetric domain*)


aAsymm = 12;
bAsymm = 104/100 bNat[aAsymm];


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* ::Subsubsection:: *)
(*System of ODES for terminal curve \[CapitalPhi] = 0*)


viContourSystem[a_, b_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{fun, p, q, slope},
      (* Scalar function \[CapitalPhi] whose contours are sought *)
      fun = vi[a, b][x, y];
      (* Components of the gradient vector *)
      p = D[fun, x];
      q = D[fun, y];
      (* Magnitude of the gradient vector *)
      slope = Sqrt[p^2 + q^2];
      (* Return system of ODEs *)
      {
        x' == q / slope,
        y' == -p / slope
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*System of ODES for T-contours*)


xyContourSystem[b_] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    With[
     {p = p[b][x, y],
      q = q[b][x, y]
     },
      Module[{slope, xDer, yDer},
        (* Magnitude of the gradient vector *)
        slope = Sqrt[p^2 + q^2];
        (* Return system of ODEs *)
        {
          x' == q / slope,
          y' == -p / slope
        } /. {
          x' -> x'[s],
          y' -> y'[s],
          x -> x[s],
          y -> y[s]
        }
      ]
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* Non-trivial critical terminal point (x_0, 0), which is hyperbolic *)
Table[
  startXYSimp[a]["hyperbolic"] = {
    {x0Simp[a], 0}
  };
, {a, aValuesSimp}];


(* Starting points along contour through (1/2 x_0, 0) *)
Table[
  startXYSimp[a]["contour"] =
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      Module[{b, nMax, yMax, sMax, xyContour},
        b = 1;
        nMax = 8;
        yMax = 2;
        (* (Probable) upper bound for arc length traversed *)
        sMax = 3/2 * yMax;
        (* Contour *)
        xyContour =
          NDSolveValue[
            {
              xyContourSystem[b],
              x[0] == 1/2 x0Simp[a], y[0] == 0,
              WhenEvent[Abs @ y[s] > yMax, "StopIntegration"]
            }, {x, y}, {s, -sMax, sMax},
            NoExtrapolation
          ];
        (* Actual arc length traversed *)
        sMax = DomainEnd[xyContour];
        (* Starting points along the contour *)
        Table[
          xyContour[s] // Through // Rationalize[#, 0] &
        , {s, Subdivide[-sMax, sMax, nMax]}]
      ]
    ];
, {a, aValuesSimp}];


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Subsubsubsection:: *)
(*Gentle regime B < B_\[Natural](A)*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "gentle";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along (probably) disconnected contour
        T = 2/3 T(x = x_\[Natural], y = 0)
      (For the current values A = 3, 10, this contour is disconnected,
      but I say "probably" since there might be some extreme choice of A & B
      for which the contour is not disconnected.)
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 2/3 tKnown[b][xNat[a], 0];
    (* Seek starting point (x == 1, y) therealong *)
    xInit = 1;
    yInit = SeekRoot[
      tKnown[b][xInit, #] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 6;
    startXYGen[a][regime]["disconnected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Join[#, yReflect /@ #] &
    );
    (*
      ----------------------------------------------------------------
      Starting points along (probably) connected contour
        T = T(x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    (* Starting point (x == x_\[Natural], y == 0) *)
    xInit = xNat[a];
    yInit = 0;
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {Abs @ y[s] > yMax},
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points within the viable domain *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Cases[{x_, y_} /; vi[a, b][x, y] > 0]
    );
    (*
      ----------------------------------------------------------------
      Critical terminal points along y-axis
        (x == 0, y == y_0)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["vertical"] = {
      {0, y0[a, b]},
      {0, -y0[a, b]}
    };
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Gentle-to-fair transition B = B_\[Natural](A)*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "gentle_fair";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along (probably) disconnected contour
        T = 2/3 T(x = x_\[Natural], y = 0)
      (For the current values A = 3, 10, this contour is disconnected,
      but I say "probably" since there might be some extreme choice of A & B
      for which the contour is not disconnected.)
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 2/3 tKnown[b][xNat[a], 0];
    (* Seek starting point (x == 1, y) therealong *)
    xInit = 1;
    yInit = SeekRoot[
      tKnown[b][xInit, #] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 6;
    startXYGen[a][regime]["disconnected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Join[#, yReflect /@ #] &
    );
    (*
      ----------------------------------------------------------------
      Starting points along connected contour
        T = T(x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    (* Starting point (x == x_\[Natural], y == 0) *)
    xInit = xNat[a];
    yInit = 0;
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {Abs @ y[s] > yMax},
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points within the viable domain *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Critical terminal points along y-axis
        (x == 0, y == y_0)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["vertical"] = {
      {0, y0[a, b]},
      {0, -y0[a, b]}
    };
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal point (x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xNat[a], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Fair regime B_\[Natural](A) < B < 1*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "fair";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along (probably) disconnected contour
        T = 2/3 T(x = x_\[Natural], y = 0)
      (For the current values A = 3, 10, this contour is disconnected,
      but I say "probably" since there might be some extreme choice of A & B
      for which the contour is not disconnected.)
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 2/3 tKnown[b][xNat[a], 0];
    (* Seek starting point (x == 1, y) therealong *)
    xInit = 1;
    yInit = SeekRoot[
      tKnown[b][xInit, #] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 6;
    startXYGen[a][regime]["disconnected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Join[#, yReflect /@ #] &
    );
    (*
      ----------------------------------------------------------------
      Starting points along (probably) connected contour
        T = T(x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    (* Starting point (x == x_\[Natural], y == 0) *)
    xInit = xNat[a];
    yInit = 0;
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {Abs @ y[s] > yMax},
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points within the viable domain *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Critical terminal points along y-axis
        (x == 0, y == y_0)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["vertical"] = {
      {0, y0[a, b]},
      {0, -y0[a, b]}
    };
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal points
        (x = x_\[Flat], y = 0)
        (x = x_\[Sharp], y = 0)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xFlat[a, b], 0},
      {xSharp[a, b], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Fair-to-steep transition B = 1*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "fair_steep";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along connected contour
        T = 1/2
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 1/2;
    (* Seek starting point (x, y == 0) therealong *)
    yInit = 0;
    xInit = SeekRoot[
      tKnown[b][#, yInit] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal point (x = x_\[Sharp], y = 0)
      (The trivial point (x = x_\[Flat], y = 0) is useless.)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xSharp[a, b], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Steep regime B > 1*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "steep";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along connected contour
        T = 1/2
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 1/2;
    (* Seek starting point (x, y == 0) therealong *)
    yInit = 0;
    xInit = SeekRoot[
      tKnown[b][#, yInit] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal point (x = x_\[Sharp], y = 0)
      (The point (x = x_\[Flat], y = 0) lies in the unphysical region T < 0.)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xSharp[a, b], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Starting point groups for each regime*)


idListGen =
  With[{d = "disconnected", c = "connected", h = "hyperbolic", v = "vertical"},
    Association[
      "gentle" -> {d, c, v},
      "gentle_fair" -> {d, c, h, v},
      "fair" -> {d, c, h, v},
      "fair_steep" -> {c, h},
      "steep" -> {c, h}
    ]
  ];


(* ::Subsubsubsection:: *)
(*List of regimes*)


regimeListGen = Keys[idListGen];


(* ::Subsection:: *)
(*Terminal curve x = x(s), y = y(s)*)


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


xyTermSimp[a_?NumericQ] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{b, sMax},
      b = 1;
      sMax = 5;
      NDSolveValue[
        {
          viContourSystem[a, b],
          x[0] == x0Simp[a], y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsection:: *)
(*Traced boundaries x = x(y)*)


(* ::Subsubsection:: *)
(*Traced boundary derivative dx/dy*)


(* ::Text:: *)
(*See (r5.18) (Page r5-3).*)
(*By default the lower branch is used since this corresponds*)
(*to the candidate boundary for y > 0.*)


xTraDer[a_, b_, upperBranch_: False] := Function[{x, y},
  With[
   {p = p[b][x, y],
    q = q[b][x, y],
    f = f[a, b][x, y],
    vi = vi[a, b][x, y],
    sign = If[upperBranch, -1, +1]
   },
    Divide[
      p q + sign f Sqrt[vi],
      q^2 - f^2
    ]
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*Traced boundary curvature*)


(* ::Text:: *)
(*See (r5.20) (Page r5-5).*)


curTra[a_, b_, upperBranch_: False] := Function[{x, y},
  Module[{d, xDer, xDer2, cur},
    (* Abbreviation for y-derivative *)
    d = Dt[#, y, Constants -> {a, b}] &;
    (* x' and x'' *)
    (*
      Note: these are by default for the lower branch,
      since xTraDer uses the lower branch
     *)
    xDer = xTraDer[a, b, upperBranch][x, y];
    xDer2 = d[xDer];
    (* Curvature evaluated *)
    cur = xDer2 /. {d[x] -> xDer};
    cur = cur /. {d @ If[_, _, _] -> 0}
    (*
      Note: ugly d @ If[_, _, _] replacement needed
      since upperBranch returns an If expression.
      d[_If] doesn't work since it evaluates.
      d @ If[___] raises the warning If::argbu.
     *)
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*System of ODES for inflection frontier \[Kappa] = 0*)


curContourSystem[a_, b_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{fun, p, q, slope},
      (* Scalar function \[Kappa] whose contours are sought *)
      (*
        Note: this is for the lower branch,
        since curTra uses the lower branch
       *)
      fun = curTra[a, b][x, y];
      (* Components of the gradient vector *)
      p = D[fun, x];
      q = D[fun, y];
      (* Magnitude of the gradient vector *)
      slope = Sqrt[p^2 + q^2];
      (* Return system of ODEs *)
      {
        x' == Re[q / slope],
        y' == Re[-p / slope]
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* ::Subsubsubsection:: *)
(*Candidate boundary*)


xTraCandSimp[a_?NumericQ, terminateAtStraightContour_: False] :=
  With[{x = \[FormalX]},
    Module[{b, yMax},
      b = 1;
      yMax = 5;
      NDSolveValue[
        {
          x'[y] == Re @ xTraDer[a, b][x[y], y],
          x[0] == x0Simp[a],
          WhenEvent[
            terminateAtStraightContour && x[y] > xStraight,
            "StopIntegration"
          ]
        }, x, {y, 0, yMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsubsection:: *)
(*Critical value y = y(A) for inflection along x = \[Pi]/2.*)


(* ::Text:: *)
(*See (r5.24) (Page r5-6).*)


yCurCritSimp[a_] :=
  With[{y = \[FormalY], s = \[FormalCapitalS]},
    ArcSinh[s] /. First @ Solve[
      {
        2 s - (1 + s^2) (a^2 s + 4 Sqrt[a^2 (1 + s^2) - 1]) == 0,
        0 < a < 1
      },
      s, Reals
    ]
  ] // Evaluate;


yCurCritSimp[\[FormalCapitalA]]


(* ::Subsubsubsection:: *)
(*A_i (inflection dimensionless group)*)


(* ::Text:: *)
(*See (r5.25) (Page r5-6).*)


(* Compute A_i using the bisection algorithm *)
(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
aInflSimp = Module[{dest, aMin, aMax, a, num},
  dest = "cosine_simple-a-inflection.txt";
  ExportIfNotExists[dest,
    (* A_i is around 0.8 *)
    aMin = 7/10;
    aMax = 9/10;
    (* Solve x(y(A)) - \[Pi]/2 == 0 *)
    SeekRootBisection[
      xTraCandSimp[#] @ yCurCritSimp[#] - Pi / 2 &,
      {aMin, aMax},
      "ReturnIterations" -> True
    ] // Compress
  ];
  (* Import *)
  {a, num} = Import[dest] // Uncompress;
  (* Print iterations used *)
  Print @ FString["Bisection algorithm: {num} iterations"];
  (* Return A_i *)
  a
]


(* ::Subsubsubsection:: *)
(*Approximate determination of A_i*)


(* ::Text:: *)
(*See (r5.29) (Page r5-7).*)


aInflSimpApprox =
  With[{a = \[FormalA]},
    a /. First @ Solve[
      {
        a^6 - a^4 + 44 a^2 - 28 == 0,
        a > 0
      },
      a, Reals
    ]
  ];
aInflSimpApprox // N


(* ::Subsubsubsection:: *)
(*Convex domain aspect ratio*)


(* ::Text:: *)
(*Ratio of height to width of the thin, lens-like convex domains.*)


aspectRatioSimp[a_] :=
  Module[{xRad, xStart, xEnd, yStart, yEnd},
    (* Traced boundary *)
    xRad = xTraCandSimp[a, True];
    xStart = xRad[0];
    xEnd = xStraight;
    yStart = 0;
    yEnd = DomainEnd[xRad];
    (* Return aspect ratio *)
    2 (yEnd - yStart) / (xEnd - xStart)
  ];


(* ::Subsubsubsection:: *)
(*Non-convex lens-like candidate self-radiation upper bound*)


(* ::Text:: *)
(*See Pages r6-4 and r6-5 of manuscripts/radiation-6-self.pdf,*)
(*in particular (r6.32) for the ultra-crude upper bound for R.*)


(* Compute A vs R{ultra} *)
(* (This is not terribly slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
candidateBoundaryRatioTable =
  Module[
    {
      dest,
      aMin, aMax,
      numValues, aValues,
      b,
      aRTable,
      xCandidate, yEnd,
      xDer, xDer2, xDer3,
      yInfl, xInfl,
      yView, xView,
      yEx, xEx,
      yXDer2Max, xDer2Max,
      yXDerMin, xDerMin,
      tMin, tMax,
      rBound,
      rCutoff,
      dummyForTrailingCommas
    },
    dest = "cosine_simple-candidate-boundary-ratio-bound-ultra.txt";
    ExportIfNotExists[dest,
      (* Values of A for sampling *)
      aMin = 0;
      aMax = aInflSimp;
      numValues = 50;
      aValues = Subdivide[aMin, aMax, numValues + 1] // Rest // Most;
      (* Value of B *)
      b = 1;
      (* For each value of A *)
      Table[
        (* Compute candidate boundary *)
        xCandidate = xTraCandSimp[a, True];
        (* Get ending y-coordinate y_e *)
        yEnd = DomainEnd[xCandidate];
        (* Analytic expressions for derivatives (better than numeric) *)
        xDer[y_] := xTraDer[a, b][xCandidate[y], y];
        xDer2[y_] := curTra[a, b][xCandidate[y], y];
        xDer3[y_] := (
          (* dF/dy == \[PartialD]F/\[PartialD]x * dx/dy + \[PartialD]F/\[PartialD]y *)
          Derivative[1, 0][curTra[a, b]][xCandidate[y], y] * xDer[y]
          + Derivative[0, 1][curTra[a, b]][xCandidate[y], y]
        );
        (* Find inflection y-coordinate y_i *)
        yInfl = SeekRoot[xDer2, {0, yEnd}, 32] // Quiet;
        xInfl = xCandidate[yInfl];
        (* Find self-viewing extremity y-coordinate y_v (see (r6.29)) *)
        yView = SeekRoot[xCandidate[#] + (yEnd - #) xDer[#] - Pi/2 &, {0, yInfl}];
        xView = xCandidate[yView];
        (* Extremum for second derivative *)
        yEx = SeekRoot[xDer3, {0, yEnd}, 32] // Quiet;
        xEx = xCandidate[yEx];
        (* Maximum absolute value for second derivative *)
        yXDer2Max =
          First @ MaximalBy[
            {yView, yEnd, If[yView < yEx < yEnd, yEx, Nothing]},
            Abs @* xDer2,
            1
          ];
        xDer2Max = Abs @ xDer2[yXDer2Max];
        (* Minimum absolute value for first derivative *)
        yXDerMin =
          First @ MinimalBy[
            {yView, yEnd, yInfl},
            Abs @* xDer,
            1
          ];
        xDerMin = Abs @ xDer[yXDerMin];
        (* Minimum and maximum temperature *)
        (* NOTE: only need endpoints since dT/dy is never zero *)
        {tMin, tMax} = MinMax @ Table[tKnown[b][xCandidate[y], y], {y, {yView, yEnd}}];
        (* Compute ultra-crude bound for self-incident boundary ratio R (see (r6.32)) *)
        rBound =
          Divide[
            (yEnd - yView)^2 * (tMax / tMin)^4 * xDer2Max^2,
            8 (1 + xDerMin^2) ^ 2
          ];
        (* Return pair of values (A, R) *)
        {a, rBound}
        , {a, aValues}
      ] // Compress
    ];
    (* Import *)
    Import[dest] // Uncompress
  ];


(* ::Subsubsubsection:: *)
(*Cutoff where non-convex lens-like candidates are yet practical*)


rPracCandSimp = 1/100; (* 1 percent *)
aPracCandSimp =
  Module[{rUltraInterpolation, aMin, aMax},
    rUltraInterpolation = Interpolation[candidateBoundaryRatioTable];
    aMin = DomainStart[rUltraInterpolation];
    aMax = DomainEnd[rUltraInterpolation];
    SeekRoot[rUltraInterpolation[#] - rPracCandSimp &, {aMin, aMax}]
  ]


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Subsubsubsection:: *)
(*Candidate natural boundary (B = B_\[Natural](A), through x = x_\[Natural])*)


xTraCandNatGen[a_?NumericQ, terminateAtStraightContour_: False] :=
  With[{x = \[FormalX]},
    Module[{b, yMax},
      b = bNat[a];
      yMax = 5;
      NDSolveValue[
        {
          x'[y] == Re @ xTraDer[a, b][x[y], y],
          x[0] == xNat[a],
          WhenEvent[
            terminateAtStraightContour && x[y] > xStraight,
            "StopIntegration"
          ]
        }, x, {y, 0, yMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsection:: *)
(*Traced boundaries x = x(s), y = y(s)*)


(* ::Subsubsection:: *)
(*System of ODES for tracing*)


xyTraSystem[a_, b_, upperBranch_: True] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    With[
     {p = p[b][x, y],
      q = q[b][x, y],
      f = f[a, b][x, y],
      vi = vi[a, b][x, y],
      sign = If[upperBranch, +1, -1]
     },
      Module[{grad2, xDer, yDer},
        (* Square of gradient *)
        grad2 = p^2 + q^2;
        (* Return system of ODEs *)
        xDer = (-q f + sign p Re @ Sqrt[vi]) / grad2;
        yDer = (+p f + sign q Re @ Sqrt[vi]) / grad2;
        {x' == xDer, y' == yDer} /. {
          x' -> x'[s],
          y' -> y'[s],
          x -> x[s],
          y -> y[s]
        }
      ]
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* ::Subsubsubsection:: *)
(*Generic boundaries (for the chosen representative values of A)*)


Module[
 {b,
  idList,
  sMax,
  xyInitList, viTol,
  xInit, yInit
 },
  b = 1;
  idList = {"contour", "hyperbolic"};
  sMax = 4;
  Table[
    Table[
      xyInitList = startXYSimp[a][id];
      viTol = 10^-6;
      (*
        NOTE: leniency is required;
        when travelling away from y == 0,
        the traced boundaries asymptote very quickly
        towards the terminal curve.
       *)
      xyTraSimp[a][id] =
        Table[
          {xInit, yInit} = xyInit;
          With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
            NDSolveValue[
              {
                xyTraSystem[a, b],
                x[0] == xInit, y[0] == yInit,
                WhenEvent[
                  {
                    vi[a, b][x[s], y[s]] < viTol,
                    tKnown[b][x[s], y[s]] < 0
                  },
                  "StopIntegration"
                ]
              }, {x, y}, {s, -sMax, sMax},
              NoExtrapolation
            ]
          ]
        , {xyInit, xyInitList}];
    , {id, idList}]
  , {a, aValuesSimp}]
];


(* ::Subsubsubsection:: *)
(*Candidate boundary*)


(*
  "Cand" means candidate,
  referring to the traced boundary through (x_0, 0)
  which is almost the same as the terminal curve.
  The upper branch is taken for y > 0;
  the lower branch is taken for y < 0.
 *)
xyTraCandSimp[a_?NumericQ, terminateAtStraightContour_: False] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{b, sMax},
      b = 1;
      sMax = 5;
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == x0Simp[a], y[0] == 0,
          WhenEvent[
            terminateAtStraightContour && x[s] > xStraight,
            "StopIntegration"
          ]
        }, {x, y}, {s, 0, sMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Subsubsubsection:: *)
(*Generic boundaries (for the chosen representative values of A)*)


Module[
 {aValues,
  b, idList, sMax,
  xyInitList, xInit, yInit
 },
  (* Representative values of A *)
  aValues = aValuesGen;
  (* For each value of A *)
  Table[
    (* For each regime *)
    Table[
      (* Value of B *)
      b = bValueGen[regime][a];
      (* Starting point groups *)
      idList = idListGen[regime];
      (* Range for s *)
      sMax = 6;
      (* For each group of starting points *)
      Table[
        (* Starting points *)
        xyInitList = startXYGen[a][regime][id];
        (* Solve for traced boundaries *)
        xyTraGen[a][regime][id] =
          Table[
            {xInit, yInit} = xyInit;
            With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
              NDSolveValue[
                {
                  xyTraSystem[a, b],
                  x[0] == xInit, y[0] == yInit,
                  WhenEvent[
                    {
                      vi[a, b][x[s], y[s]] < 0,
                      tKnown[b][x[s], y[s]] < 0
                    },
                    "StopIntegration"
                  ]
                }, {x, y}, {s, -sMax, sMax},
                NoExtrapolation
              ]
            ]
          , {xyInit, xyInitList}]
      , {id, idList}]
    , {regime, regimeListGen}]
  , {a, aValues}]
];


(* ::Subsubsubsection:: *)
(*Constructing an asymmetric domain*)


(* ::Text:: *)
(*Non-trivial intersections (i.e. not (x_\[Flat], 0) or (x_\[Sharp], 0)) for the lower branch*)
(*  between the inflection frontier and the x-axis (y = 0), "axis inflection", and*)
(*  between the inflection frontier and the straight boundary (x = \[Pi]/2), "straight inflection".*)


Module[
 {a, b,
  xFl, xSh, xMid,
  sMax,
  xyInflAxis, xyInflStraight
 },
  a = aAsymm;
  b = bAsymm;
  (* Useful constants *)
  xFl = xFlat[a, b];
  xSh = xSharp[a, b];
  xMid = Way[xFl, xSh];
  (*
    ------------------------------------------------
    Axis inflection
    ------------------------------------------------
    xInflAxis "x_i(axis)":
      non-trivial x along y == 0 where inflection occurs
    yInflAxis "y_i(axis)":
      y where traced boundary through (x_i(axis), 0) intersects x == \[Pi]/2
   *)
  xInflAxis =
    SeekFirstRootBisection[
      curTra[a, b][#, 0] &,
      {xMid, xSh}
    ];
  sMax = 5;
  xyInflAxis =
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInflAxis, y[0] == 0,
          WhenEvent[
            x[s] > xStraight,
            "StopIntegration"
          ]
        }, {x, y}, {s, 0, sMax},
        NoExtrapolation
      ]
    ];
  yInflAxis = -xyInflAxis[[2]] @ DomainEnd[xyInflAxis];
  (*
    ------------------------------------------------
    Straight inflection
    ------------------------------------------------
    yInflStraight "y_i(straight)":
      non-trivial y along x == \[Pi]/2 where inflection occurs
    xInflStraight "x_i(straight)":
      x where traced boundary through (\[Pi]/2, y_i(straight)) intersects y == 0
   *)
  yInflStraight =
    SeekFirstRootBisection[
      curTra[a, b][xStraight, #] &,
      {0.4, 0.8}
    ];
  sMax = 5;
  xyInflStraight =
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xStraight, y[0] == -yInflStraight,
          WhenEvent[
            y[s] > 0,
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, 0},
        NoExtrapolation
      ]
    ];
  xInflStraight = xyInflStraight[[1]] @ DomainStart[xyInflStraight];
]


(* ::Text:: *)
(*The asymmetric domain shall have corner at*)
(*  x-coordinate 1/10 of the way from x_i(straight) to x_i(axis), and*)
(*  y-coordinate 3/4 of the way from 0 to the negative y, along that x,*)
(*    where inflection occurs for the lower branch.*)


Module[
 {a, b,
  yInflNegative,
  idList, upperBranch,
  sMax, sLower, sUpper
 },
  a = aAsymm;
  b = bAsymm;
  (*
    ------------------------------------------------
    Asymmetric domain corner
    ------------------------------------------------
    xAsymmCorner:
      1/10 of the way from x_i(straight) to x_i(axis)
    yAsymmCorner:
      3/4 of the way from 0 to the negative y, along that x,
      where inflection occurs for the lower branch.
   *)
  xAsymmCorner = Way[xInflStraight, xInflAxis, 1/10];
  yInflNegative = SeekFirstRootBisection[
    curTra[a, b][xAsymmCorner, #] &,
    {0, -0.5}
  ];
  yAsymmCorner = Way[0, yInflNegative, 3/4];
  (*
    ------------------------------------------------
    Asymmetric domain boundaries
    ------------------------------------------------
   *)
  idList = {"upper", "lower"};
  sMax = 2;
  Table[
    upperBranch = id === "upper";
    sLower = If[upperBranch, 0, -sMax];
    sUpper = If[upperBranch, sMax, 0];
    xyTraAsymm[id] =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyTraSystem[a, b, upperBranch],
            x[0] == xAsymmCorner, y[0] == yAsymmCorner,
            WhenEvent[
              x[s] > xStraight,
              "StopIntegration"
            ]
          }, {x, y}, {s, sLower, sUpper},
          NoExtrapolation
        ]
      ];
  , {id, idList}];
];


(* ::Subsection:: *)
(*Numerical verification (finite elements)*)


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* ::Text:: *)
(*Here we consider the thin, lens-like convex domains.*)


(* ::Subsubsubsection:: *)
(*Generate finite element mesh*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[
 {dest,
  yReflect,
  xyRad, xyRadEvenExt,
  xStart, sEnd, xEnd, yEnd,
  sSpacing, bPointList, nB,
  bMesh, mesh,
  prRad, prBath
 },
  Table[
    dest = FString[
      "cosine_simple-verification-mesh-{aNamesSimpConvex[a]}.txt"
    ];
    ExportIfNotExists[dest,
      (* Reflect y coordinate *)
      yReflect = # * {1, -1} &;
      (* Radiation boundary x == x(s), y == y(s) *)
      (* NOTE: y decreases with increasing s *)
      xyRad = xyTraCandSimp[a, True];
      (* Left-hand extremity thereof *)
      xStart = xyRad[[1]][0];
      (* Upper endpoint thereof *)
      sEnd = DomainEnd[xyRad];
      {xEnd, yEnd} = xyRad[sEnd] // Through // yReflect;
      (* Check that endpoint x value equals \[Pi]/2 *)
      If[xEnd == xStraight, Null,
        Print @ FString @ StringJoin[
          "WARNING (a == {N[a]}): ",
          "xEnd == {N[xEnd]} does not equal \[Pi]/2"
        ]
      ];
      (* Make spacing of boundary points 1/5 of the width *)
      sSpacing = 1/5 (xEnd - xStart);
      (* Boundary points *)
      xyRadEvenExt[s_] /; s >= 0 := xyRad[s] // Through;
      xyRadEvenExt[s_] /; s < 0 := xyRadEvenExt[-s] // yReflect;
      bPointList = Join[
        (* Radiation boundary (x(s), y(s)) *)
        Table[
          xyRadEvenExt[s]
        , {s, UniformRange[-sEnd, sEnd, sSpacing]}],
        (* Straight boundary x == \[Pi]/2 *)
        Table[
          {xStraight, y}
        , {y, UniformRange[-yEnd, yEnd, sSpacing]}]
          // Rest // Most
      ];
      bPointList = bPointList;
      nB = Length[bPointList];
      (* Build boundary element mesh *)
      bMesh = ToBoundaryMesh[
        "Coordinates" -> bPointList,
        "BoundaryElements" -> {
          LineElement[
            Table[{n, n + 1}, {n, nB}] // Mod[#, nB, 1] &
          ]
        }
      ];
      (* Build mesh *)
      mesh = ToElementMesh[bMesh,
        "ImproveBoundaryPosition" -> True
      ];
      (* Predicate functions for radiation and bath (straight) boundaries *)
      prRad = Function[{x, y}, x < xStraight // Evaluate];
      prBath = Function[{x, y}, x == xStraight // Evaluate];
      (* Export *)
      {a, mesh, prRad, prBath} // Compress
    ]
  , {a, aValuesSimpConvex}]
]


(* ::Subsubsubsection:: *)
(*Solve boundary value problem*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Table[
  Module[
   {dest, source,
    a, mesh, prRad, prBath,
    tSol
   },
    dest = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    ExportIfNotExists[dest,
      (* Import mesh *)
      source = FString[
        "cosine_simple-verification-mesh-{aNamesSimpConvex[a]}.txt"
      ];
      {a, mesh, prRad, prBath} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* Radiation boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prRad[x, y]],
            (* Heat bath (straight edge) boundary condition *)
            DirichletCondition[t[x, y] == 1, prBath[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      ];
      tSol // Compress
    ]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Text:: *)
(*Here we consider the asymmetric convex domain constructed above.*)


(* ::Subsubsubsection:: *)
(*Generate finite element mesh*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[
 {dest,
  xyUpper, xyLower, yTop, yBottom,
  sSpacing, bPointList, nB,
  bMesh, mesh,
  prRad, prBath
 },
  dest = "cosine_general-verification-mesh-asymmetric.txt";
  ExportIfNotExists[dest,
    (* Top- and bottom-corner y-coordinates *)
    (*
      NOTE: "upper" and "lower" refer to the branch,
      not the relative vertical position of the two boundary curves
     *)
    xyUpper = xyTraAsymm["upper"];
    xyLower = xyTraAsymm["lower"];
    yTop = xyLower[[2]] @ DomainStart[xyLower];
    yBottom = xyUpper[[2]] @ DomainEnd[xyUpper];
    (* Make spacing of boundary points 1/50 of the height *)
    sSpacing = 1/50 (yTop - yBottom);
    (* Boundary points *)
    bPointList = Join[
      (* Radiation boundary: lower branch (physically higher) *)
      Table[
        xyLower[s] // Through
      , {s, UniformRange[
          DomainStart[xyLower],
          DomainEnd[xyLower],
          sSpacing
        ]
      }] // Most,
      (* Radiation boundary: upper branch (physically lower) *)
      Table[
        xyUpper[s] // Through
      , {s, UniformRange[
          DomainStart[xyUpper],
          DomainEnd[xyUpper],
          sSpacing
        ]
      }] // Most,
      (* Straight boundary x == \[Pi]/2 *)
      Table[
        {xStraight, y}
      , {y, UniformRange[yBottom, yTop, sSpacing]}] // Most
    ];
    nB = Length[bPointList];
    (* Build boundary element mesh *)
    bMesh = ToBoundaryMesh[
      "Coordinates" -> bPointList,
      "BoundaryElements" -> {
        LineElement[
          Table[{n, n + 1}, {n, nB}] // Mod[#, nB, 1] &
        ]
      }
    ];
    (* Build mesh *)
    mesh = ToElementMesh[bMesh,
      "ImproveBoundaryPosition" -> True
    ];
    (* Predicate functions for radiation and bath (straight) boundaries *)
    prRad = Function[{x, y}, x < xStraight // Evaluate];
    prBath = Function[{x, y}, x == xStraight // Evaluate];
    (* Export *)
    {mesh, prRad, prBath} // Compress
  ]
]


(* ::Subsubsubsection:: *)
(*Solve boundary value problem*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[
 {dest, source,
  mesh, prRad, prBath,
  a, tSol
 },
  dest = "cosine_general-verification-solution-asymmetric.txt";
  ExportIfNotExists[dest,
    (* Import mesh *)
    source = "cosine_general-verification-mesh-asymmetric.txt";
    {mesh, prRad, prBath} = Import[source] // Uncompress;
    (* Solve boundary value problem *)
    a = aAsymm;
    With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
      tSol = NDSolveValue[
        {
          (* Steady state heat equation *)
          -Laplacian[t[x, y], {x, y}] ==
            (* Radiation boundary condition *)
            NeumannValue[(-1 / a) t[x, y]^4, prRad[x, y]],
          (* Heat bath (straight edge) boundary condition *)
          DirichletCondition[t[x, y] == 1, prBath[x, y]]
        }, t, Element[{x, y}, mesh]
      ]
    ];
    tSol // Compress
  ]
]


(* ::Subsection:: *)
(*Italicise symbols*)


aIt = Italicise["A"];
bIt = Italicise["B"];
sIt = Italicise["s"];
xIt = Italicise["x"];
yIt = Italicise["y"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;


natStyle = Red;
flatStyle = Blue;
sharpStyle = Orange;


contStyle = LightGray;
straightStyle = Directive[Dashed, Magenta];
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;
inflStyle = Directive[Thick, Green];


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


simp0Style = Blue;


pointStyle = PointSize[Large];
startingPointStyle = Directive[Opacity[0.7], pointStyle];
glowStyle = Directive[Thick, Yellow, Opacity[0.7]];


inflDotStyle = Directive[Red, Opacity[0.7], pointStyle];


(* ::Subsection:: *)
(*Computations for figures*)


(* ::Subsubsection:: *)
(*Simple case (B = 1) traced boundaries*)


(* ::Subsubsubsection:: *)
(*Starting points along T = 0 contour*)


startXYSimpFigure =
  Module[{b, yMax, num, eps, yValues},
    (* Constants *)
    b = 1;
    yMax = 2;
    num = 16;
    eps = 10^-4;
    yValues = Subdivide[-yMax, yMax, num];
    (* Return starting points *)
    Table[
      {ArcCos[(1 - eps) Sech[y] / b], y}
    , {y, yValues}]
  ];


(* ::Subsubsubsection:: *)
(*Traced boundaries*)


xyTraSimpFigure =
  Module[{a, b, sMax, viTol, xInit, yInit},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 10^-6;
    Table[
      {xInit, yInit} = xyInit;
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyTraSystem[a, b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                vi[a, b][x[s], y[s]] < viTol,
                tKnown[b][x[s], y[s]] < 0
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax}
        ]
      ]
    , {xyInit, startXYSimpFigure}]
  ];


(* ::Subsubsection:: *)
(*Simple case (B = 1) traced boundaries, patched spiky*)


(*
  Corners (x_1, y_1), ..., (x_n, y_n) are (carefully) chosen,
  then the intersections between lower-branch(i) and upper-branch(i+1),
  for i == 1, ..., n - 1, chosen. Note that within each spike,
  the lower-branch curve is physically higher (in y).
*)


patchedCornerList = {
  {1, -0.5},
  {0.7, -0.25},
  {0.85, -0.13},
  {1.1, 0.2},
  {1.05, 0.3},
  {0.9, 0.8},
  {1.3, 1.23},
  Nothing
};


patchedCornerNum = Length[patchedCornerList];


patchedBoundaryUpperList =
  Module[{a, b, sMax, viTol, xInit, yInit},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 0;
    Table[
      {xInit, yInit} = xyInit;
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyTraSystem[a, b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                vi[a, b][x[s], y[s]] < viTol,
                tKnown[b][x[s], y[s]] < 0
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax}
        ]
      ]
    , {xyInit, patchedCornerList}]
  ];


patchedBoundaryLowerList =
  Module[{a, b, sMax, viTol, xInit, yInit},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 0;
    Table[
      {xInit, yInit} = xyInit;
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyTraSystem[a, b, False],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                vi[a, b][x[s], y[s]] < viTol,
                tKnown[b][x[s], y[s]] < 0
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax}
        ]
      ]
    , {xyInit, patchedCornerList}]
  ];


patchedIntersectionList =
  Table[
    SeekParametricIntersection[
      patchedBoundaryLowerList[[i]],
      patchedBoundaryUpperList[[i + 1]]
    ]
  , {i, patchedCornerNum - 1}];


{patchedIntersectionLowerList, patchedIntersectionUpperList} =
  Transpose @ patchedIntersectionList;
AppendTo[patchedIntersectionLowerList,
  DomainStart @ patchedBoundaryLowerList[[patchedCornerNum]]
];
PrependTo[patchedIntersectionUpperList,
  DomainEnd @ patchedBoundaryUpperList[[1]]
];


(* ::Subsubsection:: *)
(*Simple case (B = 1) traced boundaries, patched smooth*)


patchedBoundarySmooth =
  Module[{a, b, sMax, viTol},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 0;
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == x0Simp[a], y[0] == 0,
          WhenEvent[
            tKnown[b][x[s], y[s]] < 0,
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*General case (B arbitrary) asymmetric domain construction*)


(* ::Subsubsubsection:: *)
(*x_\[Flat] and x_\[Sharp]*)


xFlatAsymm = xFlat[aAsymm, bAsymm];
xSharpAsymm = xSharp[aAsymm, bAsymm];


(* ::Subsubsubsection:: *)
(*Inflection frontiers for the lower branch*)


asymmInflectionFrontierList =
  Module[{a, b, xyInitList},
    a = aAsymm;
    b = bAsymm;
    xyInitList = {
      {xInflAxis, 0},
      {xStraight, yInflStraight},
      Nothing
    };
    Table[
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            curContourSystem[a, b],
            x[0] == xyInit[[1]], y[0] == xyInit[[2]],
            WhenEvent[
              {
                x[s] < 0,
                x[s] > Pi,
                Abs @ y[s] > 2,
                vi[a, b][x[s], y[s]] < 0,
                False
              },
              "StopIntegration"
            ]
          },
          {x, y},
          {s, -6, 6}
          , NoExtrapolation
        ]
      ]
      , {xyInit, xyInitList}
    ]
  ];


(* ::Subsubsubsection:: *)
(*Convex portions of traced boundaries*)


asymmConvexPortionsList =
  Module[{a, b, xInitList, sConvexLower},
    a = aAsymm;
    b = bAsymm;
    xInitList = Subdivide[xInflStraight, xInflAxis, 4];
    xInitList = Append[xInitList, xSharpAsymm];
    Table[
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        (*
          NOTE: WhenEvent detection doesn't work at the initial point,
          hence the If hack.
        *)
        sConvexLower = If[xInit >= xInflAxis, 0, -6];
        NDSolveValue[
          {
            xyTraSystem[a, b],
            x[0] == xInit, y[0] == 0,
            WhenEvent[
              curTra[a, b, True][x[s], y[s]] < 0,
              "StopIntegration"
            ]
          },
          {x, y},
          {s, sConvexLower, 6}
          , NoExtrapolation
        ]
      ]
      , {xInit, xInitList}
    ]
  ];


(* ::Section:: *)
(*Known solution*)


(* ::Subsection:: *)
(*Algebra check for \[CapitalPhi]*)


(* ::Text:: *)
(*See (r5.17) (Page r5-2).*)


With[{a = \[FormalCapitalA], b = \[FormalCapitalB], x = \[FormalX], y = \[FormalY]},
  vi[a, b][x, y] ==
    b^2 (Sin[x]^2 + Sinh[y]^2)
    - (1 - b Cos[x] Cosh[y])^8 / a^2
    // FullSimplify
]


(* ::Subsection:: *)
(*Interactive visualiser for known solution*)


DynamicModule[
 {bInit, bMin, bMax,
  xMin, xMax, yMax
 },
  (* Values of B *)
  bInit = 1;
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2;
  yMax = 3;
  (* Plot *)
  Manipulate[
    Plot3D[
      tKnown[b][x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ClippingStyle -> unphysStyle,
      PlotLabel -> BoxedLabel[bIt == N[b]],
      PlotOptions[Axes] // Evaluate,
      PlotRange -> {0, Automatic}
    ]
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}]
]


(* ::Subsection:: *)
(*Interactive visualiser for unphysical domain (T < 0)*)


DynamicModule[
 {bInit, bMin, bMax,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinCont, xMaxCont, yMaxCont,
  numTo1, numBeyond1
 },
  (* Values of B *)
  bInit = 1;
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 3;
  (* Plot *)
  Manipulate[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[bIt == N[b]]
      ],
      (* Unphysical domain *)
      RegionPlot[
        tKnown[b][x, y] < 0,
        {x, xMin, xMax}, {y, -yMax, yMax},
        BoundaryStyle -> unphysStyle,
        PlotPoints -> 50,
        PlotStyle -> unphysStyle
      ],
      (* Known solution contours *)
      ContourPlot[
        tKnown[b][x, y],
        {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
        ContourLabels -> None,
        Contours -> numTo1 + numBeyond1,
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
      ],
      (* Straight contour *)
      Graphics @ {straightStyle,
        Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
      }
    ]
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}]
]


(* ::Section:: *)
(*Viable domain*)


(* ::Subsection:: *)
(*Interactive visualiser*)


DynamicModule[
 {aInit, aMin, aMax,
  bInit, bMin, bMax,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  xMinCont, xMaxCont, yMaxCont,
  numTo1, numBeyond1,
  resetLabel,
  critTermPoint
 },
  (* Values of A *)
  aInit = 0.2;
  aMin = 0.01;
  aMax = 2;
  (* Values of B *)
  bInit = 1.01 bNat[aInit];
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.1;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 3;
  (* Label for reset *)
  resetLabel[var_] := Row @ {"Reset ", var};
  (* Critical terminal point *)
  critTermPoint[x0_, style_] :=
    Graphics @ {Directive[style, pointStyle],
      Point @ {x0, 0}
    };
  (* Plot *)
  Manipulate[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      RegionPlot[tKnown[b][x, y] < 0,
        {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
        BoundaryStyle -> unphysStyle,
        PlotPoints -> 50,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours *)
      ContourPlot[
        tKnown[b][x, y],
        {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
        ContourLabels -> None,
        Contours -> numTo1 + numBeyond1,
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
      ],
      (* Straight contour *)
      Graphics @ {straightStyle,
        Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
      },
      (* Critical terminal points along x == 0 *)
      Which[
        (* B < 1 *)
        b < 1,
        {
          Graphics @ {Directive[Black, pointStyle],
            Point @ {{0, y0[a, b]}, {0, -y0[a, b]}}
          }
        },
        (* Otherwise don't care *)
        True,
        {}
      ],
      (* Critical terminal points along y == 0 *)
      Which[
        (* Two distinct terminal points, x_\[Flat] & x_\[Sharp] *)
        b > bNat[a],
        {
          critTermPoint[Re @ xFlat[a, b], flatStyle],
          critTermPoint[Re @ xSharp[a, b], sharpStyle]
          (* Re chops off small imaginary part in floating point arithmetic *)
        },
        (* One terminal point x_\[Natural] *)
        b == bNat[a],
        {
          critTermPoint[xNat[a], natStyle]
        },
        (* Zero critical terminal points *)
        True,
        {}
      ]
    ]
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}
  , Button[resetLabel[aIt], a = aInit]
  , Button[resetLabel[bIt], b = bInit]
  ]
]


(* ::Subsection:: *)
(*Interactive visualiser for \[CapitalPhi] along y = 0*)


DynamicModule[
 {aInit, aMin, aMax,
  bInit, bMin, bMax,
  xMin, xMax
 },
  (* Values of A *)
  aInit = 0.3;
  aMin = 0.01;
  aMax = 3;
  (* Values of B *)
  bInit = 1;
  bMin = 0.1;
  bMax = 3;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  (* Plot *)
  Manipulate[
    Show[
      EmptyFrame[{xMin, xMax}, {-0.1, 1},
        FrameLabel -> {xIt, "\[CapitalPhi]"[yIt == 0]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      Plot[{vi[a, b][x, 0], 0}, {x, xMin, xMax},
        Filling -> {1 -> {2}},
        PlotRange -> Full
      ],
      (* Location of straight contour *)
      Graphics @ {Directive[straightStyle, pointStyle],
        Point @ {xStraight, 0}
      }
    ]
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}]
]


(* ::Subsection:: *)
(*Animation for \[CapitalPhi] along y = 0, simple case (B = 1)*)


Module[
 {aMin, aMax, aStep,
  aValues, b,
  xMin, xMax
 },
  (* Values of A *)
  aMin = 1/10;
  aMax = 2;
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  (* Animation *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-0.1, 1},
        FrameLabel -> {xIt, "\[CapitalPhi]"[yIt == 0]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* \[CapitalPhi] *)
      Plot[{vi[a, b][x, 0], 0}, {x, xMin, xMax},
        Filling -> {1 -> {2}},
        PlotRange -> Full
      ],
      (* Location of straight contour *)
      Graphics @ {Directive[straightStyle, pointStyle],
        Point @ {xStraight, 0}
      },
      (* Non-trivial critical terminal point x_0 *)
      Graphics @ {Directive[Red, pointStyle = PointSize[Large]],
        Point @ {x0Simp[a], 0}
      }
    ]
  , {a, aValues}]
] // Ex["cosine_simple-phi.gif", gifOpts]


(* ::Subsection:: *)
(*Animation for viable domain, simple case (B = 1)*)


Module[
 {aMin, aMax, aStep,
  aValues, b,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  xMinCont, xMaxCont, yMaxCont,
  unphysicalDomain,
  numTo1, numBeyond1,
  genericContours, straightContour
 },
  (* Values of A *)
  aMin = 1/10;
  aMax = 2;
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.15;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (*
    NOTE: since B is fixed at 1,
    the unphysical domain and known solution contours
    do not need to be redrawn for different values of A.
   *)
  (* Unphysical domain *)
  unphysicalDomain =
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ];
  (* Known solution contours *)
  numTo1 = 5;
  numBeyond1 = 3;
  genericContours =
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
      ContourLabels -> None,
      Contours -> numTo1 + numBeyond1,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
    ];
  straightContour =
    Graphics @ {straightStyle,
      Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
    };
  (* Animation *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 180,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      unphysicalDomain,
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours (generic) *)
      genericContours,
      straightContour,
      (* Critical terminal point (x_0, 0) *)
      Graphics @ {Directive[simp0Style, pointStyle],
        Point @ {x0Simp[a], 0}
      },
      (* Known solution contour through (x_0, 0) *)
      ContourPlot[
        tKnown[b][x, y] == tKnown[b][x0Simp[a], 0],
        {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
        ContourLabels -> None,
        ContourStyle -> simp0Style,
        PlotRange -> {0, 1}
      ]
    ]
  , {a, aValues}]
] // Ex["cosine_simple-viable.gif", gifOpts]


(* ::Subsection:: *)
(*Animation for viable domain, general case (A = const, B varying)*)


Module[
 {a,
  bStep, bCrit, bMin, bMax, bValues,
  xMin, xMax, yMax,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  straightContour,
  critTermPoint
 },
  (* Value of A *)
  a = 3;
  (* Values of B *)
  bStep = 1/100;
  bCrit = bNat[a];
  bMin = Max[0, Floor[bCrit - 4 bStep, bStep]];
  bMax = Ceiling[bCrit + 8 bStep, bStep];
  bValues = Range[bMin, bMax, bStep];
  bValues = Append[bCrit] @ bValues;
  bValues = bValues // Sort[#, Less] &;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Plot ranges with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Straight contour (independent of A and B) *)
  straightContour =
    Graphics @ {straightStyle,
      Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
    };
  (* Critical terminal point *)
  critTermPoint[x0_, style_] :=
    Graphics @ {Directive[style, pointStyle],
      Point @ {x0, 0}
    };
  (* Animation *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 180,
        PlotLabel -> BoxedLabel[
          Row[
            {
              aIt == a,
              bIt == N[b, 2]
            } // If[b == bCrit,
              Style[#, natStyle] & /@ # &,
              Identity
            ],
            ","
          ],
          FrameStyle -> If[b == bCrit,
            natStyle,
            Automatic
          ]
        ]
      ],
      (* Unphysical domain *)
      RegionPlot[tKnown[b][x, y] < 0,
        {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
        BoundaryStyle -> unphysStyle,
        PlotPoints -> 50,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours *)
      ContourPlot[tKnown[b][x, y],
        {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
        ContourLabels -> None,
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, tKnown[b][xMax, 0]}
      ],
      straightContour,
      (* Critical terminal points along y == 0 *)
      Which[
        (* Two distinct terminal points, x_\[Flat] & x_\[Sharp] *)
        b > bCrit,
        {
          critTermPoint[Re @ xFlat[a, b], flatStyle],
          critTermPoint[Re @ xSharp[a, b], sharpStyle]
          (* Re chops off small imaginary part in floating point arithmetic *)
        },
        (* One terminal point x_\[Natural] *)
        b == bCrit,
        {
          critTermPoint[N @ xNat[a], natStyle]
          (* Doesn't work without applying N for some reason *)
        },
        (* Zero critical terminal points *)
        True,
        {}
      ]
    ]
  , {b, bValues}]
] // Ex["cosine_general-viable.gif", gifOpts]


(* ::Section:: *)
(*Simple case (B = 1) plots*)


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


Module[
 {b,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  xMinCont, xMaxCont, yMaxCont,
  unphysicalDomain,
  numTo1, numBeyond1, genericContours,
  xStraight, straightContour,
  idList
 },
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.1;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (*
    NOTE: since B is fixed at 1,
    the unphysical domain and known solution contours
    do not need to be redrawn for different values of A.
   *)
  (* Unphysical domain *)
  unphysicalDomain =
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ];
  (* Known solution contours (generic & straight) *)
  numTo1 = 5;
  numBeyond1 = 3;
  genericContours =
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
      ContourLabels -> None,
      Contours -> numTo1 + numBeyond1,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
    ];
  xStraight = Pi/2;
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
  };
  (* Group names *)
  idList = {"contour", "hyperbolic"};
  (* Plots for various A *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      unphysicalDomain,
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours (generic) *)
      genericContours,
      straightContour,
      (* Starting points *)
      ListPlot[
        Table[startXYSimp[a][id], {id, idList}],
        LabelingFunction -> Function @ Placed[
          #2[[2]],
          Center
        ],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
    // Ex @ FString[
      "cosine_simple-a-{aNamesSimp[a]}-traced-starting.pdf"
    ]
  , {a, aValuesSimp}]
]


(* ::Subsection:: *)
(*Traced boundaries*)


Module[
 {b,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  unphysicalDomain,
  xStraight, straightContour,
  idList
 },
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.1;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (*
    NOTE: since B is fixed at 1,
    the unphysical domain and known solution contours
    do not need to be redrawn for different values of A.
   *)
  (* Unphysical domain *)
  unphysicalDomain =
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ];
  (* Straight contour *)
  xStraight = Pi/2;
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxNon}, {xStraight, yMaxNon}}
  };
  (* Group names *)
  idList = {"contour", "hyperbolic"};
  (* Plots for various A *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      unphysicalDomain,
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Straight contour *)
      straightContour,
      (* Traced boundaries (general) *)
      Table[
        Table[
          ParametricPlot[
            xy[s]
              // Through
              // {#, {#[[1]], -#[[2]]}} &
              // Evaluate,
            {s, DomainStart[xy], DomainEnd[xy]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {xy, xyTraSimp[a][id]}]
      , {id, idList}],
      (* Traced boundaries (hyperbolic) *)
      Table[
        Table[
          ParametricPlot[
            xy[s]
              // Through
              // {#, {#[[1]], -#[[2]]}} &
              // Evaluate,
            {s, DomainStart[xy], DomainEnd[xy]},
            PlotStyle -> glowStyle
          ]
        , {xy, xyTraSimp[a][id]}]
      , {id, {"hyperbolic"}}]
    ]
    // Ex @ FString[
      "cosine_simple-a-{aNamesSimp[a]}-traced-full.pdf"
    ]
  , {a, aValuesSimp}]
]


(* ::Text:: *)
(*The terminal curve is practically a traced boundary.*)
(*For each A < 1 there will be some (x_i, y_i) where it inflects;*)
(*if x_i > \[Pi]/2, then we will have a convex domain*)
(*bounded by that traced boundary and by x = pi/2,*)
(*with thermal radiation and constant temperature along them respectively.*)
(*For a (contrived) practical situation, imagine rod at constant temperature*)
(*with a knob affixed at the end, and radiating into vacuum.*)


(* ::Subsection:: *)
(*Check boundary condition along terminal curve*)


With[{x = \[FormalX], y = \[FormalY], a = \[FormalCapitalA]},
  Module[
   {aValues, b, t, gradT,
    n, fluxL, fluxR, res,
    sMax,
    xTerm, yTerm
   },
    (* Values of A *)
    aValues = {1/10^3, 1/100, 1/10, 1/2, 3/4, 1, 5/4, 3/2};
    (* Value of B *)
    b = 1;
    (* Known solution T *)
    t = tKnown[b][x, y];
    (* Gradient of T *)
    gradT = Grad[t, {x, y}];
    (* Normal vector n *)
    n = -{y', -x'};
    (* Left hand side of boundary condition *)
    fluxL = n . gradT;
    (* Right hand side of boundary condition *)
    fluxR = -t^4 / a;
    (* Residual of boundary condition *)
    res["abs"] = fluxL - fluxR;
    res["rel"] = fluxR / fluxL - 1;
    (* Plots of residuals *)
    sMax = 5;
    Table[
      Plot[
        Table[
          (* Terminal curve x == x(s), y == y(s) *)
          {xTerm, yTerm} = xyTermSimp[a];
          (* Evaluate residual therealong *)
          res[type] /. {
            x' -> xTerm'[s],
            y' -> yTerm'[s],
            x -> xTerm[s],
            y -> yTerm[s]
          }
        , {a, aValues}] // Evaluate,
        {s, -sMax, sMax},
        AxesLabel -> {sIt, FString @ "Residual ({type})"},
        ImageSize -> 360,
        PlotLabel -> "Along terminal curve:",
        PlotLegends -> LineLegend[aValues // N, LegendLabel -> aIt],
        PlotOptions[Axes] // Evaluate
      ] // Ex @ FString[
        "cosine_simple-terminal-boundary-residual-{type}.pdf"
      ]
    , {type, {"abs", "rel"}}]
  ]
]


(* ::Subsection:: *)
(*Check boundary condition along candidate traced boundary*)


(* ::Text:: *)
(*(The candidate traced boundary is the traced boundary through (x_0, 0),*)
(*which is almost the same as the terminal curve.)*)


With[{x = \[FormalX], y = \[FormalY], a = \[FormalCapitalA]},
  Module[
   {aValues, b, t, gradT,
    n, fluxL, fluxR, res,
    sMax,
    xCand, yCand
   },
    (* Values of A *)
    aValues = {1/10^3, 1/100, 1/10, 1/2, 3/4, 1, 5/4, 3/2};
    (* Value of B *)
    b = 1;
    (* Known solution T *)
    t = tKnown[b][x, y];
    (* Gradient of T *)
    gradT = Grad[t, {x, y}];
    (* Normal vector n *)
    n = {y', -x'};
    (* Left hand side of boundary condition *)
    fluxL = n . gradT;
    (* Right hand side of boundary condition *)
    fluxR = -t^4 / a;
    (* Residual of boundary condition *)
    res["abs"] = fluxL - fluxR;
    res["rel"] = fluxR / fluxL - 1;
    (* Plots of residuals *)
    sMax = 5;
    Table[
      Plot[
        Table[
          (* Candidate boundary x == x(s), y == y(s) *)
          {xCand, yCand} = xyTraCandSimp[a];
          (* Evaluate residual therealong *)
          res[type] /. {
            x' -> xCand'[s],
            y' -> yCand'[s],
            x -> xCand[s],
            y -> yCand[s]
          } /. {
            s -> Abs[s]
          }
        , {a, aValues}] // Evaluate,
        {s, -sMax, sMax},
        AxesLabel -> {sIt, FString @ "Residual ({type})"},
        ImageSize -> 360,
        PlotLabel -> "Along terminal curve:",
        PlotLegends -> LineLegend[aValues // N, LegendLabel -> aIt],
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ] // Ex @ FString[
        "cosine_simple-candidate-boundary-residual-{type}.pdf"
      ]
    , {type, {"abs", "rel"}}]
  ]
]


(* ::Subsection:: *)
(*Curvature algebra*)


(* ::Text:: *)
(*See (r5.21) (Page r5-5).*)


With[{y = \[FormalY], a = \[FormalCapitalA]},
  Block[{$Assumptions = 0 < a <= 1 && 0 < y < ArcCosh[1/a]},
    curTra[a, 1][Pi/2, y] == (
      a^2 Cosh[y] / Sqrt[a^2 Cosh[y]^2 - 1] * (
        2 Sinh[y] - Cosh[y]^2 (a^2 Sinh[y] + 4 Sqrt[a^2 Cosh[y]^2 - 1])
      )
    )
    // FullSimplify
  ]
]


(* ::Text:: *)
(*See (r5.22) & (r5.23) (Page r5-6).*)


With[{y = \[FormalY], a = \[FormalCapitalA], s = \[FormalCapitalS]},
  2 Sinh[y] - Cosh[y]^2 (a^2 Sinh[y] + 4 Sqrt[a^2 Cosh[y]^2 - 1]) == (
    2 s - (1 + s^2) (a^2 s + 4 Sqrt[a^2 (1 + s^2) - 1])
  )
    /. {s -> Sinh[y]}
    // FullSimplify
]


(* ::Text:: *)
(*More compact version:*)


With[{y = \[FormalY], a = \[FormalCapitalA], c = \[FormalCapitalC], s = \[FormalCapitalS]},
  Block[{$Assumptions = 0 < a <= 1 && 0 < y < ArcCosh[1/a]},
    curTra[a, 1][Pi/2, y] == (
      a^2 c / Sqrt[a^2 c^2 - 1] * (
        2 s - (1 + s^2) (a^2 s + 4 Sqrt[a^2 (1 + s^2) - 1])
      )
    )
      /. {c -> Cosh[y], s -> Sinh[y]}
      // FullSimplify
  ]
]


(* ::Subsection:: *)
(*Horizontal coordinate x at critical y = y(A)*)


(* ::Text:: *)
(*A plot of x(y(A)) - \[Pi]/2, whose zero A = A_i was sought.*)
(*See (r5.25) (Page r5-6).*)


Module[
 {aMin, aMax, aInfl
 },
  (* Plot range *)
  aMin = 0.15;
  aMax = 1;
  (* Value of A_i *)
  aInfl = aInflSimp;
  (* Plot *)
  Show[
    (* A vs x(y(A)) - \[Pi]/2 *)
    Plot[
      xTraCandSimp[a] @ yCurCritSimp[a] - Pi / 2 // Evaluate,
      {a, aMin, aMax},
      AxesLabel -> {aIt, xIt @ yIt[aIt] - Pi / 2},
      ImageSize -> 360,
      PlotOptions[Axes] // Evaluate
    ],
    (* A_i *)
    Graphics @ {inflDotStyle,
      Point @ {aInfl, 0}
    },
    Graphics @ Text[
      Subscript[aIt, "i"] == N[aInfl]
        // textStyle,
      {aInfl, 0},
      {0, -1.3}
    ]
  ]
] // Ex["cosine_simple-candidate-x-minus-half-pi.pdf"]


(* ::Subsection:: *)
(*Approximate A_i (more curvature algebra)*)


(* ::Text:: *)
(*Check (r5.28) (Page r5-6).*)


With[{a = \[FormalCapitalA]},
  vi[a, 1][Pi/2, ArcSech[a]] == 0
]


(* ::Subsubsection:: *)
(*Curvature of \[CapitalPhi] = 0 at (x, y) = (\[Pi]/2, arcsech(A))*)


(* ::Text:: *)
(*Recall that the terminal curve \[CapitalPhi] = 0*)
(*is almost indistinguishable from the candidate boundary.*)
(*See (r5.29) (Page r5-7).*)


With[{x = \[FormalX], y = \[FormalY], a = \[FormalCapitalA]},
  Module[{div, grad, cur},
    (* Abbreviations *)
    div = Div[#, {x, y}] &;
    grad = Grad[#, {x, y}] &;
    (* Curvature *)
    cur = div @ Normalize @ grad @ vi[a, 1][x, y];
    (* Evaluate at crossing *)
    cur /. {x -> Pi/2, y -> ArcSech[a]}
      // FullSimplify[#, 0 < a < 1] &
  ]
]


Module[
 {aMin, aMax, aInfl
 },
  (* Plot range *)
  aMin = 0;
  aMax = 1;
  (* Value of A_i(approx) *)
  aInfl = aInflSimpApprox;
  (* Plot *)
  Show[
    (* A vs \[Kappa] *)
    Plot[
      a (a^6 - a^4 + 44 a^2 - 28) / (16 + a^2 - a^4)^(3/2),
      {a, aMin, aMax},
      AxesLabel -> {aIt, "\[Kappa]"},
      ImageSize -> 360,
      PlotLabel -> Column @ {
        "Curvature of terminal curve",
        Row @ {"at ", ""[xIt == Pi/2, yIt == ArcSech[aIt]]}
      },
      PlotOptions[Axes] // Evaluate
    ],
    (* A_i *)
    Graphics @ {inflDotStyle,
      Point @ {aInfl, 0}
    },
    Graphics @ Text[
      Subscript[aIt, "i"] == N[aInfl]
        // textStyle,
      {aInfl, 0},
      {0, -1.3}
    ]
  ]
] // Ex["cosine_simple-terminal-curvature-x-at-half-pi.pdf"]


(* ::Subsection:: *)
(*Convex domains without known solution*)


Module[
 {b,
  xRad, yEnd,
  xMin, xMax, yMax
 },
  Table[
    (* Value of B *)
    b = 1;
    (* Radiation boundary x == x(y) for convex domain *)
    xRad = xTraCandSimp[a, True];
    yEnd = DomainEnd[xRad];
    xRad = Function[{y}, xRad[Abs @ y] // Evaluate];
    (* Plot range *)
    xMin = xStraight - 0.6 yEnd;
    xMax = xStraight + 0.6 yEnd;
    yMax = yEnd;
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[aIt == N[a]]
      ],
      (* Convex domain *)
      ParametricPlot[
        {
          {xRad[y], y},
          {xStraight, y}
        }, {y, -yEnd, yEnd},
        PlotStyle -> convexStyle
      ]
    ] // Ex @ FString[
      "cosine_simple-traced-{aNamesSimpConvex[a]}.pdf"
    ]
  , {a, aValuesSimpConvex}]
]


(* ::Subsection:: *)
(*Convex domain aspect ratio*)


Module[{aMin, aMax},
  (* Plot range *)
  aMin = aInflSimp;
  aMax = 1;
  (* Plot *)
  Show[
    (* A vs \[Kappa] *)
    Plot[
      aspectRatioSimp[a],
      {a, aMin, aMax},
      AxesLabel -> {aIt, Null},
      ImageSize -> 360,
      PlotLabel -> "Aspect ratio of convex domains",
      PlotRange -> {0, Automatic},
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["cosine_simple-traced-convex-aspect-ratio.pdf"]


(* ::Subsection:: *)
(*Non-convex lens-like candidate self-radiation*)


(* ::Subsubsection:: *)
(*A vs R_crude*)


Module[
  {
    aMin, aMax,
    rMin, rMax,
    dummyForTrailingCommas
  },
  aMin = 0;
  aMax = aInflSimp;
  rMin = 10^-8;
  rMax = 20;
  Show[
    ListLogPlot[
      candidateBoundaryRatioTable
      , AxesLabel -> {aIt, Subscript[Italicise["R"], "crude"]}
      , Joined -> True
      , PlotRange -> {rMin, rMax}
      , PlotRangeClipping -> False
      , PlotOptions[Axes] // Evaluate
    ],
    Plot[
      rPracCandSimp // Log
      , {a, aMin, aMax}
      , PlotStyle -> Red
    ],
    Graphics @ {
      Directive[Red, Dashed],
      Line @ {{aPracCandSimp, rPracCandSimp // Log}, {aPracCandSimp, rMin // Log}},
      Directive[Red, PointSize[Large]],
      Point @ {aPracCandSimp, rMin // Log},
      Red,
      Text[
        Subscript[Italicise["R"], "p"] == PercentForm @ N[rPracCandSimp]
        , {aMin, rPracCandSimp // Log}
        , {-1.2, -1.2}
      ],
      Text[
        Subscript[aIt, "p"] == aPracCandSimp
        , {aPracCandSimp, rMin // Log}
        , {-1.2, -1}
      ],
      {}
    },
    {}
  ]
] // Ex["cosine_simple-candidate-boundary-ratio-bound-ultra.pdf"]


(* ::Subsubsection:: *)
(*Practical non-convex lenses*)


Module[
  {
    b,
    aValues,
    aMin, xRadAMin, yEndAMin,
    eps,
    yMaxFrame, xMinFrame, xMaxFrame,
    imageSize,
    plotList,
    xRad, yEnd,
    xRadEvenExtension,
    yInfl, xInfl,
    dummyForTrailingCommas
   },
  (* Value of B *)
  b = 1;
  (* Values of A *)
  aValues = Subdivide[aPracCandSimp, aInflSimp, 4];
  (* Plot range *)
  aMin = Min[aValues];
  xRadAMin = xTraCandSimp[aMin, True];
  yEndAMin = DomainEnd[xRadAMin];
  eps = 0.02;
  yMaxFrame = yEndAMin;
  xMinFrame = (1 - eps) x0Simp[aMin];
  xMaxFrame = (1 + eps) xStraight;
  (* List of plots *)
  imageSize = 96;
  plotList =
    Table[
      (* Radiation boundary x == x(y) for convex domain *)
      xRad = xTraCandSimp[a, True];
      yEnd = DomainEnd[xRad];
      xRadEvenExtension = Function[{y}, xRad[Abs @ y] // Evaluate];
      (* Point of inflection *)
      Which[
        a < aInflSimp,
          yInfl = SeekRoot[
            curTra[a, b][xRad[#], #] &,
            {DomainStart[xRad], DomainEnd[xRad]}
          ],
        a == aInflSimp,
          yInfl = DomainEnd[xRad],
        True,
          yInfl = Indeterminate
      ];
      xInfl = xRad[yInfl];
      (* Plot *)
      Show[
        EmptyFrame[
          {xMinFrame, xMaxFrame}, {-yMaxFrame, yMaxFrame}
          , Frame -> None
          , ImageSize -> 4 imageSize
            (* NOTE: GraphicsRow makes the plot smaller *)
          , PlotRangePadding -> None
        ],
        (* Convex domain *)
        ParametricPlot[
          {
            {xRadEvenExtension[y], y},
            {xStraight, y}
          }
          , {y, -yEnd, yEnd}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle /@ {"Traced", "Contour"}
        ],
        (* Point of inflection *)
        If[NumericQ[yInfl],
          Graphics @ {GeneralStyle["Point"],
            Point @ {{xInfl, yInfl}, {xInfl, -yInfl}}
          },
          {}
        ],
        {}
      ]
      , {a, aValues}
    ];
  (* Final figure *)
  Grid @ {{
    Column[
      {
        Row @ {
          GraphicsRow[plotList, Spacings -> {3.2 imageSize, 0}],
          Row @ {Graphics[ImageSize -> 0.5 imageSize]}
        },
        Nothing
      }
    , Center
    , Spacings -> 0
    ],
    Nothing
  }}
]


(* ::Section:: *)
(*General case (B arbitrary) plots*)


(* ::Subsection:: *)
(*Critical terminal points along y = 0*)


Module[
 {a,
  bCrit, xCrit,
  bMin, bMax,
  xMin, xMax
 },
  (* Value of A *)
  a = 3;
  (* Critical values B_\[Natural] and x_\[Natural] *)
  (* (N needed otherwise Epilog doesn't work) *)
  bCrit = bNat[a] // N;
  xCrit = xNat[a] // N;
  (* B range *)
  bMin = 0;
  bMax = 5;
  (* x range *)
  xMin = 0;
  xMax = 1.1 xSharp[a, #] & [
    With[{b = \[FormalCapitalB]},
      b /. Last @ Solve[
        D[xSharp[a, b], b] == 0,
        b
      ]
    ]
  ];
  (* Plot *)
  Plot[
    (* x_\[Sharp] & x_\[Flat]*)
    {
      xSharp[a, b],
      xFlat[a, b],
      Piecewise @ {
        {xMin, b < bCrit},
        {Indeterminate, True}
      }
    },
    {b, bMin + 10^-6, bMax},
    AxesLabel -> {bIt, xIt},
    Epilog -> {
      (* x == x_\[Natural], B == B_\[Natural] *)
      Directive[pointStyle, natStyle],
      Point @ {bCrit, xCrit},
      Text[
        Column[
          textStyle /@ {
            Subscript[xIt, "nat"] == xCrit,
            Subscript[bIt, "nat"] == bCrit
          }
        ],
        {bCrit, xCrit},
        {-1.2, 0}
      ],
    },
    Filling -> {
      1 -> Top,
      2 -> Bottom,
      3 -> Top
    },
    FillingStyle -> nonStyle,
    PlotLabel -> BoxedLabel[aIt == a],
    PlotLegends -> (
      xIt == Subscript[xIt, #] & /@
        {"sharp", "flat"}
    ),
    PlotRange -> {xMin, xMax},
    PlotStyle -> {flatStyle, sharpStyle, None},
    PlotOptions[Axes] // Evaluate
  ] // PrettyString[
    "flat" -> "\[Flat]",
    "nat" -> "\[Natural]",
    "sharp" -> "\[Sharp]"
  ] // Ex @ FString[
    "cosine_general-critical-a-{ToName[a]}.pdf"
  ]
]


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* (This is somewhat slow, as 10 files are being generated in one go.) *)
Module[
 {xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  labelFun,
  regime, b, idList
 },
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 4;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 240,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* Labelling function *)
  labelFun = Function @ Placed[
    #2[[2]],
    Center
  ];
  (* For each value of A *)
  Table[
    (* For each regime *)
    Table[
      (* Value of B *)
      b = bValueGen[regime][a];
      (* Starting point groups *)
      idList = idListGen[regime];
      (* Plot *)
      Show[
        emptyFrame[a, b],
        unphysicalDomain[b],
        nonViableDomain[a, b],
        generalContours[b],
        straightContour,
        ListPlot[
          Table[
            startXYGen[a][regime][id]
          , {id, idList}],
          LabelingFunction -> labelFun,
          PlotLegends -> idList,
          PlotStyle -> startingPointStyle
        ]
      ] // Ex @ FString[
        "cosine_general-a-{ToName[a]}-{regime}-traced-starting.pdf"
      ]
    , {regime, regimeListGen}]
  , {a, aValuesGen}]
]


(* ::Subsection:: *)
(*Traced boundaries*)


(* (This is somewhat slow, as 10 files are being generated in one go.) *)
Module[
 {xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  regime, b, idList
 },
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 4;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 240,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* For each value of A *)
  Table[
    (* For each regime *)
    Table[
      (* Value of B *)
      b = bValueGen[regime][a];
      (* Starting point groups *)
      idList = idListGen[regime];
      (* Plot *)
      Show[
        emptyFrame[a, b],
        unphysicalDomain[b],
        nonViableDomain[a, b],
        generalContours[b],
        straightContour,
        (* Traced boundaries (general) *)
        Table[
          Table[
            ParametricPlot[
              xy[s]
                // Through
                // {#, {#[[1]], -#[[2]]}} &
                // Evaluate,
              {s, DomainStart[xy], DomainEnd[xy]},
              PlotStyle -> {upperStyle, lowerStyle}
            ]
          , {xy, xyTraGen[a][regime][id]}]
        , {id, idList}]
        (* Not bothering with hyperbolic traced boundaries *)
      ] // Ex @ FString[
        "cosine_general-a-{ToName[a]}-{regime}-traced-full.pdf"
      ]
    , {regime, regimeListGen}]
  , {a, aValuesGen}]
]


(* ::Subsection:: *)
(*Traced boundaries from x-axis*)


(* ::Subsubsection:: *)
(*Visualiser (poor man's Manipulate)*)


Module[
 {a, b,
  xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  lineOfSymmetry,
  sMax,
  xInitMin, xInitMax, xInitList,
  xyList,
  tracedBoundaries,
  xyMid, sInflList,
  xInitInfl, yInitInfl,
  xyInflList,
  inflectionFrontiers
 },
  (* Values of A and B (to be set manually) *)
  a = 12;
  b = 1.05 bNat[a];
  (*
    ------------------------------------------------
    Barely an infinite family of convex domains
    ------------------------------------------------
    a = 7.5;
    b = 1.05 bNat[a];
    ------------------------------------------------
    Probably only one convex domain (lens-like)
    ------------------------------------------------
    a = 7;
    b = 1.1 bNat[a];
    ------------------------------------------------
    Weird inflection frontier
    ------------------------------------------------
    a = 7.5;
    b = 1.3 bNat[a];
    ------------------------------------------------
    Slightly asymmetric convex domains
    ------------------------------------------------
    a = 12;
    b = 1.05 bNat[a];
   *)
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* Line of symmetry (y == 0) *)
  lineOfSymmetry[a_, b_] := Graphics @ {Directive[Thin, Orange],
    Line @ N @ {{xFlat[a, b], 0}, {xSharp[a, b], 0}}
  };
  (* Range for s *)
  sMax = 6;
  (* Initial x values along axis *)
  xInitMin = xFlat[a, b];
  xInitMax = xSharp[a, b];
  xInitList = Subdivide[xInitMin, xInitMax, 4];
  (* Exclude x_\[Flat] == 0 for B = 1 *)
  If[b == 1, xInitList = xInitList // Rest];
  (* Compute traced boundaries *)
  xyList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInit, y[0] == 0,
          WhenEvent[
            {
              vi[a, b][x[s], y[s]] < 0,
              tKnown[b][x[s], y[s]] < 0
            },
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  , {xInit, xInitList}];
  tracedBoundaries =
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> {upperStyle, lowerStyle}
      ]
    , {xy, xyList}];
  (* Compute inflection frontiers for the lower branch *)
  xyMid = xyList[[3]];
  sInflList =
    Table[
      SeekFirstRootBisection[
        curTra[a, b] @@ Through @ xyMid[#] &,
        {0, end[xyMid]}
      ]
    , {end, {DomainStart, DomainEnd}}];
  xyInflList =
    Table[
      {xInitInfl, yInitInfl} = xyMid[sInit] // Through;
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            curContourSystem[a, b],
            x[0] == xInitInfl, y[0] == yInitInfl,
            WhenEvent[
              {
                x[s] < xMinMar,
                x[s] > xMaxMar,
                Abs @ y[s] > yMaxMar,
                vi[a, b][x[s], y[s]] < 0
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ]
    , {sInit, sInflList}];
  inflectionFrontiers[
    opt : OptionsPattern @ {"Mirror" -> False}
  ] :=
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // If[OptionValue["Mirror"],
            {#, {#[[1]], -#[[2]]}} &,
            Identity
          ]
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> inflStyle
      ]
    , {xy, xyInflList}];
  (* Plot *)
  Show[
    emptyFrame[a, b],
    unphysicalDomain[b],
    nonViableDomain[a, b],
    generalContours[b],
    straightContour,
    lineOfSymmetry[a, b],
    tracedBoundaries,
    inflectionFrontiers[]
  ]
]


(* ::Subsubsection:: *)
(*Animation*)


Module[
 {bMin, bMax, bValues,
  xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  sMax,
  xInitMin, xInitMax, xInitList,
  xyList
 },
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2.5;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 200,
      PlotLabel -> BoxedLabel[
        Column[
          {aIt == a, bIt == N[b, 2]},
          Alignment -> Center
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* Range for s *)
  sMax = 6;
  (* For various values of A *)
  Table[
    (* Values of B *)
    bMin = bNat[a];
    bMax = 1;
    bValues = Subdivide[bMin, bMax, 20];
    (* B-manipulation *)
    Table[
      (* Initial x values along axis *)
      xInitMin = xFlat[a, b];
      xInitMax = xSharp[a, b];
      xInitList = Subdivide[xInitMin, xInitMax, 4];
      (* Exclude x_\[Flat] == 0 for B = 1 *)
      If[b == 1, xInitList = xInitList // Rest];
      (* Compute traced boundaries *)
      xyList = Table[
        With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
          NDSolveValue[
            {
              xyTraSystem[a, b],
              x[0] == xInit, y[0] == 0,
              WhenEvent[
                {
                  vi[a, b][x[s], y[s]] < 0,
                  tKnown[b][x[s], y[s]] < 0
                },
                "StopIntegration"
              ]
            }, {x, y}, {s, -sMax, sMax},
            NoExtrapolation
          ]
        ]
      , {xInit, xInitList}];
      (* Plot *)
      Show[
        emptyFrame[a, b],
        unphysicalDomain[b],
        nonViableDomain[a, b],
        generalContours[b],
        straightContour,
        (* Traced boundaries *)
        Table[
          ParametricPlot[
            xy[s]
              // Through
              // {#, {#[[1]], -#[[2]]}} &
              // Evaluate,
            {s, DomainStart[xy], DomainEnd[xy]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {xy, xyList}]
      ]
    , {b, bValues}]
      // Ex[
        FString["cosine_general-a-{a}-traced-from-axis.gif"],
        gifOpts
      ]
  , {a, aValuesGen}]
]


(* ::Subsection:: *)
(*Curvature along x-axis*)


(* ::Subsubsection:: *)
(*Interactive visualiser*)


DynamicModule[
 {aInit, aMin, aMax,
  bInit, bMin, bMax,
  xMin, xMax,
  curMax, curMaxMar,
  resetLabel,
  critTermPoint
 },
  (* Values of A *)
  aInit = 12;
  aMin = 1;
  aMax = 20;
  (* Values of B *)
  bInit = 1.05 bNat[aInit];
  bMin = 0;
  bMax = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  curMax = 1;
  curMaxMar = 1.2 curMax;
  (* Label for reset *)
  resetLabel[var_] := Row @ {"Reset ", var};
  (* Critical terminal point *)
  critTermPoint[x0_, style_] :=
    Graphics @ {Directive[style, pointStyle],
      Point @ {x0, 0}
    };
  (* Plot *)
  Manipulate[
    Show[
      Plot[
        curTra[a, b][x, 0], {x, xMin, xMax},
        AxesLabel -> {Automatic, "\[Kappa]"},
        PlotRange -> {-curMax, curMax},
        PlotOptions[Axes] // Evaluate
      ],
      (* Straight contour *)
      Graphics @ {straightStyle,
        Line @ {{xStraight, -curMaxMar}, {xStraight, curMaxMar}}
      },
      (* Critical terminal points along y == 0 *)
      Which[
        (* Two distinct terminal points, x_\[Flat] & x_\[Sharp] *)
        b > bNat[a],
        {
          critTermPoint[Re @ xFlat[a, b], flatStyle],
          critTermPoint[Re @ xSharp[a, b], sharpStyle]
          (* Re chops off small imaginary part in floating point arithmetic *)
        },
        (* One terminal point x_\[Natural] *)
        b == bNat[a],
        {
          critTermPoint[xNat[a], natStyle]
        },
        (* Zero critical terminal points *)
        True,
        {}
      ]
    ]
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}
  , Button[resetLabel[aIt], a = aInit]
  , Button[resetLabel[bIt], b = bInit]
  ]
]


(* ::Subsection:: *)
(*Constructing an asymmetric domain*)


Module[
 {a, b,
  includeYReflection,
  xFl, xSh, sMax,
  xMin, xMax, yMax,
  mar, xMinMar, xMaxMar, yMaxMar,
  xInitGenericList, xyGenericList,
  xyInitInflectionList, xyInflectionList,
  xInitConvexList, sConvexLower, xyConvexList,
  emptyFrame,
  nonViableDomain,
  straightContour,
  lineOfSymmetry,
  inflectionFrontiers
 },
  (*
    ------------------------------------------------
    Definitions
    ------------------------------------------------
   *)
  (* A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Include reflection in y (across x-axis) *)
  includeYReflection = {#, # * {1, -1}} &;
  (* Useful constants *)
  xFl = xFlat[a, b];
  xSh = xSharp[a, b];
  sMax = 6;
  (* Plot range *)
  xMin = Way[0, xFl] // N;
  xMax = 2;
  yMax = 2/3 (xMax - xMin);
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Generic traced boundaries from axis *)
  xInitGenericList = Subdivide[xFl, xSh, 4];
  xyGenericList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInit, y[0] == 0,
          WhenEvent[
            {
              x[s] < xMinMar,
              x[s] > xMaxMar,
              Abs @ y[s] > xMaxMar,
              vi[a, b][x[s], y[s]] < 0,
              tKnown[b][x[s], y[s]] < 0
            },
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  , {xInit, xInitGenericList}];
  (* Inflection frontiers for the lower branch *)
  xyInitInflectionList = {
    {xInflAxis, 0},
    {xStraight, yInflStraight}
  };
  xyInflectionList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          curContourSystem[a, b],
          x[0] == xyInit[[1]], y[0] == xyInit[[2]],
          WhenEvent[
            {
              x[s] < xMinMar,
              x[s] > xMaxMar,
              Abs @ y[s] > xMaxMar,
              vi[a, b][x[s], y[s]] < 0
            },
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  , {xyInit, xyInitInflectionList}];
  (* Convex portions of traced boundaries *)
  xInitConvexList = Subdivide[xInflStraight, xInflAxis, 4];
  xyConvexList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      (*
        NOTE: WhenEvent detection doesn't work at the initial point,
        hence the If hack.
       *)
      sConvexLower = If[xInit == Last[xInitConvexList], 0, -sMax];
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInit, y[0] == 0,
          WhenEvent[
            curTra[a, b, True][x[s], y[s]] < 0,
            "StopIntegration"
          ]
        }, {x, y}, {s, sConvexLower, sMax},
        NoExtrapolation
      ]
    ]
  , {xInit, xInitConvexList}];
  (*
    ------------------------------------------------
    Re-used plots
    ------------------------------------------------
   *)
  (* Empty frame *)
  emptyFrame =
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ];
  (* Non-viable domain *)
  nonViableDomain =
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ];
  (* Straight contour *)
  straightContour =
    Graphics @ {straightStyle,
      Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
    };
  (* Line of symmetry (y == 0) *)
  lineOfSymmetry =
    Graphics @ {Directive[Thin, Orange],
      Line @ N @ {{xFlat[a, b], 0}, {xSharp[a, b], 0}}
    };
  (* Inflection frontiers for the lower branch *)
  inflectionFrontiers[
    opt : OptionsPattern @ {"Mirror" -> False}
  ] :=
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // If[OptionValue["Mirror"], includeYReflection, Identity]
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> inflStyle
      ]
    , {xy, xyInflectionList}];
  (*
    ------------------------------------------------
    Plots
    ------------------------------------------------
   *)
  {
    (* Generic traced boundaries from axis *)
    Show[
      emptyFrame,
      nonViableDomain,
      straightContour,
      lineOfSymmetry,
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // includeYReflection
            // Evaluate,
          {s, DomainStart[xy], DomainEnd[xy]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {xy, xyGenericList}],
      inflectionFrontiers[]
    ]
      // Ex["cosine_general-asymmetric-traced-generic.pdf"],
    (* Convex portions of traced boundaries *)
    Show[
      emptyFrame,
      nonViableDomain,
      straightContour,
      lineOfSymmetry,
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // includeYReflection
            // Evaluate,
          {s, DomainStart[xy], DomainEnd[xy]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {xy, xyConvexList}],
      inflectionFrontiers["Mirror" -> True]
    ]
      // Ex["cosine_general-asymmetric-traced-convex_portions.pdf"],
    (* Convex asymmetric domain boundaries *)
    Show[
      emptyFrame,
      nonViableDomain,
      straightContour,
      lineOfSymmetry,
      ParametricPlot[
        Table[
          xyTraAsymm[id][s] // Through
        , {id, {"upper", "lower"}}] // Evaluate,
        {s, -sMax, sMax},
        PlotStyle -> convexStyle
      ],
      inflectionFrontiers["Mirror" -> True]
    ]
      // Ex["cosine_general-asymmetric-traced-convex_domain.pdf"]
  }
]


(* ::Section:: *)
(*Numerical verification (B = 1) plots*)


(* ::Subsection:: *)
(*Finite element mesh*)


Table[
  Module[
   {source,
    a, mesh, prRad, prBath,
    bCoords, bCoordsRad, bCoordsBath,
    dest
   },
    (* Import mesh *)
    source = FString[
      "cosine_simple-verification-mesh-{aNamesSimpConvex[a]}.txt"
    ];
    {a, mesh, prRad, prBath} = Import[source] // Uncompress;
    (* Boundary coordinates *)
    bCoords = Part[
      mesh["Coordinates"],
      List @@ First @ mesh["BoundaryElements"] // Flatten // DeleteDuplicates
    ];
    bCoordsRad = Select[bCoords, prRad @@ # &];
    bCoordsBath = Select[bCoords, prBath @@ # &];
    (* Export plot *)
    dest = StringReplace[source, ".txt" -> ".pdf"];
    Show[
      mesh["Wireframe"],
      ListPlot[{bCoordsRad, bCoordsBath},
        PlotStyle -> (
          Directive[#, pointStyle, Opacity[0.7]] &
            /@ {Blue, Red}
        )
      ],
      ImageSize -> 240
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsection:: *)
(*Numerical solution*)


Table[
  Module[{source, tSol, mesh, dest},
    (* Import solution *)
    source = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Plot *)
    dest = source // StringReplace[".txt" -> ".png"];
    With[{x = \[FormalX], y = \[FormalY]},
      Plot3D[tSol[x, y], Element[{x, y}, mesh],
        AxesLabel -> Italicise /@ {"x", "y", "T"},
        PlotLabel -> Column[
          {"Numerical solution", aIt == N[a]},
          Center
        ],
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsection:: *)
(*Relative error (3D)*)


Table[
  Module[{source, tSol, mesh, dest},
    (* Import solution *)
    source = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Plot *)
    dest = FString[
      "cosine_simple-verification-rel_error-3d-{aNamesSimpConvex[a]}.png"
    ];
    With[{x = \[FormalX], y = \[FormalY]},
      Plot3D[
        tSol[x, y] / tKnown[1][x, y] - 1, Element[{x, y}, mesh],
        PlotLabel -> Column[
          {"Relative error of numerical solution", aIt == N[a]},
          Center
        ],
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsection:: *)
(*Relative error (2D)*)


Table[
  Module[{source, tSol, mesh, dest},
    (* Import solution *)
    source = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Plot *)
    (* Note the aspect ratio *)
    dest = FString[
      "cosine_simple-verification-rel_error-2d-{aNamesSimpConvex[a]}.png"
    ];
    With[{x = \[FormalX], y = \[FormalY]},
      Show[
        DensityPlot[
          tSol[x, y] / tKnown[1][x, y] - 1, Element[{x, y}, mesh],
          ColorFunction -> "Rainbow",
          ImageSize -> 320,
          PlotLabel -> Column[
            {"Relative error of numerical solution", aIt == N[a]},
            Center
          ],
          PlotRange -> Full,
          PlotLegends -> Automatic,
          PlotOptions[Axes] // Evaluate
        ],
        mesh @ "Wireframe"[
          "MeshElementStyle" -> EdgeForm[Opacity[0.1]]
        ]
      ]
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Section:: *)
(*Numerical verification (B arbitrary) plots*)


(* ::Subsection:: *)
(*Finite element mesh*)


Module[
 {source,
  mesh, prRad, prBath,
  bCoords, bCoordsRad, bCoordsBath,
  dest
 },
  (* Import mesh *)
  source = "cosine_general-verification-mesh-asymmetric.txt";
  {mesh, prRad, prBath} = Import[source] // Uncompress;
  (* Boundary coordinates *)
  bCoords = Part[
    mesh["Coordinates"],
    List @@ First @ mesh["BoundaryElements"] // Flatten // DeleteDuplicates
  ];
  bCoordsRad = Select[bCoords, prRad @@ # &];
  bCoordsBath = Select[bCoords, prBath @@ # &];
  (* Export plot *)
  dest = StringReplace[source, ".txt" -> ".pdf"];
  Show[
    mesh["Wireframe"],
    ListPlot[{bCoordsRad, bCoordsBath},
      PlotStyle -> (
        Directive[#, pointStyle, Opacity[0.7]] &
          /@ {Blue, Red}
      )
    ],
    ImageSize -> 240
  ] // Ex[dest]
]


(* ::Subsection:: *)
(*Numerical solution*)


Module[{a, b, source, tSol, mesh, dest},
  a = aAsymm;
  b = bAsymm;
  (* Import solution *)
  source = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  dest = source // StringReplace[".txt" -> ".png"];
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y], Element[{x, y}, mesh],
      AxesLabel -> Italicise /@ {"x", "y", "T"},
      PlotLabel -> Column[
        {
          "Numerical solution",
          Row[{aIt == a, bIt == N[b]}, ","]
        },
        Center
      ],
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ] // Ex[dest]
]


(* ::Subsection:: *)
(*Relative error (3D)*)


Module[{a, b, source, tSol, mesh, dest},
  a = aAsymm;
  b = bAsymm;
  (* Import solution *)
  source = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  dest = "cosine_general-verification-rel_error-3d-asymmetric.png";
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tKnown[b][x, y] - 1, Element[{x, y}, mesh],
      AxesLabel -> Italicise /@ {"x", "y", "T"},
      PlotLabel -> Column[
        {
          "Rel. error of numerical solution",
          Row[{aIt == a, bIt == N[b]}, ","]
        },
        Center
      ],
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ] // Ex[dest]
]


(* ::Subsection:: *)
(*Relative error (2D)*)


Module[{a, b, source, tSol, mesh, dest},
  a = aAsymm;
  b = bAsymm;
  (* Import solution *)
  source = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  dest = "cosine_general-verification-rel_error-2d-asymmetric.png";
  With[{x = \[FormalX], y = \[FormalY]},
    Show[
      DensityPlot[tSol[x, y] / tKnown[b][x, y] - 1, Element[{x, y}, mesh],
        ColorFunction -> "Rainbow",
        PlotLabel -> Column[
          {
            "Rel. error of numerical solution",
            Row[{aIt == a, bIt == N[b]}, ","]
          },
          Center
        ],
        PlotRange -> Full,
        PlotLegends -> Automatic,
        PlotOptions[Frame] // Evaluate
      ],
      mesh["Wireframe"],
      (* Prevent bunched-up x ticks *)
      PlotRange -> {{1., 1.8}, {Automatic, Automatic}}
    ]
  ] // Ex[dest]
]


(* ::Section:: *)
(*Figure: known solution (known-solution)*)


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    b,
    xMin, xMax, yMax,
    dummyForTrailingCommas
  },
  (* Values of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = xStraight;
  yMax = 2;
  (* Make plot *)
  Show[
    Plot3D[
      tKnown[b][x, y]
      , {x, xMin, xMax}
      , {y, -yMax, yMax}
      , AxesEdge -> {{-1, -1}, {-1, -1}, {-1, +1}}
      , AxesLabel -> {xIt, yIt, Italicise["T"] // Margined @ {{0, 5}, {0, 0}}}
      , ClippingStyle -> BoundaryTracingStyle["Unphysical"]
      , LabelStyle -> Directive[Black, 18]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , PlotPoints -> 50
      , PlotRange -> {0, Full}
      , PlotStyle -> SlidesStyle["InteriorRegion"]
      , TicksStyle -> 14
      , ViewPoint -> {-2.6, -1.4, 1}
    ],
    Graphics3D @ {
      Dashed, SlidesStyle["Source"], Thickness[0.01],
      Line @ {{xStraight, -yMax, 1}, {xStraight, +yMax, 1}},
      {}
    },
    {}
    , ImageSize -> 360
  ]
] // Ex["cosine_simple-known-solution-slides.png"
  , Background -> None
];


(* ::Section:: *)
(*Figure: (Un)physical region (cosine-physical.pdf)*)


Module[
 {bStep, bValues, bMax,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinCont, xMaxCont, yMaxCont,
  yMaxContStraight,
  numTo1, numBeyond1,
  plotList,
  textStyle, arrowStyle,
  parameterArrow,
  xGraphicsBTick,
  legendLabelStyle,
  legendRegions, legendCurves,
  dummyForTrailingCommas
 },
  (* Values of B *)
  bStep = 3/10;
  bValues = {1 - bStep, 1, 1 + bStep};
  bMax = Max[bValues];
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = SeekRoot[tKnown[bMax][#, yMax] &, {0, xStraight}] + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + 2 eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 2;
  (* List of plots *)
  plotList =
    Table[
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax},
          Frame -> None,
          ImageSize -> Automatic,
          PlotRangePadding -> None
        ],
        (* Unphysical domain *)
        RegionPlot[
          tKnown[b][x, y] < 0,
          {x, xMin, xMax}, {y, -yMax, yMax},
          BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
          PlotPoints -> 6,
          PlotStyle -> BoundaryTracingStyle["Unphysical"]
        ],
        (* Known solution contours *)
        ContourPlot[
          tKnown[b][x, y],
          {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
          ContourLabels -> None,
          Contours -> numTo1 + numBeyond1,
          ContourShading -> None,
          ContourStyle -> BoundaryTracingStyle["ContourPlain"],
          PlotPoints -> If[b <= 1, 5, 4],
          PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
        ],
        (* Straight contour *)
        ParametricPlot[{xStraight, y},
          {y, -yMaxContStraight, yMaxContStraight},
          PlotPoints -> 2,
          PlotRange -> Full,
          PlotStyle -> BoundaryTracingStyle["ContourImportant"]
        ]
      ]
    , {b, bValues}];
  (* Parameter (B) increase indicator arrow *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* B-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["B"] // textStyle,
          {1, 0},
          {-2, 0}
        ]
      },
      (* B == 1 *)
      xGraphicsBTick = 0.48;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsBTick, 0},
          {xGraphicsBTick, -0.01}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsBTick, 0}
          , {0, 1.5}
        ]
      },
      PlotRange -> All
    ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendRegions =
    RegionLegend[
      {BoundaryTracingStyle["Unphysical"]},
      {"unphysical region"}
      , LabelStyle -> legendLabelStyle
    ];
  legendCurves =
    CurveLegend[
      {BoundaryTracingStyle["ContourPlain"]},
      {Row @ {Italicise["T"], "\[Hyphen]contour"}}
      , LabelStyle -> legendLabelStyle
    ];
  (* Final figure *)
  {
    GraphicsRow[plotList
      , ImageSize -> ImageSizeTextWidth
      , Spacings -> {0.2 ImageSizeTextWidth, 0}
    ] // Ex["cosine-physical.pdf"]
    ,
    Show[parameterArrow
      , ImageSize -> ImageSizeTextWidth
    ] // Ex["cosine-physical-arrow.pdf"]
    ,
    GraphicsGrid[List @ Join[legendRegions, legendCurves]
      , AspectRatio -> 0.1
      , ImageSize -> ImageSizeTextWidth
      , Spacings -> {{0, -0.3} ImageSizeTextWidth, 0}
    ] // Ex["cosine-physical-legend.pdf"]
    ,
    Nothing
  }
]


(* ::Section:: *)
(*Figure: Simple case (un)physical region and (non-)viable domain (cosine_simple-physical-viable.pdf)*)


Module[
 {b,
  aStep, aValues,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinCont, xMaxCont, yMaxCont,
  xMinViable, xMaxViable, yMaxViable,
  xMaxTerminal,
  yMaxContStraight,
  numTo1, numBeyond1,
  plotList,
  textStyle, arrowStyle,
  parameterArrow,
  xGraphicsATick,
  legendLabelStyle,
  legendRegions, legendCurves,
  dummyForTrailingCommas
 },
  (* Value of B *)
  b = 1;
  (* Values of A *)
  aStep = 5/10;
  aValues = {1 - aStep, 1, 1 + aStep};
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin;
  xMaxUnphys = xMax;
  yMaxUnphys = yMax;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 2;
  (* List of plots *)
  plotList =
    Table[
      (* Plot range for viable domain *)
      xMinViable = x0Simp[a] - eps;
      xMaxViable = xMax + eps;
      yMaxViable = yMax + eps;
      (* Plot range for terminal curve *)
      xMaxTerminal =
        SeekFirstRootBisection[
          vi[a, b][#, yMax] &,
          {xMinViable, xMax}
        ] + eps;
      (* Plot *)
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax},
          Frame -> None,
          ImageSize -> Automatic,
          PlotRangePadding -> None
        ],
        (* Unphysical domain *)
        RegionPlot[
          tKnown[b][x, y] < 0,
          {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
          BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
          PlotPoints -> 6,
          PlotStyle -> BoundaryTracingStyle["Unphysical"]
        ],
        (* Non-viable domain *)
        RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
          {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
          BoundaryStyle -> None,
          PlotPoints -> 2,
          PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        (* Known solution contours *)
        ContourPlot[
          tKnown[b][x, y],
          {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
          ContourLabels -> None,
          Contours -> numTo1 + numBeyond1,
          ContourShading -> None,
          ContourStyle -> BoundaryTracingStyle["ContourPlain"],
          PlotPoints -> 5,
          PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
        ],
        (* Straight contour *)
        ParametricPlot[{xStraight, y},
          {y, -yMaxContStraight, yMaxContStraight},
          PlotPoints -> 2,
          PlotRange -> Full,
          PlotStyle -> BoundaryTracingStyle["ContourImportant"]
        ],
        (* Terminal curve *)
        ContourPlot[vi[a, b][x, y] == 0,
          {x, xMinViable, xMaxTerminal}, {y, -yMaxViable, yMaxViable},
          ContourStyle -> BoundaryTracingStyle["Terminal"],
          PlotPoints -> 3
        ]
      ]
    , {a, aValues}];
  (* Parameter (A) increase indicator arrow *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* A-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["A"] // textStyle,
          {1, 0},
          {-2, 0}
        ]
      },
      (* A == 1 *)
      xGraphicsATick = 0.51;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsATick, 0},
          {xGraphicsATick, -0.01}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsATick, 0}
          , {0, 1.5}
        ]
      },
      PlotRange -> All
    ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendRegions =
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable", "Unphysical"},
      {"non\[Hyphen]viable domain", "unphysical region"}
      , LabelStyle -> legendLabelStyle
    ];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Terminal", "ContourPlain"},
      {"terminal curve", Row @ {Italicise["T"], "\[Hyphen]contour"}}
      , LabelStyle -> legendLabelStyle
    ];
  (* Final figure *)
  {
    GraphicsRow[plotList
      , ImageSize -> ImageSizeTextWidth
      , Spacings -> {0.2 ImageSizeTextWidth, 0}
    ] // Ex["cosine_simple-physical-viable.pdf"]
    ,
    Show[parameterArrow
      , ImageSize -> ImageSizeTextWidth
    ] // Ex["cosine_simple-physical-viable-arrow.pdf"]
    ,
    GraphicsGrid[{legendRegions, legendCurves}
      , Alignment -> Left
      , ImageSize -> ImageSizeTextWidth
      , ItemAspectRatio -> 0.1
      , Spacings -> {{0, -0.08} ImageSizeTextWidth, -0.02 ImageSizeTextWidth}
    ] // Ex["cosine_simple-physical-viable-legend.pdf"]
    ,
    Nothing
  }
]


(* ::Section:: *)
(*Figure: Simple case traced boundaries (cosine_simple-traced-boundaries.pdf)*)


Module[
 {a, b,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinViable, xMaxViable, yMaxViable,
  yMaxContStraight
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 0.42 ImageSizeTextWidth;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin;
  xMaxUnphys = xMax;
  yMaxUnphys = yMax;
  (* Plot range for viable domain *)
  xMinViable = x0Simp[a] - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -17}},
        Italicise["y"]
      },
      FrameTicksStyle -> LabelSize["Tick"],
      ImageSize -> imageSize,
      LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"],
      PlotRangePadding -> None
    ],
    (* Unphysical domain *)
    RegionPlot[
      tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
      PlotPoints -> 5,
      PlotStyle -> BoundaryTracingStyle["Unphysical"]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> None,
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {xy, xyTraSimpFigure}]
  ]
] // Ex["cosine_simple-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: Simple case traced boundaries, patched spiky (cosine_simple-traced-boundaries-patched-spiky.pdf)*)


Module[
 {a, b,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinViable, xMaxViable, yMaxViable,
  yMaxContStraight
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 0.36 ImageSizeTextWidth;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin;
  xMaxUnphys = xMax;
  yMaxUnphys = yMax;
  (* Plot range for viable domain *)
  xMinViable = x0Simp[a] - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      Frame -> None,
      ImageSize -> imageSize,
      PlotRangePadding -> None
    ],
    (* Traced boundaries (background) *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Background"]
      ]
    , {xy, Join[patchedBoundaryUpperList, patchedBoundaryLowerList]}],
    (* Unphysical domain *)
    RegionPlot[
      tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
      PlotPoints -> 6,
      PlotStyle -> BoundaryTracingStyle["Unphysical"]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> None,
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries (lower) *)
    Table[
      ParametricPlot[
        patchedBoundaryLowerList[[i]][s]
          // Through
          // Evaluate,
        {s, patchedIntersectionLowerList[[i]], 0}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {i, patchedCornerNum}],
    (* Traced boundaries (upper) *)
    Table[
      ParametricPlot[
        patchedBoundaryUpperList[[i]][s]
          // Through
          // Evaluate,
        {s, 0, patchedIntersectionUpperList[[i]]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {i, patchedCornerNum}],
    {}
  ]
] // Ex["cosine_simple-traced-boundaries-patched-spiky.pdf"]


(* ::Section:: *)
(*Figure: Simple case traced boundaries, patched smooth (cosine_simple-traced-boundaries-patched-smooth.pdf)*)


Module[
 {a, b,
  x0,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinViable, xMaxViable, yMaxViable,
  yMaxContStraight,
  textStyle, textStyleBracket
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Critical terminal x-coordinate x_0 *)
  x0 = x0Simp[a];
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 0.36 ImageSizeTextWidth;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin;
  xMaxUnphys = xMax;
  yMaxUnphys = yMax;
  (* Plot range for viable domain *)
  xMinViable = x0Simp[a] - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      Frame -> None,
      ImageSize -> imageSize,
      PlotRangePadding -> None
    ],
    (* Traced boundaries (background) *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, DomainStart[xy], 0}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Background"]
      ]
    , {xy, {patchedBoundarySmooth}}],
    (* Unphysical domain *)
    RegionPlot[
      tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
      PlotPoints -> 6,
      PlotStyle -> BoundaryTracingStyle["Unphysical"]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> None,
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries (smooth patching) *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, 0, DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {xy, {patchedBoundarySmooth}}],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {x0, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]0"],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {x0, 0},
        {1.5, -0.2}
      ]
        // textStyle
    },
    {}
  ]
] // Ex["cosine_simple-traced-boundaries-patched-smooth.pdf"]


(* ::Section:: *)
(*Figure: Simple case terminal points (cosine_simple-terminal-points.pdf)*)


Module[
 {a, b,
  x0,
  eps,
  xMin, xMax, yMax,
  xMinCont, xMaxCont, yMaxCont,
  xMinViable, xMaxViable, yMaxViable,
  xMaxTerminal,
  contNum, tContMin, tContStep, tContValues,
  tContOrd, xOrdGuess, yOrdGuess, xOrd, yOrd,
  textStyle, textStyleBracket, textVerticalShift,
  plot,
  legendLabelStyle,
  legendCurves, legendRegions,
  dummyForTrailingCommas
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Critical terminal x-coordinate x_0 *)
  x0 = x0Simp[a];
  (* Margin *)
  eps = 0.01;
  (* Plot range *)
  xMin = x0 - 0.5 (xStraight - x0);
  xMax = xStraight;
  yMax = 1;
  (* Plot range for contours *)
  xMinCont = x0 - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Contours *)
  contNum = 4;
  tContMin = tKnown[b][x0, 0];
  tContStep = 0.03;
  tContValues = tContMin + tContStep (Range[contNum] - 1);
  (* Plot range for viable domain *)
  xMinViable = x0 - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Plot range for terminal curve *)
  xMaxTerminal =
    SeekFirstRootBisection[
      vi[a, b][#, yMax] &,
      {xMinViable, xMax}
    ] + eps;
  (* Determine ordinary terminal point *)
  tContOrd = tContMin + tContStep;
  xOrdGuess = Way[x0, xMax];
  yOrdGuess = Way[0, yMax];
  {xOrd, yOrd} = FindRoot[
    Function[{x, y},
      {tKnown[b][x, y] - tContOrd, vi[a, b][x, y]}
    ],
    {{xOrdGuess}, {yOrdGuess}}
  ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textVerticalShift = -0.25;
  (* Plot *)
  plot = Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      AspectRatio -> 2,
      Frame -> None,
      PlotRangePadding -> None
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> None,
      PlotPoints -> 3,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Known solution contours *)
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
      ContourLabels -> None,
      Contours -> tContValues,
      ContourShading -> None,
      ContourStyle -> BoundaryTracingStyle["ContourPlain"],
      PlotPoints -> 6
    ],
    (* Terminal curve *)
    ContourPlot[vi[a, b][x, y] == 0,
      {x, xMinViable, xMaxTerminal}, {y, -yMaxViable, yMaxViable},
      ContourStyle -> BoundaryTracingStyle["Terminal"],
      PlotPoints -> 5
    ],
    (* Ordinary terminal point (x_ord, y_ord) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xOrd, yOrd}
    },
    Graphics @ {
      Text[
        "ordinary",
        {xOrd, yOrd},
        {-1.4, textVerticalShift}
      ] // textStyle,
      {}
    },
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {x0, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]0"],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {x0, 0},
        {1.5, -0.2}
      ] // textStyle,
      Text[
        "critical",
        {x0, 0},
        {-1.45, textVerticalShift}
      ] // textStyle,
      {}
    },
    {}
  ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"ContourPlain", "Terminal"},
      {Row @ {Italicise["T"], "\[Hyphen]contour"}, "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ];
  legendRegions =
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ];
  (* Combined *)
  GraphicsRow[
    {
      plot,
      Column[Join[legendCurves, legendRegions]
        , Spacings -> {0, {-1.5, -1.5, -1.3}}
      ]
    }
    , ItemAspectRatio -> 2 (* took a while to figure this out *)
    , ImageSize -> 0.58 ImageSizeTextWidth
    , Spacings -> {0.1 ImageSizeTextWidth, 0}
  ]
] // Ex["cosine_simple-terminal-points.pdf"]


(* ::Section:: *)
(*Figure: Simple case candidate domains (cosine_simple-candidate-domains.pdf)*)


Module[
  {
    b,
    numA,
    aStep, aValues,
    aMin, xRadAMin, yEndAMin,
    eps,
    yMaxFrame, xMinFrame, xMaxFrame,
    numCandidates, imageSize,
    plotList,
    xRad, yEnd,
    xRadEvenExtension,
    yInfl, xInfl,
    textStyle, arrowStyle,
    parameterArrow,
    xGraphicsAInfl, xGraphicsA1,
    legendLabelStyle,
    legendCurves, legendNodes,
    textStyleBracket,
    dummyForTrailingCommas
   },
  (* Value of B *)
  b = 1;
  (* Number of A from A_i to 1 inclusive *)
  numA = 4;
  (* Values of A *)
  aStep = (1 - aInflSimp) / (numA - 1);
  aValues = aInflSimp + aStep * Range[-2, numA - 1];
  (* Plot range *)
  aMin = Min[aValues];
  xRadAMin = xTraCandSimp[aMin, True];
  yEndAMin = DomainEnd[xRadAMin];
  eps = 0.02;
  yMaxFrame = yEndAMin;
  xMinFrame = (1 - eps) x0Simp[aMin];
  xMaxFrame = 1.1 xStraight;
  (* Image width arithmetic *)
  numCandidates = Length[aValues];
  imageSize["Main"] = 0.55 ImageSizeTextWidth;
  imageSize["Legend"] = 0.4 ImageSizeTextWidth;
  imageSize["Candidate"] = imageSize["Main"] / (numCandidates + 1);
  (* List of plots *)
  plotList =
    Table[
      (* Radiation boundary x == x(y) for convex domain *)
      xRad = xTraCandSimp[a, True];
      yEnd = DomainEnd[xRad];
      xRadEvenExtension = Function[{y}, xRad[Abs @ y] // Evaluate];
      (* Point of inflection *)
      Which[
        a < aInflSimp,
          yInfl = SeekRoot[
            curTra[a, b][xRad[#], #] &,
            {DomainStart[xRad], DomainEnd[xRad]}
          ],
        a == aInflSimp,
          yInfl = DomainEnd[xRad],
        True,
          yInfl = Indeterminate
      ];
      xInfl = xRad[yInfl];
      (* Plot *)
      Show[
        EmptyFrame[
          {xMinFrame, xMaxFrame}, {-yMaxFrame, yMaxFrame}
          , Frame -> None
          , ImageSize -> imageSize["Candidate"]
          , PlotRangePadding -> None
        ],
        (* Convex domain *)
        ParametricPlot[
          {
            {xRadEvenExtension[y], y},
            {xStraight, y}
          }
          , {y, -yEnd, yEnd}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle /@ {"Traced", "Contour"}
        ],
        (* Point of inflection *)
        If[NumericQ[yInfl],
          Graphics @ {GeneralStyle["Point"],
            Point @ {{xInfl, yInfl}, {xInfl, -yInfl}}
          },
          {}
        ],
        {}
      ]
      , {a, aValues}
    ];
  (* Parameter (A) increase indicator arrow *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* A-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["A"] // textStyle
          , {1, 0}
          , {-2, 0}
        ]
      },
      (* A == A_infl *)
      xGraphicsAInfl = 0.388;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsAInfl, 0},
          {xGraphicsAInfl, -0.02}
        }
      },
      Graphics @ {
        Text[
          Subscript[Italicise["A"], "i"] == SignificantFiguresForm[5][aInflSimp]
            // textStyle
          , {xGraphicsAInfl, 0}
          , {0, 1.2}
        ]
      },
      (* A == 1 *)
      xGraphicsA1 = 0.859;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsA1, 0},
          {xGraphicsA1, -0.02}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsA1, 0}
          , {0, 1.2}
        ]
      },
      (* Convexity-interval labels *)
      Graphics @ {
        Text[
          "non\[Hyphen]convex" // textStyle
          (* NOTE: use '\\[Hyphen]' to prevent parsing as minus sign *)
          , {Way[0, xGraphicsAInfl, 0.47], 0}
          , {0, -1.2}
        ],
        Text[
          "convex" // textStyle
          , {Way[xGraphicsAInfl, xGraphicsA1], 0}
          , {0, -1.2}
        ],
        {}
      },
      {}
      , ImageSize -> imageSize["Main"]
      , PlotRange -> All
      , PlotRangePadding -> None
    ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced", "Contour"},
      {"radiation", "constant temperature"}
      , LabelStyle -> legendLabelStyle
    ];
  legendNodes =
    NodeLegend[
      {Automatic},
      {
        textStyleBracket = Style[#, LabelSize["Legend"] + 2] &;
        Row @ {
          "inflection",
          " ",
          "(" // textStyleBracket,
            "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], "\[VeryThinSpace]i"],
          ",",
          "\[VeryThinSpace]",
          Subscript[Italicise["y"], "i"],
          ")" // textStyleBracket,
          Nothing
        }
      }
      , LabelStyle -> legendLabelStyle
      , LegendMarkerSize -> Automatic
    ];
  (* Export *)
  {
    (* Main *)
    Column[
      {
        Grid[
          List @ Append[
            plotList,
            Graphics[ImageSize -> imageSize["Candidate"]]
          ]
          , Spacings -> 0
        ],
        parameterArrow
      }
      , Spacings -> 0.1
    ] // Ex["cosine_simple-candidate-domains.pdf"]
    ,
    (* Legend *)
    GraphicsColumn[
      Join[legendCurves, legendNodes]
      , Alignment -> Left
      , ImageSize -> imageSize["Legend"]
      , ItemAspectRatio -> 0.1
      , Spacings -> 0
    ]
      // Ex["cosine_simple-candidate-domains-legend.pdf"]
  }
]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    b,
    numA,
    aStep, aValues,
    aMin, xRadAMin, yEndAMin,
    eps,
    yMaxFrame, xMinFrame, xMaxFrame,
    numCandidates, imageSize,
    plotList,
    xRad, yEnd,
    xRadEvenExtension,
    yInfl, xInfl,
    textStyle, arrowStyle,
    parameterArrow,
    xGraphicsAInfl, xGraphicsA1,
    legendLabelStyle,
    legendCurves, legendNodes,
    textStyleBracket,
    dummyForTrailingCommas
   },
  (* Value of B *)
  b = 1;
  (* Number of A from A_i to 1 inclusive *)
  numA = 4;
  (* Values of A *)
  aStep = (1 - aInflSimp) / (numA - 1);
  aValues = aInflSimp + aStep * Range[-2, numA - 1];
  (* Plot range *)
  aMin = Min[aValues];
  xRadAMin = xTraCandSimp[aMin, True];
  yEndAMin = DomainEnd[xRadAMin];
  eps = 0.02;
  yMaxFrame = yEndAMin;
  xMinFrame = (1 - eps) x0Simp[aMin];
  xMaxFrame = 1.1 xStraight;
  (* Image width arithmetic *)
  numCandidates = Length[aValues];
  imageSize["Main"] = 0.5 ImageSizeTextWidthBeamer;
  imageSize["Legend"] = 0.4 ImageSizeTextWidth;
  imageSize["Candidate"] = imageSize["Main"] / (numCandidates + 1);
  (* List of plots *)
  plotList =
    Table[
      (* Radiation boundary x == x(y) for convex domain *)
      xRad = xTraCandSimp[a, True];
      yEnd = DomainEnd[xRad];
      xRadEvenExtension = Function[{y}, xRad[Abs @ y] // Evaluate];
      (* Point of inflection *)
      Which[
        a < aInflSimp,
          yInfl = SeekRoot[
            curTra[a, b][xRad[#], #] &,
            {DomainStart[xRad], DomainEnd[xRad]}
          ],
        a == aInflSimp,
          yInfl = DomainEnd[xRad],
        True,
          yInfl = Indeterminate
      ];
      xInfl = xRad[yInfl];
      (* Plot *)
      Show[
        EmptyFrame[
          {xMinFrame, xMaxFrame}, {-yMaxFrame, yMaxFrame}
          , Frame -> None
          , ImageSize -> imageSize["Candidate"]
          , PlotRangePadding -> None
        ],
        (* Convex domain *)
        ParametricPlot[
          {
            {xStraight, y},
            {xRadEvenExtension[y], y}
          }
          , {y, -yEnd, yEnd}
          , PlotPoints -> 2
          , PlotStyle -> {
              Directive[BoundaryTracingStyle["Contour"], SlidesStyle["Source"]],
              Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
            }
        ],
        (* Point of inflection *)
        If[NumericQ[yInfl],
          Graphics @ {GeneralStyle["Point"],
            Point @ {{xInfl, yInfl}, {xInfl, -yInfl}}
          },
          {}
        ],
        {}
      ]
      , {a, aValues}
    ];
  (* Parameter (A) increase indicator arrow *)
  textStyle = Style[#, 10] &;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.07]];
  parameterArrow =
    Show[
      (* A-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["A"] // textStyle
          , {1, 0}
          , {-1.5, 0}
        ]
      },
      (* A == A_infl *)
      xGraphicsAInfl = 0.38;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsAInfl, 0},
          {xGraphicsAInfl, -0.02}
        }
      },
      Graphics @ {
        Text[
          Subscript[Italicise["A"], "i"] == SignificantFiguresForm[5][aInflSimp]
            // textStyle
          , {xGraphicsAInfl, 0}
          , {0, 1.5}
        ]
      },
      (* A == 1 *)
      xGraphicsA1 = 0.84;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsA1, 0},
          {xGraphicsA1, -0.02}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsA1, 0}
          , {0, 1.5}
        ]
      },
      (* Convexity-interval labels *)
      Graphics @ {
        Text[
          "non\[Hyphen]convex" // textStyle
          (* NOTE: use '\\[Hyphen]' to prevent parsing as minus sign *)
          , {Way[0, xGraphicsAInfl, 0.47], 0}
          , {0, -1.2}
        ],
        Text[
          "convex" // textStyle
          , {Way[xGraphicsAInfl, xGraphicsA1], 0}
          , {0, -1.2}
        ],
        {}
      },
      {}
      , ImageSize -> imageSize["Main"]
      , PlotRange -> All
      , PlotRangePadding -> None
    ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced", "Contour"},
      {"radiation", "constant temperature"}
      , LabelStyle -> legendLabelStyle
    ];
  legendNodes =
    NodeLegend[
      {Automatic},
      {
        textStyleBracket = Style[#, LabelSize["Legend"] + 2] &;
        Row @ {
          "inflection",
          " ",
          "(" // textStyleBracket,
            "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], "\[VeryThinSpace]i"],
          ",",
          "\[VeryThinSpace]",
          Subscript[Italicise["y"], "i"],
          ")" // textStyleBracket,
          Nothing
        }
      }
      , LabelStyle -> legendLabelStyle
      , LegendMarkerSize -> Automatic
    ];
  (* Export *)
  {
    (* Main *)
    Column[
      {
        Grid[
          List @ Append[
            plotList,
            Graphics[ImageSize -> imageSize["Candidate"]]
          ]
          , Spacings -> 0
        ],
        parameterArrow
      }
      , Spacings -> 0.1
    ] // Ex["cosine_simple-candidate-domains-slides.pdf"]
  }
];


(* ::Section:: *)
(*Figure: Simple case self-viewing bounds (cosine_simple-self-viewing-bounds-*.pdf)*)


(* Here we use a fake traced boundary (tanh) to exaggerate the features. *)
Module[
  {
    yStart, yEnd,
    xStart, xEnd,
    fun, xTraced,
    xInflection, yInflection,
    xView, yView,
    textStyle, textStyleBracket, textStyleLocalPosition,
    coordinatePairLabel,
      yEndLabel,
      xInflectionLabel, yInflectionLabel,
      xViewLabel, yViewLabel,
      y1Label, y2Label,
    tangentStyle, goodStyle, selfStyle, rayArrowStyle,
    pointStyle,
    yTickLength, yTick,
    markLength, localPositionMark, localPositionLabel,
    xa, ya, xb, yb, xc, yc,
    commonPlot,
    xLocal, yLocal, y1, y2,
    yStarValues,
    dummyForTrailingCommas
  },
  (* Radiation boundary range *)
  {yStart, yEnd} = {-1.2, 2};
  {xStart, xEnd} = {0.2, xStraight};
  (* Fake traced boundary for a radiation boundary *)
  fun[y_] := Tanh[2.5y];
  xTraced[y_] :=
    Way[
      xStart, xEnd,
      (fun[yStart] - fun[y]) / (fun[yStart] - fun[yEnd])
    ] // Evaluate;
  (* Point of inflection *)
  yInflection = 0; (* since xTraced is a rescaled version of tanh *)
  xInflection = xTraced[yInflection];
  (* Lowest visible point (self-viewing extremity) *)
  yView = SeekRoot[xTraced[#] + (yEnd - #) xTraced'[#] - Pi/2 &, {yStart, yInflection}];
  xView = xTraced[yView];
  (* Text functions *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textStyleLocalPosition = Style[#, LabelSize["Point"] + 1] & @* LaTeXStyle;
  coordinatePairLabel[x_, y_] :=
    Row @ {
      "(" // textStyleBracket,
      "",
      x,
      ",\[ThinSpace]",
      y,
      ")" // textStyleBracket
    } // textStyle;
  (* Subscripts with invisible spacing *)
  yEndLabel = Subscript[Italicise["y"], "\[VeryThinSpace]e"];
  xInflectionLabel = Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]i"];
  yInflectionLabel = Subscript[Italicise["y"], "\[VeryThinSpace]i"];
  xViewLabel = Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]v"];
  yViewLabel = Subscript[Italicise["y"], "\[VeryThinSpace]v"];
  y1Label = Subscript[Italicise["y"], "\[VeryThinSpace]1"];
  y2Label = Subscript[Italicise["y"], "\[VeryThinSpace]2"];
  (* Plot styles *)
  tangentStyle = Directive[Thickness[Small], GeneralStyle["Dashed"]];
  goodStyle = Directive[Black, GeneralStyle["Thick"], CapForm["Butt"]];
  selfStyle = Directive[GrayLevel[0.7], GeneralStyle["Thick"], CapForm["Butt"]];
  rayArrowStyle[p_] := Arrowheads @ {{0.07, p}};
  pointStyle = PointSize[0.05];
  yTickLength = (xEnd - xStart) / 25;
  yTick[y_, xStart_: xStraight] := {
    If[xStart < xStraight,
      {BoundaryTracingStyle["Contour"], Gray, Line @ {{xStart, y}, {xStraight, y}}},
      {}
    ],
    {Black, Line @ {{xStraight, y}, {xStraight + yTickLength, y}}},
    {}
  };
  markLength = (xEnd - xStart) / 18;
  localPositionMark[y_, ang_: 0] := (
    {Directive[AbsoluteThickness[2.5], Black],
      Line[markLength / Sqrt[2] * {{-1, -1}, {+1, +1}}],
      Line[markLength / Sqrt[2] * {{+1, -1}, {-1, +1}}],
      {}
    }
      // Rotate[#, ang, {0, 0}] &
      // Translate[#, {xTraced[y], y}] &
  );
  localPositionLabel =
    coordinatePairLabel[
      Italicise["x"],
      Row @ {Italicise["y"], "\[VeryThinSpace]"}
    ] // textStyleLocalPosition;
  (* Make common plot *)
  commonPlot[endStyle_: goodStyle, inflectionStyle_: goodStyle] :=
    Show[
      (* Constant-temperature boundary *)
      Graphics @ {
        BoundaryTracingStyle["Contour"],
        Line @ {{xStraight, yStart}, {xStraight, yEnd}},
        Text[
          xIt == SeparatedRow["VeryThin"]["\[Pi]", "/", 2] // textStyle
          , {xStraight, yStart}
          , {-1.3, -1.2}
          , {0, 1}
        ]
      },
      (* Uppermost tip *)
      Graphics @ {pointStyle,
        {endStyle, Point @ {xEnd, yEnd}},
        Text[
          coordinatePairLabel[
            SeparatedRow["VeryThin"]["\[Pi]", "/", 2],
            yEndLabel
          ]
          , {xEnd, yEnd}
          , {1.3, -0.2}
        ],
        {}
      },
      (* Point of inflection *)
      Graphics @ {pointStyle,
        {inflectionStyle, Point @ {xInflection, yInflection}},
        Text[
          coordinatePairLabel[xInflectionLabel, yInflectionLabel]
          , {xInflection, yInflection}
          , {-1.25, 0.5}
        ],
        {}
      },
      (* Lowest visible point (self-viewing extremity) *)
      Graphics @ {pointStyle,
        Point @ {xView, yView},
        Text[
          coordinatePairLabel[xViewLabel, yViewLabel]
          , {xView, yView}
          , {-1.25, 0.25}
        ],
        {}
      },
      {}
      , ImageSize -> 0.32 ImageSizeTextWidth
      , PlotRange -> {{xStart, xEnd}, {yStart, yEnd}}
      , PlotRangePadding -> {
          {Scaled[0.05], Scaled[0.17]},
          {Scaled[0.01], Scaled[0.03]}
        }
    ];
  (* Local position for each case *)
  (*
    We use the intersection of the tangent line in case (b)
    as the local position for case (c), because reciprocity.
  *)
  ya = Way[yStart, yView, 0.25];
  xa = xTraced[ya];
  yb = Way[yView, yInflection, 0.45];
  xb = xTraced[yb];
  yc = SeekRootBisection[xTraced[#] + (yb + #) xTraced'[#] - Pi/2 &, {yInflection, yEnd}];
  xc = xTraced[yc];
  {
    (*
      --------------------------------
      (a) y < y_v
      --------------------------------
    *)
    {xLocal, yLocal} = {xa, ya};
    Show[commonPlot[],
      (* Radiation boundary *)
      ParametricPlot[
        {xTraced[y], y}
        , {y, yStart, yEnd}
        , PlotPoints -> 2
        , PlotStyle -> goodStyle
      ],
      (* Tangent line through (x_v, y_v) *)
      Graphics @ {tangentStyle,
        HalfLine @ {{xEnd, yEnd}, {xView, yView}}
      },
      (* Local position (x, y) *)
      Graphics @ {
        localPositionMark[yLocal],
        Text[
          localPositionLabel
          , {xLocal, yLocal}
          , {-1.65, -0.15}
        ],
        {}
      },
      (* Inflection tick *)
      Graphics @ {
        yTick[yInflection],
        Text[
          yInflectionLabel // textStyle
          , {xStraight + yTickLength, yInflection}
          , {-1.8, -0.1}
        ],
        {}
      },
      (* Convexity labels *)
      Graphics @ {
        Text[
          "convex" // textStyle
          , {xStraight, Way[yStart, yInflection]}
          , {0, 1}
          , {0, 1}
        ],
        Text[
          "concave" // textStyle
          , {xStraight, Way[yInflection, yEnd]}
          , {0, 1}
          , {0, 1}
        ],
        {}
      },
      {}
    ]
      // Ex["cosine_simple-self-viewing-bounds-none.pdf"]
    ,
    (*
      --------------------------------
      (b) y_v < y < y_i
      --------------------------------
    *)
    {y1, y2} = {yc, yEnd};
    {xLocal, yLocal} = {xb, yb};
    Show[commonPlot[selfStyle],
      (* Radiation boundary *)
      ParametricPlot[
        {xTraced[y], y}
        , {y, yStart, y1}
        , PlotPoints -> 2
        , PlotStyle -> goodStyle
      ],
      ParametricPlot[
        {xTraced[y], y}
        , {y, y1, y2}
        , PlotPoints -> 2
        , PlotStyle -> selfStyle
      ],
      (* Tangent line through (x_b, y_y) *)
      Graphics @ {tangentStyle,
        HalfLine[{xb, yb}, {xb, yb} - {xc, yc}]
      },
      (* Local position (x, y) *)
      Graphics @ {
        localPositionMark[yLocal, -20 Degree],
        Text[
          localPositionLabel
          , {xLocal, yLocal}
          , {-1.7, 0.2}
        ],
        {}
      },
      (* Bound ticks *)
      Graphics @ {
        (* y_2 *)
        Text[
          y2Label // textStyle
          , {xStraight + yTickLength, y2}
          , {-1.4, -0.1}
        ],
        (* y_1 *)
        yTick[y1, xTraced[y1]],
        Text[
          y1Label // textStyle
          , {xStraight + yTickLength, y1}
          , {-1.7, -0.1}
        ],
        {}
      },
      (* Self-viewing rays *)
      yStarValues = Subdivide[y1, y2, 3];
      Graphics @ {
        rayArrowStyle[1/2],
        Table[
          {
            If[yStar > y1, Nothing, tangentStyle],
            Arrow @ {
              {xTraced[yStar], yStar},
              {xLocal, yLocal}
            }
          }
          , {yStar, yStarValues}
        ]
      },
      {}
    ]
      // Ex["cosine_simple-self-viewing-bounds-concave.pdf"]
    ,
    (*
      --------------------------------
      (c) y_i < y < y_e
      --------------------------------
    *)
    {y1, y2} = {yb, yEnd};
    {xLocal, yLocal} = {xc, yc};
    Show[commonPlot[selfStyle, selfStyle],
      (* Radiation boundary *)
      ParametricPlot[
        {xTraced[y], y}
        , {y, yStart, y1}
        , PlotPoints -> 2
        , PlotStyle -> goodStyle
      ],
      ParametricPlot[
        {xTraced[y], y}
        , {y, y1, y2}
        , PlotPoints -> 2
        , PlotStyle -> selfStyle
      ],
      (* Tangent line through (x_b, y_y) *)
      Graphics @ {tangentStyle,
        HalfLine[{xb, yb}, {xb, yb} - {xc, yc}]
      },
      (* Local position (x, y) *)
      Graphics @ {
        localPositionMark[yLocal, -10 Degree],
        Text[
          localPositionLabel
          , {xLocal, yLocal}
          , {1.65, -0.5}
        ],
        {}
      },
      (* Bound ticks *)
      Graphics @ {
        (* y_2 *)
        Text[
          y2Label // textStyle
          , {xStraight + yTickLength, y2}
          , {-1.4, -0.1}
        ],
        (* y_1 *)
        yTick[y1, xTraced[y1]],
        Text[
          y1Label // textStyle
          , {xStraight + yTickLength, y1}
          , {-1.7, -0.1}
        ],
        {}
      },
      (* Self-viewing rays *)
      yStarValues = Join[
        Subdivide[y1, yLocal, 2] // Most,
        {Way[yLocal, yEnd, 0.25], yEnd},
        {}
      ];
      Graphics @ {
        Table[
          {
            Which[
              yStar == y1, rayArrowStyle[0.65],
              yStar < yLocal, rayArrowStyle[0.5],
              yStar < Way[yLocal, yEnd], rayArrowStyle[0.5],
              True, rayArrowStyle[0.72]
            ],
            If[yStar > y1, Nothing, tangentStyle],
            Arrow @ {
              {xTraced[yStar], yStar},
              {xLocal, yLocal}
            }
          }
          , {yStar, yStarValues}
        ]
      },
      {}
    ]
      // Ex["cosine_simple-self-viewing-bounds-both.pdf"]
  }
]


(* ::Section:: *)
(*Figure: Simple case self-viewing ratio upper bound (cosine_simple-self-viewing-ratio-bound.pdf)*)


(*
  Copy, paste, edit from "A vs R_crude",
  in "Non-convex lens-like candidate self-radiation".
*)
Module[
  {
    aMin, aMax,
    rMin, rMax,
    textStyle,
    dummyForTrailingCommas
  },
  aMin = 0;
  aMax = aInflSimp;
  rMin = 10^-6;
  rMax = 10^2.7;
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  Show[
    ListLogPlot[
      candidateBoundaryRatioTable
      , AxesLabel -> {
          Italicise["A"] // Margined @ {{2, 0}, {5, 0}},
          Subscript[Italicise["R"], Style["+", Magnification -> 1.35]]
            // Margined @ {{0, 0}, {-5, 0}}
        }
      , Epilog -> {
          BoundaryTracingStyle["Contour"],
          Line @ {{aMin, #}, {aMax, #}} & [rPracCandSimp // Log],
          Line @ {{aPracCandSimp, rMin // Log}, {aPracCandSimp, rPracCandSimp // Log}},
          PointSize[Medium], Point @ {aPracCandSimp, rMin // Log},
          Text[
            Subscript[Italicise["A"], "p"] == SignificantFiguresForm[5][aPracCandSimp]
              // textStyle
            , {aPracCandSimp, rMin // Log}
            , {-1.05, -1.08}
          ],
          {}
        }
      , Joined -> True
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , PlotRange -> {rMin, rMax}
      , PlotRangeClipping -> False
      , PlotStyle -> Black
      , TicksStyle -> LabelSize["Tick"]
    ],
    {}
    , ImageSize -> 0.52 ImageSizeTextWidth
  ]
] // Ex["cosine_simple-self-viewing-ratio-bound.pdf"]


(* ::Section:: *)
(*Figure: General case (un)physical region and (non-)viable domain (cosine_general-physical-viable.pdf)*)


Module[
  {
    a, bN, bValues,
    xMin, xMax, yMax,
    eps,
    textStylePoint,
    plotList,
    xMinUnphys, xMaxUnphys, yMaxUnphys,
    plotPointsUnphys,
    xMinViable, xMaxViable, yMaxViable,
    plotPointsViable,
    textStyle, arrowStyle,
    parameterArrow,
    xGraphicsBNat, xGraphicsB1,
    legendLabelStyle,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = 1;
  (* Value of B at gentle-to-fair transition *)
  bN = bNat[a];
  (* Values of B *)
  bValues = {
    (* Case 1 *) 0.7 bN,
    (* Case 2 *) bN,
    (* Case 3 *) Way[bN, 1, 0.03],
    (* Case 4 *) 1,
    (* Case 5 *) 2,
    Nothing
  };
  (* Plot range *)
  xMin = 0;
  xMax = 2;
  yMax = 3;
  (* Margin *)
  eps = 0.1;
  (* List of plots *)
  textStylePoint = Style[#, LabelSize["Label"] + 1] & @* LaTeXStyle;
  plotList =
    Table[
      (* Plot range for unphysical domain *)
      xMinUnphys = xMin;
      xMaxUnphys = xMax;
      yMaxUnphys = yMax;
      plotPointsUnphys = If[b == 1, 5, 7];
      (* Plot range for viable domain *)
      xMinViable = If[b < 1, xMin - eps, xSharp[a, b] - eps];
      xMaxViable = xMax + eps;
      yMaxViable = yMax + eps;
      plotPointsViable = If[b < 1, 7, 4];
      (* Plot *)
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax}
          , Frame -> None
          , ImageSize -> Automatic
          , PlotRangePadding -> None
        ],
        (* Unphysical domain *)
        RegionPlot[
          tKnown[b][x, y] < 0
          , {x, xMinUnphys, xMaxUnphys}
          , {y, -yMaxUnphys, yMaxUnphys}
          , BoundaryStyle -> BoundaryTracingStyle["Unphysical"]
          , PlotPoints -> plotPointsUnphys
          , PlotStyle -> BoundaryTracingStyle["Unphysical"]
        ],
        (* Non-viable domain *)
        RegionPlot[
          vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0
          , {x, xMinViable, xMaxViable}
          , {y, -yMaxViable, yMaxViable}
          , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
          , PlotPoints -> plotPointsViable
          , PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        (* Critical terminal points *)
        Which[
          b < bN,
          {
            Graphics @ {GeneralStyle["Point"],
              Point @ {{0, N @ y0[a, b]}, {0, N @ -y0[a, b]}}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["y"], 0] // textStylePoint
                , {0, N @ y0[a, b]}
                , {-1.2, -1.2}
             ]
            },
            {}
          },
          b == bN,
          {
            Graphics @ {GeneralStyle["Point"],
              Point @ {{0, N @ y0[a, b]}, {0, N @ -y0[a, b]}}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["y"], 0] // textStylePoint
                , {0, N @ y0[a, b]}
                , {-1, -1.3}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xNat[a], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Natural]"] // textStylePoint
                , {N @ xNat[a], 0}
                , {0.3, -1.35}
             ]
            },
            {}
          },
          bN < b < 1,
          {
            Graphics @ {GeneralStyle["Point"],
              Point @ {{0, N @ y0[a, b]}, {0, N @ -y0[a, b]}}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["y"], 0] // textStylePoint
                , {0, N @ y0[a, b]}
                , {-0.9, -1.3}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xFlat[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Flat]"] // textStylePoint
                , {N @ xFlat[a, b], 0}
                , {-0.15, -1.25}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xSharp[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Sharp]"] // textStylePoint
                , {N @ xSharp[a, b], 0}
                , {0.9, -1.25}
             ]
            },
            {}
          },
          b == 1,
          {
            Graphics @ {Lighter[Gray, 0.75], GeneralStyle["Translucent"], GeneralStyle["Point"],
              Point @ {N @ xFlat[a, b], 0}
            },
            Graphics @ {White,
              Text[
                Subscript[Italicise["y"], 0] // textStylePoint
                , {N @ xFlat[a, b], 0}
                , {-1, -1.9}
             ]
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Flat]"] // textStylePoint
                , {N @ xFlat[a, b], 0}
                , {-2.1, -0.2}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xSharp[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Sharp]"] // textStylePoint
                , {N @ xSharp[a, b], 0}
                , {1.9, -0.2}
             ]
            },
            {}
          },
          b > 1,
          {
            Graphics @ {White, GeneralStyle["Point"],
              Point @ {N @ xFlat[a, b], 0}
            },
            Graphics @ {White,
              Text[
                Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Flat]"] // textStylePoint
                , {N @ xFlat[a, b], 0}
                , {-2, -0.2}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xSharp[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Sharp]"] // textStylePoint
                , {N @ xSharp[a, b], 0}
                , {1.9, -0.2}
             ]
            },
            {}
          },
          True, {}
        ],
        {}
      ]
      , {b, bValues}
    ];
  (* Parameter (B) increase indicator arrow *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* B-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["B"] // textStyle,
          {1, 0},
          {-2, 0}
        ]
      },
      (* B == B_nat *)
      xGraphicsBNat = 0.288;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsBNat, 0},
          {xGraphicsBNat, -0.01}
        }
      },
      Graphics @ {
        Text[
          Subscript[Italicise["B"], "\[Natural]"]
            // textStyle
          , {xGraphicsBNat, 0}
          , {0, 1.4}
        ]
      },
      (* B == 1 *)
      xGraphicsB1 = 0.725;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsB1, 0},
          {xGraphicsB1, -0.01}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsB1, 0}
          , {0, 1.4}
        ]
      },
      (* Regime labels *)
      Graphics @ {
        Text[
          "gentle" // textStyle
          , {Way[0, xGraphicsBNat, 0.25], 0}
          , {0, -1.3}
        ],
        Text[
          "fair" // textStyle
          , {Way[xGraphicsBNat, xGraphicsB1, 0.48], 0}
          , {0, -1.3}
        ],
        Text[
          "steep" // textStyle
          , {Way[xGraphicsB1, 1, 0.69], 0}
          , {0, -1.3}
        ],
        {}
      },
      {}
      , PlotRange -> All
    ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  (* Final figure *)
  {
    GraphicsRow[plotList
      , ImageSize -> ImageSizeTextWidth
      , Spacings -> {0.3 ImageSizeTextWidth, 0}
    ] // Ex["cosine_general-physical-viable.pdf"]
    ,
    Show[parameterArrow
      , ImageSize -> ImageSizeTextWidth
    ] // Ex["cosine_general-physical-viable-arrow.pdf"]
    ,
    GraphicsGrid[
      Transpose @ {
        RegionLegend[
          {BoundaryTracingStyle["NonViable"]},
          {"non\[Hyphen]viable domain"}
          , LabelStyle -> legendLabelStyle
        ],
        CurveLegend[
          {BoundaryTracingStyle["Terminal"]},
          {"terminal curve"}
          , LabelStyle -> legendLabelStyle
        ],
        RegionLegend[
          {BoundaryTracingStyle["Unphysical"]},
          {"unphysical region"}
          , LabelStyle -> legendLabelStyle
        ],
        Nothing
      }
      , Alignment -> Left
      , ImageSize -> ImageSizeTextWidth
      , ItemAspectRatio -> 0.11
    ] // Ex["cosine_general-physical-viable-legend.pdf"]
  }
]


(* ::Section:: *)
(*Figure: General case critical terminal points (cosine_general-critical.pdf)*)


Module[
  {
    a,
    bN, xN,
    bMin, bMax,
    xMin, xMax, xMaxFilling,
    xFillingFunction,
    visiblePlotStyles, plot,
    bMargined, xMargined,
    legendLabelStyleCurves, legendLabelStyleRegions,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = 1;
  (* Critical values B_\[Natural] and x_\[Natural] *)
  (* (N needed otherwise Epilog doesn't work) *)
  bN = bNat[a] // N;
  xN = xNat[a] // N;
  (* Plot range *)
  bMin = 0;
  bMax = 1.7;
  xMin = 0;
  xMax = 2;
  xMaxFilling = Way[Pi/2, xMax, 2/3];
  (* Filling function for B < B_\[Natural] *)
  xFillingFunction[b_] :=
    Piecewise @ {
      {xMin, b < bN},
      {Indeterminate, True}
    };
  (* Plot *)
  visiblePlotStyles = {
    Directive[Black, Thick],
    Directive[Black, Thick, GeneralStyle["Dashed"]]
  };
  plot = Plot[
    {xSharp[a, b], xFlat[a, b], xFillingFunction[b]}
    , {b, bMin, bMax}
    , AspectRatio -> 1 / 1.5
    , AxesLabel -> {
        bIt,
        xIt // Margined @ {{0, 0}, {-5, 0}}
      }
    , Epilog -> {
        (* Guiding lines for critical values *)
        BoundaryTracingStyle["Contour"],
        Line @ {{0, xN}, {bN, xN}, {bN, 0}},
        {}
      }
    , Filling -> {1 -> xMaxFilling, 2 -> xMin, 3 -> xMaxFilling}
    , FillingStyle -> BoundaryTracingStyle["NonViable"]
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Label"]
    , PlotPoints -> 5
    , PlotRange -> {xMin, xMax}
    , PlotRangeClipping -> False
    , PlotStyle -> Append[visiblePlotStyles, None]
    , Ticks -> {
        {
          {0, 0 // Margined @ {{0, 0}, {0, -13}}},
          {bN, Subscript[bIt, "\[Natural]"] // Margined @ {{0, 0}, {0, -12}}},
          {1, 1 // Margined @ {{0, 0}, {0, -11}}}
        },
        {0, {xN, Subscript[xIt, "\[Natural]"] // Margined @ {{0, 2}, {4, 0}}}}
      }
    , TicksStyle -> LabelSize["Label"]
  ];
  (* Legend *)
  legendLabelStyleCurves = LatinModernLabelStyle @ LabelSize["Label"];
  legendLabelStyleRegions = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      visiblePlotStyles,
      xIt == Subscript[xIt, #] &
        /@ {"\[Sharp]", "\[Flat]"}
      , LabelStyle -> legendLabelStyleCurves
    ];
  legendRegions =
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable"}
      , LabelStyle -> legendLabelStyleRegions
    ];
  (* Export *)
  {
    Show[
      plot
      , ImageSize -> 0.45 ImageSizeTextWidth
    ]
      // Ex["cosine_general-critical.pdf"]
    ,
    GraphicsColumn[
      Join[legendCurves, legendRegions]
      , Alignment -> Left
      , ImageSize -> 0.2 ImageSizeTextWidth
      , ItemAspectRatio -> 0.17
    ]
      // Ex["cosine_general-critical-legend.pdf"]
  }
]


(* ::Section:: *)
(*Figure: General case asymmetric convex domain construction (cosine_general-asymmetric-construction)*)


(* ::Subsection:: *)
(*General case traced boundary convex portions (cosine_general-traced-boundaries-convex.pdf)*)


Module[
  {
    a, b,
    yReflection, includeYReflection,
    textStyle, textStylePoint,
    xMin, xMax, yMax,
    margin,
    xMinNon, xMaxNon, yMaxNon,
    yMaxCont,
    xyInflectionOuter,
    sInflectionOuterLabel, xInflectionOuterLabel, yInflectionOuterLabel,
    xyInflectionInner,
    sInflectionInnerLabel, xInflectionInnerLabel, yInflectionInnerLabel,
    xyTracedLeftmost,
    sTracedLeftmostLabel, xTracedLeftmostLabel, yTracedLeftmostLabel,
    dummyForTrailingCommas
  },
  (* Values of A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Reflection in y (across x-axis) *)
  yReflection = # * {1, -1} &;
  includeYReflection = {#, yReflection[#]} &;
  (* Text style *)
  textStyle = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  textStylePoint = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  (* Plot range *)
  xMin = Floor[0.9 xFlat[a, b], 0.2];
  xMax = Ceiling[1.1 xSharp[a, b], 0.2];
  yMax = Ceiling[0.9 (xMax - xMin), 0.2];
  (* Absolute plot range margin *)
  margin = 0.1;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - margin;
  xMaxNon = xMax + margin;
  yMaxNon = (
    Table[
      SeekRoot[vi[a, b][x, #] &, {0, yMax}]
      , {x, {xMin, xMax}}
    ]
      // Max[#] + margin &
  );
  (* Plot range for straight contour *)
  yMaxCont = yMax + margin;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , ImageSize -> 0.45 ImageSizeTextWidth
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Inflection frontiers *)
    Table[
      ParametricPlot[
        xy[s] (* lower branch *)
          // Through
          // includeYReflection (* upper branch *)
          // Evaluate
          ,
        {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> BoundaryTracingStyle["Background"]
        , PlotPoints -> 2
      ]
      , {xy, asymmInflectionFrontierList}
    ],
    (* Inflection frontier labels (outer) *)
    xyInflectionOuter = asymmInflectionFrontierList // Last;
    xInflectionOuterLabel = Way[xMin, xStraight, 0.5];
    sInflectionOuterLabel = SeekRoot[
      xyInflectionOuter[[1]][#] - xInflectionOuterLabel &,
      {DomainStart[xyInflectionOuter], DomainEnd[xyInflectionOuter]}
    ];
    yInflectionOuterLabel = xyInflectionOuter[[2]] @ sInflectionOuterLabel;
    Graphics @ {Gray,
      Text[
        "lower" // textStyle
        , {xInflectionOuterLabel, yInflectionOuterLabel}
        , {0, -1.2}
      ],
      Text[
        "upper" // textStyle
        , {xInflectionOuterLabel, yInflectionOuterLabel} // yReflection
        , {0, 0.8}
      ],
      {}
    },
    (* Inflection frontier labels (inner) *)
    xyInflectionInner = asymmInflectionFrontierList // First;
    xInflectionInnerLabel = Way[xFlatAsymm, xInflAxis, 0.28];
    sInflectionInnerLabel = SeekRoot[
      xyInflectionInner[[2]]'[#] &,
      {DomainStart[xyInflectionInner], DomainEnd[xyInflectionInner]}
    ];
    {xInflectionInnerLabel, yInflectionInnerLabel} =
      xyInflectionInner[sInflectionInnerLabel] // Through;
    Graphics @ {Gray,
      Text[
        "lower" // textStyle
        , {xInflectionInnerLabel, yInflectionInnerLabel}
        , {0, 0.8}
      ],
      Text[
        "upper" // textStyle
        , {xInflectionInnerLabel, yInflectionInnerLabel} // yReflection
        , {0, -1.2}
      ],
      {}
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[a, b][x, y] < 0
      , {x, xMinNon, xMaxNon}
      , {y, -yMaxNon, yMaxNon}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 4
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Straight contour *)
    ParametricPlot[
      {xStraight, y},
      {y, -yMaxCont, yMaxCont}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    (* Convex portions of traced boundaries *)
    Table[
      ParametricPlot[
        xy[s] (* upper branch *)
          // Through
          // includeYReflection (* lower branch *)
          // Evaluate
          ,
        {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
        , PlotPoints -> 2
      ]
      , {xy, asymmConvexPortionsList}
    ],
    (* Traced boundary labels *)
    xyTracedLeftmost = asymmConvexPortionsList // First;
    sTracedLeftmostLabel =
      Way[DomainStart[xyTracedLeftmost], DomainEnd[xyTracedLeftmost], 0.63];
    {xTracedLeftmostLabel, yTracedLeftmostLabel} =
      xyTracedLeftmost[sTracedLeftmostLabel] // Through;
    Graphics @ {
      Text[
        "upper" // textStyle
        , {xTracedLeftmostLabel, yTracedLeftmostLabel}
        , {0, 0.7}
        , ParametricTangent[xyTracedLeftmost, sTracedLeftmostLabel]
      ],
      Text[
        "lower" // textStyle
        , {xTracedLeftmostLabel, yTracedLeftmostLabel} // yReflection
        , {0, -1.2}
        , ParametricTangent[xyTracedLeftmost, sTracedLeftmostLabel] // yReflection
      ],
      {}
    },
    (* Straight contour and constant temperature labels *)
    Graphics @ {
      Text[
        xIt == SeparatedRow["VeryThin"]["\[Pi]", "/", 2] // textStyle
        , {xStraight, Way[yInflectionOuterLabel, yMax, 2/3]}
        , {0, -1.2}
        , {0, 1}
      ]
    },
    Graphics @ {
      Text[
        Italicise["T"] == 1 // textStyle
        , {xStraight, -Way[yInflectionOuterLabel, yMax, 2/3]}
        , {0, -1.2}
        , {0, 1}
      ]
    },
    (* Critical terminal point (x_\[Flat], 0) *)
    Graphics @ {GeneralStyle["Point"],
      Point @ {N @ xFlatAsymm, 0}
    },
    Graphics @ {
      Text[
        Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Flat]"] // textStylePoint
        , {N @ xFlatAsymm, 0}
        , {-1.85, -0.15}
      ]
    },
    (* Critical terminal point (x_\[Sharp], 0) *)
    Graphics @ {GeneralStyle["Point"],
      Point @ {N @ xSharpAsymm, 0}
    },
    Graphics @ {
      Text[
        Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]\[Sharp]"] // textStylePoint
        , {N @ xSharpAsymm, 0}
        , {1.7, -0.15}
      ]
    },
    {}
  ]
] // Ex["cosine_general-traced-boundaries-convex.pdf"]


(* ::Subsection:: *)
(*General case asymmetric convex domain (cosine_general-asymmetric_domain.pdf)*)


Module[
  {
    a, b,
    yReflection, includeYReflection,
    textStyle,
    xMin, xMax, yMax,
    margin,
    yTop, yBottom,
    dummyForTrailingCommas
  },
  (* Values of A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Reflection in y (across x-axis) *)
  yReflection = # * {1, -1} &;
  includeYReflection = {#, yReflection[#]} &;
  (* Text style *)
  textStyle = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  (* Plot range *)
  xMin = Floor[0.9 xFlat[a, b], 0.2];
  xMax = Ceiling[1.1 xSharp[a, b], 0.2];
  yMax = Ceiling[0.9 (xMax - xMin), 0.2];
  (* Absolute plot range margin *)
  margin = 0.1;
  (* Endpoints for constant temperature boundary *)
  yTop = xyTraAsymm["lower"][[2]] @ DomainStart @ xyTraAsymm["lower"];
  yBottom = xyTraAsymm["upper"][[2]] @ DomainEnd @ xyTraAsymm["upper"];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , ImageSize -> 0.45 ImageSizeTextWidth
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Inflection frontiers *)
    Table[
      ParametricPlot[
        xy[s] (* lower branch *)
          // Through
          // includeYReflection (* upper branch *)
          // Evaluate
          ,
        {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> BoundaryTracingStyle["Background"]
        , PlotPoints -> 2
      ]
      , {xy, asymmInflectionFrontierList}
    ],
    (* Straight contour *)
    ParametricPlot[
      {xStraight, y},
      {y, yTop, yBottom}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    (* Domain radiation boundaries *)
    Table[
      ParametricPlot[
        xyTraAsymm[id][s] // Through,
        {s, DomainStart @ xyTraAsymm[id], DomainEnd @ xyTraAsymm[id]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {id, {"upper", "lower"}}
    ],
    (* Straight contour and constant temperature labels *)
    Graphics @ {
      Text[
        xIt == SeparatedRow["VeryThin"]["\[Pi]", "/", 2] // textStyle
        , {xStraight, Way[yBottom, yTop, 7/10]}
        , {0, 1}
        , {0, 1}
      ]
    },
    Graphics @ {
      Text[
        Italicise["T"] == 1 // textStyle
        , {xStraight, Way[yBottom, yTop, 3/10]}
        , {0, 1}
        , {0, 1}
      ]
    },
    {}
  ]
] // Ex["cosine_general-asymmetric_domain.pdf"]


(* ::Subsubsection:: *)
(*Version for slides*)


Module[
  {
    a, b,
    yReflection, includeYReflection,
    textStyle,
    xMin, xMax, yMax,
    margin,
    yTop, yBottom,
    dummyForTrailingCommas
  },
  (* Values of A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Reflection in y (across x-axis) *)
  yReflection = # * {1, -1} &;
  includeYReflection = {#, yReflection[#]} &;
  (* Text style *)
  textStyle = Style[#, 9] &;
  (* Plot range *)
  xMin = Floor[0.9 xFlat[a, b], 0.2];
  xMax = Ceiling[1.1 xSharp[a, b], 0.2];
  yMax = Ceiling[0.7 (xMax - xMin), 0.2];
  (* Absolute plot range margin *)
  margin = 0.1;
  (* Endpoints for constant temperature boundary *)
  yTop = xyTraAsymm["lower"][[2]] @ DomainStart @ xyTraAsymm["lower"];
  yBottom = xyTraAsymm["upper"][[2]] @ DomainEnd @ xyTraAsymm["upper"];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> 8
      , LabelStyle -> 11
    ],
    (* Line of symmetry *)
    Graphics @ {
      Gray, Dashed,
      Line @ {{xMin, 0}, {xMax, 0}}
    },
    (* Straight contour *)
    ParametricPlot[
      {xStraight, y},
      {y, yTop, yBottom}
      , PlotPoints -> 2
      , PlotStyle -> Directive[BoundaryTracingStyle["Contour"], SlidesStyle["Source"]]
    ],
    (* Domain radiation boundaries *)
    Table[
      ParametricPlot[
        xyTraAsymm[id][s] // Through,
        {s, DomainStart @ xyTraAsymm[id], DomainEnd @ xyTraAsymm[id]}
        , PlotPoints -> 2
        , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
      ]
      , {id, {"upper", "lower"}}
    ],
    (* Straight contour and constant temperature labels *)
    Graphics @ {
      SlidesStyle["Source"],
      Text[
        xIt == SeparatedRow["VeryThin"]["\[Pi]", "/", 2] // textStyle
        , {xStraight, Way[yBottom, yTop, 7/10]}
        , {0, 1.3}
        , {0, 1}
      ]
    },
    Graphics @ {
      SlidesStyle["Source"],
      Text[
        Italicise["T"] == 1 // textStyle
        , {xStraight, Way[yBottom, yTop, 3/10]}
        , {0, 1.5}
        , {0, 1}
      ]
    },
    {}
    , ImageSize -> 0.5 * 0.8 ImageSizeTextWidthBeamer
  ]
] // Ex["cosine_general-asymmetric_domain-slides.pdf"];


(* ::Subsubsection:: *)
(*Version for slides, with mesh*)


Module[
  {
    a, b,
    mesh,
    yReflection, includeYReflection,
    textStyle,
    xMin, xMax, yMax,
    margin,
    yTop, yBottom,
    dummyForTrailingCommas
  },
  (* Values of A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Finite element mesh *)
  mesh = Import["cosine_general-verification-mesh-asymmetric.txt"] // Uncompress // First;
  (* Reflection in y (across x-axis) *)
  yReflection = # * {1, -1} &;
  includeYReflection = {#, yReflection[#]} &;
  (* Text style *)
  textStyle = Style[#, 9] &;
  (* Plot range *)
  xMin = Floor[0.9 xFlat[a, b], 0.2];
  xMax = Ceiling[1.1 xSharp[a, b], 0.2];
  yMax = Ceiling[0.7 (xMax - xMin), 0.2];
  (* Absolute plot range margin *)
  margin = 0.1;
  (* Endpoints for constant temperature boundary *)
  yTop = xyTraAsymm["lower"][[2]] @ DomainStart @ xyTraAsymm["lower"];
  yBottom = xyTraAsymm["upper"][[2]] @ DomainEnd @ xyTraAsymm["upper"];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> 8
      , LabelStyle -> 11
    ],
    (* Line of symmetry *)
    Graphics @ {
      Gray, Dashed,
      Line @ {{xMin, 0}, {xMax, 0}}
    },
    (* Straight contour *)
    ParametricPlot[
      {xStraight, y},
      {y, yTop, yBottom}
      , PlotPoints -> 2
      , PlotStyle -> Directive[BoundaryTracingStyle["Contour"], SlidesStyle["Source"]]
    ],
    (* Domain radiation boundaries *)
    Table[
      ParametricPlot[
        xyTraAsymm[id][s] // Through,
        {s, DomainStart @ xyTraAsymm[id], DomainEnd @ xyTraAsymm[id]}
        , PlotPoints -> 2
        , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
      ]
      , {id, {"upper", "lower"}}
    ],
    (* Straight contour and constant temperature labels *)
    Graphics @ {
      SlidesStyle["Source"],
      Text[
        xIt == SeparatedRow["VeryThin"]["\[Pi]", "/", 2] // textStyle
        , {xStraight, Way[yBottom, yTop, 7/10]}
        , {0, 1.3}
        , {0, 1}
      ]
    },
    Graphics @ {
      SlidesStyle["Source"],
      Text[
        Italicise["T"] == 1 // textStyle
        , {xStraight, Way[yBottom, yTop, 3/10]}
        , {0, 1.5}
        , {0, 1}
      ]
    },
    (* Finite element mesh *)
    mesh["Wireframe"],
    {}
    , ImageSize -> 0.5 * 0.8 ImageSizeTextWidthBeamer
  ]
] // Ex["cosine_general-asymmetric_mesh-slides.pdf"];


(* ::Subsection:: *)
(*Legend (cosine_general-asymmetric-construction-legend.pdf)*)


Module[
  {
    legendLabelStyle,
    topRow, bottomRow,
    dummyForTrailingCommas
  },
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  topRow =
    CurveLegend[
      BoundaryTracingStyle /@ {"Background", "Traced", "Contour"},
      {"inflection frontier", "traced boundary", "constant temperature"}
      , LabelStyle -> legendLabelStyle
    ];
  bottomRow = Join[
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ],
    CurveLegend[
      BoundaryTracingStyle /@ {"Terminal"},
      {"terminal curve"}
      , LabelStyle -> legendLabelStyle
    ],
    {}
  ];
  (* Combine *)
  GraphicsGrid[{topRow, bottomRow}
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.115
    , Spacings -> {Automatic, 0}
  ]
] // Ex["cosine_general-asymmetric-construction-legend.pdf"]


(* ::Section:: *)
(*Figure: Verification domains and meshes (cosine-verification-domain-meshes.pdf)*)


Module[
  {
    xHalfWidth, yMax,
    simpleDomain, simpleMeshList, simpleMesh,
    generalDomain, generalMesh,
    legendCurves, legend,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xHalfWidth = 0.45;
  yMax = 0.8;
  (* Simple case: lens-shaped domain *)
  Module[
    {
      a, b,
      xRad, yEnd, xRadEvenExtension,
      xCentre, xMin, xMax,
      nameSuffix, mesh,
      nMeshPortions, yMinMeshPortion, yMaxMeshPortion,
      dummyForTrailingCommas1
    },
    (* Values of A and B *)
    a = aValuesSimpConvex // First;
    b = 1;
    (* Radiation boundary x == x(y) for convex domain *)
    xRad = xTraCandSimp[a, True];
    yEnd = DomainEnd[xRad];
    xRadEvenExtension = Function[{y}, xRad[Abs @ y] // Evaluate];
    (* Plot range *)
    xCentre = Way[x0Simp[a], xStraight];
    xMin = xCentre - xHalfWidth;
    xMax = xCentre + xHalfWidth;
    (* Import mesh *)
    nameSuffix = aNamesSimpConvex[a];
    mesh =
      Import @ FString["cosine_simple-verification-mesh-{nameSuffix}.txt"]
        // Uncompress // #[[2]] &;
    (* Domain plot *)
    simpleDomain =
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax}
          , FrameLabel -> {
              Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
              Italicise["y"]
            }
          , FrameStyle -> LabelSize["Axis"]
          , FrameTicksStyle -> LabelSize["Tick"]
          , LabelStyle -> LatinModernLabelStyle[LabelSize["Tick"] - 1]
          , PlotLabel -> Column[
              {
                "Lens\[Hyphen]shaped",
                SeparatedRow[","][aIt == N[a], bIt == b]
              }
              , Alignment -> Center
            ]
        ],
        (* Boundaries *)
        ParametricPlot[
          {{xRadEvenExtension[y], y}, {xStraight, y}}
          , {y, -yEnd, yEnd}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle /@ {"Traced", "Contour"}
        ],
        {}
      ];
    (* Mesh *)
    nMeshPortions = 3;
    simpleMeshList =
      Table[
        yMinMeshPortion = Way[yEnd, -yEnd, (n - 1) / nMeshPortions];
        yMaxMeshPortion = Way[yEnd, -yEnd, n / nMeshPortions];
        Show[mesh["Wireframe"]
          , PlotRange -> {Automatic, {yMinMeshPortion, yMaxMeshPortion}}
        ]
        , {n, nMeshPortions}
      ];
    simpleMesh =
      GraphicsGrid[
        {simpleMeshList}
        , Spacings -> {{2 -> 0.25 ImageSizeTextWidth, 3 -> 0}, 0}
        (* non-uniform spacing looks better (an optical illusion) *)
      ];
  ];
  (* General case: asymmetric domain *)
  Module[
    {
      a, b,
      yTop, yBottom,
      xCentre, xMin, xMax,
      mesh,
      dummyForTrailingCommas1
    },
    (* Values of A and B *)
    a = aAsymm;
    b = bAsymm;
    (* Endpoints for constant temperature boundary *)
    yTop = xyTraAsymm["lower"][[2]] @ DomainStart @ xyTraAsymm["lower"];
    yBottom = xyTraAsymm["upper"][[2]] @ DomainEnd @ xyTraAsymm["upper"];
    (* Plot range *)
    xCentre = Way[
      xyTraAsymm["upper"][DomainStart @ xyTraAsymm["upper"]] // Through // First,
      xStraight
    ];
    xMin = xCentre - xHalfWidth;
    xMax = xCentre + xHalfWidth;
    (* Import mesh *)
    mesh =
      Import["cosine_general-verification-mesh-asymmetric.txt"]
        // Uncompress // First;
    (* Domain plot *)
    generalDomain =
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax}
          , FrameLabel -> {
              Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
              Italicise["y"]
            }
          , FrameStyle -> LabelSize["Axis"]
          , FrameTicksStyle -> LabelSize["Tick"]
          , LabelStyle -> LatinModernLabelStyle[LabelSize["Tick"] - 1]
          , PlotLabel -> Column[
              {
                "Asymmetric",
                SeparatedRow[","][aIt == a, bIt == SignificantFiguresForm[5][b]]
              }
              , Alignment -> Center
              , Spacings -> 0.27
            ]
        ],
        (* Constant-temperature boundary *)
        ParametricPlot[
          {xStraight, y},
          {y, yTop, yBottom}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Contour"]
        ],
        (* Radiation boundaries *)
        Table[
          ParametricPlot[
            xyTraAsymm[id][s] // Through,
            {s, DomainStart @ xyTraAsymm[id], DomainEnd @ xyTraAsymm[id]}
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ]
          , {id, {"upper", "lower"}}
        ],
        {}
      ];
    (* Mesh *)
    generalMesh = Show[mesh["Wireframe"]];
  ];
  (* Legend *)
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced", "Contour"},
      {"radiation", "constant temperature"}
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
    ];
  legend = Grid[{legendCurves}, Spacings -> {2.5, 0}];
  (* Export *)
  {
    Show[simpleDomain
      , ImageSize -> {Automatic, 0.6 ImageSizeTextWidth}
    ] // Ex["cosine-verification-lens-domain.pdf"]
    ,
    Show[simpleMesh
      , ImageSize -> 0.16 ImageSizeTextWidth
    ] // Ex["cosine-verification-lens-mesh.pdf"]
    ,
    Show[generalDomain
      , ImageSize -> {Automatic, 0.6 ImageSizeTextWidth}
    ] // Ex["cosine-verification-asymmetric-domain.pdf"]
    ,
    Show[generalMesh
      , ImageSize -> 0.11 ImageSizeTextWidth
    ] // Ex["cosine-verification-asymmetric-mesh.pdf"]
    ,
    GraphicsGrid[
      List @ Join[legendCurves]
      , Alignment -> Left
      , ImageSize -> 0.65 ImageSizeTextWidth
      , ItemAspectRatio -> 0.12
      , Spacings -> -0.1 ImageSizeTextWidth
    ] // Ex["cosine-verification-legend.pdf"]
    ,
    Nothing
  }
]


(* ::Subsection:: *)
(*Number of mesh elements*)


(* ::Subsubsection:: *)
(*Simple case: lens-shaped domain*)


Module[{a, nameSuffix, mesh},
  a = aValuesSimpConvex // First;
  nameSuffix = aNamesSimpConvex[a];
  mesh =
    Import @ FString["cosine_simple-verification-mesh-{nameSuffix}.txt"]
      // Uncompress // #[[2]] &;
  mesh
]


(* ::Subsubsection:: *)
(*General case: asymmetric domain*)


Module[{mesh},
  mesh =
    Import["cosine_general-verification-mesh-asymmetric.txt"]
      // Uncompress // First;
  mesh
]


(* ::Subsection:: *)
(*Maximum relative error throughout mesh*)


(* ::Subsubsection:: *)
(*Simple case: lens-shaped domain*)


Module[
  {
    a, b,
    nameSuffix, name, tSol, mesh,
    relError,
    dummyForTrailingCommas1
  },
  (* Values of A and B *)
  a = aValuesSimpConvex // First;
  b = 1;
  (* Import numerical solution *)
  nameSuffix = aNamesSimpConvex[a];
  name = FString @ "cosine_simple-verification-solution-{nameSuffix}.txt";
  tSol = Import[name] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tKnown[b][x, y] - 1;
  relError @@@ mesh["Coordinates"] // Abs // Max // PercentForm
]


(* ::Subsubsection:: *)
(*General case: asymmetric domain*)


Module[
  {
    a, b,
    name, tSol, mesh,
    relError,
    dummyForTrailingCommas1
  },
  (* Values of A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Import numerical solution *)
  name = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[name] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tKnown[b][x, y] - 1;
  relError @@@ mesh["Coordinates"] // Abs // Max // PercentForm
]


(* ::Section:: *)
(*Figure: Verification relative error (cosine-verification-*-relative-error)*)


(* ::Subsection:: *)
(*Simple case: lens-shaped domain*)


Module[
  {
    a, b,
    nameSuffix, name, tSol, mesh,
    xMin, xMax, yMin, yMax,
    meshXDensity, meshYDensity,
    relError,
    dummyForTrailingCommas1
  },
  (* Values of A and B *)
  a = aValuesSimpConvex // First;
  b = 1;
  (* Import numerical solution *)
  nameSuffix = aNamesSimpConvex[a];
  name = FString @ "cosine_simple-verification-solution-{nameSuffix}.txt";
  tSol = Import[name] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Mesh density *)
  {{xMin, xMax}, {yMin, yMax}} = mesh["Bounds"];
  meshXDensity = 2;
  meshYDensity = 0.3 (yMax - yMin) / (xMax - xMin) * meshXDensity // Floor;
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tKnown[b][x, y] - 1;
  (* Make plot *)
  Module[{x, y},
    Plot3D[relError[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{0, 5}, {0, -10}},
          Italicise["y"] // Margined @ {{20, 0}, {0, -20}},
          None
        }
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> {Automatic, Automatic, 0.08}
      , ImageSize -> 0.42 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , Mesh -> {meshXDensity, meshYDensity}
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotLabel -> (Style["Relative error", LabelSize["Axis"] - 1] // Margined @ {{0, 25}, {0, 0}})
      , PlotRange -> Full
      , PlotRangePadding -> {Scaled /@ {0.22, 0.27}, Scaled[0.02], Scaled[0.01]}
      , PlotStyle -> BoundaryTracingStyle["Solution3D"]
      , TicksStyle -> LabelSize["Tick"]
      , ViewPoint -> {+0.8, -4, 0.4}
    ]
  ]
] // Ex["cosine-verification-lens-relative-error.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Subsection:: *)
(*General case: asymmetric domain*)


Module[
  {
    a, b,
    name, tSol, mesh,
    xMin, xMax, yMin, yMax,
    meshXDensity, meshYDensity,
    relError,
    dummyForTrailingCommas1
  },
  (* Values of A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Import numerical solution *)
  name = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[name] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Mesh density *)
  {{xMin, xMax}, {yMin, yMax}} = mesh["Bounds"];
  meshXDensity = 5;
  meshYDensity = 2/3 (yMax - yMin) / (xMax - xMin) * meshXDensity // Floor;
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tKnown[b][x, y] - 1;
  (* Make plot *)
  Module[{x, y},
    Plot3D[relError[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{0, 5}, {0, -10}},
          Italicise["y"] // Margined @ {{20, 0}, {0, -20}},
          None
        }
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> {Automatic, Automatic, 0.21}
      , ImageSize -> 0.42 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , Mesh -> {meshXDensity, meshYDensity}
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotLabel -> (Style["Relative error", LabelSize["Axis"] - 1] // Margined @ {{0, 25}, {0, 0}})
      , PlotRange -> Full
      , PlotRangePadding -> {Scaled /@ {0.15, 0.07}, Scaled[0.01], Scaled[0.05]}
      , PlotStyle -> BoundaryTracingStyle["Solution3D"]
      , TicksStyle -> LabelSize["Tick"]
      , ViewPoint -> {+1, -2.7, 0.7}
    ]
  ]
] // Ex["cosine-verification-asymmetric-relative-error.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Figure: Verification relative error histogram (cosine-verification-*-relative-error-histogram)*)


(* ::Subsection:: *)
(*Simple case: lens-shaped domain*)


Module[
  {
    a, b,
    nameSuffix, name, tSol, mesh,
    relError, relErrorValues,
    dummyForTrailingCommas
  },
  (* Values of A and B *)
  a = aValuesSimpConvex // First;
  b = 1;
  (* Import numerical solution *)
  nameSuffix = aNamesSimpConvex[a];
  name = FString @ "cosine_simple-verification-solution-{nameSuffix}.txt";
  tSol = Import[name] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tKnown[b][x, y] - 1;
  relErrorValues = relError @@@ mesh["Coordinates"] // Abs;
  Histogram[relErrorValues, {"Log", 13}, "LogCount"
    , ChartStyle -> LightGray
    , Frame -> {{True, False}, {True, False}}
    , FrameLabel -> {
        Style["Relative error", Smaller] // Margined @ {{0, 0}, {0, -20}},
        None
      }
    , FrameTicks -> {
        Table[
          {
            10^p,
            Superscript[10, p] // Margined @ {{0, 0}, {0, -7}},
            {0, 0.02}
          }
          , {p, -16, -4, 4}
        ],
        Automatic
      }
    , FrameTicksStyle -> LabelSize["Tick"]
    , ImageSize -> 0.36 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> All
  ]
] // Ex["cosine-verification-lens-relative-error-histogram.pdf"]


(* ::Subsection:: *)
(*General case: asymmetric domain*)


Module[
  {
    a, b,
    name, tSol, mesh,
    relError, relErrorValues,
    dummyForTrailingCommas
  },
  (* Values of A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Import numerical solution *)
  name = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[name] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tKnown[b][x, y] - 1;
  relErrorValues = relError @@@ mesh["Coordinates"] // Abs;
  Histogram[relErrorValues, {"Log", 13}, "LogCount"
    , ChartStyle -> LightGray
    , Frame -> {{True, False}, {True, False}}
    , FrameLabel -> {
        Style["Relative error", Smaller] // Margined @ {{0, 0}, {0, -20}},
        None
      }
    , FrameTicks -> {
        Table[
          {
            10^p,
            Superscript[10, p] // Margined @ {{0, 0}, {0, -7}},
            {0, 0.02}
          }
          , {p, -16, -4, 4}
        ],
        Automatic
      }
    , FrameTicksStyle -> LabelSize["Tick"]
    , ImageSize -> 0.36 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> All
  ]
] // Ex["cosine-verification-asymmetric-relative-error-histogram.pdf"]
