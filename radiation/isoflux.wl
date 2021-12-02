(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "isoflux"};


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Conformal mapping*)


zOfZeta[zeta_] := zeta^(1/2);
zetaOfZ[z_] := z^2;


(* ::Subsection:: *)
(*Analytic function W*)


wOfZeta[zeta_] := zeta;


(* ::Subsection:: *)
(*Known solution T = Re[W]*)


tOfZeta[zeta_] := Re[wOfZeta[zeta]] // Evaluate;
tOfZ[z_] := tOfZeta[zetaOfZ[z]] // Evaluate;
tOfXY[x_, y_] := tOfZ[x + I y] // ComplexExpand // Evaluate;


(* ::Subsection:: *)
(*Slope dW/dz = dW/d\[Zeta] . d\[Zeta]/dz*)


wSlopeOfZ[z_] := wOfZeta'[zetaOfZ[z]] zetaOfZ'[z] // Evaluate;


(* ::Subsection:: *)
(*Flux function f*)


f = 1;


(* ::Subsection:: *)
(*Viability function \[CapitalPhi]*)


phiOfZ[z_] := Abs[wSlopeOfZ[z]]^2 - f // Evaluate;
phiOfZeta[zeta_] := phiOfZ[zOfZeta[zeta]] // Evaluate;


(* ::Subsection:: *)
(*Solution-positive region function*)


solutionPositiveRegionFunction =
  Function[{x, y}, tOfZ[x + I y] > 0 // Evaluate];


(* ::Subsection:: *)
(*Terminal point of given phase*)


zTerminal[phase_] := Module[{z}, z /. First @ Solve[phiOfZ[z] == 0 && Arg[z] == phase, z, Complexes]];
zetaTerminal[phase_] := zetaOfZ @ zTerminal[phase];


(* ::Subsection:: *)
(*Critical terminal point*)


zCritical = Module[{z}, z /. First @ Solve[phiOfZ[z] == 0 && z > 0, z, Reals]];
zetaCritical = zetaOfZ[zCritical];


(* ::Subsection:: *)
(*Traced boundary derivative d\[Zeta]/ds*)


zetaTracedDerivative[zeta_, sign_] :=
  Divide[
    I * f + (sign) Sqrt[phiOfZeta[zeta]],
    wOfZeta'[zeta]
  ];


(* ::Subsection:: *)
(*Traced boundary \[Zeta](s)*)


zetaTraced[zeta0_, {sStart_, sEnd_}, sign_, terminationPhi_: 0] :=
  Module[{zeta},
    NDSolveValue[
      {
        zeta'[s] == zetaTracedDerivative[zeta[s], sign],
        zeta[0] == zeta0,
        WhenEvent[phiOfZeta[zeta[s]] < terminationPhi, "StopIntegration"],
        WhenEvent[tOfZeta[zeta[s]] < 0, "StopIntegration"],
        Nothing
      }
      , zeta, {s, sStart, sEnd}
    ]
  ]


(* ::Subsubsection:: *)
(*(from walls)*)


yWallList = Join[
  Subdivide[1/2, 1, 5],
  Subdivide[1/2, 1, 5] // Minus // Reverse,
  {}
];


zetaTracedWallListUpper =
  Module[{zeta0},
    Table[
      zeta0 = zetaOfZ[Abs[y] + I y];
      zetaTraced[zeta0, {0, 3}, +1]
      , {y, yWallList // Most}
    ]
  ];


zetaTracedWallListLower =
  Module[{zeta0},
    Table[
      zeta0 = zetaOfZ[Abs[y] + I y];
      zetaTraced[zeta0, {0, -3}, -1]
      , {y, yWallList // Rest}
    ]
  ];


(* ::Subsubsection:: *)
(*(from terminal curve)*)


phiTracedTerminalList = Pi/4 {-1, -3/4, -1/2, 1/2, 3/4, 1};


zetaTracedTerminalListUpper =
  Module[{zeta0},
    Table[
      zeta0 = zetaTerminal[phase];
      zetaTraced[zeta0, Sign[phase] {0, 3}, +1]
      , {phase, phiTracedTerminalList // Rest}
    ]
  ];


zetaTracedTerminalListLower =
  Module[{zeta0},
    Table[
      zeta0 = zetaTerminal[phase];
      zetaTraced[zeta0, Sign[phase] {0, 3}, -1]
      , {phase, phiTracedTerminalList // Most}
    ]
  ];


(* ::Subsubsection:: *)
(*(from axis)*)


zetaAxisPoint = zetaOfZ[7/10];


zetaTracedAxisUpper = zetaTraced[zetaAxisPoint, {-3, 3}, +1];
zetaTracedAxisLower = zetaTraced[zetaAxisPoint, {-3, 3}, -1];


(* ::Subsubsection:: *)
(*(through critical terminal points)*)


zetaTracedCriticalUpper = zetaTraced[zetaCritical, {-3, 3}, +1];
zetaTracedCriticalLower = zetaTraced[zetaCritical, {-3, 3}, -1];


(* ::Subsection:: *)
(*Styles*)


upperStyle = Blue;
lowerStyle = Red;


criticalStyle = Directive[Yellow, PointSize[Large]];


(* ::Section:: *)
(*Checks*)


zOfZeta[\[FormalZeta]]
zetaOfZ[\[FormalZ]]


tOfZeta[\[FormalZeta]]
tOfZ[\[FormalZ]]


wSlopeOfZ[\[FormalZ]]


f


phiOfZ[\[FormalZ]]
phiOfZeta[\[FormalZeta]]


zCritical
zetaCritical


zTerminal[0] == zCritical


zetaTracedDerivative[\[FormalZeta], +1]
zetaTracedDerivative[\[FormalZeta], -1]


(* ::Section:: *)
(*Computer algebra*)


(* ::Subsection:: *)
(*Initial \[Zeta] != 1/4*)


With[{s = \[FormalS], a0 = \[FormalA]},
  Module[
    {
      nMax, a, x, lhs, rhs,
      dummyForTrailingCommas
    },
    nMax = 8;
    a[0] = a0;
    x = Sum[a[n] s^n, {n, 0, nMax}] + O[s]^(nMax + 1);
    lhs = D[x, s];
    rhs = Sqrt[4 Sqrt[x^2 + s^2] - 1];
    Table[
      {
        Subscript[Italicise["a"], n],
        a[n] =
          a[n] /. First @ Solve[
            SeriesCoefficient[lhs - rhs, {s, 0, n - 1}] == 0
              // FullSimplify[#, a0 > 0] &
              // Evaluate
            , a[n]
          ]
      }
      , {n, nMax}
    ]
  ] // ReplaceAll[a0 -> Subscript[Italicise["a"], 0]]
] // TableForm // Ex["isoflux-ordinary-series.pdf"]


(* ::Subsection:: *)
(*Initial \[Zeta] = 1/4 and s > 0*)


With[{s = \[FormalS], a0 = \[FormalA]},
  Module[
    {
      nMax, a, x, lhs, rhs,
      dummyForTrailingCommas
    },
    nMax = 16;
    a[0] = a0;
    x = 1/4 + 2 s^2 + Sum[a[n] s^n, {n, 4, nMax, 2}] + O[s]^(nMax + 2);
    lhs = D[x, s];
    rhs = Sqrt[4 Sqrt[x^2 + s^2] - 1];
    Table[
      {
        Subscript[Italicise["a"], n],
        a[n] =
          a[n] /. First @ Solve[
            SeriesCoefficient[lhs - rhs, {s, 0, n - 1}] == 0
              // FullSimplify[#, s > 0] &
              // Evaluate
            , a[n]
          ]
      }
      , {n, 4, nMax, 2}
    ]
  ] // ReplaceAll[a0 -> Subscript[Italicise["a"], 0]]
] // TableForm // Ex["isoflux-special-series-positive.pdf"]


(* ::Section:: *)
(*Visualisation*)


(* ::Subsection:: *)
(*Known solution T*)


Plot3D[
  tOfZ[x + I y] // Evaluate
  , {x, 0, 3}
  , {y, -3, 3}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
]


(* ::Subsection:: *)
(*Viability function \[CapitalPhi]*)


Plot3D[
  phiOfZ[x + I y] // Evaluate
  , {x, 0, 2}
  , {y, -2, 2}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
  , RegionFunction -> solutionPositiveRegionFunction
]


(* ::Section:: *)
(*Boundary tracing*)


(* ::Subsection:: *)
(*\[Zeta]-space*)


Show[
  EmptyFrame[
    {0, 1}, {-1, 1}
    , FrameLabel -> ({"Re", "Im"}["\[Zeta]"] // Through)
  ],
  (* Non-viable domain *)
  RegionPlot[
    phiOfZeta[xx + I yy] < 0 && tOfZeta[xx + I yy] > 0
    , {xx, 0, 1}
    , {yy, -1, 1}
    , BoundaryStyle -> None
    , PlotStyle -> BoundaryTracingStyle["NonViable"]
  ],
  (* Terminal curve *)
  ContourPlot[
    phiOfZeta[xx + I yy] == 0
    , {xx, 0, 1}
    , {yy, -1, 1}
    , ContourLabels -> None
    , ContourStyle -> BoundaryTracingStyle["Terminal"]
  ],
  (* Solution contours *)
  ContourPlot[
    tOfZeta[xx + I yy]
    , {xx, 0, 1}
    , {yy, -1, 1}
    , Contours ->
        Module[
          {
            tCriticalValue,
            divisionsUntoCritical,
            divisionsBeyondCritical,
            tSingleDivision,
            dummyForTrailingCommas
          },
          tCriticalValue = tOfZ[zCritical];
          divisionsUntoCritical = 2;
          divisionsBeyondCritical = 11;
          tSingleDivision = tCriticalValue / divisionsUntoCritical;
          tSingleDivision * Range[0, divisionsUntoCritical + divisionsBeyondCritical]
        ]
    , ContourLabels -> None
    , ContourShading -> None
  ],
  (* Traced boundaries (from walls) *)
  Table[
    ParametricPlot[
      zeta[s] // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> upperStyle
    ]
    , {zeta, zetaTracedWallListUpper}
  ],
  Table[
    ParametricPlot[
      zeta[s] // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> lowerStyle
    ]
    , {zeta, zetaTracedWallListLower}
  ],
  (* Traced boundaries (from terminal curve) *)
  Table[
    ParametricPlot[
      zeta[s] // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> upperStyle
    ]
    , {zeta, zetaTracedTerminalListUpper}
  ],
  Table[
    ParametricPlot[
      zeta[s] // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> lowerStyle
    ]
    , {zeta, zetaTracedTerminalListLower}
  ],
  (* Traced boundaries (from axis) *)
  Table[
    ParametricPlot[
      zeta[s] // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> If[zeta === zetaTracedAxisUpper, upperStyle, lowerStyle]
    ]
    , {zeta, {zetaTracedAxisUpper, zetaTracedAxisLower}}
  ],
  (* Traced boundaries (through critical terminal point) *)
  Table[
    ParametricPlot[
      zeta[s] // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> If[zeta === zetaTracedCriticalUpper, upperStyle, lowerStyle]
    ]
    , {zeta, {zetaTracedCriticalUpper, zetaTracedCriticalLower}}
  ],
  (* Critical terminal point *)
  Graphics @ {criticalStyle,
    Point[zetaCritical // ReIm]
  },
  {}
] // Ex["isoflux-traced-boundaries-zeta-space.pdf"]


(* ::Subsection:: *)
(*z-space*)


Show[
  EmptyFrame[
    {0, 1.3}, {-1, 1}
  ],
  (* Non-viable domain *)
  RegionPlot[
    phiOfZ[x + I y] < 0 && tOfZ[x + I y] > 0
    , {x, 0, 1}
    , {y, -1, 1}
    , BoundaryStyle -> None
    , PlotStyle -> BoundaryTracingStyle["NonViable"]
  ],
  (* Terminal curve *)
  ContourPlot[
    phiOfZ[x + I y] == 0
    , {x, 0, 1}
    , {y, -1, 1}
    , ContourLabels -> None
    , ContourStyle -> BoundaryTracingStyle["Terminal"]
    , RegionFunction -> solutionPositiveRegionFunction
  ],
  (* Solution contours *)
  ContourPlot[
    tOfZ[x + I y]
    , {x, 0, 1.5}
    , {y, -1, 1}
    , Contours ->
        Module[
          {
            tCriticalValue,
            divisionsUntoCritical,
            divisionsBeyondCritical,
            tSingleDivision,
            dummyForTrailingCommas
          },
          tCriticalValue = tOfZ[zCritical];
          divisionsUntoCritical = 2;
          divisionsBeyondCritical = 11;
          tSingleDivision = tCriticalValue / divisionsUntoCritical;
          tSingleDivision * Range[0, divisionsUntoCritical + divisionsBeyondCritical]
        ]
    , ContourLabels -> None
    , ContourShading -> None
  ],
  (* Traced boundaries (from walls) *)
  Table[
    ParametricPlot[
      zeta[s] // zOfZeta // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> upperStyle
    ]
    , {zeta, zetaTracedWallListUpper}
  ],
  Table[
    ParametricPlot[
      zeta[s] // zOfZeta // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> lowerStyle
    ]
    , {zeta, zetaTracedWallListLower}
  ],
  (* Traced boundaries (from terminal curve) *)
  Table[
    ParametricPlot[
      zeta[s] // zOfZeta // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> upperStyle
    ]
    , {zeta, zetaTracedTerminalListUpper}
  ],
  Table[
    ParametricPlot[
      zeta[s] // zOfZeta // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> lowerStyle
    ]
    , {zeta, zetaTracedTerminalListLower}
  ],
  (* Traced boundaries (from axis) *)
  Table[
    ParametricPlot[
      zeta[s] // zOfZeta // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> If[zeta === zetaTracedAxisUpper, upperStyle, lowerStyle]
    ]
    , {zeta, {zetaTracedAxisUpper, zetaTracedAxisLower}}
  ],
  (* Traced boundaries (through critical terminal point) *)
  Table[
    ParametricPlot[
      zeta[s] // zOfZeta // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> If[zeta === zetaTracedCriticalUpper, upperStyle, lowerStyle]
    ]
    , {zeta, {zetaTracedCriticalUpper, zetaTracedCriticalLower}}
  ],
  (* Critical terminal point *)
  Graphics @ {criticalStyle,
    Point[zCritical // ReIm]
  },
  {}
] // Ex["isoflux-traced-boundaries-z-space.pdf"]


(* ::Section:: *)
(*Sanity check*)


With[{x = \[FormalX], y = \[FormalY]},
  {
    {"Known solution T", tOfXY[x, y]},
    {"T along positive diagonal", tOfXY[x, x]},
    {"T along negative diagonal", tOfXY[x, -x]},
    {"Flux along vertical line", {1, 0}.Grad[tOfXY[x, y], {x, y}] /. {x -> Re[zCritical]}},
    Nothing
  }
] // TableForm
