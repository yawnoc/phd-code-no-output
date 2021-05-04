(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "conformal_triangle"};


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Conformal map (disk interior to triangle exterior)*)


(* ::Subsubsection:: *)
(*dz/d\[Zeta]*)


zOfZetaDerivative[zeta_, c_] := c / zeta^2 * (1 - zeta^3)^(2/3);
zOfZetaDerivative[\[FormalZeta], \[FormalCapitalC]]


(* ::Subsubsection:: *)
(*z(\[Zeta])*)


zOfZeta[zeta_, c_] := Integrate[zOfZetaDerivative[zeta, c], zeta] // Evaluate;
zOfZeta[\[FormalZeta], \[FormalCapitalC]]


(* ::Subsubsection:: *)
(*Multiplicative constant*)


constant = 1 / zOfZeta[1, 1];
-1/constant == Hypergeometric2F1[-2/3, -1/3, 2/3, 1]
{constant, constant // N}


(* ::Subsubsection:: *)
(*z(\[Zeta]) with fixed constant*)


zOfZetaDerivative[zeta_] := zOfZetaDerivative[zeta, constant] // Evaluate;
zOfZeta[zeta_] := zOfZeta[zeta, constant] // Evaluate;


xyOfZeta[zeta_] := ReIm @ zOfZeta[zeta] // Evaluate;


(* ::Subsubsection:: *)
(*Principal cubic root of unity*)


omega = Exp[I 2 Pi/3];
omega^3


(* ::Subsubsection:: *)
(*Checks*)


zOfZeta[1] == 1
zOfZeta[omega^2] == omega
zOfZeta[omega] == omega^2
zOfZeta'[\[FormalZeta]] == zOfZetaDerivative[\[FormalZeta]] // FullSimplify
zOfZeta[\[FormalZeta]] ==
  Divide[
    Hypergeometric2F1[-2/3, -1/3, 2/3, \[FormalZeta]^3],
    \[FormalZeta] Hypergeometric2F1[-2/3, -1/3, 2/3, 1]
  ]


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*Analytical function, W(\[Zeta])*)


b = 1.6;
rho0 = Exp[-b];
w[zeta_] := Log[zeta / rho0] // Evaluate;


(* ::Subsubsection:: *)
(*Physical temperature, Re{W}*)


temperature[zeta_] := Re @ w[zeta] // Evaluate;


(* ::Subsection:: *)
(*Boundary condition*)


(* ::Subsubsection:: *)
(*Flux*)


a = 1.5;
flux[zeta_] := -temperature[zeta] / a // Evaluate;


(* ::Subsubsection:: *)
(*Gradient squared*)


gradSquared[zeta_] := Abs[w'[zeta] / zOfZetaDerivative[zeta]] ^ 2 // Evaluate;
gradSquared[\[FormalZeta]] * constant^2


(* ::Subsubsection:: *)
(*Viability*)


viability[zeta_] := gradSquared[zeta] - flux[zeta]^2 // Evaluate;


(* ::Subsection:: *)
(*Hyperbolic critical terminal point*)


(* Principal point *)
viabilityTolerance = 10^-6;
angHyperbolic = Pi/3;
radHyperbolic = SeekRootBisection[
  viability[# Exp[I angHyperbolic]] - viabilityTolerance &,
  {rho0, 1}
  , viabilityTolerance / 1000
];
zetaHyperbolic = radHyperbolic * Exp[I angHyperbolic];


viability[zetaHyperbolic]


(* ::Subsection:: *)
(*Traced boundaries*)


(* ::Subsubsection:: *)
(*Solvers for ODE*)


tracedRHS[zeta_, branchSign_: 1] :=
  Divide[
    I flux[zeta] + branchSign Sqrt[viability[zeta]],
    w'[zeta]
  ] // Evaluate;


sMax = 1.5 zOfZeta[rho0] // Ceiling;


tracedBoundary[zeta0_, branchSign_: 1, terminationViability_: 0] :=
  Module[{zeta},
    NDSolveValue[
      {
        zeta'[s] == tracedRHS[zeta[s], branchSign],
        zeta[0] == zeta0,
        WhenEvent[
          {
            Abs[zeta[s]] < rho0,
            Abs[zeta[s]] > 1,
            viability[zeta[s]] < terminationViability
          },
          "StopIntegration"
        ]
      },
      zeta, {s, -sMax, sMax}
      , NoExtrapolation
    ]
  ];


(* ::Subsubsection:: *)
(*Actual boundaries*)


(* ::Subsubsubsection:: *)
(*Starting points*)


tracedTypeList = {"general", "hyperbolic"};


zetaStartList["general"] =
  Module[{rad, angValues},
    rad = Way[rho0, radHyperbolic];
    angValues = 2 Pi / 3 {0.25, 0.5, 0.75};
    Table[rad Exp[I ang], {ang, angValues}]
  ];


zetaStartList["hyperbolic"] = {zetaHyperbolic};


(* ::Subsubsubsection:: *)
(*Boundaries*)


Table[
  zetaTracedList[type, branchSign] =
    Table[tracedBoundary[zeta0, branchSign], {zeta0, zetaStartList[type]}
    ];
  , {type, tracedTypeList}
  , {branchSign, {-1, 1}}
];


(* ::Subsection:: *)
(*Numerical verification*)


(* ::Subsubsection:: *)
(*Portion of hyperbolic traced boundary*)


zetaTracedVerification = zetaTracedList["hyperbolic", -1] // First;


sVerificationStart = 0;
sVerificationEnd =
  SeekRoot[
    Arg @ zetaTracedVerification[#] &,
    {sVerificationStart, DomainEnd[zetaTracedVerification]}
  ];


(* ::Subsubsection:: *)
(*Finite element mesh*)


Module[
  {
    spacingZeta, spacingZ,
    outerZetaListOneSixth, outerZetaListOneThird, outerZetaList,
    outerZList, outerXYList, outerInradius,
    innerZListOneThird, innerZList, innerXYList, innerCircumradius,
    numPointsOuter, numPointsInner, mod,
    boundaryMesh, mesh,
    discriminatingRadius, predicateOuter, predicateInner,
    dummyForTrailingCommas
  },
  (* Spacing of points in \[Zeta]- and z-space *)
  spacingZeta = 0.1;
  spacingZ = 0.1;
  (* Outer boundary (flux condition) *)
  outerZetaListOneSixth =
    Table[
      zetaTracedVerification[s]
      , {s, UniformRange[sVerificationStart, sVerificationEnd, spacingZeta]}
    ];
  outerZetaListOneThird =
    Join[
      outerZetaListOneSixth,
      outerZetaListOneSixth // Conjugate // Most // Reverse,
      {}
    ] // Most;
  outerZetaList = Join @@ Table[outerZetaListOneThird / omega^k, {k, 0, 2}];
  outerZList = zOfZeta /@ outerZetaList;
  outerXYList = ReIm[outerZList];
  outerInradius = Min @ Abs[outerZList];
  (* Inner boundary (constant temperature) *)
  innerZListOneThird = Subdivide[1, omega, Abs[omega - 1] / spacingZ // Ceiling] // Most;
  innerZList = Join @@ Table[innerZListOneThird * omega^k, {k, 0, 2}];
  innerXYList = ReIm[innerZList];
  innerCircumradius = Max @ Abs[innerZList];
  (* Numbering *)
  numPointsOuter = Length[outerXYList];
  numPointsInner = Length[innerXYList];
  mod[n_] := Mod[#, n, 1] &;
  (* Build boundary mesh *)
  boundaryMesh = ToBoundaryMesh[
    "Coordinates" -> Join[outerXYList, innerXYList],
    "BoundaryElements" -> {
      LineElement[
        Table[{n, n+1}, {n, numPointsOuter}] // mod[numPointsOuter]
      ],
      LineElement[
        numPointsOuter
          +
        (Table[{n, n+1}, {n, numPointsInner}] // mod[numPointsInner])
      ]
    }
  ];
  (* Build mesh *)
  mesh = ToElementMesh[boundaryMesh
    , "ImproveBoundaryPosition" -> True
    , "RegionHoles" -> {0, 0}
  ];
  (* Predicate functions for boundaries *)
  discriminatingRadius = Way[innerCircumradius, outerInradius];
  If[Not[innerCircumradius < discriminatingRadius < outerInradius],
    Print["WARNING: bad discriminatingRadius for boundary predicate functions"];
  ];
  predicateOuter = Function[{x, y}, Abs[x + I y] > discriminatingRadius // Evaluate];
  predicateInner = Function[{x, y}, Abs[x + I y] < discriminatingRadius // Evaluate];
  (* Store in global variables *)
  verificationMesh = mesh;
  verificationPredicateOuter = predicateOuter;
  verificationPredicateInner = predicateInner;
]


(* ::Subsubsection:: *)
(*Solve boundary value problem*)


verificationSolution =
  With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
    NDSolveValue[
      {
        (* Steady state heat equation *)
        -Laplacian[t[x, y], {x, y}] ==
          (* Outer (flux) boundary condition *)
          NeumannValue[-t[x, y] / a, verificationPredicateOuter[x, y]],
        (* Inner (constant temperature) boundary condition *)
        DirichletCondition[t[x, y] == b, verificationPredicateInner[x, y]]
      }, t, Element[{x, y}, verificationMesh]
    ]
  ];


(* ::Section:: *)
(*Visualise transformation*)


Module[
  {
    xMin, xMax, yMin, yMax,
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    dummyForTrailingCommas
  },
  (* Plot range (z-space) *)
  {xMin, xMax} = {yMin, yMax} = Norm[xyOfZeta[rho0]] {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[rho0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {yMin, yMax}
      , ImageSize -> 360
    ],
    (* Radial contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, radValues}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {ang, angValues}]
      , {rad, radMin, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
    ],
    {}
  ]
]


(* ::Section:: *)
(*Visualise known solution*)


Module[
  {
    radMin, radMax,
    angMin, angMax,
    dummyForTrailingCommas
  },
  (* Plot range (\[Zeta]-space) *)
  {radMin, radMax} = {rho0, 1};
  {angMin, angMax} = {0, 2 Pi};
  (* Make plot *)
  ParametricPlot3D[
    Append[
      xyOfZeta[#],
      temperature[#]
    ] & [
      rad Exp[I ang]
    ]
    , {rad, radMin, radMax}
    , {ang, angMin, angMax}
    , BoxRatios -> {Automatic, Automatic, 2.5}
    , Exclusions -> None
  ]
]


(* ::Section:: *)
(*\[Zeta]-space visualisation*)


Module[
  {
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    zetaContourStyle,
    upperStyle, lowerStyle, hyperbolicStyle,
    dummyForTrailingCommas
  },
  (* Contours *)
  radValues = Subdivide[rho0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Styles *)
  zetaContourStyle = Nest[Lighter, Blue, 3];
  hyperbolicStyle = Directive[Thick, Black];
  upperStyle = Blue;
  lowerStyle = Red;
  (* Make plot *)
  Show[
    (* Radial contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {rad, radValues}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> zetaContourStyle
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {ang, angValues}]
      , {rad, radMin, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> zetaContourStyle
    ],
    (* Non-viable domain *)
    RegionPlot[
      viability[xx + I yy] < 0
      , {xx, -2 radMax, 2 radMax}
      , {yy, -2 radMax, 2 radMax}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 50
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Hyperbolic critical terminal points *)
    Graphics @ {
      Black, PointSize[Large],
      Point @ Table[zetaHyperbolic * omega^k // ReIm, {k, 0, 2}]
    },
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          zeta[s] * omega^k // ReIm // Evaluate
          , {s, DomainStart[zeta], DomainEnd[zeta]}
          , PlotStyle -> Which[
              type == "hyperbolic", hyperbolicStyle,
              branchSign == 1, upperStyle,
              branchSign == -1, lowerStyle,
              True, Automatic
            ]
        ]
        , {zeta, zetaTracedList[type, branchSign]}
      ]
      , {type, tracedTypeList}
      , {branchSign, {-1, 1}}
      , {k, 0, 2}
    ],
    {}
  ]
]


(* ::Section:: *)
(*z-space visualisation*)


Module[
  {
    xMin, xMax, yMin, yMax,
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    zetaContourStyle,
    upperStyle, lowerStyle, hyperbolicStyle,
    dummyForTrailingCommas
  },
  (* Plot range (z-space) *)
  {xMin, xMax} = {yMin, yMax} = Norm[xyOfZeta[rho0]] {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[rho0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Styles *)
  zetaContourStyle = Nest[Lighter, Blue, 3];
  hyperbolicStyle = Directive[Thick, Black];
  upperStyle = Blue;
  lowerStyle = Red;
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {yMin, yMax}
      , ImageSize -> 360
    ],
    (* Radial contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, radValues}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> zetaContourStyle
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {ang, angValues}]
      , {rad, radMin, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> zetaContourStyle
    ],
    (* Non-viable domain *)
    (*
      The inverse map \[Zeta](z) is expensive to compute,
      so we instead generate the region in \[Zeta]-space
      and map forward to z-space, which is cheap.
      Credit to Michael E2 for the implementation,
      see <https://mathematica.stackexchange.com/a/85922>.
      Re-use permission granted in comments, see archived version:
      <https://web.archive.org/web/20210407060154/https://mathematica.stackexchange.com/questions/85919/transforming-a-region-obtained-with-regionplot>
    *)
    Module[
      {
        zetaRegion, zetaBoundaryRegion,
        fun, zBoundaryMesh, zMesh, zRegion,
        dummyForTrailingCommas1
      },
      (* \[Zeta]-space *)
      zetaRegion =
        DiscretizeGraphics @ RegionPlot[
          viability[xx + I yy] < 0 && rho0 < Abs[xx + I yy] < 1
          , {xx, -2 radMax, 2 radMax}
          , {yy, -2 radMax, 2 radMax}
        ];
      zetaBoundaryRegion = BoundaryMesh[zetaRegion];
      (* Forward transformation *)
      fun = Function[{xx, yy}, xyOfZeta[xx + I yy]];
      zBoundaryMesh =
        ToBoundaryMesh[
          "Coordinates" -> fun @@@ MeshCoordinates[zetaBoundaryRegion],
          "BoundaryElements" -> {LineElement @@ Thread[MeshCells[zetaBoundaryRegion, 1], Line]}
        ];
      (* z-space *)
      zMesh =
        ToElementMesh[zBoundaryMesh
          , MaxCellMeasure -> {"Area" -> Infinity}
          , "MeshOrder" -> 1
        ];
      zRegion = MeshRegion[zMesh];
      (* Plot *)
      RegionPlot[zRegion
        , BoundaryStyle -> None(*BoundaryTracingStyle["Terminal"]*)
        , PlotStyle -> BoundaryTracingStyle["NonViable"]
      ]
    ],
    (* Hyperbolic critical terminal points *)
    Graphics @ {
      Black, PointSize[Large],
      Point @ Table[zetaHyperbolic * omega^k // xyOfZeta, {k, 0, 2}]
    },
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          zeta[s] * omega^k // xyOfZeta // Evaluate
          , {s, DomainStart[zeta], DomainEnd[zeta]}
          , PlotPoints -> 2
          , PlotStyle -> Which[
              type == "hyperbolic", hyperbolicStyle,
              branchSign == 1, upperStyle,
              branchSign == -1, lowerStyle,
              True, Automatic
            ]
        ]
        , {zeta, zetaTracedList[type, branchSign]}
      ]
      , {type, tracedTypeList}
      , {branchSign, {-1, 1}}
      , {k, 0, 2}
    ],
    {}
  ]
]


(* ::Section:: *)
(*Numerical verification*)


(* ::Subsection:: *)
(*Visualise mesh*)


verificationMesh["Wireframe"]


(* ::Subsection:: *)
(*Comparison direct plot*)


Module[
  {
    radMin, radMax,
    angMin, angMax,
    x, y,
    dummyForTrailingCommas
  },
  (* Plot range (\[Zeta]-space) for exact solution *)
  {radMin, radMax} = {rho0, 1};
  {angMin, angMax} = {0, 2 Pi};
  (* Make plot *)
  Show[
    (* Exact solution *)
    ParametricPlot3D[
      Append[
        xyOfZeta[#],
        temperature[#]
      ] & [
        rad Exp[I ang]
      ]
      , {rad, radMin, radMax}
      , {ang, angMin, angMax}
      , BoxRatios -> {Automatic, Automatic, 2.5}
      , Exclusions -> None
    ],
    (* Numerical solution *)
    Plot3D[
      verificationSolution[x, y], Element[{x, y}, verificationMesh]
      , PlotStyle -> Directive[Opacity[0.7], Blue]
    ],
    {}
  ]
]


(* ::Subsubsection:: *)
(*Relative error on mesh (by z-coordinates)*)


Module[
  {
    radValues, angValues,
    zetaValues, exactSolution,
    relativeError,
    x, y,
    dummyForTrailingCommas
  },
  (* Exact solution on mesh (z-space) *)
  radValues = Subdivide[rho0, 1, 100];
  angValues = Subdivide[0, 2 Pi, 100];
  zetaValues = Join @@ Table[
    rad Exp[I ang]
    , {rad, radValues}
    , {ang, angValues}
  ];
  exactSolution = Interpolation[
    Table[
      {xyOfZeta[zeta], temperature[zeta]}
      , {zeta, zetaValues // DeleteDuplicates}
    ]
    , InterpolationOrder -> 1
  ];
  (* Relative error *)
  relativeError[x_, y_] := verificationSolution[x, y] / exactSolution[x, y] - 1;
  Plot3D[relativeError[x, y], Element[{x, y}, verificationMesh]
    , PlotRange -> Full
  ]
]


(* ::Section:: *)
(*Figure: polar grid (conformal_triangle-grid-*-space)*)


(* ::Subsection:: *)
(*\[Zeta]-space*)


Module[
  {
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    radMarkerRad, angMarkerRad,
    radMarkerAng, angMarkerAng,
    axesLabel, textStyle,
    dummyForTrailingCommas
  },
  (* Contours *)
  radValues = Subdivide[0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Polar coordinate markers *)
  radMarkerRad = radValues[[-2]];
  angMarkerRad = radValues[[-4]];
  radMarkerAng = angMarkerAng = angValues[[2]];
  (* Make plot *)
  axesLabel[string_] := Row @ {string, " ", "\[Zeta]" // LaTeXStyle};
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    EmptyFrame[{-1, 1}, {-1, 1}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -10}},
          axesLabel["Im"] // Margined @ {{0, -4}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
      , PlotRangePadding -> Scaled[0.03]
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {ang, angValues}]
      , {rad, radMin, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Background"]
    ],
    (* Radial contours (\[Rho] < 1) *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {rad, radValues // Most}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Background"]
    ],
    (* Radial contours (\[Rho] == 1) *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {rad, radValues // Last // List}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    (* Radial coordinate marker *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {ang, {radMarkerAng}}]
      , {rad, radMin, radMarkerRad}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> Black
    ] /. {line_Line :> {Arrowheads @ {{0.063, 1}}, Arrow[line]}},
    Graphics @ {
      Text[
        "\[Rho]" // textStyle
        , 0.83 radMarkerRad * Exp[I radMarkerAng] // ReIm
        , {0, -1.15}
      ]
    },
    (* Angular coordinate marker *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {rad, {angMarkerRad}}]
      , {ang, angMin, angMarkerAng}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> Black
    ] /. {line_Line :> {Arrowheads @ {{0.055, 0.98}}, Arrow[line]}},
    Graphics @ {
      Text[
        "\[CurlyPhi]" // textStyle
        , angMarkerRad * Exp[I 3/5 angMarkerAng] // ReIm
        , {-2, 0}
      ]
    },
    {}
    , ImageSize -> 0.47 ImageSizeTextWidth
  ]
] // Ex["conformal_triangle-grid-zeta-space.pdf"]


(* ::Subsection:: *)
(*z-space*)


Module[
  {
    xMin, xMax, yMin, yMax,
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    radMarkerRad, angMarkerRad,
    radMarkerAng, angMarkerAng,
    axesLabel, textStyle,
    dummyForTrailingCommas
  },
  (* Plot range (z-space) *)
  {xMin, xMax} = {yMin, yMax} = Abs[zOfZeta[rho0]] {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[0, 1, 6] /. {0 -> rho0/2};
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Polar coordinate markers *)
  radMarkerRad = radValues[[-2]];
  angMarkerRad = radValues[[-4]];
  radMarkerAng = angMarkerAng = angValues[[2]];
  (* Make plot *)
  axesLabel[string_] := Row @ {string, "\[ThinSpace]", Italicise["z"]};
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    EmptyFrame[{xMin, xMax}, {yMin, yMax}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -15}},
          axesLabel["Im"] // Margined @ {{0, -2}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {ang, angValues}]
      , {rad, radMin/2, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Background"]
    ],
    (* Radial contours (\[Rho] < 1) *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, radValues // Most}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Background"]
    ],
    (* Radial contours (\[Rho] == 1) *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, radValues // Last // List}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    (* Radial coordinate marker *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {ang, {radMarkerAng}}]
      , {rad, radMin, radMarkerRad}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> Black
    ] /. {line_Line :> {Arrowheads @ {{0.065, 0.7}}, Arrow[line]}},
    Graphics @ {
      Text[
        "\[Rho]" // textStyle
        , 0.3 radMarkerRad * Exp[I radMarkerAng] // xyOfZeta
        , {1.2, 0.5}
      ]
    },
    (* Angular coordinate marker *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, {angMarkerRad}}]
      , {ang, angMin, angMarkerAng}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> Black
    ] /. {line_Line :> {Arrowheads @ {{0.06, 0.97}}, Arrow[line]}},
    Graphics @ {
      Text[
        "\[CurlyPhi]" // textStyle
        , angMarkerRad * Exp[I 1/2 angMarkerAng] // xyOfZeta
        , {-2.25, 0}
      ]
    },
    {}
    , ImageSize -> 0.47 ImageSizeTextWidth
  ]
] // Ex["conformal_triangle-grid-z-space.pdf"]


(* ::Section:: *)
(*Figure: known solution (conformal_triangle-known-solution)*)


Module[
  {
    radMin, radMax,
    angMin, angMax,
    axesLabel,
    dummyForTrailingCommas
  },
  (* Plot range (\[Zeta]-space) *)
  {radMin, radMax} = {rho0, 1};
  {angMin, angMax} = {0, 2 Pi};
  (* Make plot *)
  axesLabel[string_] := Row @ {string, "\[ThinSpace]", Italicise["z"]};
  ParametricPlot3D[
    Append[
      xyOfZeta[#],
      temperature[#]
    ] & [
      rad Exp[I ang]
    ]
    , {rad, radMin, radMax}
    , {ang, angMin, angMax}
    , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
    , AxesLabel -> {
        axesLabel["Re"] // Margined @ {{5, 0}, {0, -2}},
        axesLabel["Im"] // Margined @ {{3, 0}, {0, 0}},
        Italicise["T"] // Margined @ {{0, 0}, {0, 30}}
      }
    , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
    , Boxed -> {Back, Bottom, Left}
    , BoxRatios -> {Automatic, Automatic, 3.5}
    , Exclusions -> None
    , ImageSize -> 0.48 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
    , Lighting -> GeneralStyle["AmbientLighting"]
    , Mesh -> {5, 8}
    , MeshStyle -> BoundaryTracingStyle["Edge3D"]
    , PlotPoints -> 50
    , PlotStyle -> Directive[GeneralStyle["Translucent"], BoundaryTracingStyle["Solution3D"]]
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["conformal_triangle-known-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Figure: traced boundaries (conformal_triangle-traced-boundaries-*)*)


(* ::Subsection:: *)
(*\[Zeta]-space*)


Module[
  {
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    axesLabel, textStyle,
    dummyForTrailingCommas
  },
  (* Contours *)
  radValues = Subdivide[0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Make plot *)
  axesLabel[string_] := Row @ {string, " ", "\[Zeta]" // LaTeXStyle};
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    EmptyFrame[{-1, 1}, {-1, 1}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -10}},
          axesLabel["Im"] // Margined @ {{0, -4}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
      , PlotRangePadding -> Scaled[0.03]
    ],
    (* Non-viable domain *)
    RegionPlot[
      viability[xx + I yy] < 0
        && Abs[xx + I yy] < 1
        && temperature[xx + I yy] > 0
      , {xx, -radMax, radMax}
      , {yy, -radMax, radMax}
      , BoundaryStyle -> None
      , PlotPoints -> 15
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Radial contour \[Rho] == 1 *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {rad, radValues // Last // List}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          zeta[s] * omega^k // ReIm // Evaluate
          , {s, DomainStart[zeta], DomainEnd[zeta]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
        , {zeta, zetaTracedList[type, branchSign]}
      ]
      , {type, tracedTypeList}
      , {branchSign, {-1, 1}}
      , {k, 0, 2}
    ],
    (* Unphysical region *)
    Graphics @ {BoundaryTracingStyle["Unphysical"],
      Disk[{0, 0}, rho0]
    },
    {}
    , ImageSize -> 0.47 ImageSizeTextWidth
  ]
] // Ex["conformal_triangle-traced-boundaries-zeta-space.pdf"]


(* ::Subsection:: *)
(*z-space*)


Module[
  {
    xMin, xMax, yMin, yMax,
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    axesLabel, textStyle,
      rHyperbolic,
      more,
    dummyForTrailingCommas
  },
  (* Plot range (z-space) *)
  {xMin, xMax} = {yMin, yMax} = Abs[zOfZeta[rho0]] {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[0, 1, 6] /. {0 -> rho0/2};
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Make plot *)
  axesLabel[string_] := Row @ {string, "\[ThinSpace]", Italicise["z"]};
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    EmptyFrame[{xMin, xMax}, {yMin, yMax}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -15}},
          axesLabel["Im"] // Margined @ {{0, -2}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
    ],
    (* Non-viable domain *)
    (*
      The inverse map \[Zeta](z) is expensive to compute,
      so we instead generate the region in \[Zeta]-space
      and map forward to z-space, which is cheap.
      Credit to Michael E2 for the implementation,
      see <https://mathematica.stackexchange.com/a/85922>.
      Re-use permission granted in comments, see archived version:
      <https://web.archive.org/web/20210407060154/https://mathematica.stackexchange.com/questions/85919/transforming-a-region-obtained-with-regionplot>
    *)
    Module[
      {
        zetaRegion, zetaBoundaryRegion,
        fun, zBoundaryMesh, zMesh, zRegion,
        dummyForTrailingCommas1
      },
      (* \[Zeta]-space *)
      zetaRegion =
        DiscretizeGraphics @ RegionPlot[
          viability[xx + I yy] < 0
            && Abs[xx + I yy] < 1
            && temperature[xx + I yy] > 0
          , {xx, -radMax, radMax}
          , {yy, -radMax, radMax}
        ];
      zetaBoundaryRegion = BoundaryMesh[zetaRegion];
      (* Forward transformation *)
      fun = Function[{xx, yy}, xyOfZeta[xx + I yy]];
      zBoundaryMesh =
        ToBoundaryMesh[
          "Coordinates" -> fun @@@ MeshCoordinates[zetaBoundaryRegion],
          "BoundaryElements" -> {LineElement @@ Thread[MeshCells[zetaBoundaryRegion, 1], Line]}
        ];
      (* z-space *)
      zMesh =
        ToElementMesh[zBoundaryMesh
          , MaxCellMeasure -> {"Area" -> Infinity}
          , "MeshOrder" -> 1
        ];
      zRegion = DiscretizeRegion[MeshRegion[zMesh], MaxCellMeasure -> 1];
      (* Plot *)
      rHyperbolic = Abs @ zOfZeta[radHyperbolic];
      Quiet[
        RegionPlot[RegionMember[zRegion, {x, y}]
          , {x, -rHyperbolic, rHyperbolic}
          , {y, -rHyperbolic, rHyperbolic}
          , BoundaryStyle -> None
          , PlotStyle -> Directive[BoundaryTracingStyle["NonViable"], GrayLevel[0.7]]
        ]
        , {ImplicitRegion::bcond}
          (*
            Using RegionMember raises the warning ImplicitRegion::bcond,
            but reduces the file size by over 600 kB.
          *)
      ]
    ],
    (* Radial contour \[Rho] == 1 (which is the triangle in z-space) *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, radValues // Last // List}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          zeta[s] * omega^k // xyOfZeta // Evaluate
          , {s, DomainStart[zeta], DomainEnd[zeta]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
        , {zeta, zetaTracedList[type, branchSign]}
      ]
      , {type, tracedTypeList}
      , {branchSign, {-1, 1}}
      , {k, 0, 2}
    ],
    (* Unphysical region *)
    (*
      We cheat and pretend the boundary is circular.
      Close enough for practical purposes.
    *)
    more = 1.2;
    RegionPlot[
      Abs[x + I y] > Abs[zOfZeta[rho0]]
      , {x, more * xMin, more * xMax}
      , {y, more * yMin, more * yMax}
      , BoundaryStyle -> None
      , PlotStyle -> BoundaryTracingStyle["Unphysical"]
      , PlotPoints -> 5
    ],
    {}
    , ImageSize -> 0.47 ImageSizeTextWidth
  ]
] // Ex["conformal_triangle-traced-boundaries-z-space.pdf"]


(* ::Subsection:: *)
(*Legend*)


Module[{legendLabelStyle},
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  GraphicsGrid[
    Transpose @ {
      CurveLegend[
        {BoundaryTracingStyle["Traced"]},
        {"traced boundary"}
        , LabelStyle -> legendLabelStyle
      ],
      RegionLegend[
        {BoundaryTracingStyle["NonViable"]},
        {"non\[Hyphen]viable domain"}
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
  ]
] // Ex["conformal_triangle-traced-boundaries-legend.pdf"]


(* ::Section:: *)
(*Figure: domain (conformal_triangle-domain)*)


Module[
  {
    xMin, xMax, yMin, yMax,
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    axesLabel,
    dummyForTrailingCommas
  },
  (* Plot range (z-space) *)
  {xMin, xMax} = {yMin, yMax} =
    1.15 Abs[zOfZeta @ zetaTracedVerification[sVerificationEnd]] {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[0, 1, 6] /. {0 -> rho0/2};
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Make plot *)
  axesLabel[string_] := Row @ {string, "\[ThinSpace]", Italicise["z"]};
  Show[
    EmptyFrame[{xMin, xMax}, {yMin, yMax}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -15}},
          axesLabel["Im"] // Margined @ {{0, -1}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
    ],
    (* Radial contour \[Rho] == 1 (which is the triangle in z-space) *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, radValues // Last // List}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          reflect[zeta[s]] * omega^k // xyOfZeta // Evaluate
          , {s, sVerificationStart, sVerificationEnd}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
        , {zeta, {zetaTracedVerification}}
        , {reflect, {Identity, Conjugate}}
      ]
      , {k, 0, 2}
    ],
    {}
    , ImageSize -> 0.44 ImageSizeTextWidth
  ]
] // Ex["conformal_triangle-domain.pdf"]
