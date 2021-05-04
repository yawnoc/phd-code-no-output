(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< LaplaceYoung`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "half_plane"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Contact angles*)


(* gpd stands for "gamma per degree". *)
gpdValues = {1, 5, 10, 15, 30, 45, 60, 75};


(* ::Subsection:: *)
(*Finite element mesh*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["half_plane-mesh.txt",
  Module[
   {xMax, yMax, ySpacing, yValues,
    bPointList, nB, mod,
    bMesh, mesh,
    prWet
   },
    (* Length scales *)
    xMax = 10;
    yMax = 1/2;
    ySpacing = 0.005;
    (* Boundary points *)
    yValues = UniformRange[yMax, -yMax, -ySpacing];
    bPointList = Join[
      (* Fine length scale along x == 0 *)
      Table[{0, y}, {y, yValues}],
      (* Default length scale elsewhere *)
      {{xMax, -yMax}, {xMax, yMax}}
    ];
    nB = Length[bPointList];
    mod[n_] := Mod[#, n, 1] &;
    (* Boundary element mesh *)
    bMesh = ToBoundaryMesh[
      "Coordinates" -> bPointList,
      "BoundaryElements" -> {
        LineElement[
          Table[{n, n + 1}, {n, nB}] // mod[nB]
        ]
      }
    ];
    (* Build mesh *)
    mesh = ToElementMesh[bMesh,
      "ImproveBoundaryPosition" -> True
    ];
    (* Predicate function for wetting boundary (x == 0) *)
    prWet = Function[{x, y}, x == 0 // Evaluate];
    {mesh, prWet}
      // Compress
  ]
]


(* ::Subsection:: *)
(*Solve boundary value problem*)


(* (This is slow (~7 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["half_plane-solution.txt",
  Module[{mesh, prWet, gamma},
    (* Import mesh *)
    {mesh, prWet} = Import["half_plane-mesh.txt"] // Uncompress;
    (* Solve Laplace--Young equation *)
    Association @ Table[
      gamma = gpd * Degree;
      gpd -> SolveLaplaceYoung[gamma, mesh, prWet]
    , {gpd, gpdValues}]
      // Compress
  ]
]


(* ::Subsection:: *)
(*Solve PDE (fixed-point iteration)*)


(* (This is extremely slow (~10 min), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["half_plane-solution-fixed_point.txt",
  Module[{mesh, prWet, gamma},
    (* Import mesh *)
    {mesh, prWet} = Import["half_plane-mesh.txt"] // Uncompress;
    (* Solve Laplace--Young equation *)
    Association @ Table[
      gamma = gpd * Degree;
      gpd -> SolveLaplaceYoungFixedPoint[gamma, mesh, prWet]
    , {gpd, gpdValues}]
      // Compress
  ]
]


(* ::Subsection:: *)
(*Italicise symbols*)


gIt = Style["\[Gamma]"];


(* ::Subsection:: *)
(*Global styles for plots*)


pointStyle = PointSize[Large];


(* ::Section:: *)
(*Finite element mesh*)


Module[{mesh, prWet, nElem},
  (* Import mesh *)
  {mesh, prWet} = Import["half_plane-mesh.txt"] // Uncompress;
  nElem = Length @ mesh[[2, 1, 1]];
  (* Plot *)
  Show[
    mesh["Wireframe"],
    PlotLabel -> FString @ "{nElem} mesh elements",
    ImageSize -> 1440,
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["half_plane-mesh.pdf"]


(* ::Section:: *)
(*Fine length scale*)


(* ::Text:: *)
(*Test what fine length scale to use along x = 0*)
(*by benchmarking against the exact half-plane solution for gamma = 1 degree.*)


(* (This is slow (~20 sec).) *)
Module[
  {
    gamma,
    hTheory,
    xMax, yMax,
    mod,
    prWet,
    ySpacingValues,
    yValues, boundaryPointList, numPoints,
    boundaryMesh, mesh,
    tSol, hNonlinear, relErrorNonlinear,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  gamma = 1 Degree;
  (* Theoretical height rise *)
  hTheory = HHalfPlane[gamma];
  (* Length scales *)
  xMax = 10;
  yMax = 1/2;
  (* Modulo *)
  mod[n_] := Mod[#, n, 1] &;
  (* Predicate function for wetting boundary (x == 0) *)
  prWet = Function[{x, y}, x == 0];
  (* Fine length scales to test *)
  ySpacingValues = {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1};
  (* For each fine length scale *)
  Table[
    (* Boundary points *)
    yValues = UniformRange[yMax, -yMax, -ySpacing];
    boundaryPointList = Join[
      (* Fine length scale along x == 0 *)
      Table[{0, y}, {y, yValues}],
      (* Default length scale elsewhere *)
      {{xMax, -yMax}, {xMax, yMax}},
      (* Dummy for trailing commas *)
      {}
    ];
    numPoints = Length[boundaryPointList];
    (* Boundary element mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> boundaryPointList,
        "BoundaryElements" -> {
          LineElement[
            Table[{n, n + 1}, {n, numPoints}] // mod[numPoints]
          ]
        },
        {}
      ];
    (* Build mesh *)
    mesh =
      ToElementMesh[
        boundaryMesh
        , "ImproveBoundaryPosition" -> True
    ];
    (* Compute wall height *)
    tSol = SolveLaplaceYoung[gamma, mesh, prWet];
    hNonlinear = tSol[0, 0];
    (* Compare with theoretical height *)
    relErrorNonlinear = hNonlinear / hTheory - 1;
    (* Return table row *)
    {
      ySpacing,
      Length @ mesh[[2, 1, 1]],
      hNonlinear,
      relErrorNonlinear,
      Nothing
    }
    , {ySpacing, ySpacingValues}
  ]
    //
      TableForm[#
        , TableHeadings -> {
            None,
            {
              "Fine length",
              "Mesh elements",
              "Computed height",
              "rel. error",
              Nothing
            }
          }
      ] &
    //
      Column @ {
        {"Contact angle", gamma},
        {"Theoretical height", hTheory // N},
        #
      } &
    //
      Ex["half_plane-fine-length-scale.pdf"]
]


(* ::Section:: *)
(*Numerical solution*)


(* ::Subsection:: *)
(*Using built-in nonlinear solver*)


Module[{tSolAss, tSol, mesh, gamma},
  tSolAss = Import["half_plane-solution.txt"] // Uncompress;
  Table[
    tSol = tSolAss[gpd];
    mesh = tSol["ElementMesh"];
    gamma = gpd * Degree;
    Block[{x = \[FormalX], y = \[FormalY]},
    (* (Using With results in SetDelayed::wrsym and an empty plot) *)
      Plot3D[tSol[x, y], Element[{x, y}, mesh],
        AxesLabel -> Italicise /@ {"x", "y", "T"},
        PlotLabel -> Column[
          {
            "Numerical solution",
            gIt == gamma
          },
          Alignment -> Center
        ],
        PlotRange -> {0, Sqrt[2]},
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex @ FString @ "half_plane-solution-gpd-{gpd}.png"
  , {gpd, gpdValues}]
]


(* ::Subsection:: *)
(*Using fixed-point iteration*)


Module[{tSolAss, tSol, n, mesh, gamma},
  tSolAss = Import["half_plane-solution-fixed_point.txt"] // Uncompress;
  Table[
    {tSol, n} = tSolAss[gpd];
    mesh = tSol["ElementMesh"];
    gamma = gpd * Degree;
    Block[{x = \[FormalX], y = \[FormalY]},
    (* (Using With results in SetDelayed::wrsym and an empty plot) *)
      Plot3D[tSol[x, y], Element[{x, y}, mesh],
        AxesLabel -> Italicise /@ {"x", "y", "T"},
        PlotLabel -> Column[
          {
            "Numerical solution: {n} iterations" // FString,
            gIt == gamma
          },
          Alignment -> Center
        ],
        PlotRange -> {0, Sqrt[2]},
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex @ FString @ "half_plane-solution-fixed_point-gpd-{gpd}.png"
  , {gpd, gpdValues}]
]


(* ::Subsection:: *)
(*Discrepancy between nonlinear solver and fixed-point iteration*)


Module[
 {tSolAssNonlinear, tSolAssFixedPoint,
  tSolNonlinear, tSolFixedPoint,
  gamma, xMax
 },
  tSolAssNonlinear = Import["half_plane-solution.txt"] // Uncompress;
  tSolAssFixedPoint = Import["half_plane-solution-fixed_point.txt"] // Uncompress;
  Table[
    tSolNonlinear = tSolAssNonlinear[gpd];
    tSolFixedPoint = tSolAssFixedPoint[gpd] // First;
    gamma = gpd * Degree;
    xMax = 3;
    Plot[tSolNonlinear[x, 0] - tSolFixedPoint[x, 0], {x, 0, xMax},
      AxesLabel -> {Italicise @ "x", ""},
      PlotLabel -> Column[
        {
          "[Nonlinear]" - "[Fixed-point]",
          gIt == gamma
        },
        Alignment -> Center
      ],
      PlotRange -> Full,
      PlotStyle -> Red,
      PlotOptions[Axes] // Evaluate
    ] // Ex @ FString @ "half_plane-solution-discrepancy-gpd-{gpd}.pdf"
  , {gpd, gpdValues}]
]


(* ::Subsection:: *)
(*Compare with theoretical height rise values*)


Module[
 {tSolAssNonlinear, tSolAssFixedPoint,
  tSolNonlinear, tSolFixedPoint, n,
  gamma, hTheory,
  hNonlinear, relErrorNonlinear,
  hFixedPoint, relErrorFixedPoint
 },
  tSolAssNonlinear = Import["half_plane-solution.txt"] // Uncompress;
  tSolAssFixedPoint = Import["half_plane-solution-fixed_point.txt"] // Uncompress;
  Table[
    tSolNonlinear = tSolAssNonlinear[gpd];
    {tSolFixedPoint, n} = tSolAssFixedPoint[gpd];
    (* Theoretical *)
    gamma = gpd * Degree;
    hTheory = HHalfPlane[gamma];
    (* Nonlinear solver *)
    hNonlinear = tSolNonlinear[0, 0];
    relErrorNonlinear = hNonlinear / hTheory - 1;
    (* Fixed-point iteration *)
    hFixedPoint = tSolFixedPoint[0, 0];
    relErrorFixedPoint = hFixedPoint / hTheory - 1;
    {
      gamma, hTheory // N,
      hNonlinear, relErrorNonlinear,
      hFixedPoint, relErrorFixedPoint,
      Abs[relErrorNonlinear] < Abs[relErrorFixedPoint]
    }
  , {gpd, gpdValues}]
] // TableForm[#,
  TableHeadings -> {
    None,
    {
      "\[Gamma]", "Theory",
      "Nonlinear", "rel. error",
      "Fixed-point", "rel. error",
      "Nonlinear better"
    }
  }
] & // Ex["half_plane-height-comparison.pdf"]


(* ::Subsection:: *)
(*Numerical error for nonlinear solver (\[Gamma] = 1\[Degree])*)


Module[
 {gpd, gamma,
  tSolAss, tSol,
  xFun, xMax, tMax, tMin, tStep, tInterp
 },
  (* Contact angle *)
  gpd = 1;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tSolAss = Import["half_plane-solution.txt"] // Uncompress;
  tSol = tSolAss[gpd];
  (* Implicit exact half-plane solution x == x(T) *)
  xFun[t_] := XHalfPlane[gamma][t] // Evaluate;
  (* Interpolation to obtain T == T(x) *)
  tMax = HHalfPlane[gamma];
  xMax = 3;
  tMin = Quiet[
    1/2 * SeekRoot[xFun[#] - xMax &, {0, tMax}]
  , {Power::infy}];
  tStep = 0.001; (* Recall ySpacing = 0.005 above *)
  tInterp = Interpolation @ (
    Table[{xFun[t], t}, {t, tMax, tMin, -tStep}]
  );
  (* Plot *)
  LogPlot[tSol[x, 0] / tInterp[x] - 1 // Abs,
    {x, 0, xMax},
    AxesLabel -> {Italicise @ "x", "Rel. error"},
    AxesOrigin -> {-0.125, Automatic},
    Exclusions -> None,
    PlotRange -> Full,
    PlotOptions[Axes] // Evaluate
  ] // Ex @ FString["half_plane-solution-rel_error-gpd-{gpd}.pdf"]
]


(* ::Section:: *)
(*Figure: exact half-plane solution (half-plane-solution)*)


Module[
  {
    gamma,
    xMin, xMax,
    tStart, h, tEnd,
    plotRangePadding,
    angleMarkerRadius, labelSize, textStyle, textStyleZero,
    verticalTicksHorizontalOffset,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  gamma = 30 Degree;
  (* Plot range *)
  {xMin, xMax} = {0, 2};
  tStart = h = HHalfPlane[gamma];
  tEnd = SeekRoot[XHalfPlane[gamma][#] - xMax &, {0, tStart}, 10] // Quiet;
  (* Plot *)
  plotRangePadding = 0.15 h;
  angleMarkerRadius = 0.2 h;
  labelSize = LabelSize["Label"];
  textStyle = Style[#, labelSize] & @* LaTeXStyle;
  textStyleZero = Style[#, labelSize - 1] & @* LaTeXStyle;
  Show[
    EmptyAxes[{xMin, xMax}, {tStart, tEnd}
      , AxesLabel -> {
          Italicise["x"],
          Italicise["T"] // Margined @ {{0, 0}, {-4, 0}}
        }
      , ImageSize -> 0.48 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle[labelSize]
      , PlotRange -> All
      , PlotRangeClipping -> False
      , PlotRangePadding -> {{0, plotRangePadding}, {0, plotRangePadding}}
      , Ticks -> None
    ],
    (* Contact angle *)
    Graphics @ {GeneralStyle["DefaultThick"],
      Circle[
        {xMin, tStart},
        angleMarkerRadius,
        {-Pi/2, -Pi/2 + gamma}
      ]
    },
    Graphics @ {
      Text[
        "\[Gamma]" // textStyle
        , {xMin, tStart} + XYPolar[angleMarkerRadius, -Pi/2 + gamma/2]
        , {-0.55, 0.7}
      ]
    },
    (* Liquid *)
    ParametricPlot[
      {XHalfPlane[gamma][t], t}
      , {t, tStart, tEnd}
      , AspectRatio -> Automatic
      , ImageSize -> 480
      , PlotPoints -> 2
      , PlotRange -> {0, All}
      , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
    ],
    (* Manual ticks (for better spacing) *)
    Graphics @ {
      Text[
        "0" // textStyleZero
        , {0, 0}
        , {-0.25, 0.9}
      ],
      {}
    },
    verticalTicksHorizontalOffset = 2.5;
    Graphics @ {
      Text[
        "0" // textStyleZero
        , {0, 0}
        , {verticalTicksHorizontalOffset, -0.25}
      ],
      Text[
        Italicise["h"] // textStyle
        , {0, h}
        , {verticalTicksHorizontalOffset, -0.25}
      ],
      {}
    },
    {}
  ]
] // Ex["half_plane-solution.pdf"]


(* ::Section:: *)
(*Figure: half-plane mesh (half-plane-mesh-*)*)


Module[
  {
    mesh, meshWireframe,
    xMin, xMax, yMin, yMax,
    xSimplified, rectangleSimplified,
    posSimplified, meshWireframeSimplified,
    xHalf, havlesOptions, plotLeftHalf, plotRightHalf,
    xFineDetail, yFineDetail, plotFineDetail,
    dummyForTrailingCommas
  },
  (* Import mesh *)
  mesh = Import["half_plane-mesh.txt"] // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  {{xMin, xMax}, {yMin, yMax}} = mesh["Bounds"];
  (*
    Simplify mesh to reduce file size
    by replacing detailed portion with black rectangle
  *)
  xSimplified = 0.03;
  rectangleSimplified = Graphics @ Rectangle[{xMin, yMin}, {xSimplified, yMax}];
  posSimplified = MeshWireframePositions[mesh, {x_, y_} /; x > xSimplified];
  meshWireframeSimplified = mesh["Wireframe"[posSimplified]];
  (* Plot halves *)
  xHalf = Way[xMin, xMax];
  havlesOptions = {
    FrameLabel -> {
      Italicise["x"] // Margined @ {{0, 0}, {0, -0.06 ImageSizeTextWidth}},
      Italicise["y"]
    },
    Frame -> {{True, False}, {True, True}},
    FrameTicks -> {Automatic,
      Table[
        {y // N, y // NumberForm[#, {2, 1}] &}
        , {y, Subdivide[yMin, yMax, 2]}
      ]
    },
    FrameTicksStyle -> LabelSize["Tick"],
    ImagePadding -> {{Automatic, 0.02 ImageSizeTextWidth}, {Automatic, Automatic}},
    LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"],
    PlotRangeClipping -> True,
    PlotOptions[Frame],
    Nothing
  };
  plotLeftHalf =
    Show[
      meshWireframeSimplified,
      rectangleSimplified,
      {}
      , PlotRange -> {{xMin, xHalf}, {yMin, yMax}}
      , havlesOptions
    ];
  plotRightHalf =
    Show[
      meshWireframe
      , PlotRange -> {{xHalf, xMax}, {yMin, yMax}}
      , havlesOptions
    ];
  (* Plot fine detail *)
  xFineDetail = 0.2;
  yFineDetail = 0.1;
  plotFineDetail =
    Show[
      meshWireframe
      , Frame -> True
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -0.05 ImageSizeTextWidth}},
          Italicise["y"]
        }
      , ImagePadding -> Automatic
      , ImageSize -> 0.4 ImageSizeTextWidth
      , FrameTicks -> Automatic
      , PlotRange -> {{xMin, xFineDetail}, {-yFineDetail, yFineDetail}}
      , havlesOptions
    ];
  (* Export *)
  {
    GraphicsColumn[{plotLeftHalf, plotRightHalf}
      , ImageSize -> 0.58 ImageSizeTextWidth
      , ItemAspectRatio -> 0.4
      , Spacings -> -0.04 ImageSizeTextWidth
    ]
      // Ex["half_plane-mesh-halves.pdf"]
    ,
    Show[plotFineDetail
      , ImageSize -> 0.38 ImageSizeTextWidth
    ]
      // Ex["half_plane-mesh-detail.pdf"]
  }
]


(* ::Section:: *)
(*Figure: numerical solution relative error (half_plane-relative-error)*)


Module[
  {
    gpd, gamma,
    tComputedAssoc, tComputed,
    xFun, h,
    xSamplingValues,
    tExactValues, xExactValues,
    tComputedValues, relErrorValues,
    relErrorTable,
    xFirst, relErrorFirst,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  gpd = 1;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tComputedAssoc = Import["half_plane-solution.txt"] // Uncompress;
  tComputed = tComputedAssoc[gpd];
  (* Implicit exact half-plane solution x == x(T) *)
  xFun = XHalfPlane[gamma];
  h = HHalfPlane[gamma];
  (* Exact values *)
  xSamplingValues = Subdivide[0, 3, 30];
  tExactValues =
    Table[
      SeekRoot[xFun[#] - x &, {0, h}, 10] // Quiet
      , {x, xSamplingValues}
    ];
  xExactValues = xFun /@ tExactValues;
  (* Computed values and relative error *)
  tComputedValues = Table[tComputed[x, 0], {x, xExactValues}];
  relErrorValues = tComputedValues / tExactValues - 1 // Abs;
  (* Plot *)
  relErrorTable = {xExactValues, relErrorValues} // Transpose;
  ListLogPlot[
    relErrorTable
    , AxesLabel -> {
        Italicise["x"] // Margined @ {{0, 1}, {5, 0}},
        Style["Relative error", Smaller]
      }
    , ImageSize -> 0.6 ImageSizeTextWidth
    , Joined -> True
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotMarkers -> {Automatic, 5}
    , PlotStyle -> Black
    , PlotOptions[Axes] // Evaluate
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["half_plane-relative-error.pdf"]
