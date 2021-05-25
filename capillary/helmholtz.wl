(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "helmholtz"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Wedge half-angle*)


(* ::Subsubsection:: *)
(*\[Alpha]*)


alpha = Pi/4;


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*U*)


uKnown[x_, y_] := Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]];


(* ::Subsubsection:: *)
(*P = \[PartialD]U/\[PartialD]x, Q = \[PartialD]U/\[PartialD]y*)


p = Derivative[1, 0] @ uKnown;
q = Derivative[0, 1] @ uKnown;


(* ::Subsection:: *)
(*Flux function F*)


f = 1;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


vi[x_, y_] := p[x, y]^2 + q[x, y]^2 - f^2 // Expand // Evaluate


(* ::Subsection:: *)
(*Terminal curve x = x(y)*)


xTerm[y_] := 1/Sqrt[2] * Log[2 Cosh[Sqrt[2] y]];


(* ::Subsection:: *)
(*Traced boundaries (arc-length parametrisation)*)


xDerivative[x_, y_, branchSign_: 1] :=
  Divide[
    -q[x, y] f + branchSign p[x, y] Re @ Sqrt @ vi[x, y],
    p[x, y]^2 + q[x, y]^2
  ];
yDerivative[x_, y_, branchSign_: 1] :=
  Divide[
    +p[x, y] f + branchSign q[x, y] Re @ Sqrt @ vi[x, y],
    p[x, y]^2 + q[x, y]^2
  ];


(* Arc-length parametrisation x == x(s), y == y(s). *)
xyTraced[
  {xInitial_, yInitial_}, sInitial_,
  {sStart_, sEnd_},
  sSign_: 1,
  branchSign_: 1,
  terminationVi_: 0, terminateAtWalls_: False
] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    NDSolveValue[
      {
        x'[s] == sSign * xDerivative[x[s], y[s], branchSign],
        y'[s] == sSign * yDerivative[x[s], y[s], branchSign],
        x[sInitial] == xInitial,
        y[sInitial] == yInitial,
        WhenEvent[
          {
            vi[x[s], y[s]] < terminationVi,
            If[terminateAtWalls, Abs @ y[s] > x[s], False],
            False
          },
          "StopIntegration"
        ]
      }
      , {x, y}
      , {s, sStart, sEnd}
      , NoExtrapolation
    ]
  ];


(* ::Subsection:: *)
(*Contours (arc-length parametrisation)*)


(* Arc-length parametrisation x == x(s), y == y(s). *)
xyContour[
  {xInitial_, yInitial_}, sInitial_,
  {sStart_, sEnd_},
  sSign_: 1,
  terminateAtWalls_: False
] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    NDSolveValue[
      {
        x'[s] == sSign * +q[x[s], y[s]] / Sqrt[p[x[s], y[s]]^2 + q[x[s], y[s]]^2],
        y'[s] == sSign * -p[x[s], y[s]] / Sqrt[p[x[s], y[s]]^2 + q[x[s], y[s]]^2],
        x[sInitial] == xInitial,
        y[sInitial] == yInitial,
        WhenEvent[
          If[terminateAtWalls, Abs @ y[s] > x[s], False],
          "StopIntegration"
        ]
      }
      , {x, y}
      , {s, sStart, sEnd}
      , NoExtrapolation
    ]
  ];


(* ::Subsection:: *)
(*Critical terminal point x_0*)


x0 = Log[2] / Sqrt[2];
p @@ {x0, 0} == -1


(* ::Section:: *)
(*Algebra*)


(* ::Subsection:: *)
(*Check known solution*)


With[{x = \[FormalX], y = \[FormalY]},
  {
    (* Helmholtz equation *)
    Laplacian[uKnown[x, y], {x, y}] == uKnown[x, y],
    (* Boundary conditions *)
    XYPolar[1, alpha + Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> x Tan[alpha]},
    XYPolar[1, -alpha - Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> -x Tan[alpha]},
    (* Derivatives *)
    p[x, y] == -1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]]),
    q[x, y] ==  1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] - Exp[(-x-y)/Sqrt[2]]),
    Nothing
  } // FullSimplify
]


(* ::Subsection:: *)
(*Gradient squared*)


With[{x = \[FormalX], y = \[FormalY]},
  p[x, y]^2 + q[x, y]^2 // Expand
]


(* ::Subsection:: *)
(*Check viability*)


With[{x = \[FormalX], y = \[FormalY]},
  vi[x, y] == Exp[-Sqrt[2] (x-y)] + Exp[-Sqrt[2] (x+y)] - 1
    // FullSimplify
]


(* ::Subsection:: *)
(*Check terminal curve*)


With[{x = \[FormalX], y = \[FormalY]},
  vi[x, y] == 0 /. {x -> xTerm[y]}
    // FullSimplify
]


(* ::Section:: *)
(*Figure: wedge domain (helmholtz-wedge-domain)*)


(* RUN ON WINDOWS FOR ITALICISED OMEGA *)
Module[
  {
    wedgeXMax, wedgeYMax,
    clearanceProportion, xMin, xMax, yMax,
    alpha, angleMarkerRadius,
    normalVectorPosition, normalVectorLength,
    orthogonalityMarkerLength, orthogonalityMarkerStyle,
    normalThicknessCorrection,
    textStyle, testStyleOmega,
    fixSpacing,
    dummyForTrailingCommas
  },
  (* Wedge dimensions *)
  wedgeXMax = 1;
  wedgeYMax = 1;
  (* Plot range *)
  clearanceProportion = 1/10;
  xMin = Way[0, wedgeXMax, -clearanceProportion];
  xMax = Way[0, wedgeXMax, 1 + clearanceProportion];
  yMax = Way[0, wedgeYMax, 1 + clearanceProportion];
  (* Wedge half-angle *)
  alpha = Pi/4;
  angleMarkerRadius = 1/5 wedgeXMax;
  (* Normal vector *)
  normalVectorPosition = 1/2 {wedgeXMax, wedgeYMax};
  normalVectorLength = Sqrt[2]/3 wedgeXMax;
  orthogonalityMarkerLength = 1/14 wedgeXMax;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[None]];
  normalThicknessCorrection = -0.06;
  (* Diagram *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  testStyleOmega = Style[#, LabelSize["LabelOmega"]] & @* LaTeXStyle;
  fixSpacing[string_, number_, extra_: ""] :=
    Row @ {
      string,
      If[$OperatingSystem == "Windows", StringRepeat["\[NegativeThickSpace]", number] <> extra, ""]
    };
  (* Windows spacing fix *)
  Show[
    EmptyAxes[{xMin, xMax}, {-yMax, yMax}
      , AspectRatio -> Automatic
      , AxesStyle -> Darker[Gray]
      , ImageSize -> 0.4 * 0.8 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Method -> {"AxesInFront" -> False}
      , Ticks -> None
    ],
    (* Wedge half-angle *)
    Graphics @ {
      Circle[{0, 0}, angleMarkerRadius, {0, alpha}]
    },
    Graphics @ {
      Text[
        Equal[
          fixSpacing["\[Alpha]" // textStyle, 7],
          SeparatedRow[Style["\[NegativeVeryThinSpace]/\[NegativeVeryThinSpace]", Magnification -> 1]][
            fixSpacing["\[Pi]" // textStyle, 7],
            4
          ]
        ]
          // textStyle
        , XYPolar[angleMarkerRadius, alpha / 2]
        , {-1.2, -0.5}
      ]
    },
    (* Normal vector *)
    Graphics @ {Directive[GeneralStyle["Thick"], Arrowheads[Medium]],
      Arrow @ {
        normalVectorPosition,
        normalVectorPosition + XYPolar[normalVectorLength, alpha + Pi/2]
      }
    },
    Graphics @ {
      Text[
        Embolden["n"] // textStyle
        , normalVectorPosition + XYPolar[normalVectorLength, alpha + Pi/2]
        , {0, -1}
      ]
    },
    Graphics @ {orthogonalityMarkerStyle,
      Rectangle[
        normalVectorPosition,
        normalVectorPosition
          + orthogonalityMarkerLength {1 - normalThicknessCorrection, 1}
      ] // Rotate[#, alpha, normalVectorPosition] &
    },
    (* Wedge domain *)
    Graphics @ {
      Text[
        Italicise["\[CapitalOmega]"] // testStyleOmega
        , {3/4 wedgeXMax, 1/3 wedgeYMax}
        , If[$OperatingSystem == "Windows", {-0.7, -0.1}, {-1.2, -0.5}]
      ]
    },
    (* Wedge domain boundary *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{wedgeXMax, wedgeYMax}, {0, 0}, {wedgeXMax, -wedgeYMax}}
    },
    Graphics @ {
      Text[
        SeparatedRow[
          If[$OperatingSystem == "Windows",
            StringRepeat["\[NegativeThickSpace]", 9] <> "\[NegativeVeryThinSpace]" // textStyle,
            ""
          ]
        ] @@ {
          "\[PartialD]" // testStyleOmega,
          Italicise["\[CapitalOmega]"] // testStyleOmega
        }
        , 4/5 {wedgeXMax, wedgeYMax}
        , If[$OperatingSystem == "Windows", {0.5, -1}, {1, -1}]
      ]
    },
    {}
  ]
] // Ex["helmholtz-wedge-domain.pdf"]


(* ::Section:: *)
(*Figure: known solution (helmholtz-wedge-solution)*)


Module[
  {
    wallHorizontalThickness, wallHeight,
    xMinWall, xMax, yMax, yMaxWall,
    verticalEdge, wedgeBaseVertices,
    solutionBoundingBoxVertices,
    makeBaseVertex3D,
    dummyForTrailingCommas
  },
  (* Wall properties *)
  wallHorizontalThickness = 0.25;
  wallHeight = 1.21 uKnown[0, 0];
  (* Plot range *)
  xMinWall = -wallHorizontalThickness;
  xMax = 4;
  yMax = xMax/2;
  yMaxWall = yMax + wallHorizontalThickness;
  (* Manually drawn wedge wall wedges *)
  verticalEdge @ {x_, y_} := Line @ {{x, y, 0}, {x, y, wallHeight}};
  wedgeBaseVertices = {
    {yMax - wallHorizontalThickness / 2, -yMax - wallHorizontalThickness / 2},
    {yMax, -yMax},
    {yMax, +yMax},
    {yMax - wallHorizontalThickness / 2, +yMax + wallHorizontalThickness / 2},
    Nothing
  };
  solutionBoundingBoxVertices = {
    {yMax, -yMax},
    {xMax, 0},
    {yMax, +yMax},
    Nothing
  };
  makeBaseVertex3D @ {x_, y_} := {x, y, 0};
  (* Plot *)
  Show[
    (* Known solution *)
    Plot3D[
      uKnown[x, y], {x, xMinWall, xMax}, {y, -yMaxWall, yMaxWall}
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> Italicise @ {"x", "y", "U"}
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> Automatic
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Solution3D"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Label"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotPoints -> 20
      , PlotStyle -> BoundaryTracingStyle["Solution3D"]
      , RegionFunction -> Function[{x, y},
          And[
            Abs[y] < x,
            Abs[y] < xMax - x
          ]
        ]
      , TicksStyle -> LabelSize["Tick"]
      , ViewPoint -> {2.5, -1.1, 1.4}
    ],
    (* Wedge walls *)
    Plot3D[
      wallHeight, {x, xMinWall, xMax}, {y, -yMaxWall, yMaxWall}
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Wall3D"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , Mesh -> None
      , PlotPoints -> 20
      , PlotStyle -> BoundaryTracingStyle["Wall3D"]
      , RegionFunction -> Function[{x, y},
          And[
            Abs[y] - wallHorizontalThickness < x < Abs[y],
            Abs[y] < xMax - x
          ]
        ]
    ],
    (* Manually drawn wedge wall edges *)
    Graphics3D @ {
      BoundaryTracingStyle["Edge3D"],
      verticalEdge /@ wedgeBaseVertices,
      Line[makeBaseVertex3D /@ wedgeBaseVertices[[{1, 2}]]],
      Line[makeBaseVertex3D /@ wedgeBaseVertices[[{-1, -2}]]],
      Line @ {{0, 0, uKnown[0, 0]}, {0, 0, wallHeight}},
      {}
    },
    (* Manually drawn solution bounding box edges *)
    Graphics3D @ {
      BoundaryTracingStyle["Edge3D"],
      Line[makeBaseVertex3D /@ solutionBoundingBoxVertices],
      Line @ {{xMax, 0, 0}, {xMax, 0, uKnown[xMax, 0]}},
      {}
    },
    {}
    , ImageSize -> 0.55 ImageSizeTextWidth
  ]
] // Ex["helmholtz-wedge-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Figure: non-viable domain (helmholtz-viable)*)


Module[
  {
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 1;
  yMax = xMax;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , ImageSize -> 0.5 * 0.7 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{xMaxMore, yMaxMore}, {0, 0}, {xMaxMore, -yMaxMore}}
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0 && Abs[y] < x
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
  ]
] // Ex["helmholtz-viable.pdf"]


(* ::Section:: *)
(*Figure: generic traced boundaries (helmholtz-traced-boundaries)*)


Module[
  {
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    xyStartList, xyTracedList,
    sStartPlotting,
    sEndPlotting,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 1;
  yMax = xMax;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundaries (upper branch) *)
  xyStartList = Table[XYPolar[r, -alpha], {r, Subdivide[rMax, 10]}];
  xyTracedList =
    Table[
      xyTraced[xyStart, 0, {0, 2 rMax}, -1]
      , {xyStart, xyStartList}
    ];
  (* Arc-length to start and end plotting at *)
  sStartPlotting[xy_] := DomainStart[xy];
  sEndPlotting[xy_] :=
    Module[{sStart, sEnd},
      sStart = DomainStart[xy];
      sEnd = DomainEnd[xy];
      If[sEnd < xMaxMore,
        sEnd,
        SeekRoot[
          xy[[1]][#] - xMaxMore &,
          {sStart, sEnd},
          5
        ]
      ]
    ];
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , ImageSize -> 0.5 * 0.7 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{xMaxMore, yMaxMore}, {0, 0}, {xMaxMore, -yMaxMore}}
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0 && Abs[y] < x
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 3
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // IncludeYReflection
          // Evaluate
        , {s, sStartPlotting[xy], sEndPlotting[xy]}
        , PlotPoints -> 25
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedList}
    ],
    {}
  ]
] // Ex["helmholtz-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: legend (helmholtz-legend)*)


Module[
  {
    legendLabelStyle,
    row,
    dummyForTrailingCommas
  },
  legendLabelStyle = LatinModernLabelStyle[LabelSize["Legend"]];
  row = Join[
    CurveLegend[
      BoundaryTracingStyle /@ {"Wall", "Terminal"},
      {"wedge wall", "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ],
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ] // Margined @ {{0, 0}, {0, -3}} /@ # &,
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced"},
      {"traced boundary"}
      , LabelStyle -> legendLabelStyle
    ],
    {}
  ];
  GraphicsGrid[{row}
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.15
    , Spacings -> {{0, -0.1, 0.03, 0.06} ImageSizeTextWidth, 0}
  ]
] // Ex["helmholtz-legend.pdf"]


(* ::Section:: *)
(*Figure: traced boundaries, patched (helmholtz-traced-boundaries-patched)*)


Module[
  {
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    normalVectorLength,
    commonPlot,
    cornerList, numPlots,
    sRange,
    xyTracedUpperList, xyTracedLowerList,
    numCorners,
    sIntersectionLower, sIntersectionUpper,
    xTracedForNormal, yTracedForNormal, sNormal, xyNormal, normalVector,
    xyLower, xyUpper,
    textStyle,
    plotList,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 1;
  yMax = xMax;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  normalVectorLength = 1/3 xMax;
  (* Location of corners (increasing in y-coordinate) *)
  cornerList = Association[
    1 -> {{0.45, -0.3}},
    2 -> {{0.4, 0}, {0.45, 0.35}},
    3 -> {{0.3, -0.12}, {0.4, 0.25}, {0.5, 0.45}},
    4 -> {{0.45, -0.4}, {0.31, 0.1}, {0.5, 0.45}, {0.6, 0.6}},
    Nothing
  ];
  numPlots = Length[cornerList];
  (* Compute traced boundaries therethrough *)
  sRange = 2 rMax {-1, 1};
  Table[
    (* Upper branch *)
    xyTracedUpperList[n] = xyTraced[#, 0, sRange, -1, 1, 0, True] & /@ cornerList[n];
    (* Lower branch *)
    xyTracedLowerList[n] = xyTraced[#, 0, sRange, 1, -1, 0, True] & /@ cornerList[n];
    (* Intersections *)
    numCorners = Length @ cornerList[n];
    Table[
      (* Starting arc length for lower branch *)
      sIntersectionLower[n][k] =
        If[k == 1
          ,
          xyLower = xyTracedLowerList[n][[k]];
          SeekRoot[
            xyLower[[1]][#] - xMaxMore &,
            {DomainStart[xyLower], DomainEnd[xyLower]}, 5
          ]
          ,
          xyLower = xyTracedLowerList[n][[k]];
          xyUpper = xyTracedUpperList[n][[k - 1]];
          SeekParametricIntersection[xyLower, xyUpper] // First
        ];
      (* Ending arc length for upper branch *)
      sIntersectionUpper[n][k] =
        If[k == numCorners
          ,
          xyUpper = xyTracedUpperList[n][[k]];
          SeekRoot[
            xyUpper[[1]][#] - xMaxMore &,
            {DomainStart[xyUpper], DomainEnd[xyUpper]}, 5
          ]
          ,
          xyUpper = xyTracedUpperList[n][[k]];
          xyLower = xyTracedLowerList[n][[k + 1]];
          SeekParametricIntersection[xyUpper, xyLower] // First
        ];
      , {k, numCorners}
    ];
    (* Normal vector *)
    (* xyTracedForNormal, sNormal, xyNormal, normal *)
    {xTracedForNormal, yTracedForNormal} = xyTracedUpperList[n][[1]];
    sNormal = Way[0, sIntersectionUpper[n][1], If[n == 1, 1/3, 1/2]];
    xyNormal[n] = {xTracedForNormal, yTracedForNormal}[sNormal] // Through;
    normalVector[n] = (
      {xTracedForNormal', yTracedForNormal'}[sNormal]
        // Through
        // Normalize
        // Cross
    );
    , {n, numPlots}
  ];
  (* Common plot *)
  commonPlot = Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
    ],
    (* Wedge walls *)(*
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{xMaxMore, yMaxMore}, {0, 0}, {xMaxMore, -yMaxMore}}
    },*)
    {}
  ];
  (* Plots *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  plotList = Table[
    Show[
      (* Common plot *)
      commonPlot,
      (* Traced boundaries (upper) *)
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Background"]
        ]
        , {xy, xyTracedUpperList[n]}
      ],
      (* Traced boundaries (background) *)
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Background"]
        ]
        , {xy, xyTracedLowerList[n]}
      ],
      (* Normal vector *)
      Graphics @ {Directive[GeneralStyle["Thick"], Arrowheads[Medium]],
        Arrow @ {
          xyNormal[n],
          xyNormal[n] + normalVectorLength * normalVector[n]
        }
      },
      Graphics @ {
        Text[
          Embolden["n"] // textStyle
          , xyNormal[n] + normalVectorLength * normalVector[n]
          , {0.25, -0.9}
        ]
      },
      (* Traced boundaries (patched) *)
      numCorners = Length @ cornerList[n];
      Table[
        {
          (* Lower branch portions *)
          ParametricPlot[
            xyTracedLowerList[n][[k]][s]
              // Through
              // Evaluate
            , {s, sIntersectionLower[n][k], 0}
            , PlotPoints -> 2
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ],
          (* Upper branch portions *)
          ParametricPlot[
            xyTracedUpperList[n][[k]][s]
              // Through
              // Evaluate
            , {s, 0, sIntersectionUpper[n][k]}
            , PlotPoints -> 2
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ],
          {}
        }
        , {k, numCorners}
      ],
      {}
    ]
  , {n, numPlots}];
  GraphicsGrid[{plotList}
    , ImageSize -> ImageSizeTextWidth
  ]
] // Ex["helmholtz-traced-boundaries-patched.pdf"]


(* ::Section:: *)
(*Figure: terminal-points (helmholtz-terminal-points)*)


Module[
  {
    xMin, xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    xyStartList, xyContourList,
    xyContourOrdinary, sOrdinary, xyOrdinary,
    textStyle, textStyleBracket, textVerticalShift,
    plot,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMin = x0 / 2;
  xMax = 1;
  yMax = 1;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Contours *)
  xyStartList = Table[{x, 0}, {x, x0 + 0.06 Range[0, 3]}];
  xyContourList =
    Table[
      xyContour[xyStart, 0, 2 rMax {-1, 1}, 1, True]
      , {xyStart, xyStartList}
    ];
  (* Orindary terminal point to be shown *)
  xyContourOrdinary = xyContourList[[2]];
  sOrdinary = SeekRoot[
    xyContourOrdinary[[1]][#] - xTerm @ xyContourOrdinary[[2]][#] &,
    {0, xMaxMore}, 5
  ];
  xyOrdinary = xyContourOrdinary[sOrdinary] // Through;
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textVerticalShift = -0.23;
  (* Plot *)
  plot = Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , Frame -> None
      , PlotRangePadding -> {{Scaled[0.05], Scaled[0.04]}, None}
    ],
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0 && Abs[y] < x
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> None
      , PlotPoints -> 3
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Contours *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["ContourPlain"]
      ]
      , {xy, xyContourList}
    ],
    (* Terminal curve *)
    ContourPlot[
      vi[x, y] == 0
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , ContourStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {x0, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]0"],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {x0, 0},
        {1.5, textVerticalShift}
      ] // textStyle,
      Text[
        "critical",
        {x0, 0},
        {-1.35, textVerticalShift}
      ] // textStyle,
      {}
    },
    (* Ordinary terminal point *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ xyOrdinary
    },
    Graphics @ {
      Text[
        "ordinary",
        xyOrdinary,
        {-1.3, textVerticalShift}
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
      {Row @ {Italicise["U"], "\[Hyphen]contour"}, "terminal curve"}
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
        , Spacings -> {0, {-1.5, -1.5, -1.4}}
      ]
    }
    , ItemAspectRatio -> 1.75
    , ImageSize -> 0.63 ImageSizeTextWidth
    , Spacings -> {-0.1 ImageSizeTextWidth, 0}
  ]
] // Ex["helmholtz-terminal-points.pdf"]


(* ::Section:: *)
(*Figure: traced boundaries through (x_0, 0) (helmholtz-traced-boundaries-hyperbolic-*)*)


(* ::Subsection:: *)
(*Both*)


Module[
  {
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    xy,
    sStartPlotting, sEndPlotting,
    textStyle, textStyleBracket, textVerticalShift,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 1;
  yMax = xMax;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary (upper branch) *)
  xy = xyTraced[{x0, 0}, 0, 2 rMax {-1, 1}, -1, 1, 0, True];
  (* Arc-length to start and end plotting at *)
  sStartPlotting[xy_] := DomainStart[xy];
  sEndPlotting[xy_] :=
    Module[{sStart, sEnd},
      sStart = DomainStart[xy];
      sEnd = DomainEnd[xy];
      If[sEnd < xMaxMore,
        sEnd,
        SeekRoot[
          xy[[1]][#] - xMaxMore &,
          {sStart, sEnd},
          5
        ]
      ]
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textVerticalShift = -0.1;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
      , ImageSize -> 0.45 * 0.55 ImageSizeTextWidth
      , PlotRangePadding -> {{0.04, Automatic}, Automatic}
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{xMaxMore, yMaxMore}, {0, 0}, {xMaxMore, -yMaxMore}}
    },
    (* Traced boundaries *)
    ParametricPlot[
      xy[s]
        // Through
        // IncludeYReflection
        // Evaluate
      , {s, sStartPlotting[xy], sEndPlotting[xy]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {x0, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]0"],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {x0, 0},
        {1.45, textVerticalShift}
      ] // textStyle,
      {}
    },
    {}
  ]
] // Ex["helmholtz-traced-boundaries-hyperbolic-both.pdf"]


(* ::Subsection:: *)
(*Rounding*)


Module[
  {
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    xy,
    sStartPlotting, sEndPlotting,
    textStyle, textStyleBracket, textVerticalShift,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 1;
  yMax = xMax;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary (upper branch) *)
  xy = xyTraced[{x0, 0}, 0, 2 rMax {0, 1}, -1];
  (* Arc-length to start and end plotting at *)
  sStartPlotting[xy_] := DomainStart[xy];
  sEndPlotting[xy_] :=
    Module[{sStart, sEnd},
      sStart = DomainStart[xy];
      sEnd = DomainEnd[xy];
      If[sEnd < xMaxMore,
        sEnd,
        SeekRoot[
          xy[[1]][#] - xMaxMore &,
          {sStart, sEnd},
          5
        ]
      ]
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textVerticalShift = -0.1;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
      , ImageSize -> 0.45 * 0.55 ImageSizeTextWidth
      , PlotRangePadding -> {{0.04, Automatic}, Automatic}
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{xMaxMore, yMaxMore}, {0, 0}, {xMaxMore, -yMaxMore}}
    },
    (* Traced boundaries *)
    ParametricPlot[
      xy[s]
        // Through
        // IncludeYReflection
        // Evaluate
      , {s, sStartPlotting[xy], sEndPlotting[xy]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {x0, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]0"],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {x0, 0},
        {1.45, textVerticalShift}
      ] // textStyle,
      {}
    },
    {}
  ]
] // Ex["helmholtz-traced-boundaries-hyperbolic-rounding.pdf"]


(* ::Section:: *)
(*Figure: height-rise profiles (helmholtz-height-rise-profiles)*)


Module[
  {
    sMax,
    xySharp,
    xyRoundUpperHalf, xyRound,
    dummyForTrailingCommas
  },
  (* Plot range *)
  sMax = 3;
  (* Sharp corner (original walls) *)
  xySharp[s_] := XYPolar[Abs[s], Sign[s] alpha];
  (* Rounded corner (upper portion of upper branch) *)
  xyRoundUpperHalf = xyTraced[{x0, 0}, 0, {0, sMax}, -1];
  xyRound[s_] := xyRoundUpperHalf @ Abs[s] // Through;
  (* Plot *)
  Plot[
    {uKnown @@ xySharp[s], uKnown @@ xyRound[s]} // Evaluate
    , {s, -sMax, sMax}
    , AspectRatio -> Automatic
    , AxesLabel -> Italicise /@ {"s", "U"}
    , ImageSize -> 0.65 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {0, All}
    , PlotStyle -> {
        BoundaryTracingStyle["Wall"],
        Directive[BoundaryTracingStyle["Traced"], GeneralStyle["DefaultThick"]]
      }
    , PlotOptions[Axes] // Evaluate
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["helmholtz-height-rise-profiles.pdf"]
