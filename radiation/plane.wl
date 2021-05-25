(* ::Package:: *)

(* ::Text:: *)
(*See r1 (manuscripts/radiation-1-plane.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< Curvilinear`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "plane"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Traced boundary y = y(x)*)


(* ::Text:: *)
(*See (r1.17) to (r1.20) (Page r1-2 of manuscripts/radiation-1-plane.pdf).*)


(* NOTE: making yTraDer negative was a really dumb decision. Too late to change it though. *)
yTraDer = Function[{x}, -x^4 / Sqrt[1 - x^8]];
yTra = Function[{x}, Integrate[yTraDer[x], x] // Evaluate];
yTraMax = Abs @ yTra[1];


With[{x = \[FormalX]},
  {yTraDer[x], yTra[x], yTraMax}
]


(* ::Subsection:: *)
(*Self-incident boundary ratio*)


(* ::Subsubsection:: *)
(*Exact expression*)


(* ::Text:: *)
(*See (r6.15) (Page r6-2 of manuscripts/radiation-6-self.pdf).*)


boundaryRatioIntegrand[x_, xx_] :=
  Module[{y, yDer},
    y = -yTra[#] &;
    yDer = -yTraDer[#] &;
    Divide[
      Times[
        xx^4,
        -(x - xx) yDer[xx] + (y[x] - y[xx]),
        +(x - xx) yDer[x]  - (y[x] - y[xx])
      ],
      Times[
        2,
        ((x - xx)^2 + (y[x] - y[xx])^2) ^ (3/2),
        Sqrt[1 + yDer[x]^2]
      ]
    ]
  ] // Evaluate;


boundaryRatio[x1_, x2_][x_] :=
  1/x^4 * NIntegrate[boundaryRatioIntegrand[x, xx], {xx, x1, x2}];


(* ::Subsubsection:: *)
(*Crude upper bound*)


(* ::Text:: *)
(*See (r6.22) (Page r6-3 of manuscripts/radiation-6-self.pdf).*)


boundaryRatioBoundIntegrand[x_, xx_] :=
  Module[{y},
    y = -yTra[#] &;
    Divide[
      xx^4 * (x - xx)^4,
      ((x - xx)^2 + (y[x] - y[xx])^2) ^ (3/2)
    ]
  ] // Evaluate;


boundaryRatioBoundPrefactor[x1_, x2_][x_] :=
  Module[{y, yDer, yDer2},
    y = -yTra[#] &;
    yDer = -yTraDer[#] &;
    yDer2 = yDer';
    Divide[
      (yDer2 @ Max[x1, x2, x])^2,
      8 x^4 * Sqrt[1 + yDer[x]^2]
    ]
  ] // Evaluate;


boundaryRatioBound[x1_, x2_][x_] :=
  Times[
    boundaryRatioBoundPrefactor[x1, x2][x],
    NIntegrate[boundaryRatioBoundIntegrand[x, xx], {xx, x1, x2}]
  ]


(* ::Subsubsection:: *)
(*Ultra-crude upper bound*)


(* ::Text:: *)
(*See (r6.26) (Page r6-4 of manuscripts/radiation-6-self.pdf).*)


boundaryRatioBoundUltra[x1_, x2_] /; x1 < x2 :=
  Module[{y, yDer, yDer2},
    y = -yTra[#] &;
    yDer = -yTraDer[#] &;
    yDer2 = yDer';
    Divide[
      (x2 - x1)^2 * (x2/x1)^4 * yDer2[x2]^2,
      8 (1 + yDer[x1]^2) ^ 2
    ] // N
  ] // Evaluate;


(* ::Subsection:: *)
(*Parameters for figures*)


(* Avoid confusion (note yTra == -2F1(x; ...)) *)
yTraUpper[c_][x_] := c + yTra[x];
yTraLower[c_][x_] := c - yTra[x];


(* ::Subsubsection:: *)
(*Patched portions of traced boundaries*)


(*
  REMEMBER:
    the positive signed branch is the lower;
    the negative signed branch is the upper.
  This is because the sign in the tracing equation is \[MinusPlus].
  Therefore in a spike (shaped like <),
  the lower-branch curve is actually higher.
  However, writing each curve as y == c \[MinusPlus] 2F1(x; ...),
  note that the lower-branch curve has lower c.
  Now, to generate assorted patchings of traced boundaries,
  we choose the corners (x_1, y_1), ..., (x_n, y_n)
  with y_1 < ... < y_n
  and determine the intersections
  between lower-branch(i) and upper-branch(i+1)
  for i == 1, ..., n - 1.
  The corners (x_i, y_i) cannot be chosen arbitrarily,
  so don't be dumb.
 *)


(* Lists of corners *)
patchedCornerList = Association[
  "regular-long" -> Table[{0.4, y}, {y, Subdivide[-0.6, 0.6, 3]}],
  "regular-short" -> Table[{0.8, y}, {y, Subdivide[-0.7, 0.7, 8]}],
  "irregular" -> {
    {0.6, -0.54},
    {0.81, -0.435},
    {0.77, -0.15},
    {0.23, -0.06},
    {0.45, 0.},
    {0.7, 0.15},
    {0.3, 0.6}
  }
];


patchedIdList = Keys[patchedCornerList];


Table[
  {patchedCornerXList[id], patchedCornerYList} =
    Transpose @ patchedCornerList[id]
, {id, patchedIdList}];


(*
  Lists of constants {cUpper, cLower}
  for the pair of traced boundaries through each corner point
 *)
Table[
  Module[{xCorner, yCorner, cUpper, cLower},
    {patchedCUpperList[id], patchedCLowerList[id]} =
      Transpose @ Table[
        {xCorner, yCorner} = corner;
        cUpper = yCorner - yTra[xCorner];
        cLower = yCorner + yTra[xCorner];
        {cUpper, cLower}
      , {corner, patchedCornerList[id]}]
  ]
, {id, patchedIdList}];


(*
  Lists of intersections between lower-branch(i) and upper-branch(i+1)
 *)
Table[
  patchedCornerNum[id] = Length @ patchedCornerList[id];
  Module[{cUpper, cLower, xInt, yInt},
    {patchedIntXList[id], patchedIntYList[id]} =
      Transpose @ Table[
        cLower = patchedCLowerList[id][[i]];
        cUpper = patchedCUpperList[id][[i + 1]];
        xInt = SeekRoot[
          yTraUpper[cUpper][#] - yTraLower[cLower][#] &,
          {0, 1}
        ];
        yInt = yTraUpper[cUpper][xInt];
        {xInt, yInt}
      , {i, patchedCornerNum[id] - 1}]
  ]
, {id, patchedIdList}];


(* ::Subsubsection:: *)
(*Corners for constructing domains*)


(*
  Corners used to build domains (specified by index i)
  among each list of corners
*)
domainCornerRangeList["regular-long"] = {{2, 3}};
domainCornerRangeList["regular-short"] = {{3, 7}};
domainCornerRangeList["irregular"] = {{1, 2}, {3, 6}, {7, 7}};
(*
  Location of constant-temperature boundaries
  among each list of corners
*)
domainXBathList["regular-long"] = {1};
domainXBathList["regular-short"] = {Way[patchedIntXList["regular-short"][[3]], 1]};
domainXBathList["irregular"] = {
  Way[patchedIntXList["irregular"][[1]], patchedIntXList["irregular"][[2]], 1/3],
  Way[patchedIntXList["irregular"][[2]], patchedIntXList["irregular"][[-2]], 1/5],
  Way[patchedCornerXList["irregular"][[-2]], patchedIntXList["irregular"][[-2]]]
};


(* ::Subsection:: *)
(*Fin example*)


(* ::Subsubsection:: *)
(*Physical constants*)


stefanBoltzmann = 5.67 * 10^-8 (* W /m^2 /K^4 *);


(*
  ## Emissivity ##
  Figure 4, in:
  Wade & Preedy (2003).
  Fujihokka: A high-emissivity approach to aluminum anodizing.
  Metal Finishing, 101(12), 8--13.
  <https://doi.org/10.1016/S0026-0576(03)90093-6>
*)
exampleEmissivity = 0.9;


(*
  ## Conductivity ##
  Figure 1, in:
  Cook, Moore, Matsumura, & Van der Meer. (1975)
  The Thermal and Electrical Conductivity of Aluminum.
  ORNL-5079.
  Oak Ridge National Lab.
  <https://doi.org/10.2172/5066461>
*)
exampleConductivity = 236 (* W /m /K *);


(*
  ## Radiation constant ##
  c == \[Epsilon] \[Sigma] / k
*)
With[
  {
    eps = exampleEmissivity,
    sigma = stefanBoltzmann,
    k = exampleConductivity,
    dummyForTrailingCommas = Null
  },
  exampleConstant = eps sigma / k (* /m /K^3 *);
];


(* ::Subsubsection:: *)
(*Interval*)


{exampleX1, exampleX2} = {0.5, 0.6};


(* ::Subsubsection:: *)
(*Profile*)


exampleYTraced[x_] := yTraLower[0][x] - yTraLower[0][exampleX1] // Evaluate;
exampleY2 = exampleYTraced[exampleX2];


(* ::Subsubsection:: *)
(*Tip temperature T_1*)


With[
  {
    c = exampleConstant,
    x1Hat = exampleX1,
    x2Hat = exampleX2,
    dummyForTrailingCommas = Null
  },
  exampleT1[length_] :=
    x1Hat
    Divide[
      x2Hat - x1Hat,
      c * length
    ]^(1/3);
];


(* ::Subsubsection:: *)
(*Base temperature T_2*)


With[
  {
    c = exampleConstant,
    x1Hat = exampleX1,
    x2Hat = exampleX2,
    dummyForTrailingCommas = Null
  },
  exampleT2[length_] :=
    x2Hat
    Divide[
      x2Hat - x1Hat,
      c * length
    ]^(1/3);
];


(* ::Subsubsection:: *)
(*Power per length p*)


With[
  {
    k = exampleConductivity,
    c = exampleConstant,
    x1Hat = exampleX1,
    x2Hat = exampleX2,
    y2Hat = exampleY2,
    dummyForTrailingCommas = Null
  },
  exampleP[length_] :=
    2 y2Hat k
    Divide[
      x2Hat - x1Hat,
      c * length
    ]^(1/3);
];


(* ::Subsubsection:: *)
(*Celsius, ugh*)


celsiusOffset = 273.15 (* K *);


(* ::Subsection:: *)
(*Directional dependence of output*)


directionalDependenceIntegrand[xx_][phi_] :=
  Module[{yDer},
    yDer = -yTraDer[#] &;
    xx^4 * (yDer[xx] Cos[phi] + Sin[phi])
  ] // Evaluate;


directionalDependence[x1_, x2_][phi_?NumericQ] :=
  NIntegrate[directionalDependenceIntegrand[xx][phi], {xx, x1, x2}];


(* ::Section:: *)
(*Traced boundary algebra*)


(* ::Subsection:: *)
(*General*)


With[{x = \[FormalX]},
  {
    {"dy/dx", yTraDer[x]},
    {"y", yTra[x]},
    {"y(1)", yTraMax},
    {"y(1)", yTraMax // N}
  } // TableForm
] // Ex["plane-traced-general.pdf"]


(* ::Subsection:: *)
(*Near x = 0*)


With[{x = \[FormalX]},
  Module[{yDerAsy, yAsy},
    yDerAsy = yTraDer[x] + O[x]^12;
    yAsy = Integrate[yDerAsy, x];
    {
      {"dy/dx", yDerAsy},
      {"y", yAsy}
    } // TableForm
  ]
] // Ex["plane-traced-x-near-0.pdf"]


(* ::Subsection:: *)
(*Near x = 1*)


With[{xi = \[FormalXi]},
  Module[{yDerAsy, yAsy},
    yDerAsy = -yTraDer[1 - xi] + O[xi]^(1/2);
    yAsy = Integrate[yDerAsy, xi];
    {
      {"dy/d\[Xi]", yDerAsy},
      {"y", yAsy}
    } // TableForm
  ]
] // Ex["plane-traced-x-near-1.pdf"]


(* ::Section:: *)
(*Traced boundary plot*)


(* ::Subsection:: *)
(*Without known solution*)


Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
  AspectRatio -> Automatic,
  PlotRange -> Full,
  PlotStyle -> Black,
  PlotOptions[Axes] // Evaluate
] // Ex["plane-traced.pdf"]


(* ::Subsection:: *)
(*With known solution*)


Show[
  EmptyAxes[{0, 1}, {-yTraMax, yTraMax}],
  (* Known solution *)
  DensityPlot[x, {x, 0, 1}, {y, -yTraMax, yTraMax},
    ColorFunction -> "TemperatureMap",
    RegionFunction -> Function[{x, y}, Abs[y] < -yTra[x]]
  ],
  (* Traced boundaries *)
  Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
    AspectRatio -> Automatic,
    PlotRange -> Full,
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["plane-traced-with-known.png"]


(* ::Section:: *)
(*Self-incident boundary ratio plots*)


(* ::Subsection:: *)
(*Exact expression*)


Module[
  {
    intervalList,
    numPoints,
    xSub, flatFractions, stringFractions,
    x1, x2,
    x1String, x2String, fileName,
    xValues, xRatioTable,
    dummyForTrailingCommas
  },
  (* List of intervals {x1, x2} to plot the ratio for *)
  intervalList = {
    {0, 1},
    {0, 1/2}, {1/2, 1},
    {0, 1/3}, {1/3, 2/3}, {2/3, 1},
    {0, 1/4}, {1/4, 1/2}, {1/2, 3/4}, {3/4, 1},
    Nothing
  };
  (* Number of points for sampling *)
  numPoints = 50;
  (* Nice labels *)
  xSub[n_] := Subscript[Italicise["x"], n];
  flatFractions = Rational[a_, b_] :> SeparatedRow["/"][a, b];
  stringFractions = {
    Rational[a_, b_] :> StringForm["``o``", a, b],
    n_Integer :> ToString[n]
  };
  (* Make plots *)
  Table[
    (* Get interval endpoints *)
    {x1, x2} = interval;
    (* Build table of values *)
    xValues = Subdivide[x1, x2, numPoints];
    xRatioTable = Table[{x, boundaryRatio[x1, x2][x]}, {x, xValues}];
    (* Plot and export *)
    {x1String, x2String} = {x1, x2} /. stringFractions;
    fileName = StringTemplate["plane-boundary-ratio-x1-``-x2-``.png"][x1String, x2String];
    ListPlot[xRatioTable
      , Joined -> True
      , AxesLabel -> Italicise /@ {"x", "R"}
      , PlotLabel -> {xSub[1], xSub[2]} == (interval /. flatFractions)
      , PlotOptions[Axes] // Evaluate
    ] // Ex[fileName]
    , {interval, intervalList}
  ] // Quiet
]


(* ::Subsection:: *)
(*Crude upper bound*)


Module[
  {
    intervalList,
    numPoints,
    xSub, flatFractions, stringFractions,
    x1, x2,
    x1String, x2String, fileName,
    xValues, xRatioTable,
    dummyForTrailingCommas
  },
  (* List of intervals {x1, x2} to plot the ratio for *)
  intervalList = {
    {0, 1},
    {0, 1/2}, {1/2, 1},
    {0, 1/3}, {1/3, 2/3}, {2/3, 1},
    {0, 1/4}, {1/4, 1/2}, {1/2, 3/4}, {3/4, 1},
    Nothing
  };
  (* Number of points for sampling *)
  numPoints = 50;
  (* Nice labels *)
  xSub[n_] := Subscript[Italicise["x"], n];
  flatFractions = Rational[a_, b_] :> SeparatedRow["/"][a, b];
  stringFractions = {
    Rational[a_, b_] :> StringForm["``o``", a, b],
    n_Integer :> ToString[n]
  };
  (* Make plots *)
  Table[
    (* Get interval endpoints *)
    {x1, x2} = interval;
    (* Build table of values *)
    xValues = Subdivide[x1, x2, numPoints];
    xRatioTable = Table[{x, boundaryRatio[x1, x2][x]}, {x, xValues}];
    (* Plot and export *)
    {x1String, x2String} = {x1, x2} /. stringFractions;
    fileName = StringTemplate["plane-boundary-ratio-with-bound-x1-``-x2-``.png"][x1String, x2String];
    Show[
      (* Crude upper bound *)
      Plot[boundaryRatioBound[x1, x2][x]
        , {x, x1, x2}
        , AxesLabel -> Italicise /@ {"x", "R"}
        , PlotLabel -> Column @ {
            {xSub[1], xSub[2]} == (interval /. flatFractions),
            Subscript[Italicise["R"], "ultra"] == boundaryRatioBoundUltra[x1, x2]
          }
        , PlotRange -> {{x1, x2}, {0, Automatic}}
        , PlotStyle -> Directive[Red, Dashed]
        , PlotOptions[Axes] // Evaluate
      ],
      (* Exact expression *)
      ListPlot[xRatioTable
        , Joined -> True
        , PlotRange -> All
      ],
      {}
    ] // Ex[fileName]
    , {interval, intervalList}
  ] // Quiet
]


(* ::Section:: *)
(*Fin example (length = 1 metre)*)


Module[
  {
    k, c, x1Hat, x2Hat, y2Hat,
    l, h, temp1, temp2, p,
    dummyForTrailingCommas
  },
  (* Abbreviations *)
  k = exampleConductivity;
  c = exampleConstant;
  {x1Hat, x2Hat} = {exampleX1, exampleX2};
  y2Hat = exampleY2;
  (* Computed quantities *)
  l = 1;
  h = 2 y2Hat l / (x2Hat - x1Hat);
  temp1 = x1Hat ((x2Hat - x1Hat) / (c l))^(1/3);
  temp2 = x2Hat ((x2Hat - x1Hat) / (c l))^(1/3);
  p = 2 y2Hat k ((x2Hat - x1Hat) / (c l))^(1/3);
  (* Make table *)
  {
    {"Fin length", "L",
      l // Quantity[#, "Meters"] &
    },
    {"Fin thickness at base", "H",
      h // Quantity[#, "Meters"] &
    },
    {"Temperature at tip", "T_1",
      temp1 // Quantity[#, "Kelvins"] &,
      temp1 - celsiusOffset // Quantity[#, "DegreesCelsius"] &,
      {}
    },
    {"Temperature at base", "T_2",
      temp2 // Quantity[#, "Kelvins"] &,
      temp2 - celsiusOffset // Quantity[#, "DegreesCelsius"] &,
      temp2 == exampleT2[l] (* Check *)
    },
    {"Power per length", "p",
      p // Quantity[#, "Watts" / "Meters"] &,
      {},
      p == exampleP[l] (* Check *)
    },
    {"Max self-viewing ratio", "R",
      Table[
        boundaryRatio[x1Hat, x2Hat][xHat]
        , {xHat, Subdivide[x1Hat, x2Hat, 10]}
      ] // Max // ScientificForm
    },
    Nothing
  } // TableForm
] // Ex["plane-fin-1-metre.pdf"]


(* ::Section:: *)
(*Figure: Slab BVP (plane-slab-bvp.pdf)*)


Module[
  {
    xSlabMin, xSlabMax,
    ySlabMin, ySlabMax, ySlabMid,
    xRadiation, yRadiationMargin,
    numRadiation, yRadiationValues,
    radiationArrowStyle, radiationArrowSize,
    xAxisMin, xAxisMax, yAxis,
    axisStyle,
    textStyle,
    dummyForTrailingCommas
  },
  (* Slab dimensions *)
  {xSlabMin, xSlabMax} = {0, 1};
  {ySlabMin, ySlabMax} = {0, 1.7};
  ySlabMid = Way[ySlabMin, ySlabMax];
  (* Radiation arrow specs *)
  xRadiation = Way[xSlabMin, xSlabMax, -0.15];
  yRadiationMargin = Way[ySlabMin, ySlabMax, 0.05];
  numRadiation = 6;
  yRadiationValues = Subdivide[
    ySlabMin + yRadiationMargin,
    ySlabMax - yRadiationMargin,
    numRadiation - 1
  ];
  radiationArrowStyle = Arrowheads[0.025];
  radiationArrowSize = 2.2 Abs[xRadiation];
  (* x-axis specs *)
  xAxisMin = xRadiation - 2 radiationArrowSize;
  xAxisMax = xSlabMax + 2 radiationArrowSize;
  yAxis = ySlabMin - 0.7 radiationArrowSize;
  axisStyle = Directive[Arrowheads[0.05]];
  (* Make diagram *)
  textStyle = Style[#, LabelSize["Label"] - 1] & @* LaTeXStyle;
  Show[
    (* Slab horizontal boundaries *)
    Graphics @ {GeneralStyle["DefaultThick"],
      Line @ {{xSlabMin, ySlabMax}, {xSlabMax, ySlabMax}},
      Line @ {{xSlabMin, ySlabMin}, {xSlabMax, ySlabMin}},
      {}
    },
    (* Slab radiation boundary *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Line @ {{xSlabMin, ySlabMax}, {xSlabMin, ySlabMin}},
      {}
    },
    (* Radiation arrows *)
    Graphics @ {radiationArrowStyle,
      Table[
        SquigglyArrow[{xRadiation, yRadiation}, Pi, radiationArrowSize]
        , {yRadiation, yRadiationValues}
      ],
      Text[
        "radiation" // textStyle
        , {xRadiation - radiationArrowSize, ySlabMid}
        , {1.3, -0.2}
      ],
      {}
    },
    (* Constant temperature boundary *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {{xSlabMax, ySlabMax}, {xSlabMax, ySlabMin}},
      Text[
        Column[{"constant", "temperature"}
          , Alignment -> Center
          , Spacings -> 0
        ] // textStyle
        , {xSlabMax, ySlabMid}
        , {-1.2, 0}
      ],
      {}
    },
    (* x-axis *)
    Graphics @ {axisStyle,
      Arrow @ {{xAxisMin, yAxis}, {xAxisMax, yAxis}},
      Text[
        Italicise["x"] // textStyle
        , {xAxisMax, yAxis}
        , {-2.5, -0.2}
      ],
      {}
    },
    {}
    , ImageSize -> 0.45 ImageSizeTextWidth
  ]
] // Ex["plane-slab-bvp.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    xSlabMin, xSlabMax,
    ySlabMin, ySlabMax, ySlabMid,
    xRadiation, yRadiationMargin,
    numRadiation, yRadiationValues,
    radiationArrowStyle, radiationArrowSize,
    xAxisMin, xAxisMax, yAxis,
    axisStyle,
    textStyle,
    dummyForTrailingCommas
  },
  (* Slab dimensions *)
  {xSlabMin, xSlabMax} = {0, 1};
  {ySlabMin, ySlabMax} = {0, 1.7};
  ySlabMid = Way[ySlabMin, ySlabMax];
  (* Radiation arrow specs *)
  xRadiation = Way[xSlabMin, xSlabMax, -0.15];
  yRadiationMargin = Way[ySlabMin, ySlabMax, 0.05];
  numRadiation = 6;
  yRadiationValues = Subdivide[
    ySlabMin + yRadiationMargin,
    ySlabMax - yRadiationMargin,
    numRadiation - 1
  ];
  radiationArrowStyle = Arrowheads[0.025];
  radiationArrowSize = 2.2 Abs[xRadiation];
  (* x-axis specs *)
  xAxisMin = xRadiation - 2 radiationArrowSize;
  xAxisMax = xSlabMax + 2 radiationArrowSize;
  yAxis = ySlabMin - 0.7 radiationArrowSize;
  axisStyle = Directive[Arrowheads[0.05]];
  (* Make diagram *)
  textStyle = Style[#, 9] &;
  Show[
    (* Slab *)
    Graphics @ {
      SlidesStyle["InteriorRegion"],
      Rectangle[{xSlabMin, ySlabMin}, {xSlabMax, ySlabMax}]
    },
    (* Slab horizontal boundaries *)
    Graphics @ {GeneralStyle["DefaultThick"],
      Line @ {{xSlabMin, ySlabMax}, {xSlabMax, ySlabMax}},
      Line @ {{xSlabMin, ySlabMin}, {xSlabMax, ySlabMin}},
      {}
    },
    (* Slab radiation boundary *)
    Graphics @ {
      BoundaryTracingStyle["Traced"],
      SlidesStyle["Boundary"],
      Line @ {{xSlabMin, ySlabMax}, {xSlabMin, ySlabMin}},
      {}
    },
    (* Radiation arrows *)
    Graphics @ {radiationArrowStyle,
      Table[
        SquigglyArrow[{xRadiation, yRadiation}, Pi, radiationArrowSize]
        , {yRadiation, yRadiationValues}
      ],
      Text[
        "radiation" // textStyle
        , {xRadiation - radiationArrowSize, ySlabMid}
        , {1.3, -0.2}
      ],
      {}
    },
    (* Constant temperature boundary *)
    Graphics @ {
      {
        BoundaryTracingStyle["Contour"], SlidesStyle["Source"],
        Line @ {{xSlabMax, ySlabMax}, {xSlabMax, ySlabMin}}
      },
      Text[
        Column[{"constant", "temperature"}
          , Alignment -> Center
          , Spacings -> 0
        ] // textStyle
        , {xSlabMax, ySlabMid}
        , {-1.2, 0}
      ],
      {}
    },
    (* x-axis *)
    Graphics @ {axisStyle,
      Arrow @ {{xAxisMin, yAxis}, {xAxisMax, yAxis}},
      Text[
        Italicise["x"] // textStyle
        , {xAxisMax, yAxis}
        , {-2.5, -0.2}
      ],
      {}
    },
    {}
    , ImageSize -> 0.55 ImageSizeTextWidthBeamer
  ]
] // Ex["plane-slab-bvp-slides.pdf"];


(* ::Section:: *)
(*Figure: Traced boundaries, single spike (plane-traced-boundary-spike.pdf)*)


Module[{yMax, yMaxMore, xTerm},
  yMax = Ceiling[yTraMax, 0.1];
  yMaxMore = 1.2 yTraMax;
  xTerm = 1;
  Show[
    (* Pair of traced boundaries forming a spike *)
    Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
      AspectRatio -> Automatic,
      ImageSize -> 240,
      LabelStyle -> LatinModernLabelStyle[12],
      PlotRange -> {-yMax, yMax},
      PlotStyle -> BoundaryTracingStyle["Traced"],
      PlotOptions[Axes] // Evaluate
    ],
    (* Critical terminal curve *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {{xTerm, -yMaxMore}, {xTerm, yMaxMore}}
    }
  ]
] // Ex["plane-traced-boundary-spike.pdf"]


(* ::Section:: *)
(*Figure: Traced boundaries, various (plane-traced-boundary-various.pdf)*)


(*
  REMEMBER:
    the positive signed branch is the lower;
    the negative signed branch is the upper.
  This is because the sign in the tracing equation is \[MinusPlus].
  Therefore in a spike (shaped like <),
  the lower-branch curve is actually higher.
  However, writing each curve as y == c \[MinusPlus] 2F1(x; ...),
  note that the lower-branch curve has lower c.
  Now, to generate assorted patchings of traced boundaries,
  we choose the corners (x_1, y_1), ..., (x_n, y_n)
  with y_1 < ... < y_n
  and determine the intersections
  between lower-branch(i) and upper-branch(i+1)
  for i == 1, ..., n - 1.
  The corners (x_i, y_i) cannot be chosen arbitrarily,
  so don't be dumb.
 *)
Module[
 {yTraUpper, yTraLower,
  cornerListList,
  xMin, xMax, yMin, yMax,
  imageSize,
  plotList,
  xCornerList, yCornerList,
  cUpperList, cLowerList,
  xCorner, yCorner,
  cUpper, cLower,
  n, xIntList, yIntList,
  restrict,
  xInt, yInt,
  xLeft, xRight,
  xTerm
 },
  (* Avoid confusion (note yTra == -2F1(x; ...)) *)
  yTraUpper[c_][x_] := c + yTra[x];
  yTraLower[c_][x_] := c - yTra[x];
  (* List of lists of corners *)
  cornerListList = List[
    Table[{0.4, y}, {y, Subdivide[-0.6, 0.6, 4]}],
    Table[{0.8, y}, {y, Subdivide[-0.7, 0.7, 8]}],
    {0.1, 0.3} # & /@ {
      {4, -2.2},
      {9, -1.3},
      {7.7, -0.5},
      {2.3, -0.2},
      {4.5, 0},
      {7, 0.5},
      {3, 2}
    }
  ];
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.05 xTerm;
  yMax = 0.7;
  imageSize = 240;
  (* Build a plot for each of these lists *)
  plotList = Table[
    {xCornerList, yCornerList} = Transpose[cornerList];
    (*
      List of constants {cUpper, cLower}
      for the pair of traced boundaries
      through each corner point
     *)
    {cUpperList, cLowerList} =
      Transpose @ Table[
        {xCorner, yCorner} = corner;
        cUpper = yCorner - yTra[xCorner];
        cLower = yCorner + yTra[xCorner];
        {cUpper, cLower}
      , {corner, cornerList}];
    (*
      List of intersections between
      lower-branch(i) and upper-branch(i+1)
     *)
    n = Length[cornerList];
    {xIntList, yIntList} =
      Transpose @ Table[
        cLower = cLowerList[[i]];
        cUpper = cUpperList[[i + 1]];
        xInt = SeekRoot[
          yTraUpper[cUpper][#] - yTraLower[cLower][#] &,
          {0, 1}
        ];
        yInt = yTraUpper[cUpper][xInt];
        {xInt, yInt}
      , {i, n - 1}];
    (* Restrict domain *)
    restrict[x0_, x1_] :=
      Piecewise[
        {{1, x0 <= # <= x1}},
        Indeterminate
      ] &;
    (* Plot *)
    Show[
      (* Traced boundaries *)
      Plot[
        Table[
          {
            (* Upper-branch(i) *)
            xLeft = xCornerList[[i]];
            xRight = If[i > 1,
              xIntList[[i - 1]],
              Max[xIntList]
            ];
            cUpper = cUpperList[[i]];
            (
              yTraUpper[cUpper][x] restrict[xLeft, xRight][x]
            ),
            (* Lower-branch(i) *)
            xLeft = xCornerList[[i]];
            xRight = If[i < n,
              xIntList[[i]],
              Max[xIntList]
            ];
            cLower = cLowerList[[i]];
            (
              yTraLower[cLower][x] restrict[xLeft, xRight][x]
            )
          }
        , {i, n}] // Evaluate,
        {x, 0, 1},
        AspectRatio -> Automatic,
        Axes -> None,
        ImageSize -> 240,
        PlotRange -> {{xMin, xMax}, {-yMax, yMax}},
        PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
      (* Critical terminal curve *)
      Graphics @ {BoundaryTracingStyle["Contour"],
        Line @ {{xTerm, -yMax}, {xTerm, yMax}}
      }
    ]
  , {cornerList, cornerListList}]
  // GraphicsRow[#,
    Spacings -> {
      {0, -0.2, 0.55} imageSize,
      Automatic
    }
  ] &
] // Ex["plane-traced-boundary-various.pdf"]


(* ::Section:: *)
(*Figure: Traced boundaries (plane-traced-boundaries.pdf)*)


Module[
 {xMin, xMax, yMax,
  xTerm,
  cMax, cList
 },
  (* Plot range *)
  xMin = 0;
  xMax = 1;
  yMax = Ceiling[2 yTraMax, 0.1];
  (* Critical terminal curve *)
  xTerm = 1;
  (* Values of integration constant *)
  cMax = 0.95 yMax;
  cList = Subdivide[-yMax, yMax, 8];
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      AspectRatio -> Automatic,
      FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
        Italicise["y"]
      },
      FrameTicksStyle -> LabelSize["Tick"],
      ImageSize -> 0.45 ImageSizeTextWidth,
      LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Traced boundaries *)
    Table[
      Plot[c + yTra[x] {1, -1}, {x, xMin, xMax},
        PlotPoints -> 3,
        PlotRange -> {-yMax, yMax},
        PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {c, cList}],
    (* Critical terminal curve *)
    ParametricPlot[
      {xTerm, y}, {y, -yMax, yMax},
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["Traced"]
    ]
  ]
] // Ex["plane-traced-boundaries.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


(* ::Subsubsection:: *)
(*Everything*)


Module[
 {xMin, xMax, yMax,
  xMinPlot, xMaxPlot, yMaxPlot,
  xTerm,
  cMax, cList
 },
  (* Viable & physical range *)
  xMin = 0;
  xMax = 1;
  yMax = Ceiling[2 yTraMax, 0.2];
  (* Actual plot range *)
  xMinPlot = -0.267;
  xMaxPlot = +1.27;
  yMaxPlot = 0.95 yMax;
  (* Critical terminal curve *)
  xTerm = 1;
  (* Values of integration constant *)
  cMax = 0.95 yMax;
  cList = Subdivide[-yMax, yMax, 8];
  Show[
    EmptyFrame[{xMinPlot, xMaxPlot}, {-yMaxPlot, yMaxPlot},
      AspectRatio -> Automatic,
      FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -10}},
        Italicise["y"]
      },
      FrameTicksStyle -> 8,
      ImageSize -> 0.85 * 0.5 ImageSizeTextWidthBeamer,
      LabelStyle -> 10
    ],
    (* Unphysical region *)
    Graphics @ {BoundaryTracingStyle["Unphysical"],
      Rectangle[
        {Way[xMinPlot, xMin, -1], -2 yMaxPlot},
        {xMin, +2 yMaxPlot}
      ]
    },
    (* Non-viable domain *)
    Graphics @ {BoundaryTracingStyle["NonViable"],
      Rectangle[
        {xMax, -2 yMaxPlot},
        {Way[xMax, xMaxPlot, 2], +2 yMaxPlot}
      ]
    },
    (* Traced boundaries *)
    Table[
      Plot[c + yTra[x] {1, -1}, {x, xMin, xMax},
        PlotPoints -> 3,
        PlotRange -> {-yMax, yMax},
        PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
      ]
    , {c, cList}],
    (* Critical terminal curve *)
    ParametricPlot[
      {xTerm, y}, {y, -yMax, yMax},
      PlotPoints -> 2,
      PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
    ]
  ]
] // Ex["plane-traced-boundaries-slides_everything.pdf"];


(* ::Subsubsection:: *)
(*Both regions (unphysical and non-viable)*)


Module[
 {xMin, xMax, yMax,
  xMinPlot, xMaxPlot, yMaxPlot,
  xTerm,
  cMax, cList
 },
  (* Viable & physical range *)
  xMin = 0;
  xMax = 1;
  yMax = Ceiling[2 yTraMax, 0.2];
  (* Actual plot range *)
  xMinPlot = -0.267;
  xMaxPlot = +1.27;
  yMaxPlot = 0.95 yMax;
  (* Critical terminal curve *)
  xTerm = 1;
  (* Values of integration constant *)
  cMax = 0.95 yMax;
  cList = Subdivide[-yMax, yMax, 8];
  Show[
    EmptyFrame[{xMinPlot, xMaxPlot}, {-yMaxPlot, yMaxPlot},
      AspectRatio -> Automatic,
      FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -10}},
        Italicise["y"]
      },
      FrameTicksStyle -> 8,
      ImageSize -> 0.85 * 0.5 ImageSizeTextWidthBeamer,
      LabelStyle -> 10
    ],
    (* Unphysical region *)
    Graphics @ {BoundaryTracingStyle["Unphysical"],
      Rectangle[
        {Way[xMinPlot, xMin, -1], -2 yMaxPlot},
        {xMin, +2 yMaxPlot}
      ]
    },
    (* Non-viable domain *)
    Graphics @ {BoundaryTracingStyle["NonViable"],
      Rectangle[
        {xMax, -2 yMaxPlot},
        {Way[xMax, xMaxPlot, 2], +2 yMaxPlot}
      ]
    },
    (* Critical terminal curve *)
    ParametricPlot[
      {xTerm, y}, {y, -yMax, yMax},
      PlotPoints -> 2,
      PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
    ]
  ]
] // Ex["plane-traced-boundaries-slides_regions.pdf"];


(* ::Subsubsection:: *)
(*Unphysical region*)


Module[
 {xMin, xMax, yMax,
  xMinPlot, xMaxPlot, yMaxPlot,
  xTerm,
  cMax, cList
 },
  (* Viable & physical range *)
  xMin = 0;
  xMax = 1;
  yMax = Ceiling[2 yTraMax, 0.2];
  (* Actual plot range *)
  xMinPlot = -0.267;
  xMaxPlot = +1.27;
  yMaxPlot = 0.95 yMax;
  (* Critical terminal curve *)
  xTerm = 1;
  (* Values of integration constant *)
  cMax = 0.95 yMax;
  cList = Subdivide[-yMax, yMax, 8];
  Show[
    EmptyFrame[{xMinPlot, xMaxPlot}, {-yMaxPlot, yMaxPlot},
      AspectRatio -> Automatic,
      FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -10}},
        Italicise["y"]
      },
      FrameTicksStyle -> 8,
      ImageSize -> 0.85 * 0.5 ImageSizeTextWidthBeamer,
      LabelStyle -> 10
    ],
    (* Unphysical region *)
    Graphics @ {BoundaryTracingStyle["Unphysical"],
      Rectangle[
        {Way[xMinPlot, xMin, -1], -2 yMaxPlot},
        {xMin, +2 yMaxPlot}
      ]
    },
    {}
  ]
] // Ex["plane-traced-boundaries-slides_unphysical.pdf"];


(* ::Section:: *)
(*Figure: Traced boundaries, patched (plane-traced-boundaries-patched.pdf)*)


Module[
 {xTerm,
  xMin, xMax, yMin, yMax,
  plotList,
  plotPointsGeneral, plotPointsPatched,
  textStyle,
  n, cUpperList, cLowerList, xCornerList, xIntList,
  cUpper, cLower,
  xLeft, xRight
 },
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.1 xTerm;
  yMax = 0.7;
  (* Plot points *)
  plotPointsGeneral = 3;
  plotPointsPatched = 2;
  (* Styles *)
  textStyle = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  (* Build a plot for each list of corners *)
  plotList = Table[
    (* *)
    n = patchedCornerNum[id];
    cUpperList = patchedCUpperList[id];
    cLowerList = patchedCLowerList[id];
    xCornerList = patchedCornerXList[id];
    xIntList = patchedIntXList[id];
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        AspectRatio -> Automatic,
        Frame -> None,
        ImageSize -> Automatic
      ],
      (* General boundaries: upper-branch(i) *)
      Table[
        cUpper = cUpperList[[i]];
        Plot[
          yTraUpper[cUpper][x],
          {x, 0, 1},
          PlotPoints -> plotPointsGeneral,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Background"]
        ]
      , {i, n}],
      (* General boundaries: lower-branch(i) *)
      Table[
        cLower = cLowerList[[i]];
        Plot[
          yTraLower[cLower][x],
          {x, 0, 1},
          PlotPoints -> plotPointsGeneral,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Background"]
        ]
      , {i, n}],
      (* Patched portions: upper-branch(i) *)
      Table[
        xLeft = xCornerList[[i]];
        xRight = If[i > 1,
          xIntList[[i - 1]],
          Max[xIntList]
        ];
        cUpper = cUpperList[[i]];
        Plot[
          yTraUpper[cUpper][x],
          {x, xLeft, xRight},
          PlotPoints -> plotPointsPatched,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
      , {i, n}],
      (* Patched portions: lower-branch(i) *)
      Table[
        xLeft = xCornerList[[i]];
        xRight = If[i < n,
          xIntList[[i]],
          Max[xIntList]
        ];
        cLower = cLowerList[[i]];
        Plot[
          yTraLower[cLower][x],
          {x, xLeft, xRight},
          PlotPoints -> plotPointsPatched,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
      , {i, n}],
      (* Critical terminal curve *)
      ParametricPlot[
        {xTerm, y}, {y, -yMax, yMax},
        PlotPoints -> 2,
        PlotStyle -> BoundaryTracingStyle["Terminal"]
      ],
      Graphics @ {
        Text[
          Italicise["x"] == 1
          , {1, 0}
          , {0, 1}
          , {0, 1}
        ] // textStyle
      },
      {}
    ]
  , {id, patchedIdList}]
  // GraphicsRow[#,
    Spacings -> {
      0.3 ImageSizeTextWidth,
      Automatic
    }
  ] &
  // Show[#, ImageSize -> ImageSizeTextWidth] &
] // Ex["plane-traced-boundaries-patched.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[
 {xTerm,
  xMin, xMax, yMin, yMax,
  plotList,
  plotPointsGeneral, plotPointsPatched,
  textStyle,
  n, cUpperList, cLowerList, xCornerList, xIntList,
  cUpper, cLower,
  xLeft, xRight
 },
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.1 xTerm;
  yMax = 0.7;
  (* Plot points *)
  plotPointsGeneral = 3;
  plotPointsPatched = 2;
  (* Styles *)
  textStyle = Style[#, 9] &;
  (* Build a plot for each list of corners *)
  plotList = Table[
    (* *)
    n = patchedCornerNum[id];
    cUpperList = patchedCUpperList[id];
    cLowerList = patchedCLowerList[id];
    xCornerList = patchedCornerXList[id];
    xIntList = patchedIntXList[id];
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        AspectRatio -> Automatic,
        Frame -> None,
        ImageSize -> Automatic
      ],
      (* General boundaries: upper-branch(i) *)
      Table[
        cUpper = cUpperList[[i]];
        Plot[
          yTraUpper[cUpper][x],
          {x, 0, 1},
          PlotPoints -> plotPointsGeneral,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Background"]
        ]
      , {i, n}],
      (* General boundaries: lower-branch(i) *)
      Table[
        cLower = cLowerList[[i]];
        Plot[
          yTraLower[cLower][x],
          {x, 0, 1},
          PlotPoints -> plotPointsGeneral,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Background"]
        ]
      , {i, n}],
      (* Patched portions: upper-branch(i) *)
      Table[
        xLeft = xCornerList[[i]];
        xRight = If[i > 1,
          xIntList[[i - 1]],
          Max[xIntList]
        ];
        cUpper = cUpperList[[i]];
        Plot[
          yTraUpper[cUpper][x],
          {x, xLeft, xRight},
          PlotPoints -> plotPointsPatched,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
        ]
      , {i, n}],
      (* Patched portions: lower-branch(i) *)
      Table[
        xLeft = xCornerList[[i]];
        xRight = If[i < n,
          xIntList[[i]],
          Max[xIntList]
        ];
        cLower = cLowerList[[i]];
        Plot[
          yTraLower[cLower][x],
          {x, xLeft, xRight},
          PlotPoints -> plotPointsPatched,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
        ]
      , {i, n}],
      (* Critical terminal curve *)
      ParametricPlot[
        {xTerm, y}, {y, -yMax, yMax},
        PlotPoints -> 2,
        PlotStyle -> BoundaryTracingStyle["Terminal"]
      ],
      Graphics @ {
        Text[
          Italicise["x"] == 1
          , {1, 0}
          , {0, 1.2}
          , {0, 1}
        ] // textStyle
      },
      {}
    ]
  , {id, patchedIdList}]
  // GraphicsRow[#,
    Spacings -> {
      0.2 ImageSizeTextWidthBeamer,
      Automatic
    }
  ] &
  // Show[#, ImageSize -> ImageSizeTextWidthBeamer] &
] // Ex["plane-traced-boundaries-patched-slides.pdf"];


(* ::Section:: *)
(*Figure: Domains (plane-domains.pdf)*)


(* ::Subsection:: *)
(*Domains*)


Module[
 {xTerm,
  xMin, xMax, yMin, yMax,
  plotPointsPatched,
  textStyle, textStyleLabel, textStyleBracket,
  plotList,
  n, cUpperList, cLowerList, xCornerList, xIntList,
  iRangeList, xBathList,
  iMin, iMax, xBath, yBathBottom, yBathTop,
  cUpper, cLower,
  xLeft, xRight,
  domainLabel,
  dummyForTrailingCommas
 },
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.1 xTerm;
  yMax = 0.7;
  (* Plot points *)
  plotPointsPatched = 2;
  (* Styles *)
  textStyle = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Build a plot for each list of corners *)
  plotList = Table[
    (* *)
    n = patchedCornerNum[id];
    cUpperList = patchedCUpperList[id];
    cLowerList = patchedCLowerList[id];
    xCornerList = patchedCornerXList[id];
    xIntList = patchedIntXList[id];
    iRangeList = domainCornerRangeList[id];
    xBathList = domainXBathList[id];
    (* Plot *)
    Show[
      EmptyAxes[{xMin, xMax}, {-yMax, yMax},
        AspectRatio -> Automatic,
        Axes -> None,
        ImageSize -> Automatic
      ],
      (* Critical terminal curve *)
      ParametricPlot[
        {xTerm, y}, {y, -yMax, yMax},
        PlotPoints -> 2,
        PlotStyle -> BoundaryTracingStyle["Terminal", "Background"]
      ],
      Graphics @ {Gray,
        Text[
          Italicise["x"] == 1
          , {1, 0}
          , {0, 1}
          , {0, 1}
        ] // textStyle
      },
      (* Domains *)
      Table[
        {iMin, iMax} = iRangeList[[j]];
        {
          (* Constant-temperature (Dirichlet) boundary *)
          xBath = xBathList[[j]];
          yBathBottom = yTraUpper[cUpperList[[iMin]]] @ xBath;
          yBathTop = yTraLower[cLowerList[[iMax]]] @ xBath;
          ParametricPlot[
            {xBath, y}, {y, yBathBottom, yBathTop},
            PlotPoints -> 2,
            PlotStyle -> BoundaryTracingStyle["Contour"]
          ],
          (* Patched portions: upper-branch(i) *)
          Table[
            xLeft = xCornerList[[i]];
            xRight = If[i > iMin, xIntList[[i - 1]], xBath];
            cUpper = cUpperList[[i]];
            Plot[
              yTraUpper[cUpper][x],
              {x, xLeft, xRight},
              PlotPoints -> plotPointsPatched,
              PlotRange -> {-yMax, yMax},
              PlotStyle -> BoundaryTracingStyle["Traced"]
            ]
          , {i, iMin, iMax}],
          (* Patched portions: lower-branch(i) *)
          Table[
            xLeft = xCornerList[[i]];
            xRight = If[i < iMax, xIntList[[i]], xBath];
            cLower = cLowerList[[i]];
            Plot[
              yTraLower[cLower][x],
              {x, xLeft, xRight},
              PlotPoints -> plotPointsPatched,
              PlotRange -> {-yMax, yMax},
              PlotStyle -> BoundaryTracingStyle["Traced"]
            ]
          , {i, iMin, iMax}]
        }
      , {j, Length[iRangeList]}],
      (* Domain labels *)
      domainLabel[str_][scaledX_, scaledY_] :=
        Text[
          Row @ {
            "(" // textStyleBracket,
            str,
            ")" // textStyleBracket,
          } // textStyleLabel
          , {Way[xMin, xMax, scaledX], scaledY * yMax}
        ];
      Graphics @ {
        Which[
          id == patchedIdList[[1]],
            domainLabel["a"][0.4, 0],
          id == patchedIdList[[2]],
            domainLabel["b"][0.5, 0],
          id == patchedIdList[[3]],
            {
              domainLabel["c"][0.3, 0.75],
              domainLabel["d"][0.11, 0.05],
              domainLabel["e"][0.5, -0.62],
              Nothing
            },
          True, {}
        ]
      },
      {}
    ]
  , {id, patchedIdList}];
  (* Combine *)
  GraphicsRow[plotList
    , ImageSize -> ImageSizeTextWidth
    , Spacings -> {0.3 ImageSizeTextWidth, Automatic}
  ]
] // Ex["plane-domains.pdf"]


(* ::Subsubsection:: *)
(*Version for slides*)


Module[
 {xTerm,
  xMin, xMax, yMin, yMax,
  plotPointsPatched,
  textStyle, textStyleLabel, textStyleBracket,
  plotList,
  n, cUpperList, cLowerList, xCornerList, xIntList,
  iRangeList, xBathList,
  iMin, iMax, xBath, yBathBottom, yBathTop,
  cUpper, cLower,
  xLeft, xRight,
  domainLabel,
  dummyForTrailingCommas
 },
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.1 xTerm;
  yMax = 0.7;
  (* Plot points *)
  plotPointsPatched = 2;
  (* Styles *)
  textStyle = Style[#, 9] &;
  textStyleLabel = Style[#, 10] &;
  textStyleBracket = Style[#, 10] &;
  (* Build a plot for each list of corners *)
  plotList = Table[
    (* *)
    n = patchedCornerNum[id];
    cUpperList = patchedCUpperList[id];
    cLowerList = patchedCLowerList[id];
    xCornerList = patchedCornerXList[id];
    xIntList = patchedIntXList[id];
    iRangeList = domainCornerRangeList[id];
    xBathList = domainXBathList[id];
    (* Plot *)
    Show[
      EmptyAxes[{xMin, xMax}, {-yMax, yMax},
        AspectRatio -> Automatic,
        Axes -> None,
        ImageSize -> Automatic
      ],
      (* Domains *)
      Table[
        {iMin, iMax} = iRangeList[[j]];
        {
          (* Constant-temperature (Dirichlet) boundary *)
          xBath = xBathList[[j]];
          yBathBottom = yTraUpper[cUpperList[[iMin]]] @ xBath;
          yBathTop = yTraLower[cLowerList[[iMax]]] @ xBath;
          ParametricPlot[
            {xBath, y}, {y, yBathBottom, yBathTop},
            PlotPoints -> 2,
            PlotStyle -> Directive[BoundaryTracingStyle["Contour"], SlidesStyle["Source"]]
          ],
          (* Patched portions: upper-branch(i) *)
          Table[
            xLeft = xCornerList[[i]];
            xRight = If[i > iMin, xIntList[[i - 1]], xBath];
            cUpper = cUpperList[[i]];
            Plot[
              yTraUpper[cUpper][x],
              {x, xLeft, xRight},
              PlotPoints -> plotPointsPatched,
              PlotRange -> {-yMax, yMax},
              PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
            ]
          , {i, iMin, iMax}],
          (* Patched portions: lower-branch(i) *)
          Table[
            xLeft = xCornerList[[i]];
            xRight = If[i < iMax, xIntList[[i]], xBath];
            cLower = cLowerList[[i]];
            Plot[
              yTraLower[cLower][x],
              {x, xLeft, xRight},
              PlotPoints -> plotPointsPatched,
              PlotRange -> {-yMax, yMax},
              PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
            ]
          , {i, iMin, iMax}]
        }
      , {j, Length[iRangeList]}],
      (* Domain labels *)
      domainLabel[str_][scaledX_, scaledY_] :=
        Text[
          Row @ {
            "(" // textStyleBracket,
            str,
            ")" // textStyleBracket,
          } // textStyleLabel
          , {Way[xMin, xMax, scaledX], scaledY * yMax}
        ];
      Graphics @ {
        Which[
          id == patchedIdList[[1]],
            domainLabel["a"][0.4, 0],
          id == patchedIdList[[2]],
            domainLabel["b"][0.5, 0],
          id == patchedIdList[[3]],
            {
              domainLabel["c"][0.3, 0.75],
              domainLabel["d"][0.11, 0.05],
              domainLabel["e"][0.5, -0.62],
              Nothing
            },
          True, {}
        ]
      },
      {}
    ]
  , {id, patchedIdList}]
  // GraphicsRow[#,
    Spacings -> {
      0.2 ImageSizeTextWidthBeamer,
      Automatic
    }
  ] &
  // Show[#, ImageSize -> ImageSizeTextWidthBeamer] &
] // Ex["plane-domains-slides.pdf"];


(* ::Subsection:: *)
(*Legend*)


Module[{legendCurves},
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle @* ReleaseHold /@
        {"Traced", "Contour", Hold @ Sequence["Terminal", "Background"]},
      {"radiation", "constant temperature", "critical terminal curve"}
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
    ];
  GraphicsGrid[{legendCurves}
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.11
    , Spacings -> {{0, -0.24, 0.1} ImageSizeTextWidth, 0}
  ]
] // Ex["plane-domains-legend.pdf"]


(* ::Section:: *)
(*Figure: Self-viewing radiation ratio (plane-self-viewing-ratio.pdf)*)


Module[
  {
    labelPlacementFromInterval, intervalList,
    numPoints, data,
    x1, x2, xValues,
    textStyle, textStyleBracket,
    intervalText,
    proportion, offset, xLabel, rLabel,
    dummyForTrailingCommas
  },
  (* List of intervals {x1, x2} to plot the ratio for *)
  labelPlacementFromInterval = Association[
    {0, 0.3} -> {0.8, {1.35, 0}},
    {0.2, 0.5} -> {0.1, {-1.2, -0.7}},
    {0.5, 0.8} -> {0, {-0.55, -0.9}},
    {0.7, 1} -> {0.3, {0, -1.4}},
    {0.3, 0.4} -> {1, {-1, 0.6}},
    {0.5, 0.6} -> {0.7, {-0.8, 0.8}},
    {0.65, 0.75} -> {0.7, {-1, 0.65}},
    {0.8, 0.9} -> {0.8, {-1, 0.6}},
    Nothing
  ];
  intervalList = Keys[labelPlacementFromInterval];
  (* Number of points for sampling *)
  numPoints = 50;
  (* Compute data points to plot *)
  Table[
    {x1, x2} = interval;
    xValues = Subdivide[x1, x2, numPoints];
    If[x1 == 0,
      xValues = Join[
        Subdivide[x1, xValues[[4]], 20],
        xValues[[5 ;;]]
      ];
    ];
    If[x2 == 1,
      xValues = Join[
        xValues[[;; -5]],
        Subdivide[xValues[[-4]], x2, 20]
      ];
    ];
    data[interval] =
      Cases[
        Table[{x, boundaryRatio[x1, x2][x]}, {x, xValues}],
        {_?NumericQ, _?NumericQ}
      ] // Quiet;
    , {interval, intervalList}
  ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  intervalText[interval_] :=
    Text[
      {x1, x2} = interval;
      {proportion, offset} = labelPlacementFromInterval[interval];
      Row @ {
        "(" // textStyleBracket,
        "\[NegativeVeryThinSpace]",
        x1,
        ",\[ThinSpace]",
        x2,
        ")" // textStyleBracket
      } // textStyle
      ,
        xLabel = Way[x1, x2, proportion];
        rLabel = Interpolation[data[interval], xLabel];
        {xLabel, Log[rLabel]}
      , offset
    ];
  (* Make plot *)
  ListLogPlot[
    Table[data[interval], {interval, intervalList}]
    , AxesLabel -> {
        Italicise["x"] // Margined @ {{0, 1}, {5, 0}},
        Italicise["R"]
      }
    , Epilog -> {
        {
          BoundaryTracingStyle["Contour"],
          Line @ {{0, #}, {1, #}} & [1/100 // Log]
        },
        Table[intervalText[interval], {interval, intervalList}],
        {}
      }
    , ImageSize -> 0.55 ImageSizeTextWidth
    , Joined -> True
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {{0, 1}, {All, 10}}
    , PlotRangeClipping -> False
    , PlotStyle -> Black
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["plane-self-viewing-ratio.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    labelPlacementFromInterval, intervalList,
    numPoints, data,
    x1, x2, xValues,
    textStyle, textStyleBracket,
    intervalText,
    proportion, offset, xLabel, rLabel,
    dummyForTrailingCommas
  },
  (* List of intervals {x1, x2} to plot the ratio for *)
  labelPlacementFromInterval = Association[
    {0, 0.3} -> {0.8, {1.35, 0}},
    {0.2, 0.5} -> {0.1, {-1.2, -0.7}},
    {0.5, 0.8} -> {0, {-0.55, -0.9}},
    {0.7, 1} -> {0.3, {0, -1.4}},
    {0.3, 0.4} -> {1, {-1, 0.6}},
    {0.5, 0.6} -> {0.7, {-0.8, 0.8}},
    {0.65, 0.75} -> {0.7, {-1, 0.65}},
    {0.8, 0.9} -> {0.8, {-1, 0.6}},
    Nothing
  ];
  intervalList = Keys[labelPlacementFromInterval];
  (* Number of points for sampling *)
  numPoints = 50;
  (* Compute data points to plot *)
  Table[
    {x1, x2} = interval;
    xValues = Subdivide[x1, x2, numPoints];
    If[x1 == 0,
      xValues = Join[
        Subdivide[x1, xValues[[4]], 20],
        xValues[[5 ;;]]
      ];
    ];
    If[x2 == 1,
      xValues = Join[
        xValues[[;; -5]],
        Subdivide[xValues[[-4]], x2, 20]
      ];
    ];
    data[interval] =
      Cases[
        Table[{x, boundaryRatio[x1, x2][x]}, {x, xValues}],
        {_?NumericQ, _?NumericQ}
      ] // Quiet;
    , {interval, intervalList}
  ];
  (* Text style *)
  textStyle = Style[#, 6] &;
  textStyleBracket = Style[#, 6] &;
  intervalText[interval_] :=
    Text[
      {x1, x2} = interval;
      {proportion, offset} = labelPlacementFromInterval[interval];
      Row @ {
        "(" // textStyleBracket,
        "\[NegativeVeryThinSpace]",
        x1,
        ",\[ThinSpace]",
        x2,
        ")" // textStyleBracket
      } // textStyle
      ,
        xLabel = Way[x1, x2, proportion];
        rLabel = Interpolation[data[interval], xLabel];
        {xLabel, Log[rLabel]}
      , offset
    ];
  (* Make plot *)
  ListLogPlot[
    Table[data[interval], {interval, intervalList}]
    , AxesLabel -> {
        Italicise["x"] // Margined @ {{0, 1}, {5, 0}},
        Italicise["R"]
      }
    , Epilog -> {
        {
          BoundaryTracingStyle["Contour"],
          Line @ {{0, #}, {1, #}} & [1/100 // Log]
        },
        Table[intervalText[interval], {interval, intervalList}],
        {}
      }
    , ImageSize -> 0.5 ImageSizeTextWidthBeamer
    , Joined -> True
    , LabelStyle -> 10
    , PlotRange -> {{0, 1}, {All, 10}}
    , PlotRangeClipping -> False
    , PlotStyle -> SlidesStyle["Boundary"]
    , TicksStyle -> 7
  ]
] // Ex["plane-self-viewing-ratio-slides.pdf"];


(* ::Section:: *)
(*Figure: example fin dimensions (plane-fin-dimensions)*)


Module[
  {
    x1, x2, y, y2,
    dimensionMarkerStyle,
    lengthMarkerY, thicknessMarkerX,
    textStyle, textStyleBracket,
    dummyForTrailingCommas
  },
  (* Abbreviations *)
  {x1, x2} = {exampleX1, exampleX2};
  y = exampleYTraced;
  y2 = exampleY2;
  (* Dimension markers *)
  dimensionMarkerStyle = Arrowheads @ {-Small, Small};
  lengthMarkerY = -1.5 y2;
  thicknessMarkerX = Way[x1, x2, 1.045];
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] & @* LaTeXStyle;
  (* Make diagram *)
  Show[
    (* Radiation boundaries *)
    Plot[{-1, 1} y[x], {x, x1, x2}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Constant-temperature boundary *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {{x2, -y2}, {x2, +y2}}
    },
    (* Length L marker *)
    Graphics @ {dimensionMarkerStyle,
      Arrow @ {
        {x1, lengthMarkerY},
        {x2, lengthMarkerY}
      },
      Text[
        Italicise["L"] // textStyle
        , {Way[x1, x2], lengthMarkerY}
        , {0, 1}
      ]
    },
    (* Thickness H marker *)
    Graphics @ {dimensionMarkerStyle,
      Arrow @ {
        {thicknessMarkerX, -y2},
        {thicknessMarkerX, y2}
      },
      Text[
        Italicise["H"] // textStyle
        , {thicknessMarkerX, 0}
        , {-1.8, -0.15}
      ]
    },
    (* Tip coordinate *)
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], 1],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        } // textStyle
        , {x1, 0}
        , {1.5, -0.15}
      ]
    },
    (* Base coordinate *)
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], 2],
          ",\[VeryThinSpace]\[VeryThinSpace]",
          Subscript[Italicise["y"], 2],
          ")" // textStyleBracket
        } // textStyle
        , {x2, y2}
        , {0, -1.4}
      ]
    },
    {}
    , AspectRatio -> Automatic
    , Axes -> False
    , ImageSize -> 0.72 ImageSizeTextWidth
    , PlotRange -> All
    , PlotRangeClipping -> False
    , PlotRangePadding -> {Automatic, {Automatic, Scaled[0.4]}}
  ]
] // Ex["plane-fin-dimensions.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    x1, x2, y, y2,
    dimensionMarkerStyle,
    lengthMarkerY, thicknessMarkerX,
    textStyle, textStyleBracket,
    dummyForTrailingCommas
  },
  (* Abbreviations *)
  {x1, x2} = {exampleX1, exampleX2};
  y = exampleYTraced;
  y2 = exampleY2;
  (* Dimension markers *)
  dimensionMarkerStyle = Arrowheads @ {-Small, Small};
  lengthMarkerY = -1.5 y2;
  thicknessMarkerX = Way[x1, x2, 1.045];
  textStyle = Style[#, 8] &;
  textStyleBracket = Style[#, 8] &;
  (* Make diagram *)
  Show[
    (* Constant-temperature boundary *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      SlidesStyle["Source"],
      Line @ {{x2, -y2}, {x2, +y2}}
    },
    (* Radiation boundaries *)
    Plot[{-1, 1} y[x], {x, x1, x2}
      , PlotPoints -> 2
      , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
    ],
    (* Length L marker *)
    Graphics @ {dimensionMarkerStyle,
      Arrow @ {
        {x1, lengthMarkerY},
        {x2, lengthMarkerY}
      },
      Text[
        Italicise["L"] // textStyle
        , {Way[x1, x2], lengthMarkerY}
        , {0, 1}
      ]
    },
    (* Thickness H marker *)
    Graphics @ {dimensionMarkerStyle,
      Arrow @ {
        {thicknessMarkerX, -y2},
        {thicknessMarkerX, y2}
      },
      Text[
        Italicise["H"] // textStyle
        , {thicknessMarkerX, 0}
        , {-1.85, -0.15}
      ]
    },
    (* Tip coordinate *)
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], 1],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        } // textStyle
        , {x1, 0}
        , {1.5, -0.15}
      ]
    },
    (* Base coordinate *)
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], 2],
          ",\[VeryThinSpace]\[VeryThinSpace]",
          Subscript[Italicise["y"], 2],
          ")" // textStyleBracket
        } // textStyle
        , {x2, y2}
        , {0, -1.4}
      ]
    },
    {}
    , AspectRatio -> Automatic
    , Axes -> False
    , ImageSize -> 0.5 ImageSizeTextWidthBeamer
    , PlotRange -> All
    , PlotRangeClipping -> False
    , PlotRangePadding -> Automatic
  ]
] // Ex["plane-fin-dimensions-slides.pdf"];


(* ::Section:: *)
(*Figure: example fin temperature (plane-fin-temperature)*)


Module[{textStyle},
textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
Plot[
  {exampleT1[l], exampleT2[l]} - celsiusOffset
  , {l, 0, 4}
  , AxesLabel -> {
      (* L / m *)
      SeparatedRow["VeryThin"][
        Italicise["L"],
        Style["/", Magnification -> 1.25],
        "m"
      ] // Margined @ {{-2, 0}, {2, 0}},
      (* t / \[Degree]C *)
      SeparatedRow["VeryThin"][
        Italicise["t"],
        Style["/", Magnification -> 1.25],
        SeparatedRow["VeryThin"][
          AdjustmentBox[
            Style["\[Degree]", Magnification -> 1.2]
            , BoxBaselineShift -> -0.2
          ],
          "C"
        ] // DisplayForm
      ] // Margined @ {{0, 0}, {-5, -5}}
    }
  , Epilog -> {
      Text[
        "base" // textStyle
        , {#, exampleT2[#] - celsiusOffset} & [1]
        , {-1, -1}
      ],
      Text[
        "tip" // textStyle
        , {#, exampleT1[#] - celsiusOffset} & [0.9]
        , {1.3, 0.4}
      ],
      {}
    }
  , ImageSize -> 0.48 ImageSizeTextWidth
  , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
  , PlotStyle -> Black
  , PlotPoints -> 3
  , TicksStyle -> LabelSize["Tick"]
]
] // Ex["plane-fin-temperature.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[{textStyle},
textStyle = Style[#, 10] &;
Plot[
  {exampleT1[l], exampleT2[l]} - celsiusOffset
  , {l, 0, 4}
  , AspectRatio -> 0.8
  , AxesLabel -> {
      (* L / m *)
      SeparatedRow["Thin"][
        Italicise["L"],
        Style["/", Magnification -> 1],
        "m"
      ] // Margined @ {{-2, 0}, {2, 0}},
      (* t / \[Degree]C *)
      SeparatedRow["Thin"][
        Italicise["t\[ThinSpace]"],
        Style["/", Magnification -> 1],
        SeparatedRow[][Style["\[Degree]", Magnification -> 1], "C"]
      ] // Margined @ {{0, 0}, {-10, -5}}
    }
  , Epilog -> {
      Text[
        "base" // textStyle
        , {#, exampleT2[#] - celsiusOffset} & [1]
        , {-1, -1}
      ],
      Text[
        "tip" // textStyle
        , {#, exampleT1[#] - celsiusOffset} & [0.9]
        , {1.3, 0.4}
      ],
      {}
    }
  , ImageSize -> 0.5 ImageSizeTextWidthBeamer
  , LabelStyle -> 10
  , PlotLabel -> "Aluminium"
  , PlotStyle -> Black
  , PlotPoints -> 3
  , TicksStyle -> 8
]
] // Ex["plane-fin-temperature-slides.pdf"];


(* ::Section:: *)
(*Figure: example fin power per length (plane-fin-power-per-length)*)


Plot[
  exampleP[l] / 10^3
  , {l, 0, 6}
  , AxesLabel -> {
      (* L / m *)
      SeparatedRow["VeryThin"][
        Italicise["L"],
        Style["/", Magnification -> 1.25],
        "m"
      ] // Margined @ {{-2, 0}, {2, 0}},
      (* p / (W/m) *)
      SeparatedRow["VeryThin"][
        Italicise["p"],
        Style["/", Magnification -> 1.25],
        Row @ {
          Style["(", Magnification -> 1.4],
          SeparatedRow["Thin"][
            "k\[NegativeThinSpace]\[VeryThinSpace]W",
            Superscript[
              "m",
              RowBox @ List @ AdjustmentBox[
                Row @ {"-", Style[1, Magnification -> 0.9]}
                , BoxBaselineShift -> 0.12
              ] // DisplayForm
            ]
          ],
          Style[")", Magnification -> 1.4]
        }
      ] // Margined @ {{0, 0}, {-3, -5}}
    }
  , ImageSize -> 0.48 ImageSizeTextWidth
  , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
  , PlotPoints -> 3
  , PlotRange -> {0, Automatic}
  , PlotStyle -> Black
  , TicksStyle -> LabelSize["Tick"]
] // Ex["plane-fin-power-per-length.pdf"]


(* ::Section:: *)
(*Figure: directional dependence example fin vs strip (plane-directional-dependence-fin-strip)*)


Module[
  {
    x1, x2,
    yDer,
    yDerMaxFin, phiMaxFin,
    finDependence, finConstant, finDependenceNormalised,
    stripDependence, stripConstant, stripDependenceNormalised,
    textStyle,
    horizontalTicks, verticalTicks,
    dummyForTrailingCommas
  },
  (* Abbreviations *)
  {x1, x2} = {exampleX1, exampleX2};
  (* Fin angular range *)
  yDer = -yTraDer[#] &;
  yDerMaxFin = yDer[x2];
  phiMaxFin = Pi - ArcTan[yDerMaxFin];
  (* Fin directional dependence *)
  finDependence[phi_] := directionalDependence[x1, x2] @ Abs[phi];
  finConstant = NIntegrate[finDependence[phi], {phi, -phiMaxFin, phiMaxFin}];
  finDependenceNormalised[phi_] := finDependence[phi] / finConstant;
  (* Strip directional dependence *)
  stripDependence[phi_] := Cos[phi];
  stripConstant = Integrate[stripDependence[phi], {phi, -Pi/2, Pi/2}];
  stripDependenceNormalised[phi_] := stripDependence[phi] / stripConstant;
  (* Etc. *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  horizontalTicks =
    Table[
      {
        gpd * Degree,
        SeparatedRow["VeryThin"][
          gpd,
          AdjustmentBox[
            Style["\[Degree]", Magnification -> 1.2]
            , BoxBaselineShift -> -0.1
          ]
        ]
          // DisplayForm
          // Margined @ {{0, 0}, {0, -7}}
      }
      , {gpd, -180, 180, 90}
    ];
  verticalTicks = Automatic;
  (* Make plot *)
  Plot[
    {stripDependenceNormalised, finDependenceNormalised}[phi]
      // Through
      // Evaluate
    , {phi, -Pi, Pi}
    , AspectRatio -> 0.4
    , AxesLabel -> {
        "\[CurlyPhi]" // LaTeXStyle // Margined @ {{0, 0}, {4, 0}},
        SeparatedRow["Thick"]["Normalised", Italicise["I"]]
          // Margined @ {{0, 0}, {-4, 0}}
      }
    , Epilog -> {
        Text[
          "fin" // textStyle
          , {#, finDependenceNormalised[#]} & [-125 Degree]
          , {2.4, 0}
        ],
        Text[
          "strip" // textStyle
          , {#, stripDependenceNormalised[#]} & [20 Degree]
          , {-1.8, 0}
        ],
        {}
      }
    , ImageSize -> 0.65 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotPoints -> 2
    , PlotRange -> {0, Automatic}
    , PlotStyle -> Black
    , Ticks -> {horizontalTicks, verticalTicks}
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["plane-directional-dependence-fin-strip.pdf"]


(* ::Section:: *)
(*Figure: example fin directional dependence geometry (plane-directional-dependence-geometry-*)*)


Module[
  {
    x1, x2, y, y2,
    xFinCentre,
    yDer,
    yDerMaxFin, phiMaxFin,
    finDependence, finConstant,
    finDependenceNormalised,
    yMaxStrip,
    stripDependence, stripConstant, stripDependenceNormalised,
    rho,
    gammaCorrect,
    dependenceColourFunction,
    dependenceThicknessFunction, replaceData,
    thicken,
    commonPlot,
      rMinDirectionReference, rMaxDirectionReference,
      rMidDirectionReference,
      phiMarker,
      textStyle,
    directionalOptions,
    radiationArrowheadSize,
    radiationArrowSize, radiationArrowClearance,
    radiationXValuesFin, radiationYValuesStrip,
    dummyForTrailingCommas
  },
  (* Fin *)
  {x1, x2} = {exampleX1, exampleX2};
  y = exampleYTraced;
  y2 = exampleY2;
  xFinCentre = Way[x1, x2];
  yDer = -yTraDer[#] &;
  yDerMaxFin = yDer[x2];
  phiMaxFin = Pi - ArcTan[yDerMaxFin];
  finDependence[phi_] := directionalDependence[x1, x2] @ Abs[phi];
  finConstant = NIntegrate[finDependence[phi], {phi, -phiMaxFin, phiMaxFin}];
  finDependenceNormalised[phi_] := finDependence[phi] / finConstant;
  (* Strip *)
  yMaxStrip = 1.5 y2;
  stripDependence[phi_] := Cos[phi];
  stripConstant = Integrate[stripDependence[phi], {phi, -Pi/2, Pi/2}];
  stripDependenceNormalised[phi_] := stripDependence[phi] / stripConstant;
  (* Cylinder at infinity *)
  rho = 0.97 (x2 - x1);
  (* Dependence colouring *)
  gammaCorrect[level_] := Clip[level, {0, 1}] ^ 0.3;
  dependenceColourFunction[dependence_, phi_, {minPhi_, maxPhi_}] :=
    GrayLevel @ gammaCorrect[1 - Rescale[dependence[phi], dependence /@ {minPhi, maxPhi}]];
  (* Dependence thickening *)
  (*
    `thicken` below is a modified version of Kuba's `thick$color` function
    in <https://mathematica.stackexchange.com/a/28207>.
    Re-use permission granted in comments, see archived version:
    <https://web.archive.org/web/20210315151008/https://mathematica.stackexchange.com/questions/28202/vary-the-thickness-of-a-plotted-function/>
  *)
  dependenceThicknessFunction[dependence_][x_, y_] := dependence[Abs @ ArcTan[-x, y]];
  thicken[referenceThickness_][dependence_] := ReplaceAll[
    GraphicsComplex[points_List, data_List, options : OptionsPattern[__]] :>
    GraphicsComplex[
      points,
      data /. {
       Line[linePoints_List, OptionsPattern[]] :>
         Sequence @@ (
           {
             AbsoluteThickness[
               referenceThickness *
               dependenceThicknessFunction[dependence] @@
                 (Mean /@ Transpose @ points[[#]])
             ],
             Line[#, VertexColors -> Automatic]
           } & /@
             Partition[linePoints, 2, 1]
         )
      },
      options
    ]
  ];
  (* Make plots *)
  rMinDirectionReference = 1.15 rho;
  rMaxDirectionReference = 1.25 rho;
  rMidDirectionReference = Way[rMinDirectionReference, rMaxDirectionReference];
  phiMarker = 17 Degree;
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  commonPlot = EmptyFrame[{-rho, rho}, {-rho, rho}
    , Epilog -> {
        (* Reference direction (minus x) *)
        Line @ {{-rMinDirectionReference, 0}, {-rMaxDirectionReference, 0}},
        (* Arc *)
        Circle[{0, 0}, rMidDirectionReference, {Pi, Pi - phiMarker}],
        {Arrowheads[0.04],
          Arrow @ Table[
            XYPolar[rMidDirectionReference, Pi - phi]
            , {phi, Subdivide[0.7, 1.1, 2] phiMarker}
          ]
        },
        (* Label *)
        Text[
          "\[CurlyPhi]" // textStyle
          , XYPolar[rMidDirectionReference, Pi - phiMarker]
          , {-0.5, -1.3}
        ],
        {}
      }
    , Frame -> None
    , ImageSize -> 0.33 ImageSizeTextWidth
    , PlotRangePadding -> {Scaled /@ {0.15, 0.03}, Scaled[0.03]}
  ];
  directionalOptions = {Nothing
    , ColorFunctionScaling -> False
    , PlotPoints -> 3
  };
  radiationArrowheadSize = 0.025;
  radiationArrowSize = 0.25 (x2 - x1);
  radiationArrowClearance = 0.35 radiationArrowSize;
  radiationXValuesFin = Way[x1, x2, Subdivide[0.1, 0.9, 3]];
  radiationYValuesStrip = {-1, 1} 0.7 yMaxStrip;
  {
    (* Fin *)
    Show[
      commonPlot,
      (* Directional dependence *)
      ParametricPlot[
        XYPolar[rho, Pi - phi]
        , {phi, -phiMaxFin, phiMaxFin}
        , ColorFunction -> Function[{x, y, phi},
            dependenceColourFunction[
              finDependenceNormalised, phi,
              {phiMaxFin, Pi/2}
            ]
          ]
        , directionalOptions // Evaluate
      ] // thicken[13][finDependenceNormalised],
      (* Radiation boundaries *)
      ParametricPlot[
        {{x - xFinCentre, -y[x]}, {x - xFinCentre, y[x]}}
        , {x, x1, x2}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
      (* Radiation arrows *)
      Graphics @ {Arrowheads[radiationArrowheadSize],
        Table[
          SquigglyArrow[
            {x - xFinCentre, sign * y[x]} + XYPolar[radiationArrowClearance, #],
            #,
            radiationArrowSize
          ] & [
            sign * (ArcTan[yDer[x]] + Pi/2)
          ]
          , {x, radiationXValuesFin}
          , {sign, {-1, 1}}
        ]
      },
      {}
    ] // Ex["plane-directional-dependence-geometry-fin.pdf"]
    ,
    (* Strip *)
    Show[
      commonPlot,
      (* Directional dependence *)
      ParametricPlot[
        XYPolar[rho, Pi - phi]
        , {phi, -Pi/2, Pi/2}
        , ColorFunction -> Function[{x, y, phi},
            dependenceColourFunction[
              stripDependenceNormalised, phi,
              {Pi/2, 0}
            ]
          ]
        , directionalOptions // Evaluate
      ] // thicken[8][stripDependenceNormalised],
      (* Radiation strip *)
      Graphics @ {BoundaryTracingStyle["Traced"],
        Line @ {{0, -yMaxStrip}, {0, +yMaxStrip}}
      },
      (* Radiation arrows *)
      Graphics @ {Arrowheads[1.15 radiationArrowheadSize],
        Table[
          SquigglyArrow[
            {-radiationArrowClearance, y},
            -Pi,
            1.15 radiationArrowSize
          ]
          , {y, radiationYValuesStrip}
        ]
      },
      {}
    ] // Ex["plane-directional-dependence-geometry-strip.pdf"]
  }
]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    x1, x2, y, y2,
    xFinCentre,
    yDer,
    yDerMaxFin, phiMaxFin,
    finDependence, finConstant,
    finDependenceNormalised,
    yMaxStrip,
    stripDependence, stripConstant, stripDependenceNormalised,
    rho,
    gammaCorrect,
    dependenceColourFunction,
    dependenceThicknessFunction, replaceData,
    thicken,
    commonPlot,
      rMinDirectionReference, rMaxDirectionReference,
      rMidDirectionReference,
      phiMarker,
      textStyle,
    directionalOptions,
    radiationArrowheadSize,
    radiationArrowSize, radiationArrowClearance,
    radiationXValuesFin, radiationYValuesStrip,
    dummyForTrailingCommas
  },
  (* Fin *)
  {x1, x2} = {exampleX1, exampleX2};
  y = exampleYTraced;
  y2 = exampleY2;
  xFinCentre = Way[x1, x2];
  yDer = -yTraDer[#] &;
  yDerMaxFin = yDer[x2];
  phiMaxFin = Pi - ArcTan[yDerMaxFin];
  finDependence[phi_] := directionalDependence[x1, x2] @ Abs[phi];
  finConstant = NIntegrate[finDependence[phi], {phi, -phiMaxFin, phiMaxFin}];
  finDependenceNormalised[phi_] := finDependence[phi] / finConstant;
  (* Strip *)
  yMaxStrip = 1.5 y2;
  stripDependence[phi_] := Cos[phi];
  stripConstant = Integrate[stripDependence[phi], {phi, -Pi/2, Pi/2}];
  stripDependenceNormalised[phi_] := stripDependence[phi] / stripConstant;
  (* Cylinder at infinity *)
  rho = 0.97 (x2 - x1);
  (* Dependence colouring *)
  gammaCorrect[level_] := Clip[level, {0, 1}] ^ 0.3;
  dependenceColourFunction[dependence_, phi_, {minPhi_, maxPhi_}] :=
    GrayLevel @ gammaCorrect[1 - Rescale[dependence[phi], dependence /@ {minPhi, maxPhi}]];
  (* Dependence thickening *)
  (*
    `thicken` below is a modified version of Kuba's `thick$color` function
    in <https://mathematica.stackexchange.com/a/28207>.
    Re-use permission granted in comments, see archived version:
    <https://web.archive.org/web/20210315151008/https://mathematica.stackexchange.com/questions/28202/vary-the-thickness-of-a-plotted-function/>
  *)
  dependenceThicknessFunction[dependence_][x_, y_] := dependence[Abs @ ArcTan[-x, y]];
  thicken[referenceThickness_][dependence_] := ReplaceAll[
    GraphicsComplex[points_List, data_List, options : OptionsPattern[__]] :>
    GraphicsComplex[
      points,
      data /. {
       Line[linePoints_List, OptionsPattern[]] :>
         Sequence @@ (
           {
             AbsoluteThickness[
               referenceThickness *
               dependenceThicknessFunction[dependence] @@
                 (Mean /@ Transpose @ points[[#]])
             ],
             Line[#, VertexColors -> Automatic]
           } & /@
             Partition[linePoints, 2, 1]
         )
      },
      options
    ]
  ];
  (* Make plots *)
  rMinDirectionReference = 1.15 rho;
  rMaxDirectionReference = 1.25 rho;
  rMidDirectionReference = Way[rMinDirectionReference, rMaxDirectionReference];
  phiMarker = 17 Degree;
  textStyle = Style[#, 10] &;
  commonPlot = EmptyFrame[{-rho, rho}, {-rho, rho}
    , Epilog -> {
        (* Reference direction (minus x) *)
        Line @ {{-rMinDirectionReference, 0}, {-rMaxDirectionReference, 0}},
        (* Arc *)
        Circle[{0, 0}, rMidDirectionReference, {Pi, Pi - phiMarker}],
        {Arrowheads[0.04],
          Arrow @ Table[
            XYPolar[rMidDirectionReference, Pi - phi]
            , {phi, Subdivide[0.7, 1.1, 2] phiMarker}
          ]
        },
        (* Label *)
        Text[
          "\[CurlyPhi]" // textStyle
          , XYPolar[rMidDirectionReference, Pi - phiMarker]
          , {-0.5, -1.7}
        ],
        {}
      }
    , Frame -> None
    , ImageSize -> 0.8 * 0.5 ImageSizeTextWidthBeamer
    , LabelStyle -> Directive[11, SlidesStyle["Boundary"]]
    , PlotRangePadding -> {Scaled /@ {0.15, 0.03}, Scaled[0.03]}
  ];
  directionalOptions = {Nothing
    , ColorFunctionScaling -> False
    , PlotPoints -> 3
  };
  radiationArrowheadSize = 0.025;
  radiationArrowSize = 0.25 (x2 - x1);
  radiationArrowClearance = 0.35 radiationArrowSize;
  radiationXValuesFin = Way[x1, x2, Subdivide[0.1, 0.9, 3]];
  radiationYValuesStrip = {-1, 1} 0.7 yMaxStrip;
  {
    (* Fin *)
    Show[
      commonPlot,
      (* Directional dependence *)
      ParametricPlot[
        XYPolar[rho, Pi - phi]
        , {phi, -phiMaxFin, phiMaxFin}
        , ColorFunction -> Function[{x, y, phi},
            dependenceColourFunction[
              finDependenceNormalised, phi,
              {phiMaxFin, Pi/2}
            ]
          ]
        , directionalOptions // Evaluate
      ] // thicken[13][finDependenceNormalised],
      (* Radiation boundaries *)
      ParametricPlot[
        {{x - xFinCentre, -y[x]}, {x - xFinCentre, y[x]}}
        , {x, x1, x2}
        , PlotPoints -> 2
        , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
      ],
      (* Radiation arrows *)
      Graphics @ {Arrowheads[radiationArrowheadSize],
        Table[
          SquigglyArrow[
            {x - xFinCentre, sign * y[x]} + XYPolar[radiationArrowClearance, #],
            #,
            radiationArrowSize
          ] & [
            sign * (ArcTan[yDer[x]] + Pi/2)
          ]
          , {x, radiationXValuesFin}
          , {sign, {-1, 1}}
        ]
      },
      {}
      , PlotLabel -> ("Fin" // Margined @ {{0, 0}, {5, 0}})
    ] // Ex["plane-directional-dependence-geometry-fin-slides.pdf"]
    ,
    (* Strip *)
    Show[
      commonPlot,
      (* Directional dependence *)
      ParametricPlot[
        XYPolar[rho, Pi - phi]
        , {phi, -Pi/2, Pi/2}
        , ColorFunction -> Function[{x, y, phi},
            dependenceColourFunction[
              stripDependenceNormalised, phi,
              {Pi/2, 0}
            ]
          ]
        , directionalOptions // Evaluate
      ] // thicken[8][stripDependenceNormalised],
      (* Radiation strip *)
      Graphics @ {Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]],
        Line @ {{0, -yMaxStrip}, {0, +yMaxStrip}}
      },
      (* Radiation arrows *)
      Graphics @ {Arrowheads[1.15 radiationArrowheadSize],
        Table[
          SquigglyArrow[
            {-radiationArrowClearance, y},
            -Pi,
            1.15 radiationArrowSize
          ]
          , {y, radiationYValuesStrip}
        ]
      },
      {}
      , PlotLabel -> ("Strip" // Margined @ {{0, 0}, {5, 0}})
    ] // Ex["plane-directional-dependence-geometry-strip-slides.pdf"]
  }
];
