(* ::Package:: *)

SetDirectory @ ParentDirectory @ NotebookDirectory[]
<< Conway`
<< Curvilinear`
SetDirectory @ NotebookDirectory[]


(* ::Section:: *)
(*Testing Conway.wl*)


Plot[x^2, {x, 0, 1},
  PlotLabel -> BoxedLabel["This is a boxed label"],
  PlotOptions[Axes] // Evaluate
]


DefaultColours


Module[{y, iFun, iFunReInt},
  iFun = NDSolveValue[
    {y'[x] == 0, y[0] == 1}, y, {x, 0, 1},
    MaxStepFraction -> 1/2000,
    NoExtrapolation
  ];
  iFunReInt = ReInterpolate[iFun];
  {
    DomainStart[iFun],
    DomainEnd[iFun],
    iFun[1 + DomainEnd[iFun]],
    iFun["Grid"] // Length,
    iFunReInt["Grid"] // Length
  }
]


EmptyAxes[{0, 1}, {0, 1}]


EmptyFrame[{0, 1}, {0, 1}]


Module[{var = 0},
  {
    "{{1 + 1}} formats as {1 + 1}.",
    "{{var}} formats as {var}."
  } // FString // TableForm
]


Block[{blockVar = 0},
  "Block variables also work: {blockVar}" // FString
]


Module[{outerVar = 0},
  Module[{innerVar = 1},
    "Nested Modules: outerVar is {outerVar} and innerVar is {innerVar}."
      // FString
  ]
]


"a + beta" // PrettyString["beta" -> "\[Beta]"]


(* Compare displayed and exported LaTeXStyle: *)
Module[{expressions},
  expressions = LaTeXStyle @ Flatten @ {
    $Version,
    Join[
      CharacterRange["0", "9"],
      CharacterRange["A", "Z"],
      CharacterRange["a", "z"]
    ]
      // {Identity, Italicise, Embolden}
      // Through,
    Join[
      CharacterRange["\[CapitalAlpha]", "\[CapitalOmega]"],
      CharacterRange["\[Alpha]", "\[Omega]"]
    ],
    "e" ^ ("i" "\[Pi]") + 1 == 0 // TraditionalForm,
    {"\[PartialD]", "\[Epsilon]", "\[CurlyTheta]", "\[CurlyKappa]", "\[Phi]", "\[CurlyRho]", "\[CurlyPi]", Infinity}
  } // Style[#, 24] &;
  {
    expressions,
    expressions // Ex @ FString[
      "test-latex-{ToLowerCase[$SystemID]}-v{Floor[$VersionNumber]}.pdf"
    ]
  }
]


(* Tests for tracing gamma *)
(
  {
    (* U+0303 COMBINING TILDE \tilde *)
    LaTeXStyle /@ {"\[Gamma]" , "\:0303"} // SeparatedRow["\[NegativeMediumSpace]\[NegativeMediumSpace]"] @@ # &,
    (* U+2022 BULLET \bullet *)
    (* (but bullet is too small) *)
    LaTeXStyle @ Subscript["\[Gamma]", "\[Bullet]"],
    LaTeXStyle @ Subscript["\[Gamma]", Style["\[NegativeVeryThinSpace]\[Bullet]", Magnification -> 1.1]],
    (* U+22C6 STAR OPERATOR \star becomes an asterisk in Mathematica *)
    (* U+2605 BLACK STAR made smaller to approximate \star *)
    (* (but this approach does not scale with font size) *)
    LaTeXStyle @ Subscript["\[Gamma]", Style["\[NegativeVeryThinSpace]\[NegativeThinSpace]\[FivePointedStar]", Magnification -> 0.6]],
    (* Dummy for trailing commas *)
    Nothing
  }
    // SeparatedRow[" "] @@ # &
    // Style[#, 24] &
    // Ex["gamma-tracing.pdf"]
)


3


Module[{fun, x},
  fun = Cos[#] - # &;
  x = SeekRoot[fun, {0, 1}];
  {x, fun[x]}
]


Module[{fun, x},
  fun = Cos[#] - # &;
  x = SeekRoot[fun, 0.7];
  {x, fun[x]}
]


SeekRootBisection[Cos[#] - # &, {0, 1}]
SeekRootBisection[Cos[#] - # &, {0, 1}, "ReturnIterations" -> True]


Module[{fun, xMin, xMax, xAny, xFirst},
  fun = Function[{x}, (x - 1/2) Sinc[x]];
  (*
    The roots are 1/2, \[Pi], 2\[Pi], ....
    SeekRootBisection does not always find the first root.
    To that end, use SeekFirstRootBisection.
  *)
  xMin = 0;
  xMax = 5/2 Pi;
  xAny = SeekRootBisection[fun, {xMin, xMax}];
  xFirst = SeekFirstRootBisection[fun, {xMin, xMax}];
  Plot[fun[x], {x, xMin, xMax},
    Epilog -> {
      PointSize[Large],
      Red, Point @ {xAny, 0}, Text["Any", {xAny, 0}, {-1, 2}],
      Blue, Point @ {xFirst, 0}, Text["First", {xFirst, 0}, {-1, 2}]
    }
  ]
]


UniformRange[0, 1, 0.099]
UniformRange[0, 2 Pi, 5 Degree]


Table[
  {x, ToName[x]}
, {x, {0, 0., 1, Pi, 1.*^-10, 3 + 0. I, 3 + 4. I, Infinity}}] // TableForm


Way[0, 1]


Way["x", "y", "p"]


Module[{points},
  points = RandomPoint[Circle[], 1000];
  points = points // DeleteNearbyPoints[0.2] // SortByPhi;
  ListPlot[points,
    AspectRatio -> Automatic,
    LabelingFunction -> Function @ Placed[#2[[2]], Center],
    PlotRange -> All,
    PlotStyle -> Yellow
  ]
]


(* Perspective test *)
Module[
  {
    points, cube3d,
    eye, canvasNormal, canvasPoint, vertical,
    cube3dprojected, cube2d,
    dummyForTrailingCommas
  },
  points = Tuples[{0, 1}, 3] - 1/2;
  cube3d = {
    Red, Line @ points[[{2, 4, 8, 6, 2}]],
    Black,
      Line @ points[[{1, 2}]],
      Line @ points[[{3, 4}]],
      Line @ points[[{5, 6}]],
      Line @ points[[{7, 8}]],
    Blue, Line @ points[[{1, 3, 7, 5, 1}]],
    Nothing
  };
  eye = 2 {1, 2, 3};
  canvasNormal = {1, 2, 3};
  canvasPoint = {1, 2, 3};
  vertical = {0, 0, 1};
  cube3dprojected = cube3d /. {
    {x_?NumericQ, y_?NumericQ, z_?NumericQ} :>
      OnePointPerspective[3][eye, canvasNormal, canvasPoint, vertical]
        @ {x, y, z}
  };
  cube2d = cube3d /. {
    {x_?NumericQ, y_?NumericQ, z_?NumericQ} :>
      OnePointPerspective[2][eye, canvasNormal, canvasPoint, vertical]
        @ {x, y, z}
  };
  {
    Show[
      Graphics3D @ {
        cube3d,
        Opacity[0.5], Yellow, Hyperplane[canvasNormal, canvasPoint],
        cube3dprojected,
        PointSize[Large], Black, Point @ {eye},
        {}
      },
      {}
      , Boxed -> False
      , ViewPoint -> {-1.6, 2, 3}
    ],
    Graphics[cube2d]
  }
]


(* ::Section:: *)
(*Testing Curvilinear.wl*)


(* NOTE: For some reason \[FormalPhi] is not protected in Version 11.0. *)
With[{r = \[FormalR], phi = \[FormalPhi], x = \[FormalX], y = \[FormalY]},
  {
    XPolar[r, phi],
    YPolar[r, phi],
    XYPolar[r, phi],
    RPolar[x, y],
    PhiPolar[x, y],
    RPhiPolar[x, y],
    HRPolar[r, phi],
    HPhiPolar[r, phi],
    ARPolar[r, phi],
    APhiPolar[r, phi],
    APhiPolar[r, phi] == Cross @ ARPolar[r, phi]
  } // TableForm[#, TableDepth -> 1] &
]


With[{u = \[FormalU], v = \[FormalV], x = \[FormalX], y = \[FormalY]},
  {
    XBipolar[u, v],
    YBipolar[u, v],
    XYBipolar[u, v],
    UBipolar[x, y],
    VBipolar[x, y],
    UVBipolar[x, y],
    HBipolar[u, v],
    HUBipolar[u, v],
    HVBipolar[u, v],
    AUBipolar[u, v],
    AVBipolar[u, v],
    AVBipolar[u, v] == Cross @ AUBipolar[u, v] // FullSimplify
  } // TableForm[#, TableDepth -> 1] &
]


(* ::Section:: *)
(*Testing greys for printing*)


Module[
  {
    fun1, fun2,
    pValues,
    greyList, blackList,
    intensityValues, thicknessValues,
    dummyForTrailingCommas
  },
  (* Functions to be plotted *)
  fun1[p_][x_] := (1 - p) - x;
  fun2[p_][x_] := x^p;
  pValues = Subdivide[4];
  greyList[x_] := Table[{fun1[p][x], fun2[p][x]}, {p, pValues}];
  blackList[x_] := Table[fun2[p][x], {p, pValues}];
  (* List of plots *)
  intensityValues = Range[0.3, 0.8, 0.1];
  thicknessValues = {Automatic, Thick};
  Table[
    Show[
      (* Background greys *)
      Plot[greyList[x], {x, 0, 1}
        , AxesLabel -> {None, None}
        , AspectRatio -> 1/2
        , ImageSize -> 360
        , PlotLabel -> StringTemplate["GrayLevel[``], ``"][intensity, thickness]
        , PlotRange -> {-0.2, 1}
        , PlotStyle -> GrayLevel[intensity]
        , PlotOptions[Axes] // Evaluate
      ],
      (* Main blacks *)
      Plot[blackList[x], {x, 0.2, 0.4}
        , PlotStyle -> Directive[Black, thickness]
        , PlotOptions[Axes] // Evaluate
      ],
      {}
    ]
    , {intensity, intensityValues}
    , {thickness, thicknessValues}
  ]
    // TableForm[#, TableHeadings -> {intensityValues, thicknessValues}] &
] // Ex["test-grey-figure.pdf"]
