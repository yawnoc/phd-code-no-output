(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ NotebookDirectory[];


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Section:: *)
(*Figure: Terminal points (terminal-*.pdf)*)


(* ::Subsection:: *)
(*Main (sub)figures*)


(*
  In all cases we put \[CapitalPhi] == y + 3/2 x^2,
  with a terminal point (x, y) == (0, 0).
  The terminal curve is y == -3/2 x^2.
  The non-viable domain is y < -3/2 x^2.
  We choose three known solutions
  with local T-contour T == 0
  to achieve the three (non-degenerate) cases
  for critical terminal points:
    "ordinary"
      T == (x - 1)^2 + (y + 1)^2 - 2
      T == 0 (y ~ x) crosses the terminal curve at a nonzero-angle.
    "critical_hyperbolic"
      T == y + 1/2 x^2
      T == 0 (y ~ 1/2 x^2) is tangential on the viable side.
    "critical_elliptic"
      T == y + 4 x^2
      T == 0 (y ~ -4 x^2) is tangential on the non-viable side.
  Since \[CapitalPhi] == (grad T)^2 - F^2,
  we have F == sqrt((grad T)^2 - \[CapitalPhi]).
 *)


Module[
 {phi,
  caseList, tAss,
  imageSize, xMax, yMax,
  xMaxLess, yMaxLess,
  xMaxMore, yMaxMore,
  vi, t, p, q, grad2, f,
  xyTraSystem,
  sMax, sLow, sHigh, xyTra
 },
  (* Viability function \[CapitalPhi] *)
  phi = Function[{x, y}, y + 3/2 x^2];
  (* List of cases to be plotted *)
  caseList = {
    "ordinary",
    "critical_hyperbolic",
    "critical_elliptic"
  };
  (* Known solutions for the cases *)
  tAss = AssociationThread[
    caseList -> {
      Function[{x, y}, (x - 1)^2 + (y + 1)^2 - 2],
      Function[{x, y}, y + 1/2 x^2],
      Function[{x, y}, y + 4 x^2]
    }
  ];
  (* Plotting constants *)
  imageSize = 0.325 ImageSizeTextWidth;
  xMax = 0.65;
  yMax = 0.35;
  xMaxLess = 0.95 xMax;
  yMaxLess = 0.95 yMax;
  xMaxMore = 1.2 xMax;
  yMaxMore = 1.2 yMax;
  (* For the three cases *)
  Table[
    Module[{x, y, s},
      (* Viability function \[CapitalPhi] *)
      vi = phi[x, y];
      (* Known solution T *)
      t = tAss[case][x, y];
      (* Derivatives of T *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient *)
      grad2 = p^2 + q^2;
      (* Flux function *)
      f = Sqrt[grad2 - vi];
      (* System of tracing equations *)
      xyTraSystem[upperBranch_] :=
        Module[{sign, xDer, yDer},
          sign = If[upperBranch, +1, -1];
          (* Return system of ODEs *)
          xDer = (-q f + sign p Re @ Sqrt[vi]) / grad2;
          yDer = (+p f + sign q Re @ Sqrt[vi]) / grad2;
          {x' == xDer, y' == yDer } /. {
            x' -> x'[s],
            y' -> y'[s],
            x -> x[s],
            y -> y[s]
          }
        ];
      (* Solve for traced boundaries *)
      If[case != "critical_elliptic",
        (* Not elliptic case: traced boundaries exist *)
        sMax = 1;
        sLow = -sMax;
        sHigh = If[case === "ordinary", 0, sMax];
        xyTra = Table[
          NDSolveValue[
            {
              xyTraSystem[upperBranch],
              x[0] == 0, y[0] == 0
            }, {x, y}, {s, sLow, sHigh},
            NoExtrapolation
          ]
        , {upperBranch, {True, False}}],
        (* Elliptic case: traced boundaries do not exist *)
        xyTra = {}
      ];
      (* Plot *)
      Show[
        EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
          Frame -> None,
          ImageSize -> imageSize
        ],
        (* Non-viable domain *)
        RegionPlot[vi < 0,
          {x, -xMaxMore, xMaxMore}, {y, -yMaxMore, 0},
          BoundaryStyle -> None,
          PlotPoints -> 4,
          PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        (* Terminal curve *)
        ContourPlot[vi == 0,
          {x, -xMaxMore, xMaxMore}, {y, -yMaxMore, 0},
          ContourStyle -> BoundaryTracingStyle["Terminal"],
          PlotPoints -> 8
        ],
        (* Local T-contour *)
        ContourPlot[t == 0,
          {x, -xMaxLess, xMaxLess}, {y, -yMaxLess, yMaxLess},
          ContourStyle -> BoundaryTracingStyle["ContourPlain"],
          PlotPoints -> 10
        ],
        (* Traced boundaries *)
        Table[
          ParametricPlot[xy[s] // Through,
            {s, 1/2 DomainStart[xy], 1/2 DomainEnd[xy]},
            PlotPoints -> 2,
            PlotStyle -> BoundaryTracingStyle["Traced"]
          ]
        , {xy, xyTra}],
        (* Terminal point *)
        Graphics @ {
          GeneralStyle["Point"],
          GeneralStyle["Translucent"],
          Point @ {0, 0}
        }
      ]
    ] // Ex @ FString["terminal-{case}.pdf"]
  , {case, caseList}]
]


(* ::Subsection:: *)
(*Legend*)


Module[
  {
    curves, regions,
    dummyForTrailingCommas
  },
  curves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Terminal", "ContourPlain", "Traced"},
      {"terminal curve", Row @ {Italicise["T"], "\[Hyphen]contour"}, "traced boundary"}
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
    ];
  regions = (
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable", "Viable"},
      {"non\[Hyphen]viable domain", "viable domain"}
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
    ]
  );
  (* Combine *)
  GraphicsGrid[{curves, regions // Insert[Null, 2]}
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.11
    , Spacings -> {{0, 0.19, 0.05} ImageSizeTextWidth, 0}
  ]
] // Ex["terminal-legend.pdf"]


(* ::Section:: *)
(*Figure: Curvilinear coordinates (orthogonal-curvilinear-coordinates.pdf)*)


Module[
 {f, fInv, u, v, x, y,
  hu, hv,
  au, av,
  x0, y0, u0, v0,
  du, dv,
  xLeft, xRight, yTop, yBottom,
  originalStyle, displacedStyle,
  orthogonalityMarkerLength, orthogonalityMarkerStyle,
  basisVectorLength, basisVectorStyle,
  textStyle,
  auTipXY, avTipXY,
  huduUV, huduXY,
  hvdvUV, hvdvXY,
  uConstUV, uConstXY,
  vConstUV, vConstXY,
  constText, vectorText, displacedLengthText
 },
  (*
    Any orthogonal curvilinear coordinate system (u, v)
    can be written u + i v == f(x + i y)
    for some analytic function f.
    We take f == exp, f^(-1) == log
    so we have "reversed" polar coordinates.
  *)
  (* Transformation *)
  f = Exp;
  fInv = InverseFunction[f];
  u[x_, y_] := Re @ f[x + I y] // ComplexExpand // Evaluate;
  v[x_, y_] := Im @ f[x + I y] // ComplexExpand // Evaluate;
  x[u_, v_] := Re @ fInv[u + I v] // ComplexExpand // Evaluate;
  y[u_, v_] := (
    Im @ fInv[u + I v]
    // ComplexExpand
    // # /. {Arg[u + I v] -> ArcTan[u, v]} & (* ComplexExpand is dumb *)
    // Evaluate
  );
  (* Local basis *)
  hu[u_, v_] := D[{x, y} @@ {u, v} // Through, u] // Evaluate;
  hv[u_, v_] := D[{x, y} @@ {u, v} // Through, v] // Evaluate;
  (* Local orthonormal basis *)
  au[u_, v_] := Normalize @ hu[u, v];
  av[u_, v_] := Normalize @ hv[u, v];
  (* Original coordinates *)
  {x0, y0} = {0.35, 0.55};
  {u0, v0} = {u, v} @@ {x0, y0} // Through;
  (* Displacements *)
  du = 0.3;
  dv = 0.25;
  (* Plot range *)
  xLeft = 0.25;
  xRight = 0.35;
  yTop = 0.3;
  yBottom = 0.3;
  (* Styles *)
  originalStyle = Black;
  displacedStyle = BoundaryTracingStyle["ContourPlain"];
  orthogonalityMarkerLength = 0.02;
  orthogonalityMarkerStyle = Directive[
    EdgeForm[Black],
    FaceForm[None]
  ];
  basisVectorLength = 0.15;
  basisVectorStyle = Directive[AbsoluteThickness[2.75], Arrowheads[0.08]];
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Coordinates of tips of vectors *)
  auTipXY = {x0, y0} + basisVectorLength au[u0, v0];
  avTipXY = {x0, y0} + basisVectorLength av[u0, v0];
  (* Coordinates of displaced length labels *)
  huduUV = {u0 + du / 2, v0};
  huduXY = {x, y} @@ huduUV // Through;
  hvdvUV = {u0, v0 + dv / 2};
  hvdvXY = {x, y} @@ hvdvUV // Through;
  (* Coordinates of contour labels *)
  uConstUV = {u0, v0 - dv};
  uConstXY = {x, y} @@ uConstUV // Through;
  vConstUV = {u0 - 2/3 du, v0};
  vConstXY = {x, y} @@ vConstUV // Through;
  (* Text functions *)
  constText[coord_] := Italicise[coord] == "const";
  vectorText[coord_] :=
    Subscript[
      Embolden["a"],
      Row @ {"\[NegativeVeryThinSpace]", Spacer[0.2], Italicise[coord]}
    ];
  displacedLengthText[coord_] :=
    Row @ {
      Subscript[
        Italicise["h"],
        Row @ {"\[NegativeThinSpace]", Spacer[0.8], Italicise[coord]}
      ],
      "\[ThinSpace]",
      Row @ {
        "d", "\[NegativeThinSpace]",
        Spacer[0.85], Italicise[coord]
      }
    };
  (* Plot *)
  Show[
    (* Contours *)
    ContourPlot[
      {
        u[x, y] == u0 + du,
        v[x, y] == v0 + dv,
        u[x, y] == u0,
        v[x, y] == v0,
        Nothing
      },
      {x, x0 - xLeft, x0 + xRight},
      {y, y0 - yBottom, y0 + yTop},
      AspectRatio -> Automatic,
      ContourLabels -> None,
      ContourStyle -> {displacedStyle, displacedStyle, originalStyle, originalStyle},
      Frame -> None,
      ImageSize -> 0.4 ImageSizeTextWidth,
      PlotPoints -> 10
    ],
    Graphics @ {
      (* u == const *)
      Text[
        constText["u"] // textStyle,
        uConstXY,
        {0, -1.2},
        av @@ uConstUV
      ],
      (* v == const *)
      Text[
        constText["v"] // textStyle,
        vConstXY,
        {0, 0.7},
        au @@ vConstUV
      ]
    },
    (* Orthogonality marker *)
    Graphics @ {orthogonalityMarkerStyle,
      Rotate[
        Rectangle[
          {x0, y0},
          {x0, y0} + orthogonalityMarkerLength
        ],
        ArcTan @@ au[u0, v0],
        {x0, y0}
      ]
    },
    (* Basis vectors *)
    Graphics @ {basisVectorStyle,
      (* a_u *)
      Arrow @ {{x0, y0}, auTipXY},
      Text[
        vectorText["u"] // textStyle,
        auTipXY,
        {-0.3, -1.6}
      ],
      (* a_v *)
      Arrow @ {{x0, y0}, avTipXY},
      Text[
        vectorText["v"] // textStyle,
        avTipXY,
        {-0.9, 1.6}
      ]
    },
    (* Displaced length labels *)
    Graphics @ {
      (* h_u du *)
      Text[
        displacedLengthText["u"] // textStyle,
        huduXY,
        {0, 1.8},
        au @@ huduUV
      ],
      (* h_v du *)
      Text[
        displacedLengthText["v"] // textStyle,
        hvdvXY,
        {0, -2.2},
        av @@ hvdvUV
      ]
    }
  ]
] // Ex @ "orthogonal-curvilinear-coordinates.pdf"


(* ::Section:: *)
(*Figure: Cartesian normal vector (cartesian-normal)*)


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    dx, dy, ds, slope,
    positionStart, positionEnd,
    normalVector, normalVectorLength, normalVectorBase, normalVectorTip,
    textStyle,
    dummyForTrailingCommas
  },
  (* Differential lengths *)
  dx = 7;
  dy = 2;
  ds = Sqrt[dx^2 + dy^2];
  slope = dy/dx;
  (* Differential positions *)
  positionStart = {-dx, 0};
  positionEnd = {0, dy};
  (* Normal vector *)
  normalVector = {-dy, dx} / ds;
  normalVectorLength = 1.3 ds;
  normalVectorBase = Way[positionStart, positionEnd, 0.55];
  normalVectorTip = normalVectorBase + normalVectorLength * normalVector;
  (* Make diagram *)
  textStyle = Style[#, 10] &;
  Show[
    (* Displacements *)
    Graphics @ {
      Line @ {positionStart, {0, 0}, positionEnd},
      Text[
        Row @ {"d", Italicise["x"]} // textStyle
        , 1/2.2 positionStart
        , {0, 0.9}
      ],
      Text[
        Row @ {"d", Italicise["y"]} // textStyle
        , 1/2 positionEnd
        , {-1.7, 0}
      ],
      {}
    },
    (* Boundary *)
    Graphics @ {
      BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"],
      Line @ {
        Way[positionStart, positionEnd, -0.4],
        Way[positionStart, positionEnd, +1.5],
        Nothing
      },
      {}
    },
    (* Normal *)
    Graphics @ {
      Arrowheads[0.2], Thickness[0.04],
      Arrow @ {normalVectorBase, normalVectorTip},
      Text[
        Embolden["n"] // textStyle
        , normalVectorTip
        , {-3, 0}
      ],
      {}
    },
    {}
    , ImageSize -> 0.2 ImageSizeTextWidthBeamer
  ]
] // Ex["cartesian-normal-slides.pdf"];
