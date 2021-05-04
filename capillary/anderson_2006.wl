(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(*
  This be an independent implementation of
    Anderson M. L., Bassom A. P., & Fowkes N. (2006).
    Exact solutions of the Laplace--Young equation.
    Proc. R. Soc. Lond. Ser. A Math. Phys. Eng. Sci., 462 (2076), 3645--3656.
    <https://doi.org/10.1098/rspa.2006.1744>
  for the purposes of independently reproducing its Figure 3
  (traced boundaries with various indentations for the half-plane solution).
*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "anderson_2006"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Constants for gamma in Manipulate*)


gammaInitial = 30 Degree;
gammaClearance = 10^-3;


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*Half-plane wall height (2.4)*)


etaWall[gamma_] := Sqrt[2 (1 - Sin @ gamma)];


(* ::Subsubsection:: *)
(*Implicit known half-plane solution (2.6)*)


yKnown[eta_] := ArcCosh[2 / eta] - Sqrt[4 - eta^2] - ArcCosh @ Sqrt[2] + Sqrt[2];


(* ::Subsubsection:: *)
(*Solution derivative \[PartialD]\[Eta]/\[PartialD]y*)


etaYDerivative[eta_] := 1 / yKnown'[eta];


(* ::Subsubsubsection:: *)
(*Verify solution derivative satisfies (2.3) with A = 1*)


With[{eta = \[FormalEta], a = 1},
  (1 + etaYDerivative[eta]^2) ^ (-1/2) == a - 1/2 eta^2
    // FullSimplify[#, 0 < eta < Sqrt[2]] &
]


(* ::Subsection:: *)
(*Boundary tracing*)


(* ::Subsubsection:: *)
(*Boundary tracing system of ODEs (3.2)*)


etaWallPlus[gamma_] := Sqrt[2 (1 + Sin @ gamma)];


yDerivative[eta_] :=
  Divide[
    eta^2 - 2,
    eta Sqrt[4 - eta^2]
  ];


xDerivative[gamma_][eta_] :=
  Divide[
    2 yDerivative[eta] Cos[gamma],
    Sqrt[(eta^2 - etaWall[gamma]^2) (etaWallPlus[gamma]^2 - eta^2)]
  ];


(* ::Subsubsubsection:: *)
(*Verify tracing system satisfies boundary condition (3.1)*)


(* Mathematica needs a bit of help *)
assumedNegativeQuantity[gamma_][eta_] := 2 - 4 eta^2 + eta^4 + 2 Cos[2 gamma];


With[{eta = \[FormalEta], gamma = \[FormalGamma]},
  FullSimplify[
    Divide[
      xDerivative[gamma][eta] etaYDerivative[eta],
      Sqrt[xDerivative[gamma][eta]^2 + yDerivative[eta]^2]
    ]
      ==
    Cos[gamma] Sqrt[1 + etaYDerivative[eta]^2]
    , And[
        etaWall[gamma] < eta < Sqrt[2],
        0 < gamma < Pi/2,
        assumedNegativeQuantity[gamma][eta] < 0,
        True
      ]
  ]
]


(* ::Subsubsubsection:: *)
(*Verify assumedPositiveQuantity is actually positive*)


Manipulate[
  Plot[
    assumedNegativeQuantity[gamma][eta]
    , {eta, etaWall[gamma], Sqrt[2]}
  ]
  , {{gamma, gammaInitial}, gammaClearance, Pi/2 - gammaClearance}
]


(* ::Subsubsection:: *)
(*Elliptic integral arguments (unnumbered equation after (3.4))*)


curlyPhi[gamma_][eta_] :=
  ArcSin @ Sqrt @ Divide[
    2 + 2 Sin[gamma] - eta^2,
    4 Sin[gamma]
  ];


muPlus[gamma_] := 2 Sin[gamma] / (1 + Sin[gamma]);
muMinus[gamma_] := 2 Sin[gamma] / (1 - Sin[gamma]);


(* ::Subsubsection:: *)
(*Traced boundaries (3.3) & (3.4)*)


(*
  **IMPORTANT**
  1. There is a typo in (3.3). The very last mu_plus should be a mu_minus.
  2. F and Pi are defined as per Mathematica's EllipticF and EllipticPi,
     with the last parameter being m == k^2,
     rather than k used in Gradsteyn & Ryzhik (1980).
*)
xTraced[gamma_][eta_] := (
  Sqrt[2] Cos[gamma] / Sqrt[1 - Sin[gamma]]
  Subtract[
    EllipticF[curlyPhi[gamma][eta], -muMinus[gamma]],
    1 / (1 + Sin[gamma]) EllipticPi[muPlus[gamma], curlyPhi[gamma][eta], -muMinus[gamma]]
  ]
);


(* Same as known solution. *)
yTraced[eta_] := yKnown[eta];


(* ::Subsubsubsection:: *)
(*Verify traced boundaries*)


(* Note the introduced negative sign for the x-component. *)
With[{eta = \[FormalEta], gamma = \[FormalGamma]},
  FullSimplify[
    {
      -xTraced[gamma]'[eta] == xDerivative[gamma][eta],
      yTraced'[eta] == yDerivative[eta]
    }
    , And[
        etaWall[gamma] < eta < Sqrt[2],
        0 < gamma < Pi/2,
        assumedNegativeQuantity[gamma][eta] < 0,
        True
      ]
  ]
]


(* ::Subsubsubsection:: *)
(*Visualise traced boundaries*)


Manipulate[
  ParametricPlot[
    {
      {xTraced[gamma][eta], yTraced[eta]},
      {-xTraced[gamma][eta], yTraced[eta]}
    }
    , {eta, etaWall[gamma], Sqrt[2]}
    , Epilog -> {
        Dashed,
        Line @ {
          {-2, yTraced @ etaWall[gamma]},
          {+2, yTraced @ etaWall[gamma]}
        },
        {}
      }
    , PlotRange -> All
  ]
  , {{gamma, gammaInitial}, gammaClearance, Pi/2 - gammaClearance}
]


(* ::Subsubsection:: *)
(*Canonically-translated copy*)


(*
  Defined so that:
  1. The corner lies at the origin; and
  2. The boundary extends to the positive-x half.
*)
xCanonical[gamma_][eta_] := xTraced[gamma][Sqrt[2]] - xTraced[gamma][eta];
yCanonical[eta_] := yTraced[eta];


(* ::Subsubsubsection:: *)
(*Visualise traced boundaries*)


Manipulate[
  ParametricPlot[
    {
      {xCanonical[gamma][eta], yCanonical[eta]},
      {-xCanonical[gamma][eta], yCanonical[eta]}
    }
    , {eta, etaWall[gamma], Sqrt[2]}
    , Epilog -> {
        Dashed,
        Line @ {
          {-2, yTraced @ etaWall[gamma]},
          {+2, yTraced @ etaWall[gamma]}
        },
        {}
      }
    , PlotRange -> All
  ]
  , {{gamma, gammaInitial}, gammaClearance, Pi/2 - gammaClearance}
]


(* ::Section:: *)
(*Reproducing Figure 3*)


(*
  1. From the caption of Figure 2, we guess that Figure 3
     also has gamma == pi/4.
  2. We guess that the topmost plot of Figure 3
     shows the canonically translated copy.
  3. We guess that all plots have the full vertical range 0 < y < y_wall
     corresponding to eta_wall < eta < sqrt(2).
*)


(*
  4. We note that the slope dy/dx at y == 0 (the corner),
     equal to y-derivative / x-derivative at eta == sqrt(2)),
     is equal to tan(gamma):
*)
With[{eta = \[FormalEta], gamma = \[FormalGamma]},
  Limit[yDerivative[eta] / xDerivative[gamma][eta], eta -> Sqrt[2]]
    // FullSimplify[#, 0 < gamma < Pi/2] &
]


(*
  5. The top plot of Figure 3 appears to have slope 1 at the corner,
     which is indeed consistent with gamma == pi/4.
  6. Measurements indicate that all indentations in Figure 3 are full-depth,
     i.e. extend all the way down to y == 0 (equivalently eta == sqrt(2)).
*)


Module[
  {
    gamma,
    h, sqrt2,
    xC, yC,
    reflectHoriztonal,
    tracedOptions,
    plotSpacing,
    maxDepth, maxHalfWidth,
    xMax, yMax, yMin,
    commonPlot,
    plot,
    dummyForTrailingCommas
  },
  (* Constants *)
  gamma = Pi/4;
  (* Abbreviations *)
  h = etaWall[gamma];
  sqrt2 = Sqrt[2];
  xC = xCanonical[gamma];
  yC = yCanonical;
  reflectHoriztonal[{x_, y_}] := {{x, y}, {-x, y}};
  tracedOptions = {
    PlotPoints -> 2,
    PlotRange -> Full,
    PlotStyle -> BoundaryTracingStyle["Traced"],
    Nothing
  };
  plotSpacing = 32;
  (* Maximum depth and half-width of indentations *)
  maxDepth = yC[h];
  maxHalfWidth = xC[h];
  (* Plot range *)
  xMax = 2.1;
  yMax = +1.1 maxDepth;
  yMin = -0.1 maxDepth;
  (* Common plot (terminal curve) *)
  commonPlot =
    Show[
      EmptyFrame[{-xMax, xMax}, {yMin, yMax}
        , Frame -> None
        , ImageSize -> 240
      ],
      Graphics @ {
        BoundaryTracingStyle["Contour"],
        Line @ {{-xMax, maxDepth}, {+xMax, maxDepth}}
      },
      {}
    ];
  (* Plot 1: V *)
  plot[1] =
    Show[
      commonPlot,
      (* V portion *)
      ParametricPlot[
        {xC[eta], yC[eta]} // reflectHoriztonal
        , {eta, sqrt2, h}
        , tracedOptions // Evaluate
      ],
      (* Terminal curve portion *)
      ParametricPlot[
        {x, yC[h]} // reflectHoriztonal
        , {x, xC[h], xMax}
        , tracedOptions // Evaluate
      ],
      {}
    ];
  (* Plot 2: W *)
  plot[2] =
    Module[{xCorner, etaIntersection},
      xCorner = 0.15;
      etaIntersection = SeekRoot[xC[#] - xCorner &, {h, sqrt2}, 5];
      Show[
        commonPlot,
        (* W outer portions *)
        ParametricPlot[
          {xCorner + xC[eta], yC[eta]} // reflectHoriztonal
          , {eta, sqrt2, h}
          , tracedOptions // Evaluate
        ],
        (* W inner portion *)
        ParametricPlot[
          {xCorner - xC[eta], yC[eta]} // reflectHoriztonal
          , {eta, sqrt2, etaIntersection}
          , tracedOptions // Evaluate
        ],
        (* Terminal curve portion *)
        ParametricPlot[
          {x, yC[h]} // reflectHoriztonal
          , {x, xC[h], xMax}
          , tracedOptions // Evaluate
        ],
        {}
      ]
    ];
  (* Plot 3: four bumps *)
  plot[3] =
    Module[{xBump, etaBump},
      xBump = xMax / 4;
      etaBump = SeekRoot[xC[#] - xBump &, {h, sqrt2}, 5];
      Show[
        commonPlot,
        (* Bumps *)
        ParametricPlot[
          Table[
            {
              {xC[eta] + xCopy, yC[eta]} // reflectHoriztonal,
              {2 xBump - xC[eta] + xCopy, yC[eta]} // reflectHoriztonal
            }
            , {xCopy, {0, 2 xBump}}
          ]
          , {eta, sqrt2, etaBump}
          , tracedOptions // Evaluate
        ],
        {}
      ]
    ];
  (* Plot 4: three Ws and two half-Ws *)
  plot[4] =
    Module[{xCorner, etaIntersection, xBump, xBumpDiff, etaBumpDiff},
      xCorner = 0.15;
      etaIntersection = SeekRoot[xC[#] - xCorner &, {h, sqrt2}, 5];
      xBump = xMax / 4 ;
      xBumpDiff = xBump - xCorner;
      etaBumpDiff = SeekRoot[xC[#] - xBumpDiff &, {h, sqrt2}, 5];
      Show[
        commonPlot,
        (* W outer portions *)
        ParametricPlot[
          Table[
            {
              {xC[eta] + xCorner + xCopy, yC[eta]} // reflectHoriztonal,
              {2 xBump - xCorner - xC[eta] + xCopy, yC[eta]} // reflectHoriztonal
            }
            , {xCopy, {0, 2 xBump}}
          ]
          , {eta, sqrt2, etaBumpDiff}
          , tracedOptions // Evaluate
        ],
        (* W inner portion *)
        ParametricPlot[
          Table[
            {
              {xCorner - xC[eta] + xCopy, yC[eta]} // reflectHoriztonal,
              {2 xBump - xCorner + xC[eta] + xCopy, yC[eta]} // reflectHoriztonal
            }
            , {xCopy, {0, 2 xBump}}
          ]
          , {eta, sqrt2, etaIntersection}
          , tracedOptions // Evaluate
        ],
        {}
      ]
    ];
  (* Plot 5: Lots of upside-down Vs *)
  plot[5] =
    Module[{xCorner, etaIntersection},
      xCorner = 0.15;
      etaIntersection = SeekRoot[xC[#] - xCorner &, {h, sqrt2}, 5];
      Show[
        commonPlot,
        ParametricPlot[
          Table[
            {
              {xCorner - xC[eta] + xCopy, yC[eta]} // reflectHoriztonal,
              {xCorner + xC[eta] + xCopy, yC[eta]} // reflectHoriztonal
            }
            , {xCopy, 0, xMax - xCorner, 2 xCorner}
          ]
          , {eta, sqrt2, etaIntersection}
          , tracedOptions // Evaluate
        ],
        {}
      ]
    ];
  plot[All] =
    GraphicsColumn[
      {
        plot[1],
        GraphicsRow[plot /@ {2, 3}, Spacings -> plotSpacing],
        GraphicsRow[plot /@ {4, 5}, Spacings -> plotSpacing],
        Nothing
      }
      , Alignment -> Center
      , ImageSize -> 0.85 ImageSizeTextWidth
      , Spacings -> 0.05 ImageSizeTextWidth
    ];
  plot[All]
] // Ex["half_plane-traced-boundaries.pdf"]
