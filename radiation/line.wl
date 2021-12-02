(* ::Package:: *)

(* ::Text:: *)
(*See r2 (manuscripts/radiation-2-line.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "line"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*\[Psi] = \[Psi](r)*)


(* ::Text:: *)
(*See (r2.16).*)


psi[r_] := r Log[r]^4;


(* ::Subsection:: *)
(*r_\[Natural] (transition terminal radius)*)


(* ::Text:: *)
(*See (r2.19).*)


rNat = Exp[-4]


(* ::Subsection:: *)
(*A_\[Natural] (transition dimensionless group)*)


(* ::Text:: *)
(*See (r2.20).*)


aNat = psi[rNat]


(* ::Subsection:: *)
(*r_\[Flat] (inner terminal radius)*)


(* InverseFunction can't invert psi ... *)(*
InverseFunction @ Function[{r},
  ConditionalExpression[
    psi[r],
    0 < r < rNat
  ] // Evaluate
]*)


(* ... so spot pattern from special values in Solve: *)(*
With[{r = \[FormalR]},
  Table[
    Select[
      r /. Solve[psi[r] == a, r, Reals],
      0 < # < rNat &
    ]
  , {a, {1, E, 3, Pi}}]
]*)


rFlat[a_] := Exp[4 ProductLog[-1, -a^(1/4) / 4]];


(* Check this: *)(*
With[{a = \[FormalCapitalA]},
  Block[{$Assumptions = 0 < a < aNat},
    psi @ rFlat[a] == a // FullSimplify
  ]
]*)


(* ::Subsection:: *)
(*r_\[Sharp] (outer terminal radius)*)


(* InverseFunction can't invert psi ... *)(*
InverseFunction @ Function[{r},
  ConditionalExpression[
    psi[r],
    rNat < r < 1
  ] // Evaluate
]*)


(* ... so spot pattern from special values in Solve: *)(*
With[{r = \[FormalR]},
  Table[
    Select[
      r /. Solve[psi[r] == a, r, Reals],
      rNat < # < 1 &
    ]
  , {a, {1, E, 3, Pi}}]
]*)


rSharp[a_] := Exp[4 ProductLog[-a^(1/4) / 4]];


(* Check this: *)(*
With[{a = \[FormalCapitalA]},
  Block[{$Assumptions = 0 < a < aNat},
    psi @ rSharp[a] == a // FullSimplify
  ]
]*)


(* ::Subsection:: *)
(*r_i (inflection radius)*)


(* ::Text:: *)
(*See (r2.32).*)


rInfl = Exp[-2];


(* ::Subsection:: *)
(*\[CapitalPhi] (viability)*)


(* phi I am reserving for azimuthal angle *)
vi[a_] := Function[{r},
  (a^2 - psi[r]^2) / (a^2 r^2) // Evaluate
] // Evaluate;


(* ::Subsection:: *)
(*Representative values of A*)


aHot = 3;
aColdHot = aNat;
aCold = 5;
aHot < aColdHot < aCold


(* ::Subsection:: *)
(*Traced boundary \[Phi] = \[Phi](r)*)


phiTraDer[a_] := Function[{r},
  -psi[r] / (r Sqrt[a^2 - psi[r]^2]) // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]*)


phiTraHot["outer"] =
  With[{a = aHot, phi = \[FormalPhi]},
    Module[{rMin, rMax, rInit},
      rMin = rSharp[a];
      rMax = 1;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


phiTraHot["inner"] =
  With[{a = aHot, phi = \[FormalPhi]},
    Module[{rMin, rMax, rInit},
      rMin = 0;
      rMax = rFlat[a];
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*Cold-to-hot transition A = A_\[Natural]*)


phiTraColdHot["outer"] =
  With[{a = aColdHot, phi = \[FormalPhi]},
    Module[{rMin, rMax, rInit},
      rMin = rNat;
      rMax = 1;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


phiTraColdHot["inner"] =
  With[{a = aColdHot, phi = \[FormalPhi]},
    Module[{rMin, rMax, rInit},
      rMin = 0;
      rMax = rNat;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*Cold regime A > A_\[Natural]*)


phiTraCold =
  With[{a = aCold, phi = \[FormalPhi]},
    Module[{rMin, rMax, rInit},
      rMin = 0;
      rMax = 1;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsection:: *)
(*Numerical verification (finite elements)*)


(* ::Subsubsection:: *)
(*Generate finite element mesh*)


(* ::Text:: *)
(*This is the "generic" convex boundary plotted in*)
(*"line-hot-traced-convex-generic-with-known.png".*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[{dest},
  dest = "line-verification-mesh.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {a, rSh,
      nSpike, modNSpike, mod2Pi,
      rSpike, phiSpike, phiCentre, phiHalfWidth,
      rBath, tBath, phSpacing, sSpacing,
      phi, phValues,
      extPointList, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      (* Dimensionless group *)
      a = aHot;
      (* Critical terminal curve *)
      rSh = rSharp[a];
      (* Spikes *)
      nSpike = 2;
      modNSpike = Mod[#, nSpike, 1] &;
      mod2Pi = Mod[#, 2 Pi] &;
      (* Spike radii *)
      rSpike[1] = rInfl;
      rSpike[2] = Way[rSh, rInfl, 3/4];
      Table[
        phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
      , {n, nSpike}];
      (* Spike centres *)
      phiCentre[1] = 0;
      phiCentre[2] = 135 Degree;
      phiCentre[n_] := phiCentre[n // modNSpike];
      (* Spike angular half-widths *)
      Table[
        phiHalfWidth[n] =
          phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
      , {n, nSpike}];
      phiHalfWidth[n_] := phiHalfWidth[n // modNSpike];
      (* Heat bath *)
      rBath = rSh / 2;
      tBath = Log[1 / rBath];
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 5 degrees along the terminal curve) *)
      phSpacing = 5 Degree;
      sSpacing = rSh phSpacing;
      (* External (radiation) boundary points *)
      (*
        I am lazy, so I am simply listing the boundary points from
        "line-hot-traced-convex-generic-with-known.png",
        then removing duplicates and sorting them by azimuthal angle.
       *)
      extPointList = Join[
        (* Spikes *)
        Join @@ Table[
          phi[r_] := phiTraHot["outer"][r] - phiSpike[n];
          Join @@ Table[
            (* (here, # is r such that phi[r] equals ph) *)
            SeekRoot[phi[#] - ph &, {rSh, rSpike[n]}, 5]
            // {
              XYPolar[#, phiCentre[n] +  phi[#]],
              XYPolar[#, phiCentre[n] -  phi[#]]
            } &
          , {ph, UniformRange[0, phiHalfWidth[n], phSpacing]}]
        , {n, nSpike}],
        (* Terminal curve *)
        Join @@ Table[
          Table[
            XYPolar[rSh, ph]
          , {ph, UniformRange[
            phiCentre[n] + phiHalfWidth[n] // mod2Pi,
            phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi,
            phSpacing
          ]}]
        , {n, nSpike}]
      ];
      extPointList = DeleteDuplicates[
        extPointList,
        Norm[#2 - #1] < sSpacing / 4 &
      ];
      extPointList = SortBy[
        extPointList,
        mod2Pi[ArcTan @@ #] &
      ];
      (* Internal (heat bath) boundary *)
      intPointList = Table[
        XYPolar[rBath, ph]
      , {ph, UniformRange[0, 2 Pi, sSpacing / rBath]}] // Most;
      (* Numbering *)
      nExt = Length[extPointList];
      nInt = Length[intPointList];
      mod[n_] := Mod[#, n, 1] &;
      (* Build boundary element mesh *)
      bMesh = ToBoundaryMesh[
        "Coordinates" -> Join[extPointList, intPointList],
        "BoundaryElements" -> {
          (* External *)
          LineElement[
            Table[{n, n + 1}, {n, nExt}] // mod[nExt]
          ],
          (* Internal *)
          LineElement[
            nExt + (Table[{n, n + 1}, {n, nInt}] // mod[nInt])
          ]
        }
      ];
      (* Build mesh *)
      mesh = ToElementMesh[bMesh,
        "ImproveBoundaryPosition" -> True,
        "RegionHoles" -> {0, 0}
      ];
      (* Predicate functions for exterior and interior boundaries *)
      prExt = Function[{x, y},
        RPolar[x, y] > Way[rBath, rSh] // Evaluate
      ];
      prInt = Function[{x, y},
        RPolar[x, y] < Way[rBath, rSh] // Evaluate
      ];
      {a, tBath, mesh, prExt, prInt}
       // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Solve boundary value problem*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[{dest},
  dest = "line-verification-solution.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, tBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "line-verification-mesh.txt";
      {a, tBath, mesh, prExt, prInt} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* External (radiation) boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prExt[x, y]],
            (* Internal (heat bath) boundary condition *)
            DirichletCondition[t[x, y] == tBath, prInt[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      ];
      tSol // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsection:: *)
(*Italicise symbols*)


aIt = Italicise["A"];
rIt = Italicise["r"];
tIt = Italicise["T"];


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
sharpStyle = Magenta;
flatStyle = Darker[Green];


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;
inflStyle = Directive[Dashed, Green];


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


pointStyle = PointSize[Large];


transStyle = Function[{a}, If[a == aNat, Style[#, Red] &, Identity]];


(* ::Section:: *)
(*\[Psi]*)


(* ::Subsection:: *)
(*Properties*)


With[{r = \[FormalR]},
  {
    {"psi", psi[r]},
    {"psi(0)", Limit[psi[r], r -> 0]},
    {"psi(1)", psi[1]},
    {"psi'", psi'[r] // FullSimplify},
    {"psi'(0)", Limit[psi'[r], r -> 0]},
    {"psi'(rNat)", psi'[rNat]},
    {"psi'(1)", psi'[1]},
    {"psi''", psi''[r] // FullSimplify},
    {"psi''(rNat)", psi''[rNat]},
    {"psi''(rNat)", psi''[rNat] // N},
    {"rNat", rNat},
    {"rNat", rNat // N},
    {"ANat", aNat},
    {"ANat", aNat // N}
  } // PrettyString["Nat" -> "_\[Natural]", "psi" -> "\[Psi]"]
    // TableForm
] // Ex["line-psi-properties.pdf"]


(* ::Subsection:: *)
(*Plot*)


Show[
  (* psi *)
  Plot[psi[r], {r, 0, 1},
    AxesLabel -> {rIt, aIt},
    PlotRange -> All,
    PlotRangeClipping -> False,
    PlotStyle -> Blue,
    PlotOptions[Axes] // Evaluate
  ],
  (* rNat *)
  Graphics @ {Directive[Dashed, natStyle],
    Line @ {{rNat, 0}, {rNat, aNat}}
  },
  Graphics @ {natStyle,
    Text[
      Subscript[rIt, "Nat"] // textStyle,
      {rNat, 0},
      {-0.7, 0.8}
    ]
  },
  (* aNat *)
  Graphics @ {Directive[Dashed, natStyle],
    Line @ {{rNat, aNat}, {0, aNat}}
  },
  Graphics @ {natStyle,
    Text[
      Subscript[aIt, "Nat"] // textStyle,
      {0, aNat},
      {1.3, 0}
    ]
  }
] // PrettyString["Nat" -> "\[Natural]"] // Ex["line-psi.pdf"]


(* ::Section:: *)
(*Viable domain*)


(* ::Subsection:: *)
(*Plot (radial space)*)


Module[{eps},
  eps = 0.07;
  Show[
    EmptyFrame[{0, (1 + eps) aNat}, {0, 1 + eps},
      AspectRatio -> 1,
      FrameLabel -> {aIt, rIt}
    ],
    (* Viable domain *)
    Graphics @ {
      Text[
        "Viable domain" // textStyle,
        {(1 + eps) aNat / 2, (1 + eps) / 2}
      ]
    },
    (* Terminal radii *)
    Plot[{rSharp[a], rFlat[a]}, {a, 0, aNat},
    (* Non-viable domain *)
      Filling -> {1 -> {2}},
      FillingStyle -> nonStyle,
      PlotLegends -> textStyle /@ {
        rIt == Subscript[rIt, "Sharp"],
        rIt == Subscript[rIt, "Flat"]
      },
      PlotRange -> Full,
      PlotStyle -> {sharpStyle, flatStyle}
    ],
    Graphics @ {
      Text[
        "Non-viable domain" // textStyle,
        {(1 + eps) aNat / 3, (1 + eps) / 20}
      ]
    },
    (* Unphysical domain *)
    Plot[1, {a, 0, (1 + 2 eps) aNat},
      Filling -> Top,
      FillingStyle -> unphysStyle,
      PlotStyle -> unphysStyle
    ],
    Graphics @ {White,
      Text[
        "Unphysical domain" // textStyle,
        {(1 + eps) aNat / 2, 1 + eps / 2}
      ]
    }
  ] // PrettyString["Sharp" -> "\[Sharp]", "Flat" -> "\[Flat]"]
] // Ex["line-viable-r.pdf"]


(* ::Subsection:: *)
(*Animations*)


Module[
 {aMin, aMax, aStep,
  aValues,
  eps, rMax, rMaxUnphys, rMaxNon
 },
  (* Values of A *)
  aMin = 1/5;
  aMax = Ceiling[aNat];
  aStep = 1/5;
  aValues = Range[aMin, aMax, aStep];
  aValues = Append[aValues, aNat] // Sort[#, Less] &;
  (* Animation *)
  Table[
    eps = 0.1;
    rMax = 1;
    rMaxUnphys = 1 + eps;
    rMaxNon = If[a < aNat, rSharp[a], rNat];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Transition terminal curve *)
      If[a == aNat,
        ContourPlot[RPolar[x, y] == rNat,
          {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
          ContourLabels -> None,
          ContourStyle -> termStyle
        ],
        {}
      ]
    ]
  , {a, aValues}]
] // Ex["line-viable-full.gif", gifOpts]


Module[
 {aMin, aMax, aStep,
  aValues,
  eps, rMax, rMaxUnphys, rMaxNon
 },
  (* Values of A *)
  aMin = Floor[aNat];
  aMax = Ceiling[aNat];
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  aValues = Append[aValues, aNat] // Sort[#, Less] &;
  (* Animation *)
  Table[
    eps = 0.1;
    rMax = 2 rNat;
    rMaxUnphys = 1 + eps;
    rMaxNon = If[a < aNat, rSharp[a], rNat];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Transition terminal curve *)
      If[a == aNat,
        ContourPlot[RPolar[x, y] == rNat,
          {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
          ContourLabels -> None,
          ContourStyle -> termStyle
        ],
        {}
      ]
    ]
  , {a, aValues}]
] // Ex["line-viable-zoom.gif", gifOpts]


(* ::Section:: *)
(*Traced boundary algebra*)


(* ::Subsection:: *)
(*Near x = 1 (upper branch)*)


With[{a = \[FormalCapitalA], r = \[FormalR]},
  Module[{phiDerivativeAsymptotic, phiAsymptotic},
    phiDerivativeAsymptotic = Series[phiTraDer[a][r], {r, 1, 4}];
    phiAsymptotic = Integrate[phiDerivativeAsymptotic, r];
    {
      {"d\[Phi]/dr", phiDerivativeAsymptotic},
      {"\[Phi]", phiAsymptotic},
      Nothing
    }
      // FullSimplify[#, a > 0] &
      // TableForm
  ]
] // Ex["line-traced-r-near-1.pdf"]


(* ::Subsection:: *)
(*Curvature and radius of inflection*)


(* ::Text:: *)
(*See Page r2-8.*)


With[{r = \[FormalR], phi = \[FormalPhi], a = \[FormalCapitalA]},
  Module[
   {ar, aphi, d, dSimp, vel, acc, comp, cur,
    phiDer, phiDer2, curEval, rInfl,
    brevity
   },
    (* Local orthogonal basis *)
    {ar, aphi};
    (* Abbreviation for r-derivative *)
    d = Dt[#, r, Constants -> a] &;
    (* Simplification rules *)
    dSimp = {
      d[ar] -> phi'[r] aphi,
      d[aphi] -> -phi'[r] ar
    };
    (* Velocity and acceleration *)
    vel = ar + r phi'[r] aphi;
    acc = d[vel] /. dSimp;
    (* Components for local (r, \[Phi], z) basis *)
    comp = {
      ar -> {1, 0, 0},
      aphi -> {0, 1, 0}
    };
    {vel, acc} = {vel, acc} /. comp;
    (* Curvature *)
    cur = Last @ Cross[vel, acc];
    (* \[Phi]' and \[Phi]'' *)
    phiDer = phiTraDer[a][r];
    phiDer2 = d[phiDer] // FullSimplify;
    (* Curvature evaluated *)
    curEval = cur /. {
      phi'[r] -> phiDer,
      phi''[r] -> phiDer2
    } // FullSimplify;
    (* Radius of inflection *)
    rInfl = r /. Solve[curEval == 0 && 0 < r < 1, r];
    (* For brevity *)
    brevity = {Log[r] -> \[FormalCapitalL]};
    (* Table of results *)
    {
      {"Curvature", cur},
      {"phi'", phiDer},
      {"phi''", phiDer2},
      {"Curvature", curEval},
      {"r_i", rInfl}
    } /. brevity
      // PrettyString["phi" -> "\[Phi]"]
      // TableForm
  ]
] // Ex["line-curvature-inflection.pdf"]


(* ::Subsection:: *)
(*Series solution for hot regime outer boundaries*)


(* ::Text:: *)
(*See (r2.33) (Page r2-9).*)


With[{r = \[FormalR], xi = \[FormalXi], l = \[FormalCapitalL]},
  Module[{nMax, psiSeries, a, sqrtSeries, phiDerSeries, phiSeries, negLog},
    nMax = 8;
    (* Series for \[Psi](r_\[Sharp] + \[Xi]) *)
    psiSeries = psi[r + xi] + O[xi] ^ nMax;
    (* By construction, \[Psi](r_\[Sharp]) == A *)
    a = psi[r];
    (* Series for sqrt(A^2 - \[Psi]^2(r_\[Sharp] + \[Xi])) *)
    sqrtSeries = Sqrt[a^2 - psiSeries^2];
    (* Series for d\[Phi]/d\[Xi] (see (r2.23a)) *)
    phiDerSeries = psiSeries / ((r + xi) sqrtSeries);
    (* Series for \[Phi] *)
    phiSeries = Integrate[phiDerSeries, xi];
    (* Simplification *)
    phiSeries = phiSeries /. {Log[r] -> -negLog};
    phiSeries =
      Block[{$Assumptions = rNat < r < 1 && 0 < negLog < 4},
        FullSimplify /@ phiSeries
      ];
    phiSeries = phiSeries /. {negLog -> -Log[r]};
    phiSeries = Factor /@ phiSeries;
    phiSeries /. {
      Log[r] -> Subscript[l, "sharp"],
      r -> Subscript[r, "sharp"]
    } // PrettyString["sharp" -> "\[Sharp]"]
  ]
] // Ex["line-hot-traced-outer-series.pdf"]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsection:: *)
(*Hot regime A < A_\[Natural]*)


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  phiOffsetValues,
  eps, rMaxAss, rMaxUnphys, rMaxNonAss,
  rMax, rMaxNon, phi
 },
  a = aHot;
  phiOffsetValues = Range[0, 360, 30] Degree // Most;
  eps = 0.1;
  rMaxAss = Association[
    "full" -> 1,
    "zoom_1" -> rInfl,
    "zoom_2" -> rFlat[a]
  ];
  rMaxNonAss = Association[
    "full" -> rSharp[a],
    "zoom_1" -> rInfl (1 + eps),
    "zoom_2" -> rFlat[a] (1 + eps)
  ];
  rMaxUnphys = 1 + eps;
  Table[
    rMax = rMaxAss[zoom];
    rMaxNon = rMaxNonAss[zoom];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Traced boundaries *)
      Table[
        phi = phiTraHot[id];
        ParametricPlot[
          {
            XYPolar[r, phi[r] + phiOffset],
            XYPolar[r, -phi[r] + phiOffset]
          } // Evaluate, {r, DomainStart[phi], DomainEnd[phi]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {phiOffset, phiOffsetValues}
      , {id, {"outer", "inner"}}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourLabels -> None,
        ContourStyle -> inflStyle
      ]
    ] // Ex["line-hot-traced-" <> zoom <> ".pdf"]
  , {zoom, Keys[rMaxAss]}]
]


(* ::Subsubsection:: *)
(*Convex (generic) without known solution*)


Module[
 {a, rSh,
  nSpike, modNSpike, mod2Pi,
  rSpike, phiSpike, phiCentre, phiHalfWidth,
  rMax, eps, rMaxNon, phi
 },
  a = aHot;
  rSh = rSharp[a];
  nSpike = 2;
  modNSpike = Mod[#, nSpike, 1] &;
  mod2Pi = Mod[#, 2 Pi] &;
  rSpike[1] = rInfl;
  rSpike[2] = Way[rSh, rInfl, 3/4];
  Table[
    phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
  , {n, nSpike}];
  phiCentre[1] = 0;
  phiCentre[2] = 135 Degree;
  phiCentre[n_] := phiCentre[n // modNSpike];
  Table[
    phiHalfWidth[n] =
      phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
  , {n, nSpike}];
  phiHalfWidth[n_] := phiHalfWidth[n // modNSpike];
  rMax = rInfl;
  eps = 0.1;
  rMaxNon = rInfl (1 + eps);
  Show[
    EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
      PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a] @ RPolar[x, y] < 0,
      {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries (spikes) *)
    Table[
      phi[r_] := phiTraHot["outer"][r] - phiSpike[n];
      ParametricPlot[
        {
          XYPolar[r, phiCentre[n] +  phi[r]],
          XYPolar[r, phiCentre[n] -  phi[r]]
        } // Evaluate, {r, rSh, rSpike[n]},
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Traced boundaries (terminal curve portions) *)
    Table[
      ParametricPlot[
        XYPolar[rSh, ph] // Evaluate,
        {
          ph,
          phiCentre[n] + phiHalfWidth[n] // mod2Pi,
          phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
        },
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Radius of inflection *)
    ContourPlot[RPolar[x, y] == rInfl,
      {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
      ContourLabels -> None,
      ContourStyle -> inflStyle
    ]
  ]
] // Ex["line-hot-traced-convex-generic.pdf"]


(* ::Subsubsection:: *)
(*Convex (generic) with known solution*)


Module[
 {a, rSh,
  nSpike, modNSpike, mod2Pi,
  rSpike, phiSpike, phiCentre, phiHalfWidth,
  rMax, eps, rMaxNon, rBath, phi
 },
  a = aHot;
  rSh = rSharp[a];
  nSpike = 2;
  modNSpike = Mod[#, nSpike, 1] &;
  mod2Pi = Mod[#, 2 Pi] &;
  rSpike[1] = rInfl;
  rSpike[2] = Way[rSh, rInfl, 3/4];
  Table[
    phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
  , {n, nSpike}];
  phiCentre[1] = 0;
  phiCentre[2] = 135 Degree;
  phiCentre[n_] := phiCentre[n // modNSpike];
  Table[
    phiHalfWidth[n] =
      phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
  , {n, nSpike}];
  phiHalfWidth[n_] := phiHalfWidth[n // modNSpike];
  rMax = rInfl;
  eps = 0.1;
  rMaxNon = rInfl (1 + eps);
  rBath = rSh / 2;
  Show[
    EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
      PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a] @ RPolar[x, y] < 0,
      {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[-Log @ RPolar[x, y],
      {x, -rMax, rMax}, {y, -rMax, rMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        RPolar[x, y] > rBath && Or[
          RPolar[x, y] < rSh,
          Sequence @@ Table[
            phi[r_] := phiTraHot["outer"][r] - phiSpike[n];
            Less[
              phiCentre[n] - phi @ RPolar[x, y],
              PhiPolar[x, y] // If[n == 3, mod2Pi, Identity],
              phiCentre[n] + phi @ RPolar[x, y]
            ]
          , {n, nSpike}]
        ] // Evaluate
      ]
    ],
    (* Traced boundaries (spikes) *)
    Table[
      phi[r_] := phiTraHot["outer"][r] - phiSpike[n];
      ParametricPlot[
        {
          XYPolar[r, phiCentre[n] +  phi[r]],
          XYPolar[r, phiCentre[n] -  phi[r]]
        } // Evaluate, {r, rSh, rSpike[n]},
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Traced boundaries (terminal curve portions) *)
    Table[
      ParametricPlot[
        XYPolar[rSh, ph] // Evaluate,
        {
          ph,
          phiCentre[n] + phiHalfWidth[n] // mod2Pi,
          phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
        },
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Radius of inflection *)
    ContourPlot[RPolar[x, y] == rInfl,
      {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
      ContourLabels -> None,
      ContourStyle -> inflStyle
    ]
  ]
] // Ex["line-hot-traced-convex-generic-with-known.png"]


(* ::Subsubsection:: *)
(*Convex (regular) without known solution*)


Module[
 {a, rSh,
  mod2Pi,
  phi,
  rSpike, phiSpike, phiCentre, phiHalfWidth,
  rMax, eps, rMaxNon
 },
  a = aHot;
  rSh = rSharp[a];
  Table[
    mod2Pi = Mod[#, 2 Pi] &;
    phi[r_] := phiTraHot["outer"][r];
    rSpike = SeekRoot[
      phi[(1 + 10^-6) rSh] -  phi[#] - Pi / nSpike &,
      {rSh, rInfl}
    ];
    phiSpike = phi[rSpike];
    phiHalfWidth = Pi / nSpike;
    phiCentre[n_] := n * 2 phiHalfWidth;
    rMax = rInfl;
    eps = 0.1;
    rMaxNon = rInfl (1 + eps);
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Traced boundaries (spikes) *)
      Table[
        phi[r_] := phiTraHot["outer"][r] - phiSpike;
        ParametricPlot[
          {
            XYPolar[r, phiCentre[n] +  phi[r]],
            XYPolar[r, phiCentre[n] -  phi[r]]
          } // Evaluate, {r, rSh, rSpike},
          PlotStyle -> convexStyle
        ]
      , {n, 0, nSpike - 1}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourLabels -> None,
        ContourStyle -> inflStyle
      ]
    ] // Ex[
      "line-hot-traced-convex-regular-"
      <> ToString[nSpike]
      <> ".pdf"
    ]
  , {nSpike, 4, 8}]
]


(* ::Subsubsection:: *)
(*Convex (regular) with known solution*)


Module[
 {a, rSh,
  mod2Pi,
  phi,
  rSpike, phiSpike, phiCentre, phiHalfWidth,
  rMax, eps, rMaxNon, rBath
 },
  a = aHot;
  rSh = rSharp[a];
  Table[
    mod2Pi = Mod[#, 2 Pi] &;
    phi[r_] := phiTraHot["outer"][r];
    rSpike = SeekRoot[
      phi[(1 + 10^-6) rSh] -  phi[#] - Pi / nSpike &,
      {rSh, rInfl}
    ];
    phiSpike = phi[rSpike];
    phiHalfWidth = Pi / nSpike;
    phiCentre[n_] := n * 2 phiHalfWidth;
    rMax = rInfl;
    eps = 0.1;
    rMaxNon = rInfl (1 + eps);
    rBath = rSh / 2;
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
    (* Known solution *)
    DensityPlot[-Log @ RPolar[x, y],
      {x, -rMax, rMax}, {y, -rMax, rMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        RPolar[x, y] > rBath && Or[
          RPolar[x, y] < rSh,
          Sequence @@ Table[
            phi[r_] := phiTraHot["outer"][r] - phiSpike;
            Less[
              phiCentre[n] - phi @ RPolar[x, y],
              PhiPolar[x, y] // If[n > (nSpike - 1) / 2, mod2Pi, Identity],
              phiCentre[n] + phi @ RPolar[x, y]
            ]
          , {n, 0, nSpike - 1}]
        ] // Evaluate
      ]
    ],
      (* Traced boundaries (spikes) *)
      Table[
        phi[r_] := phiTraHot["outer"][r] - phiSpike;
        ParametricPlot[
          {
            XYPolar[r, phiCentre[n] +  phi[r]],
            XYPolar[r, phiCentre[n] -  phi[r]]
          } // Evaluate, {r, rSh, rSpike},
          PlotStyle -> convexStyle
        ]
      , {n, 0, nSpike - 1}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourLabels -> None,
        ContourStyle -> inflStyle
      ]
    ] // Ex[
      "line-hot-traced-convex-regular-"
      <> ToString[nSpike]
      <> "-with-known.png"
    ]
  , {nSpike, 4, 8}]
]


(* ::Subsection:: *)
(*Cold-to-hot transition A = A_\[Natural]*)


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  phiOffsetValues,
  eps, rMaxAss, rMaxNonAss, rMaxUnphys,
  rMax, rMaxNon, phi
 },
  a = aColdHot;
  phiOffsetValues = Range[0, 360, 30] Degree // Most;
  eps = 0.1;
  rMaxAss = Association[
    "full" -> 1,
    "zoom_1" -> rInfl,
    "zoom_2" -> rFlat[a]
  ];
  rMaxNonAss = Association[
    "full" -> rSharp[a],
    "zoom_1" -> rInfl (1 + eps),
    "zoom_2" -> rFlat[a] (1 + eps)
  ];
  rMaxUnphys = 1 + eps;
  Table[
    rMax = rMaxAss[zoom];
    rMaxNon = rMaxNonAss[zoom];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == N[a] // transStyle[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Traced boundaries *)
      Table[
        phi = phiTraColdHot[id];
        ParametricPlot[
          {
            XYPolar[r, phi[r] + phiOffset],
            XYPolar[r, -phi[r] + phiOffset]
          } // Evaluate, {r, DomainStart[phi], DomainEnd[phi]},
          PlotStyle -> {upperStyle, lowerStyle},
          RegionFunction -> Function[{x, y},
            Abs[RPolar[x, y] - rNat] > 10^-3 // Evaluate
          ]
        ]
      , {phiOffset, phiOffsetValues}
      , {id, {"outer", "inner"}}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourLabels -> None,
        ContourStyle -> inflStyle
      ],
      (* Transition terminal curve *)
      ContourPlot[RPolar[x, y] == rNat,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        ContourLabels -> None,
        ContourStyle -> termStyle
      ]
    ] // Ex["line-cold_hot-traced-" <> zoom <> ".pdf"]
  , {zoom, Keys[rMaxAss]}]
]


(* ::Subsection:: *)
(*Cold regime A > A_\[Natural]*)


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  phiOffsetValues,
  eps, rMaxAss, rMaxUnphys,
  rMax, phi
 },
  a = aCold;
  phiOffsetValues = Range[0, 360, 30] Degree // Most;
  eps = 0.1;
  rMaxAss = Association[
    "full" -> 1,
    "zoom_1" -> rInfl,
    "zoom_2" -> 2 rNat
  ];
  rMaxUnphys = 1 + eps;
  Table[
    rMax = rMaxAss[zoom];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Traced boundaries *)
      Table[
        phi = phiTraCold;
        ParametricPlot[
          {
            XYPolar[r, phi[r] + phiOffset],
            XYPolar[r, -phi[r] + phiOffset]
          } // Evaluate, {r, DomainStart[phi], DomainEnd[phi]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {phiOffset, phiOffsetValues}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourLabels -> None,
        ContourStyle -> inflStyle
      ]
    ] // Ex["line-cold-traced-" <> zoom <> ".pdf"]
  , {zoom, Keys[rMaxAss]}]
]


(* ::Section:: *)
(*Numerical verification plots*)


(* ::Subsection:: *)
(*Finite element mesh*)


Module[
 {source,
  a, tBath, mesh, prExt, prInt,
  bCoords, bCoordsExt, bCoordsInt
 },
  (* Import mesh *)
  source = "line-verification-mesh.txt";
  {a, tBath, mesh, prExt, prInt} = Import[source] // Uncompress;
  (* Boundary coordinates *)
  bCoords = Part[
    mesh["Coordinates"],
    List @@ First @ mesh["BoundaryElements"] // Flatten // DeleteDuplicates
  ];
  bCoordsExt = Select[bCoords, prExt @@ # &];
  bCoordsInt = Select[bCoords, prInt @@ # &];
  (* Plot *)
  Show[
    mesh["Wireframe"],
    ListPlot[{bCoordsExt, bCoordsInt},
      PlotStyle -> (Directive[#, pointStyle] & /@ {Blue, Red})
    ],
    ImageSize -> 480
  ]
] // Ex["line-verification-mesh.pdf"]


(* ::Subsection:: *)
(*Numerical solution*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "line-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y], Element[{x, y}, mesh],
      AxesLabel -> Italicise /@ {"x", "y", "T"},
      PlotLabel -> "Numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["line-verification-solution.png"]


(* ::Subsection:: *)
(*Relative error (3D)*)


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "line-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Log[1 / RPolar[x, y]];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["line-verification-rel_error-3d.png"]


(* ::Subsection:: *)
(*Relative error (2D)*)


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "line-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Log[1 / RPolar[x, y]];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Show[
      DensityPlot[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
        ColorFunction -> "Rainbow",
        PlotLabel -> "Relative error of numerical solution",
        PlotRange -> Full,
        PlotLegends -> Automatic,
        PlotOptions[Frame] // Evaluate
      ],
      mesh["Wireframe"]
    ]
  ]
] // Ex["line-verification-rel_error-2d.png"]


(* ::Section:: *)
(*Figure: known solution (line-known-solution)*)


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    rMax, tMax,
    dummyForTrailingCommas
  },
  (* Plot range *)
  rMax = 1.1;
  tMax = 2.5;
  (* Make plot *)
  Plot3D[
    Log[1 / RPolar[x, y]]
    , {x, -rMax, rMax}
    , {y, -rMax, rMax}
    , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
    , AxesLabel -> Italicise /@ {"x", "y", "T"}
    , Boxed -> {Back, Bottom, Left}
    , BoxRatios -> {Automatic, Automatic, Automatic}
    , ClippingStyle -> {BoundaryTracingStyle["Unphysical"], None}
    , ImageSize -> 360
    , LabelStyle -> Directive[Black, 16]
    , Lighting -> GeneralStyle["AmbientLighting"]
    , PlotPoints -> 50
    , PlotRange -> {0, tMax}
    , PlotStyle -> SlidesStyle["InteriorRegion"]
    , TicksStyle -> 12
  ]
] // Ex["line-known-solution-slides.png"];


(* ::Section:: *)
(*Figure: Auxiliary function (line-auxiliary-function)*)


Module[{textStyle},
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    Plot[psi[r], {r, 0, 1}
      , AspectRatio -> 0.5
      , AxesLabel -> {rIt // Margined @ {{3, 0}, {5, 0}}, LaTeXStyle["\[Psi]"]}
      , ImageSize -> 0.75 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , PlotRange -> All
      , PlotRangeClipping -> False
      , PlotRangePadding -> {None, {0.05, Automatic}}
      , PlotStyle -> Directive[Black, Thick]
      , PlotOptions[Axes] // Evaluate
      , TicksStyle -> LabelSize["Tick"]
    ],
    Graphics @ {
      (* Guiding lines for critical values *)
      BoundaryTracingStyle["Contour"],
      Line @ {{0, aNat}, {rNat, aNat}, {rNat, 0}},
      (* rNat *)
      Text[
        Subscript[rIt, "Nat"] // textStyle,
        {rNat, 0},
        {-0.2, 0.9}
      ],
      (* aNat *)
      Text[
        Subscript[aIt, "Nat"] // textStyle,
        {0, aNat},
        {1.5, 0}
      ],
      {}
    },
    {}
  ] // PrettyString["Nat" -> "\[Natural]"]
] // Ex["line-auxiliary-function.pdf"]


(* ::Section:: *)
(*Figure: viable domain (line-viable)*)


Module[
  {
    aStep, aValues, aMin,
    xMax, yMin, yMax,
    plotList,
    rMaxNon,
    textStyle, arrowStyle,
    parameterArrow,
    xGraphicsANat,
    legendLabelStyle,
    dummyForTrailingCommas
  },
  (* Values of A *)
  aStep = 3/100;
  aValues = aNat + aStep * Range[-1, 1];
  aMin = Min[aValues];
  (* Plot range *)
  xMax = 1.25 rSharp[aMin];
  yMin = -rSharp[aMin];
  yMax = 1.2 rSharp[aMin];
  (* List of plots *)
  plotList =
    Table[
      rMaxNon = If[a < aNat, rSharp[a], rNat];
      Show[
        EmptyFrame[{-xMax, xMax}, {yMin, yMax}
          , Frame -> None
          , FrameLabel -> None
          , FrameTicks -> None
          , PlotRangePadding -> {Automatic, {None, Automatic}}
        ],
        Which[
          (* Hot regime *)
          a < aNat,
          {
            Graphics @ {
              FaceForm @ BoundaryTracingStyle["NonViable"],
              EdgeForm @ BoundaryTracingStyle["Terminal"],
              Annulus @ {rFlat[a], rSharp[a]}
            },
            Graphics @ {
              Text[
                rIt == Subscript[rIt, "\[Flat]"] // textStyle
                , {0, rFlat[a]}
                , {0, 1.2}
             ]
            },
            Graphics @ {
              Text[
                rIt == Subscript[rIt, "\[Sharp]"] // textStyle
                , {0, rSharp[a]}
                , {0, -1}
             ]
            },
            {}
          },
          (* Transition *)
          a == aNat,
          {
            Graphics @ {
              BoundaryTracingStyle["Terminal"],
              Circle[{0, 0}, rNat]
            },
            Graphics @ {
              Text[
                rIt == Subscript[rIt, "\[Natural]"] // textStyle
                , {0, rNat}
                , {0, -1}
             ]
            },
            {}
          },
          (* (Nothing to plot for cold regime) *)
          True, {}
        ],
        {}
      ]
      , {a, aValues}
    ];
  (* Parameter (B) increase indicator arrow *)
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
          aIt // textStyle,
          {1, 0},
          {-2, 0}
        ]
      },
      (* A == A_nat *)
      xGraphicsANat = 0.515;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsANat, 0},
          {xGraphicsANat, -0.01}
        }
      },
      Graphics @ {
        Text[
          Subscript[Italicise["A"], "\[Natural]"]
            // textStyle
          , {xGraphicsANat, 0}
          , {0, 1.4}
        ]
      },
      (* Regime labels *)
      Graphics @ {
        Text[
          "hot" // textStyle
          , {Way[0, xGraphicsANat, 0.31], 0}
          , {0, -1.3}
        ],
        Text[
          "cold" // textStyle
          , {Way[xGraphicsANat, 1, 0.5], 0}
          , {0, -1.3}
        ],
        {}
      },
      {}
      , PlotRange -> All
    ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  (* Export *)
  {
    GraphicsRow[plotList
      , ImageSize -> 0.8 ImageSizeTextWidth
      , Spacings -> 0
    ]
      // Ex["line-viable.pdf"]
    ,
    Show[parameterArrow
      , ImageSize -> 0.8 ImageSizeTextWidth
    ]
      // Ex["line-viable-arrow.pdf"]
    ,
    GraphicsGrid[
      List @ Join[
        RegionLegend[
          {BoundaryTracingStyle["NonViable"]},
          {"non\[Hyphen]viable domain"}
          , LabelStyle -> legendLabelStyle
        ],
        RegionLegend[
          {BoundaryTracingStyle["Viable"]},
          {"viable domain"}
          , LabelStyle -> legendLabelStyle
        ],
        CurveLegend[
          {BoundaryTracingStyle["Terminal"]},
          {"terminal curve"}
          , LabelStyle -> legendLabelStyle
        ],
        {}
      ]
      , Alignment -> Left
      , ImageSize -> ImageSizeTextWidth
      , ItemAspectRatio -> 0.1
      , Spacings -> {{0, 0, -0.12} ImageSizeTextWidth, 0}
    ]
      // Ex["line-viable-legend.pdf"]
    ,
    Nothing
  }
]


(* ::Section:: *)
(*Figure: traced boundaries (traced-boundaries)*)


Module[
  {
    aStep, aValues, aMin,
    imageSize, rMaxShow,
    tracedStyle,
    rMaxNon,
    inequality,
    regimeName,
    dummyForTrailingCommas
  },
  (* Values of A *)
  aStep = 0.13;
  aValues = aNat + aStep * Range[-1, 1];
  aMin = Min[aValues];
  (* Plot range *)
  imageSize = 0.325 ImageSizeTextWidth;
  rMaxShow = 1.35 rSharp[aMin];
  (* Slightly less thick traced boundaries *)
  tracedStyle = (
    BoundaryTracingStyle["Traced"]
      /. {AbsoluteThickness[d_] :> AbsoluteThickness[0.9 d]}
  );
  (* Plots *)
  Table[
    rMaxNon = If[a < aNat, rSharp[a], rNat];
    inequality = Which[a < aNat, Less, a == aNat, Equal, a > aNat, Greater];
    regimeName = Which[a < aNat, "hot", a == aNat, "cold_hot", a > aNat, "cold"];
    Show[
      EmptyFrame[{-rMaxShow, rMaxShow}, {-rMaxShow, rMaxShow}
        , FrameLabel -> {
            Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
            Italicise["y"]
          }
        , FrameStyle -> LabelSize["Axis"]
        , FrameTicksStyle -> LabelSize["Tick"]
        , LabelStyle -> LatinModernLabelStyle[LabelSize["Label"] - 1]
        , ImageSize -> imageSize
        , PlotLabel -> inequality[aIt, Subscript[aIt, "\[Natural]"]]
      ],
      (* Non-viable domain *)
      If[a < aNat,
        Graphics @ {
          EdgeForm @ BoundaryTracingStyle["Terminal"],
          FaceForm @ BoundaryTracingStyle["NonViable"],
          Annulus[{0, 0}, {rFlat[a], rSharp[a]}]
        },
        {}
      ],
      (* Traced boundaries *)
      Which[
        (* Hot regime || Cold-to-hot transition *)
        a <= aNat,
        {
          (* Inner viable island *)
          With[{r = \[FormalR]},
            Module[{eps, rMin, rMax, rInit, phiMin, phiMax, rTraced, phiOffsetValues},
              eps = 10.^-4;
              rMin = eps;
              rMax = (1 - eps) rFlat[a];
              rInit = Way[rMin, rMax];
              phiMin = -2 Pi;
              phiMax = 2 Pi;
              rTraced = NDSolveValue[
                {
                  r'[phi] == 1 / phiTraDer[a][r[phi]],
                  r[0] == rInit,
                  WhenEvent[{r[phi] < rMin, r[phi] > rMax}, "StopIntegration"]
                },
                r, {phi, phiMin, phiMax},
                NoExtrapolation
              ];
              phiOffsetValues = Subdivide[0, 2 Pi, 3] // Most;
              (* Traced boundaries *)
              Table[
                ParametricPlot[
                  {
                    XYPolar[rTraced[phi], phi + phiOffset],
                    XYPolar[rTraced[phi], -phi + phiOffset]
                  } // Evaluate, {phi, phiMin, phiMax}
                  , PlotPoints -> 3
                  , PlotStyle -> tracedStyle
                ]
              , {phiOffset, phiOffsetValues}]
            ]
          ],
          (* Outer viable mainland *)
          With[{r = \[FormalR]},
            Module[{eps, rMin, rMax, rInit, phiMin, phiMax, rTraced, phiOffsetValues},
              eps = 10.^-4;
              rMin = (1 + eps) rSharp[a];
              rMax = 1.5 rMaxShow;
              rInit = Way[rMin, rMax];
              phiMin = -Pi;
              phiMax = 1.5 Pi;
              rTraced = NDSolveValue[
                {
                  r'[phi] == 1 / phiTraDer[a][r[phi]],
                  r[0] == rInit,
                  WhenEvent[{r[phi] < rMin, r[phi] > rMax}, "StopIntegration"]
                },
                r, {phi, phiMin, phiMax},
                NoExtrapolation
              ];
              phiOffsetValues = Subdivide[0, 2 Pi, 6] // Most;
              (* Traced boundaries *)
              Table[
                ParametricPlot[
                  {
                    XYPolar[rTraced[phi], phi + phiOffset],
                    XYPolar[rTraced[phi], -phi + phiOffset]
                  } // Evaluate, {phi, phiMin, phiMax}
                  , PlotPoints -> 3
                  , PlotStyle -> tracedStyle
                ]
              , {phiOffset, phiOffsetValues}]
            ]
          ],
          {}
        },
        (* Cold regime *)
        a > aNat,
        {
          With[{r = \[FormalR]},
            Module[{eps, rMin, rMax, rInit, phiMin, phiMax, rTraced, phiOffsetValues},
              eps = 10.^-4;
              rMin = eps;
              rMax = 1.5 rMaxShow;
              rInit = Way[rMin, rMax];
              phiMin = -Pi;
              phiMax = 2.5 Pi;
              rTraced = NDSolveValue[
                {
                  r'[phi] == 1 / phiTraDer[a][r[phi]],
                  r[0] == rInit,
                  WhenEvent[{r[phi] < rMin, r[phi] > rMax}, "StopIntegration"]
                },
                r, {phi, phiMin, phiMax},
                NoExtrapolation
              ];
              phiOffsetValues = Subdivide[0, 2 Pi, 3] // Most;
              (* Traced boundaries *)
              Table[
                ParametricPlot[
                  {
                    XYPolar[rTraced[phi], phi + phiOffset],
                    XYPolar[rTraced[phi], -phi + phiOffset]
                  } // Evaluate, {phi, phiMin, phiMax}
                  , PlotPoints -> 3
                  , PlotStyle -> tracedStyle
                ]
              , {phiOffset, phiOffsetValues}]
            ]
          ]
        },
        True, {}
      ],
      {}
    ] // Ex @ FString["line-traced-boundaries-{regimeName}.pdf"]
    , {a, aValues}
  ]
]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    aStep, aValues, aMin,
    imageSize, rMaxShow,
    tracedStyle,
    rMaxNon,
    inequality,
    regimeName,
    dummyForTrailingCommas
  },
  (* Values of A *)
  aStep = 0.13;
  aValues = aNat + aStep * Range[-1, 1];
  aMin = Min[aValues];
  (* Plot range *)
  imageSize = 0.33 ImageSizeTextWidthBeamer;
  rMaxShow = 1.35 rSharp[aMin];
  (* Slightly less thick traced boundaries *)
  tracedStyle = (
    Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
      /. {AbsoluteThickness[d_] :> AbsoluteThickness[0.9 d]}
  );
  (* Plots *)
  Table[
    rMaxNon = If[a < aNat, rSharp[a], rNat];
    inequality = Which[a < aNat, Less, a == aNat, Equal, a > aNat, Greater];
    regimeName = Which[a < aNat, "hot", a == aNat, "cold_hot", a > aNat, "cold"];
    Show[
      EmptyFrame[{-rMaxShow, rMaxShow}, {-rMaxShow, rMaxShow}
        , Frame -> None
        , ImageSize -> imageSize
      ],
      (* Non-viable domain *)
      If[a < aNat,
        Graphics @ {
          EdgeForm @ tracedStyle,
          FaceForm @ BoundaryTracingStyle["NonViable"],
          Annulus[{0, 0}, {rFlat[a], rSharp[a]}]
        },
        {}
      ],
      (* Traced boundaries *)
      Which[
        (* Hot regime || Cold-to-hot transition *)
        a <= aNat,
        {
          (* Inner viable island *)
          With[{r = \[FormalR]},
            Module[{eps, rMin, rMax, rInit, phiMin, phiMax, rTraced, phiOffsetValues},
              eps = 10.^-4;
              rMin = eps;
              rMax = (1 - eps) rFlat[a];
              rInit = Way[rMin, rMax];
              phiMin = -2 Pi;
              phiMax = 2 Pi;
              rTraced = NDSolveValue[
                {
                  r'[phi] == 1 / phiTraDer[a][r[phi]],
                  r[0] == rInit,
                  WhenEvent[{r[phi] < rMin, r[phi] > rMax}, "StopIntegration"]
                },
                r, {phi, phiMin, phiMax},
                NoExtrapolation
              ];
              phiOffsetValues = Subdivide[0, 2 Pi, 3] // Most;
              (* Traced boundaries *)
              Table[
                ParametricPlot[
                  {
                    XYPolar[rTraced[phi], phi + phiOffset],
                    XYPolar[rTraced[phi], -phi + phiOffset]
                  } // Evaluate, {phi, phiMin, phiMax}
                  , PlotPoints -> 3
                  , PlotStyle -> tracedStyle
                ]
              , {phiOffset, phiOffsetValues}]
            ]
          ],
          (* Outer viable mainland *)
          With[{r = \[FormalR]},
            Module[{eps, rMin, rMax, rInit, phiMin, phiMax, rTraced, phiOffsetValues},
              eps = 10.^-4;
              rMin = (1 + eps) rSharp[a];
              rMax = 1.5 rMaxShow;
              rInit = Way[rMin, rMax];
              phiMin = -Pi;
              phiMax = 1.5 Pi;
              rTraced = NDSolveValue[
                {
                  r'[phi] == 1 / phiTraDer[a][r[phi]],
                  r[0] == rInit,
                  WhenEvent[{r[phi] < rMin, r[phi] > rMax}, "StopIntegration"]
                },
                r, {phi, phiMin, phiMax},
                NoExtrapolation
              ];
              phiOffsetValues = Subdivide[0, 2 Pi, 6] // Most;
              (* Traced boundaries *)
              Table[
                ParametricPlot[
                  {
                    XYPolar[rTraced[phi], phi + phiOffset],
                    XYPolar[rTraced[phi], -phi + phiOffset]
                  } // Evaluate, {phi, phiMin, phiMax}
                  , PlotPoints -> 3
                  , PlotStyle -> tracedStyle
                ]
              , {phiOffset, phiOffsetValues}]
            ]
          ],
          {}
        },
        (* Cold regime *)
        a > aNat,
        {
          With[{r = \[FormalR]},
            Module[{eps, rMin, rMax, rInit, phiMin, phiMax, rTraced, phiOffsetValues},
              eps = 10.^-4;
              rMin = eps;
              rMax = 1.5 rMaxShow;
              rInit = Way[rMin, rMax];
              phiMin = -Pi;
              phiMax = 2.5 Pi;
              rTraced = NDSolveValue[
                {
                  r'[phi] == 1 / phiTraDer[a][r[phi]],
                  r[0] == rInit,
                  WhenEvent[{r[phi] < rMin, r[phi] > rMax}, "StopIntegration"]
                },
                r, {phi, phiMin, phiMax},
                NoExtrapolation
              ];
              phiOffsetValues = Subdivide[0, 2 Pi, 3] // Most;
              (* Traced boundaries *)
              Table[
                ParametricPlot[
                  {
                    XYPolar[rTraced[phi], phi + phiOffset],
                    XYPolar[rTraced[phi], -phi + phiOffset]
                  } // Evaluate, {phi, phiMin, phiMax}
                  , PlotPoints -> 3
                  , PlotStyle -> tracedStyle
                ]
              , {phiOffset, phiOffsetValues}]
            ]
          ]
        },
        True, {}
      ],
      {}
    ] // Ex @ FString["line-traced-boundaries-{regimeName}-slides.pdf"]
    , {a, aValues}
  ]
];


(* ::Section:: *)
(*Figure: critical terminal radii (line-critical)*)


Module[
  {
    aMax,
    textStyle,
    plotStyles, plot,
    rPatch, aPatch,
    legendLabelStyleCurves, legendLabelStyleRegions,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  aMax = Ceiling[aNat] - 0.1;
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Plot *)
  plotStyles = {
    Directive[Black, Thick],
    Directive[Black, Thick, GeneralStyle["Dashed"]]
  };
  plot = Show[
    (* Terminal radii *)
    Plot[{rSharp[a], rFlat[a]}, {a, 0, aNat}
      (* Non-viable domain *)
      , AspectRatio -> 1.05
      , AxesLabel -> {aIt, rIt}
      , Filling -> {1 -> {2}}
      , FillingStyle -> BoundaryTracingStyle["NonViable"]
      , ImageSize -> 0.7 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , PlotRange -> {{0, aMax}, Full}
      , PlotRangeClipping -> False
      , PlotRangePadding -> {{0.05, Automatic}, {0, Automatic}}
      , PlotStyle -> plotStyles
      , TicksStyle -> LabelSize["Tick"]
    ],
    (* Patch missing piece near r == 1 *)
    rPatch = 0.9877;
    aPatch = psi[rPatch];
    ParametricPlot[{a, rSharp[a]}, {a, 0, aPatch}
      , PlotPoints -> 2
      , PlotStyle -> plotStyles[[1]]
    ],
    (* Critical radii *)
    Graphics @ {
      (* Guiding lines *)
      BoundaryTracingStyle["Contour"],
      Line @ {{0, rNat}, {aNat, rNat}, {aNat, 0}},
      (* rNat *)
      Text[
        Subscript[rIt, "\[Natural]"] // textStyle,
        {0, rNat},
        {2, -0.15}
      ],
      (* rNat *)
      Text[
        Subscript[aIt, "\[Natural]"] // textStyle,
        {aNat, 0},
        {0, 1.1}
      ],
      {}
    },
    {}
  ];
  (* Legend *)
  legendLabelStyleCurves = LatinModernLabelStyle @ LabelSize["Legend"];
  legendLabelStyleRegions = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      plotStyles,
      Row @ {rIt, "\[VeryThinSpace]"} == Subscript[rIt, #] &
        /@ {"\[VeryThinSpace]\[VeryThinSpace]\[Sharp]", "\[VeryThinSpace]\[VeryThinSpace]\[VeryThinSpace]\[Flat]"}
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
    Show[plot
      , ImageSize -> 0.5 ImageSizeTextWidth
    ] // Ex["line-critical.pdf"]
    ,
    GraphicsColumn[
      Join[legendCurves, legendRegions]
      , Alignment -> Left
      , ImageSize -> 0.2 ImageSizeTextWidth
      , ItemAspectRatio -> 0.2
      , Spacings -> -0.02 ImageSizeTextWidth
    ] // Ex["line-critical-legend.pdf"]
  }
]


(* ::Section:: *)
(*Figure: cold regime concave corner (line-traced-boundaries-cold-concave)*)


Module[
  {
    a,
    imageSize, rMaxShow,
    rCorner, phiCorner,
    phiMin, phiMax, rTraced,
    xyUpper, xyLower,
    textStyle,
    pUpperLabel, pLowerLabel,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = 1.1 aNat;
  (* Plot range *)
  imageSize = 240;
  rMaxShow = 2 rNat;
  (* Location of concave corner *)
  rCorner = 0.1 rMaxShow;
  phiCorner = 213 Degree;
  (* Traced boundaries *)
  phiMin = -0.83 Pi;
  phiMax = 0;
  With[{r = \[FormalR]},
    rTraced = NDSolveValue[
      {r'[phi] == 1 / phiTraDer[a][r[phi]], r[0] == rCorner},
      r, {phi, phiMin, phiMax},
      NoExtrapolation
    ];
  ];
  xyUpper[p_] := XYPolar[rTraced[#], # + phiCorner] & @ Way[phiMin, phiMax, p];
  xyLower[p_] := XYPolar[rTraced[#], -# + phiCorner] & @ Way[phiMin, phiMax, 1 - p];
  (* Text style *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Branch label proportions of traced boundary *)
  pUpperLabel = 0.05;
  pLowerLabel = 0.59;
  (* Plot *)
  Show[
    (* Traced boundaries *)
    ParametricPlot[
      {xyUpper[p], xyLower[p]} // Evaluate
      , {p, 0, 1}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{-1, 0}, {3, 0}},
          Italicise["y"]
        }
      , AxesStyle -> Darker[Gray]
      , ImageSize -> 0.4 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
        (* NOTE: LatinModernLabelStyle[16] results in cut-off x-label *)
      , Method -> {"AxesInFront" -> False}
        (* https://mathematica.stackexchange.com/a/26918 *)
      , PlotStyle -> BoundaryTracingStyle["Traced"]
      , Ticks -> None
    ]
      /. line_Line :> {Arrowheads[{0, 0.08, 0}], Arrow[line]}
      (* https://mathematica.stackexchange.com/a/102182 *)
    ,
    (* Branch labels *)
    Graphics @ {
      Text[
        "upper" // textStyle
        , xyUpper[pUpperLabel]
        , {0, 1.1}
      ],
      Text[
        "lower" // textStyle
        , xyLower[pLowerLabel]
        , {0, -1.6}
      ],
      {}
    },
    {}
  ]
] // Ex["line-traced-boundaries-cold-concave.pdf"]


(* ::Section:: *)
(*Figure: hot regime concave corner on inner viable island (line-traced-boundaries-hot-inner-concave.pdf)*)


Module[
  {
    a,
    imageSize, rMaxShow, rMaxNon,
    rCorner, phiCorner,
    phiMin, phiMax, rTraced,
    xyUpper, xyLower,
    textStyle, textStyleBigger,
    pUpperLabel, pLowerLabel,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = (1 - 4/1000) aNat;
  (* Plot range *)
  imageSize = 0.4 ImageSizeTextWidth;
  rMaxShow = 1 rSharp[a];
  rMaxNon = rSharp[a];
  (* Location of concave corner *)
  rCorner = 0.27 rFlat[a];
  phiCorner = 222 Degree;
  (* Traced boundaries *)
  phiMin = -0.75 Pi;
  phiMax = 0;
  With[{r = \[FormalR]},
    rTraced = NDSolveValue[
      {r'[phi] == 1 / phiTraDer[a][r[phi]], r[0] == rCorner},
      r, {phi, phiMin, phiMax},
      NoExtrapolation
    ];
  ];
  xyUpper[p_] := XYPolar[rTraced[#], # + phiCorner] & @ Way[phiMin, phiMax, p];
  xyLower[p_] := XYPolar[rTraced[#], -# + phiCorner] & @ Way[phiMin, phiMax, 1 - p];
  (* Text style *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleBigger = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Branch label proportions of traced boundary *)
  pUpperLabel = 0.28;
  pLowerLabel = 0.55;
  (* Plot *)
  Show[
    EmptyFrame[{-rMaxShow, rMaxShow}, {-rMaxShow, rMaxShow}
      , Frame -> None
      , ImageSize -> imageSize
      , PlotRange -> All
      , PlotRangeClipping -> False
    ],
    (* Non-viable domain *)
    Graphics @ {
      FaceForm @ BoundaryTracingStyle["NonViable"],
      EdgeForm @ BoundaryTracingStyle["Terminal"],
      Annulus @ {rFlat[a], rSharp[a]}
    },
    (* Traced boundaries *)
    ParametricPlot[
      {xyUpper[p], xyLower[p]} // Evaluate
      , {p, 0, 1}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ] /. line_Line :> {Arrowheads[{0, {0.07, 0.52}, 0}], Arrow[line]},
    (* Branch labels *)
    Graphics @ {
      Text[
        "upper" // textStyle
        , xyUpper[pUpperLabel]
        , {-1.4, 0}
      ],
      Text[
        "lower" // textStyle
        , xyLower[pLowerLabel]
        , {0, -1.5}
      ],
      {}
    },
    (* Terminal curve labels *)
    Graphics @ {
      Text[
        rIt == Subscript[rIt, "\[Flat]"] // textStyleBigger
        , {0, -rFlat[a]}
        , {0, -1.3}
     ]
    },
    Graphics @ {
      Text[
        rIt == Subscript[rIt, "\[Sharp]"] // textStyleBigger
        , {0, -rSharp[a]}
        , {0, 0.65}
     ]
    },
    {}
  ]
] // Ex["line-traced-boundaries-hot-inner-concave.pdf"]


(* ::Section:: *)
(*Figure: hot regime outer viable mainland tangent line (line-hot-outer-tangent-line.pdf)*)


Module[
  {
    rSharp, phi, xi,
    yMax,
    angleMarkerLength,
    orthogonalityMarkerLength, orthogonalityMarkerStyle,
    textStyle,
    dummyForTrailingCommas
  },
  (* Coordinates *)
  rSharp = 1;
  phi = 40 Degree;
  xi = rSharp (Sec[phi] - 1); (* see (r2.24) *)
  (* Vertical extent of tangent line *)
  yMax = 1.2 rSharp;
  (* Styles *)
  angleMarkerLength = 0.25 rSharp;
  orthogonalityMarkerLength = 0.1 rSharp;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[None]];
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Plot *)
  Show[
    Graphics @ {
      (* Circle r == r_sharp *)
      {
        BoundaryTracingStyle["Terminal"],
        Circle[{0, 0}, rSharp]
      },
      (* Radius unto tangent line *)
      Line @ {{0, 0}, {rSharp, 0}},
      Text[
        Subscript[rIt, "\[Sharp]"] // textStyle
        , {rSharp/1.9, 0}
        , {0, 0.6}
      ],
      (* Radius at azimuthal angle *)
      Line @ {{0, 0}, XYPolar[rSharp + xi, phi]},
      Text[
        Subscript[rIt, "\[Sharp]"] // textStyle
        , XYPolar[rSharp/1.9, phi]
        , XYPolar[1.5, 0.8 phi - Pi/2]
      ],
      Text[
        "\[Xi]" // textStyle
        , XYPolar[rSharp + xi/2, phi]
        , XYPolar[1.4, phi - Pi/2]
      ],
      (* Azimuthal angle marker *)
      Circle[{0, 0}, angleMarkerLength, {0, phi}],
      Text[
        "\[Phi]" // textStyle
        , XYPolar[angleMarkerLength, phi/2]
        , XYPolar[2, 0.62 phi/2 + Pi]
      ],
      (* Tangent line *)
      {
        BoundaryTracingStyle["Traced"],
        Line @ {{rSharp, -yMax}, {rSharp, yMax}}
      },
      (* Orthogonality marker *)
      Line[
        {rSharp, 0} + orthogonalityMarkerLength * # & /@{
          {0, 1}, {-1, 1}, {-1, 0}
        }
      ],
      {}
    }
    , ImageSize -> 0.33 ImageSizeTextWidth
  ]
] // Ex["line-hot-outer-tangent-line.pdf"]


(* ::Section:: *)
(*Figure: concave spike near r = 1 (line-traced-boundaries-hot-spike)*)


Module[
  {
    a, rFl, rSh,
    xMin, xMax, yMax,
    phiTracedRaw, rRawStart, rRawEnd,
    phiRawEnd, phiTraced,
    phiTerminal,
    rInflection, phiInflection,
    textStyle,
    dummyForTrailingCommas
  },
  (* Constants *)
  a = aHot;
  rFl = rFlat[a];
  rSh = rSharp[a];
  (* Plot range *)
  xMax = 1;
  yMax = 1.7 rSh;
  xMin = -1.5 rSh;
  (* Traced boundary *)
  phiTracedRaw = phiTraHot["outer"];
  rRawStart = DomainStart[phiTracedRaw];
  rRawEnd = DomainEnd[phiTracedRaw];
  (* Enforce phi == 0 at r == 1 *)
  phiRawEnd = phiTracedRaw[rRawEnd];
  phiTraced[r_] := phiTracedRaw[r] - phiRawEnd // Evaluate;
  (* Endpoint phi for terminal curve traced boundary portion *)
  phiTerminal = phiTraced[rRawStart] // Abs;
  (* Inflection coordinates *)
  rInflection = rInfl;
  phiInflection = phiTraced[rInflection] // Abs;
  (* Text style *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Make plot *)
  Show[
    EmptyAxes[{xMin, xMax}, {-yMax, yMax}
      , AspectRatio -> Automatic
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , ImageSize -> 0.7 ImageSizeTextWidth
      , PlotRange -> All
      , PlotRangePadding -> {{Scaled[0.1], Automatic}, Scaled[0.2]}
      , PlotRangeClipping -> None
      , Ticks -> {Automatic, None}
      , TicksStyle -> LabelSize["Tick"]
    ],
    (* Non-viable domain *)
    Graphics @ {
      FaceForm @ BoundaryTracingStyle["NonViable"],
      EdgeForm @ BoundaryTracingStyle["Terminal"],
      Annulus @ {rFl, rSh}
    },
    (* Traced boundaries (spike) *)
    ParametricPlot[
      XYPolar[r, phiTraced[r]]
        // IncludeYReflection
        // Evaluate
      , {r, rRawStart, rRawEnd}
      , PlotPoints -> Automatic
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    Graphics @ {
      Text[
        "concave" // textStyle
        , {0.5, 0}
        , {0, -1.4}
      ],
      {}
    },
    (* Traced boundaries (terminal curve portion) *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Circle[{0, 0}, rSh, {phiTerminal, 2 Pi - phiTerminal}]
    },
    Graphics @ {
      Text[
        "convex" // textStyle
        , XYPolar[rSh, Way[phiTerminal, 2 Pi - phiTerminal, 0.2]]
        , {1, -1.1}
      ],
      {}
    },
    (* Point of inflection *)
    Graphics @ {GeneralStyle["Point"],
      Point @ IncludeYReflection @ XYPolar[rInflection, phiInflection]
    },
    Graphics @ {
      Text[
        "inflection" // textStyle
        , XYPolar[rInflection, phiInflection]
        , {-0.7, -1.4}
      ],
      {}
    },
    {}
  ] // Ex["line-traced-boundaries-hot-spike.pdf"]
]


(* ::Section:: *)
(*Figure: hot regime protrusion*)
(*(line-traced-boundaries-hot-protrusion.pdf)*)


Module[
  {
    a,
    imageSize, rMaxShow, rSh,
    rCorner, phiCorner,
    eps, rMin, phiTraced,
    xyUpper, xyLower,
    textStyle, textStyleBigger,
    pUpperLabel, pLowerLabel,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = 0.75 aNat;
  (* Plot range *)
  imageSize = 0.33 ImageSizeTextWidth;
  rMaxShow = rInfl;
  rSh = rSharp[a];
  (* Location of (convex) protrusion corner *)
  rCorner = Way[rSh, rInfl, 2/3];
  phiCorner = 90 Degree;
  (* Traced boundaries *)
  eps = 10.^-6;
  rMin = (1 + eps) rSh;
  With[{phi = \[FormalPhi]},
    phiTraced = NDSolveValue[
      {phi'[r] == phiTraDer[a][r], phi[rCorner] == 0},
      phi, {r, rMin, rCorner},
      NoExtrapolation
    ];
  ];
  xyUpper[p_] := XYPolar[#, phiTraced[#] + phiCorner] & @ Way[rMin, rCorner, 1 - p];
  xyLower[p_] := XYPolar[#, -phiTraced[#] + phiCorner] & @ Way[rMin, rCorner, p];
  (* Text style *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleBigger = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Branch label proportions of traced boundary *)
  pUpperLabel = 0.7;
  pLowerLabel = 0.27;
  (* Plot *)
  Show[
    EmptyFrame[{-rMaxShow, rMaxShow}, {-rMaxShow, rMaxShow}
      , Frame -> None
      , ImageSize -> imageSize
      , PlotRange -> All
      , PlotRangeClipping -> False
    ],
    (* Outer terminal curve r == r_sharp *)
    Graphics @ {
      BoundaryTracingStyle["Terminal"],
      Circle[{0, 0}, rSh]
    },
    (* Frontier of inflection r == r_infl *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      Circle[{0, 0}, rInfl]
    },
    (* Traced boundaries (protrusion) *)
    ParametricPlot[
      xyLower[p]
      , {p, 0, 1}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ] /. line_Line :> {Arrowheads[{0, {0.08, 0.75}, 0}], Arrow[line]},
    ParametricPlot[
      xyUpper[p]
      , {p, 0, 1}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ] /. line_Line :> {Arrowheads[{0, {0.08, 0.52}, 0}], Arrow[line]},
    (* Branch labels *)
    Graphics @ {
      Text[
        "upper" // textStyle
        , xyUpper[pUpperLabel]
        , {0, -1.3}
        , -xyUpper'[pUpperLabel]
      ],
      Text[
        "lower" // textStyle
        , xyLower[pLowerLabel]
        , {0, -1.1}
        , -xyLower'[pLowerLabel]
      ],
      {}
    },
    (* Radii labels *)
    Graphics @ {
      Text[
        rIt == Subscript[rIt, "\[Sharp]"] // textStyleBigger
        , {0, -rSh}
        , {0, -1.5}
     ]
    },
    Graphics @ {
      Text[
        rIt == Subscript[rIt, "i"] // textStyleBigger
        , {0, -rInfl}
        , {0, -1.4}
     ]
    },
    {}
  ]
] // Ex["line-traced-boundaries-hot-protrusion.pdf"]


(* ::Section:: *)
(*Figure: convex domains (line-domains.pdf)*)


(*
  Basically the following:
  - "line-hot-traced-convex-regular-5.pdf"
  - "line-hot-traced-convex-regular-4.pdf"
  - "line-hot-traced-convex-generic.pdf"
*)


Module[
  {
    a, rSh,
    caseList,
    mod2Pi, phi,
    xMax, yMax,
    numSpikes, modNumSpikes,
    rSpike, phiSpike, phiCentre, phiHalfWidth,
    plotList, caseIsRegular,
    phiPlotted,
    dummyForTrailingCommas
  },
  (* Constants *)
  a = aHot;
  rSh = rSharp[a];
  (* Cases to be plotted *)
  caseList = {"pentagon", "square", "generic"};
  (* Abbreviations *)
  mod2Pi = Mod[#, 2 Pi] &;
  phi[r_] := phiTraHot["outer"][r];
  (* List of plots *)
  xMax = rInfl;
  yMax = 0.9 rInfl;
  plotList = Table[
    (* Define spikes *)
    ClearAll[numSpikes, rSpike, phiSpike, phiHalfWidth, phiCentre];
    caseIsRegular = case != "generic";
    Which[
      caseIsRegular,
        numSpikes = <|"pentagon" -> 5, "square" -> 4|> @ case;
        rSpike[n_] := SeekRoot[
          phi[(1 + 10^-6) rSh] - phi[#] - Pi / numSpikes &,
          {rSh, rInfl}
        ] // Evaluate;
        phiSpike[n_] = phi[rSpike[1]];
        phiHalfWidth[n_] = Pi / numSpikes;
        phiCentre[n_] := (n - 1) * 2 phiHalfWidth[n];
      ,
      Not[caseIsRegular],
        numSpikes = 2;
        modNumSpikes = Mod[#, numSpikes, 1] &;
        rSpike[1] = rInfl;
        rSpike[2] = Way[rSh, rInfl, 3/4];
        Table[
          phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
          , {n, numSpikes}
        ];
        phiCentre[1] = 0;
        phiCentre[2] = 135 Degree;
        phiCentre[n_] := phiCentre[n // modNumSpikes];
        Table[
          phiHalfWidth[n] =
            phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
          , {n, numSpikes}
        ];
        phiHalfWidth[n_] := phiHalfWidth[n // modNumSpikes];
      ,
      True, {}
    ];
    (* Plot of constructed domain *)
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax}
        , Frame -> None
        , PlotRangePadding -> {{Scaled[-0.07], Scaled[0.03]}, None}
      ],
      (* Outer terminal curve (effective incircle) *)
      Graphics @ {
        BoundaryTracingStyle["Terminal", "Background"],
        Circle[{0, 0}, rSh]
      },
      (* Traced boundaries (spikes) *)
      Table[
        phiPlotted[r_] := phiTraHot["outer"][r] - phiSpike[n];
        ParametricPlot[
          {
            XYPolar[r, phiCentre[n] + phiPlotted[r]],
            XYPolar[r, phiCentre[n] - phiPlotted[r]]
          } // Evaluate
          , {r, rSh, rSpike[n]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
        , {n, numSpikes}
      ],
      (* Traced boundaries (terminal curve portions) *)
      If[Not[caseIsRegular],
        Table[
          ParametricPlot[
            XYPolar[rSh, ph] // Evaluate,
            {
              ph,
              phiCentre[n] + phiHalfWidth[n] // mod2Pi,
              phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
            }
            , PlotPoints -> 2
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ]
          , {n, numSpikes}
        ]
        , {}
      ],
      {}
    ]
    , {case, caseList}
  ];
  (* Final figure *)
  GraphicsRow[plotList
    , Alignment -> Left
    , ImageSize -> 0.95 ImageSizeTextWidth
    , PlotRangePadding -> 0
    , Spacings -> 0
  ]
] // Ex["line-domains.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    a, rSh,
    caseList,
    mod2Pi, phi,
    xMax, yMax,
    numSpikes, modNumSpikes,
    rSpike, phiSpike, phiCentre, phiHalfWidth,
    plotList, caseIsRegular,
    phiPlotted,
    dummyForTrailingCommas
  },
  (* Constants *)
  a = aHot;
  rSh = rSharp[a];
  (* Cases to be plotted *)
  caseList = {"pentagon", "square", "generic"};
  (* Abbreviations *)
  mod2Pi = Mod[#, 2 Pi] &;
  phi[r_] := phiTraHot["outer"][r];
  (* List of plots *)
  xMax = rInfl;
  yMax = 0.9 rInfl;
  plotList = Table[
    (* Define spikes *)
    ClearAll[numSpikes, rSpike, phiSpike, phiHalfWidth, phiCentre];
    caseIsRegular = case != "generic";
    Which[
      caseIsRegular,
        numSpikes = <|"pentagon" -> 5, "square" -> 4|> @ case;
        rSpike[n_] := SeekRoot[
          phi[(1 + 10^-6) rSh] - phi[#] - Pi / numSpikes &,
          {rSh, rInfl}
        ] // Evaluate;
        phiSpike[n_] = phi[rSpike[1]];
        phiHalfWidth[n_] = Pi / numSpikes;
        phiCentre[n_] := (n - 1) * 2 phiHalfWidth[n];
      ,
      Not[caseIsRegular],
        numSpikes = 2;
        modNumSpikes = Mod[#, numSpikes, 1] &;
        rSpike[1] = rInfl;
        rSpike[2] = Way[rSh, rInfl, 3/4];
        Table[
          phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
          , {n, numSpikes}
        ];
        phiCentre[1] = 0;
        phiCentre[2] = 135 Degree;
        phiCentre[n_] := phiCentre[n // modNumSpikes];
        Table[
          phiHalfWidth[n] =
            phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
          , {n, numSpikes}
        ];
        phiHalfWidth[n_] := phiHalfWidth[n // modNumSpikes];
      ,
      True, {}
    ];
    (* Plot of constructed domain *)
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax}
        , Frame -> None
        , PlotRangePadding -> {{Scaled[-0.07], Scaled[0.03]}, None}
      ],
      (* Outer terminal curve (effective incircle) *)
      Graphics @ {
        BoundaryTracingStyle["Terminal", "Background"],
        Circle[{0, 0}, rSh]
      },
      (* Traced boundaries (spikes) *)
      Table[
        phiPlotted[r_] := phiTraHot["outer"][r] - phiSpike[n];
        ParametricPlot[
          {
            XYPolar[r, phiCentre[n] + phiPlotted[r]],
            XYPolar[r, phiCentre[n] - phiPlotted[r]]
          } // Evaluate
          , {r, rSh, rSpike[n]}
          , PlotPoints -> 2
          , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
        ]
        , {n, numSpikes}
      ],
      (* Traced boundaries (terminal curve portions) *)
      If[Not[caseIsRegular],
        Table[
          ParametricPlot[
            XYPolar[rSh, ph] // Evaluate,
            {
              ph,
              phiCentre[n] + phiHalfWidth[n] // mod2Pi,
              phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
            }
            , PlotPoints -> 2
            , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
          ]
          , {n, numSpikes}
        ]
        , {}
      ],
      (* Line-source singularity *)
      Graphics @ {
        SlidesStyle["Source"],
        Point @ {0, 0}
      },
      {}
    ]
    , {case, caseList}
  ];
  (* Final figure *)
  GraphicsRow[plotList
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidthBeamer
    , PlotRangePadding -> 0
    , Spacings -> 0
  ]
] // Ex["line-domains-slides.pdf"]


(* ::Subsection:: *)
(*Version for slides, without terminal curve*)


Module[
  {
    a, rSh,
    caseList,
    mod2Pi, phi,
    xMax, yMax,
    numSpikes, modNumSpikes,
    rSpike, phiSpike, phiCentre, phiHalfWidth,
    plotList, caseIsRegular,
    phiPlotted,
    dummyForTrailingCommas
  },
  (* Constants *)
  a = aHot;
  rSh = rSharp[a];
  (* Cases to be plotted *)
  caseList = {"pentagon", "square", "generic"};
  (* Abbreviations *)
  mod2Pi = Mod[#, 2 Pi] &;
  phi[r_] := phiTraHot["outer"][r];
  (* List of plots *)
  xMax = rInfl;
  yMax = 0.9 rInfl;
  plotList = Table[
    (* Define spikes *)
    ClearAll[numSpikes, rSpike, phiSpike, phiHalfWidth, phiCentre];
    caseIsRegular = case != "generic";
    Which[
      caseIsRegular,
        numSpikes = <|"pentagon" -> 5, "square" -> 4|> @ case;
        rSpike[n_] := SeekRoot[
          phi[(1 + 10^-6) rSh] - phi[#] - Pi / numSpikes &,
          {rSh, rInfl}
        ] // Evaluate;
        phiSpike[n_] = phi[rSpike[1]];
        phiHalfWidth[n_] = Pi / numSpikes;
        phiCentre[n_] := (n - 1) * 2 phiHalfWidth[n];
      ,
      Not[caseIsRegular],
        numSpikes = 2;
        modNumSpikes = Mod[#, numSpikes, 1] &;
        rSpike[1] = rInfl;
        rSpike[2] = Way[rSh, rInfl, 3/4];
        Table[
          phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
          , {n, numSpikes}
        ];
        phiCentre[1] = 0;
        phiCentre[2] = 135 Degree;
        phiCentre[n_] := phiCentre[n // modNumSpikes];
        Table[
          phiHalfWidth[n] =
            phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
          , {n, numSpikes}
        ];
        phiHalfWidth[n_] := phiHalfWidth[n // modNumSpikes];
      ,
      True, {}
    ];
    (* Plot of constructed domain *)
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax}
        , Frame -> None
        , PlotRangePadding -> {{Scaled[-0.07], Scaled[0.03]}, None}
      ],
      (* Traced boundaries (spikes) *)
      Table[
        phiPlotted[r_] := phiTraHot["outer"][r] - phiSpike[n];
        ParametricPlot[
          {
            XYPolar[r, phiCentre[n] + phiPlotted[r]],
            XYPolar[r, phiCentre[n] - phiPlotted[r]]
          } // Evaluate
          , {r, rSh, rSpike[n]}
          , PlotPoints -> 2
          , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
        ]
        , {n, numSpikes}
      ],
      (* Traced boundaries (terminal curve portions) *)
      If[Not[caseIsRegular],
        Table[
          ParametricPlot[
            XYPolar[rSh, ph] // Evaluate,
            {
              ph,
              phiCentre[n] + phiHalfWidth[n] // mod2Pi,
              phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
            }
            , PlotPoints -> 2
            , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
          ]
          , {n, numSpikes}
        ]
        , {}
      ],
      (* Line-source singularity *)
      Graphics @ {
        SlidesStyle["Source"],
        Point @ {0, 0}
      },
      {}
    ]
    , {case, caseList}
  ];
  (* Final figure *)
  GraphicsRow[plotList
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidthBeamer
    , PlotRangePadding -> 0
    , Spacings -> 0
  ]
] // Ex["line-domains-slides_plain.pdf"]


(* ::Subsection:: *)
(*Version for slides, with non-singular sources*)


Module[
  {
    a, rSh,
    caseList,
    mod2Pi, phi,
    xMax, yMax,
    numSpikes, modNumSpikes,
    rSpike, phiSpike, phiCentre, phiHalfWidth,
    plotList, caseIsRegular,
    phiPlotted,
    dummyForTrailingCommas
  },
  (* Constants *)
  a = aHot;
  rSh = rSharp[a];
  (* Cases to be plotted *)
  caseList = {"pentagon", "square", "generic"};
  (* Abbreviations *)
  mod2Pi = Mod[#, 2 Pi] &;
  phi[r_] := phiTraHot["outer"][r];
  (* List of plots *)
  xMax = rInfl;
  yMax = 0.9 rInfl;
  plotList = Table[
    (* Define spikes *)
    ClearAll[numSpikes, rSpike, phiSpike, phiHalfWidth, phiCentre];
    caseIsRegular = case != "generic";
    Which[
      caseIsRegular,
        numSpikes = <|"pentagon" -> 5, "square" -> 4|> @ case;
        rSpike[n_] := SeekRoot[
          phi[(1 + 10^-6) rSh] - phi[#] - Pi / numSpikes &,
          {rSh, rInfl}
        ] // Evaluate;
        phiSpike[n_] = phi[rSpike[1]];
        phiHalfWidth[n_] = Pi / numSpikes;
        phiCentre[n_] := (n - 1) * 2 phiHalfWidth[n];
      ,
      Not[caseIsRegular],
        numSpikes = 2;
        modNumSpikes = Mod[#, numSpikes, 1] &;
        rSpike[1] = rInfl;
        rSpike[2] = Way[rSh, rInfl, 3/4];
        Table[
          phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
          , {n, numSpikes}
        ];
        phiCentre[1] = 0;
        phiCentre[2] = 135 Degree;
        phiCentre[n_] := phiCentre[n // modNumSpikes];
        Table[
          phiHalfWidth[n] =
            phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
          , {n, numSpikes}
        ];
        phiHalfWidth[n_] := phiHalfWidth[n // modNumSpikes];
      ,
      True, {}
    ];
    (* Plot of constructed domain *)
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax}
        , Frame -> None
        , PlotRangePadding -> {{Scaled[-0.07], Scaled[0.03]}, None}
      ],
      (* Traced boundaries (spikes) *)
      Table[
        phiPlotted[r_] := phiTraHot["outer"][r] - phiSpike[n];
        ParametricPlot[
          {
            XYPolar[r, phiCentre[n] + phiPlotted[r]],
            XYPolar[r, phiCentre[n] - phiPlotted[r]]
          } // Evaluate
          , {r, rSh, rSpike[n]}
          , PlotPoints -> 2
          , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
        ]
        , {n, numSpikes}
      ],
      (* Traced boundaries (terminal curve portions) *)
      If[Not[caseIsRegular],
        Table[
          ParametricPlot[
            XYPolar[rSh, ph] // Evaluate,
            {
              ph,
              phiCentre[n] + phiHalfWidth[n] // mod2Pi,
              phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
            }
            , PlotPoints -> 2
            , PlotStyle -> Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]]
          ]
          , {n, numSpikes}
        ]
        , {}
      ],
      (* Constant-temperature boundary *)
      Graphics @ {
        BoundaryTracingStyle["Contour"], SlidesStyle["Source"],
        Circle[{0, 0}, 0.5 rSh]
      },
      {}
    ]
    , {case, caseList}
  ];
  (* Final figure *)
  GraphicsRow[plotList
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidthBeamer
    , PlotRangePadding -> 0
    , Spacings -> 0
  ]
] // Ex["line-domains-slides_sources.pdf"];


(* ::Section:: *)
(*Figure: verification domain and mesh (line-verification-domain-mesh.pdf)*)


Module[
  {
    a, rSh,
    rBath, tBath,
    mod2Pi, phi,
    xMax, yMax, imageSize,
    textStyle,
    numSpikes, modNumSpikes,
    domainPlot,
    rSpike, phiSpike, phiCentre, phiHalfWidth,
    phiPlotted,
    mesh,
    meshGraphic,
    legendCurves, legend,
    dummyForTrailingCommas
  },
  (* Constants *)
  a = aHot;
  rSh = rSharp[a];
  (* Heat bath (Dirichlet condition) *)
  rBath = rSh / 2;
  tBath = Log[1 / rBath];
  (* Abbreviations *)
  mod2Pi = Mod[#, 2 Pi] &;
  phi[r_] := phiTraHot["outer"][r];
  (* Plot range *)
  xMax = rInfl;
  yMax = 0.9 rInfl;
  imageSize = 240;
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Define spikes *)
  numSpikes = 2;
  modNumSpikes = Mod[#, numSpikes, 1] &;
  rSpike[1] = rInfl;
  rSpike[2] = Way[rSh, rInfl, 3/4];
  Table[
    phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
    , {n, numSpikes}
  ];
  phiCentre[1] = 0;
  phiCentre[2] = 135 Degree;
  phiCentre[n_] := phiCentre[n // modNumSpikes];
  Table[
    phiHalfWidth[n] =
      phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
    , {n, numSpikes}
  ];
  phiHalfWidth[n_] := phiHalfWidth[n // modNumSpikes];
  (* Domain plot *)
  domainPlot = Show[
    EmptyFrame[{-xMax, xMax}, {-yMax, yMax}
      , Frame -> None
      , ImageSize -> imageSize
      , PlotRange -> All
    ],
    (* Heat bath (Dirichlet condition) *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      Circle[{0, 0}, rBath]
    },
    (* Traced boundaries (spikes) *)
    Table[
      phiPlotted[r_] := phiTraHot["outer"][r] - phiSpike[n];
      ParametricPlot[
        {
          XYPolar[r, phiCentre[n] + phiPlotted[r]],
          XYPolar[r, phiCentre[n] - phiPlotted[r]]
        } // Evaluate
        , {r, rSh, rSpike[n]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {n, numSpikes}
    ],
    (* Traced boundaries (terminal curve portions) *)
    Table[
      ParametricPlot[
        XYPolar[rSh, ph] // Evaluate,
        {
          ph,
          phiCentre[n] + phiHalfWidth[n] // mod2Pi,
          phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
        }
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {n, numSpikes}
    ],
    (* Heat bath (Dirichlet condition) labels *)
    Graphics @ {
      Text[
        tIt == Subscript[tIt, "\[NegativeVeryThinSpace]d"] // textStyle
        (* NOTE: negative thick space for better kerning *)
        , {0, rBath}
        , {0, 1.6}
     ]
    },
    Graphics @ {
      Text[
        rIt == Subscript[rIt, "\[VeryThinSpace]d"] // textStyle
        , {0, -rBath}
        , {0, -1.6}
     ]
    },
    {}
  ];
  (* Import mesh *)
  mesh = Import["line-verification-mesh.txt"] // Uncompress // #[[3]]&;
  (* Mesh graphic *)
  meshGraphic =
    Show[
      mesh["Wireframe"]
      , ImageSize -> imageSize
    ];
  (* Legend *)
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced", "Contour"},
      {"radiation", "constant temperature"}
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
    ];
  {
    GraphicsRow[{domainPlot, meshGraphic}
      , ImageSize -> 0.8 ImageSizeTextWidth
    ] // Ex["line-verification-domain-mesh.pdf"]
    ,
    GraphicsRow[legendCurves
      , Alignment -> Left
      , ImageSize -> ImageSizeTextWidth
      , ItemAspectRatio -> 0.05
      , Spacings -> -0.5 ImageSizeTextWidth
    ] // Ex["line-verification-domain-mesh-legend.pdf"]
  }
]


(* ::Subsection:: *)
(*Number of mesh elements*)


Module[{mesh},
  mesh = Import["line-verification-mesh.txt"] // Uncompress // #[[3]]&;
  mesh
]


(* ::Subsection:: *)
(*Maximum relative error throughout mesh*)


Module[{source, tSol, mesh, tExact, relError},
  (* Import solution *)
  source = "line-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Log[1 / RPolar[x, y]];
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tExact[x, y] - 1;
  relError @@@ mesh["Coordinates"] // Abs // Max // PercentForm
]


(* ::Section:: *)
(*Figure: verification solution and relative error (line-verification-*)*)


(* ::Subsection:: *)
(*Solution*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "line-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  Module[{x, y},
    Plot3D[tSol[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{0, 30}, {0, -20}},
          Italicise["y"] // Margined @ {{5, 0}, {0, 0}},
          None
        }
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> {Automatic, Automatic, 0.11}
      , ImageSize -> 0.42 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotLabel -> (Style["Numerical solution", LabelSize["Axis"] - 1] // Margined @ {{0, 0}, {0, 0}})
      , PlotRange -> Full
      , PlotRangePadding -> {Scaled /@ {0.04, 0.01}, Scaled[0.05], Scaled /@ {0.1, 0}}
      , PlotStyle -> Directive[GeneralStyle["Translucent"], BoundaryTracingStyle["Solution3D"]]
      , TicksStyle -> LabelSize["Tick"]
      , ViewPoint -> {1.6, -2.4, 1.5}
    ]
  ]
] // Ex["line-verification-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Subsection:: *)
(*Relative error*)


Module[{source, tSol, mesh, tExact, relError},
  (* Import solution *)
  source = "line-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Log[1 / RPolar[x, y]];
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tExact[x, y] - 1 // Abs;
  Module[{x, y},
    Plot3D[relError[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{0, 30}, {0, -20}},
          Italicise["y"] // Margined @ {{5, 0}, {0, 0}},
          None
        }
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> {Automatic, Automatic, 0.13}
      , ImageSize -> 0.42 {1, 594/691} ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotLabel -> (Style["Relative error", LabelSize["Axis"] - 1] // Margined @ {{0, 10}, {0, 0}})
      , PlotRange -> Full
      , PlotRangePadding -> {Scaled /@ {0.04, 0.01}, Scaled[0.05], Scaled /@ {0, 0.1}}
      , PlotStyle -> Directive[GeneralStyle["Translucent"], BoundaryTracingStyle["Solution3D"]]
      , TicksStyle -> LabelSize["Tick"]
      , ViewPoint -> {1.46, -2.4, 1.5}
    ]
  ]
] // Ex["line-verification-relative-error.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Physical range*)


(* ::Text:: *)
(*See (r2.45) and (r2.47).*)


Module[
  {
    r, eps, k, sigma,
    tMin, tMax,
    pMin, pMax,
    dummyForTrailingCommas
  },
  (* Inradius *)
  r = Quantity[4, "Centimeters"];
  (* PVC emissivity and conductivity *)
  (*
    See Table 2 of
    M. Lucchi, M. Lorenzini (2018).
    Transient analysis of the radiative heating of rotating PVC pipes
    in a oven for end-forming process.
    Applied Thermal Engineering 129, 84--92.
  *)
  eps = 0.93;
  k = Quantity[0.18, "Watts" / "Meters" / "Kelvins"];
  (* Stefan--Boltzmann constant *)
  sigma = Quantity[5.67 10^-8, "Watts" / "Meters"^2 / "Kelvins"^4];
  (* Temperature range *)
  tMin = 1/2^(2/3) * (k / (eps sigma r))^(1/3);
  tMax = 1/2^(1/3) * (k / (eps sigma r))^(1/3);
  (* Power per length *)
  pMin = Pi k^(4/3) / (2^(5/3) (eps sigma r)^(1/3));
  pMax = Pi k^(4/3) / (2^(1/3) (eps sigma r)^(1/3));
  (* Pretty table *)
  {
    {"Inradius", r},
    {"Emissivity", eps},
    {"Conductivity", k},
    {"Temperature", tMin, tMax},
    {"", UnitConvert[tMin, "DegreesCelsius"], UnitConvert[tMax, "DegreesCelsius"]},
    {"Power per length", pMin, pMax}
  } // TableForm
] // Ex["line-physical-range.pdf"]
