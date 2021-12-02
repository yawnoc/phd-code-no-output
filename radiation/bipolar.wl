(* ::Package:: *)

(* ::Text:: *)
(*See r3 (manuscripts/radiation-3-bipolar.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "bipolar"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Critical terminal points along u = \[Pi]*)


(* ::Subsubsection:: *)
(*(v_\[Natural]\[Pi], A_\[Natural]\[Pi])*)
(*(Warm-to-hot transition hyperbolic coordinate and dimensionless group)*)


(* ::Text:: *)
(*See Section 4.2.1 of r3 (Page r3-8), particularly (r3.43) and (r3.44).*)


vNatPi =
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v Sinh[v] / 4,
        v > 0
      ]
    ]
  ][0];


aNatPi = Function[{v}, 4 v^3 / Sinh[v]] @ vNatPi;


{vNatPi, aNatPi} // N


(* ::Subsubsection:: *)
(*v_\[Flat]\[Pi] and v_\[Sharp]\[Pi]*)


(* ::Text:: *)
(*See Sections 4.1.1 and 4.1.2 of r3 (Page r3-6).*)


(* The argument a to vFlatPi and vSharpPi must be EXACT. *)


vFlatPi[a_] /; a <= aNatPi :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v^4 / a,
        0 < v <= vNatPi
      ] // Evaluate
    ]
  ][0] // Evaluate;


vSharpPi[a_] /; a <= aNatPi :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v^4 / a,
        v >= vNatPi
      ] // Evaluate
    ]
  ][0] // Evaluate;


(* Check warm-to-hot transition (u = \[Pi]) *)
vFlatPi[aNatPi] == vNatPi == vSharpPi[aNatPi]


(* ::Subsection:: *)
(*Critical terminal points along u = 0*)


(* ::Subsubsection:: *)
(*(v_\[Natural]0, A_\[Natural]0)*)
(*(Cold-to-warm transition hyperbolic coordinate and dimensionless group)*)


(* ::Text:: *)
(*See Section 4.2.2 of r3 (Page r3-8), particularly (r3.49) and (r3.50).*)


vNat0 =
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v Sinh[v] / 4,
        v > 0
      ]
    ]
  ][0];


aNat0 = Function[{v}, 4 v^3 / Sinh[v]] @ vNat0;


{vNat0, aNat0} // N


(* ::Subsubsection:: *)
(*v_\[Flat]0 and v_\[Sharp]0*)


(* ::Text:: *)
(*See Sections 4.1.1 to 4.1.4 of r3 (Pages r3-6 and r3-7).*)


(* The argument a to vFlat0 and vSharp0 must be EXACT. *)


vFlat0[a_] /; a <= aNat0 :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v^4 / a,
        0 < v <= vNat0
      ] // Evaluate
    ]
  ][0] // Evaluate;


vSharp0[a_] /; a <= aNat0 :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v^4 / a,
        v >= vNat0
      ] // Evaluate
    ]
  ][0] // Evaluate;


(* Check cold-to-warm transition (u = 0) *)
vFlat0[aNat0] == vNat0 == vSharp0[aNat0]


(* ::Subsection:: *)
(*\[CapitalPhi] (viability)*)


(* ::Text:: *)
(*See (r3.39) (Page r3-5).*)


vi[a_] := Function[{u, v},
  (Cosh[v] - Cos[u])^2 - v^8 / a^2
];


(* ::Subsection:: *)
(*\[Omega] = \[Omega](v_\[Sharp]0)*)


(* ::Text:: *)
(*See (r3.86) (Page r3-14).*)


omega[v_] := (Cosh[v] - 1) / (v Sinh[v]);


(* ::Subsection:: *)
(*Representative values of A*)


aHot = 9;
aWarmHot = aNatPi;
aWarm = 9 + 1/2;
aColdWarm = aNat0;
aCold = 10;
aHot < aWarmHot < aWarm < aColdWarm < aCold


(* ::Subsection:: *)
(*Traced boundary derivative dv/du*)


(* ::Text:: *)
(*See (r3.53) (Page r3-9).*)


vTraDer[a_] := Function[{u, v},
  (a / v^4) Sqrt @ vi[a][u, v] // Evaluate
] // Evaluate;


(* ::Subsection:: *)
(*Traced boundary curvature*)


(* ::Text:: *)
(*See (r3.63) (Page r3-10).*)


curTra[a_] := Function[{u, v},
  Module[{d, dd, h, cur, vDer, vDer2},
    (* Abbreviation for u-derivative *)
    d = Dt[#, u, Constants -> a] &;
    (* Abbreviations *)
    h = HBipolar[u, v];
    dd = Sinh[v] - d[v] Sin[u];
    (* Curvature *)
    cur = h^3 (dd (1 + d[v]^2) + d @ d[v] / h);
    (* v' and v'' *)
    vDer = vTraDer[a][u, v];
    vDer2 = d[vDer];
    (* Curvature evaluated *)
    cur
      /. {d @ d[v] -> vDer2}
      /. {d[v] -> vDer}
  ] // Evaluate
] // Evaluate;


(* Version for lower branch *)
curTraLower[a_] := Function[{u, v},
  Module[{d, dd, h, cur, vDer, vDer2},
    (* Abbreviation for u-derivative *)
    d = Dt[#, u, Constants -> a] &;
    (* Abbreviations *)
    h = HBipolar[u, v];
    dd = Sinh[v] - d[v] Sin[u];
    (* Curvature *)
    cur = h^3 (dd (1 + d[v]^2) + d @ d[v] / h);
    (* v' and v'' *)
    vDer = -vTraDer[a][u, v];
    vDer2 = d[vDer];
    (* Curvature evaluated *)
    cur
      /. {d @ d[v] -> vDer2}
      /. {d[v] -> vDer}
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*v_i (inflection hyperbolic coordinate along u = -\[Pi])*)


(* ::Text:: *)
(*See (r3.65) (Page r3-11).*)


vInfl = InverseFunction[
    Function[{v},
      ConditionalExpression[
        v/2 Tanh[v/2] - 1,
        v > 0
      ] // Evaluate
    ]
  ][0];


(* ::Subsubsection:: *)
(*Inner candidate boundary corner coordinate (at u = -\[Pi])*)


(* ::Text:: *)
(*See the bottom of Page r3-9.*)


(* The hyperbolic coordinate (v) at u == -\[Pi] *)
vTraCandCorn[a_] /; 0 < a <= aNat0 :=
  With[{v = \[FormalV]},
    NDSolveValue[
      {
        v'[u] == Re @ vTraDer[a][u, v[u]],
        v[0] == vSharp0[a],
        WhenEvent[vi[a][u, v[u]] < 0, "StopIntegration"]
      }, v[-Pi], {u, -Pi, 0},
      NoExtrapolation,
      PreciseOptions[]
    ]
  ];


(* ::Subsubsection:: *)
(*Outer candidate boundary corner coordinate (at u = -\[Pi])*)


(* ::Text:: *)
(*From v = v_\[Flat]0, not considered in radiation-3-bipolar.pdf.*)


(* The hyperbolic coordinate (v) at u == -\[Pi] *)
vTraCandCornOuter[a_] /; 0 < a <= aNat0 :=
  With[{v = \[FormalV]},
    NDSolveValue[
      {
        v'[u] == Re @ vTraDer[a][u, v[u]],
        v[0] == vFlat0[a],
        WhenEvent[vi[a][u, v[u]] < 0, "StopIntegration"]
      }, v[-Pi], {u, -Pi, 0},
      NoExtrapolation,
      PreciseOptions[]
    ]
  ];


(* ::Subsubsection:: *)
(*A_i (inflection dimensionless group)*)


(* ::Text:: *)
(*See (r3.68) and (r3.69) (Page r3-11).*)


(* Compute A_i using the bisection algorithm *)
aInfl = Module[{dest, a, aMin, aMax, num},
  (* (This is really slow (~30 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "bipolar-a-inflection.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {a, num} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    (* (argument to vTraCandCorn must be EXACT) *)
    aMin = aNat0 - 4/1000;
    aMax = aNat0;
    {a, num} = SeekRootBisection[
      vTraCandCorn[#] - vInfl &,
      {aMin, aMax},
      10^-7,
      N -> False,
      "ReturnIterations" -> True
    ];
    {a, num} // Compress // Ex[dest]
  ];
  (* Print iterations used *)
  Print["Bisection algorithm: " <> ToString[num] <> " iterations"];
  (* Return value of A_i *)
  a
];
aInfl // N


(* ::Subsubsection:: *)
(*Asymmetry*)


(* ::Text:: *)
(*See (r3.92) (Page r3-16).*)


asymm[vSharp0_, vCorn_] :=
  Divide[
    1 - XBipolar[-Pi, vCorn],
    XBipolar[0, vSharp0] - 1
  ] // FullSimplify // Evaluate;


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* ::Subsubsection:: *)
(*System of ODES for terminal curve \[CapitalPhi] = 0*)


(* ::Text:: *)
(*See (57) and (58) (Page 13 of manuscripts/boundary-tracing.pdf)*)
(*for arc length parametrisation of a contour.*)


viContourSystem[a_, sSign_: 1] :=
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{h, fun, p, q, slope},
      (* Scale factor for bipolar coordinates *)
      h = HBipolar[u, v];
      (* Scalar function \[CapitalPhi] whose contours are sought *)
      fun = vi[a][u, v];
      (* Components of the gradient vector *)
      p = (1 / h) D[fun, u];
      q = (1 / h) D[fun, v];
      (* Magnitude of the gradient vector *)
      slope = Sqrt[p^2 + q^2];
      (* Return system of ODEs *)
      {
        u' == q / (h slope),
        v' == -p / (h slope)
      } /. {
        u' -> Sign[sSign] u'[s],
        v' -> Sign[sSign] v'[s],
        u -> u[s],
        v -> v[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


(* Starting points which are hyperbolic critical terminal points *)
(* (the two critical terminal points along u == 0) *)
startUVHot["hyperbolic"] =
  Module[{a},
    a = aHot;
    {
      {0, vSharp0[a]},
      {0, vFlat0[a]}
    } // N // Rationalize[#, 0] &
  ];


(* Starting points along outer terminal curve *)
startUVHot["outer"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, uvTerm},
      a = aHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      smax = 2 Pi (XBipolar[0, vFlat0[a]] - 1);
      (* Outer terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a, -1],
            u[0] == 0, v[0] == vFlat0[a],
            WhenEvent[u[s] == 2 Pi, "StopIntegration"]
          }, {u, v}, {s, 0, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[0, smax, nmax]}] // Most
    ]
  ];


(* Starting points along inner terminal curve *)
startUVHot["inner"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, uvTerm},
      a = aHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      smax = 2 Pi (XBipolar[0, vSharpPi[a]] - 1);
      (* Inner terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a],
            u[0] == 0, v[0] == vSharp0[a],
            WhenEvent[u[s] == 2 Pi, "StopIntegration"]
          }, {u, v}, {s, 0, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[0, smax, nmax]}] // Most
    ]
  ];


(* Starting points along axis *)
startUVHot["axis"] =
  Module[{a, xOuter, xValues},
    a = aHot;
    xOuter = XBipolar[0, vFlat0[a]];
    xValues =
      Join[
        Subdivide[1, 17/16, 4] // Rest,
        {9/8, 5/4, 3/2, 2}
      ] xOuter;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


(* Starting points which are hyperbolic critical terminal points *)
(* (the two critical terminal points along u == 0) *)
startUVWarmHot["hyperbolic"] =
  Module[{a},
    a = aWarmHot;
    {
      {0, vSharp0[a]},
      {0, vFlat0[a]}
    } // N // Rationalize[#, 0] &
  ];


(* Starting points along outer terminal curve *)
startUVWarmHot["outer"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, vPert, uvTerm},
      a = aWarmHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      (* (between -smax and smax to avoid pincer at vNatPi) *)
      smax = Pi (XBipolar[0, vFlat0[a]] - 1);
      (* Prevent sharp turn at vNatPi *)
      vPert = 10^-6;
      (* Outer terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a, -1],
            u[0] == 0, v[0] == vFlat0[a] - vPert,
            WhenEvent[Abs @ u[s] == Pi, "StopIntegration"]
          }, {u, v}, {s, -smax, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[-smax, smax, nmax]}]
        // Most
        // RotateLeft[#, nmax / 2] &
    ]
  ];


(* Starting points along inner terminal curve *)
startUVWarmHot["inner"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, vPert, uvTerm},
      a = aWarmHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      (* (between -smax and smax to avoid pincer at vNatPi) *)
      smax = Pi (XBipolar[0, vSharpPi[a]] - 1);
      (* Prevent sharp turn at vNatPi *)
      vPert = 10^-6;
      (* Inner terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a],
            u[0] == 0, v[0] == vSharp0[a] + vPert,
            WhenEvent[Abs @ u[s] == Pi, "StopIntegration"]
          }, {u, v}, {s, -smax, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[-smax, smax, nmax]}]
        // Most
        // RotateLeft[#, nmax / 2] &
    ]
  ];


(* Starting points along axis *)
startUVWarmHot["axis"] =
  Module[{a, xOuter, xValues},
    a = aWarmHot;
    xOuter = XBipolar[0, vFlat0[a]];
    xValues =
      Join[
        Subdivide[1, 17/16, 4] // Rest,
        {9/8, 5/4, 3/2, 2}
      ] xOuter;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


(* Starting points which are hyperbolic critical terminal points *)
(* (the two critical terminal points along u == 0) *)
startUVWarm["hyperbolic"] =
  Module[{a},
    a = aWarm;
    {
      {0, vSharp0[a]},
      {0, vFlat0[a]}
    } // N // Rationalize[#, 0] &
  ];


(* Starting points along terminal curve *)
startUVWarm["terminal"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, uvTerm},
      a = aWarm;
      nmax = 16;
      (* (Probable) upper bound for arc length traversed *)
      smax = 2 Pi (XBipolar[0, vFlat0[a]] - 1);
      (* Outer terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a, -1],
            u[0] == 0, v[0] == vFlat0[a],
            WhenEvent[u[s] == 0 && u'[s] > 0, "StopIntegration"]
          }, {u, v}, {s, 0, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[0, smax, nmax]}] // Most
    ]
  ];


(* Starting points along axis *)
startUVWarm["axis"] =
  Module[{a, xOuter, xValues},
    a = aWarm;
    xOuter = XBipolar[0, vFlat0[a]];
    xValues =
      Join[
        Subdivide[1, 17/16, 4] // Rest,
        {9/8, 5/4, 3/2, 2}
      ] xOuter;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


(* Starting point which is a hyperbolic critical terminal point *)
(* (the lone critical terminal point along u == 0) *)
startUVColdWarm["hyperbolic"] = {{0, vNat0}} // N // Rationalize[#, 0] &;


(* Starting points along terminal curve *)
(* (the terminal curve is degenerate: it is simply the point {0, vNat0}) *)


(* Starting points along axis *)
startUVColdWarm["axis"] =
  Module[{nmax, xTerminal, xValues},
    nmax = 16;
    xTerminal = XBipolar[0, vNat0];
    xValues = Subdivide[1, 2 xTerminal - 1, nmax];
    xValues = xValues // Rest;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Cold regime A > A_\[Natural]0*)


(* There are no critical terminal points; there is no terminal curve *)


(* Starting points along axis *)
startUVCold["axis"] =
  Module[{nmax, xValues},
    nmax = 16;
    xValues = Subdivide[1, 2, nmax] // Rest;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* Starting points along axis, left of the positive singularity *)
startUVCold["left"] =
  Module[{nmax, xValues},
    nmax = 16;
    xValues = Subdivide[0.6, 1, nmax] // Most;
    Table[UVBipolar[x, 0], {x, xValues}]
  ];


(* Starting points along a v-contour *)
startUVCold["contour"] =
  Module[{uValues, v},
    uValues = Subdivide[0, 2 Pi, 16] // Most;
    v = 3;
    Table[{u, v}, {u, uValues}]
  ];


(* ::Subsection:: *)
(*Traced boundary v = v(u)*)


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aHot;
    Table[
      uvInitList = startUVHot[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraHot[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "outer", "inner", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aWarmHot;
    Table[
      uvInitList = startUVWarmHot[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraWarmHot[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "outer", "inner", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aWarm;
    Table[
      uvInitList = startUVWarm[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraWarm[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "terminal", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aColdWarm;
    Table[
      uvInitList = startUVColdWarm[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraColdWarm[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Cold regime A > A_\[Natural]0*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aCold;
    Table[
      uvInitList = startUVCold[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraCold[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 4 Pi, uInit + 4 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"axis", "left", "contour"}}]
  ];
];


(* ::Subsection:: *)
(*Numerical verification for candidate boundaries (finite elements)*)


(* ::Text:: *)
(*These are not slow, nevertheless compute once and store.*)
(*Delete the corresponding file manually to compute from scratch.*)


(* ::Subsubsection:: *)
(*System of ODES for tracing by arc length parametrisation*)
(*u = u(s), v = v(s)*)


(* ::Text:: *)
(*See (48) & (49) (Page 12) and (r3.34) to (r3.39) (Page 3-5).*)


uvTraSystem[a_] :=
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{h, p, q, f, phi, grad2, uDer, vDer},
      (* Scale factor for bipolar coordinates *)
      h = HBipolar[u, v];
      (* Properties of known solution T == v *)
      p = 0;
      q = 1 / h;
      f = -v^4 / a;
      phi = vi[a][u, v];
      grad2 = p^2 + q^2;
      (* Return system of ODEs *)
      uDer = (-q f + p Re @ Sqrt[phi]) / (h grad2) // FullSimplify;
      vDer = (+p f + q Re @ Sqrt[phi]) / (h grad2);
      {u' == uDer, v' == vDer} /. {
        u' -> u'[s],
        v' -> v'[s],
        u -> u[s],
        v -> v[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Generate mesh for nice case*)


aNice = aInfl // Floor[#, 1/100] &


Module[{dest},
  dest = "bipolar-nice-verification-mesh.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {a, vSh0,
      sMax, uv,
      vBath,
      radSh0, sSpacing,
      sValues, reflect, extPointList,
      cenBath, radBath, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      (* Dimensionless group *)
      a = aNice;
      (* Hyperbolic critical terminal point *)
      vSh0 = vSharp0[a];
      (* Upper bound for arc length traversed *)
      sMax = Pi (1 - XBipolar[-Pi, vTraCandCorn[a]]);
      (* Traced boundaries *)
      uv = With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
        NDSolveValue[
          {
            uvTraSystem[a],
            u[0] == 0, v[0] == vSh0,
            WhenEvent[u[s] == -Pi, "StopIntegration"]
          }, {u, v}, {s, -sMax, 0},
          NoExtrapolation
        ]
      ];
      (* Actual arc length traversed *)
      sMax = -DomainStart[uv];
      (* Heat bath coordinate *)
      vBath = 1.1 vSh0;
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 7.5 degrees along v == vSh0) *)
      sSpacing = Csch[vSh0] * 7.5 Degree;
      (* External (radiation) boundary points *)
      sValues = UniformRange[-sMax, 0, sSpacing];
      reflect[{x_, y_}] := {x, -y};
      extPointList = Join[
        Table[
          XYBipolar @@ Through @ uv[s] // reflect
        , {s, sValues}] // Reverse // Most,
        Table[
          XYBipolar @@ Through @ uv[s]
        , {s, sValues}] // Most
      ];
      (* Internal (heat bath) boundary *)
      cenBath = {Coth[vBath], 0};
      radBath = Csch[vBath];
      intPointList = Table[
        cenBath + XYPolar[radBath, ph]
      , {ph, UniformRange[0, 2 Pi, sSpacing / radBath]}] // Most;
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
        "RegionHoles" -> {1, 0}
      ];
      (* Predicate functions for exterior and interior boundaries *)
      prExt = Function[{x, y}, VBipolar[x, y] < Way[vSh0, vBath] // Evaluate];
      prInt = Function[{x, y}, VBipolar[x, y] > Way[vSh0, vBath] // Evaluate];
      {a, vBath, mesh, prExt, prInt}
        // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Generate mesh for hot regime A < A_\[Natural]\[Pi]*)


(* ::Text:: *)
(*See "bipolar-hot-traced-convex-with-known.png".*)


Module[{dest},
  dest = "bipolar-hot-verification-mesh.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {a, vSh0,
      sMax, uv,
      vBath,
      radSh0, sSpacing,
      sValues, reflect, extPointList,
      cenBath, radBath, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      (* Dimensionless group *)
      a = aHot;
      (* Hyperbolic critical terminal point *)
      vSh0 = vSharp0[a];
      (* Upper bound for arc length traversed *)
      sMax = Pi (1 - XBipolar[-Pi, vTraCandCorn[a]]);
      (* Traced boundaries *)
      uv = With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
        NDSolveValue[
          {
            uvTraSystem[a],
            u[0] == 0, v[0] == vSh0,
            WhenEvent[u[s] == -Pi, "StopIntegration"]
          }, {u, v}, {s, -sMax, 0},
          NoExtrapolation
        ]
      ];
      (* Actual arc length traversed *)
      sMax = -DomainStart[uv];
      (* Heat bath coordinate *)
      vBath = 1.2 vSh0;
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 5 degrees along v == vSh0) *)
      sSpacing = Csch[vSh0] * 5 Degree;
      (* External (radiation) boundary points *)
      sValues = UniformRange[-sMax, 0, sSpacing];
      reflect[{x_, y_}] := {x, -y};
      extPointList = Join[
        Table[
          XYBipolar @@ Through @ uv[s] // reflect
        , {s, sValues}] // Reverse // Most,
        Table[
          XYBipolar @@ Through @ uv[s]
        , {s, sValues}] // Most
      ];
      (* Internal (heat bath) boundary *)
      cenBath = {Coth[vBath], 0};
      radBath = Csch[vBath];
      intPointList = Table[
        cenBath + XYPolar[radBath, ph]
      , {ph, UniformRange[0, 2 Pi, sSpacing / radBath]}] // Most;
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
        "RegionHoles" -> {1, 0}
      ];
      (* Predicate functions for exterior and interior boundaries *)
      prExt = Function[{x, y}, VBipolar[x, y] < Way[vSh0, vBath] // Evaluate];
      prInt = Function[{x, y}, VBipolar[x, y] > Way[vSh0, vBath] // Evaluate];
      {a, vBath, mesh, prExt, prInt}
        // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Generate mesh for warm-to-hot transition A = A_\[Natural]\[Pi]*)


(* ::Text:: *)
(*See "bipolar-warm_hot-traced-convex-with-known.png".*)


Module[{dest},
  dest = "bipolar-warm_hot-verification-mesh.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {a, vSh0,
      sMax, uv,
      vBath,
      radSh0, sSpacing,
      sValues, reflect, extPointList,
      cenBath, radBath, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      (* Dimensionless group *)
      a = aWarmHot;
      (* Hyperbolic critical terminal point *)
      vSh0 = vSharp0[a];
      (* Upper bound for arc length traversed *)
      sMax = Pi (1 - XBipolar[-Pi, vTraCandCorn[a]]);
      (* Traced boundaries *)
      uv = With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
        NDSolveValue[
          {
            uvTraSystem[a],
            u[0] == 0, v[0] == vSh0,
            WhenEvent[u[s] == -Pi, "StopIntegration"]
          }, {u, v}, {s, -sMax, 0},
          NoExtrapolation
        ]
      ];
      (* Actual arc length traversed *)
      sMax = -DomainStart[uv];
      (* Heat bath coordinate *)
      vBath = 1.2 vSh0;
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 5 degrees along v == vSh0) *)
      sSpacing = Csch[vSh0] * 5 Degree;
      (* External (radiation) boundary points *)
      sValues = UniformRange[-sMax, 0, sSpacing];
      reflect[{x_, y_}] := {x, -y};
      extPointList = Join[
        Table[
          XYBipolar @@ Through @ uv[s] // reflect
        , {s, sValues}] // Reverse // Most,
        Table[
          XYBipolar @@ Through @ uv[s]
        , {s, sValues}] // Most
      ];
      (* Internal (heat bath) boundary *)
      cenBath = {Coth[vBath], 0};
      radBath = Csch[vBath];
      intPointList = Table[
        cenBath + XYPolar[radBath, ph]
      , {ph, UniformRange[0, 2 Pi, sSpacing / radBath]}] // Most;
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
        "RegionHoles" -> {1, 0}
      ];
      (* Predicate functions for exterior and interior boundaries *)
      prExt = Function[{x, y}, VBipolar[x, y] < Way[vSh0, vBath] // Evaluate];
      prInt = Function[{x, y}, VBipolar[x, y] > Way[vSh0, vBath] // Evaluate];
      {a, vBath, mesh, prExt, prInt}
        // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Generate mesh for warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


(* ::Text:: *)
(*See "bipolar-warm-traced-convex-with-known.png".*)


Module[{dest},
  dest = "bipolar-warm-verification-mesh.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {a, vSh0,
      sMax, uv,
      vBath,
      radSh0, sSpacing,
      sValues, reflect, extPointList,
      cenBath, radBath, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      (* Dimensionless group *)
      a = aWarm;
      (* Hyperbolic critical terminal point *)
      vSh0 = vSharp0[a];
      (* Upper bound for arc length traversed *)
      sMax = Pi (1 - XBipolar[-Pi, vTraCandCorn[a]]);
      (* Traced boundaries *)
      uv = With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
        NDSolveValue[
          {
            uvTraSystem[a],
            u[0] == 0, v[0] == vSh0,
            WhenEvent[u[s] == -Pi, "StopIntegration"]
          }, {u, v}, {s, -sMax, 0},
          NoExtrapolation
        ]
      ];
      (* Actual arc length traversed *)
      sMax = -DomainStart[uv];
      (* Heat bath coordinate *)
      vBath = 1.2 vSh0;
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 5 degrees along v == vSh0) *)
      sSpacing = Csch[vSh0] * 5 Degree;
      (* External (radiation) boundary points *)
      sValues = UniformRange[-sMax, 0, sSpacing];
      reflect[{x_, y_}] := {x, -y};
      extPointList = Join[
        Table[
          XYBipolar @@ Through @ uv[s] // reflect
        , {s, sValues}] // Reverse // Most,
        Table[
          XYBipolar @@ Through @ uv[s]
        , {s, sValues}] // Most
      ];
      (* Internal (heat bath) boundary *)
      cenBath = {Coth[vBath], 0};
      radBath = Csch[vBath];
      intPointList = Table[
        cenBath + XYPolar[radBath, ph]
      , {ph, UniformRange[0, 2 Pi, sSpacing / radBath]}] // Most;
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
        "RegionHoles" -> {1, 0}
      ];
      (* Predicate functions for exterior and interior boundaries *)
      prExt = Function[{x, y}, VBipolar[x, y] < Way[vSh0, vBath] // Evaluate];
      prInt = Function[{x, y}, VBipolar[x, y] > Way[vSh0, vBath] // Evaluate];
      {a, vBath, mesh, prExt, prInt}
        // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Generate mesh for cold-to-warm transition A = A_\[Natural]0*)


(* ::Text:: *)
(*See "bipolar-cold_warm-traced-pseudo_convex-with-known.png".*)
(*NOTE: this boundary is actually slightly concave near the corner at u = -\[Pi].*)


Module[{dest},
  dest = "bipolar-cold_warm-verification-mesh.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {a, vSh0,
      sMax, uv,
      vBath,
      radSh0, sSpacing,
      sValues, reflect, extPointList,
      cenBath, radBath, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      (* Dimensionless group *)
      a = aColdWarm;
      (* Hyperbolic critical terminal point *)
      vSh0 = vSharp0[a];
      (* Upper bound for arc length traversed *)
      sMax = Pi (1 - XBipolar[-Pi, vTraCandCorn[a]]);
      (* Traced boundaries *)
      uv = With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
        NDSolveValue[
          {
            uvTraSystem[a],
            u[0] == 0, v[0] == vSh0,
            WhenEvent[u[s] == -Pi, "StopIntegration"]
          }, {u, v}, {s, -sMax, 0},
          NoExtrapolation
        ]
      ];
      (* Actual arc length traversed *)
      sMax = -DomainStart[uv];
      (* Heat bath coordinate *)
      vBath = 1.2 vSh0;
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 5 degrees along v == vSh0) *)
      sSpacing = Csch[vSh0] * 5 Degree;
      (* External (radiation) boundary points *)
      sValues = UniformRange[-sMax, 0, sSpacing];
      reflect[{x_, y_}] := {x, -y};
      extPointList = Join[
        Table[
          XYBipolar @@ Through @ uv[s] // reflect
        , {s, sValues}] // Reverse // Most,
        Table[
          XYBipolar @@ Through @ uv[s]
        , {s, sValues}] // Most
      ];
      (* Internal (heat bath) boundary *)
      cenBath = {Coth[vBath], 0};
      radBath = Csch[vBath];
      intPointList = Table[
        cenBath + XYPolar[radBath, ph]
      , {ph, UniformRange[0, 2 Pi, sSpacing / radBath]}] // Most;
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
        "RegionHoles" -> {1, 0}
      ];
      (* Predicate functions for exterior and interior boundaries *)
      prExt = Function[{x, y}, VBipolar[x, y] < Way[vSh0, vBath] // Evaluate];
      prInt = Function[{x, y}, VBipolar[x, y] > Way[vSh0, vBath] // Evaluate];
      {a, vBath, mesh, prExt, prInt}
        // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Solve BVP for nice case*)


Module[{dest},
  dest = "bipolar-nice-verification-solution.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, vBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "bipolar-nice-verification-mesh.txt";
      {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* External (radiation) boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prExt[x, y]],
            (* Internal (heat bath) boundary condition *)
            DirichletCondition[t[x, y] == vBath, prInt[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      ];
      tSol // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Solve BVP for hot regime A < A_\[Natural]\[Pi]*)


Module[{dest},
  dest = "bipolar-hot-verification-solution.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, vBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "bipolar-hot-verification-mesh.txt";
      {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* External (radiation) boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prExt[x, y]],
            (* Internal (heat bath) boundary condition *)
            DirichletCondition[t[x, y] == vBath, prInt[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      ];
      tSol // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Solve BVP for warm-to-hot transition A = A_\[Natural]\[Pi]*)


Module[{dest},
  dest = "bipolar-warm_hot-verification-solution.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, vBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "bipolar-warm_hot-verification-mesh.txt";
      {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* External (radiation) boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prExt[x, y]],
            (* Internal (heat bath) boundary condition *)
            DirichletCondition[t[x, y] == vBath, prInt[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      ];
      tSol // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Solve BVP for warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


Module[{dest},
  dest = "bipolar-warm-verification-solution.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, vBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "bipolar-warm-verification-mesh.txt";
      {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* External (radiation) boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prExt[x, y]],
            (* Internal (heat bath) boundary condition *)
            DirichletCondition[t[x, y] == vBath, prInt[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      ];
      tSol // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Solve BVP for cold-to-warm transition A = A_\[Natural]0*)


Module[{dest},
  dest = "bipolar-cold_warm-verification-solution.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, vBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "bipolar-cold_warm-verification-mesh.txt";
      {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* External (radiation) boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prExt[x, y]],
            (* Internal (heat bath) boundary condition *)
            DirichletCondition[t[x, y] == vBath, prInt[x, y]]
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
uIt = Italicise["u"];
vIt = Italicise["v"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Bipolar coordinates*)


(* ::Subsubsection:: *)
(*Effectively v = Infinity*)


vInfinity = 10;


(* ::Subsubsection:: *)
(*Contours*)


defContStyle = LightGray;


uContour[u_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {v, -vInfinity, vInfinity},
    Axes -> None,
    Exclusions -> If[Mod[u, 2 Pi] == 0, 0, Automatic],
    PlotStyle -> style
  ];


vContour[v_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {u, 0, 2 Pi},
    Axes -> None,
    PlotStyle -> style
  ];


uContours = uContour /@ Most[Range[0, 360, 30] Degree];
vContours = vContour /@ {-2, -1, -1/2, 0, 1/2, 1, 2};


(* ::Subsubsection:: *)
(*Local orthonormal basis vectors*)


uVector[u_, v_, style_] :=
  Graphics @ {Directive[Thick, style],
    Arrow @ {
      XYBipolar[u, v],
      XYBipolar[u, v] + AUBipolar[u, v]
    }
  };


vVector[u_, v_, style_] :=
  Graphics @ {Directive[Thick, style],
    Arrow @ {
      XYBipolar[u, v],
      XYBipolar[u, v] + AVBipolar[u, v]
    }
  };


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;
textStyleBigger = Style[#, 20] &;


natStyle = Red;


flat0Style = Darker[Green];
sharp0Style = Magenta;


flatPiStyle = Blue;
sharpPiStyle = Orange;


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;
inflStyle = Directive[Dashed, Green];


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


guideStyle = Dashed;
pointStyle = PointSize[Large];
glowStyle = Directive[Thick, Yellow, Opacity[0.7]];


uStyle = Blue;
vStyle = Red;


lighter = Directive[Opacity[0.5], #] &;


transStyle =
  Function[{a},
    If[a == aNat0 || a == aNatPi,
      Style[#, natStyle] &,
      Identity
    ]
  ];


(* ::Section:: *)
(*Bipolar coordinates*)


(* ::Subsection:: *)
(*Algebra*)


(* Check v == log(r_- / r_+) (Page r3-2) *)
With[{x = \[FormalX], y = \[FormalY]},
  Module[{rMinus, rPlus},
    rMinus = Sqrt[(x + 1)^2 + y^2];
    rPlus = Sqrt[(x - 1)^2 + y^2];
    VBipolar[x, y] == Log[rMinus / rPlus]
      // FullSimplify[#, Element[{x, y}, Reals]] &
  ]
]


(* ::Subsection:: *)
(*Interactive visualiser*)


DynamicModule[
 {uInit, vInit, xyInit,
  vMax, xMax, yMax
 },
  uInit = 90 Degree;
  vInit = 3/10;
  xyInit = XYBipolar[uInit, vInit];
  vMax = 5;
  xMax = 4;
  yMax = 3;
  Manipulate[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
    (* NOTE: In Version 11, TrackingFunction won't work *)
    (* unless endpoint 2 Pi is numericised *)
    (* https://mathematica.stackexchange.com/a/101695 *)
  , {{u, uInit, uIt}, 0, 2 Pi,
      Appearance -> "Open",
      TrackingFunction -> Function[{uNew},
        u = uNew;
        xy = XYBipolar[u, v];
      ]
  }
  , {{v, vInit, vIt}, -vMax, vMax,
      Appearance -> "Open",
      TrackingFunction -> Function[{vNew},
        v = vNew;
        xy = XYBipolar[u, v];
      ]
  }
  , {{xy, xyInit}, Locator,
      TrackingFunction -> Function[{xyNew},
        xy = xyNew;
        u = UBipolar @@ xy // Mod[#, 2 Pi] &;
        v = VBipolar @@ xy;
      ]
  }]
]


(* ::Subsection:: *)
(*Animation (u varying, v fixed)*)


Module[
 {uValues, v,
  xMax, yMax
 },
  uValues = Subdivide[360, 36] Degree // Most;
  v = 0.3;
  xMax = 4;
  yMax = 3;
  Table[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
  , {u, uValues}]
] // Ex["bipolar-sweep-u.gif", gifOpts]


(* ::Subsection:: *)
(*Animation (u fixed, v varying)*)


Module[
 {u, vValues,
  xMax, yMax
 },
  u = 90 Degree;
  vValues = Subdivide[-4, 4, 20];
  xMax = 4;
  yMax = 3;
  Table[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
  , {v, vValues}]
] // Ex["bipolar-sweep-v.gif", gifOpts]


(* ::Section:: *)
(*Critical terminal points*)


Module[
 {aCold, aWarm, aHot,
  vMin, vMax,
  regList,
  aAss,
  a, nString,
  pointAlongCurve, lineToAxis, textBelowAxis
 },
  (* Non-borderline representative values of A *)
  aCold = Ceiling[1.1 aNat0, 1/2];
  aWarm = Way[aNatPi, aNat0];
  aHot = Floor[0.95 aNatPi, 1/2];
  (* Horizontal plot range *)
  vMin = Floor[vFlat0[aHot], 1/2];
  vMax = Ceiling[vSharp0[aHot], 1/2];
  (* List of regimes *)
  regList = {"hot", "warm_hot", "warm", "cold_warm", "cold"};
  (* Associations for various quantities *)
  aAss = AssociationThread[regList -> {aHot, aNatPi, aWarm, aNat0, aCold}];
  (* Plots for all five regimes *)
  Table[
    a = aAss[reg];
    nString = ToString @ First @ FirstPosition[regList, reg];
    pointAlongCurve = Point @ {#, #^4 / a} &;
    lineToAxis = Line @ {{#, #^4 / a}, {#, 0}} &;
    textBelowAxis = Function[{v, sub},
      Text[
        Subscript[
          vIt,
          sub // PrettyString[
            "Sharp" -> "\[Sharp]",
            "Nat" -> "\[Natural]",
            "Flat" -> "\[Flat]",
            "Pi" -> "\[Pi]"
          ]
        ] // textStyle,
        {v, 0},
        {0, 2.5}
      ]
    ];
    Show[
      (* Functions of v *)
      Plot[{Cosh[v] + 1, Cosh[v] - 1, v^4 / a}, {v, vMin, vMax},
        AxesLabel -> {vIt, Null},
        ImageSize -> 640,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], reg // PrettyString["_" -> "-to-"] // ""},
            "  "
          ] // transStyle[a]
        ],
        PlotLegends -> {
          Row[{Cosh[vIt] + 1, uIt == Pi // ""}, "  "],
          Row[{Cosh[vIt] - 1, uIt == 0 // ""}, "  "],
          Row[{vIt^4, aIt}, "/"]
        },
        PlotRange -> All,
        PlotRangeClipping -> False,
        PlotStyle -> {lighter[Blue], lighter[Red], Black},
        PlotOptions[Axes] // Evaluate
      ],
      (* Critical terminal points along u == 0 *)
      Which[
        (* Two distinct terminal points *)
        a < aNat0,
        {
          (* v_\[Flat]0 *)
          Graphics @ {Directive[flat0Style, pointStyle],
            pointAlongCurve @ vFlat0[a]
          },
          Graphics @ {Directive[flat0Style, guideStyle],
            lineToAxis @ vFlat0[a]
          },
          Graphics @ {flat0Style,
            textBelowAxis[vFlat0[a], "Flat0"]
          },
          (* v_\[Sharp]0 *)
          Graphics @ {Directive[sharp0Style, pointStyle],
            pointAlongCurve @ vSharp0[a]
          },
          Graphics @ {Directive[sharp0Style, guideStyle],
            lineToAxis @ vSharp0[a]
          },
          Graphics @ {sharp0Style,
            textBelowAxis[vSharp0[a], "Sharp0"]
          }
        },
        (* One terminal point *)
        a == aNat0,
        {
          (* v_\[Natural]0 *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongCurve @ vNat0
          },
          Graphics @ {Directive[natStyle, guideStyle],
            lineToAxis @ vNat0
          },
          Graphics @ {natStyle,
            textBelowAxis[vNat0, "Nat0"]
          }
        },
        (* Zero terminal points *)
        True,
        {}
      ],
      (* Critical terminal points along u == \[Pi] *)
      Which[
        (* Two distinct terminal points *)
        a < aNatPi,
        {
          (* v_\[Flat]\[Pi] *)
          Graphics @ {Directive[flatPiStyle, pointStyle],
            pointAlongCurve @ vFlatPi[a]
          },
          Graphics @ {Directive[flatPiStyle, guideStyle],
            lineToAxis @ vFlatPi[a]
          },
          Graphics @ {flatPiStyle,
            textBelowAxis[vFlatPi[a], "FlatPi"]
          },
          (* v_\[Sharp]\[Pi] *)
          Graphics @ {Directive[sharpPiStyle, pointStyle],
            pointAlongCurve @ vSharpPi[a]
          },
          Graphics @ {Directive[sharpPiStyle, guideStyle],
            lineToAxis @ vSharpPi[a]
          },
          Graphics @ {sharpPiStyle,
            textBelowAxis[vSharpPi[a], "SharpPi"]
          }
        },
        (* One terminal point *)
        a == aNatPi,
        {
          (* v_\[Natural]\[Pi] *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongCurve @ vNatPi
          },
          Graphics @ {Directive[natStyle, guideStyle],
            lineToAxis @ vNatPi
          },
          Graphics @ {natStyle,
            textBelowAxis[vNatPi, "NatPi"]
          }
        },
        (* Zero terminal points *)
        True,
        {
          Graphics @ {Transparent (* for consistent plot height *),
            textBelowAxis[vNatPi, "NatPi"]
          }
        }
      ]
    ] // Ex @ StringJoin["bipolar-critical-", nString, "-", reg, ".pdf"]
  , {reg, regList}]
]


(* ::Section:: *)
(*Viable domain animations*)


Module[
 {aMin, aMax, aStep,
  aValues,
  xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  (* Values of A *)
  aMin = 1;
  aMax = 10;
  aStep = 1/4;
  aValues = Range[aMin, aMax, aStep];
  aValues = Join[aValues, {aNat0, aNatPi}] // Sort[#, Less] &;
  (* Animation *)
  Table[
    (* Plot range *)
    xMax = 3.5;
    yMax = 2.5;
    (* Plot range for unphysical domain *)
    eps = 0.1;
    xMaxUnphys = xMax (1 + eps);
    yMaxUnphys = yMax (1 + eps);
    (* Plot range for non-viable domain *)
    (* (these bounds are merely guesses which work) *)
    xMinNon = If[a < aNatPi, vFlatPi[a], vNatPi] // 0.9 XBipolar[Pi, #] &;
    xMaxNon = If[a < aNat0, vFlat0[a], vNat0] // 1.01 XBipolar[0, #] &;
    yMaxNon = (xMaxNon - xMinNon) / 2;
    (* Plot *)
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[x < 0,
        {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Non-viable domain *)
      RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ]
    ]
  , {a, aValues}]
] // Ex["bipolar-viable-full.gif", gifOpts]


Module[
 {aMin, aMax, aStep,
  aValues,
  pointAlongAxis, textAboveBelowAxis,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  (* Values of A *)
  aMin = 8;
  aMax = 10;
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  aValues = Join[aValues, {aNat0, aNatPi}] // Sort[#, Less] &;
  (* Graphics abbreviations *)
  pointAlongAxis[u_, v_] := Point @ {XYBipolar[u, v]};
  textAboveBelowAxis[u_, v_, sub_, ysign_: 1] :=
    Text[
      Subscript[
        vIt,
        sub // PrettyString[
          "Sharp" -> "\[Sharp]",
          "Nat" -> "\[Natural]",
          "Flat" -> "\[Flat]",
          "Pi" -> "\[Pi]"
        ]
      ] // textStyle // BoxedLabel[#, Background -> White] &,
      XYBipolar[u, v],
      {0, -Sign[ysign] 1.7}
    ];
  (* Animation *)
  Table[
    (* Plot range *)
    xMin = vFlatPi[aMin] // 0.99 XBipolar[Pi, #] &;
    xMax = vFlat0[aMin] // 1.01 XBipolar[0, #] &;
    yMax = (xMax - xMin) / 2;
    (* Plot range for non-viable domain *)
    (* (these bounds are merely guesses which work) *)
    xMinNon = If[a < aNatPi, vFlatPi[a], vNatPi] // 0.9 XBipolar[Pi, #] &;
    xMaxNon = If[a < aNat0, vFlat0[a], vNat0] // 1.01 XBipolar[0, #] &;
    yMaxNon = (xMaxNon - xMinNon) / 2;
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Non-viable domain *)
      RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Critical terminal points along u == 0 *)
      Which[
        (* Two distinct terminal points *)
        a < aNat0,
        {
          (* v_\[Flat]0 *)
          Graphics @ {Directive[flat0Style, pointStyle],
            pointAlongAxis[0, vFlat0[a]]
          },
          Graphics @ {flat0Style,
            textAboveBelowAxis[0, vFlat0[a], "Flat0", -1]
          },
          (* v_\[Sharp]0 *)
          Graphics @ {Directive[sharp0Style, pointStyle],
            pointAlongAxis[0, vSharp0[a]]
          },
          Graphics @ {sharp0Style,
            textAboveBelowAxis[0, vSharp0[a], "Sharp0"]
          }
        },
        (* One terminal point *)
        a == aNat0,
        {
          (* v_\[Natural]0 *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongAxis[0, vNat0]
          },
          Graphics @ {natStyle,
            textAboveBelowAxis[0, vNat0, "Nat0"]
          }
        },
        (* Zero terminal points *)
        True,
        {}
      ],
      (* Critical terminal points along u == \[Pi] *)
      Which[
        (* Two distinct terminal points *)
        a < aNatPi,
        {
          (* v_\[Flat]\[Pi] *)
          Graphics @ {Directive[flatPiStyle, pointStyle],
            pointAlongAxis[Pi, vFlatPi[a]]
          },
          Graphics @ {flatPiStyle,
            textAboveBelowAxis[Pi, vFlatPi[a], "FlatPi"]
          },
          (* v_\[Sharp]\[Pi] *)
          Graphics @ {Directive[sharpPiStyle, pointStyle],
            pointAlongAxis[Pi, vSharpPi[a]]
          },
          Graphics @ {sharpPiStyle,
            textAboveBelowAxis[Pi, vSharpPi[a], "SharpPi", -1]
          }
        },
        (* One terminal point *)
        a == aNatPi,
        {
          (* v_\[Natural]\[Pi] *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongAxis[Pi, vNatPi]
          },
          Graphics @ {natStyle,
            textAboveBelowAxis[Pi, vNatPi, "NatPi"]
          }
        },
        (* Zero terminal points *)
        True,
        {}
      ]
    ]
  , {a, aValues}]
] // Ex["bipolar-viable-zoom.gif", gifOpts]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  idList
 },
  a = aHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0.99 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.07 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vFlatPi[a] // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Group names *)
  idList = {"outer", "inner", "axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVHot[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 4, Left, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-hot-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot range for non-viable domain *)
  xMinNon = vFlatPi[a] // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-hot-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0.95 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.05 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vFlatPi[a] // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-hot-traced-zoom_1.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aHot;
  (* Plot range *)
  xMin = vSharpPi[a] // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vSharpPi[a] // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraHot[id]}]
    , {id, {"inner"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-hot-traced-zoom_2.pdf"]


(* ::Subsubsection:: *)
(*Convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v
 },
  a = aHot;
  (* Plot range *)
  xMin = vSharpPi[a] // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vSharpPi[a] // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-hot-traced-convex.pdf"]


(* ::Subsubsection:: *)
(*Convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v, vBath
 },
  a = aHot;
  (* Plot range *)
  xMin = vSharpPi[a] // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vSharpPi[a] // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-hot-traced-convex-with-known.png"]


(* ::Subsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  idList
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0.99 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.07 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Group names *)
  idList = {"outer", "inner", "axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVWarmHot[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 4, Left, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-warm_hot-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm_hot-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.95 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.05 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm_hot-traced-zoom_1.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"outer", "inner"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm_hot-traced-zoom_2.pdf"]


(* ::Subsubsection:: *)
(*Convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarmHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm_hot-traced-convex.pdf"]


(* ::Subsubsection:: *)
(*Convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v, vBath
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarmHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm_hot-traced-convex-with-known.png"]


(* ::Subsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  idList
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.07 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Group names *)
  idList = {"terminal", "axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVWarm[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 3, Left, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-warm-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarm[id]}]
    , {id, {"terminal", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.97 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.03 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarm[id]}]
    , {id, {"terminal", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm-traced-zoom.pdf"]


(* ::Subsubsection:: *)
(*Convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.97 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.1 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm-traced-convex.pdf"]


(* ::Subsubsection:: *)
(*Convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v, vBath
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.97 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.1 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm-traced-convex-with-known.png"]


(* ::Subsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  idList
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 0.99;
  xMax = vNat0 // 2.01 XBipolar[0, #] - 1 &;
  yMax = (xMax - xMin) / 2;
  (* Group names *)
  idList = {"axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVColdWarm[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 2, Right, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-cold_warm-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vNat0 // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-cold_warm-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 9/8;
  yMax = (xMax - xMin) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-cold_warm-traced-zoom.pdf"]


(* ::Subsubsection:: *)
(*Pseudo-convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  v
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 17/16;
  yMax = (xMax - xMin) / 2;
  (* Traced boundaries *)
  v = vTraColdWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-cold_warm-traced-pseudo_convex.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  v, uStart,
  reflect, tanVec,
  uInfl, vInfl, xyInfl,
  tanVecInfl, tanLineInfl,
  xyCorn, xyInt,
  mark,
  newLine = "\n",
  concaveStyle = Red,
  inflDotStyle = Directive[Yellow, Opacity[0.5], pointStyle],
  inflMarkStyle = Directive[Thin, Darker[Green]],
  tanStyle = Directive[Blue, Thin],
  cornDotStyle = Directive[Black, Opacity[0.7], pointStyle],
  cornMarkStyle = Directive[Red, Thin],
  intDotStyle = Directive[Cyan, Opacity[0.5], pointStyle],
  intMarkStyle = Directive[Blue, Opacity[0.7], Thin]
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 17/16;
  yMax = (xMax - xMin) / 2;
  (* Traced boundaries (NOTE: uStart is NEGATIVE) *)
  v = vTraColdWarm["hyperbolic"] // First;
  uStart = DomainStart[v];
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Reflection across x-axis *)
  reflect[{x_, y_}] := {x, -y};
  (* Cartesian components of tangent vector to traced boundary *)
  tanVec[u_, v_] :=
    HBipolar[u, v] Plus[
      AUBipolar[u, v],
      vTraDer[a][u, v] AVBipolar[u, v]
    ] // Evaluate;
  (* Points of inflection (NOTE: uInfl is NEGATIVE) *)
  uInfl = SeekRoot[curTra[a][#, v[#]] &, {-Pi, 0}];
  vInfl = v[uInfl];
  xyInfl = XYBipolar[uInfl, vInfl];
  tanVecInfl = tanVec[uInfl, vInfl];
  tanLineInfl = InfiniteLine[xyInfl, tanVecInfl];
  (* Corner formed by concave portions *)
  xyCorn = XYBipolar[-Pi, v[-Pi]];
  (* Intersection of tangent lines *)
  xyInt =
    With[{x = \[FormalX], y = \[FormalY]},
      {x, y} /. First @ FindInstance[
        And[
          Element[{x, y}, tanLineInfl],
          Element[{x, y}, reflect /@ tanLineInfl]
        ]
      , {x, y}]
    ];
  (* Marker *)
  mark[{x_, y_}, rot_, len_: 0.01] := (
    Line @ {{0, -len}, {0, len}}
      // Rotate[#, rot] &
      // Translate[#, {x, y}] &
  );
  (* Plots *)
  {
    (* Traced boundary convexity as a function of u *)
    Show[
      Plot[
        curTra[a][u, v[u]] {
          If[u < uInfl, 1, Indeterminate],
          If[u > uInfl, 1, Indeterminate]
        } // Evaluate, {u, -Pi, 0},
        AxesLabel -> {uIt, "Cur."},
        PlotRange -> Full,
        PlotStyle -> {concaveStyle, convexStyle},
        PlotOptions[Axes] // Evaluate
      ],
      (* Inflection *)
      Graphics @ {inflDotStyle,
        Point @ {uInfl, 0}
      },
      Graphics @ Text[
        Subscript[uIt, "i"] == uInfl // textStyle,
        {uInfl, 0},
        {-1.1, -1.5}
      ]
    ] // Ex["bipolar-cold_warm-traced-pseudo_convex-curvature.pdf"],
    (* Detailed plot demonstrating inflection and concavity *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Traced boundaries (convex portion) *)
      ParametricPlot[
        XYBipolar[u, v[u]],
        {u, uInfl, -uInfl},
        PlotStyle -> convexStyle
      ],
      (* Traced boundaries (concave portions) *)
      ParametricPlot[
        XYBipolar[u, v[u]],
        {u, uStart, uInfl},
        PlotStyle -> concaveStyle
      ],
      ParametricPlot[
        XYBipolar[u, v[u]],
        {u, -uInfl, -uStart},
        PlotStyle -> concaveStyle
      ],
      (* Points of inflection *)
      Graphics @ {inflDotStyle,
        Point @ {xyInfl, xyInfl // reflect}
      },
      (* Mark at points of inflection *)
      Graphics @ {inflMarkStyle,
        mark[xyInfl, {{1, 0}, tanVecInfl}],
        mark[xyInfl // reflect, {{1, 0}, tanVecInfl // reflect}]
      },
      (* Tangent line at points of inflection *)
      Graphics @ {tanStyle,
        tanLineInfl,
        reflect /@ tanLineInfl
      },
      (* Corner formed by concave portions *)
      Graphics @ {cornDotStyle,
        Point @ {xyCorn}
      },
      (* Marker at corner *)
      Graphics @ {cornMarkStyle,
        mark[xyCorn, 0]
      },
      (* Intersection of tangent lines *)
      Graphics @ {intDotStyle,
        Point @ {xyInt}
      },
      (* Mark at intersection of tangent lines *)
      Graphics @ {intMarkStyle,
        mark[xyInt, 0]
      },
      (* Explanatory text *)
      Graphics @ {
        Text[
          StringJoin[
            "Black: convex portion of boundary", newLine,
            "Red: concave portions of boundary", newLine,
            "Green marks (yellow dots): points of inflection", newLine,
            "Blue: tangent line at points of inflection", newLine,
            "Red mark (black dot): corner formed by concave portions", newLine,
            "Blue mark (cyan dot): intersection of tangent lines", newLine,
            newLine,
            "This PDF consists of vector graphics.", newLine,
            "Zoom in to see the small gap between the red and blue marks."
          ] // Style[#, 7, Background -> White] &,
          {Way[xMin, xMax], Way[0, yMax, 7/10]}
        ]
      }
    ] // Ex["bipolar-cold_warm-traced-pseudo_convex-detailed.pdf"]
  }
]


(* ::Subsubsection:: *)
(*Pseudo-convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  v, vBath
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 17/16;
  yMax = (xMax - xMin) / 2;
  (* Traced boundaries *)
  v = vTraColdWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-cold_warm-traced-pseudo_convex-with-known.png"]


(* ::Subsection:: *)
(*Cold regime A > A_\[Natural]0*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  idList
 },
  a = aCold;
  (* Plot range *)
  xMin = 0.6;
  xMax = 2;
  yMax = (xMax - xMin) / 2;
  (* Group names *)
  idList = {"axis", "left", "contour"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVCold[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 2, Top, Bottom]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-cold-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys
 },
  a = aCold;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vNat0 // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraCold[id]}]
    , {id, {"axis", "left", "contour"}}]
  ]
] // Ex["bipolar-cold-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax
 },
  a = aCold;
  (* Plot range *)
  xMin = 0.6;
  xMax = 1.25;
  yMax = (xMax - xMin) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraCold[id]}]
    , {id, {"axis", "left", "contour"}}]
  ]
] // Ex["bipolar-cold-traced-zoom.pdf"]


(* ::Section:: *)
(*Traced boundary algebra*)


(* ::Text:: *)
(*See Section 4.4 (Page r3-10).*)


(* ::Subsection:: *)
(*Check derivatives*)


(* ::Text:: *)
(*See (r3.56) to (r3.59)(Page r3-10).*)


With[{u = \[FormalU], v = \[FormalV], a = \[FormalCapitalA]},
  Module[
   {au, av, d, ss, cc, h, dd, ee
   },
    (* Local orthogonal basis *)
    au = AUBipolar[u, v];
    av = AVBipolar[u, v];
    (* Abbreviation for u-derivative *)
    d = Dt[#, u, Constants -> a] &;
    (* Abbreviations *)
    ss = Sin[u] Sinh[v];
    cc = Cos[u] Cosh[v] - 1;
    h = HBipolar[u, v];
    dd = Sinh[v] - d[v] Sin[u];
    ee = d[v] Sinh[v] + Sin[u];
    (* Table of results *)
    {
      {
        "(h S)' == h^2 C D",
        "(r3.56)",
        d[h ss] == h^2 cc dd // FullSimplify
      },
      {
        "(h C)' == -h^2 S D",
        "(r3.57)",
        d[h cc] == -h^2 ss dd // FullSimplify
      },
      {
        "(a_u)' == h D a_v",
        "(r3.58)",
        d[au] == h dd av // FullSimplify
      },
      {
        "(a_v)' == -h D a_u",
        "(r3.59)",
        d[av] == -h dd au // FullSimplify
      }
    } // TableForm
  ]
] // Ex["bipolar-check-derivatives.pdf"]


(* ::Subsection:: *)
(*Check curvature*)


(* ::Text:: *)
(*See (r3.62) and (r3.63) (Page r3-10).*)


With[{u = \[FormalU], v = \[FormalV], a = \[FormalCapitalA]},
  Module[
   {au, av, d, ss, cc, h, dd, ee,
    dSimp, vel, acc, comp, cur,
    newLine = "\n  "
   },
    (* Local orthogonal basis *)
    {au, av};
    (* Abbreviation for u-derivative *)
    d = Dt[#, u, Constants -> a] &;
    (* Abbreviations *)
    ss = Sin[u] Sinh[v];
    cc = Cos[u] Cosh[v] - 1;
    h = HBipolar[u, v];
    dd = Sinh[v] - d[v] Sin[u];
    ee = d[v] Sinh[v] + Sin[u];
    (* Simplification rules *)
    dSimp = {
      d[au] -> h dd av,
      d[av] -> -h dd au
    };
    (* Velocity and acceleration *)
    vel = h au + h d[v] av;
    acc = d[vel] /. dSimp;
    (* Components for local (u, v, z) basis *)
    comp = {
      au -> {1, 0, 0},
      av -> {0, 1, 0}
    };
    (* Curvature *)
    cur = Last @ Cross[
      vel /. comp,
      acc /. comp
    ];
    (* Table of results *)
    {
      {
        StringJoin[
          "Acc. ==", newLine,
          "- h^2 (D v' + E) a_u", newLine,
          "+ h (v'' - h E v' + h D) a_v"
        ],
        "(r3.62)",
        acc ==
          - h^2 (dd d[v] + ee) au
          + h (d @ d[v] - h ee d[v] + h dd) av
          /. comp
          // FullSimplify
      },
      {
        "Cur. == h^3 (D (1 + v'^2) + v'' / h)",
        "(r3.63)",
        cur == h^3 (dd (1 + d[v]^2) + d @ d[v] / h)
          /. comp 
          // FullSimplify
      }
    } // TableForm
  ]
] // Ex["bipolar-check-curvature.pdf"]


(* ::Subsection:: *)
(*Curvature at u = -\[Pi]*)


(* ::Text:: *)
(*See (r3.64) to (r3.66) (Page r3-10).*)


With[{v = \[FormalV], a = \[FormalCapitalA]},
  {
    {
      "Cur. (u = -\[Pi])",
      "(r3.65)",
      curTra[a][-Pi, v] // FullSimplify
    },
    {
      "v_i",
      "(r3.66)",
      vInfl // N
    },
    {
      "",
      "(check)",
      curTra[a][-Pi, vInfl] == 0 // FullSimplify
    }
  } // TableForm
] // Ex["bipolar-curvature-corner.pdf"]


(* Verify same result for lower branch too *)
With[{v = \[FormalV], a = \[FormalCapitalA]},
  Module[{upper, lower, equal},
    upper = curTra[a][-Pi, v] // FullSimplify;
    lower = curTraLower[a][Pi, v] // FullSimplify;
    equal = upper == lower;
    {
      {"Upper", "Lower", "Equal"},
      {upper, lower, equal}
    } // Transpose // TableForm
  ]
] // Ex["bipolar-curvature-corner-branches.pdf"]


(* ::Subsection:: *)
(*Hyperbolic coordinate at u = -\[Pi] for inner candidate boundaries*)


Module[
 {aMin, aMax, num, aValues,
  dest, aVTable,
  inflDotStyle = Directive[Red, Opacity[0.7], pointStyle]
 },
  (* Values of dimensionless group A *)
  (* (argument to vTraCandCorn must be EXACT) *)
  aMin = aNat0 - 4/1000;
  aMax = aNat0;
  num = 100;
  aValues = Subdivide[aMin, aMax, num];
  (* Compute table of (A, v(-\[Pi])) values *)
  (* (This is really slow (~25 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "bipolar-candidate-v-corner.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    aVTable = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aVTable = Table[
      {a // N, vTraCandCorn[a]}
    , {a, aValues}];
    aVTable // Compress // Ex[dest]
  ];
  (* Plot *)
  Show[
    (* A vs v(-\[Pi]) *)
    ListPlot[
      aVTable,
      AxesLabel -> {aIt, vIt[-Pi]},
      Joined -> True,
      PlotRange -> All,
      PlotStyle -> Black,
      PlotOptions[Axes] // Evaluate
    ],
    (* v_i *)
    Plot[vInfl, {a, aMin, aMax},
      PlotStyle -> inflStyle
    ],
    (* A_i *)
    Graphics @ {inflDotStyle,
      Point @ {aInfl, vInfl}
    },
    Graphics @ Text[
      Column @ {
        Subscript[aIt, "i"] == N[aInfl],
        Subscript[vIt, "i"] == N[vInfl]
      } // textStyle,
      {aInfl, vInfl},
      {-1.2, -1.2}
    ]
  ]
] // Ex["bipolar-candidate-v-corner.pdf"]


(* ::Subsection:: *)
(*Hyperbolic coordinate at u = -\[Pi] for outer candidate boundaries*)


(* ::Subsubsection:: *)
(*Visualisation for both outer and inner*)


(* ::Text:: *)
(*Conclusion: outer candidate boundary is never convex.*)


DynamicModule[
  {
    aMin, aMax, aInit,
    xMin, xMax, yMax,
    aInflRationalised,
    aRationalised,
    vInnerInfl,
    vInner, vOuter,
    dummyForTrailingCommas
  },
  (* A values *)
  {aMin, aMax} = {0.9986, 1} aNat0;
  aInit = 0.999 aInfl // N;
  (* Plot range *)
  xMin = 0.5 XBipolar[0, vNat0];
  xMax = 1.03 XBipolar[0, vNat0];
  yMax = 0.1 xMax;
  (* Visualise *)
  Manipulate[
    (* Rationalise values of A *)
    aInflRationalised = aInfl // N // Rationalize[#, 0] &;
    aRationalised = a // N // Rationalize[#, 0] &;
    (* Compute inflection (borderline convex) inner boundary *)
    vInnerInfl =
      With[{v = \[FormalV]},
        NDSolveValue[
          {
            v'[u] == Re @ vTraDer[aInflRationalised][u, v[u]],
            v[0] == vSharp0[aInflRationalised]
          }, v, {u, -Pi, 0}
          , NoExtrapolation
          , PreciseOptions[]
        ]
      ];
    (* Compute inner and outer candidate boundaries *)
    vInner =
      With[{v = \[FormalV]},
        NDSolveValue[
          {
            v'[u] == Re @ vTraDer[aRationalised][u, v[u]],
            v[0] == vSharp0[aRationalised]
          }, v, {u, -Pi, 0}
          , NoExtrapolation
          , PreciseOptions[]
        ]
      ];
    vOuter =
      With[{v = \[FormalV]},
        NDSolveValue[
          {
            v'[u] == Re @ vTraDer[aRationalised][u, v[u]],
            v[0] == vFlat0[aRationalised]
          }, v, {u, -Pi, 0}
          , NoExtrapolation
          , PreciseOptions[]
        ]
      ];
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax}
        , ImageSize -> Full
        , PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]}
          , {u, DomainStart[v], DomainEnd[v]}
          , PlotStyle -> LightGray
        ]
        , {v, {vInnerInfl}}
      ],
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]}
          , {u, DomainStart[v], DomainEnd[v]}
          , PlotStyle -> {upperStyle, lowerStyle}
        ]
        , {v, {vOuter, vInner}}
      ],
      {}
    ]
    , {{a, aInit}, aMin, aMax}
  ]
]


(* ::Section:: *)
(*Physical temperature range algebra*)


(* ::Text:: *)
(*See Section 5.1 (Page r3-14).*)


(* ::Subsection:: *)
(*\[Omega] properties*)


With[{v = \[FormalV]},
  Module[{omegaAsyInfinity},
    omegaAsyInfinity[v_] := 1 / v - 2 Exp[-v] / v;
    {
      {"omega", omega[v]},
      {"Small v series", omega[v] + O[v]^3},
      {"Large v asymptotic", omegaAsyInfinity[v]},
      {
        "(check)",
        And @@ (
          Limit[#, v -> Infinity] == 1 & /@
          {
            (1 / v) / omega[v],
            (omega[v] - (1 / v)) / (- 2 Exp[-v] / v)
          }
        )
      },
      {"domega/dv", omega'[v] // FullSimplify}
    } /. {v -> Subscript[v, "sharp0"]}
      // PrettyString["omega" -> "\[Omega]", "v" -> ToString[v] <> "_sharp0"]
      // PrettyString["sharp" -> "\[Sharp]"]
      // TableForm
  ]
] // Ex["bipolar-omega-properties.pdf"]


(* ::Subsection:: *)
(*\[Omega] plot*)


Module[{omegaAsyOrigin, omegaAsyInfinity, vItSharp0},
  omegaAsyOrigin[v_] := 1 / 2 - v^2 / 24;
  omegaAsyInfinity[v_] := 1 / v - 2 Exp[-v] / v;
  vItSharp0 = Subscript[vIt, "sharp0"];
  Show[
    Plot[
      {
        omega[v],
        Piecewise @ {
          {omegaAsyOrigin[v], v < 3.4},
          {Indeterminate, True}
        },
        Piecewise @ {
          {omegaAsyInfinity[v], v > 0.7},
          {Indeterminate, True}
        }
      }, {v, 0, 10},
      AxesLabel -> {vItSharp0, "omega"},
      PlotLegends -> {
        "omega"[vItSharp0],
        omegaAsyOrigin[vItSharp0],
        omegaAsyInfinity[vItSharp0]
      },
      PlotRange -> All,
      PlotRangeClipping -> False,
      PlotStyle -> {
        Blue,
        Directive[guideStyle, Darker[Green]],
        Directive[guideStyle, Magenta]
      },
      PlotOptions[Axes] // Evaluate
    ],
    (* {vNat0, 1/4} *)
    Graphics @ {Directive[guideStyle, natStyle],
      Line @ {{0, 1/4}, {vNat0, 1/4}, {vNat0, 0}}
    },
    Graphics @ {natStyle,
      Text[Row[{1, 4}, "/"] // textStyle, {0, 1/4}, {1.3, 0}],
      Text[Subscript[vIt, "nat0"] // textStyle, {vNat0, 0}, {0, 2}]
    }
  ]
    // PrettyString[
      "omega" -> "\[Omega]",
      "nat" -> "\[Natural]",
      "sharp" -> "\[Sharp]"
    ]
] // Ex["bipolar-omega.pdf"]


(* ::Section:: *)
(*Physical asymmetry*)


(* ::Subsection:: *)
(*Check asymmetry*)


(* ::Text:: *)
(*See (r3.92) (Page r3-16).*)


With[{v = \[FormalV]},
  Module[{vSharp0, vCorn},
    vSharp0 = Subscript[v, "sharp0"];
    vCorn = v[-Pi];
    {
      {"Asymmetry", asymm[vSharp0, vCorn]},
      {
        "(r3.92)",
        Equal[
          asymm[vSharp0, vCorn],
          Divide[
            1 - Sinh[vCorn] / (Cosh[vCorn] + 1),
            Sinh[vSharp0] / (Cosh[vSharp0] - 1) - 1
          ]
        ] // FullSimplify
      }
    }
      // PrettyString["sharp" -> "\[Sharp]"]
      // TableForm
  ]
] // Ex["bipolar-check-asymmetry.pdf"]


(* ::Subsection:: *)
(*Plot asymmetry*)


Module[
 {aMin, aMax, num, aValues,
  dest, aAsymmTable, asymmNatPi, asymmNat0
 },
  (* Values of dimensionless group A *)
  (* (argument to vTraCandCorn must be EXACT) *)
  aMin = 8;
  aMax = aNat0;
  num = 100;
  aValues = Subdivide[aMin, aMax, num];
  (* Compute table of (A, asymmetry) values *)
  (* (This is really slow (~35 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "bipolar-candidate-asymmetry.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {aAsymmTable, asymmNatPi, asymmNat0} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aAsymmTable = Table[
      {a // N, asymm[vSharp0[a], vTraCandCorn[a]]}
    , {a, aValues}];
    asymmNatPi = asymm[vSharp0[aNatPi], vTraCandCorn[aNatPi]];
    asymmNat0 = Last[aAsymmTable][[2]];
    {aAsymmTable, asymmNatPi, asymmNat0} // Compress // Ex[dest]
  ];
  (* Plot *)
  Show[
    (* A vs asymmetry *)
    ListPlot[
      aAsymmTable,
      AxesLabel -> {aIt, "Asymm."},
      Joined -> True,
      PlotRange -> All,
      PlotStyle -> Black,
      PlotOptions[Axes] // Evaluate
    ],
    (* aNatPi *)
    Graphics @ {Directive[natStyle, guideStyle],
      Line @ {{aMin, asymmNatPi}, {aNatPi, asymmNatPi}, {aNatPi, 0}}
    },
    Graphics @ {natStyle,
      Text[N[asymmNatPi], {aMin, asymmNatPi}, {1.2, 0}],
      Text[Subscript[aIt, "natpi"] // textStyle, {aNatPi, 0}, {0, 2.3}]
    },
    (* aNat0 *)
    Graphics @ {Directive[natStyle, guideStyle],
      Line @ {{aMin, asymmNat0}, {aNat0, asymmNat0}, {aNat0, 0}}
    },
    Graphics @ {natStyle,
      Text[N[asymmNat0], {aMin, asymmNat0}, {1.2, 0}],
      Text[Subscript[aIt, "nat0"] // textStyle, {aNat0, 0}, {0, 2.3}]
    },
    (* Regime labels *)
    Graphics @ {
      Text["hot" // textStyle, {Way[aMin, aNatPi], 3/5}],
      Text["warm" // textStyle, {Way[aNatPi, aNat0], 3/4}]
    },
    (* Don't chop labels *)
    ImageMargins -> None,
    PlotRange -> All,
    PlotRangeClipping -> False
  ] // PrettyString["nat" -> "\[Natural]", "pi" -> "\[Pi]"]
] // Ex["bipolar-candidate-asymmetry.pdf"]


(* ::Section:: *)
(*Numerical verification plots*)


(* ::Subsection:: *)
(*Finite element mesh*)


(* ::Subsubsection:: *)
(*Nice case*)


Module[
  {
    source,
    a, vBath, mesh, prExt, prInt,
    bCoords, bCoordsExt, bCoordsInt,
    dummyForTrailingCommas
  },
  (* Import mesh *)
  source = "bipolar-nice-verification-mesh.txt";
  {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
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
    ListPlot[{bCoordsExt, bCoordsInt}
      , PlotStyle -> {Blue, Red}
    ],
    {}
    , ImageSize -> 480
  ]
] // Ex["bipolar-nice-verification-mesh.pdf"]


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


Module[
 {source,
  a, vBath, mesh, prExt, prInt,
  bCoords, bCoordsExt, bCoordsInt
 },
  (* Import mesh *)
  source = "bipolar-hot-verification-mesh.txt";
  {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
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
] // Ex["bipolar-hot-verification-mesh.pdf"]


(* ::Subsubsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


Module[
 {source,
  a, vBath, mesh, prExt, prInt,
  bCoords, bCoordsExt, bCoordsInt
 },
  (* Import mesh *)
  source = "bipolar-warm_hot-verification-mesh.txt";
  {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
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
] // Ex["bipolar-warm_hot-verification-mesh.pdf"]


(* ::Subsubsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


Module[
 {source,
  a, vBath, mesh, prExt, prInt,
  bCoords, bCoordsExt, bCoordsInt
 },
  (* Import mesh *)
  source = "bipolar-warm-verification-mesh.txt";
  {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
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
] // Ex["bipolar-warm-verification-mesh.pdf"]


(* ::Subsubsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


Module[
 {source,
  a, vBath, mesh, prExt, prInt,
  bCoords, bCoordsExt, bCoordsInt
 },
  (* Import mesh *)
  source = "bipolar-cold_warm-verification-mesh.txt";
  {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
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
] // Ex["bipolar-cold_warm-verification-mesh.pdf"]


(* ::Subsection:: *)
(*Numerical solution*)


(* ::Subsubsection:: *)
(*Nice case*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "bipolar-nice-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y], Element[{x, y}, mesh]
      , AxesLabel -> Italicise /@ {"x", "y", "T"}
      , PlotLabel -> "Numerical solution"
      , PlotRange -> Full
      , PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["bipolar-nice-verification-solution.png"]


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "bipolar-hot-verification-solution.txt";
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
] // Ex["bipolar-hot-verification-solution.png"]


(* ::Subsubsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "bipolar-warm_hot-verification-solution.txt";
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
] // Ex["bipolar-warm_hot-verification-solution.png"]


(* ::Subsubsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "bipolar-warm-verification-solution.txt";
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
] // Ex["bipolar-warm-verification-solution.png"]


(* ::Subsubsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "bipolar-cold_warm-verification-solution.txt";
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
] // Ex["bipolar-cold_warm-verification-solution.png"]


(* ::Subsection:: *)
(*Relative errors (3D and 2D)*)


(* ::Subsubsection:: *)
(*Nice case*)


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-nice-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["bipolar-nice-verification-rel_error-3d.png"]


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-nice-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
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
] // Ex["bipolar-nice-verification-rel_error-2d.png"]


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-hot-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["bipolar-hot-verification-rel_error-3d.png"]


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-hot-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
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
] // Ex["bipolar-hot-verification-rel_error-2d.png"]


(* ::Subsubsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-warm_hot-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["bipolar-warm_hot-verification-rel_error-3d.png"]


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-warm_hot-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
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
] // Ex["bipolar-warm_hot-verification-rel_error-2d.png"]


(* ::Subsubsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-warm-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["bipolar-warm-verification-rel_error-3d.png"]


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-warm-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
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
] // Ex["bipolar-warm-verification-rel_error-2d.png"]


(* ::Subsubsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-cold_warm-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["bipolar-cold_warm-verification-rel_error-3d.png"]


Module[{source, tSol, mesh, tExact},
  (* Import solution *)
  source = "bipolar-cold_warm-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
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
] // Ex["bipolar-cold_warm-verification-rel_error-2d.png"]


(* ::Section:: *)
(*Figure: bipolar coordinates (bipolar-coordinates)*)


Module[
  {
    xMax, yMax,
    textStyleLabel, textStyleContour,
    uContourText, vContourText, basisVectorText,
    curvilinearStyle,
    basisVectorStyle, basisVectorLength,
    orthogonalityMarkerLength,
    updValues, uValues, uValuesLessThanPi, vValues,
    uForLabel, vForLabel, xShiftForLabel, yShiftForLabel,
    uvBase, xyBase, xyTipU, xyTipV,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 2.5;
  yMax = 1.8;
  (* Text functions *)
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleContour = Style[#, LabelSize["Label"] - 2] & @* LaTeXStyle;
  uContourText[value_] :=
    Italicise["u"] ==
      SeparatedRow["VeryThin"][
        value / Degree,
        AdjustmentBox[
          Magnify["\[Degree]", 1.1]
          , BoxBaselineShift -> -0.1
        ]
      ]
      // DisplayForm
      // textStyleContour;
  vContourText[value_] := Italicise["v"] == value // textStyleContour;
  basisVectorText[coord_] := Subscript[Embolden["a"], Italicise[coord]] // textStyleLabel;
  (* Styles *)
  curvilinearStyle = GeneralStyle["DefaultThick"];
  basisVectorStyle = Directive[Thickness[0.008], Arrowheads[0.03]];
  basisVectorLength = 0.25 xMax;
  orthogonalityMarkerLength = basisVectorLength / 6.5;
  (* Make plot *)
  Show[
    EmptyAxes[{-xMax, xMax}, {-yMax, yMax}
      , AspectRatio -> Automatic
      , AxesStyle -> curvilinearStyle
      , Ticks -> None
    ],
    (* v-contours *)
    vValues = {0.5, 1, 2};
    vValues = Flatten @ {-vValues // Reverse, vValues};
    Graphics @ {
      curvilinearStyle,
      Table[
        {
          Circle[{Coth[v], 0}, Csch[v] // Abs],
          uForLabel =
            Which[
              v > 1, 55.5,
              v < -1, 42,
              True, 63
            ] Degree;
          Text[
            vContourText[v]
            , XYBipolar[uForLabel, v]
            , {0, -1}
            , Sign[-v] AUBipolar[uForLabel, v]
          ]
        }
        , {v, vValues}
      ],
      Text[
        vContourText[0]
        , XYBipolar[67 Degree, 0]
        , {0, 0.7}
        , {0, 1}
      ],
      {}
    },
    (* u-contours *)
    updValues = Range[0, 360, 45] // Select[!Divisible[#, 180] &];
    uValues = updValues * Degree;
    (* (avoid doubled u-contours) *)
    uValuesLessThanPi = uValues // Select[# < Pi &];
    Graphics @ {
      curvilinearStyle,
      Table[
        Circle[{0, Cot[u]}, Csc[u] // Abs]
        , {u, uValuesLessThanPi}
      ],
      (* Labels along y-axis *)
      Table[
        vForLabel = 0;
        xShiftForLabel = If[u == Pi/2, 0.03, -0.105];
        yShiftForLabel = If[u < Pi, -1, -1.15];
        Text[
          uContourText[u]
          , XYBipolar[u, vForLabel]
          , {xShiftForLabel, yShiftForLabel}
          , AVBipolar[u, vForLabel]
        ]
        , {u, uValues // Select[90 <= # / Degree <= 270 &]}
      ],
      {
        Text[
          uContourText[180 Degree]
          , {0, 0}
          , {-0.105, -1}
        ]
      },
      (* Labels along x-axis *)
      {
        Text[
          uContourText[0]
          , XYBipolar[0, 1.32]
          , {0, -1}
        ],
        Text[
          uContourText[0]
          , XYBipolar[0, -1.24]
          , {0, -1}
        ],
        {}
      },
      (* Labels which need doubling *)
      Table[
        xShiftForLabel = If[u > Pi/2, -0.15, -0.1];
        yShiftForLabel = If[u > Pi/2, -1.2, -1];
        {
          vForLabel = 0.68;
          Text[
            uContourText[u]
            , XYBipolar[u, vForLabel]
            , {xShiftForLabel, yShiftForLabel}
            , AVBipolar[u, vForLabel]
          ],
          vForLabel = 0.73;
          Text[
            uContourText[u]
            , XYBipolar[u, -vForLabel]
            , {xShiftForLabel, yShiftForLabel}
            , AVBipolar[u, -vForLabel]
          ]
        }
        , {u, uValues // Select[! 90 <= # / Degree <= 270 &]}
      ],
      {}
    },
    (* Basis vectors with orthogonality marker *)
    uvBase = {270 Degree, 0.5};
    xyBase = XYBipolar @@ uvBase;
    xyTipU = xyBase + basisVectorLength * AUBipolar @@ uvBase;
    xyTipV = xyBase + basisVectorLength * AVBipolar @@ uvBase;
    Graphics @ {
      Line[orthogonalityMarkerLength {{1, 0}, {1, 1}, {0, 1}}]
        // Rotate[#, PhiPolar @@ AUBipolar @@ uvBase, {0, 0}] &
        // Translate[#, xyBase] &
    },
    Graphics @ {basisVectorStyle,
      Arrow @ {xyBase, xyTipU},
      Text[
        basisVectorText["u"]
        , xyTipU
        , {-0.8, 0.5}
      ],
      Arrow @ {xyBase, xyTipV},
      Text[
        basisVectorText["v"]
        , xyTipV
        , {-1.4, -0.5}
      ],
      {}
    },
    {}
    , ImageSize -> 0.87 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
  ]
] // Ex["bipolar-coordinates.pdf"]


(* ::Section:: *)
(*Figure: bipolar arcs of constant u (bipolar-u-*-than-pi)*)


Module[
  {
    uValue,
    singularityNegative, singularityPositive,
    xMax, yMin, yMax,
    textStyle, textStyleBracket, coordinatePairLabel,
    curvilinearStyle,
    u, arcCentre, arcRadius, angStart, angEnd,
    vSign, v, chordVertex,
    angleMarkerRadius, angStartMarker, angEndMarker,
    angleMarkerOffset,
    singularityLabelOffset,
    dummyForTrailingCommas
  },
  (* Constants *)
  uValue["less"] = 65 Degree;
  uValue["more"] = 180 Degree + uValue["less"];
  singularityNegative = {-1, 0};
  singularityPositive = {+1, 0};
  (* Plot range (assumes uValue["less"] < 90 Degree) *)
  xMax = Csc @ uValue["less"] // Abs;
  yMin = YBipolar[uValue["more"], 0];
  yMax = YBipolar[uValue["less"], 0];
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  coordinatePairLabel[x_, y_] :=
    Row @ {
      "(" // textStyleBracket,
      "",
      x,
      ",\[ThinSpace]",
      y,
      ")" // textStyleBracket
    } // textStyle;
  (* Styles *)
  curvilinearStyle = GeneralStyle["DefaultThick"];
  (* Make plots *)
  Table[
    Show[
      u = uValue[case];
      (* u == constant arc *)
      arcCentre = {0, Cot[u]};
      arcRadius = Csc[u] // Abs;
      angStart = PhiPolar @@ (singularityPositive - arcCentre);
      angEnd = PhiPolar @@ (singularityNegative - arcCentre);
        If[u < Pi, angEnd += 2 Pi];
      Graphics @ {curvilinearStyle,
        Circle[arcCentre, arcRadius, {angStart, angEnd}]
      },
      (* Vertex of chords *)
      vSign = If[u < Pi, +1, -1];
      v = 0.4 vSign;
      chordVertex = XYBipolar[u, v];
      Graphics @ {
        Line @ {singularityNegative, chordVertex, singularityPositive}
      },
      (* Angle marker *)
      angStartMarker = PhiPolar @@ (singularityNegative - chordVertex);
      angEndMarker = PhiPolar @@ (singularityPositive - chordVertex);
      If[u < Pi,
        angleMarkerRadius = 0.3;
        angleMarkerOffset = {0.8, 0.55};
        ,
        angleMarkerRadius = 0.2;
        angEndMarker += 2 Pi;
        angleMarkerOffset = {0.5, 0.5};
      ];
      Graphics @ {
        Circle[chordVertex, angleMarkerRadius, {angStartMarker, angEndMarker}],
        Text[
          Italicise["u"] // textStyle
          , chordVertex + XYPolar[angleMarkerRadius, Way[angStartMarker, angEndMarker]]
          , angleMarkerOffset
        ],
        {}
      },
      (* Singularities *)
      singularityLabelOffset = If[u < Pi, {0.18, 1}, {0.22, -1.5}];
      Graphics @ {PointSize[Medium],
        Point @ {singularityNegative, singularityPositive},
        Text[
          coordinatePairLabel[Row @ {"-", Italicise["a"]}, 0]
          , singularityNegative
          , singularityLabelOffset
        ],
        Text[
          coordinatePairLabel[Row @ {"+", Italicise["a"]}, 0]
          , singularityPositive
          , singularityLabelOffset
        ],
        {}
      },
      {}
      , ImageSize -> 0.35 ImageSizeTextWidth
      , PlotRange -> {{-xMax, xMax}, {yMin, yMax}}
      , PlotRangePadding -> {{Scaled[0.12], Scaled[0.07]}, {Scaled[0.15], Scaled[0.03]}}
    ]
      // Ex @ FString["bipolar-u-{case}-than-pi.pdf"]
    , {case, {"less", "more"}}
  ]
]


(* ::Section:: *)
(*Figure: critical terminal points & viable domain (bipolar-viable)*)


(* ::Subsection:: *)
(*Main figure*)


Module[
  {
    regimeList, realA,
    xMinRaw, xMaxRaw, yMaxRaw,
    xMin, xMax, yMax,
    fakeCoshPlus, fakeCoshMinus, fakeQuartic,
    vMinFake, vMaxFake, fMinFake, fMaxFake,
    fakeV, fakeA,
    textStyle, textStyleAxis, textStyleContour,
    terminalStyle, guideStyleFake,
    coshStyleFake, quarticStyleFake,
    curvilinearStyle,
    criticalPlot,
      vFlat0Fake, vSharp0Fake,
      vFlatPiFake, vSharpPiFake,
      fakeVLabelVerticalOffset,
    viablePlot,
      xMinNon, xMaxNon, yMaxNon,
      viablePlotPoints,
    dummyForTrailingCommas
  },
  (* Five real regimes *)
  regimeList = {"hot", "warm_hot", "warm", "cold_warm", "cold"};
  realA = AssociationThread[regimeList ->
    {aHot, aNatPi, aWarm, aNat0, aCold}
  ];
  (* Plot range for non-viable domain plot *)
  xMinRaw = vFlatPi @ realA["hot"] // XBipolar[Pi, #] &;
  xMaxRaw = vFlat0 @ realA["hot"] // XBipolar[0, #] &;
  yMaxRaw = (xMaxRaw - xMinRaw) / 2;
  xMin = Way[xMinRaw, xMaxRaw, -0.2];
  xMax = Way[xMinRaw, xMaxRaw, +1.2];
  yMax = 1.05 yMaxRaw;
  (*
    For the critical terminal points plot,
    we use fake versions of cosh and v^4 because
    the real versions are too bunched up to comprehend.
  *)
  fakeCoshPlus[v_] := v^4 + 1;
  fakeCoshMinus[v_] := v^4 + 1/2;
  fakeQuartic[a_][v_] := v / a;
  {vMinFake, vMaxFake} = {0, 1.3};
  {fMinFake, fMaxFake} = {0, fakeCoshPlus[vMaxFake]};
  (* Five fake regimes *)
  {fakeV["warm_hot"], fakeA["warm_hot"]} =
    FindRoot[
      Function[{v, a},
        {
          fakeCoshPlus[v] - fakeQuartic[a][v],
          fakeCoshPlus'[v] - fakeQuartic[a]'[v]
        }
      ],
      {{0.75}, {0.55}}
    ];
  {fakeV["cold_warm"], fakeA["cold_warm"]} =
    FindRoot[
      Function[{v, a},
        {
          fakeCoshMinus[v] - fakeQuartic[a][v],
          fakeCoshMinus'[v] - fakeQuartic[a]'[v]
        }
      ],
      {{0.6}, {0.95}}
    ];
  fakeA["hot"] = 0.8 fakeA["warm_hot"];
  fakeA["warm"] = Way[fakeA["warm_hot"], fakeA["cold_warm"], 1/3];
  fakeA["cold"] = 2.2 fakeA["cold_warm"];
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleAxis = Style[#, LabelSize["Axis"]] & @* LaTeXStyle;
  textStyleContour = Style[#, LabelSize["Label"] - 2] & @* LaTeXStyle;
  (* Styles *)
  terminalStyle = AbsolutePointSize[6];
  guideStyleFake = Directive[Gray];
  coshStyleFake = Directive[Black, GeneralStyle["DefaultThick"]];
  quarticStyleFake = Directive[Gray, GeneralStyle["DefaultThick"], GeneralStyle["Dashed"]];
  curvilinearStyle = BoundaryTracingStyle["Background"];
  (* Make all plots *)
  Table[
    (* BEGIN Table content *)
    criticalPlot = Show[
      (* Functions of v *)
      Plot[
        {fakeCoshPlus[v], fakeCoshMinus[v], fakeQuartic[fakeA[regime]][v]}
        , {v, vMinFake, vMaxFake}
        , Axes -> None
        , PlotPoints -> 2
        , PlotRange -> {All, {fMinFake, fMaxFake}}
        , PlotRangeClipping -> False
        , PlotRangePadding -> {Automatic, {Scaled[0.15], Automatic}}
        , PlotStyle -> {coshStyleFake, coshStyleFake, quarticStyleFake}
      ],
      (* Manual axes (for better spacing than Plot's axes) *)
      Graphics @ {
        Line @ {
          {vMinFake, fMaxFake},
          {vMinFake, fMinFake},
          {vMaxFake, fMinFake},
          Nothing
        },
        Text[
          vIt // textStyleAxis
          , {vMaxFake, fMinFake}
          , {-3, -0.2}
        ],
        {}
      },
      (* Critical terminal points along u == 0 *)
      fakeVLabelVerticalOffset = 0.55;
      Which[
        (* Two distinct terminal points *)
        realA[regime] < realA["cold_warm"],
        Graphics @ {
          (* v_\[Flat]0 *)
          vFlat0Fake = SeekRoot[
            fakeCoshMinus[#] - fakeQuartic[fakeA[regime]][#] &,
            {vMinFake, fakeV["cold_warm"]}, 4
          ];
          {guideStyleFake,
            Line @ {
              {vFlat0Fake, fakeCoshMinus[vFlat0Fake]},
              {vFlat0Fake, fMinFake}
            }
          },
          {terminalStyle,
            Point @ {vFlat0Fake, fakeCoshMinus[vFlat0Fake]}
          },
          Text[
            Subscript[vIt, Row @ {"\[Flat]\[VeryThinSpace]", 0}] // textStyle
            , {vFlat0Fake, fMinFake}
            , {0, fakeVLabelVerticalOffset}
          ],
          (* v_\[Sharp]0 *)
          vSharp0Fake = SeekRoot[
            fakeCoshMinus[#] - fakeQuartic[fakeA[regime]][#] &,
            {fakeV["cold_warm"], vMaxFake}, 4
          ];
          {guideStyleFake,
            Line @ {
              {vSharp0Fake, fakeCoshMinus[vSharp0Fake]},
              {vSharp0Fake, fMinFake}
            }
          },
          {terminalStyle,
            Point @ {vSharp0Fake, fakeCoshMinus[vSharp0Fake]}
          },
          Text[
            Subscript[vIt, Row @ {"\[Sharp]\[VeryThinSpace]", 0}] // textStyle
            , {vSharp0Fake, fMinFake}
            , {If[regime == "hot", -0.2, 0], fakeVLabelVerticalOffset}
          ],
          {}
        },
        (* One terminal point *)
        realA[regime] == realA["cold_warm"],
        Graphics @ {
          (* v_\[Natural]0 *)
          vFlat0Fake = fakeV["cold_warm"];
          {guideStyleFake,
            Line @ {
              {vFlat0Fake, fakeCoshMinus[vFlat0Fake]},
              {vFlat0Fake, fMinFake}
            }
          },
          {terminalStyle,
            Point @ {vFlat0Fake, fakeCoshMinus[vFlat0Fake]}
          },
          Text[
            Subscript[vIt, Row @ {"\[Natural]\[VeryThinSpace]", 0}] // textStyle
            , {vFlat0Fake, fMinFake}
            , {0, fakeVLabelVerticalOffset}
          ],
          {}
        },
        (* Zero terminal points *)
        True, {}
      ],
      (* Critical terminal points along u == Pi *)
      Which[
        (* Two distinct terminal points *)
        realA[regime] < realA["warm_hot"],
        Graphics @ {
          (* v_\[Flat]\[Pi] *)
          vFlatPiFake = SeekRoot[
            fakeCoshPlus[#] - fakeQuartic[fakeA[regime]][#] &,
            {vMinFake, fakeV["warm_hot"]}, 4
          ];
          {guideStyleFake,
            Line @ {
              {vFlatPiFake, fakeCoshPlus[vFlatPiFake]},
              {vFlatPiFake, fMinFake}
            }
          },
          {terminalStyle,
            Point @ {vFlatPiFake, fakeCoshPlus[vFlatPiFake]}
          },
          Text[
            Subscript[vIt, Row @ {"\[Flat]", "\[Pi]"}] // textStyle
            , {vFlatPiFake, fMinFake}
            , {0, fakeVLabelVerticalOffset}
          ],
          (* v_\[Sharp]\[Pi] *)
          vSharpPiFake = SeekRoot[
            fakeCoshPlus[#] - fakeQuartic[fakeA[regime]][#] &,
            {fakeV["warm_hot"], vMaxFake}, 4
          ];
          {guideStyleFake,
            Line @ {
              {vSharpPiFake, fakeCoshPlus[vSharpPiFake]},
              {vSharpPiFake, fMinFake}
            }
          },
          {terminalStyle,
            Point @ {vSharpPiFake, fakeCoshPlus[vSharpPiFake]}
          },
          Text[
            Subscript[vIt, Row @ {"\[Sharp]", "\[Pi]"}] // textStyle
            , {vSharpPiFake, fMinFake}
            , {0.2, fakeVLabelVerticalOffset}
          ],
          {}
        },
        (* One terminal point *)
        realA[regime] == realA["warm_hot"],
        Graphics @ {
          (* v_\[Natural]\[Pi] *)
          vFlatPiFake = fakeV["warm_hot"];
          {guideStyleFake,
            Line @ {
              {vFlatPiFake, fakeCoshPlus[vFlatPiFake]},
              {vFlatPiFake, fMinFake}
            }
          },
          {terminalStyle,
            Point @ {vFlatPiFake, fakeCoshPlus[vFlatPiFake]}
          },
          Text[
            Subscript[vIt, Row @ {"\[Natural]", "\[Pi]"}] // textStyle
            , {vFlatPiFake, fMinFake}
            , {0, fakeVLabelVerticalOffset}
          ],
          {}
        },
        (* Zero terminal points *)
        True, {}
      ],
      {}
    ];
    viablePlot = Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax}
        , Frame -> False
      ],
      (* Bipolar u-contours through positive singularity *)
      Module[{updValues, uValues, uValuesLessThanPi},
        updValues = Range[0, 360, 45] // Select[!Divisible[#, 180] &];
        uValues = updValues * Degree;
        uValuesLessThanPi = uValues // Select[# < Pi &];
        Graphics @ {
          curvilinearStyle,
          Table[
            Circle[{0, Cot[u]}, Csc[u] // Abs]
            , {u, uValuesLessThanPi}
          ],
          Line @ {{0, 0}, {2 xMax, 0}},
          {}
        }
      ],
      (* Non-viable domain *)
      If[realA[regime] < realA["cold_warm"],
        Which[
          regime == "hot",
            viablePlotPoints = 5;
            xMinNon = vFlatPi @ realA[regime] // XBipolar[Pi, #] &;
            xMaxNon = vFlat0 @ realA[regime] // XBipolar[0, #] &;
            yMaxNon = 1.1 (xMaxNon - xMinNon) / 2;
          ,
          regime == "warm_hot",
            viablePlotPoints = 15;
            xMinNon = vNatPi // 0.97 XBipolar[Pi, #] &;
            xMaxNon = vFlat0 @ realA[regime] // XBipolar[0, #] &;
            yMaxNon = (xMaxNon - xMinNon) / 2;
          ,
          regime == "warm",
            viablePlotPoints = 5;
            xMinNon = vSharp0 @ realA[regime] // 0.97 XBipolar[0, #] &;
            xMaxNon = vFlat0 @ realA[regime] // XBipolar[0, #] &;
            yMaxNon = xMaxNon - xMinNon;
          ,
          True,
            viablePlotPoints = Automatic;
            {xMinNon, xMaxNon, yMaxNon} = {xMin, xMax, yMax};
        ];
        RegionPlot[
          vi @ realA[regime] @@ UVBipolar[x, y] < 0
          , {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon}
          , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
          , PlotPoints -> viablePlotPoints
          , PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        {}
      ],
      (* Critical terminal points along u == 0 *)
      Which[
        (* Two distinct terminal points *)
        realA[regime] < realA["cold_warm"],
        Graphics @ {
          (* v_\[Flat]0 *)
          {terminalStyle,
            Point @ XYBipolar[0, vFlat0 @ realA[regime]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Flat]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vFlat0 @ realA[regime]]
            , {-1.2, -0.85}
          ],
          (* v_\[Sharp]0 *)
          {terminalStyle,
            Point @ XYBipolar[0, vSharp0 @ realA[regime]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Sharp]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vSharp0 @ realA[regime]]
            , {-1.1, -0.9}
          ],
          {}
        },
        (* One terminal point *)
        realA[regime] == realA["cold_warm"],
        Graphics @ {
          (* v_\[Natural]0 *)
          {terminalStyle,
            Point @ XYBipolar[0, vNat0]
          },
          Text[
            Subscript[vIt, Row @ {"\[Natural]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vNat0]
            , {-0.85, -1}
          ],
          {}
        },
        (* Zero terminal points *)
        True, {}
      ],
      (* Critical terminal points along u == pi *)
      Which[
        (* Two distinct terminal points *)
        realA[regime] < realA["warm_hot"],
        Graphics @ {
          (* v_\[Flat]\[Pi] *)
          {terminalStyle,
            Point @ XYBipolar[Pi, vFlatPi @ realA[regime]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Flat]", "\[Pi]"}] // textStyle
            , XYBipolar[Pi, vFlatPi @ realA[regime]]
            , {1.3, -0.9}
          ],
          (* v_\[Sharp]\[Pi] *)
          {terminalStyle,
            Point @ XYBipolar[Pi, vSharpPi @ realA[regime]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Sharp]", "\[Pi]"}] // textStyle
            , XYBipolar[Pi, vSharpPi @ realA[regime]]
            , {-1.5, -0.75}
          ],
          {}
        },
        (* One terminal point *)
        realA[regime] == realA["warm_hot"],
        Graphics @ {
          (* v_\[Natural]\[Pi] *)
          {terminalStyle,
            Point @ XYBipolar[Pi, vNatPi]
          },
          Text[
            Subscript[vIt, Row @ {"\[Natural]", "\[Pi]"}] // textStyle
            , XYBipolar[Pi, vNatPi]
            , {1.4, -1}
          ],
          {}
        },
        (* Zero terminal points *)
        True, {}
      ],
      (* u == const labels *)
      (* (no point trying to automate this *)
      If[realA[regime] == realA["cold"],
        Graphics @ {Gray,
          (* u == 0 *)
          Text[
            uIt == 0 // textStyleContour
            , XYBipolar[0, 0.95 vFlat0 @ realA["warm"]]
            , {0, -1.1}
          ],
          (* u == pi/4 *)
          Text[
            uIt == SeparatedRow["VeryThin"]["\[Pi]", "/", 4] // textStyleContour
            , XYBipolar[Pi/4, vFlat0 @ realA["warm"]]
            , {0, 0.6}
            , -AVBipolar[Pi/4, vFlat0 @ realA["warm"]]
          ],
          (* u == pi/2 *)
          Text[
            uIt == SeparatedRow["VeryThin"]["\[Pi]", "/", 2] // textStyleContour
            , XYBipolar[Pi/2, 0.95 vFlatPi @ realA["hot"]]
            , {0, 0.6}
            , {0, 1}
          ],
          (* u == 3 pi/4 *)
          Text[
            uIt == Row @ {"3\[VeryThinSpace]", SeparatedRow["VeryThin"]["\[Pi]", "/", 4]}
              // textStyleContour
            , XYBipolar[3 Pi/4, 1.03 vFlat0 @ realA["warm"]]
            , {0, -1.1}
            , AVBipolar[3 Pi/4, vFlat0 @ realA["warm"]]
          ],
          (* u == pi *)
          Text[
            uIt == "\[Pi]" // textStyleContour
            , XYBipolar[Pi, 0.93 vFlatPi @ realA["hot"]]
            , {0, -1.1}
            , AVBipolar[Pi, vFlatPi @ realA["hot"]]
          ],
          (* u == 5 pi/4 *)
          Text[
            uIt == Row @ {"5\[VeryThinSpace]", SeparatedRow["VeryThin"]["\[Pi]", "/", 4]}
              // textStyleContour
            , XYBipolar[5 Pi/4, 1.03 vFlat0 @ realA["warm"]]
            , {0, -1.1}
            , AVBipolar[5 Pi/4, vFlat0 @ realA["warm"]]
          ],
          (* u == 3 pi/2 *)
          Text[
            uIt == Row @ {"3\[VeryThinSpace]", SeparatedRow["VeryThin"]["\[Pi]", "/", 2]} // textStyleContour
            , XYBipolar[3 Pi/2, 0.95 vFlatPi @ realA["hot"]]
            , {0, -1.1}
            , {0, -1}
          ],
          (* u == 7 pi/4 *)
          Text[
            uIt == Row @ {"7", SeparatedRow["VeryThin"]["\[Pi]", "/", 4]}
              // textStyleContour
            , XYBipolar[7 Pi/4, 1.03 vFlat0 @ realA["warm"]]
            , {0, -1.1}
            , -AVBipolar[7 Pi/4, vFlat0 @ realA["warm"]]
          ],
          {}
        },
        {}
      ],
      {}
    ];
    {criticalPlot, viablePlot}
    (* END Table content *)
    , {regime, regimeList}
  ]
    // GraphicsGrid[#
      , Alignment -> Center
      , ImageSize -> 0.8 {1, 1.7} ImageSizeTextWidth
    ] &
    // Ex["bipolar-viable.pdf"]
]


(* ::Subsection:: *)
(*Legend*)


Module[
  {
    textStyleLegend,
    coshStyleFake, quarticStyleFake,
    terminalLegend, viableLegend,
    dummyForTrailingCommas
  },
  textStyleLegend = LatinModernLabelStyle @ LabelSize["Legend"];
  (* Terminal points legend (first column) *)
  coshStyleFake = Directive[Black, GeneralStyle["DefaultThick"]];
  quarticStyleFake = Directive[Gray, GeneralStyle["DefaultThick"], GeneralStyle["Dashed"]];
  terminalLegend =
    CurveLegend[
      {coshStyleFake, quarticStyleFake},
      {
        SeparatedRow["VeryThin"][
          SeparatedRow["VeryThin"]["cosh", vIt],
          Style["\[ThinSpace]\[PlusMinus]\[VeryThinSpace]\[VeryThinSpace]", Magnification -> 1.4],
          1
        ],
        SeparatedRow[][vIt^4, Style["/", Magnification -> 1.2], aIt]
      }
      , LabelStyle -> textStyleLegend
    ];
  (* Non-viable domain legend (second column) *)
  viableLegend = Join[
    CurveLegend[
      {BoundaryTracingStyle["Terminal"]},
      {"terminal curve"}
      , LabelStyle -> textStyleLegend
    ],
    RegionLegend[
      {BoundaryTracingStyle["NonViable"]},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> textStyleLegend
    ],
    {}
  ];
  (* Make legend *)
  GraphicsGrid[
    {terminalLegend, viableLegend} // Transpose
    , Alignment -> {Left, Center}
    , ImageSize -> 0.8 ImageSizeTextWidth
    , ItemAspectRatio -> 0.12
    , Spacings -> {0, -0.02 ImageSizeTextWidth}
  ]
    // Ex["bipolar-viable-legend.pdf"]
]


(* ::Subsection:: *)
(*Parameter arrow*)


Module[
  {
    textStyle, arrowStyle,
    width, height,
    aTop, aBottom,
    xWidth, proportionRight, xMin, xMax,
    tickLength,
    aWarmToHot, aColdToWarm,
    dummyForTrailingCommas
  },
  (* Plot styles *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.053], Arrowheads[0.33], CapForm["Butt"]];
  (* Fixed dimensions of exported PDF (to match "bipolar-viable.pdf") *)
  width = 0.1 ImageSizeTextWidth;
  height = 0.8 * 1.7 ImageSizeTextWidth;
  (* Fake A values for vertical coordinates *)
  aTop = 1;
  aBottom = 0;
  aWarmToHot = Way[aTop, aBottom, 1.44 / 5];
  aColdToWarm = Way[aTop, aBottom, 3.46 / 5];
  (* Horizontal coordinates *)
  xWidth = width / height * (aTop - aBottom);
  proportionRight = 0.8;
  xMax = (proportionRight) xWidth;
  xMin = -(1 - proportionRight) xWidth;
  tickLength = 0.15 xWidth;
  (* Make plot *)
  Show[
    (* Parameter (A) increase indicator arrow *)
    Graphics @ {arrowStyle,
      Arrow @ {{0, aTop}, {0, aBottom}}
    },
    Graphics @ {
      Text[
        aIt // textStyle
        , {0, aBottom}
        , {-3, -1}
      ]
    },
    (* A < A_\[Natural]\[Pi] (hot regime) *)
    Graphics @ {
      Text[
        "hot" // textStyle
        , {tickLength, Way[aTop, aBottom, 0.43 / 5]}
        , {0, -0.6}
        , {0, -1.1}
      ]
    },
    (* A == A_\[Natural]\[Pi] (warm-to-hot transition) *)
    Graphics @ {arrowStyle,
      Line @ {{0, aWarmToHot}, {tickLength, aWarmToHot}}
    },
    Graphics @ {
      Text[
        Subscript[aIt, Row @ {"\[VeryThinSpace]\[Natural]", "\[Pi]"}] // textStyle
        , {tickLength, aWarmToHot}
        , {-1.5, 0}
      ]
    },
    (* A_\[Natural]\[Pi] < A < A_\[Natural]0 (warm regime) *)
    Graphics @ {
      Text[
        "warm" // textStyle
        , {tickLength, Way[aTop, aBottom, 2.43 / 5]}
        , {0, -0.6}
        , {0, -1}
      ]
    },
    (* A == A_\[Natural]0 (cold-to-warm transition) *)
    Graphics @ {arrowStyle,
      Line @ {{0, aColdToWarm}, {tickLength, aColdToWarm}}
    },
    Graphics @ {
      Text[
        Subscript[aIt, Row @ {"\[VeryThinSpace]\[Natural]", 0}] // textStyle
        , {tickLength, aColdToWarm}
        , {-1.5, 0}
      ]
    },
    (* A > A_\[Natural]0 (cold regime) *)
    Graphics @ {
      Text[
        "cold" // textStyle
        , {tickLength, Way[aTop, aBottom, 4.45 / 5]}
        , {0, -0.6}
        , {0, -1}
      ]
    },
    {}
    , ImageSize -> {width, height}
    , PlotRange -> {{xMin, xMax}, {aBottom, aTop}}
  ]
    // Ex["bipolar-viable-arrow.pdf"]
]


(* ::Section:: *)
(*Figure: generic traced boundaries (bipolar-traced-boundaries)*)


Module[
  {
    regimeList, aValue,
    xMinRaw, xMaxRaw, yMaxRaw,
    xMin, xMax, yMax,
    tracedStyle,
    a,
      xMinNon, xMaxNon, yMaxNon,
      viablePlotPoints,
      plottedCurvesSpan,
    dummyForTrailingCommas
  },
  (* Five regimes *)
  regimeList = {"hot", "warm", "cold"};
  aValue = AssociationThread[regimeList ->
    {aHot, aWarm, aCold}
  ];
  (* Plot range *)
  xMinRaw = vFlatPi[aHot] // XBipolar[Pi, #] &;
  xMaxRaw = vFlat0[aHot] // XBipolar[0, #] &;
  yMaxRaw = (xMaxRaw - xMinRaw) / 2;
  xMin = Way[xMinRaw, xMaxRaw, -0.75];
  xMax = Way[xMinRaw, xMaxRaw, +1.75];
  yMax = 0.8 (xMax - xMin) / 2;
  (* Slightly less thick traced boundaries *)
  tracedStyle = (
    BoundaryTracingStyle["Traced"]
      /. {AbsoluteThickness[d_] :> AbsoluteThickness[0.9 d]}
  );
  (* Make plots *)
  Table[
    a = aValue[regime];
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax}
        , FrameLabel -> {
            Italicise["x"] // Margined @ {{0, 0}, {0, -20}},
            Italicise["y"]
          }
        , FrameStyle -> LabelSize["Axis"]
        , FrameTicksStyle -> LabelSize["Tick"]
        , LabelStyle -> LatinModernLabelStyle[LabelSize["Label"] - 1]
        , FrameTicks -> {
            {Automatic, None},
            {FindDivisions[{xMin, xMax}, 4] // N, None}
          }
      ],
      (* Non-viable domain *)
      If[a <= aWarm,
        Which[
          regime == "hot",
            viablePlotPoints = 5;
            xMinNon = vFlatPi[a] // XBipolar[Pi, #] &;
            xMaxNon = vFlat0[a] // XBipolar[0, #] &;
            yMaxNon = 1.1 (xMaxNon - xMinNon) / 2;
          ,
          regime == "warm",
            viablePlotPoints = 5;
            xMinNon = vSharp0[a] // 0.97 XBipolar[0, #] &;
            xMaxNon = vFlat0[a] // XBipolar[0, #] &;
            yMaxNon = xMaxNon - xMinNon;
          ,
          True,
            viablePlotPoints = Automatic;
            {xMinNon, xMaxNon, yMaxNon} = {xMin, xMax, yMax};
        ];
        RegionPlot[
          vi[a] @@ UVBipolar[x, y] < 0
          , {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon}
          , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
          , PlotPoints -> viablePlotPoints
          , PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        {}
      ],
      (* Traced boundaries *)
      Which[
        (* Hot regime *)
        regime == "hot",
        plottedCurvesSpan = Association[
          "outer" -> Join[{4}, Range[-6, -3]],
          "inner" -> {1, 5},
          "axis" -> {2, 4, 5, 8},
          "hyperbolic" -> {-1},
          Nothing
        ];
        {
          Table[
            Table[
              ParametricPlot[
                XYBipolar[u, v[u]] // IncludeYReflection
                , {u, DomainStart[v], DomainEnd[v]}
                , PlotPoints -> 2
                , PlotStyle -> tracedStyle
              ]
              , {v, vTraHot[id][[plottedCurvesSpan @ id]]}
            ]
            , {id, {"outer", "inner", "axis", "hyperbolic"}}
          ],
          {}
        },
        (* Warm regime *)
        regime == "warm",
        plottedCurvesSpan = Association[
          "terminal" -> Range[-6, -2],
          "axis" -> {1, 2, 4, 5},
          "hyperbolic" -> All,
          Nothing
        ];
        {
          Table[
            Table[
              ParametricPlot[
                XYBipolar[u, v[u]] // IncludeYReflection
                , {u, DomainStart[v], DomainEnd[v]}
                , PlotPoints -> 2
                , PlotStyle -> tracedStyle
              ]
              , {v, vTraWarm[id][[plottedCurvesSpan @ id]]}
            ]
            , {id, {"terminal", "axis", "hyperbolic"}}
          ],
          {}
        },
        (* Cold regime *)
        regime == "cold",
        plottedCurvesSpan = Association[
          "axis" -> {2, 3},
          "contour" -> {5, 7, 9, 13, 16},
          Nothing
        ];
        {
          Table[
            Table[
              ParametricPlot[
                XYBipolar[u, v[u]] // IncludeYReflection
                , {u, DomainStart[v], DomainEnd[v]}
                , PlotPoints -> 2
                , PlotStyle -> tracedStyle
              ]
              , {v, vTraCold[id][[plottedCurvesSpan @ id]]}
            ]
            , {id, {"axis", "contour"}}
          ],
          {}
        },
        True, {}
      ],
      {}
      , ImageSize -> 0.45 ImageSizeTextWidth
    ]
      // Ex @ FString["bipolar-traced-boundaries-{regime}.pdf"]
    , {regime, regimeList}
  ]
]


(* ::Section:: *)
(*Figure: critical terminal points (bipolar-critical-terminal-points)*)


(* ::Subsection:: *)
(*Main*)


Module[
  {
    regimeList, aValue,
    xMinRaw, xMaxRaw, yMaxRaw,
    xMin, xMax, yMax,
    contour,
    a,
      xMinNon, xMaxNon, yMaxNon,
      viablePlotPoints,
      plottedCurvesSpan,
    textStyle,
    terminalStyle, contourStyle,
    dummyForTrailingCommas
  },
  (* Five regimes *)
  regimeList = {"hot", "warm_hot", "warm", "cold_warm"};
  aValue = AssociationThread[regimeList ->
    {aHot, aWarmHot, aWarm, aColdWarm}
  ];
  (* Plot range *)
  xMinRaw = vFlatPi[aHot] // XBipolar[Pi, #] &;
  xMaxRaw = vFlat0[aHot] // XBipolar[0, #] &;
  yMaxRaw = (xMaxRaw - xMinRaw) / 2;
  xMin = Way[xMinRaw, xMaxRaw, -0.12];
  xMax = Way[xMinRaw, xMaxRaw, +1.12];
  yMax = 0.8 (xMax - xMin) / 2;
  (* Centre and radius of a v-contour *)
  contour[v_, phiRange_: Nothing] :=
    Circle @@ {{Coth[v], 0}, Csch[v] // Abs, phiRange};
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Styles *)
  terminalStyle = GeneralStyle["Point"];
  contourStyle = BoundaryTracingStyle["ContourPlain"];
  (* Make plots *)
  Table[
    a = aValue[regime];
    Show[
      (* Non-viable domain *)
      If[a <= aWarm,
        Which[
          regime == "hot",
            viablePlotPoints = 5;
            xMinNon = vFlatPi[a] // XBipolar[Pi, #] &;
            xMaxNon = vFlat0[a] // XBipolar[0, #] &;
            yMaxNon = 1.1 (xMaxNon - xMinNon) / 2;
          ,
          regime == "warm_hot",
            viablePlotPoints = 15;
            xMinNon = vNatPi // 0.97 XBipolar[Pi, #] &;
            xMaxNon = vFlat0[a] // XBipolar[0, #] &;
            yMaxNon = (xMaxNon - xMinNon) / 2;
          ,
          regime == "warm",
            viablePlotPoints = 5;
            xMinNon = vSharp0[a] // 0.97 XBipolar[0, #] &;
            xMaxNon = vFlat0[a] // XBipolar[0, #] &;
            yMaxNon = xMaxNon - xMinNon;
          ,
          True,
            viablePlotPoints = Automatic;
            {xMinNon, xMaxNon, yMaxNon} = {xMin, xMax, yMax};
        ];
        RegionPlot[
          vi[a] @@ UVBipolar[x, y] < 0
          , {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon}
          , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
          , PlotPoints -> viablePlotPoints
          , PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        {}
      ],
      (* Critical terminal points along u == 0 *)
      Which[
        (* Two distinct terminal points *)
        a < aColdWarm,
        Graphics @ {
          (* v_\[Flat]0 *)
          {contourStyle,
            contour[vFlat0[a], {-1, 1}]
          },
          {terminalStyle,
            Point @ XYBipolar[0, vFlat0[a]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Flat]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vFlat0[a]]
            , {-1.6, -0.1}
          ],
          (* v_\[Sharp]0 *)
          {contourStyle,
            contour[vSharp0[a], If[regime == "hot", 2.4 {-1 ,1.03}, Nothing]]
          },
          {terminalStyle,
            Point @ XYBipolar[0, vSharp0[a]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Sharp]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vSharp0[a]]
            , {-1.6, -0.125}
          ],
          {}
        },
        (* One terminal point *)
        a == aColdWarm,
        Graphics @ {
          (* v_\[Natural]0 *)
          {contourStyle,
            contour[vSharp0[a]]
          },
          {terminalStyle,
            Point @ XYBipolar[0, vNat0]
          },
          Text[
            Subscript[vIt, Row @ {"\[Natural]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vNat0]
            , {-1.65, -0.15}
          ],
          {}
        },
        (* Zero terminal points *)
        True, {}
      ],
      (* Critical terminal points along u == pi *)
      Which[
        (* Two distinct terminal points *)
        a < aWarmHot,
        Graphics @ {
          (* v_\[Flat]\[Pi] *)
          {contourStyle,
            contour[vFlatPi[a], Pi + 2 {-1, 1}]
          },
          {terminalStyle,
            Point @ XYBipolar[Pi, vFlatPi[a]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Flat]", "\[Pi]"}] // textStyle
            , XYBipolar[Pi, vFlatPi[a]]
            , {1.7, -0.2}
          ],
          (* v_\[Sharp]\[Pi] *)
          {contourStyle,
            contour[vSharpPi[a], Pi + 2 {-1, 1}]
          },
          {terminalStyle,
            Point @ XYBipolar[Pi, vSharpPi[a]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Sharp]", "\[Pi]"}] // textStyle
            , XYBipolar[Pi, vSharpPi[a]]
            , {-1.6, -0.2}
          ],
          {}
        },
        (* One terminal point *)
        a == aWarmHot,
        Graphics @ {
          (* v_\[Natural]\[Pi] *)
          {contourStyle,
            contour[vNatPi, Pi + 2 {-1, 1}]
          },
          {terminalStyle,
            Point @ XYBipolar[Pi, vNatPi]
          },
          Text[
            Subscript[vIt, Row @ {"\[Natural]", "\[Pi]"}] // textStyle
            , XYBipolar[Pi, vNatPi]
            , {1.65, -0.2}
          ],
          {}
        },
        (* Zero terminal points *)
        True, {}
      ],
      {}
      , AspectRatio -> Automatic
      , Frame -> None
      , ImageSize -> {1, 0.8} 0.35 ImageSizeTextWidth
      , PlotRange -> {All, {-yMax, yMax}}
      , PlotRangeClipping -> False
    ]
      // Ex @ FString["bipolar-critical-terminal-points-{regime}.pdf"]
    , {regime, regimeList}
  ]
]


(* ::Subsection:: *)
(*Legend*)


Module[{legendLabelStyle},
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  GraphicsGrid[
    List @ Join[
      CurveLegend[
        BoundaryTracingStyle /@ {"ContourPlain", "Terminal"},
        {Row @ {Italicise["T"], "\[Hyphen]contour"}, "terminal curve"}
        , LabelStyle -> legendLabelStyle
      ],
      RegionLegend[
        BoundaryTracingStyle /@ {"NonViable"},
        {"non\[Hyphen]viable domain"}
        , LabelStyle -> legendLabelStyle
      ],
      {}
    ]
    , ImageSize -> 0.85 ImageSizeTextWidth
    , ItemAspectRatio -> 0.15
    , Spacings -> {{0, -0.04} ImageSizeTextWidth, 0}
  ]
    // Ex["bipolar-critical-terminal-points-legend.pdf"]
]


(* ::Section:: *)
(*Figure: candidate boundaries (bipolar-candidates)*)


(* ::Subsection:: *)
(*Main*)


Module[
  {
    aValues, aMin, aMid, aMax,
    xMin, yMax,
    textStyle,
    terminalStyle,
    vInner, vOuter,
    closed, merged,
    dummyForTrailingCommas
  },
  (* A values *)
  aValues = {aMin, aMid, aMax} = {0.993, 0.9992, 1} aNat0;
  (* Plot range *)
  xMin = 0.75;
  yMax = 0.102;
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Styles *)
  terminalStyle = GeneralStyle["Point"];
  (* Make plots *)
  Table[
    (* Compute traced boundary portions *)
    If[a < aMax,
      a = a // N // Rationalize[#, 0] &;
    ];
    With[{v = \[FormalV]},
      vInner =
        NDSolveValue[
          {
            v'[u] == Re @ vTraDer[a][u, v[u]],
            v[0] == vSharp0[a]
          }, v, {u, -Pi, 0}
          , NoExtrapolation
          , PreciseOptions[]
        ];
      vOuter =
        NDSolveValue[
          {
            v'[u] == Re @ vTraDer[a][u, v[u]],
            v[0] == vFlat0[a]
          }, v, {u, -Pi, 0}
          , NoExtrapolation
          , PreciseOptions[]
        ];
    ];
    (* Make plots *)
    Show[
      closed = a > aMin;
      merged = a == aMax;
      (* Candidate boundaries *)
      Table[
        ParametricPlot[
          XYBipolar[u, v[u]]
            // IncludeYReflection
            // Evaluate
          , {u, DomainStart[v], DomainEnd[v]}
          , PlotPoints -> 2
          , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
        ]
        , {v, {vInner, If[merged, Nothing, vOuter]}}
      ],
      (* Critical terminal points *)
      If[merged,
        (* v_\[Natural]0 *)
        Graphics @ {
          {terminalStyle,
            Point @ XYBipolar[0, vNat0]
          },
          Text[
            Subscript[vIt, Row @ {"\[Natural]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vNat0]
            , {1.6, -0.13}
          ],
          {}
        },
        (* v_\[Flat]0 and v_\[Sharp]0 *)
        Graphics @ {
          {terminalStyle,
            Point @ XYBipolar[0, vFlat0[a]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Flat]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vFlat0[a]]
            , {-1.65, -0.13}
          ],
          {terminalStyle,
            Point @ XYBipolar[0, vSharp0[a]]
          },
          Text[
            Subscript[vIt, Row @ {"\[Sharp]\[VeryThinSpace]", 0}] // textStyle
            , XYBipolar[0, vSharp0[a]]
            , {1.6, -0.13}
          ],
          {}
        }
      ],
      {}
      , Axes -> False
      , ImageSize -> {Automatic, 0.25 ImageSizeTextWidth}
      , PlotRange -> {{If[closed, All, xMin], All}, {-yMax, yMax}}
      , PlotRangePadding -> {None, {Scaled[0.01], Scaled[0.25]}}
      , PlotRangeClipping -> False
    ]
    , {a, aValues}
  ]
    // Quiet
    //
      GraphicsGrid[# // List
        , ImageSize -> ImageSizeTextWidth
        , Spacings -> {{0, 0, -0.035} ImageSizeTextWidth, 0}
      ] &
    // Ex["bipolar-candidates.pdf"]
]


(* ::Subsection:: *)
(*Parameter arrow*)


Module[
  {
    textStyle, arrowStyle,
    aLeft, aRight, aTransition,
    tickLength,
    dummyForTrailingCommas
  },
  (* Plot styles *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  (* Fake A values for horizontal coordinate *)
  aLeft = 0;
  aRight = 1;
  aTransition = 0.91;
  (* Make parameter arrow *)
  tickLength = 0.01;
  Show[
    (* A-axis *)
    Graphics @ {arrowStyle,
      Arrow @ {{aLeft, 0}, {aRight, 0}}
    },
    Graphics @ {
      Text[
        aIt // textStyle,
        {1, 0},
        {-2, 0}
      ]
    },
    (* A_\[Natural]0 *)
    Graphics @ {arrowStyle,
      Line @ {
        {aTransition, 0},
        {aTransition, -tickLength}
      }
    },
    Graphics @ {
      Text[
        Subscript[aIt, Row @ {"\[Natural]\[VeryThinSpace]", 0}] // textStyle
        , {aTransition, -tickLength}
        , {0, 1}
      ]
    },
    {}
    , ImageSize -> ImageSizeTextWidth
    , PlotRangePadding -> {{None, Automatic}, Automatic}
  ] // Ex["bipolar-candidates-arrow.pdf"]
]


(* ::Section:: *)
(*Figure: inner candidate boundary (bipolar-inner-candidate)*)


Module[
  {
    a,
    textStyle, textStyleContour,
    textStyleBracket, coordinatePairLabel,
    curvilinearStyle, terminalStyle,
    vInner, vEnd,
    xMinRaw, xMaxRaw,
    xMin, xMax, yMax,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = 9992/10000 aNat0;
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleContour = Style[#, LabelSize["Label"] - 2] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  coordinatePairLabel[u_, v_] :=
    Row @ {
      "(" // textStyleBracket,
      "",
      u,
      ",\[ThinSpace]",
      v,
      "\[NegativeVeryThinSpace])" // textStyleBracket
    } // textStyle;
  (* Styles *)
  curvilinearStyle = BoundaryTracingStyle["Background"];
  terminalStyle = GeneralStyle["Point"];
  (* Compute traced boundary portions *)
  With[{v = \[FormalV]},
    vInner =
      NDSolveValue[
        {
          v'[u] == Re @ vTraDer[a][u, v[u]],
          v[0] == vSharp0[a]
        }, v, {u, -Pi, 0}
        , NoExtrapolation
        , PreciseOptions[]
      ];
    vEnd = vInner[-Pi];
  ];
  (* Plot range *)
  xMinRaw = XBipolar[-Pi, vEnd];
  xMaxRaw = XBipolar[0, vSharp0[a]];
  xMin = Way[xMinRaw, xMaxRaw, -0.43];
  xMax = Way[xMinRaw, xMaxRaw, +1.4];
  yMax = 0.23 (xMax - xMin);
  (* Make plot *)
  Show[
    (* Bipolar u-contours through positive singularity *)
    Module[{updValues, uValues, uValuesLessThanPi},
      updValues = Range[0, 360, 45] // Select[!Divisible[#, 180] &];
      uValues = updValues * Degree;
      uValuesLessThanPi = uValues // Select[# < Pi &];
      Graphics @ {
        curvilinearStyle,
        Table[
          Circle[{0, Cot[u]}, Csc[u] // Abs]
          , {u, uValuesLessThanPi}
        ],
        Line @ {{0, 0}, {2, 0}},
        {}
      }
    ],
    Graphics @ {Gray,
      (* u == 0 *)
      Text[
        uIt == 0 // textStyleContour
        , {Way[XBipolar[0, vSharp0[a]], xMax, 0.55], 0}
        , {0, -1}
      ],
      (* u == pi *)
      Text[
        uIt == "\[Pi]" // textStyleContour
        , {Way[xMinRaw, 1], 0}
        , {0, -1}
      ],
      (* u == -pi *)
      Text[
        uIt == Row @ {-"\[NegativeVeryThinSpace]", "\[Pi]"} // textStyleContour
        , {Way[xMinRaw, 1], 0}
        , {0, +0.4}
      ],
      {}
    },
    (* Candidate boundaries *)
    Table[
      ParametricPlot[
        XYBipolar[u, v[u]]
          // IncludeYReflection
          // Evaluate
        , {u, DomainStart[v], DomainEnd[v]}
        , PlotPoints -> 2
        , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
      ]
      , {v, {vInner}}
    ],
    (* Right-hand end (critical terminal point) *)
    Graphics @ {
      {GeneralStyle["Point"],
        Point @ XYBipolar[0, vSharp0[a] // N]
      },
      Text[
        coordinatePairLabel[
          0,
          Subscript[vIt, Row @ {"\[Sharp]\[VeryThinSpace]", 0}]
        ] // textStyle
        , XYBipolar[0, vSharp0[a]]
        , {-1.2, +0.75}
      ],
      {}
    },
    (* Left-hand end (corner/tip) *)
    Graphics @ {
      {GeneralStyle["Point"],
        Point @ XYBipolar[-Pi, vEnd]
      },
      Text[
        coordinatePairLabel[
          \[MinusPlus] Row @ {"\[NegativeVeryThinSpace]","\[Pi]"},
          Subscript[vIt, "e\[ThinSpace]"]
        ] // textStyle
        , XYBipolar[-Pi, vEnd]
        , {1.15, +0.75}
      ],
      {}
    },
    {}
    , AspectRatio -> Automatic
    , Axes -> None
    , ImageSize -> 0.45 ImageSizeTextWidth
    , PlotRange -> {{xMin, xMax}, {-yMax, yMax}}
  ]
    // Ex["bipolar-inner-candidate.pdf"]
]


(* ::Section:: *)
(*Comparing v values for A = A_\[Natural]0*)


Module[{a, vi, ve},
  a = aNat0;
  vi = vInfl;
  ve = vTraCandCorn[aNat0];
  {
    {"A_\[Natural]0", a},
    {"v_i", vi},
    {"v_e at A = A_\[Natural]0", ve},
    {"absolute difference", ve - vi},
    {"curvature", 2 a^2 / ve^9 (-2 + ve Tanh[ve/2])},
    Nothing
  }
   // N
   // TableForm
] // Ex["bipolar-comparing-v-values-cold_warm.pdf"]


(* ::Section:: *)
(*Figure: verification domain and mesh (bipolar-verification-domain-mesh)*)


Module[
  {
    source,
    a, vBath, mesh, prExt, prInt,
    vSh0,
    sMax, uvTraced,
    radBath, cenBath,
    textStyle,
    plotRangePadding,
    meshGraphic, domainGraphic,
    dummyForTrailingCommas
  },
  (* Import mesh *)
  source = "bipolar-nice-verification-mesh.txt";
  {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
  (* Radiation (traced) boundary *)
  vSh0 = vSharp0[a];
  sMax = Pi (1 - XBipolar[-Pi, vTraCandCorn[a]]);
  uvTraced = With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    NDSolveValue[
      {
        uvTraSystem[a],
        u[0] == 0, v[0] == vSh0,
        WhenEvent[u[s] == -Pi, "StopIntegration"]
      }, {u, v}, {s, -sMax, 0},
      NoExtrapolation
    ]
  ];
  (* Internal (heat bath) boundary *)
  cenBath = {Coth[vBath], 0};
  radBath = Csch[vBath];
  (* Make graphics *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  plotRangePadding = Scaled[0.01];
  meshGraphic = Show[
    mesh["Wireframe"],
    {}
    , PlotRangePadding -> plotRangePadding
  ];
  domainGraphic = Show[
    (* Radiation (traced) boundary *)
    ParametricPlot[
      XYBipolar @@ EvaluatePair[uvTraced, s]
        // IncludeYReflection
      , {s, DomainStart[uvTraced], DomainEnd[uvTraced]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Internal (heat bath) boundary *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      Circle[cenBath, radBath]
    },
    Graphics @ {
      Text[
        Italicise["T"] == Subscript[Italicise["T"], "\[NegativeThinSpace]d"] // textStyle
        , cenBath + {0, radBath}
        , {0, 1.6}
      ]
    },
    Graphics @ {
      Text[
        Italicise["v"] == Subscript[Italicise["v"], "\[VeryThinSpace]d"] // textStyle
        , cenBath + {0, -radBath}
        , {0, -1.7}
      ]
    },
    {}
    , Axes -> None
    , PlotRangePadding -> plotRangePadding
  ];
  (* Export *)
  GraphicsRow[
    {domainGraphic, meshGraphic}
    , ImageSize -> 0.9 ImageSizeTextWidth
    , Spacings -> Scaled[0.05]
  ]
    // Ex["bipolar-verification-domain-mesh.pdf"]
]


(* ::Subsection:: *)
(*Domain details*)


Module[
  {
    source,
    a, vBath, mesh, prExt, prInt,
    radBath, cenBath,
    dummyForTrailingCommas
  },
  (* Import mesh *)
  source = "bipolar-nice-verification-mesh.txt";
  {a, vBath, mesh, prExt, prInt} = Import[source] // Uncompress;
  (* Internal (heat bath) boundary *)
  cenBath = {Coth[vBath], 0};
  radBath = Csch[vBath];
  (* Return *)
  {
    {"A", a // N},
    {"v_d", vBath},
    {"v_d / v_\[Sharp]0", vBath / vSharp0[a]},
    {"radius_d", radBath},
    Nothing
  } // TableForm
]


(* ::Subsection:: *)
(*Number of mesh elements*)


Module[{mesh},
  mesh = Import["bipolar-nice-verification-mesh.txt"] // Uncompress // #[[3]] &;
  mesh
]


(* ::Subsection:: *)
(*Maximum relative error throughout mesh*)


Module[{source, tSol, mesh, tExact, relError},
  (* Import solution *)
  source = "bipolar-nice-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
  relError[x_, y_] := tSol[x, y] / tExact[x, y] - 1;
  relError @@@ mesh["Coordinates"] // Abs // Max // PercentForm
]


(* ::Section:: *)
(*Figure: verification solution and relative error (bipolar-verification-*)*)


(* ::Subsection:: *)
(*Solution*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "bipolar-nice-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  Module[{x, y},
    Plot3D[tSol[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{0, 25}, {0, -20}},
          Italicise["y"],
          None
        }
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> {Automatic, Automatic, 0.07}
      , ImageSize -> 0.42 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotLabel -> (Style["Numerical solution", LabelSize["Axis"] - 1] // Margined @ {{0, 0}, {0, 0}})
      , PlotRange -> Full
      , PlotStyle -> Directive[GeneralStyle["Translucent"], BoundaryTracingStyle["Solution3D"]]
      , TicksStyle -> LabelSize["Tick"]
      , ViewPoint -> {1.6, -2.4, 1.5}
    ]
  ]
] // Ex["bipolar-verification-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Subsection:: *)
(*Relative error*)


Module[{source, tSol, mesh, tExact, relError},
  (* Import solution *)
  source = "bipolar-nice-verification-solution.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := VBipolar[x, y];
  (* Relative error *)
  relError[x_, y_] := tSol[x, y] / tExact[x, y] - 1 // Abs;
  Module[{x, y},
    Plot3D[relError[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{0, 25}, {0, -20}},
          Italicise["y"],
          None
        }
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> {Automatic, Automatic, 0.09}
      , ImageSize -> 0.42 {1, 532/691} ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotLabel -> (Style["Relative error", LabelSize["Axis"] - 1] // Margined @ {{0, 20}, {0, 0}})
      , PlotRange -> Full
      , PlotRangePadding -> {Scaled /@ {0.04, 0.01}, Scaled[0.05], Scaled /@ {0, 0.25}}
      , PlotStyle -> Directive[GeneralStyle["Translucent"], BoundaryTracingStyle["Solution3D"]]
      , TicksStyle -> LabelSize["Tick"]
      , ViewPoint -> {1.6, -2.4, 1.5}
    ]
  ]
] // Ex["bipolar-verification-relative-error.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Figure: inner candidate boundary asymmetry (bipolar-inner-candidate-asymmetry)*)


Module[
  {
    aAsymmTable, asymmNatPi, asymmNat0,
    aMinTable, asymmMinTable,
    aMin, aMax, asymmMin,
    asymmInterpolation,
    textStyle, textStyleTick,
      aTicks, tickProcess,
      asymmLabelHorizontalOffset,
      aLabelVerticalOffset,
      aHotLabel, asymmHotLabel,
      aWarmLabel, asymmWarmLabel,
      aColdLabel, asymmColdLabel,
    dummyForTrailingCommas
  },
  (* Import pre-computed asymmetries *)
  {aAsymmTable, asymmNatPi, asymmNat0} =
    Import["bipolar-candidate-asymmetry.txt"] // Uncompress;
  aAsymmTable = aAsymmTable // DeleteCases[{_, Indeterminate}];
  asymmInterpolation = Interpolation[aAsymmTable] @* N;
  (* Axis ranges *)
  {aMinTable, asymmMinTable} = aAsymmTable // First;
  aMin = aMinTable;
  aMax = aNat0 + (aNatPi - aMin);
  asymmMin = 0;
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleTick = Style[#, LabelSize["Tick"]] & @* LaTeXStyle;
  (* Make plot *)
  tickProcess[value_] := If[IntegerQ[value], value, N[value]];
  SetAttributes[tickProcess, Listable];
  aTicks = (
    FindDivisions[{aMin, aMax}, 6]
      // DeleteCases[9] (* prevent clash with A_\[Natural]\[Pi] *)
      // tickProcess (* make non-integers numeric *)
  );
  asymmLabelHorizontalOffset = 1.25;
  aLabelVerticalOffset = 1.2;
  ListPlot[
    aAsymmTable
    , AspectRatio -> 0.55
    , AxesLabel -> {
        aIt // Margined @ {{-3, 1}, {5, 0}},
        "Asymmetry" // Margined @ {{0, 0}, {-5, 0}}
      }
    , Epilog -> {
        (* Hot regime *)
        aHotLabel = Way[aMin, aNatPi, 0.49];
        asymmHotLabel = 0.53 asymmInterpolation[aHotLabel];
        Text[
          "hot" // textStyle
          , {aHotLabel, asymmHotLabel}
        ],
        (* Hot-to-warm transition *)
        {BoundaryTracingStyle["Contour"],
          Line @ {{aMin, asymmNatPi}, {aNatPi, asymmNatPi}},
          Line @ {{aNatPi, asymmMin}, {aNatPi, asymmNatPi}}
        },
        Text[
          asymmNatPi // SignificantFiguresForm[5] // textStyleTick
          , {aMin, asymmNatPi}
          , {asymmLabelHorizontalOffset, -0.2}
        ],
        Text[
           Subscript[aIt, Row @ {"\[VeryThinSpace]\[Natural]", "\[Pi]"}] // textStyle
          , {aNatPi, asymmMin}
          , {0, aLabelVerticalOffset}
        ],
        (* Warm regime *)
        aWarmLabel = Way[aNatPi, aNat0];
        asymmWarmLabel = 1/2 asymmInterpolation[aWarmLabel];
        Text[
          "warm" // textStyle
          , {aWarmLabel, asymmWarmLabel}
        ],
        (* Warm-to-cold transition *)
        {BoundaryTracingStyle["Contour"],
          Line @ {{aMin, asymmNat0}, {aNat0, asymmNat0}},
          Line @ {{aNat0, asymmMin}, {aNat0, asymmNat0}}
        },
        Text[
          asymmNat0 // SignificantFiguresForm[5] // textStyleTick
          , {aMin, asymmNat0}
          , {asymmLabelHorizontalOffset, 0}
        ],
        Text[
          Subscript[aIt, Row @ {"\[Natural]\[VeryThinSpace]", 0}] // textStyle
          , {aNat0, asymmMin}
          , {0, aLabelVerticalOffset}
        ],
        (* Cold regime *)
        {aColdLabel, asymmColdLabel} =
          Way[
            {aHotLabel, asymmHotLabel},
            {aWarmLabel, asymmWarmLabel}
            , 1.96
          ];
        Text[
          "cold" // textStyle
          , {aColdLabel, asymmColdLabel}
        ],
        {}
      }
    , ImagePadding -> {Scaled /@ {0.075, 0.04}, Scaled /@ {0.05, 0.05}}
    , ImageSize -> 0.6 ImageSizeTextWidth
    , Joined -> True
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {{aMin, aMax}, {asymmMin, All}}
    , PlotRangeClipping -> False
    , PlotRangePadding -> {{Automatic, Scaled[0.02]}, {0, Automatic}}
    , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
    , Ticks -> {aTicks, Automatic}
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes] // Evaluate
  ]
] // Ex["bipolar-inner-candidate-asymmetry.pdf"]


(* ::Section:: *)
(*Figure: inner candidate boundary with v_\[Sharp]0 circle (bipolar-inner-candidate-with-circle)*)


Module[
  {
    a,
    textStyle, textStyleContour,
    textStyleBracket, coordinatePairLabel,
    curvilinearStyle, terminalStyle,
    vInner,
      cenCircle, radCircle,
      phiRadiusMarker,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = Floor[aNat0, 1/100];
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleContour = Style[#, LabelSize["Label"] - 2] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  coordinatePairLabel[u_, v_] :=
    Row @ {
      "(" // textStyleBracket,
      "",
      u,
      ",\[ThinSpace]",
      v,
      "\[NegativeVeryThinSpace])" // textStyleBracket
    } // textStyle;
  (* Styles *)
  curvilinearStyle = BoundaryTracingStyle["Background"];
  terminalStyle = GeneralStyle["Point"];
  (* Compute traced boundary portions *)
  With[{v = \[FormalV]},
    vInner =
      NDSolveValue[
        {
          v'[u] == Re @ vTraDer[a][u, v[u]],
          v[0] == vSharp0[a]
        }, v, {u, -Pi, 0}
        , NoExtrapolation
        , PreciseOptions[]
      ];
  ];
  (* Make plot *)
  Show[
    (* Circle (v-contour for v_\[Sharp]0) *)
    cenCircle = {Coth[vSharp0[a]], 0};
    radCircle = Csch[vSharp0[a]];
    phiRadiusMarker = 135 Degree;
    Graphics @ {
      BoundaryTracingStyle["ContourPlain"],
        Circle[cenCircle, radCircle],
        Line @ {cenCircle, cenCircle + XYPolar[radCircle, phiRadiusMarker]},
      Darker[Gray],
        (* v == v_\[Sharp]0 *)
        Text[
          vIt == Subscript[vIt, Row @ {"\[Sharp]\[VeryThinSpace]", 0}] // textStyle
          , cenCircle + {0, -radCircle}
          , {0.1, -1.55}
        ],
        (* r_\[Sharp]0 *)
        Text[
          Subscript[Italicise["r"], Row @ {"\[NegativeVeryThinSpace]\[Sharp]\[VeryThinSpace]", 0}] // textStyle
          , cenCircle + 1/2 XYPolar[radCircle, phiRadiusMarker]
          , {-1.3, -0.57}
        ],
      {}
    },
    (* Candidate boundaries *)
    Table[
      ParametricPlot[
        XYBipolar[u, v[u]]
          // IncludeYReflection
          // Evaluate
        , {u, DomainStart[v], DomainEnd[v]}
        , PlotPoints -> 2
        , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
      ]
      , {v, {vInner}}
    ],
    (* Right-hand end (critical terminal point) *)
    Graphics @ {
      {GeneralStyle["Point"],
        Point @ XYBipolar[0, vSharp0[a] // N]
      },
      Text[
        coordinatePairLabel[
          0,
          Subscript[vIt, Row @ {"\[Sharp]\[VeryThinSpace]", 0}]
        ] // textStyle
        , XYBipolar[0, vSharp0[a]]
        , {-1.35, -0.2}
      ],
      {}
    },
    {}
    , AspectRatio -> Automatic
    , Axes -> None
    , ImageSize -> 0.52 ImageSizeTextWidth
    , PlotRange -> All
    , PlotRangeClipping -> False
  ]
    // Ex["bipolar-inner-candidate-with-circle.pdf"]
]


(* ::Section:: *)
(*Figure: auxiliary function (bipolar-auxiliary-function)*)


Module[
  {
    vSpecial, omegaSpecial,
    textStyle, textStyleSmaller,
    dummyForTrailingCommas
  },
  (* Special values *)
  vSpecial = vNat0;
  omegaSpecial = omega[vSpecial];
  (* Make plot *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleSmaller = Style[#, LabelSize["Label"] - 2] & @* LaTeXStyle;
  Plot[omega[v], {v, 0, 20}
    , AxesLabel -> {
        vIt // Margined @ {{-2, 0}, {3, 0}},
        "\[Omega]" // LaTeXStyle // Margined @ {{0, 0}, {-2, 0}}
      }
    , Epilog -> {
        {BoundaryTracingStyle["Contour"],
          Line @ {{0, omegaSpecial}, {vSpecial, omegaSpecial}},
          Line @ {{vSpecial, 0}, {vSpecial, omegaSpecial}}
        },
        Text[
          SeparatedRow["VeryThin"][1, Style["/", Larger], 4] // textStyleSmaller
          , {0, omegaSpecial}
          , {1.5, -0.07}
        ],
        Text[
          Subscript[vIt, Row @ {"\[Natural]\[VeryThinSpace]", 0}] // textStyle
          , {vSpecial, 0}
          , {0.25, 1}
        ],
        {}
      }
    , ImagePadding -> {Scaled /@ {0.04, 0.04}, Scaled /@ {0.04, 0.04}}
    , ImageSize -> 0.55 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {0, All}
    , PlotRangeClipping -> False
    , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes] // Evaluate
  ]
] // Ex["bipolar-auxiliary-function.pdf"]


(* ::Section:: *)
(*Physical range*)


Module[
  {
    r, eps, k, sigma,
    tMax, pMax,
    bipolarPrefactor, linePrefactor,
    dummyForTrailingCommas
  },
  (* Values taken from "Physical range" in line.wl. *)
  r = Quantity[4, "Centimeters"];
  eps = 0.93;
  k = Quantity[0.18, "Watts" / "Meters" / "Kelvins"];
  sigma = Quantity[5.67 10^-8, "Watts" / "Meters"^2 / "Kelvins"^4];
  (* Temperature range *)
  tMax = 1/2^(2/3) * (k / (eps sigma r))^(1/3);
  (* Power per length *)
  pMax = 2^(1/3) / vNat0 * Pi k^(4/3) / (eps sigma r)^(1/3);
  (* Pretty table *)
  {
    {"Terminal radius", r},
    {"Emissivity", eps},
    {"Conductivity", k},
    {"Temperature", tMax},
    {"Power per length", pMax},
    {"Bipolar prefactor (2^(1/3) / v_\[Natural]0)",
      bipolarPrefactor = 2^(1/3) / vNat0;
      bipolarPrefactor // N
    },
    {"Line prefactor (1 / 2^(5/3))",
      linePrefactor = 1 / 2^(5/3);
      linePrefactor // N
    },
    {"Prefactor discrepancy",
      bipolarPrefactor / linePrefactor - 1
        // N // PercentForm
    },
    Nothing
  } // TableForm
] // Ex["bipolar-physical-range.pdf"]
