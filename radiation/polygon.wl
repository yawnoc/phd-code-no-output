(* ::Package:: *)

(* ::Text:: *)
(*See r4 (manuscripts/radiation-4-polygon.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "polygon"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsubsection:: *)
(*H_n(\[Zeta])*)


(* ::Text:: *)
(*See (r4.8) (Page r4-2).*)


h[n_][zeta_] := Hypergeometric2F1[1/n, 2/n, 1 + 1/n, zeta ^ n];


(* ::Subsubsection:: *)
(*dz/d\[Zeta]*)


(* ::Text:: *)
(*See (r4.11) (Page r4-3).*)


zMapDer[n_][zeta_] := 1 / (h[n][1] (1 - zeta ^ n) ^ (2/n)) // Evaluate;


(* ::Subsubsection:: *)
(*z = z(\[Zeta])*)


(* ::Text:: *)
(*See (r4.12) (Page r4-3).*)


zMap[n_][zeta_] := zeta h[n][zeta] / h[n][1] // Evaluate;


(* ::Subsubsection:: *)
(*d\[Zeta]/dz*)


(* ::Text:: *)
(*See (r4.13) (Page r4-3).*)


zetaMapDer[n_][zeta_] := h[n][1] (1 - zeta ^ n) ^ (2/n) // Evaluate;


(* ::Subsubsection:: *)
(*(Check these)*)


With[{n = \[FormalN], zeta = \[FormalZeta]},
  {
    zMap[n]'[zeta] == zMapDer[n][zeta] // FullSimplify,
    zMapDer[n][zeta] zetaMapDer[n][zeta] == 1
  }
]


(* ::Subsection:: *)
(*Backward Schwarz--Christoffel transformation*)


(* ::Text:: *)
(*Solves z(\[Zeta]) = z by calling FindRoot with initial guess \[Zeta] = z.*)


zetaMap[n_][z_] := SeekRoot[zMap[n][#] - z &, z];


(* ::Subsection:: *)
(*Known solution (in \[Zeta]-space)*)


(* ::Subsubsection:: *)
(*G = -log(\[Zeta])*)


(* ::Text:: *)
(*See (r4.14) (Page r4-3).*)


g[zeta_] := -Log[zeta];


(* ::Subsubsection:: *)
(*G = \[Gamma] - log(\[Zeta]) (offset version)*)


(* ::Text:: *)
(*See (r4.50) (Page r4-14).*)
(*Derivatives of G are unaffected by the offset.*)


gOffset[gamma_][zeta_] := gamma - Log[zeta];


(* ::Subsubsection:: *)
(*dG/dz*)


(* ::Text:: *)
(*See (r4.23) (Page r4-4).*)


gDer[n_][zeta_] := -h[n][1] (1 - zeta ^ n) ^ (2/n) / zeta // Evaluate;


(* ::Subsubsection:: *)
(*dG/d\[Zeta]*)


gDerZeta[zeta_] := D[g[zeta], zeta] // Evaluate;


(* ::Subsection:: *)
(*Flux*)


(* ::Subsubsection:: *)
(*F*)


(* ::Text:: *)
(*See (r4.21) (Page r4-4).*)


f[a_][zeta_] := - (Re @ g[zeta]) ^ 4 / a // Evaluate;


(* ::Subsubsection:: *)
(*F (offset version)*)


fOffset[gamma_][a_][zeta_] := - (Re @ gOffset[gamma][zeta]) ^ 4 / a // Evaluate;


(* ::Subsection:: *)
(*Viable domain*)


(* ::Subsubsection:: *)
(*\[CapitalPhi] (viability)*)


(* ::Text:: *)
(*See (r4.22) (Page r4-4).*)


vi[n_][a_][zeta_] :=
  Subtract[
    (Abs @ gDer[n][zeta]) ^ 2,
    (Re @ g[zeta]) ^ 8 / a^2
  ] // Evaluate;


(* ::Subsubsection:: *)
(*\[CapitalPhi] (viability) (offset version)*)


viOffset[gamma_][n_][a_][zeta_] :=
  Subtract[
    (Abs @ gDer[n][zeta]) ^ 2,
    (Re @ gOffset[gamma][zeta]) ^ 8 / a^2
  ] // Evaluate;


(* ::Subsubsection:: *)
(*\[Psi]*)


(* ::Text:: *)
(*See (r4.24) (Page r4-4).*)


psi[n_][zeta_] := (Re @ g[zeta]) ^ 4 / (Abs @ gDer[n][zeta]);


(* ::Subsubsection:: *)
(*\[Psi] (offset version)*)


psiOffset[gamma_][n_][zeta_] :=
  (Re @ gOffset[gamma][zeta]) ^ 4 / (Abs @ gDer[n][zeta]) // Evaluate;


(* ::Subsubsection:: *)
(*\[Rho]_\[Natural]*)


(* ::Text:: *)
(*The radius \[Rho] = |\[Zeta]| achieving the maximum \[Psi] along arg(\[Zeta]) = \[CurlyPhi].*)
(*See (r4.30) (Page r4-6) and "polygon-psi-algebra.pdf".*)


rhoNat[n_?NumericQ][ph_?NumericQ] :=
  SeekRoot[
    -4 + 8 # ^ n Cos[n ph] + # ^ (2 n) (-4 + Log[#]) - Log[#] &,
    {0, 0.1},
    100
  ]


(* ::Subsubsection:: *)
(*A_\[Natural]*)


(* ::Text:: *)
(*Maximum \[Psi] along arg(\[Zeta]) = \[CurlyPhi].*)


aNat[n_?NumericQ][ph_?NumericQ] := psi[n][rhoNat[n][ph] Exp[I ph]];


(* ::Subsubsection:: *)
(*\[Rho]_\[Sharp]*)


(* ::Text:: *)
(*The larger solution to \[Psi] = A along arg(\[Zeta]) = \[CurlyPhi].*)


rhoSharp[n_?NumericQ][a_?NumericQ][ph_?NumericQ] :=
  Quiet[
    SeekRoot[
      psi[n][# Exp[I ph]] - a &,
      {1, rhoNat[n][ph]},
      100
    ]
  , {Power::infy, Infinity::indet}];


(* ::Subsection:: *)
(*Representative values of A*)


(* ::Text:: *)
(*"Hot" means A < A_\[Natural] along \[CurlyPhi] = 0.*)
(*We take A = A_\[Natural] / 2 as the representative value.*)


aHot[n_] := aNat[n][0] / 2;


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* ::Text:: *)
(*We choose points with \[Rho] between \[Rho]_\[Sharp] and 1,*)
(*and \[CurlyPhi] between 0 and 2 \[Pi] / n (since there is n-fold symmetry).*)


(* ::Subsubsection:: *)
(*Hot regime*)


Table[
  startZetaHot[n]["general"] =
    Module[
     {a,
      rhoMin, rhoMax, rho,
      phMax, phSpacing, phValues
     },
      a = aHot[n];
      (* \[Rho] *)
      rhoMin = rhoSharp[n][a][0];
      rhoMax = 1;
      rho = Way[rhoMin, rhoMax, 1/2];
      (* \[CurlyPhi] values *)
      phMax = 2 Pi / n;
      phSpacing = 30 Degree;
      phValues = UniformRange[0, phMax, phSpacing] // Most;
      (* Build list of values *)
      Table[rho Exp[I ph], {ph, phValues}]
    ]
, {n, 3, 5}];


Table[
  startZetaHot[n]["hyperbolic"] =
    Module[{a, rhoSh},
      a = aHot[n];
      (* Build list of values *)
      (* (only need singleton due to rotational symmetry) *)
      {rhoSharp[n][a][0]}
    ]
, {n, 3, 5}];


(* ::Subsection:: *)
(*Traced boundaries \[Zeta] = \[Zeta](s)*)


(* ::Subsubsection:: *)
(*\[Zeta]' = d\[Zeta]/ds*)


(* ::Text:: *)
(*See (r4.38) (Page r4-8).*)


zetaVel[n_][a_][zeta_] :=
  Divide[
    I f[a][zeta] + Sqrt @ vi[n][a][zeta],
    gDerZeta[zeta]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Solver for traced boundaries*)


(* ::Text:: *)
(*See (r4.38) (Page r4-8).*)


With[{zeta = \[FormalZeta]},
  zetaTrace[n_][a_][zetaInit_] :=
    NDSolveValue[
      {
        zeta'[s] == zetaVel[n][a] @ zeta[s],
        zeta[0] == zetaInit,
        WhenEvent[
          vi[n][a] @ zeta[s] < 0 || Abs @ zeta[s] > 1,
          "StopIntegration"
        ]
      }, zeta, {s, -Pi, Pi},
      NoExtrapolation
    ]
];


(* ::Subsubsection:: *)
(*Solver for hyperbolic traced boundaries*)


zetaTraceHyperbolic[n_][a_] :=
  Module[{rhoClearance},
    rhoClearance = 10^-8;
    zetaTrace[n][a] @ (rhoSharp[n][a][0] + rhoClearance)
  ];


(* ::Subsubsection:: *)
(*Hot regime*)


With[{zeta = \[FormalZeta]},
  Table[
    Module[{a, idList, zetaInitList},
      (* A *)
      a = aHot[n];
      (* Group names *)
      idList = {"general", "hyperbolic"};
      (* Solve for traced boundaries *)
      Table[
        zetaInitList = startZetaHot[n][id];
        zetaTraHot[n][id] =
          Table[
            zetaTrace[n][a][zetaInit]
          , {zetaInit, zetaInitList}];
      , {id, idList}];
    ]
  , {n, 3, 5}];
];


(* ::Subsection:: *)
(*Traced boundary curvature*)


(* ::Subsubsection:: *)
(*\[CurlyPhi](\[Rho] = 1)*)


phEnd[n_][a_] :=
  Module[{zeta},
    zeta = zetaTraceHyperbolic[n][a];
    zeta @ DomainStart[zeta] // Conjugate // Arg
  ];


(* ::Subsubsection:: *)
(*A_m (meeting dimensionless group)*)


(* ::Text:: *)
(*This is the smallest A > 0 for which \[CurlyPhi] = \[Pi]/n at \[Rho] = 1.*)
(*See Page r4-10.*)


(* Compute A_m using the bisection algorithm *)
aMeet = Module[
 {dest, nValues,
  aMin, aMax,
  a, num,
  aAss, numAss
 },
  (* (This is not slow, nevertheless compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "polygon-a-meeting.txt";
  nValues = Range[3, 5];
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {aAss, numAss} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aAss = Association[];
    numAss = Association[];
    aMin = 5/10;
    aMax = 12/10;
    (* NOTE: these may not work for n > 5 *)
    Table[
      {a, num} = SeekRootBisection[
        phEnd[n][#] - Pi / n &,
        {aMin, aMax},
        "ReturnIterations" -> True
      ];
      aAss = aAss // Append[n -> a];
      numAss = numAss // Append[n -> num];
    , {n, nValues}];
    {aAss, numAss} // Compress // Ex[dest]
  ];
  (* Print iterations used *)
  Table[
    Print[
      "n == {n} bisection algorithm: {num} iterations"
        // PrettyString[
          "{n}" -> ToString[n],
          "{num}" -> ToString @ numAss[n]
        ]
    ];
  , {n, nValues}];
  (* Return A_m values *)
  aAss
];
aMeet // N


(* ::Subsubsection:: *)
(*Solver for candidate traced boundaries*)


(* ::Text:: *)
(*This is essentially the same as zetaTraceHyperbolic,*)
(*but with an added termination condition \[CurlyPhi] at = \[Pi]/n*)
(*so that this need not be solved for later.*)


With[{zeta = \[FormalZeta]},
  zetaTraceCand[n_][a_] :=
    Module[{rhoClearance},
      rhoClearance = 10^-8;
      NDSolveValue[
        {
          zeta'[s] == zetaVel[n][a] @ zeta[s],
          zeta[0] == rhoSharp[n][a][0] + rhoClearance,
          WhenEvent[
            Or[
              Arg @ Conjugate @ zeta[s] > Pi / n,
              vi[n][a] @ zeta[s] < 0,
              Abs @ zeta[s] > 1
            ],
            "StopIntegration"
          ]
        }, zeta, {s, -Pi, 0},
        NoExtrapolation
      ]
    ];
];


(* ::Subsubsection:: *)
(*z' = dz/ds*)


(* ::Text:: *)
(*See (r4.39) (Page r4-11).*)


zVel[n_][a_][zeta_] :=
  Divide[
    I f[a][zeta] + Sqrt @ vi[n][a][zeta],
    gDer[n][zeta]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*z'' = d^2(z)/ds^2*)


(* ::Text:: *)
(*See (r4.40) to (r4.43) (Pages r4-11 & r4-12).*)


zAcc[n_][a_][zeta_] := Plus[
  (* 1st term of (r4.40) *)
  Times[
    (* i F + sqrt(\[CapitalPhi]) *)
    I f[a][zeta] + Sqrt @ vi[n][a][zeta],
    (* d\[Zeta]/ds *)
    zetaVel[n][a][zeta],
    (* d/d\[Zeta] (1 / (dG/dz)) *)
    D[1 / gDer[n][zeta], zeta]
  ],
  (* 2nd term of (r4.40) *)
  Divide[
    Plus[
      (* i F' *)
      Times[
        I,
        -4 (Re @ g[zeta])^3 / a,
        Re @ (zetaVel[n][a][zeta] gDerZeta[zeta])
      ],
      (* (sqrt(\[CapitalPhi]))' == \[CapitalPhi]' / (2 sqrt(\[CapitalPhi]))  *)
      Divide[
        Plus[
          (* 1st term of (r4.43) *)
          Times[
            2,
            Re @ gDer[n][zeta],
            Re @ (zetaVel[n][a][zeta] D[gDer[n][zeta], zeta])
          ],
          (* 2nd term of (r4.43) *)
          Times[
            2,
            Im @ gDer[n][zeta],
            Im @ (zetaVel[n][a][zeta] D[gDer[n][zeta], zeta])
          ],
          (* 3rd term of (r4.43) *)
          Times[
            -8 (Re @ g[zeta])^7 / a^2,
            Re @ (zetaVel[n][a][zeta] gDerZeta[zeta])
          ]
        ],
        2 Sqrt @ vi[n][a][zeta]
      ]
    ],
    (* dG/dz *)
    gDer[n][zeta]
  ]
] // Evaluate;


(* ::Subsubsection:: *)
(*Curvature x' y'' - y' x''*)


(* ::Text:: *)
(*See (r4.39a) (Page r4-11).*)


curTra[n_][a_][zeta_] :=
  Module[{xVel, yVel, xAcc, yAcc},
    {xVel, yVel} = zVel[n][a][zeta] // ReIm;
    {xAcc, yAcc} = zAcc[n][a][zeta] // ReIm;
    xVel yAcc - yVel xAcc
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Curvature at \[CurlyPhi] = \[Pi]/n*)


(* ::Text:: *)
(*The curvature evaluated at the corners formed by the candidate boundaries.*)


curCandCorner[n_][a_] /; aMeet[n] <= a < aNat[n][0] :=
  Module[{zeta},
    zeta = zetaTraceCand[n][a];
    zeta @ DomainStart[zeta] // curTra[n][a]
  ];


(* ::Subsubsection:: *)
(*A_i (inflection dimensionless group)*)


(* ::Text:: *)
(*This is the smallest A > 0 for which \[CurlyPhi] = \[Pi]/n at \[Rho] = 1.*)
(*See Page r4-10.*)


(* Compute A_i using the bisection algorithm *)
aInfl = Module[
 {dest, nValues,
  aMin, aMax,
  a, num,
  aAss, numAss
 },
  (* (This is not slow, nevertheless compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "polygon-a-inflection.txt";
  nValues = Range[3, 5];
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {aAss, numAss} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aAss = Association[];
    numAss = Association[];
    aMin = 15/10;
    aMax = 25/10;
    (* NOTE: these may not work for n > 5 *)
    Table[
      {a, num} = SeekRootBisection[
        curCandCorner[n][#] &,
        {aMin, aMax},
        "ReturnIterations" -> True
      ];
      aAss = aAss // Append[n -> a];
      numAss = numAss // Append[n -> num];
    , {n, nValues}];
    {aAss, numAss} // Compress // Ex[dest]
  ];
  (* Print iterations used *)
  Table[
    Print[
      "n == {n} bisection algorithm: {num} iterations"
        // PrettyString[
          "{n}" -> ToString[n],
          "{num}" -> ToString @ numAss[n]
        ]
    ];
  , {n, nValues}];
  (* Return A_i values *)
  aAss
];
aInfl // N


(* ::Subsection:: *)
(*Representative values of A for a convex domain*)


aConvex[n_] := Ceiling[aInfl[n] + 0.25, 0.1];
Table[
  aConvex[n] < aNat[n][0]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Representative values for offset version*)


(* ::Subsubsection:: *)
(*n and \[Gamma]*)


nOffset = 3;
gammaOffset = 1;


(* ::Subsubsection:: *)
(*A = 1.5 (joined)*)


(* ::Text:: *)
(*Non-viable neighbourhoods at \[Rho] = 1, \[CurlyPhi] = 2\[Pi]k/n*)
(*are still joined to the main non-viable moat.*)


aOffsetJoined = 15/10;


(* ::Text:: *)
(*\[Rho]_\[Sharp] is the largest critical terminal \[Rho] along \[CurlyPhi] = 0.*)


rhoOffsetJoinedSharp =
  Module[{n, gamma, a},
    n = nOffset;
    gamma = gammaOffset;
    a = aOffsetJoined;
    (* Compute \[Rho]_\[Sharp] *)
    SeekRoot[viOffset[gamma][n][a], {Exp[gamma], 1}]
  ];
rhoOffsetJoinedSharp > 1


(* ::Subsubsection:: *)
(*A = 1.7 (split)*)


(* ::Text:: *)
(*Non-viable neighbourhoods at \[Rho] = 1, \[CurlyPhi] = 2\[Pi]k/n pincer off into lakes.*)
(*Along \[CurlyPhi] = 0:*)
(*  \[Rho]_a is on the tip of the moat*)
(*  \[Rho]_b is on the tip of the lake which is nearest to the origin*)
(*  \[Rho]_\[Sharp] is on the tip of the lake which is furthest from the origin*)


aOffsetSplit = 17/10;


(* ::Text:: *)
(*\[Rho]_\[Sharp] is the largest critical terminal \[Rho] along \[CurlyPhi] = 0.*)


rhoOffsetSplitSharp =
  Module[{n, gamma, a},
    n = nOffset;
    gamma = gammaOffset;
    a = aOffsetSplit;
    (* Compute \[Rho]_\[Sharp] *)
    SeekRoot[viOffset[gamma][n][a], {Exp[gamma], 1}]
  ];


(* ::Text:: *)
(*\[Rho]_b is the 2nd largest critical terminal \[Rho] along \[CurlyPhi] = 0.*)


rhoOffsetSplitB =
  Module[{n, gamma, a},
    n = nOffset;
    gamma = gammaOffset;
    a = aOffsetSplit;
    (* Compute \[Rho]_\[Sharp] *)
    SeekRoot[viOffset[gamma][n][a], {1, 8/10}]
  ];


(* ::Text:: *)
(*\[Rho]_a is the 3rd largest critical terminal \[Rho] along \[CurlyPhi] = 0.*)


rhoOffsetSplitA =
  Module[{n, gamma, a},
    n = nOffset;
    gamma = gammaOffset;
    a = aOffsetSplit;
    (* Compute \[Rho]_\[Sharp] *)
    SeekRoot[viOffset[gamma][n][a], {8/10, 1/2}]
  ];


(* ::Text:: *)
(*Check:*)


rhoOffsetSplitA < rhoOffsetSplitB < 1 < rhoOffsetSplitSharp


(* ::Subsection:: *)
(*Starting points for boundary tracing (offset version)*)


(* ::Text:: *)
(*We choose points with \[CurlyPhi] between 0 and 2 \[Pi] / n (since there is n-fold symmetry).*)


(* ::Subsubsection:: *)
(*A = 1.5 (joined)*)


(* Starting points along terminal curve *)
(* (i.e. outer boundary of non-viable domain) *)
startZetaOffsetJoined["terminal"] =
  Module[
   {n, gamma, a, rhoSh,
    phMax, phSpacing, phValues,
    rho
   },
    n = nOffset;
    gamma = gammaOffset;
    a = aOffsetJoined;
    rhoSh = rhoOffsetJoinedSharp;
    (* \[CurlyPhi] values *)
    phMax = 2 Pi / n;
    phSpacing = 30 Degree;
    phValues = UniformRange[0, phMax, phSpacing] // Most;
    (* Build list of values *)
    Table[
      rho = SeekRoot[
        viOffset[gamma][n][a][# Exp[I ph]] &,
        {rhoSh, 1/2}
      ];
      rho Exp[I ph]
    , {ph, phValues}]
  ];


(* Starting point which is hyperbolic critical terminal point *)
startZetaOffsetJoined["hyperbolic"] = {rhoOffsetJoinedSharp};


(* ::Subsubsection:: *)
(*A = 1.7 (split)*)


(* Starting points along moat outer terminal curve *)
startZetaOffsetSplit["terminal-moat"] =
  Module[
   {n, gamma, a, rhoSh, rhoB, rhoA,
    phMax, phSpacing, phValues,
    rho
   },
    n = nOffset;
    gamma = gammaOffset;
    a = aOffsetSplit;
    rhoSh = rhoOffsetSplitSharp;
    rhoB = rhoOffsetSplitB;
    rhoA = rhoOffsetSplitA;
    (* \[CurlyPhi] values *)
    phMax = 2 Pi / n;
    phSpacing = 30 Degree;
    phValues = UniformRange[0, phMax, phSpacing] // Most;
    (* Build list of values *)
    Table[
      rho = SeekRoot[
        viOffset[gamma][n][a][# Exp[I ph]] &,
        {rhoA, 1/2}
      ];
      rho Exp[I ph]
    , {ph, phValues}]
  ];


(* Starting points along lake terminal curve *)
startZetaOffsetSplit["terminal-lake"] =
  Module[
   {n, gamma, a, rhoSh, rhoB, rhoA,
    rhoBetween, phBetween
   },
    n = nOffset;
    gamma = gammaOffset;
    a = aOffsetSplit;
    rhoSh = rhoOffsetSplitSharp;
    rhoB = rhoOffsetSplitB;
    rhoA = rhoOffsetSplitA;
    (* Build list of values *)
    {
      (* Lake point furthest from origin *)
      rhoSh,
      (* Lake point in-between *)
      rhoBetween = Way[rhoB, rhoSh];
      phBetween = SeekRoot[
        viOffset[gamma][n][a][rhoBetween Exp[I #]] &,
        {0, Pi / n}
      ];
      rhoBetween Exp[I phBetween],
      (* Lake point closest to origin *)
      rhoB
    }
  ];


(* Starting point along moat outer terminal curve which is hyperbolic critical *)
startZetaOffsetSplit["hyperbolic-moat"] = {rhoOffsetSplitA};


(* Starting points along lake terminal curve which are hyperbolic critical *)
startZetaOffsetSplit["hyperbolic-lake"] = {
  rhoOffsetSplitB,
  rhoOffsetSplitSharp
};


(* ::Subsection:: *)
(*Traced boundaries \[Zeta] = \[Zeta](s) (offset version)*)


(* ::Subsubsection:: *)
(*\[Zeta]' = d\[Zeta]/ds*)


(* ::Text:: *)
(*See (r4.38) (Page r4-8).*)


zetaVelOffset[gamma_][n_][a_][zeta_] :=
  Divide[
    I fOffset[gamma][a][zeta] + Sqrt @ viOffset[gamma][n][a][zeta],
    gDerZeta[zeta]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Solver for traced boundaries*)


(* ::Text:: *)
(*See (r4.38) (Page r4-8).*)


With[{zeta = \[FormalZeta]},
  zetaTraceOffset[gamma_][n_][a_][zetaInit_] :=
    NDSolveValue[
      {
        zeta'[s] == zetaVelOffset[gamma][n][a] @ zeta[s],
        zeta[0] == zetaInit,
        WhenEvent[
          Or[
            viOffset[gamma][n][a] @ zeta[s] < -10^-6,
            Abs @ zeta[s] > Exp[gamma]
          ],
          "StopIntegration"
        ]
      }, zeta, {s, -Pi Exp[gamma], Pi Exp[gamma]},
      NoExtrapolation
    ]
];


(* ::Subsubsection:: *)
(*A = 1.5 (joined)*)


Module[{n, gamma, a, rhoSharp, idList, zetaInitList},
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetJoined;
  (* Group names *)
  idList = {"terminal", "hyperbolic"};
  (* Solve for traced boundaries *)
  Table[
    zetaInitList = startZetaOffsetJoined[id];
    zetaTraOffsetJoined[id] =
      Table[
        zetaTraceOffset[gamma][n][a][zetaInit]
      , {zetaInit, zetaInitList}];
  , {id, idList}];
];


(* ::Subsubsection:: *)
(*A = 1.7 (split)*)


Module[{n, gamma, a, rhoSharp, idList, zetaInitList},
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetSplit;
  (* Group names *)
  idList = {
    "terminal-moat",
    "terminal-lake",
    "hyperbolic-moat",
    "hyperbolic-lake"
  };
  (* Solve for traced boundaries *)
  Table[
    zetaInitList = startZetaOffsetSplit[id];
    zetaTraOffsetSplit[id] =
      Table[
        zetaTraceOffset[gamma][n][a][zetaInit]
      , {zetaInit, zetaInitList}];
  , {id, idList}];
];


(* ::Subsection:: *)
(*Traced boundary curvature (offset version)*)


(* ::Subsubsection:: *)
(*z' = dz/ds*)


(* ::Text:: *)
(*See (r4.39) (Page r4-11).*)


zVelOffset[gamma_][n_][a_][zeta_] :=
  Divide[
    I fOffset[gamma][a][zeta] + Sqrt @ viOffset[gamma][n][a][zeta],
    gDer[n][zeta]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*z'' = d^2(z)/ds^2*)


(* ::Text:: *)
(*See (r4.40) to (r4.43) (Pages r4-11 & r4-12).*)


zAccOffset[gamma_][n_][a_][zeta_] := Plus[
  (* 1st term of (r4.40) *)
  Times[
    (* i F + sqrt(\[CapitalPhi]) *)
    I fOffset[gamma][a][zeta] + Sqrt @ viOffset[gamma][n][a][zeta],
    (* d\[Zeta]/ds *)
    zetaVelOffset[gamma][n][a][zeta],
    (* d/d\[Zeta] (1 / (dG/dz)) *)
    D[1 / gDer[n][zeta], zeta]
  ],
  (* 2nd term of (r4.40) *)
  Divide[
    Plus[
      (* i F' *)
      Times[
        I,
        -4 (Re @ gOffset[gamma][zeta])^3 / a,
        Re @ (zetaVelOffset[gamma][n][a][zeta] gDerZeta[zeta])
      ],
      (* (sqrt(\[CapitalPhi]))' == \[CapitalPhi]' / (2 sqrt(\[CapitalPhi]))  *)
      Divide[
        Plus[
          (* 1st term of (r4.43) *)
          Times[
            2,
            Re @ gDer[n][zeta],
            Re @ (zetaVelOffset[gamma][n][a][zeta] D[gDer[n][zeta], zeta])
          ],
          (* 2nd term of (r4.43) *)
          Times[
            2,
            Im @ gDer[n][zeta],
            Im @ (zetaVelOffset[gamma][n][a][zeta] D[gDer[n][zeta], zeta])
          ],
          (* 3rd term of (r4.43) *)
          Times[
            -8 (Re @ gOffset[gamma][zeta])^7 / a^2,
            Re @ (zetaVelOffset[gamma][n][a][zeta] gDerZeta[zeta])
          ]
        ],
        2 Sqrt @ viOffset[gamma][n][a][zeta]
      ]
    ],
    (* dG/dz *)
    gDer[n][zeta]
  ]
] // Evaluate;


(* ::Subsubsection:: *)
(*Curvature x' y'' - y' x''*)


(* ::Text:: *)
(*See (r4.39a) (Page r4-11).*)


curTraOffset[gamma_][n_][a_][zeta_] :=
  Module[{xVel, yVel, xAcc, yAcc},
    {xVel, yVel} = zVelOffset[gamma][n][a][zeta] // ReIm;
    {xAcc, yAcc} = zAccOffset[gamma][n][a][zeta] // ReIm;
    xVel yAcc - yVel xAcc
  ] // Evaluate;


(* ::Subsection:: *)
(*Representative values for convex offset version*)


(* ::Text:: *)
(*Chosen so that convex domains result.*)


(* ::Subsubsection:: *)
(*n, \[Gamma] and A*)


nOffsetConvex = 3;
gammaOffsetConvex = 1.6;
aOffsetConvex = 12;


(* ::Subsubsection:: *)
(*\[Rho]_\[Sharp], \[Rho]_b and \[Rho]_a*)


(* ::Text:: *)
(*Non-viable neighbourhoods at \[Rho] = 1, \[CurlyPhi] = 2\[Pi]k/n pincer off into lakes.*)
(*Along \[CurlyPhi] = 0:*)
(*  \[Rho]_a is on the tip of the moat*)
(*  \[Rho]_b is on the tip of the lake which is nearest to the origin*)
(*  \[Rho]_\[Sharp] is on the tip of the lake which is furthest from the origin*)


(* ::Text:: *)
(*\[Rho]_\[Sharp] is the largest critical terminal \[Rho] along \[CurlyPhi] = 0.*)


rhoOffsetConvexSharp =
  Module[{n, gamma, a},
    n = nOffsetConvex;
    gamma = gammaOffsetConvex;
    a = aOffsetConvex;
    (* Compute \[Rho]_\[Sharp] *)
    SeekRoot[viOffset[gamma][n][a], {Exp[gamma], 1}]
  ];


(* ::Text:: *)
(*\[Rho]_b is the 2nd largest critical terminal \[Rho] along \[CurlyPhi] = 0.*)


rhoOffsetConvexB =
  Module[{n, gamma, a},
    n = nOffsetConvex;
    gamma = gammaOffsetConvex;
    a = aOffsetConvex;
    (* Compute \[Rho]_\[Sharp] *)
    SeekRoot[viOffset[gamma][n][a], {1, 1/2}]
  ];


(* ::Text:: *)
(*\[Rho]_a is the 3rd largest critical terminal \[Rho] along \[CurlyPhi] = 0.*)


rhoOffsetConvexA =
  Module[{n, gamma, a},
    n = nOffsetConvex;
    gamma = gammaOffsetConvex;
    a = aOffsetConvex;
    (* Compute \[Rho]_\[Sharp] *)
    SeekRoot[viOffset[gamma][n][a], {1/2, 1/5}]
  ];


(* ::Text:: *)
(*Check:*)


rhoOffsetConvexA < rhoOffsetConvexB < 1 < rhoOffsetConvexSharp


(* ::Subsection:: *)
(*Starting points for boundary tracing (convex offset version)*)


(* ::Text:: *)
(*We choose points with \[CurlyPhi] between 0 and 2 \[Pi] / n (since there is n-fold symmetry).*)


(* Starting points which are not terminal *)
startZetaOffsetConvex["general"] =
  Module[
   {n, gamma, a,
    rhoA, rhoB, rho,
    phMax, phSpacing, phValues
   },
    n = nOffsetConvex;
    gamma = gammaOffsetConvex;
    a = aOffsetConvex;
    (* \[Rho] *)
    rhoB = rhoOffsetConvexB;
    rhoA = rhoOffsetConvexA;
    rho = Way[rhoA, rhoB];
    (* \[CurlyPhi] values *)
    phMax = 2 Pi / n;
    phSpacing = 30 Degree;
    phValues = UniformRange[0, phMax, phSpacing] // Most;
    (* Build list of values *)
    Table[rho Exp[I ph], {ph, phValues}]
  ];


(* Starting points along moat outer terminal curve *)
startZetaOffsetConvex["terminal-moat"] =
  Module[
   {n, gamma, a, rhoSh, rhoB, rhoA,
    phMax, phSpacing, phValues,
    rho
   },
    n = nOffsetConvex;
    gamma = gammaOffsetConvex;
    a = aOffsetConvex;
    rhoSh = rhoOffsetConvexSharp;
    rhoB = rhoOffsetConvexB;
    rhoA = rhoOffsetConvexA;
    (* \[CurlyPhi] values *)
    phMax = 2 Pi / n;
    phSpacing = 30 Degree;
    phValues = UniformRange[0, phMax, phSpacing] // Most;
    (* Build list of values *)
    Table[
      rho = SeekRoot[
        viOffset[gamma][n][a][# Exp[I ph]] &,
        {rhoA, 1/2}
      ];
      rho Exp[I ph]
    , {ph, phValues}]
  ];


(* Starting points along lake terminal curve *)
startZetaOffsetConvex["terminal-lake"] =
  Module[
   {n, gamma, a, rhoSh, rhoB, rhoA,
    rhoBetween, phBetween
   },
    n = nOffsetConvex;
    gamma = gammaOffsetConvex;
    a = aOffsetConvex;
    rhoSh = rhoOffsetConvexSharp;
    rhoB = rhoOffsetConvexB;
    rhoA = rhoOffsetConvexA;
    (* Build list of values *)
    {
      (* Lake point furthest from origin *)
      rhoSh,
      (* Lake point in-between *)
      rhoBetween = Way[rhoB, rhoSh];
      phBetween = SeekRoot[
        viOffset[gamma][n][a][rhoBetween Exp[I #]] &,
        {0, Pi / n}
      ];
      rhoBetween Exp[I phBetween],
      (* Lake point closest to origin *)
      rhoB
    }
  ];


(* Starting point along moat outer terminal curve which is hyperbolic critical *)
startZetaOffsetConvex["hyperbolic-moat"] = {rhoOffsetConvexA};


(* Starting points along lake terminal curve which are hyperbolic critical *)
startZetaOffsetConvex["hyperbolic-lake"] = {
  rhoOffsetConvexB,
  rhoOffsetConvexSharp
};


(* ::Subsection:: *)
(*Traced boundaries \[Zeta] = \[Zeta](s) (convex offset version)*)


Module[{n, gamma, a, rhoSharp, idList, zetaInitList},
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  a = aOffsetConvex;
  (* Group names *)
  idList = {
    "general",
    "terminal-moat",
    "terminal-lake",
    "hyperbolic-moat",
    "hyperbolic-lake"
  };
  (* Solve for traced boundaries *)
  Table[
    zetaInitList = startZetaOffsetConvex[id];
    zetaTraOffsetConvex[id] =
      Table[
        zetaTraceOffset[gamma][n][a][zetaInit]
      , {zetaInit, zetaInitList}];
  , {id, idList}];
];


(* ::Subsection:: *)
(*Convex offset joining arc lengths*)


(* ::Text:: *)
(*Values of s for constructing the "rotund" and "elongated" convex offset domains.*)
(*See:*)
(*  polygon_offset_z-convex-traced-hyperbolic_moat.pdf*)
(*  polygon_offset_z-convex-traced-rotund.pdf*)
(*  polygon_offset_z-convex-traced-elongated.pdf*)


(* (Assumes traced boundaries does not make a full turn) *)
Module[{n, zeta, fun, sInterval},
  n = nOffsetConvex;
  zeta = zetaTraOffsetConvex["hyperbolic-moat"] // First;
  fun[phi_] := Arg @ Conjugate @ zeta[#] - phi &;
  sInterval = {0, DomainStart[zeta]};
  (* Rotund, or disk-like (like a rounded regular polygon) *)
  sOffsetConvex["rotund"] = SeekRoot[fun[Pi / n], sInterval];
  (* Elongated, or lemon-like (like an eye) *)
  sOffsetConvex["elongated"] = SeekRoot[fun[2 Pi / n], sInterval];
];


(* ::Subsection:: *)
(*Numerical verification for convex domains (finite elements)*)


(* ::Text:: *)
(*These are not slow, nevertheless compute once and store.*)
(*Delete the corresponding file manually to compute from scratch.*)


(* ::Subsubsection:: *)
(*Generate meshes*)


Table[
  Module[{dest},
    dest = "polygon-convex-verification-mesh-" <> ToString[n] <> ".txt";
    If[Not @ FileExistsQ[dest],
      Module[
       {a, rhoSh, rSh, zeta, sMax, rhoBath, rBath, tBath,
        sSpacing, rotate, sValues,
        extPointList, intPointList,
        nExt, nInt, mod,
        bMesh, mesh,
        prExt, prInt
       },
        (* Dimensionless group *)
        a = aConvex[n];
        (* Hyperbolic critical terminal radius (\[Zeta]-space and z-space) *)
        rhoSh = rhoSharp[n][a][0];
        rSh = rhoSh // zMap[n];
        (* Traced boundary \[Zeta] == \[Zeta](s) *)
        zeta = zetaTraceCand[n][a];
        (* Half of arc length along one candidate boundary *)
        (* (which is one of the n boundaries of the domain) *)
        sMax = DomainStart[zeta] // Abs;
        (* Heat bath radius (\[Zeta]-space and z-space) *)
        rhoBath = 0.5 rhoSh;
        rBath = rhoBath // zMap[n];
        (* Heat bath temperature *)
        tBath = -Log @ rhoBath;
        (* Spacing of boundary points *)
        (* (roughly 1 boundary point every 5 degrees along rho == rhoSh) *)
        sSpacing = rhoSh * 5 Degree;
        (* External (radiation) boundary points *)
        sValues = UniformRange[-sMax, 0, sSpacing];
        rotate[k_] := Exp[I 2 Pi k / n] # &;
        extPointList = Join @@ Table[
          Join[
            Table[
              zeta[-Abs @ s] // zMap[n] // rotate[k] // ReIm
            , {s, sValues}] // Most,
            Table[
              zeta[-Abs @ s] // zMap[n] // Conjugate // rotate[k] // ReIm
            , {s, sValues}] // Reverse // Most
          ]
        , {k, 0, n - 1}];
        (* Internal (heat bath) boundary *)
        intPointList = Table[
          rhoBath Exp[I ph] // zMap[n] // ReIm
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
        prExt = Function[{x, y}, RPolar[x, y] > Way[rBath, rSh] // Evaluate];
        prInt = Function[{x, y}, RPolar[x, y] < Way[rBath, rSh] // Evaluate];
        {a, tBath, mesh, prExt, prInt}
          // Compress // Ex[dest]
      ]
    ]
  ]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Solve BVPs*)


Table[
  Module[{dest},
    dest = "polygon-convex-verification-solution-" <> ToString[n] <> ".txt";
    If[Not @ FileExistsQ[dest],
      Module[
       {source,
        a, tBath, mesh, prExt, prInt,
        tSol
       },
        (* Import mesh *)
        source = "polygon-convex-verification-mesh-" <> ToString[n] <> ".txt";
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
, {n, 3, 5}]


(* ::Subsection:: *)
(*Numerical verification for convex offset domains (finite elements)*)


(* ::Text:: *)
(*These are not slow, nevertheless compute once and store.*)
(*Delete the corresponding file manually to compute from scratch.*)


(* ::Subsubsection:: *)
(*Mesh (rotund)*)


Module[{dest},
  dest = "polygon_offset-convex-verification-mesh-rotund.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {n, gamma, a,
      rhoA, rA, rhoBath, rBath, tBath,
      zeta, sSpacing, sValues,
      extPointList, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      n = nOffsetConvex;
      gamma = gammaOffsetConvex;
      a = aOffsetConvex;
      (* Hyperbolic critical terminal radius (tip of moat) *)
      (* (both \[Zeta]-space and z-space) *)
      rhoA = rhoOffsetConvexA;
      rA = rhoA // zMap[n];
      (* Heat bath radius (\[Zeta]-space and z-space) *)
      rhoBath = 0.5 rhoA;
      rBath = rhoBath // zMap[n];
      (* Heat bath temperature *)
      tBath = gOffset[gamma][rhoBath];
      (* Traced boundary \[Zeta] == \[Zeta](s) *)
      zeta = zetaTraOffsetConvex["hyperbolic-moat"] // First;
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 5 degrees along r == rA) *)
      sSpacing = rA * 5 Degree;
      (* External (radiation) boundary points *)
      (*
        I am lazy, so I am simply listing the boundary points from
        "polygon_offset_z-convex-traced-rotund.pdf",
        then removing duplicates and sorting them by azimuthal angle.
       *)
      sValues = UniformRange[sOffsetConvex["rotund"], 0, sSpacing];
      extPointList = (
        Join @@ Table[
          Join @@ Table[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
          , {s, sValues}]
        , {k, 0, n - 1}]
          // DeleteNearbyPoints[sSpacing / 4]
          // SortByPhi
      );
      (* Internal (heat bath) boundary points *)
      intPointList = Table[
        rhoBath Exp[I ph] // zMap[n] // ReIm
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
      prExt = Function[{x, y}, RPolar[x, y] > Way[rBath, rA] // Evaluate];
      prInt = Function[{x, y}, RPolar[x, y] < Way[rBath, rA] // Evaluate];
      {a, tBath, mesh, prExt, prInt}
        // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Mesh (elongated)*)


Module[{dest},
  dest = "polygon_offset-convex-verification-mesh-elongated.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {n, gamma, a,
      rhoA, rA, rhoBath, rBath, tBath,
      zeta, sSpacing, k, sValues,
      extPointList, intPointList,
      nExt, nInt, mod,
      bMesh, mesh,
      prExt, prInt
     },
      n = nOffsetConvex;
      gamma = gammaOffsetConvex;
      a = aOffsetConvex;
      (* Hyperbolic critical terminal radius (tip of moat) *)
      (* (both \[Zeta]-space and z-space) *)
      rhoA = rhoOffsetConvexA;
      rA = rhoA // zMap[n];
      (* Heat bath radius (\[Zeta]-space and z-space) *)
      rhoBath = 0.5 rhoA;
      rBath = rhoBath // zMap[n];
      (* Heat bath temperature *)
      tBath = gOffset[gamma][rhoBath];
      (* Traced boundary \[Zeta] == \[Zeta](s) *)
      zeta = zetaTraOffsetConvex["hyperbolic-moat"] // First;
      (* Spacing of boundary points *)
      (* (roughly 1 boundary point every 5 degrees along r == rA) *)
      sSpacing = rA * 5 Degree;
      (* External (radiation) boundary points *)
      (*
        I am lazy, so I am simply listing the boundary points from
        "polygon_offset_z-convex-traced-elongated.pdf",
        then removing duplicates and sorting them by azimuthal angle.
       *)
      extPointList = (
        Join[
          (* Boundaries to the right *)
          k = 1;
          sValues = UniformRange[sOffsetConvex["elongated"], 0, sSpacing];
          Join @@ Table[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
          , {s, sValues}],
          (* Boundaries to the left *)
          k = 2;
          sValues = UniformRange[sOffsetConvex["rotund"], 0, sSpacing];
          Join @@ Table[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
          , {s, sValues}]
        ]
          // DeleteNearbyPoints[sSpacing / 4]
          // SortByPhi
      );
      (* Internal (heat bath) boundary points *)
      intPointList = Table[
        rhoBath Exp[I ph] // zMap[n] // ReIm
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
      prExt = Function[{x, y}, RPolar[x, y] > Way[rBath, rA] // Evaluate];
      prInt = Function[{x, y}, RPolar[x, y] < Way[rBath, rA] // Evaluate];
      {a, tBath, mesh, prExt, prInt}
        // Compress // Ex[dest]
    ]
  ]
]


(* ::Subsubsection:: *)
(*Solve BVP (rotund)*)


Module[{dest},
  dest = "polygon_offset-convex-verification-solution-rotund.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, tBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "polygon_offset-convex-verification-mesh-rotund.txt";
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


(* ::Subsubsection:: *)
(*Solve BVP (elongated)*)


Module[{dest},
  dest = "polygon_offset-convex-verification-solution-elongated.txt";
  If[Not @ FileExistsQ[dest],
    Module[
     {source,
      a, tBath, mesh, prExt, prInt,
      tSol
     },
      (* Import mesh *)
      source = "polygon_offset-convex-verification-mesh-elongated.txt";
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
(*Geometric regions*)


(* ::Subsubsection:: *)
(*Unit disk (\[Zeta]-space)*)


diskPoints[num_Integer] := diskPoints[num] = (
  SeedRandom[0];
  RandomPoint[Disk[], num]
);


(* ::Subsubsection:: *)
(*Regular n-gon (z-space)*)


poly[n_Integer] := RegularPolygon[{1, 0}, n];


polyVertices[n_Integer] := CirclePoints[{1, 0}, n];


(* ::Subsubsection:: *)
(*Complement of regular n-gon (z-space)*)


(* ::Text:: *)
(*This is quicker to plot than using RegionDifference with FullRegion.*)


polyComplement[n_Integer] :=
  Module[{mid},
    Table[
      (* Midpoint of k-th edge *)
      mid = Way[
        AngleVector[2 Pi (k - 1) / n],
        AngleVector[2 Pi k / n]
      ];
      (* Corresponding half space *)
      HalfSpace[-mid, mid]
    , {k, n}]
  ];


(* ::Subsection:: *)
(*Italicise symbols*)


aIt = Italicise["A"];
gIt = Style["\[Gamma]"];
nIt = Italicise["n"];
sIt = Italicise["s"];
zIt = Italicise["z"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;


contStyle = LightGray;
streamStyle = LightGray;
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


psiStyle = Blue;
aStyle = Purple;


critStyle = Red;


ellStyle = Darker[Green];
hypStyle = Magenta;


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


guideStyle = Dashed;
pointStyle = PointSize[Large];
glowStyle = Directive[Thick, Yellow, Opacity[0.7]];


(* ::Subsection:: *)
(*Repeated plots*)


(* ::Subsubsection:: *)
(*Equipotentials and streamlines (z-space)*)


equipStream[n_, opts : OptionsPattern[RegionPlot]] :=
  Module[{rMax, rhoNum, rhoValues, phNum, phValues},
    rMax = 1;
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax}, opts],
      (* Equipotentials (T == const) *)
      rhoNum = 4;
      rhoValues = Subdivide[0, 1, rhoNum] // Most;
      ParametricPlot[
        Table[
          rho Exp[I ph] // zMap[n] // ReIm
        , {rho, rhoValues}] // Evaluate,
        {ph, 0, 2 Pi},
        PlotStyle -> contStyle
      ],
      (* Streamlines (parallel to \[Del]T) *)
      phNum = 4;
      phValues = Subdivide[0, 2 Pi, n phNum] // Most;
      ParametricPlot[
        Table[
          rho Exp[I ph] // zMap[n] // ReIm
        , {ph, phValues}] // Evaluate,
        {rho, 0, 1},
        PlotStyle -> streamStyle
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*Unphysical domain (z-space)*)


unphysDomain[n_] := Graphics @ {unphysStyle, polyComplement[n]};


(* ::Section:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsection:: *)
(*Forward transformation: disk (\[Zeta]-space) to polygon (z-space)*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/159086*)


Table[
  Show[
    (* Equipotentials and streamlines *)
    equipStream[n],
    (* Unphysical domain *)
    unphysDomain[n]
  ] // Ex @ StringJoin["schwarz-christoffel-forward-", ToString[n], ".pdf"]
, {n, 3, 7}]


(* ::Section:: *)
(*\[Psi]*)


(* ::Subsection:: *)
(*Algebra*)


With[{n = \[FormalN], rho = \[FormalRho], ph = \[FormalCurlyPhi]},
  Block[{$Assumptions = n > 2 && 0 < rho < 1 && -Pi < ph < Pi},
    Module[{psiExpr, psiDerExpr},
      psiExpr = psi[n][rho Exp[I ph]] // ComplexExpand // FullSimplify;
      psiDerExpr = psiExpr // D[#, rho] & // FullSimplify;
      {
        {"psi", "(r4.29)", psiExpr},
        {"dpsi/dr", "(r4.30)", psiDerExpr}
      } /. {Gamma -> "\[CapitalGamma]"}
        // PrettyString["d" -> "\[PartialD]", "psi" -> "\[Psi]"]
        // TableForm
    ]
  ]
] // Ex["polygon-psi-algebra.pdf"]


(* ::Subsection:: *)
(*Algebra (offset version)*)


With[
 {n = \[FormalN],
  rho = \[FormalRho], ph = \[FormalCurlyPhi],
  gamma = \[FormalGamma],
  xi = \[FormalXi]
 },
  Block[{$Assumptions = n > 2 && 0 < rho < Exp[gamma] && -Pi < ph < Pi},
    Module[{psiExpr, pert},
      psiExpr = (
        psiOffset[gamma][n][rho Exp[I ph]]
          // ComplexExpand
          // FullSimplify
      );
      pert = (# /. {rho -> 1 + xi}) + O[xi]^2 &;
      {
        {"psi", "(r4.54)", psiExpr},
        {Log["rho"], "(r4.56)", Log[rho] // pert},
        {"rho" ^ "n", "(r4.56)", rho ^ n // pert}
      } /. {Gamma -> "\[CapitalGamma]"}
        // PrettyString[
          "d" -> "\[PartialD]",
          "psi" -> "\[Psi]",
          "gamma" -> "\[Gamma]",
          "rho" -> "\[Rho]"
        ]
        // TableForm
    ]
  ]
] // Ex["polygon_offset-psi-algebra.pdf"]


(* ::Subsection:: *)
(*\[Psi] plot*)


Module[{phValues},
  Table[
    (* \[CurlyPhi] values *)
    phValues = Subdivide[0, Pi / n, 4];
    (* Plot *)
    Plot[
      Table[
        psi[n][rho Exp[I ph]]
      , {ph, phValues}] // Evaluate,
      {rho, 0, 1},
      AxesLabel -> {"rho", "psi"},
      ImageSize -> 360,
      PlotLabel -> nIt == n,
      PlotLegends -> LineLegend[
        phValues,
        LegendLabel -> "ph"
      ],
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ] // PrettyString[
      "rho" -> "\[Rho]",
      "ph" -> "\[CurlyPhi]",
      "psi" -> "\[Psi]"
    ] // Ex @ StringJoin["polygon-psi-", ToString[n], ".pdf"]
  , {n, 3, 5}]
]


(* ::Subsection:: *)
(*\[Psi] plot (offset version)*)


Module[{gammaValues, phValues, rhoMax},
  (* \[Gamma] values *)
  gammaValues = {1, 2, 4};
  Table[
    (* \[CurlyPhi] values *)
    phValues = Subdivide[0, Pi / n, 4];
    (* Plot *)
    rhoMax = Min[Exp[gamma], 3];
    Plot[
      Table[
        psiOffset[gamma][n][rho Exp[I ph]]
      , {ph, phValues}] // Evaluate,
      {rho, 0, rhoMax},
      AxesLabel -> {"rho", "psi"},
      ImageSize -> 360,
      PlotLabel -> Row[{nIt == n, gIt == gamma}, ","],
      PlotLegends -> LineLegend[
        phValues,
        LegendLabel -> "ph"
      ],
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ] // PrettyString[
      "rho" -> "\[Rho]",
      "ph" -> "\[CurlyPhi]",
      "psi" -> "\[Psi]"
    ] // Ex @ StringJoin[
      "polygon_offset-psi-", ToString[n],
      "-gamma-", ToString[gamma],
      ".pdf"
    ]
  , {n, 3, 5}, {gamma, gammaValues}]
]


(* ::Subsubsection:: *)
(*A = 1.5 (joined)*)


Module[{n, gamma, a, rhoSh, rhoMax},
  (* Values of n, \[Gamma] and A *)
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetJoined;
  (* \[Rho]_\[Sharp] *)
  rhoSh = rhoOffsetJoinedSharp;
  (* Plot *)
  rhoMax = 2;
  Show[
    Plot[
      {psiOffset[gamma][n][rho], a} // Evaluate,
      {rho, 0, rhoMax},
      AxesLabel -> {"rho", "psi"},
      PlotRange -> {0, 7.5},
      PlotLabel -> Row[
        {nIt == n, gIt == gamma, aIt == N[a]},
        ","
      ],
      PlotStyle -> {psiStyle, aStyle},
      PlotOptions[Axes] // Evaluate
    ],
    (* A *)
    Graphics @ {aStyle,
      Text[aIt // textStyle, {0, a}, {3, 0}]
    },
    (* \[Rho]_\[Sharp] *)
    Graphics @ {Directive[critStyle, pointStyle],
      Point @ {rhoSh, a}
    },
    Graphics @ {Directive[critStyle, guideStyle],
      Line @ {{rhoSh, a}, {rhoSh, 0}}
    },
    Graphics @ {critStyle,
      Text[
        Subscript["rho", "sharp"] // textStyle,
        {rhoSh, 0},
        {0, 2.2}
      ]
    },
    (* Plot range *)
    PlotRange -> All,
    PlotRangeClipping -> False
  ] // PrettyString[
    "rho" -> "\[Rho]",
    "psi" -> "\[Psi]",
    "sharp" -> "\[Sharp]"
  ]
] // Ex["polygon_offset-psi-joined.pdf"]


(* ::Subsubsection:: *)
(*A = 1.7 (split)*)


Module[{n, gamma, a, rhoSh, rhoB, rhoA, rhoMax},
  (* Values of n, \[Gamma] and A *)
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetSplit;
  (* \[Rho]_\[Sharp], \[Rho]_b and \[Rho]_a *)
  rhoSh = rhoOffsetSplitSharp;
  rhoB = rhoOffsetSplitB;
  rhoA = rhoOffsetSplitA;
  (* Plot *)
  rhoMax = 2;
  Show[
    Plot[
      {psiOffset[gamma][n][rho], a} // Evaluate,
      {rho, 0, rhoMax},
      AxesLabel -> {"rho", "psi"},
      PlotRange -> {0, 7.5},
      PlotLabel -> Row[
        {nIt == n, gIt == gamma, aIt == N[a]},
        ","
      ],
      PlotStyle -> {psiStyle, aStyle},
      PlotOptions[Axes] // Evaluate
    ],
    (* A *)
    Graphics @ {aStyle,
      Text[aIt // textStyle, {0, a}, {3, 0}]
    },
    (* \[Rho]_\[Sharp] *)
    Graphics @ {Directive[critStyle, pointStyle],
      Point @ {rhoSh, a}
    },
    Graphics @ {Directive[critStyle, guideStyle],
      Line @ {{rhoSh, a}, {rhoSh, 0}}
    },
    Graphics @ {critStyle,
      Text[
        Subscript["rho", "sharp"] // textStyle,
        {rhoSh, 0},
        {0, 2.2}
      ]
    },
    (* \[Rho]_b *)
    Graphics @ {Directive[critStyle, pointStyle],
      Point @ {rhoB, a}
    },
    Graphics @ {Directive[critStyle, guideStyle],
      Line @ {{rhoB, a}, {rhoB, 0}}
    },
    Graphics @ {critStyle,
      Text[
        Subscript["rho", "b"] // textStyle,
        {rhoB, 0},
        {0, 2.2}
      ]
    },
    (* \[Rho]_a *)
    Graphics @ {Directive[critStyle, pointStyle],
      Point @ {rhoA, a}
    },
    Graphics @ {Directive[critStyle, guideStyle],
      Line @ {{rhoA, a}, {rhoA, 0}}
    },
    Graphics @ {critStyle,
      Text[
        Subscript["rho", "a"] // textStyle,
        {rhoA, 0},
        {0, 2.2}
      ]
    },
    (* Plot range *)
    PlotRange -> All,
    PlotRangeClipping -> False
  ] // PrettyString[
    "rho" -> "\[Rho]",
    "psi" -> "\[Psi]",
    "sharp" -> "\[Sharp]"
  ]
] // Ex["polygon_offset-psi-split.pdf"]


(* ::Subsection:: *)
(*\[Psi] plot (offset version) interactive visualiser*)


DynamicModule[
 {nValues, nInit,
  gammaMin, gammaMax, gammaInit,
  aMin, aMax, aInit,
  rhoMax
 },
  (* Values of n *)
  nValues = Range[3, 5];
  nInit = First[nValues];
  (* Values of \[Gamma] *)
  gammaMin = 0.1;
  gammaMax = 4;
  gammaInit = 1.6;
  (* Values of A *)
  aMin = 0.01;
  aMax = 100;
  aInit = 12;
  (* Plot range *)
  rhoMax = 2;
  Manipulate[
    Plot[
      {psiOffset[gamma][n][rho], a} // Evaluate,
      {rho, 0, rhoMax},
      AxesLabel -> {"rho", "psi"},
      PlotRange -> Automatic,
      PlotLabel -> Row[
        {nIt == n, gIt == N[gamma], aIt == N[a]},
        ","
      ],
      PlotStyle -> {psiStyle, aStyle},
      PlotOptions[Axes] // Evaluate
    ] // PrettyString[
      "rho" -> "\[Rho]",
      "psi" -> "\[Psi]"
    ]
  , {{n, nInit, nIt}, nValues}
  , {{gamma, gammaInit, gIt}, gammaMin, gammaMax, Appearance -> "Open"}
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}]
]


(* ::Subsection:: *)
(*\[Psi] plot (convex offset version)*)


Module[{n, gamma, a, rhoSh, rhoB, rhoA, rhoMax},
  (* Values of n, \[Gamma] and A *)
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  a = aOffsetConvex;
  (* \[Rho]_\[Sharp], \[Rho]_b and \[Rho]_a *)
  rhoSh = rhoOffsetConvexSharp;
  rhoB = rhoOffsetConvexB;
  rhoA = rhoOffsetConvexA;
  (* Plot *)
  rhoMax = 2;
  Show[
    Plot[
      {psiOffset[gamma][n][rho], a} // Evaluate,
      {rho, 0, rhoMax},
      AxesLabel -> {"rho", "psi"},
      PlotRange -> Automatic,
      PlotLabel -> Row[
        {nIt == n, gIt == N[gamma], aIt == N[a]},
        ","
      ],
      PlotStyle -> {psiStyle, aStyle},
      PlotOptions[Axes] // Evaluate
    ],
    (* A *)
    Graphics @ {aStyle,
      Text[aIt // textStyle, {0, a}, {4.5, 0}]
    },
    (* \[Rho]_\[Sharp] *)
    Graphics @ {Directive[critStyle, pointStyle],
      Point @ {rhoSh, a}
    },
    Graphics @ {Directive[critStyle, guideStyle],
      Line @ {{rhoSh, a}, {rhoSh, 0}}
    },
    Graphics @ {critStyle,
      Text[
        Subscript["rho", "sharp"] // textStyle,
        {rhoSh, 0},
        {0, 2.2}
      ]
    },
    (* \[Rho]_b *)
    Graphics @ {Directive[critStyle, pointStyle],
      Point @ {rhoB, a}
    },
    Graphics @ {Directive[critStyle, guideStyle],
      Line @ {{rhoB, a}, {rhoB, 0}}
    },
    Graphics @ {critStyle,
      Text[
        Subscript["rho", "b"] // textStyle,
        {rhoB, 0},
        {0, 2.2}
      ]
    },
    (* \[Rho]_a *)
    Graphics @ {Directive[critStyle, pointStyle],
      Point @ {rhoA, a}
    },
    Graphics @ {Directive[critStyle, guideStyle],
      Line @ {{rhoA, a}, {rhoA, 0}}
    },
    Graphics @ {critStyle,
      Text[
        Subscript["rho", "a"] // textStyle,
        {rhoA, 0},
        {0, 2.2}
      ]
    },
    (* Plot range *)
    PlotRange -> All,
    PlotRangeClipping -> False
  ] // PrettyString[
    "rho" -> "\[Rho]",
    "psi" -> "\[Psi]",
    "sharp" -> "\[Sharp]"
  ]
] // Ex["polygon_offset-psi-convex.pdf"]


(* ::Subsection:: *)
(*\[Rho]_\[Natural] angular dependence plot*)


Table[
  Plot[rhoNat[n][ph], {ph, 0, 2 Pi},
    AxesLabel -> {"ph", Subscript["rho", "Nat"]},
    ImageSize -> 360,
    PlotLabel -> Column[
      {
        nIt == n,
        Equal[
          Row[{Max, Min}, "/"] - 1,
          rhoNat[n][0] / rhoNat[n][Pi / n] - 1
        ]
      },
      Alignment -> Center
    ],
    PlotRange -> Full,
    PlotOptions[Axes] // Evaluate
  ] // PrettyString[
    "rho" -> "\[Rho]",
    "ph" -> "\[CurlyPhi]",
    "Nat" -> "\[Natural]"
  ] // Ex @ StringJoin["polygon-rho-nat-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Subsection:: *)
(*A_\[Natural] angular dependence plot*)


Table[
  Plot[aNat[n][ph], {ph, 0, 2 Pi},
    AxesLabel -> {"ph", Subscript[aIt, "Nat"]},
    ImageSize -> 360,
    PlotLabel -> Column[
      {
        nIt == n,
        Equal[
          Row[{Max, Min}, "/"] - 1,
          aNat[n][0] / aNat[n][Pi / n] - 1
        ]
      },
      Alignment -> Center
    ],
    PlotRange -> Full,
    PlotOptions[Axes] // Evaluate
  ] // PrettyString[
    "ph" -> "\[CurlyPhi]",
    "Nat" -> "\[Natural]"
  ] // Ex @ StringJoin["polygon-a-nat-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Section:: *)
(*Viable domain*)


(* ::Subsection:: *)
(*Interactive visualiser (\[Zeta]-space)*)


DynamicModule[
 {nValues, nInit,
  ph, aN, aMin, aMax, aInit,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  nValues = Range[3, 5];
  nInit = First[nValues];
  Manipulate[
    (* Azimuthal angle in \[Zeta]-space *)
    ph = 0;
    (* Values of A *)
    aN = aNat[n][ph];
    aMin = 1/1000;
    aMax = aN;
    aInit = aMin // N;
    Manipulate[
      eps = 0.1;
      rhoMax = 1;
      rhoMaxUnphys = 1 + eps;
      rhoMaxNon = If[a < aN, rhoSharp[n][a][ph], rhoNat[n][ph]];
      Show[
        EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
          FrameLabel -> {Re["zeta"], Im["zeta"]},
          ImageSize -> 360,
          PlotLabel -> BoxedLabel[aIt == N[a]]
        ] // PrettyString["zeta" -> "\[Zeta]"],
        (* Unphysical domain *)
        RegionPlot[RPolar[reZeta, imZeta] > 1,
          {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
          {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
          BoundaryStyle -> None,
          PlotStyle -> unphysStyle
        ],
        (* Non-viable domain *)
        RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
          {reZeta, -rhoMaxNon, rhoMaxNon},
          {imZeta, -rhoMaxNon, rhoMaxNon},
          BoundaryStyle -> termStyle,
          PlotStyle -> nonStyle
        ]
      ]
    , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}]
  , {{n, nInit, nIt}, nValues}]
]


(* ::Subsection:: *)
(*Interactive visualiser (\[Zeta]-space) (offset version)*)


DynamicModule[
 {nValues, nInit,
  gammaMin, gammaMax, gammaInit,
  aMin, aMax, aInit,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  (* Values of n *)
  nValues = Range[3, 5];
  nInit = First[nValues];
  (* Values of \[Gamma] *)
  gammaMin = 0.1;
  gammaMax = 4;
  gammaInit = 0.5;
  (* Values of A *)
  aMin = 0.01;
  aMax = 3;
  aInit = 0.17;
    (* non-viable lakes at vertices of polygon in z-space *)
  Manipulate[
    eps = 0.1;
    rhoMax = Exp[gamma];
    rhoMaxUnphys = rhoMax + eps;
    rhoMaxNon = SeekRoot[
      viOffset[gamma][n][a], (* viability \[CapitalPhi] *)
      {rhoMax, 1} (* seek largest \[Rho] *)
    ] + eps;
    Show[
      EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
        FrameLabel -> {Re["zeta"], Im["zeta"]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel @ Row[
          {gIt == N[gamma], aIt == N[a]},
          ","
        ]
      ] // PrettyString["zeta" -> "\[Zeta]"],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ]
    ]
  , {{n, nInit, nIt}, nValues}
  , {{gamma, gammaInit, gIt}, gammaMin, gammaMax, Appearance -> "Open"}
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}]
]


(* ::Subsection:: *)
(*Animations (\[Zeta]-space)*)


Table[
  Module[
   {ph,
    aN, aMin, aMax, aStep, aValues,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
   (* Azimuthal angle in \[Zeta]-space *)
   ph = 0;
   (* Values of A *)
   aN = aNat[n][ph];
   aMin = 1/10;
   aMax = Ceiling[aN];
   aStep = 1/10;
   aValues = Range[aMin, aMax, aStep];
   aValues = Append[aValues, aN] // Sort[#, Less] &;
   (* Animation *)
   Table[
     eps = 0.1;
     rhoMax = 1;
     rhoMaxUnphys = 1 + eps;
     rhoMaxNon = If[a < aN, rhoSharp[n][a][ph], rhoNat[n][ph]];
     Show[
       EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
         FrameLabel -> {Re["zeta"], Im["zeta"]},
         ImageSize -> 360,
         PlotLabel -> BoxedLabel[aIt == N[a]]
       ] // PrettyString["zeta" -> "\[Zeta]"],
       (* Unphysical domain *)
       RegionPlot[RPolar[reZeta, imZeta] > 1,
         {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
         {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
         BoundaryStyle -> None,
         PlotStyle -> unphysStyle
       ],
       (* Non-viable domain *)
       RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
         {reZeta, -rhoMaxNon, rhoMaxNon},
         {imZeta, -rhoMaxNon, rhoMaxNon},
         BoundaryStyle -> termStyle,
         PlotStyle -> nonStyle
       ]
     ]
   , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_zeta-viable-full-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 5}]


Table[
  Module[
   {ph,
    aN, aStep, aMin, aMax, aValues,
    eps, rhoN, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
   (* Azimuthal angle in \[Zeta]-space *)
   ph = 0;
   (* Values of A *)
   aN = aNat[n][ph];
   aStep = 1/20;
   aMin = Floor[aN - 5 aStep, aStep];
   aMax = Ceiling[aN, aStep];
   aValues = Range[aMin, aMax, aStep];
   aValues = Append[aValues, aN] // Sort[#, Less] &;
   (* Animation *)
   Table[
     eps = 0.1;
     rhoN = rhoNat[n][ph];
     rhoMax = 2 rhoN;
     rhoMaxUnphys = 1 + eps;
     rhoMaxNon = If[a < aN, rhoSharp[n][a][ph], rhoN];
     Show[
       EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
         FrameLabel -> {Re["zeta"], Im["zeta"]},
         ImageSize -> 360,
         PlotLabel -> BoxedLabel[aIt == N[a]]
       ] // PrettyString["zeta" -> "\[Zeta]"],
       (* Non-viable domain *)
       RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
         {reZeta, -rhoMaxNon, rhoMaxNon},
         {imZeta, -rhoMaxNon, rhoMaxNon},
         BoundaryStyle -> termStyle,
         PlotStyle -> nonStyle
       ]
     ]
   , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_zeta-viable-zoom-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Animation (\[Zeta]-space) (offset version)*)


Module[
 {n, gamma,
  aMin, aMax, aStep, aValues,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  n = nOffset;
  gamma = gammaOffset;
  (* Values of A *)
  aMin = 1/10;
  aMax = 3;
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  (* Animation *)
  eps = 0.1;
  rhoMax = Exp[gamma];
  rhoMaxUnphys = rhoMax + eps;
  Table[
    rhoMaxNon = SeekRoot[
      viOffset[gamma][n][a], (* viability \[CapitalPhi] *)
      {rhoMax, 1} (* seek largest \[Rho] *)
    ] + eps;
    Show[
      EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
        FrameLabel -> {Re["zeta"], Im["zeta"]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel[aIt == N[a]]
      ] // PrettyString["zeta" -> "\[Zeta]"],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ]
    ]
  , {a, aValues}]
] // Ex["polygon_offset_zeta-viable.gif", gifOpts]


(* ::Subsection:: *)
(*Small A plot (\[Zeta]-space)*)


(* ::Text:: *)
(*Plots to exaggerate \[CurlyPhi]-dependence.*)
(*We see that hyperbolic critical terminal points occur along \[CurlyPhi] = 2\[Pi] k/n,*)
(*whereas elliptic critical terminal points occur along \[CurlyPhi] = 2\[Pi] (k - 1/2)/n.*)


Table[
  Module[
   {a,
    eps, rhoN,
    rhoMax, rhoMaxUnphys,
    rhoMaxNon, rhoMinNon
   },
    (* Small value of A *)
    a = aNat[n][0] / 10^3;
    (* Plot *)
    eps = 0.1;
    rhoMax = 1;
    rhoMaxUnphys = 1 + eps;
    rhoMaxNon = rhoSharp[n][a][0];
    rhoMinNon = rhoSharp[n][a][Pi / n];
    Show[
      EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
        FrameLabel -> {Re["zeta"], Im["zeta"]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel[aIt == N[a]]
      ] // PrettyString["zeta" -> "\[Zeta]"],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > 1,
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Hyperbolic critical terminal points *)      
      (* (\[CurlyPhi] == 2 Pi k / n) *)
      Graphics @ {Directive[hypStyle, pointStyle],
        Point @ CirclePoints[{rhoMaxNon, 0}, n]
      },
      (* (corresponding contour) *)
      Graphics @ {hypStyle,
        Circle[{0, 0}, rhoMaxNon]
      },
      (* Elliptic critical terminal points *)
      (* (\[CurlyPhi] == 2 Pi (k - 1/2) / n) *)
      Graphics @ {Directive[ellStyle, pointStyle],
        Point @ CirclePoints[{rhoMinNon, Pi / n}, n]
      },
      (* (corresponding contour) *)
      Graphics @ {ellStyle,
        Circle[{0, 0}, rhoMinNon]
      }
    ]
  ] // Ex @ StringJoin["polygon_zeta-viable-small-a-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsection:: *)
(*Hot regime*)


(* ::Subsubsection:: *)
(*Starting points (\[Zeta]-space)*)


Table[
  Module[
   {a, idList,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
    a = aHot[n];
    (* Group names *)
    idList = {"general", "hyperbolic"};
    (* Plot ranges *)
    eps = 0.1;
    rhoMax = 1;
    rhoMaxUnphys = 1 + eps;
    rhoMaxNon = rhoSharp[n][a][0];
    (* Plot *)
    Show[
      EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
        FrameLabel -> {Re["zeta"], Im["zeta"]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel @ Row[
          {aIt == N[aHot[n]], "hot" // ""},
          "  "
        ]
      ] // PrettyString["zeta" -> "\[Zeta]"],
      (* Non-viable domain *)
      RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > 1,
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Starting points *)
      ListPlot[
        Table[startZetaHot[n][id] // ReIm, {id, idList}],
        LabelingFunction -> Function @ Placed[#2[[2]], Center],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
  ] // Ex @ StringJoin["polygon_zeta-hot-traced-starting-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Starting points (z-space)*)


Table[
  Module[{idList},
    (* Group names *)
    idList = {"general", "hyperbolic"};
    (* Plot *)
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n],
      (* Unphysical domain *)
      unphysDomain[n],
      (* Starting points *)
      ListPlot[
        Table[startZetaHot[n][id] // zMap[n] // ReIm, {id, idList}],
        LabelingFunction -> Function @ Placed[#2[[2]], Center],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
  ] // Ex @ StringJoin["polygon_z-hot-traced-starting-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*General (\[Zeta]-space)*)


Table[
  Module[
   {a,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
    a = aHot[n];
    (* Plot ranges *)
    eps = 0.1;
    rhoMax = 1;
    rhoMaxUnphys = 1 + eps;
    rhoMaxNon = rhoSharp[n][a][0];
    (* Plot *)
    Show[
      EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
        FrameLabel -> {Re["zeta"], Im["zeta"]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel @ Row[
          {aIt == N[aHot[n]], "hot" // ""},
          "  "
        ]
      ] // PrettyString["zeta" -> "\[Zeta]"],
      (* Non-viable domain *)
      RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Traced boundaries *)
      Table[
        Table[
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> {upperStyle, lowerStyle}
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"general"}}],
      Table[
        Table[
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> glowStyle
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"hyperbolic"}}],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > 1,
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ]
    ]
  ] // Ex @ StringJoin["polygon_zeta-hot-traced-full-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*General (z-space)*)


Table[
  Module[{},
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n],
      (* Traced boundaries *)
      Table[
        Table[
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // zMap[n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> {upperStyle, lowerStyle}
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"general"}}],
      Table[
        Table[
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // zMap[n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> glowStyle
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"hyperbolic"}}],
      (* Unphysical domain *)
      unphysDomain[n]
    ]
  ] // Ex @ StringJoin["polygon_z-hot-traced-full-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Hyperbolic-only animations (\[Zeta]-space)*)


Table[
  Module[
   {ph,
    aN, aMin, aMax, aStep, aValues,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon,
    zeta
   },
    (* Azimuthal angle in \[Zeta]-space *)
    ph = 0;
    (* Values of A *)
    aN = aNat[n][ph];
    aMin = 1/10;
    aStep = 1/10;
    aMax = Floor[aN, aStep];
    aValues = Range[aMin, aMax, aStep];
    aValues = Append[aValues, aN] // Sort[#, Less] &;
    (* Animation *)
    Table[
      eps = 0.1;
      rhoMax = 1;
      rhoMaxUnphys = 1 + eps;
      rhoMaxNon = If[a < aN, rhoSharp[n][a][ph], rhoNat[n][ph]];
      Show[
        EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
          FrameLabel -> {Re["zeta"], Im["zeta"]},
          ImageSize -> 360,
          PlotLabel -> BoxedLabel[aIt == N[a]]
        ] // PrettyString["zeta" -> "\[Zeta]"],
        (* Unphysical domain *)
        RegionPlot[RPolar[reZeta, imZeta] > 1,
          {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
          {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
          BoundaryStyle -> None,
          PlotStyle -> unphysStyle
        ],
        (* Non-viable domain *)
        RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
          {reZeta, -rhoMaxNon, rhoMaxNon},
          {imZeta, -rhoMaxNon, rhoMaxNon},
          BoundaryStyle -> termStyle,
          PlotStyle -> nonStyle
        ],
        (* Traced boundaries *)
        If[a < aN,
          (* If still in hot regime, plot *)
          zeta = zetaTraceHyperbolic[n][a];
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> hypStyle
            ]
          , {k, 0, n - 1}],
          (* Otherwise don't plot *)
          {}
        ]
      ]
    , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_zeta-traced-hyperbolic-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Hyperbolic-only animations (z-space)*)


Table[
  Module[
   {ph,
    aN, aMin, aMax, aStep, aValues,
    zeta
   },
    (* Azimuthal angle in \[Zeta]-space *)
    ph = 0;
    (* Values of A *)
    aN = aNat[n][ph];
    aMin = 1/10;
    aStep = 1/10;
    aMax = Floor[aN, aStep];
    aValues = Range[aMin, aMax, aStep];
    aValues = Append[aValues, aN] // Sort[#, Less] &;
    (* Animation *)
    Table[
      Show[
        (* Equipotentials and streamlines *)
        equipStream[n,
          ImageSize -> 360,
          PlotLabel -> BoxedLabel[aIt == N[a]]
        ],
        (* Traced boundaries *)
        If[a < aN,
          (* If still in hot regime, plot *)
          zeta = zetaTraceHyperbolic[n][a];
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // zMap[n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> hypStyle
            ]
          , {k, 0, n - 1}],
          (* Otherwise don't plot *)
          {}
        ],
        (* Unphysical domain *)
        unphysDomain[n]
      ]
    , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_z-traced-hyperbolic-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Convex domains*)


(* ::Subsubsection:: *)
(*Without known solution (z-space)*)


Table[
  Module[{a, zeta, rMax},
    a = aConvex[n];
    zeta = zetaTraceCand[n][a];
    rMax = zeta @ DomainStart[zeta] // zMap[n] // 1.2 Abs[#] &;
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n,
        PlotLabel -> BoxedLabel[aIt == N[a]],
        PlotRange -> {{-rMax, rMax}, {-rMax, rMax}}
      ],
      (* Traced boundaries *)
      Table[
        ParametricPlot[
          zeta[s] Exp[I 2 Pi k / n]
            // zMap[n]
            // {#, Conjugate[#]} &
            // ReIm
            // Evaluate,
          {s, DomainStart[zeta], DomainEnd[zeta]},
          PlotStyle -> convexStyle
        ]
      , {k, 0, n - 1}]
    ]
  ] // Ex @ StringJoin["polygon_z-convex-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Offset version*)


(* ::Subsubsection:: *)
(*Starting points (\[Zeta]-space) for A = 1.5 (joined)*)


Module[
 {n, gamma, a, rhoSh, idList,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetJoined;
  rhoSh = rhoOffsetJoinedSharp;
  (* Group names *)
  idList = {"terminal", "hyperbolic"};
  (* Plot ranges *)
  eps = 0.1;
  rhoMax = Exp[gamma];
  rhoMaxUnphys = rhoMax + eps;
  rhoMaxNon = rhoSh + eps;
  (* Plot *)
  Show[
    EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
      FrameLabel -> {Re["zeta"], Im["zeta"]},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ] // PrettyString["zeta" -> "\[Zeta]"],
    (* Unphysical domain *)
    RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
      {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
      {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* Non-viable domain *)
    RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
      {reZeta, -rhoMaxNon, rhoMaxNon},
      {imZeta, -rhoMaxNon, rhoMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[startZetaOffsetJoined[id] // ReIm, {id, idList}],
      LabelingFunction -> Function @ Placed[#2[[2]], Center],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["polygon_offset_zeta-joined-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*Starting points (\[Zeta]-space) for A = 1.7 (split)*)


Module[
 {n, gamma, a, rhoSh, rhoB, rhoA, idList,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetSplit;
  rhoSh = rhoOffsetSplitSharp;
  rhoB = rhoOffsetSplitB;
  rhoA = rhoOffsetSplitA;
  (* Group names *)
  idList = {
    "terminal-moat",
    "terminal-lake",
    "hyperbolic-moat",
    "hyperbolic-lake"
  };
  (* Plot ranges *)
  eps = 0.1;
  rhoMax = Exp[gamma];
  rhoMaxUnphys = rhoMax + eps;
  rhoMaxNon = rhoSh + eps;
  (* Plot *)
  Show[
    EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
      FrameLabel -> {Re["zeta"], Im["zeta"]},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ] // PrettyString["zeta" -> "\[Zeta]"],
    (* Unphysical domain *)
    RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
      {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
      {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* Non-viable domain *)
    RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
      {reZeta, -rhoMaxNon, rhoMaxNon},
      {imZeta, -rhoMaxNon, rhoMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[startZetaOffsetSplit[id] // ReIm, {id, idList}],
      LabelingFunction -> Function @ Placed[#2[[2]], Center],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["polygon_offset_zeta-split-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General (\[Zeta]-space) for A = 1.5 (joined)*)


Module[
 {n, gamma, a, rhoSh,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetJoined;
  rhoSh = rhoOffsetJoinedSharp;
  (* Plot ranges *)
  eps = 0.1;
  rhoMax = Exp[gamma];
  rhoMaxUnphys = rhoMax + eps;
  rhoMaxNon = rhoSh + eps;
  (* Plot *)
  Show[
    EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
      FrameLabel -> {Re["zeta"], Im["zeta"]},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ] // PrettyString["zeta" -> "\[Zeta]"],
    (* Unphysical domain *)
    RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
      {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
      {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* Non-viable domain *)
    RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
      {reZeta, -rhoMaxNon, rhoMaxNon},
      {imZeta, -rhoMaxNon, rhoMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetJoined[id]}]
    , {id, {"terminal"}}],
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> glowStyle
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetJoined[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["polygon_offset_zeta-joined-traced-full.pdf"]


(* ::Subsubsection:: *)
(*General (z-space) for A = 1.5 (joined)*)


Module[{n, a},
  n = nOffset;
  a = aOffsetJoined;
  (* Plot *)
  Show[
    equipStream[n,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ],
    (* \[Rho] == 1 *)
    Graphics @ {Directive[EdgeForm[contStyle], FaceForm[None]],
      poly[3]
    },
    (* Traced boundaries *)
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetJoined[id]}]
    , {id, {"terminal"}}],
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> glowStyle
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetJoined[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["polygon_offset_z-joined-traced-full.pdf"]


(* ::Subsubsection:: *)
(*General (\[Zeta]-space) for A = 1.7 (split)*)


Module[
 {n, gamma, a, rhoSh,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  n = nOffset;
  gamma = gammaOffset;
  a = aOffsetSplit;
  rhoSh = rhoOffsetSplitSharp;
  (* Plot ranges *)
  eps = 0.1;
  rhoMax = Exp[gamma];
  rhoMaxUnphys = rhoMax + eps;
  rhoMaxNon = rhoSh + eps;
  (* Plot *)
  Show[
    EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
      FrameLabel -> {Re["zeta"], Im["zeta"]},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ] // PrettyString["zeta" -> "\[Zeta]"],
    (* Unphysical domain *)
    RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
      {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
      {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* Non-viable domain *)
    RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
      {reZeta, -rhoMaxNon, rhoMaxNon},
      {imZeta, -rhoMaxNon, rhoMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetSplit[id]}]
    , {id, {"terminal-moat", "terminal-lake"}}],
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> glowStyle
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetSplit[id]}]
    , {id, {"hyperbolic-moat", "hyperbolic-lake"}}]
  ]
] // Ex["polygon_offset_zeta-split-traced-full.pdf"]


(* ::Subsubsection:: *)
(*General (z-space) for A = 1.7 (split)*)


Module[{n, a},
  n = nOffset;
  a = aOffsetSplit;
  (* Plot *)
  Show[
    equipStream[n,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ],
    (* \[Rho] == 1 *)
    Graphics @ {Directive[EdgeForm[contStyle], FaceForm[None]],
      poly[3]
    },
    (* Traced boundaries *)
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetSplit[id]}]
    , {id, {"terminal-moat", "terminal-lake"}}],
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> glowStyle
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetSplit[id]}]
    , {id, {"hyperbolic-moat", "hyperbolic-lake"}}]
  ]
] // Ex["polygon_offset_z-split-traced-full.pdf"]


(* ::Subsection:: *)
(*Convex offset version*)


(* ::Subsubsection:: *)
(*Starting points (\[Zeta]-space)*)


Module[
 {n, gamma, a, rhoSh, rhoB, rhoA, idList,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  a = aOffsetConvex;
  rhoSh = rhoOffsetConvexSharp;
  rhoB = rhoOffsetConvexB;
  rhoA = rhoOffsetConvexA;
  (* Group names *)
  idList = {
    "general",
    "terminal-moat",
    "terminal-lake",
    "hyperbolic-moat",
    "hyperbolic-lake"
  };
  (* Plot ranges *)
  eps = 0.1;
  rhoMax = rhoSh;
  rhoMaxUnphys = rhoMax + eps;
  rhoMaxNon = rhoSh + eps;
  (* Plot *)
  Show[
    EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
      FrameLabel -> {Re["zeta"], Im["zeta"]},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ] // PrettyString["zeta" -> "\[Zeta]"],
    (* Unphysical domain *)
    RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
      {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
      {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* Non-viable domain *)
    RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
      {reZeta, -rhoMaxNon, rhoMaxNon},
      {imZeta, -rhoMaxNon, rhoMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[startZetaOffsetConvex[id] // ReIm, {id, idList}],
      LabelingFunction -> Function @ Placed[#2[[2]], Center],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["polygon_offset_zeta-convex-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General (\[Zeta]-space)*)


Module[
 {n, gamma, a, rhoSh,
  eps, rhoMax, rhoMaxUnphys, rhoMaxNon
 },
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  a = aOffsetConvex;
  rhoSh = rhoOffsetConvexSharp;
  (* Plot ranges *)
  eps = 0.1;
  rhoMax = 1.2 rhoSh;
  rhoMaxUnphys = rhoMax + eps;
  rhoMaxNon = rhoSh + eps;
  (* Plot *)
  Show[
    EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
      FrameLabel -> {Re["zeta"], Im["zeta"]},
      ImageSize -> 720,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ] // PrettyString["zeta" -> "\[Zeta]"],
    (* Unphysical domain *)
    RegionPlot[RPolar[reZeta, imZeta] > Exp[gamma],
      {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
      {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* Non-viable domain *)
    RegionPlot[viOffset[gamma][n][a][reZeta + I imZeta] < 0,
      {reZeta, -rhoMaxNon, rhoMaxNon},
      {imZeta, -rhoMaxNon, rhoMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetConvex[id]}]
    , {id, {"general", "terminal-moat", "terminal-lake"}}],
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> glowStyle
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetConvex[id]}]
    , {id, {"hyperbolic-moat", "hyperbolic-lake"}}]
  ]
] // Ex["polygon_offset_zeta-convex-traced-full.pdf"]


(* ::Subsubsection:: *)
(*General (z-space)*)


Module[{n, a},
  n = nOffsetConvex;
  a = aOffsetConvex;
  (* Plot *)
  Show[
    equipStream[n,
      ImageSize -> 720,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ],
    (* \[Rho] == 1 *)
    Graphics @ {Directive[EdgeForm[contStyle], FaceForm[None]],
      poly[3]
    },
    (* Traced boundaries *)
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetConvex[id]}]
    , {id, {"general", "terminal-moat", "terminal-lake"}}],
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], DomainEnd[zeta]},
            PlotStyle -> glowStyle
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetConvex[id]}]
    , {id, {"hyperbolic-moat", "hyperbolic-lake"}}]
  ]
] // Ex["polygon_offset_z-convex-traced-full.pdf"]


(* ::Subsubsection:: *)
(*Hyperbolic from moat (z-space)*)


Module[{n, a},
  n = nOffsetConvex;
  a = aOffsetConvex;
  (* Plot *)
  Show[
    equipStream[n,
      PlotLabel -> BoxedLabel[aIt == N[a]]
    ],
    (* \[Rho] == 1 *)
    Graphics @ {Directive[EdgeForm[contStyle], FaceForm[None]],
      poly[3]
    },
    (* Traced boundaries *)
    Table[
      Table[
        Table[
          ParametricPlot[
            zeta[s] Exp[I 2 Pi k / n]
              // zMap[n]
              // {#, Conjugate[#]} &
              // ReIm
              // Evaluate,
            {s, DomainStart[zeta], 0},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {k, 0, n - 1}]
      , {zeta, zetaTraOffsetConvex[id]}]
    , {id, {"hyperbolic-moat"}}]
  ]
] // Ex["polygon_offset_z-convex-traced-hyperbolic_moat.pdf"]


(* ::Subsubsection:: *)
(*Rotund convex domain*)


Module[
 {n, a, rhoA,
  rA, rMax,
  zeta
 },
  n = nOffsetConvex;
  a = aOffsetConvex;
  rhoA = rhoOffsetConvexA;
  (* Plot range *)
  rA = rhoA // zMap[n];
  rMax = 2 rA;
  (* Plot *)
  Show[
    equipStream[n,
      PlotLabel -> BoxedLabel[aIt == N[a]],
      PlotRange -> {{-rMax, rMax}, {-rMax, rMax}}
    ],
    (* \[Rho] == 1 *)
    Graphics @ {Directive[EdgeForm[contStyle], FaceForm[None]],
      poly[3]
    },
    (* Traced boundaries *)
    zeta = zetaTraOffsetConvex["hyperbolic-moat"] // First;
    Table[
      ParametricPlot[
        zeta[s] Exp[I 2 Pi k / n]
          // zMap[n]
          // {#, Conjugate[#]} &
          // ReIm
          // Evaluate,
        {s, sOffsetConvex["rotund"], 0},
        PlotStyle -> convexStyle
      ]
    , {k, 0, n - 1}]
  ]
] // Ex["polygon_offset_z-convex-traced-rotund.pdf"]


(* ::Subsubsection:: *)
(*Elongated convex domain*)


Module[
 {n, a, rhoA,
  rA, rMax,
  zeta
 },
  n = nOffsetConvex;
  a = aOffsetConvex;
  rhoA = rhoOffsetConvexA;
  (* Plot range *)
  rA = rhoA // zMap[n];
  rMax = 3.5 rA;
  (* Plot *)
  Show[
    equipStream[n,
      PlotLabel -> BoxedLabel[aIt == N[a]],
      PlotRange -> {{-rMax, rMax}, {-rMax, rMax}}
    ],
    (* \[Rho] == 1 *)
    Graphics @ {Directive[EdgeForm[contStyle], FaceForm[None]],
      poly[3]
    },
    (* Traced boundaries *)
    zeta = zetaTraOffsetConvex["hyperbolic-moat"] // First;
    Table[
      ParametricPlot[
        zeta[s] Exp[I 2 Pi k / n]
          // zMap[n]
          // {#, Conjugate[#]} &
          // ReIm
          // Evaluate,
        {s, sOffsetConvex["elongated"], 0},
        PlotStyle -> convexStyle
      ]
    , {k, {1}}],
    Table[
      ParametricPlot[
        zeta[s] Exp[I 2 Pi k / n]
          // zMap[n]
          // {#, Conjugate[#]} &
          // ReIm
          // Evaluate,
        {s, sOffsetConvex["rotund"], 0},
        PlotStyle -> convexStyle
      ]
    , {k, {2}}]
  ]
] // Ex["polygon_offset_z-convex-traced-elongated.pdf"]


(* ::Section:: *)
(*Traced boundary curvature*)


(* ::Subsection:: *)
(*Candidate boundaries*)


(* ::Subsubsection:: *)
(*Meeting (A = A_m) (z-space)*)


Table[
  Module[{a, zeta},
    a = aMeet[n];
    zeta = zetaTraceCand[n][a];
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n,
        PlotLabel -> BoxedLabel[aIt == N[a]]
      ],
      (* Candidate traced boundaries *)
      Table[
        ParametricPlot[
          zeta[s] Exp[I 2 Pi k / n]
            // zMap[n]
            // {#, Conjugate[#]} &
            // ReIm
            // Evaluate,
          {s, DomainStart[zeta], DomainEnd[zeta]},
          PlotStyle -> hypStyle
        ]
      , {k, 0, n - 1}],
      (* Unphysical domain *)
      unphysDomain[n]
    ]
  ] // Ex @ StringJoin["polygon_z-candidate-meeting-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Curvature at \[CurlyPhi] = \[Pi]/n*)


Module[
 {nValues, dest, aCurTables,
  aMin, aMax, num, aValues,
  inflDotStyle = Directive[Red, Opacity[0.7], pointStyle]
 },
  nValues = Range[3, 5];
  (* Compute tables of (A, cur.) values *)
  (* (This is slightly slow (~4 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "polygon-candidate-curvature-corner.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    aCurTables = Import[dest] // Uncompress,
    (* Otherwise compute and store for next time *)
    aCurTables = Table[
      aMin = aMeet[n];
      aMax = Min[aNat[n][0], 3];
      num = 100;
      aValues = Subdivide[aMin, aMax, num] // Most;
      Table[
        {a, curCandCorner[n][a]}
      , {a, aValues}]
    , {n, nValues}];
    aCurTables // Compress // Ex[dest]
  ];
  (* Plot *)
  Show[
    (* A vs corner curvature *)
    ListPlot[
      aCurTables,
      AxesLabel -> {aIt, "Corner cur."},
      Joined -> True,
      PlotLegends -> LineLegend[nValues, LegendLabel -> nIt],
      PlotRange -> {{All, 3}, {All, 5}},
      PlotOptions[Axes] // Evaluate
    ],
    (* A_i (inflection) *)
    Graphics @ {inflDotStyle,
      Point @ Table[{aInfl[n], 0}, {n, nValues}]
    },
    Graphics @ Text[
      Table[{n, aInfl[n]}, {n, nValues}]
        // TableForm[#, TableHeadings -> {{}, {nIt, Subscript[aIt, "i"]}}] &
        // Style[#, 17] &,
      {aMeet[3] + 0.15, 3}
    ]
  ] // Ex["polygon-candidate-curvature-corner.pdf"]
]


(* ::Subsubsection:: *)
(*Inflection (A = A_i) (z-space)*)


Table[
  Module[{a, zeta, rMax},
    a = aInfl[n];
    zeta = zetaTraceCand[n][a];
    rMax = zeta @ DomainStart[zeta] // zMap[n] // 1.2 Abs[#] &;
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n,
        PlotLabel -> BoxedLabel[aIt == N[a]],
        PlotRange -> {{-rMax, rMax}, {-rMax, rMax}}
      ],
      (* Candidate traced boundaries *)
      Table[
        ParametricPlot[
          zeta[s] Exp[I 2 Pi k / n]
            // zMap[n]
            // {#, Conjugate[#]} &
            // ReIm
            // Evaluate,
          {s, DomainStart[zeta], DomainEnd[zeta]},
          PlotStyle -> hypStyle
        ]
      , {k, 0, n - 1}]
    ]
  ] // Ex @ StringJoin["polygon_z-candidate-inflection-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Table of critical A values*)


(* ::Text:: *)
(*This is the table of values on Page r4-12.*)


Table[
  {n, aMeet[n], aInfl[n], aNat[n][0]}
, {n, 3, 5}] // TableForm[#,
  TableHeadings -> {None, {"n", "A_m", "A_i", "A_\[Natural]"}}
] & // Ex["polygon-a-critical.pdf"]


(* ::Subsection:: *)
(*Check convexity of rotund and elongated domains*)


Module[{n, gamma, a, zeta},
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  a = aOffsetConvex;
  zeta = zetaTraOffsetConvex["hyperbolic-moat"] // First;
  Plot[
    curTraOffset[gamma][n][a] @ zeta[s],
    {s, sOffsetConvex["elongated"], 0},
    AxesLabel -> {sIt, "Cur."},
    PlotRange -> {0, Automatic},
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["polygon_offset-convex-curvature.pdf"]


(* ::Section:: *)
(*Numerical verification plots*)


(* ::Subsection:: *)
(*Finite element mesh*)


Table[
  Module[
   {source,
    a, tBath, mesh, prExt, prInt,
    bCoords, bCoordsExt, bCoordsInt
   },
    (* Import mesh *)
    source = "polygon-convex-verification-mesh-" <> ToString[n] <> ".txt";
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
  ] // Ex @ StringJoin["polygon-convex-verification-mesh-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Numerical solution*)


Table[
  Module[{source, tSol, mesh},
    (* Import solution *)
    source = "polygon-convex-verification-solution-" <> ToString[n] <> ".txt";
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
  ] // Ex @ StringJoin["polygon-convex-verification-solution-", ToString[n], ".png"]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Relative errors (3D)*)


(* (These are slow (~12 sec).) *)
Table[
  Module[{source, tSol, mesh, tExact},
    (* Import solution *)
    source = "polygon-convex-verification-solution-" <> ToString[n] <> ".txt";
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Known exact solution *)
    tExact[x_, y_] := Re @ g @ zetaMap[n][x + I y];
    (* Plot *)
    With[{x = \[FormalX], y = \[FormalY]},
      Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
        Exclusions -> None,
        PlotLabel -> "Relative error of numerical solution",
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ]
    ]
  ] // Ex @ StringJoin[
    "polygon-convex-verification-rel_error-3d-", ToString[n], ".png"
  ]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Relative errors (2D)*)


Table[
  Module[{source, tSol, mesh, tExact},
    (* Import solution *)
    source = "polygon-convex-verification-solution-" <> ToString[n] <> ".txt";
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Known exact solution *)
    tExact[x_, y_] := Re @ g @ zetaMap[n][x + I y];
    (* Plot *)
    With[{x = \[FormalX], y = \[FormalY]},
      Show[
        DensityPlot[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
          ColorFunction -> "Rainbow",
          Exclusions -> None,
          PlotLabel -> "Relative error of numerical solution",
          PlotRange -> Full,
          PlotLegends -> Automatic,
          PlotOptions[Frame] // Evaluate
        ],
        mesh["Wireframe"]
      ]
    ]
  ] // Ex @ StringJoin[
    "polygon-convex-verification-rel_error-2d-", ToString[n], ".png"
  ]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Finite element mesh (convex offset version)*)


(* ::Subsubsection:: *)
(*Rotund convex domain*)


Module[
 {source,
  a, tBath, mesh, prExt, prInt,
  bCoords, bCoordsExt, bCoordsInt
 },
  (* Import mesh *)
  source = "polygon_offset-convex-verification-mesh-rotund.txt";
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
] // Ex["polygon_offset-convex-verification-mesh-rotund.pdf"]


(* ::Subsubsection:: *)
(*Elongated convex domain*)


Module[
 {source,
  a, tBath, mesh, prExt, prInt,
  bCoords, bCoordsExt, bCoordsInt
 },
  (* Import mesh *)
  source = "polygon_offset-convex-verification-mesh-elongated.txt";
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
] // Ex["polygon_offset-convex-verification-mesh-elongated.pdf"]


(* ::Subsection:: *)
(*Numerical solution*)


(* ::Subsubsection:: *)
(*Rotund convex domain*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "polygon_offset-convex-verification-solution-rotund.txt";
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
] // Ex["polygon_offset-convex-verification-solution-rotund.png"]


(* ::Subsubsection:: *)
(*Elongated convex domain*)


Module[{source, tSol, mesh},
  (* Import solution *)
  source = "polygon_offset-convex-verification-solution-elongated.txt";
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
] // Ex["polygon_offset-convex-verification-solution-elongated.png"]


(* ::Subsection:: *)
(*Relative errors (3D)*)


(* ::Subsubsection:: *)
(*Rotund convex domain*)


Module[{n, gamma, source, tSol, mesh, tExact},
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  (* Import solution *)
  source = "polygon_offset-convex-verification-solution-rotund.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Re @ gOffset[gamma] @ zetaMap[n][x + I y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      Exclusions -> None,
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["polygon_offset-convex-verification-rel_error-3d-rotund.png"]


(* ::Subsubsection:: *)
(*Elongated convex domain*)


Module[{n, gamma, source, tSol, mesh, tExact},
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  (* Import solution *)
  source = "polygon_offset-convex-verification-solution-elongated.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Re @ gOffset[gamma] @ zetaMap[n][x + I y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
      Exclusions -> None,
      PlotLabel -> "Relative error of numerical solution",
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["polygon_offset-convex-verification-rel_error-3d-elongated.png"]


(* ::Subsection:: *)
(*Relative errors (2D)*)


(* ::Subsubsection:: *)
(*Rotund convex domain*)


Module[{n, gamma, source, tSol, mesh, tExact},
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  (* Import solution *)
  source = "polygon_offset-convex-verification-solution-rotund.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Re @ gOffset[gamma] @ zetaMap[n][x + I y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Show[
      DensityPlot[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
        ColorFunction -> "Rainbow",
        Exclusions -> None,
        PlotLabel -> "Relative error of numerical solution",
        PlotLegends -> Automatic,
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ],
      mesh["Wireframe"]
    ]
  ]
] // Ex["polygon_offset-convex-verification-rel_error-2d-rotund.png"]


(* ::Subsubsection:: *)
(*Elongated convex domain*)


Module[{n, gamma, source, tSol, mesh, tExact},
  n = nOffsetConvex;
  gamma = gammaOffsetConvex;
  (* Import solution *)
  source = "polygon_offset-convex-verification-solution-elongated.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Known exact solution *)
  tExact[x_, y_] := Re @ gOffset[gamma] @ zetaMap[n][x + I y];
  (* Plot *)
  With[{x = \[FormalX], y = \[FormalY]},
    Show[
      DensityPlot[tSol[x, y] / tExact[x, y] - 1, Element[{x, y}, mesh],
        AspectRatio -> Automatic,
        ColorFunction -> "Rainbow",
        Exclusions -> None,
        PlotLabel -> "Relative error of numerical solution",
        PlotLegends -> Automatic,
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ],
      mesh["Wireframe"]
    ]
  ]
] // Ex["polygon_offset-convex-verification-rel_error-2d-elongated.png"]
