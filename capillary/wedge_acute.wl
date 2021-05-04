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
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "wedge_acute"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Wedge half-angles*)


(* apd stands for "alpha per degree". *)
apdValues = {10, 15, 30, 40, 45, 50, 60, 75};


(* ::Subsection:: *)
(*Contact angles*)


(* gpd stands for "gamma per degree". *)
gpdValues = Join[{1}, Range[5, 85, 5]];


(* ::Subsection:: *)
(*Finite element mesh*)


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt",
    Module[
     {rMax, alpha,
      lenCoarse,
      rSpacingWet, rValuesWet,
      phiSpacingFar, phiValuesFar,
      bPointList, nB, mod,
      bMesh, mesh,
      prWet
     },
      (* Geometry *)
      rMax = 10;
      alpha = apd * Degree;
      (* Mesh coarse length scale *)
      lenCoarse = 0.2;
      (* Wet boundary (|\[Phi]| == \[Alpha]) points *)
      (* (fine length scale used here) *)
      rSpacingWet = 0.01;
      rValuesWet = UniformRange[0, rMax, rSpacingWet];
      (* Far boundary (R == R_max) points *)
      phiSpacingFar = lenCoarse / rMax;
      phiValuesFar = UniformRange[-alpha, alpha, phiSpacingFar];
      (* Boundary points *)
      bPointList = Join[
        (* Lower wall: \[Phi] == -\[Alpha], 0 < r < R_max *)
        Table[XYPolar[r, -alpha], {r, rValuesWet}] // Most,
        (* Far arc: r == R_max, -\[Alpha] < \[Phi] < \[Alpha] *)
        Table[XYPolar[rMax, phi], {phi, phiValuesFar}] // Most,
        (* Upper wall: \[Phi] == \[Alpha], 0 < r < R_max *)
        Table[XYPolar[r, alpha], {r, rValuesWet // Reverse}] // Most
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
        "ImproveBoundaryPosition" -> True,
        MeshRefinementFunction -> MeshRefinementUniform[lenCoarse]
      ];
      (* Predicate function for wetting boundaries (|\[Phi]| == \[Alpha]) *)
      prWet = Function[{x, y}, x <= rMax Cos[alpha] // Evaluate];
      {mesh, prWet}
        // Compress
    ]
  ]
, {apd, apdValues}]


(* ::Subsection:: *)
(*Solve boundary value problem*)


(* (These are slow (~1 hr), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
 *)
(* For each value of alpha: *)
Table[
  Module[{mesh, prWet, gamma, eval},
    (* If all values of gamma have been solved for *)
    If[
      Table[
        FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"
      , {gpd, gpdValues}]
        // AllTrue[FileExistsQ],
      (* Do nothing *)
      Null,
      (* Otherwise solve for solutions: *)
      (* Import mesh *)
      {mesh, prWet} =
        Import[
          FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt"
        ] // Uncompress;
      (* For each value of gamma: *)
      Table[
        gamma = gpd * Degree;
        ExportIfNotExists[
          FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt",
          eval = EvaluationData @ SolveLaplaceYoung[gamma, mesh, prWet];
          eval /@ {"Result", "Success", "FailureType", "MessagesText", "AbsoluteTiming"}
            // Compress
        ]
      , {gpd, gpdValues}]
    ]
  ]
, {apd, apdValues}]


(* ::Subsection::Closed:: *)
(*System of ODEs for T-contour (old)*)


(* ::Text:: *)
(*See (57) and (58) (Page 13 of manuscripts/boundary-tracing.pdf)*)
(*for arc length parametrisation of a contour.*)


tContourSystem[tSol_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{t, p, q, slope},
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Magnitude of gradient vector *)
      slope = Sqrt[p ^ 2 + q ^ 2];
      (* Return system of ODEs *)
      {
        x' == q / slope,
        y' == -p / slope
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsection::Closed:: *)
(*System of ODEs for terminal curve \[CapitalPhi] = 0 (old)*)


(* ::Text:: *)
(*Note that \[CapitalPhi] = 0 is equivalent to (\[Del]T)^2 = cot(\[Gamma])^2,*)
(*so a \[CapitalPhi]-contour is equivalent to a (\[Del]T)^2-contour.*)


(* vi means \[CapitalPhi]. *)
viContourSystem[tSol_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[
     {t, p, q,
      grad2, grad2P, grad2Q,
      grad2Slope
     },
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient vector *)
      grad2 = p^2 + q^2;
      (* Components of gradient of gradient squared *)
      grad2P = D[grad2, x];
      grad2Q = D[grad2, y];
      (* Magnitude of gradient of gradient squared *)
      grad2Slope = Sqrt[grad2P ^ 2 + grad2Q ^ 2];
      (* Return system of ODEs *)
      {
        x' == grad2Q / grad2Slope,
        y' == -grad2P / grad2Slope
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsection::Closed:: *)
(*System of ODEs for tracing (old)*)


(* ::Text:: *)
(*See (55) and (56) (Page 13 of manuscripts/boundary-tracing.pdf)*)
(*for arc length parametrisation of a traced boundary.*)


xyTraSystem[tSol_, gammaTra_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{t, p, q, grad2, f, vi},
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient vector *)
      grad2 = p^2 + q^2;
      (* Flux function F *)
      f = Cos[gammaTra] Sqrt[1 + grad2];
      (* Viability function \[CapitalPhi] *)
      vi = Sin[gammaTra]^2 * grad2 - Cos[gammaTra]^2;
      (* Return system of ODEs *)
      {
        x' == (-q f + p * Re @ Sqrt[vi]) / grad2,
        y' == (+p f + q * Re @ Sqrt[vi]) / grad2
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Text:: *)
(*Termination at terminal curve*)


xyTraSystemTerm[tSol_, gammaTra_] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{t, p, q, grad2, vi},
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient vector *)
      grad2 = p^2 + q^2;
      (* Viability function \[CapitalPhi] *)
      vi = Sin[gammaTra]^2 * grad2 - Cos[gammaTra]^2;
      (* Return WhenEvent for termination *)
      WhenEvent[vi < 0 // Evaluate, "StopIntegration"]
        /. {x -> x[s], y -> y[s]}
    ]
  ];


(* ::Subsection:: *)
(*Hyperbolic critical terminal point*)


x0[tSol_, gammaTra_] :=
  Module[{p, xMax},
    (* \[PartialD]T/\[PartialD]x *)
    p = Derivative[1, 0][tSol];
    (* Find x == x_0 for which -\[PartialD]T/\[PartialD]x == cot(gamma-tilde) *)
    xMax = 10;
    SeekRoot[p[#, 0] + Cot[gammaTra] &, {0, xMax}]
  ];


(* ::Subsection::Closed:: *)
(*Tracing (old)*)


(* ::Subsubsection:: *)
(*Representative angles*)


(* alpha, gamma and gamma-tilde (tracing gamma) per degree *)
apdRep = 40;
gpdRep = 60;
gpdTraRep = 70;


alphaRep = apdRep * Degree;
gammaRep = gpdRep * Degree;
gammaTraRep = gpdTraRep * Degree;


(* ::Subsubsection:: *)
(*Import known solution*)


tSolRep = Import[
  FString @ "solution/wedge_acute-solution-apd-{apdRep}-gpd-{gpdRep}.txt"
] // Uncompress // First;


meshRep = tSolRep["ElementMesh"];


(* ::Subsubsection:: *)
(*Same contact angle for tracing*)


(* ::Text:: *)
(*Hyperbolic critical terminal point (x_0)*)


x0Same = x0[tSolRep, gammaRep];


(* ::Text:: *)
(*Traced boundaries through (x, y) = (x_0, 0)*)


(* (This is somewhat slow (~2 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-traced-hyperbolic.txt",
  List @ With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{tSol, gamma, xInit, sMax},
      tSol = tSolRep;
      gamma = gammaRep;
      xInit = x0Same;
      (* Numerical integration *)
      sMax = 3 x0Same;
      NDSolveValue[
        {
          xyTraSystem[tSol, gamma],
          x[0] == xInit,
          y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ] // ReInterpolate
    ]
  ] // Compress
]


(* ::Text:: *)
(*Local T-contour through (x, y) = (x_0, 0)*)


(* (This is somewhat slow (~2 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-contour-hyperbolic.txt",
  List @ With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{tSol, xInit, sMax},
      tSol = tSolRep;
      xInit = x0Same;
      (* Numerical integration *)
      sMax = 3 x0Same;
      NDSolveValue[
        {
          tContourSystem[tSol],
          x[0] == xInit,
          y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ] // ReInterpolate
    ]
  ] // Compress
]


(* ::Text:: *)
(*Terminal curve*)


(* (This is somewhat slow (~2 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-terminal.txt",
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{tSol, xInit, sMax},
      tSol = tSolRep;
      xInit = x0Same;
      (* Numerical integration *)
      sMax = 3 x0Same;
      NDSolveValue[
        {
          viContourSystem[tSol],
          x[0] == xInit,
          y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ] // ReInterpolate
    ]
  ] // Compress
]


(* ::Text:: *)
(*Generic traced boundaries*)


(* (This is slow (~20 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-traced-general.txt",
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[
     {tSol, alpha, gamma, sMax,
      xMax, yMax,
      num, xyInitList, xInit, yInit,
      xyTraWall, xyTraTerm
     },
      tSol = tSolRep;
      alpha = alphaRep;
      gamma = gammaRep;
      (* Numerical integration *)
      sMax = 3 x0Same;
      (* Ranges *)
      xMax = 2 x0Same;
      yMax = xMax Tan[alpha];
      (* Traced boundaries from the lower wall *)
      num = 16;
      xyInitList = Subdivide[{0, 0}, {xMax, -yMax}, num] // Rest;
      xyTraWall =
        Table[
          {xInit, yInit} = xyInit;
          NDSolveValue[
            {
              xyTraSystem[tSol, gamma, -1],
              x[0] == xInit,
              y[0] == yInit,
              xyTraSystemTerm[tSol, gamma]
            }, {x, y}, {s, 0, sMax},
            NoExtrapolation
          ] // ReInterpolate
        , {xyInit, xyInitList}];
      (* Traced boundaries from terminal curve *)
      xyInitList =
        Cases[
          (* Ending points of xyTraWall *)
          Table[
            xy @ DomainEnd[xy] // Through
          , {xy, xyTraWall}],
          (* Reflect those which are on the terminal curve *)
          (* (due to traced boundary hitting terminal curve and terminating) *)
          {x_, y_} /; y < 0
            :>
          {x, -y}
        ];
      xyTraTerm =
        Table[
          {xInit, yInit} = xyInit;
          NDSolveValue[
            {
              xyTraSystem[tSol, gamma, -1],
              x[0] == xInit,
              y[0] == yInit,
              xyTraSystemTerm[tSol, gamma]
            }, {x, y}, {s, 0, sMax},
            NoExtrapolation
          ] // ReInterpolate
        , {xyInit, xyInitList}];
      (* Join lists *)
      Join[xyTraWall, xyTraTerm]
    ]
  ] // Compress
]


(* ::Subsection::Closed:: *)
(*Both branches of a curve (old)*)


(* Adds reflection about x == 0. *)
bothBranches[{x_, y_}] := {{x, y}, {x, -y}};


(* ::Subsection:: *)
(*Rounded wedges (circular arc)*)


(* ::Subsubsection:: *)
(*Finite element mesh*)


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString @ "mesh/wedge_acute-mesh-apd-{apd}-rad-{rad}.txt",
    Module[
      {
        rMax, alpha,
        xCentre, yCentre,
        rJoin,
        angMin, angMax,
        fineLengthScale, coarseLengthScale,
        boundaryPoints, numBoundaryPoints, boundaryElements,
        boundaryMesh, mesh,
        predicateWet,
        tipTheory, tipMesh,
        dummyForTrailingCommas
      },
      (* Geometry *)
      rMax = 10;
      alpha = apd * Degree;
      (* Centre of circular arc *)
      xCentre = rad / Sin[alpha];
      yCentre = 0;
      (* Radius (distance from origin) of arc--wall join *)
      (* (where the arc is tangential to the wall) *)
      rJoin = rad / Tan[alpha];
      (* Angular range of circular arc *)
      angMin = Pi/2 + alpha;
      angMax = 2 Pi - angMin;
      (* Refinement length scales along boundaries *)
      (* (refinement not needed in interior, since not boundary tracing afterward) *)
      fineLengthScale = 0.01;
      coarseLengthScale = 0.5;
      (* Boundary points and elements *)
      boundaryPoints =
        Join[
          (* Circular rounding arc *)
          (* (split into halves to ensure a point along y == 0) *)
          Table[
            {xCentre, yCentre} + XYPolar[rad, ang]
            , {ang, UniformRange[angMin, Pi, fineLengthScale / rad]}
          ] // Rest,
          Table[
            {xCentre, yCentre} + XYPolar[rad, ang]
            , {ang, UniformRange[Pi, angMax, fineLengthScale / rad]}
          ] // Rest,
          (* Lower straight wall *)
          Table[
            XYPolar[r, -alpha]
            , {r, UniformRange[rJoin, rMax, fineLengthScale]}
          ] // Rest,
          (* Far arc (numerical infinity) *)
          Table[
            XYPolar[rMax, phi]
            , {phi, UniformRange[-alpha, alpha, coarseLengthScale / rMax]}
          ] // Rest,
          (* Upper straight wall *)
          Table[
            XYPolar[r, +alpha]
            , {r, UniformRange[rMax, rJoin, -fineLengthScale]}
          ] // Rest,
          {}
        ];
      numBoundaryPoints = Length[boundaryPoints];
      boundaryElements =
        LineElement /@ {
          Table[{n, n + 1}, {n, numBoundaryPoints}]
            // Mod[#, numBoundaryPoints, 1] &
        };
      (* Build mesh *)
      boundaryMesh =
        ToBoundaryMesh[
          "Coordinates" -> boundaryPoints,
          "BoundaryElements" -> boundaryElements,
          {}
        ];
      mesh = ToElementMesh[boundaryMesh, "ImproveBoundaryPosition" -> True];
      (* Predicate function for wetting boundaries *)
      predicateWet = Function[{x, y}, x <= rMax Cos[alpha] // Evaluate];
      (* Rounding tip location *)
      tipTheory = {xCentre, yCentre} - {rad, 0};
      tipMesh = MinimalBy[mesh["Coordinates"], First] // First;
      (* Return *)
      {mesh, predicateWet, tipTheory, tipMesh} // Compress
    ]
  ]
  , {apd, {60}}
  , {rad, Range[0, 3, 0.2] // Rest}
]


(* ::Subsubsection:: *)
(*Solve for numerical rounding heights (conventional approach)*)


(* (This is slow (~30 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["wedge_acute-radius-height-rise-conventional.csv",
  Module[
    {
      apd, gpd, radValues, gamma,
      mesh, predicateWet, tipTheory, tipMesh,
      tNumerical, height,
      table, tSharp, heightSharp,
      dummyForTrailingCommas
    },
    (* Chosen parameters *)
    apd = 60;
    gpd = 75;
    radValues = Range[0, 3, 0.2] // Rest;
    gamma = gpd * Degree;
    (* For each chosen rounding radius: *)
    Table[
      (* Import mesh *)
      {mesh, predicateWet, tipTheory, tipMesh} =
        FString @ "mesh/wedge_acute-mesh-apd-{apd}-rad-{rad}.txt"
          // Import // Uncompress;
      (* Compute numerical solution to capillary BVP *)
      tNumerical = SolveLaplaceYoung[gamma, mesh, predicateWet];
      (* Obtain height at rounding tip *)
      height[rad] = tNumerical @@ tipMesh
      , {rad, radValues}
    ];
    (* Build radius vs height table *)
    table = Table[{rad, height[rad]}, {rad, radValues}];
    (* Prepend zero-radius data point (using sharp-cornered solution) *)
    tSharp =
      Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    heightSharp = tSharp[0, 0];
    table = Prepend[table, {0, heightSharp}];
    (* Return radius vs height table for export *)
    table
  ]
]


(* ::Subsubsection:: *)
(*Well-approximating radius vs rounding height (boundary tracing approach)*)


(* (This is slow (~6 sec).) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["wedge_acute-radius-height-rise-tracing.csv",
  Module[
    {
      apd, gpdTracing,
      alpha, gammaTracing,
      gpdValues,
      tNumerical,
      gamma, xCritical,
      d, xCriticalOffset,
      rad, height,
      dummyForTrailingCommas
    },
    (* Prescribed angles *)
    apd = 60;
    gpdTracing = 75;
    alpha = apd * Degree;
    gammaTracing = gpdTracing * Degree;
    (* For various known solution angles: *)
    gpdValues = Range[90 - apd, gpdTracing, 5];
    Table[
      gamma = gpd * Degree;
      (* Import numerical solutions *)
      tNumerical[gpd] =
        Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Uncompress // First;
      (* Compute critical terminal point x_0 *)
      xCritical[gpd] =
        Quiet[
          x0[tNumerical[gpd], gammaTracing]
          , {InterpolatingFunction::femdmval}
        ];
      (* Offset critical terminal point x'_0  *)
      d[gpd] = DHalfPlane[gamma, gammaTracing];
      xCriticalOffset[gpd] = xCritical[gpd] - d[gpd] / Sin[alpha];
      (* Well-approximating rounding radius *)
      rad[gpd] = xCriticalOffset[gpd] * Sin[alpha] / (1 - Sin[alpha]);
      (* Height rise *)
      height[gpd] = tNumerical[gpd][xCritical[gpd], 0];
      , {gpd, gpdValues}
    ];
    (* Return radius vs height table for export *)
    Table[{rad[gpd], height[gpd]}, {gpd, gpdValues}]
  ]
]


(* ::Subsection:: *)
(*Modifications (viable contours)*)


(* ::Subsubsection:: *)
(*Angular parameters*)


{apdMod, gpdMod} = {45, 60};
numContoursViable = 4;


(* ::Subsubsection:: *)
(*Contours*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["modification/wedge_acute-modification-contours.txt",
  Module[
    {
      gamma, tNumerical,
      xCritical,
      hWall, xCentreWall,
      numberUpToCritical, xStep, numberUpToWall,
      xCentreValues, hValues,
      sMax, xyContourList,
      dummyForTrailingCommas
    },
    (* Numerical solution *)
    gamma = gpdMod * Degree;
    tNumerical =
      Import @ FString["solution/wedge_acute-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
        // Uncompress // First;
    (* Critical terminal point x_0 *)
    xCritical = x0[tNumerical, gamma];
    (* Point along centreline which has wall height *)
    hWall = HHalfPlane[gamma];
    xCentreWall =
      SeekRoot[
        tNumerical[#, 0] - hWall &,
        {0, 10}
        , 5
      ];
    (* Centreline x-values and height rises for contours *)
    numberUpToCritical = numContoursViable;
    xStep = xCritical / numberUpToCritical;
    numberUpToWall = Floor[xCentreWall / xStep];
    xCentreValues = xCritical * Range[numberUpToWall] / numberUpToCritical;
    hValues = Table[tNumerical[x, 0], {x, xCentreValues}];
    (* Compute T-contours *)
    sMax = 4; (* probably enough *)
    xyContourList =
      Table[
        Quiet[
          ContourByArcLength[tNumerical][{x, 0}, 0, sMax {-1, 1}]
          , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
        ]
        , {x, xCentreValues}
      ];
    (* Return lists to export *)
    {xCentreValues, hValues, xyContourList} // Compress
  ]
]


(* ::Subsubsection:: *)
(*Finite element meshes*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["modification/wedge_acute-modification-meshes.txt",
  Module[
    {
      rMax, alpha,
      xCentreValues, hValues, xyContourList,
      fineLengthScale,
      coarseLengthScale,
        sContourUpper, sContourLower, sContourValues, xyContourPoints,
        rUpper1, phiUpper1, rUpper2, phiUpper2, rUpperJoin,
        rLower1, phiLower1, rLower2, phiLower2, rLowerJoin,
        boundaryPoints, numBoundaryPoints, boundaryElements,
        boundaryMesh, mesh,
        predicateWet,
      dummyForTrailingCommas
    },
    (* Geometry *)
    rMax = 10;
    alpha = apdMod * Degree;
    (* Import contours *)
    {xCentreValues, hValues, xyContourList} =
      Import["modification/wedge_acute-modification-contours.txt"] // Uncompress;
    (* Mesh length scales *)
    fineLengthScale = 0.01;
    coarseLengthScale = 0.5;
    (* For each contour *)
    Table[
      (* Raw T-contour points *)
      sContourUpper = DomainEnd[xy];
      sContourLower = DomainStart[xy];
      sContourValues = Join[
        UniformRange[sContourUpper, 0, -fineLengthScale],
        UniformRange[0, sContourLower, -fineLengthScale] // Rest,
        {}
      ];
      xyContourPoints = Table[EvaluatePair[xy, s], {s, sContourValues}];
      (* Interpolate linearly in polar coordinates unto upper wall *)
      {rUpper1, phiUpper1} = RPhiPolar @@ xyContourPoints[[1]];
      {rUpper2, phiUpper2} = RPhiPolar @@ xyContourPoints[[2]];
      rUpperJoin = Rescale[alpha, {phiUpper2, phiUpper1}, {rUpper2, rUpper1}];
      (* Interpolate linearly in polar coordinates unto lower wall *)
      {rLower1, phiLower1} = RPhiPolar @@ xyContourPoints[[-1]];
      {rLower2, phiLower2} = RPhiPolar @@ xyContourPoints[[-2]];
      rLowerJoin = Rescale[-alpha, {phiLower2, phiLower1}, {rLower2, rLower1}];
      (* Boundary points and elements *)
      boundaryPoints =
        Join[
          (* T-contour portion *)
          xyContourPoints // Rest // Most,
          (* Lower straight wall *)
          Table[
            XYPolar[r, -alpha]
            , {r, UniformRange[rLowerJoin, rMax, fineLengthScale]}
          ] // Most,
          (* Far arc (numerical infinity) *)
          Table[
            XYPolar[rMax, phi]
            , {phi, UniformRange[-alpha, alpha, coarseLengthScale / rMax]}
          ] // Most,
          (* Upper straight wall *)
          Table[
            XYPolar[r, +alpha]
            , {r, UniformRange[rMax, rUpperJoin, -fineLengthScale]}
          ],
          {}
        ];
      numBoundaryPoints = Length[boundaryPoints];
      boundaryElements =
        LineElement /@ {
          Table[{n, n + 1}, {n, numBoundaryPoints}]
            // Mod[#, numBoundaryPoints, 1] &
        };
      (* Build mesh *)
      boundaryMesh =
        ToBoundaryMesh[
          "Coordinates" -> boundaryPoints,
          "BoundaryElements" -> boundaryElements,
          {}
        ];
      mesh = ToElementMesh[boundaryMesh, "ImproveBoundaryPosition" -> True];
      (* Predicate function for wetting boundaries *)
      predicateWet = Function[{x, y}, x <= rMax Cos[alpha] // Evaluate];
      (* Return mesh *)
      {mesh, predicateWet, rUpperJoin, rLowerJoin}
      , {xy, xyContourList}
    ] // Compress
  ]
]


(* ::Subsubsection:: *)
(*Solve capillary BVP and extract height rise profiles*)


(* (This is slow (~15 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["modification/wedge_acute-modification-profiles.txt",
  Module[
    {
      alpha, gamma,
      meshStuffList,
      mesh, predicateWet, rUpperJoin, rLowerJoin,
      meshCoordinates, boundaryElements, boundaryCoordinates,
      tNumerical, heightValues,
      dummyForTrailingCommas
    },
    (* Chosen parameters *)
    {alpha, gamma} = {apdMod, gpdMod} Degree;
    (* Import meshes *)
    meshStuffList = Import["modification/wedge_acute-modification-meshes.txt"] // Uncompress;
    (* For each mesh: *)
    Table[
      (* Mesh stuff *)
      {mesh, predicateWet, rUpperJoin, rLowerJoin} = meshStuff;
      (* Extract T-contour boundary coordinates, sorted by y *)
      meshCoordinates = mesh["Coordinates"];
      boundaryElements = List @@ First @ mesh["BoundaryElements"] // Flatten;
      boundaryCoordinates = meshCoordinates[[boundaryElements]];
      boundaryCoordinates = boundaryCoordinates // Select[
        Function[
          rLowerJoin Sin[-alpha] <= #[[2]] <= rUpperJoin Sin[alpha]
            && predicateWet[#[[1]], #[[2]]]
        ]
      ];
      boundaryCoordinates = boundaryCoordinates // SortBy[Last];
      (* Solve capillary BVP *)
      tNumerical = SolveLaplaceYoung[gamma, mesh, predicateWet];
      (* Compute height values *)
      heightValues = tNumerical @@@ boundaryCoordinates;
      (* Return *)
      {boundaryCoordinates, heightValues}
      , {meshStuff, meshStuffList}
    ] // Compress
  ]
]


(* ::Subsection:: *)
(*Geometric regions*)


(* ::Subsubsection:: *)
(*Wedge walls*)


wedgeWalls[alpha_] := {
  HalfLine[{0, 0}, XYPolar[1, alpha]],
  HalfLine[{0, 0}, XYPolar[1, -alpha]]
};


(* ::Subsection:: *)
(*Global styles for plots*)


contStyle = Directive[Thin, Black];
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
wallStyle = Directive[Thick, Purple];


upperStyle = Blue;
lowerStyle = Red;


(* ::Section:: *)
(*Parsing of boundary value problem*)


(* ::Text:: *)
(*Inner workings of SolveLaplaceYoung in LaplaceYoung.wl.*)


Module[
  {
    gamma = 1 Degree,
    mesh = ToElementMesh @ Disk[],
    prWet = True &,
    dummyForTrailingCommas
  },
  With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
    Module[{grad, iGrad, iDiv, k, stateData, propList},
      grad = Grad[#, {x, y}] &;
      iGrad = Inactive[Grad][#, {x, y}] &;
      iDiv = Inactive[Div][#, {x, y}] &;
      k = 1 / Sqrt[1 + # . #] & @ grad @ t[x, y];
      (* Parsing *)
      stateData =
        Quiet[
          NDSolve`ProcessEquations[
            {
              iDiv[-k * IdentityMatrix[2] . iGrad @ t[x, y]] ==
                - t[x, y]
                + NeumannValue[Cos[gamma], prWet[x, y]]
            }
            , t
            , Element[{x, y}, mesh]
          ]
          , {NDSolve`ProcessEquations::femibcnd}
        ];
      propList = {
        "ReactionCoefficients",
        "DiffusionCoefficients",
        Nothing
      };
      Table[
        {
          prop,
          stateData[[1]]["FiniteElementData"]["PDECoefficientData"][prop]
        }
        , {prop, propList}
      ] // TableForm
    ]
  ]
]


(* ::Section:: *)
(*Finite element mesh*)


(* ::Subsection:: *)
(*Sharp-cornered wedges*)


(* For each value of alpha: *)
Table[
  Module[{mesh, prWet, nElem},
    (* Import mesh *)
    {mesh, prWet} =
      Import[
        FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt"
      ] // Uncompress;
    nElem = Length @ mesh[[2, 1, 1]];
    (* Plot *)
    Show[
      mesh["Wireframe"],
      PlotLabel -> FString @ "{nElem} mesh elements",
      PlotOptions[Axes] // Evaluate
    ]
  ] // Ex @ FString @ "mesh/wedge_acute-mesh-apd-{apd}.png"
, {apd, apdValues}]


(* ::Subsection:: *)
(*Rounded wedges*)


(* (This is slow (~10 sec).) *)
Table[
  Module[
    {
      mesh, predicateWet, tipTheory, tipMesh,
      numElements,
      boundaryCoordinates, wetBoundaryCoordinates,
      dummyForTrailingCommas
    },
    (* Import mesh *)
    {mesh, predicateWet, tipTheory, tipMesh} =
      FString @ "mesh/wedge_acute-mesh-apd-{apd}-rad-{rad}.txt"
        // Import // Uncompress;
    numElements = Length @ mesh[[2, 1, 1]];
    (* Determine predicate-satisfying boundary coordinates *)
    boundaryCoordinates =
      Part[
        mesh["Coordinates"],
        mesh["BoundaryElements"][[1, 1]] // Flatten // Union
      ];
    wetBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateWet @@ # &
      ];
    (* Plot *)
    Show[
      mesh["Wireframe"],
      ListPlot[boundaryCoordinates, PlotStyle -> Directive[Red, PointSize[Medium]]],
      ListPlot[wetBoundaryCoordinates, PlotStyle -> Directive[Blue, PointSize[Medium]]],
      {}
      , PlotLabel -> Column @ {
          FString @ "{numElements} mesh elements",
          {"Theoretical tip location", tipTheory},
          {"Actual (mesh) tip location", tipMesh},
          {"Equal?", tipTheory == tipMesh},
          Nothing
        }
    ]
  ] // Ex @ FString @ "mesh/wedge_acute-mesh-apd-{apd}-rad-{rad}.png"
  , {apd, {60}}
  , {rad, Range[0, 3, 0.2] // Rest}
]


(* ::Section:: *)
(*Numerical solution*)


(* ::Subsection:: *)
(*Using built-in nonlinear solver*)


(* ::Subsubsection:: *)
(*Table*)


(* (This is slow (~1 min) from all the calls to Import.) *)
Join @@ Table[
  Table[
    Module[{tSol, messages, timing, tCorner},
      (* Import numerical solution, error messages and abs. time *)
      {tSol, messages, timing} = Import[
        FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"
      ] // Uncompress // Part[#, {1, 4, 5}] &;
      (* Filter out uninformative messages and make pretty *)
      messages = DeleteCases[
        messages,
        str_ /; StringContainsQ[str, "StringForm"]
      ];
      messages = StringRiffle[messages, "\n"];
      messages = messages /. {"" -> "\[HappySmiley]"};
      (* Compute corner height *)
      tCorner = tSol[0, 0];
      (* Build row of table *)
      {apd * Degree, gpd * Degree, tCorner, timing, messages}
    ]
  , {gpd, gpdValues}]
, {apd, apdValues}] // TableForm[#,
  TableAlignments -> {Left, Center},
  TableDepth -> 2,
  TableHeadings -> {
    None,
    {"\[Alpha]", "\[Gamma]", "T_corner", "Abs. time", "Messages"}
  },
  TableSpacing -> {3, 3}
] & // Ex["wedge_acute-solution-table.pdf"]


(* ::Subsubsection:: *)
(*Table for theoretical vs computed values*)


(* (This is slow (~6 sec) from all the calls to Import.) *)
Module[
  {
    apd, alpha,
    gamma, k,
    infiniteHeightQ, infiniteSlopeQ,
    theoreticalHeight, theoreticalSlope,
    tSol, tXDerivative, tYDerivative,
    computedHeight, computedSlope, slopeRelError,
    dummyForTrailingCommas
  },
  (* Representative value of wedge half-angle *)
  apd = apdRep;
  alpha = alphaRep;
  (* For each contact angle *)
  Table[
    (* Theoretical corner height and slope *)
    gamma = gpd Degree;
    k = Sin[alpha] / Cos[gamma];
    infiniteHeightQ = alpha < Pi/2 - gamma;
    infiniteSlopeQ = alpha <= Pi/2 - gamma;
    theoreticalHeight = If[infiniteHeightQ, Infinity, "finite"];
    theoreticalSlope = If[infiniteSlopeQ, Infinity, 1 / Sqrt[k^2 - 1] // N];
    (* Import numerical solution *)
    tSol =
      StringTemplate["solution/wedge_acute-solution-apd-``-gpd-``.txt"][apd, gpd]
        // Import // Uncompress // First;
    tXDerivative = Derivative[1, 0][tSol];
    tYDerivative = Derivative[0, 1][tSol];
    (* Computed corner height and slope *)
    computedHeight = tSol[0, 0];
    computedSlope = Sqrt[tXDerivative[0, 0]^2 + tYDerivative[0, 0]^2];
    slopeRelError = If[infiniteSlopeQ, "N/A", computedSlope / theoreticalSlope - 1];
    (* Build row of table *)
    {
      gamma,
      theoreticalHeight,
      computedHeight,
      theoreticalSlope,
      computedSlope,
      slopeRelError,
      Nothing
    }
    , {gpd, gpdValues}
  ]
    //
      TableForm[#
        , TableHeadings -> {
            None,
            {
              "\[Gamma]",
              "height",
              "(comp.)",
              "slope",
              "(comp.)",
              "(rel. error)",
              Nothing
            }
          }
      ] &
    // Ex["wedge_acute-height-slope-comparison.pdf"]
]


(* ::Subsubsection:: *)
(*Borderline case comparison with asymptotic*)


(* ::Text:: *)
(*From (34) of King et al. (1999) Quart. J. Mech. Appl. Math., 52 (1), 73\[Dash]97:*)
(*T ~ T_0 - 2 sqrt(x/T_0) + y^2/4 sqrt(T_0/x)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tSol, t0,
    tAsy,
    phiValues,
    dummyForTrailingCommas
  },
  (* Borderline case *)
  apd = 40;
  gpd = 90 - apd;
  (* Angles *)
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tSol =
    Import @ StringTemplate["solution/wedge_acute-solution-apd-``-gpd-``.txt"][apd, gpd]
      // Uncompress // First;
  (* Corner height *)
  t0 = tSol[0, 0];
  (* Asymptotic solution *)
  tAsy[x_, y_] := t0 - 2 Sqrt[x / t0] + y^2/4 Sqrt[t0/x];
  (* Plot discrepancy *)
  phiValues = Subdivide[0, alpha, 4];
  Table[
    Plot[
      {tSol[x, x Tan[phi]], tAsy[x, x Tan[phi]]} // Evaluate
      , {x, 0, 3}
      , ImageSize -> 240
      , PlotLabel -> HoldForm @ "\[Phi]" == phi
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Dashed}
      , PlotOptions[Axes] // Evaluate
    ]
    , {phi, phiValues}
  ]
] // Ex["wedge_acute-borderline-asymptotic-comparison-test.pdf"]


(* ::Section::Closed:: *)
(*Traced boundary plots (old)*)


(* ::Subsection:: *)
(*Same contact angle for tracing*)


(* ::Subsubsection:: *)
(*Hyperbolic traced boundaries*)


Module[
 {alpha, gamma,
  xMax, yMax,
  xyTerm, xyCont, xyTra
 },
  (* Constants *)
  alpha = alphaRep;
  gamma = gammaRep;
  (* Plot range *)
  xMax = 2 x0Same;
  yMax = xMax Tan[alpha];
  (* Import curves *)
  xyTerm = Import @ "same/wedge_acute_same-terminal.txt" // Uncompress;
  xyCont = Import @ "same/wedge_acute_same-contour-hyperbolic.txt" // Uncompress;
  xyTra = Import @ "same/wedge_acute_same-traced-hyperbolic.txt" // Uncompress;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}],
    (* Walls *)
    Graphics @ {wallStyle, wedgeWalls[alpha]},
    (* Non-viable domain and terminal curve *)
    Table[
      ParametricPlot[
        {
          Way[xy[[1]][s], 2 xMax, p],
          xy[[2]][s]
        },
        {s, DomainStart[xy], DomainEnd[xy]},
        {p, 0, 1},
        PlotStyle -> nonStyle,
        BoundaryStyle -> termStyle
      ]
    , {xy, {xyTerm}}],
    (* Local T-contour *)
    Table[
      ParametricPlot[
        xy[s] // Through,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> contStyle
      ]
    , {xy, xyCont}],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s] // Through // bothBranches // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> {upperStyle, lowerStyle}
      ]
    , {xy, xyTra}]
  ]
] // Ex @ FString @ StringJoin[
  "same/wedge_acute_same-traced-hyperbolic",
  "-apd-{apdRep}-gpd-{gpdRep}-gpdt-{gpdRep}.pdf"
]


(* ::Subsubsection:: *)
(*Generic traced boundaries*)


Module[
 {alpha, gamma,
  xMax, yMax,
  xyTerm, xyCont, xyTra
 },
  (* Constants *)
  alpha = alphaRep;
  gamma = gammaRep;
  (* Plot range *)
  xMax = 2 x0Same;
  yMax = xMax Tan[alpha];
  (* Import curves *)
  xyTerm = Import @ "same/wedge_acute_same-terminal.txt" // Uncompress;
  xyCont = Import @ "same/wedge_acute_same-contour-hyperbolic.txt" // Uncompress;
  xyTra = Import @ "same/wedge_acute_same-traced-general.txt" // Uncompress;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}],
    (* Walls *)
    Graphics @ {wallStyle, wedgeWalls[alpha]},
    (* Non-viable domain and terminal curve *)
    Table[
      ParametricPlot[
        {
          Way[xy[[1]][s], 2 xMax, p],
          xy[[2]][s]
        },
        {s, DomainStart[xy], DomainEnd[xy]},
        {p, 0, 1},
        PlotStyle -> nonStyle,
        BoundaryStyle -> termStyle
      ]
    , {xy, {xyTerm}}],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s] // Through // bothBranches // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> {upperStyle, lowerStyle}
      ]
    , {xy, xyTra}]
  ]
] // Ex @ FString @ StringJoin[
  "same/wedge_acute_same-traced-general",
  "-apd-{apdRep}-gpd-{gpdRep}-gpdt-{gpdRep}.pdf"
]


(* ::Section:: *)
(*Modifications (viable contours)*)


(* ::Subsection:: *)
(*Visualise contours*)


Module[
  {
    alpha, gamma,
    xCentreValues, hValues, xyContourList,
    tNumerical, vi,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {alpha, gamma} = {apdMod, gpdMod} Degree;
  (* Import contours *)
  {xCentreValues, hValues, xyContourList} =
    Import["modification/wedge_acute-modification-contours.txt"] // Uncompress;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
      // Uncompress // First;
  vi = Last @ ContactDerivativeList[tNumerical, gamma];
  (* Make plot *)
  Show[
    (* Wedge walls *)
    Graphics @ {wallStyle,
      HalfLine[{0, 0}, XYPolar[1, alpha]],
      HalfLine[{0, 0}, XYPolar[1, -alpha]]
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, 0, 3}
      , {y, -3, 3}
      , BoundaryStyle -> termStyle
      , PlotStyle -> nonStyle
    ],
    (* Contours *)
    Table[
      ParametricPlot[
        EvaluatePair[xy, s]
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> contStyle
      ]
      , {xy, xyContourList}
    ],
    (* Centreline starting points *)
    ListPlot[
      Table[{x, 0}, {x, xCentreValues}]
    ],
    {}
    , Frame -> True
    , PlotRange -> {{0, 2.5}, {-2, 2}}
  ]
] // Ex["modification/wedge_acute-modification-contours.pdf"]


(* ::Subsection:: *)
(*Check contours are actually contours*)


Module[
  {
    alpha, gamma,
    xCentreValues, hValues, xyContourList,
    tNumerical,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {alpha, gamma} = {apdMod, gpdMod} Degree;
  (* Import contours *)
  {xCentreValues, hValues, xyContourList} =
    Import["modification/wedge_acute-modification-contours.txt"] // Uncompress;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
      // Uncompress // First;
  (* Make plot *)
  Show[
    (* Pre-computed heights *)
    Plot[
      hValues
      , {s, -2, 2}
      , PlotStyle -> Black
      , PlotOptions[Axes] // Evaluate
    ],
    (* Evaluated heights along curves *)
    Table[
      Plot[
        tNumerical @@ EvaluatePair[xy, s]
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> Directive[Dotted, Green]
      ]
      , {xy, xyContourList}
    ],
    {}
    , AxesLabel -> {Italicise["s"], "Height"}
    , AxesOrigin -> {0, 0}
    , PlotRange -> Full
  ]
] // Ex["modification/wedge_acute-modification-contours-check.pdf"]


(* ::Subsection:: *)
(*Visualise finite element meshes*)


Module[
  {
    xCentreValues, hValues, xyContourList,
    numContours,
    meshes,
      xy, mesh,
      rUpperJoin, rLowerJoin,
    dummyForTrailingCommas
  },
  (* Import contours *)
  {xCentreValues, hValues, xyContourList} =
    Import["modification/wedge_acute-modification-contours.txt"] // Uncompress;
  numContours = Length[xyContourList];
  (* Import meshes *)
  meshes = Import["modification/wedge_acute-modification-meshes.txt"] // Uncompress;
  (* For each contour: *)
  Table[
    (* Make plot *)
    Show[
      (* Mesh *)
      mesh = meshes[[n, 1]];
      mesh["Wireframe"],
      (* Contour *)
      xy = xyContourList[[n]];
      ParametricPlot[
        EvaluatePair[xy, s]
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> Directive[Red, Opacity[0.5], AbsoluteThickness[3]]
      ],
      (* Joining points *)
      {rUpperJoin, rLowerJoin} = meshes[[n, {3, 4}]];
      ListPlot[
        {
          XYPolar[rUpperJoin, apdMod * Degree],
          XYPolar[rLowerJoin, -apdMod * Degree]
        }
        (*, PlotMarkers -> {Automatic, Medium}*)
      ],
      {}
      , ImageSize -> 480
    ]
      // Ex @ FString["modification/wedge_acute-modification-mesh-{n}.png"]
    , {n, numContours}
  ]
]


(* ::Subsection:: *)
(*Height comparison*)


Module[
  {
    xCentreValues, hValues, xyContourList,
    profilesStuff, profiles,
      xyList, heightList, yList,
    dummyForTrailingCommas
  },
  (* Import contours *)
  hValues =
    Import["modification/wedge_acute-modification-contours.txt"]
      // Uncompress // #[[2]] &;
  (* Import profiles *)
  profilesStuff =
    Import["modification/wedge_acute-modification-profiles.txt"]
      // Uncompress;
  (* Convert to (y, T) pairs *)
  profiles = Table[
    {xyList, heightList} = coordinatesAndHeights;
    yList = xyList[[All, 2]];
    {yList, heightList} // Transpose
    , {coordinatesAndHeights, profilesStuff}
  ];
  (* Make plot *)
  Show[
    (* Contour heights *)
    Plot[hValues, {y, -2, 2}
      , PlotStyle -> Dashed
    ],
    (* Profile heights *)
    ListPlot[profiles
      , Joined -> True
    ],
    {}
    , AxesLabel -> {Italicise["y"], "Height"}
    , AxesOrigin -> {-2, 0.5}
    , PlotRange -> Full
    , PlotOptions[Axes] // Evaluate
  ]
] // Ex["modification/wedge_acute-modification-contours-comparison.pdf"]


(* ::Section:: *)
(*Figure: finite element mesh (wedge_acute-mesh-*)*)


(* ::Subsection:: *)
(*Full*)


(* (This is slow (~6 sec).) *)
Module[
  {
    apd, alpha,
    mesh, meshWireframe,
    wallDistanceSimplified, polygonSimplified,
    posSimplified, meshWireframeSimplified,
    xMin, xMax, yMin, yMax, rMax,
    rTickValuesMajor, rTickValuesMinor, phiTickValues,
    tickLengthMajor, tickLengthMinor,
    tickTextStyle, axesTextStyle,
    plotFull,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  apd = 40;
  alpha = apd * Degree;
  (* Import mesh *)
  mesh = Import @ FString["mesh/wedge_acute-mesh-apd-{apd}.txt"] // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  {{xMin, xMax}, {yMin, yMax}} = mesh["Bounds"];
  rMax = xMax;
  (*
    Simplify mesh to reduce file size
    by replacing detailed portions with black polygons
  *)
  wallDistanceSimplified = 0.1;
  polygonSimplified = Graphics @ Polygon @ {
    {wallDistanceSimplified / Sin[alpha], 0},
    XYPolar[rMax, alpha] + XYPolar[wallDistanceSimplified, alpha - Pi/2],
    XYPolar[rMax, alpha],
    {0, 0},
    XYPolar[rMax, -alpha],
    XYPolar[rMax, -alpha] + XYPolar[wallDistanceSimplified, -alpha + Pi/2],
    Nothing
  };
  posSimplified =
    MeshWireframePositions[mesh,
      {x_, y_} /;
        x Sin[alpha] - y Cos[alpha] > wallDistanceSimplified
          &&
        x Sin[alpha] + y Cos[alpha] > wallDistanceSimplified
    ];
  meshWireframeSimplified = mesh["Wireframe"[posSimplified]];
  (* Ticks *)
  rTickValuesMajor = FindDivisions[{0, rMax}, 5];
  rTickValuesMinor = FindDivisions[{0, rMax}, 20] // Complement[#, rTickValuesMajor] &;
  phiTickValues = FindDivisions[{-apd, apd}, 8] Degree;
  tickLengthMajor = 0.01 rMax;
  tickLengthMinor = 1/2 tickLengthMajor;
  (* Plot full mesh *)
  tickTextStyle = Style[#, LabelSize["Tick"]] & @* LaTeXStyle;
  axesTextStyle = Style[#, LabelSize["Axis"]] & @* LaTeXStyle;
  plotFull = Show[
    meshWireframeSimplified, polygonSimplified,
    (* Radial ticks *)
    Graphics @ {
      Table[
        {
          Line @ {
            XYPolar[r, -alpha],
            XYPolar[r, -alpha] + XYPolar[tickLengthMajor, -alpha - Pi/2]
          },
          Text[
            r // tickTextStyle
            , XYPolar[r, -alpha] + XYPolar[3 tickLengthMajor, -alpha - Pi/2]
            , {0, 0.15} + If[r == 10, {0.5, 0}, {0, 0}]
          ],
          {}
        }
        , {r, rTickValuesMajor}
      ]
    },
    Graphics @ {
      Table[
        {
          Line @ {
            XYPolar[r, -alpha],
            XYPolar[r, -alpha] + XYPolar[tickLengthMinor, -alpha - Pi/2]
          },
          {}
        }
        , {r, rTickValuesMinor}
      ]
    },
    (* Azimuthal ticks *)
    Graphics @ {
      Table[
        {
          Line @ {
            XYPolar[rMax, phi],
            XYPolar[rMax + tickLengthMajor, phi]
          },
          Text[
            SeparatedRow["VeryThin"][phi / Degree, Magnify["\[Degree]", 1.2]] // tickTextStyle
            , XYPolar[rMax + 5.5 tickLengthMajor, phi]
            , Which[
                phi == 0, {0.3, -0.15},
                phi < 0, {-0.15, -0.3 + 0.15 (phi/alpha)},
                True, {0, 0}
              ]
          ],
          {}
        }
        , {phi, phiTickValues}
      ]
    },
    (* Radial axis label *)
    Graphics @ {
      Text[
        Italicise["r"] // axesTextStyle
        , XYPolar[rMax / 2, -alpha]
        , {4, 1.5}
      ]
    },
    (* Azimuthal axis label *)
    Graphics @ {
      Text[
        "\[Phi]" // axesTextStyle // Margined[{{0, 1}, {0, 0}}]
        , {rMax, 0}
        , {-6.5, -0.05}
      ],
      {}
    },
    {}
    , ImageSize -> 0.67 ImageSizeTextWidth
    , PlotRange -> All
    , PlotRangeClipping -> False
  ];
  plotFull // Ex["wedge_acute-mesh-full.pdf"]
]


(* ::Subsection:: *)
(*Detail*)


Module[
  {
    apd, alpha,
    mesh, meshWireframe,
    rMax, xMax, yMax,
    rTickValues, tickLength,
    tickTextStyle, axesTextStyle,
    plotDetail,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  apd = 40;
  alpha = apd * Degree;
  (* Import mesh *)
  mesh = Import @ FString["mesh/wedge_acute-mesh-apd-{apd}.txt"] // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  rMax = 0.25;
  {xMax, yMax} = XYPolar[rMax, alpha];
  (* Ticks *)
  rTickValues = FindDivisions[{0, rMax}, 4] // N;
  tickLength = 0.02 rMax;
  (* Plot mesh detail *)
  tickTextStyle = Style[#, LabelSize["Tick"]] & @* LaTeXStyle;
  axesTextStyle = Style[#, LabelSize["Axis"]] & @* LaTeXStyle;
  plotDetail = Show[
    meshWireframe,
    (* Radial ticks *)
    Graphics @ {
      Table[
        {
          Line @ {
            XYPolar[r, -alpha],
            XYPolar[r, -alpha] + XYPolar[tickLength, -alpha - Pi/2]
          },
          Text[
            r // tickTextStyle
            , XYPolar[r, -alpha] + XYPolar[5 tickLength, -alpha - Pi/2]
            , If[r == 0, {-0.5, -0.1}, {0.2, -0.15}]
          ],
          {}
        }
        , {r, rTickValues}
      ]
    },
    (* Radial axis label *)
    Graphics @ {
      Text[
        Italicise["r"] // axesTextStyle
        , XYPolar[rMax / 2, -alpha]
        , {2.7, 2.5}
      ],
      {}
    },
    {}
    , ImageSize -> 0.27 ImageSizeTextWidth
    , PlotRange -> {{0, xMax}, {-yMax, yMax}}
    , PlotRangePadding -> {0.2, {0.3, 0.2}} rMax
  ];
  plotDetail
] // Ex["wedge_acute-mesh-detail.pdf"]


(* ::Section:: *)
(*Figure: borderline case comparison with asymptotic  (wedge_acute-borderline-asymptotic-comparison-*)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, t0,
    tAsymptotic,
    styleCommon, styleNumerical, styleAsymptotic,
    phiValues, phiCases, case,
    dummyForTrailingCommas
  },
  (* Borderline case *)
  apd = 40;
  gpd = 90 - apd;
  (* Angles *)
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Corner height *)
  t0 = tNumerical[0, 0];
  (* Asymptotic solution *)
  tAsymptotic[x_, y_] := t0 - 2 Sqrt[x / t0] + y^2/4 Sqrt[t0/x];
  (* Plot styles *)
  styleCommon = Directive[Black, GeneralStyle["DefaultThick"]];
  styleNumerical = styleCommon;
  styleAsymptotic = Directive[styleCommon, GeneralStyle["Dashed"]];
  (* Plot discrepancy *)
  phiValues = {0, alpha};
  phiCases = {"symmetry", "wall"};
  Table[
    case = AssociationThread[phiValues -> phiCases][phi];
    Plot[
      {tNumerical[x, x Tan[phi]], tAsymptotic[x, x Tan[phi]]} // Evaluate
      , {x, 0, 3}
      , AxesLabel -> {
          Italicise["x"] // Margined @ {{0, 1}, {5, 0}},
          Italicise["T"]
        }
      , ImageSize -> 0.45 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , PlotRange -> {0, All}
      , PlotStyle -> {styleNumerical, styleAsymptotic}
      , PlotOptions[Axes] // Evaluate
      , TicksStyle -> LabelSize["Tick"]
    ] // Ex @ FString["wedge_acute-borderline-asymptotic-comparison-{case}.pdf"]
    , {phi, phiValues}
  ]
]


(* ::Section:: *)
(*Figure: known solution  (wedge_acute-solution-*)*)


Module[
  {
    apd, alpha,
    gpdValues, gpdCases,
    rMax, xArcCorner, yArcCorner,
    case,
    tNumerical,
    x, y, mesh,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 40;
  alpha = apd * Degree;
  (* Contact angles *)
  gpdValues = {50, 55};
  gpdCases = {"borderline", "moderate"};
  (* Plot range *)
  rMax = 4;
  {xArcCorner, yArcCorner} = XYPolar[rMax, alpha];
  (* Make plots *)
  Table[
    case = AssociationThread[gpdValues -> gpdCases][gpd];
    (* Import numerical solution *)
    tNumerical =
      Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    (* Plot *)
    mesh = tNumerical["ElementMesh"];
    Show[
      (* Numerical solution *)
      Plot3D[
        tNumerical[x, y], Element[{x, y}, mesh]
        , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
        , AxesLabel -> {
            Italicise["x"],
            Italicise["y"],
            Italicise["T"] // Margined[{{0, 0.2}, {0, 1}} ImageSizeCentimetre]
          }
        , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
        , Boxed -> {Back, Bottom, Left}
        , BoxRatios -> Automatic
        , Filling -> 0
        , FillingStyle -> BoundaryTracingStyle["Solution3D"]
        , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
        , Lighting -> GeneralStyle["AmbientLighting"]
        , MeshStyle -> BoundaryTracingStyle["Edge3D"]
        , PlotPoints -> 20
        , PlotRange -> {0, 2.2}
        , PlotStyle -> BoundaryTracingStyle["Solution3D"]
        , TicksStyle -> LabelSize["Tick"]
        , RegionFunction -> Function[{x, y},
            RPolar[x, y] < rMax
          ]
        , ViewPoint -> {1.5, -2.5, 1.4}
      ],
      (* Manually drawn solution base edges *)
      Graphics3D @ {
        BoundaryTracingStyle["Edge3D"],
        Line @ {{0, 0, 0}, {xArcCorner, -yArcCorner, 0}},
        {}
      },
      ParametricPlot3D[
        rMax {Cos[phi], Sin[phi], 0}
        , {phi, -alpha, alpha}
        , PlotStyle -> Directive[Black, BoundaryTracingStyle["Edge3D"]]
      ],
      (* Manually drawn solution vertical edges *)
      Graphics3D @ {
        BoundaryTracingStyle["Edge3D"],
        Line @ {
          {0, 0, tNumerical[0, 0]},
          {0, 0, 0}
        },
        Line @ {
          {xArcCorner, -yArcCorner, tNumerical[xArcCorner, -yArcCorner]},
          {xArcCorner, -yArcCorner, 0}
        },
        {}
      },
      {}
      , ImageSize -> 0.45 ImageSizeTextWidth
    ]
      // Ex[FString["wedge_acute-solution-{case}.png"]
        , Background -> None
        , ImageResolution -> 4 BasicImageResolution
      ]
    , {gpd, gpdValues}
  ]
]


(* ::Section:: *)
(*Figure: generic traced boundaries  (wedge_acute-traced-boundaries)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    rStartList, xyStartList,
    xStartTerminal, yStartTerminal,
    xyTracedList,
    plot, legend,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.5 xCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundaries (upper branch) *)
  (*
    NOTE: wall integration doesn't work
    because the traced boundary ventures slightly outside the wedge domain
    and Mathematica complains about Indeterminate values.
    So the wall boundaries are drawn manually at the end.
  *)
  rStartList = Subdivide[rMax, 10] // Rest;
  xyStartList = Table[XYPolar[r, -alpha], {r, rStartList}];
    (* Append a starting point along the terminal curve *)
    xStartTerminal = Way[xCritical, xMax, 0.6];
    yStartTerminal =
      SeekRoot[vi[xStartTerminal, #] &, {0, xStartTerminal Tan[alpha]}, 5];
    xyStartList = xyStartList // Append @ {xStartTerminal, yStartTerminal};
  xyTracedList =
    Table[
      ContactTracedBoundary[derList][xyStart, 0, {0, 2 rMax}
        , -1, 1
        , 0, Function[{x, y}, x > xMaxMore]
      ]
      , {xyStart, xyStartList}
    ];
  (* Plot *)
  plot = Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      HalfLine[{0, 0}, XYPolar[1, alpha]],
      HalfLine[{0, 0}, XYPolar[1, -alpha]],
      {}
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xCritical, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // IncludeYReflection
          // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedList}
    ],
    (* Traced boundary walls *)
    ParametricPlot[
      {x, x Tan[alpha]}
        // IncludeYReflection
        // Evaluate
      , {x, 0, xMaxMore}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    {}
  ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced", "Terminal"},
      {"traced boundary", "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ];
  legendRegions =
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ];
  (* Export *)
  {
    Show[plot
      , ImageSize -> 0.4 ImageSizeTextWidth
    ] // Ex["wedge_acute-traced-boundaries.pdf"]
    ,
    GraphicsColumn[Join[legendCurves, legendRegions]
      , Alignment -> Left
      , ImageSize -> 0.3 ImageSizeTextWidth
      , ItemAspectRatio -> 0.13
      , Spacings -> {0, 0}
    ] // Ex["wedge_acute-traced-boundaries-legend.pdf"]
  }
]


(* ::Section:: *)
(*Figure: generic traced boundary branches  (wedge_acute-traced-boundaries-*-branch)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    rStartList, xyStartList,
    xStartTerminal, yStartTerminal,
    xyTracedList,
    branchCases, branchYSigns, ySign,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.5 xCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundaries (upper branch) *)
  (*
    NOTE: wall integration doesn't work
    because the traced boundary ventures slightly outside the wedge domain
    and Mathematica complains about Indeterminate values.
    So the wall boundaries are drawn manually at the end.
  *)
  rStartList = Subdivide[rMax, 12] // Rest;
  xyStartList = Table[XYPolar[r, -alpha], {r, rStartList}];
    (* Append a starting point along the terminal curve *)
    xStartTerminal = Way[xCritical, xMax, 0.55];
    yStartTerminal =
      SeekRoot[vi[xStartTerminal, #] &, {0, xStartTerminal Tan[alpha]}, 5];
    xyStartList = xyStartList // Append @ {xStartTerminal, yStartTerminal};
  xyTracedList =
    Table[
      ContactTracedBoundary[derList][xyStart, 0, {0, 2 rMax}
        , -1, 1
        , 0, Function[{x, y}, x > xMaxMore]
      ]
      , {xyStart, xyStartList}
    ];
  (* Plot *)
  branchCases = {"upper", "lower"};
  branchYSigns = {1, -1};
  Table[
    ySign = AssociationThread[branchCases -> branchYSigns][case];
    Show[
      EmptyFrame[{0, xMax}, {-yMax, yMax}
        , Frame -> None
      ],
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        HalfLine[{0, 0}, XYPolar[1, alpha]],
        HalfLine[{0, 0}, XYPolar[1, -alpha]],
        {}
      },
      (* Non-viable domain *)
      RegionPlot[
        vi[x, y] < 0
        , {x, xCritical, xMaxMore}
        , {y, -yMaxMore, yMaxMore}
        , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
        , PlotPoints -> 9
        , PlotStyle -> BoundaryTracingStyle["NonViable"]
      ],
      (* Traced boundaries *)
      Table[
        ParametricPlot[
          EvaluatePair[xy, s, ySign] // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
        , {xy, xyTracedList}
      ],
      (* Traced boundary walls *)
      ParametricPlot[
        {x, ySign * x Tan[alpha]} // Evaluate
        , {x, 0, xMaxMore}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
      {}
      , ImageSize -> 0.45 * 0.7 ImageSizeTextWidth
    ] // Ex @ FString["wedge_acute-traced-boundaries-{case}-branch.pdf"]
    , {case, branchCases}
  ]
]


(* ::Section:: *)
(*Figure: traced boundaries, patched (wedge_acute-traced-boundaries-patched)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    derList, p, q, grad2, f, vi,
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
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.5 xCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  normalVectorLength = 1/3 xMax;
  (* Location of corners (increasing in y-coordinate) *)
  cornerList = Association[
    1 -> {{0.45, 0.2}},
    2 -> {{0.4, 0}, {0.45, 0.2}},
    3 -> {{0.6, -0.5}, {0.4, 0.25}, {0.5, 0.4}},
    4 -> {{0.5, -0.38}, {0.35, 0}, {0.4, 0.25}, {0.5, 0.4}},
    Nothing
  ];
  numPlots = Length[cornerList];
  (* Compute traced boundaries therethrough *)
  sRange = 2 rMax {-1, 1};
  Table[
    (* Upper branch *)
    xyTracedUpperList[n] =
      Table[
        ContactTracedBoundary[derList][corner, 0, sRange
          , -1, 1
          , 0, Function[{x, y}, x > xMaxMore]
        ]
        , {corner, cornerList[n]}
      ];
    xyTracedLowerList[n] =
      Table[
        ContactTracedBoundary[derList][corner, 0, sRange
          , 1, -1
          , 0, Function[{x, y}, x > xMaxMore]
        ]
        , {corner, cornerList[n]}
      ];
    (* Intersections *)
    numCorners = Length @ cornerList[n];
    Table[
      (* Starting arc length for lower branch *)
      sIntersectionLower[n][k] =
        If[k == 1
          ,
          xyLower = xyTracedLowerList[n][[k]];
          DomainEnd[xyLower]
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
          DomainEnd[xyUpper]
          ,
          xyUpper = xyTracedUpperList[n][[k]];
          xyLower = xyTracedLowerList[n][[k + 1]];
          SeekParametricIntersection[xyUpper, xyLower] // First
        ];
      , {k, numCorners}
    ];
    (* Normal vector *)
    {xTracedForNormal, yTracedForNormal} = xyTracedUpperList[n][[1]];
    sNormal = Way[0, sIntersectionUpper[n][1], If[n == 1, 1/4, 1/2]];
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
          , {1., -0.8}
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
    , Spacings -> 0
  ]
] // Ex["wedge_acute-traced-boundaries-patched.pdf"]


(* ::Section:: *)
(*Figure: terminal-points (wedge_acute-terminal-points)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    p, q, grad2, f, vi,
    xMin, xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    xyStartList, xyContourList,
    xyTerminal,
    xyContourOrdinary, sOrdinary, xyOrdinary,
    textStyle, textStyleBracket, textVerticalShift,
    plot,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMin = xCritical / 2;
  xMax = 2 xCritical;
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* T-contours *)
  xyStartList = Table[{x, 0}, {x, xCritical + 0.07 Range[0, 3]}];
  xyContourList =
    Table[
      ContourByArcLength[tNumerical][xyStart, 0, 2 rMax {-1, 1}, 1, #1 > xMaxMore &]
      , {xyStart, xyStartList}
    ];
  (* Terminal curve *)
  xyTerminal =
    ContourByArcLength[vi][{xCritical, 0}, 0, 2 rMax {-1, 1}, 1, #1 > xMaxMore &];
  (* Orindary terminal point to be shown *)
  xyContourOrdinary = xyContourList[[2]];
  sOrdinary =
    First @ SeekParametricIntersection[
      xyContourOrdinary, xyTerminal,
      {0, DomainEnd[xyContourOrdinary]}, {0, DomainEnd[xyTerminal]}
    ];
  xyOrdinary = xyContourOrdinary[sOrdinary] // Through;
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textVerticalShift = -0.25;
  (* Plot *)
  plot = Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , Frame -> None
      , PlotRangePadding -> {Automatic, None}
    ],
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> None
      , PlotPoints -> 7
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
      , PlotPoints -> 7
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xCritical, 0}
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
        {xCritical, 0},
        {1.5, textVerticalShift}
      ] // textStyle,
      Text[
        "critical",
        {xCritical, 0},
        {-1.4, textVerticalShift}
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
        {-1.33, textVerticalShift}
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
      {Row @ {Italicise["T"], "\[Hyphen]contour"}, "terminal curve"}
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
    , ItemAspectRatio -> 1.7
    , ImageSize -> 0.65 ImageSizeTextWidth
    , Spacings -> {-0.1 ImageSizeTextWidth, 0}
  ]
] // Ex["wedge_acute-terminal-points.pdf"]


(* ::Section:: *)
(*Figure: traced boundaries through (x_0, 0) (wedge_acute-traced-boundaries-hyperbolic-*)*)


(* ::Subsection:: *)
(*Both*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    xy,
    textStyle, textStyleBracket,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.7 xCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary (upper branch) *)
  xy =
    ContactTracedBoundary[derList][
      {xCritical, 0}, 0, 2 rMax {-1, 1}
      , -1, 1
      , -Infinity, Function[{x, y}, x > xMaxMore]
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
      , ImageSize -> 0.45 * 0.7 ImageSizeTextWidth
      , PlotRangePadding -> {{0.04, Automatic}, Automatic}
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMaxMore {1, Tan[alpha]},
        {0, 0},
        xMaxMore {1, -Tan[alpha]}
      }
    },
    (* Traced boundaries *)
    ParametricPlot[
      xy[s]
        // Through
        // IncludeYReflection
        // Evaluate
      , {s, DomainStart[xy], DomainEnd[xy]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xCritical, 0}
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
        {xCritical, 0},
        {1.45, -0.1}
      ] // textStyle,
      {}
    },
    {}
  ]
] // Ex["wedge_acute-traced-boundaries-hyperbolic-both.pdf"]


(* ::Subsection:: *)
(*Rounding*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    xy,
    textStyle, textStyleBracket,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.7 xCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary (upper branch) *)
  xy =
    ContactTracedBoundary[derList][
      {xCritical, 0}, 0, 2 rMax {0, 1}
      , -1, 1
      , -Infinity, Function[{x, y}, x > xMaxMore]
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
      , ImageSize -> 0.45 * 0.7 ImageSizeTextWidth
      , PlotRangePadding -> {{0.04, Automatic}, Automatic}
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMaxMore {1, Tan[alpha]},
        {0, 0},
        xMaxMore {1, -Tan[alpha]}
      }
    },
    (* Traced boundaries *)
    ParametricPlot[
      xy[s]
        // Through
        // IncludeYReflection
        // Evaluate
      , {s, DomainStart[xy], DomainEnd[xy]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xCritical, 0}
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
        {xCritical, 0},
        {1.45, -0.1}
      ] // textStyle,
      {}
    },
    {}
  ]
] // Ex["wedge_acute-traced-boundaries-hyperbolic-rounding.pdf"]


(* ::Section:: *)
(*Figure: wall coordinates (wedge_acute-wall-coordinates)*)


Module[
  {
    alpha,
    wedgeRMax, wedgeXMax, wedgeYMax,
    clearanceProportion, xMin, xMax, yMax,
    xArrowBase, yArrowBase, arrowLength,
    orthogonalityMarkerLength, orthogonalityMarkerStyle,
    angleMarkerRadius, textSize, textStyle,
    dummyForTrailingCommas
  },
  (* Wedge *)
  alpha = 40 Degree;
  wedgeRMax = 1;
  {wedgeXMax, wedgeYMax} = XYPolar[wedgeRMax, alpha];
  (* Plot range *)
  clearanceProportion = 1/10;
  xMin = Way[0, wedgeXMax, -clearanceProportion];
  xMax = Way[0, wedgeXMax, 1 + clearanceProportion];
  yMax = Way[0, wedgeYMax, 1 + clearanceProportion];
  (* Reference point (base) of wall coordinate arrows *)
  {xArrowBase, yArrowBase} = XYPolar[0.6 wedgeRMax, alpha];
  arrowLength = 0.3 wedgeRMax;
  (* Orthogonality marker *)
  orthogonalityMarkerLength = 1/14 wedgeXMax;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[None]];
  (* Etc. *)
  angleMarkerRadius = 0.15 wedgeRMax;
  textSize = LabelSize["Label"];
  textStyle = Style[#, textSize] & @* LaTeXStyle;
  (* Diagram *)
  Show[
    EmptyAxes[{xMin, xMax}, {-yMax, yMax}
      , AspectRatio -> Automatic
      , ImageSize -> 0.3 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ textSize
      , Method -> {"AxesInFront" -> False}
      , Ticks -> None
    ],
    (* Wedge walls *)
    Graphics @ {
      BoundaryTracingStyle["Wall"],
      Line @ {{wedgeXMax, wedgeYMax}, {0, 0}, {wedgeXMax, -wedgeYMax}}
    },
    (* Wedge half-angle *)
    Graphics @ {
      Circle[{0, 0}, angleMarkerRadius, {0, alpha}]
    },
    Graphics @ {
      Text[
        "\[Alpha]" // textStyle
        , XYPolar[angleMarkerRadius, alpha / 2]
        , {-2.5, -0.3}
      ]
    },
    (* Coordinate arrows *)
    Graphics @ {
      Directive[GeneralStyle["Thick"], Arrowheads[0.1]],
      Arrow @ {
        {xArrowBase, yArrowBase},
        {xArrowBase, yArrowBase} + XYPolar[arrowLength, alpha - Pi/2]
      },
      Text[
        "\[Xi]" // textStyle
        , {xArrowBase, yArrowBase} + XYPolar[arrowLength, alpha - Pi/2]
        , {-2.25, -1}
      ],
      Arrow @ {
        {xArrowBase, yArrowBase},
        {xArrowBase, yArrowBase} + XYPolar[arrowLength, alpha]
      },
      Text[
        "\[Eta]" // textStyle
        , {xArrowBase, yArrowBase} + XYPolar[arrowLength, alpha]
        , {-1.25, 0.75}
      ],
      {}
    },
    Graphics @ {orthogonalityMarkerStyle,
      Rectangle[
        {xArrowBase, yArrowBase},
        {xArrowBase, yArrowBase} + orthogonalityMarkerLength {1, 1}
      ] // Rotate[#, -(Pi/2 - alpha), {xArrowBase, yArrowBase}] &
    },
    {}
  ]
] // Ex["wedge_acute-wall-coordinates.pdf"]


(* ::Section:: *)
(*Figure: traced boundaries, different contact angle   (wedge_acute-traced-boundaries-different-angle-*)*)


(* ::Subsection:: *)
(*Main plots*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracingValues, caseNameList,
    gammaTracing,
    xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    rMaxViable,
    rStartList, xyStartList,
    xyTracedList,
    plot,
    caseName,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 40;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Different tracing contact angle *)
  gpdTracingValues = {50, 70};
  caseNameList = {"less", "more"};
  Table[
    gammaTracing = gpdTracing * Degree;
    (* Critical terminal point x_0 *)
    xCritical = x0[tNumerical, gammaTracing];
    (* Plot range *)
    xMax = 2 xCritical;
    yMax = xMax Tan[alpha];
    rMax = RPolar[xMax, yMax];
    (* Plot range but more *)
    more = 0.1;
    xMaxMore = xMax + more;
    yMaxMore = yMax + more;
    (* Derivative list for boundary tracing *)
    {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
    (* Starting points *)
    If[
      gammaTracing < gamma,
        rStartList = Subdivide[xCritical, 4] // Rest;
        xyStartList = Table[XYPolar[r, -alpha], {r, rStartList}];
      ,
        rStartList = Subdivide[rMax, 9] // Rest;
        xyStartList = Join[
          Table[XYPolar[r, +alpha], {r, rStartList}],
          Table[XYPolar[r, -alpha], {r, rStartList}],
          {{10^-6, 0}},
          {}
        ];
    ];
    (* Traced boundaries (upper branch) *)
    xyTracedList =
      Table[
        ContactTracedBoundary[derList][xyStart, 0, {0, 2 rMax}
          , -1, 1
          , 0, Function[{x, y}, x > xMaxMore]
        ]
        , {xyStart, xyStartList}
      ];
    (* Plot *)
    plot =
      Show[
        EmptyFrame[{0, xMax}, {-yMax, yMax}
          , FrameLabel -> {
              Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
              Italicise["y"]
            }
          , FrameTicksStyle -> LabelSize["Tick"]
          , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
        ],
        (* Non-viable domain *)
        If[
          gammaTracing < gamma,
          {
            RegionPlot[
               vi[x, y] < 0
               , {x, xCritical, xMaxMore}
               , {y, -yMaxMore, yMaxMore}
               , BoundaryStyle -> None
               , PlotPoints -> 7
               , PlotStyle -> BoundaryTracingStyle["NonViable"]
             ],
            ContourPlot[
               vi[x, y] == 0
               , {x, xCritical, xMaxMore}
               , {y, -yMaxMore, yMaxMore}
               , ContourLabels -> None
               , ContourStyle -> BoundaryTracingStyle["Terminal"]
               , PlotPoints -> 7
             ]
           },
           {
             RegionPlot[
               vi[x, y] < 0
               , {x, xCritical, xMaxMore}
               , {y, -yMaxMore, yMaxMore}
               , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
               , PlotPoints -> 7
               , PlotStyle -> BoundaryTracingStyle["NonViable"]
             ]
           }
        ],
        (* Wedge walls *)
        Graphics @ {BoundaryTracingStyle["Wall"],
          Line @ {
            xMaxMore {1, Tan[alpha]},
            {0, 0},
            xMaxMore {1, -Tan[alpha]}
          }
        },
        (* Traced boundaries *)
        Table[
          ParametricPlot[
            xy[s]
              // Through
              // IncludeYReflection
              // Evaluate
            , {s, DomainStart[xy], DomainEnd[xy]}
            , PlotPoints -> 2
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ]
          , {xy, xyTracedList}
        ],
        {}
      ];
    (* Export *)
    caseName = AssociationThread[gpdTracingValues, caseNameList][gpdTracing];
    Show[plot
      , ImageSize -> 0.4 ImageSizeTextWidth
    ] // Ex @ FString["wedge_acute-traced-boundaries-different-angle-{caseName}.pdf"]
    , {gpdTracing, gpdTracingValues}
  ]
]


(* ::Subsection:: *)
(*Legend*)


Module[
  {
    legendLabelStyle,
    row,
    dummyForTrailingCommas
  },
  legendLabelStyle = LatinModernLabelStyle[LabelSize["Legend"]];
  row = Join[
    CurveLegend[
      BoundaryTracingStyle /@ {"Wall", "Traced", "Terminal"},
      {"wedge wall", "traced boundary", "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ],
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ] // Margined @ {{0, 0}, {0, -3}} /@ # &,
    {}
  ];
  GraphicsGrid[{row}
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.15
    , Spacings -> {{0, -0.125, 0.1, 0.} ImageSizeTextWidth, 0}
  ]
] // Ex["wedge_acute-traced-boundaries-different-angle-legend.pdf"]


(* ::Section:: *)
(*Figure: terminal , hyperbolic, etc., different contact angle*)
(*(wedge_acute-terminal-points-different-angle)*)
(*(wedge_acute-traced-boundaries-hyperbolic-different-angle)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracing, gammaTracing,
    xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    xOrdinary, yOrdinary,
    xyContourList,
    textStyle, textStyleBracket, textVerticalShift,
    plot,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 40;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Different tracing contact angle *)
  gpdTracing = 70;
  gammaTracing = gpdTracing * Degree;
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gammaTracing];
  (* Plot range *)
  xMax = 2 xCritical;
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.1;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
  (* Ordinary terminal point *)
  yOrdinary = Way[0, yMax, 0.36];
  xOrdinary = SeekRoot[vi[#, yOrdinary] &, {xCritical, xMax}, 5];
  (* T-contours *)
  xyContourList =
    Table[
      ContourByArcLength[tNumerical][
        xyInitial, 0, 2 rMax {-1, 1}
        , 1, Function[{x, y}, x > xMaxMore]
      ]
      , {xyInitial, {{xCritical, 0}, {xOrdinary, yOrdinary}}}
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textVerticalShift = -0.2;
  (* Plot *)
  plot =
    Show[
      EmptyFrame[{0, xMax}, {-yMax, yMax}
        , Frame -> None
        , PlotRangePadding -> {{0.1, Automatic}, Automatic}
      ],
      (* Non-viable domain *)
      RegionPlot[
        vi[x, y] < 0
        , {x, xCritical, xMaxMore}
        , {y, -yMaxMore, yMaxMore}
        , BoundaryStyle -> None
        , PlotPoints -> 8
        , PlotStyle -> BoundaryTracingStyle["NonViable"]
      ],
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          xMaxMore {1, Tan[alpha]},
          {0, 0},
          xMaxMore {1, -Tan[alpha]}
        }
      },
      (* Local T-contours *)
      Table[
        ParametricPlot[
          xy[s] // Through
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["ContourPlain"]
        ]
        , {xy, xyContourList}
      ],
      (* Terminal curve *)
      ContourPlot[
        vi[x, y] == 0
        , {x, xCritical, xMaxMore}
        , {y, -yMaxMore, yMaxMore}
        , ContourStyle -> BoundaryTracingStyle["Terminal"]
        , PlotPoints -> 14
      ],
      (* Critical terminal point (x_0, 0) *)
      Graphics @ {
        GeneralStyle["Point"],
        Point @ {xCritical, 0}
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
          {xCritical, 0},
          {1.5, textVerticalShift}
        ] // textStyle,
        Text[
          "critical",
          {xCritical, 0},
          {-1.4, textVerticalShift}
        ] // textStyle,
        {}
      },
      (* Ordinary terminal point *)
      Graphics @ {
        GeneralStyle["Point"],
        Point @ {xOrdinary, yOrdinary}
      },
      Graphics @ {
        Text[
          "ordinary",
          {xOrdinary, yOrdinary},
          {-1.15, 0.6}
        ] // textStyle,
        {}
      },
      {}
    ];
  Show[plot
    , ImageSize -> 0.5 * 0.63 ImageSizeTextWidth
  ] // Ex["wedge_acute-terminal-points-different-angle.pdf"]
]


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracing, gammaTracing,
    xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    xy,
    textStyle, textStyleBracket, textVerticalShift,
    plot,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 40;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Different tracing contact angle *)
  gpdTracing = 70;
  gammaTracing = gpdTracing * Degree;
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gammaTracing];
  (* Plot range *)
  xMax = 2 xCritical;
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.1;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
  (* Traced boundary (upper branch) *)
  xy =
    ContactTracedBoundary[derList][
      {xCritical, 0}, 0, 2 rMax {-1, 1}
      , -1, 1
      , -Infinity, Function[{x, y}, x > xMaxMore]
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  textVerticalShift = -0.2;
  (* Plot *)
  plot =
    Show[
      EmptyFrame[{0, xMax}, {-yMax, yMax}
        , Frame -> None
        , PlotRangePadding -> {{0.1, Automatic}, Automatic}
      ],
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          xMaxMore {1, Tan[alpha]},
          {0, 0},
          xMaxMore {1, -Tan[alpha]}
        }
      },
      (* Traced boundaries *)
      ParametricPlot[
        xy[s]
          // Through
          // IncludeYReflection
          // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
      (* Critical terminal point (x_0, 0) *)
      Graphics @ {
        GeneralStyle["Point"],
        Point @ {xCritical, 0}
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
          {xCritical, 0},
          {1.5, textVerticalShift}
        ] // textStyle,
        {}
      },
      {}
    ];
  Show[plot
    , ImageSize -> 0.5 * 0.63 ImageSizeTextWidth
  ] // Ex["wedge_acute-traced-boundaries-hyperbolic-different-angle.pdf"]
]


(* ::Section:: *)
(*Figure: different contact angle offset*)
(*(wedge_acute-different-angle-offset-*)*)


(* ::Subsection:: *)
(*Top view*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracing, gammaTracing,
    xCritical,
    xEnd, yEnd, sMax,
    derList, p, q, grad2, f, vi,
    xy,
    d, xArrowBase, yArrowBase,
    xMaxWall, yMaxWall,
    xMaxOffset, yMaxOffset,
    arrowLength,
    orthogonalityMarkerLength, orthogonalityMarkerStyle,
    wallThicknessCorrection,
    textStyleLabel, textStylePoint, textStylePointBracket,
    plot,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 40;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Different tracing contact angle *)
  gpdTracing = 70;
  gammaTracing = gpdTracing * Degree;
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gammaTracing];
  (* Ending x for traced boundary *)
  xEnd = 1.9 xCritical;
  sMax = 2 xEnd;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
  (* Traced boundary (upper branch) *)
  xy =
    ContactTracedBoundary[derList][
      {xCritical, 0}, 0, {0, sMax}
      , -1, 1
      , -Infinity, Function[{x, y}, x > xEnd]
    ];
  {xEnd, yEnd} = xy @ DomainEnd[xy] // Through;
  (* Half-plane offset distance *)
  d = DHalfPlane[gamma, gammaTracing];
  (* Base of wall coordinate (xi) arrow *)
  xArrowBase = xEnd + d Cos[alpha + Pi/2];
  yArrowBase = xArrowBase Tan[alpha];
  (* End point of wall to be shown *)
  {xMaxWall, yMaxWall} = 1.15 {xArrowBase, yArrowBase};
  (* End point of offset wedge to be shown *)
  {xMaxOffset, yMaxOffset} = {xMaxWall, yMaxWall} + XYPolar[d, alpha - Pi/2];
  (* Etc. *)
  arrowLength = 3 d;
  orthogonalityMarkerLength = 0.35 d;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[None]];
  wallThicknessCorrection = 0.23;
  (* Text style *)
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStylePoint = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStylePointBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Plot *)
  plot =
    Show[
      (* Coordinate arrow *)
      Graphics @ {
        Directive[GeneralStyle["Thick"], Arrowheads[0.1]],
        Arrow @ {
          {xArrowBase, yArrowBase},
          {xArrowBase, yArrowBase} + XYPolar[arrowLength, alpha - Pi/2]
        },
        Text[
          "\[Xi]" // textStyleLabel
          , {xArrowBase, yArrowBase} + XYPolar[arrowLength, alpha - Pi/2]
          , {-2.25, -1}
        ],
        {}
      },
      Graphics @ {orthogonalityMarkerStyle,
        Rectangle[
          {xArrowBase, yArrowBase},
          {xArrowBase, yArrowBase}
            + orthogonalityMarkerLength {1 + wallThicknessCorrection, -1}
        ] // Rotate[#, -(Pi/2 - alpha), {xArrowBase, yArrowBase}] &
      },
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          {xMaxWall, +yMaxWall},
          {0, 0},
          {xMaxWall, -yMaxWall}
        }
      },
      Graphics @ {
        Text[
          "\[Xi]" == SeparatedRow[]["\[VeryThinSpace]", 0] // textStyleLabel
          , {xMaxWall, +yMaxWall}
          , {-1.7, -0.2}
          , XYPolar[1, alpha]
        ]
      },
      Graphics @ {
        Text[
          "\[Xi]" ==
            Row @ {
              Italicise["d"],
              "(" // textStylePointBracket,
              "\[NegativeVeryThinSpace]",
              "\[Gamma]",
              ",\[ThinSpace]",
              Subscript["\[Gamma]", Style["\[NegativeVeryThinSpace]\[Bullet]", Magnification -> 1.8]],
              ")" // textStylePointBracket
            }
            // textStyleLabel
          , {xMaxOffset, +yMaxOffset}
          , {-1.35, -0.2}
          , XYPolar[1, alpha]
        ]
      },
      (* Offset wedge *)
      (* (separate lines to ensure dotting at the corner) *)
      Graphics @ {BoundaryTracingStyle["Contour"],
        Line @ {
          {d / Sin[alpha], 0},
          {xMaxOffset, +yMaxOffset} * 1.02 (* optical correction *)
        },
        Line @ {
          {d / Sin[alpha], 0},
          {xMaxOffset, -yMaxOffset}
        },
        {}
      },
      (* Traced boundaries *)
      ParametricPlot[
        xy[s]
          // Through
          // IncludeYReflection
          // Evaluate
        , {s, DomainStart[xy], 0.98 DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
      (* Critical terminal point (x_0, 0) *)
      Graphics @ {
        GeneralStyle["Point"],
        Point @ {xCritical, 0}
      },
      Graphics @ {
        Text[
          Row @ {
            "(" // textStylePointBracket,
            "\[NegativeVeryThinSpace]",
            Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]0"],
            ",\[ThinSpace]",
            0,
            ")" // textStylePointBracket
          } // textStylePoint
          , {xCritical, 0}
          , {-1.5, -0.2}
        ],
        {}
      },
      {}
      , PlotRange -> All
      , PlotRangePadding -> 0.1
    ];
  Show[plot
    , ImageSize -> 0.5 * 0.8 ImageSizeTextWidth
  ] // Ex["wedge_acute-different-angle-offset-top.pdf"]
]


(* ::Subsection:: *)
(*Side view*)


Module[
  {
    gamma, gammaTracing,
    d, h, hTracing,
    xiEnd, tStart, tEnd,
    xiMax, tMin, tMaxWall, tMaxAxis,
    wallThickness, angleMarkerRadius, angleMarkerRadiusTracing,
    textStyleLabel, textStylePointBracket,
    dummyForTrailingCommas
  },
  (* Constants *)
  gamma = 25 Degree;
  gammaTracing = 60 Degree;
  (* Half-plane offset distance & wall-height *)
  d = DHalfPlane[gamma, gammaTracing];
  h = HHalfPlane[gamma];
  hTracing = HHalfPlane[gammaTracing];
  (* Endpoints for parameter T *)
  xiEnd = 2.8 d;
  tStart = h;
  tEnd = SeekRoot[XHalfPlane[gamma][#] - xiEnd &, {0, tStart}, 10] // Quiet;
  (* Plot range *)
  xiMax = 1.1 xiEnd;
  tMin = 0;
  tMaxWall = 1.2 tStart;
  tMaxAxis = 1.08 tMaxWall;
  (* Plotting constants *)
  wallThickness = 0.03 d;
  angleMarkerRadius = 0.18 h;
  angleMarkerRadiusTracing = 0.13 h;
  (* Text style *)
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStylePointBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Diagram *)
  Show[
    (* Liquid surface *)
    ParametricPlot[
      {XHalfPlane[gamma][t], t}
      , {t, tStart, tEnd}
      , Axes -> None
      , ImageSize -> 0.5 ImageSizeTextWidth
      , PlotPoints -> 2
      , PlotRange -> All
      , PlotRangeClipping -> False
      , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
    ],
    (* Original contact angle *)
    Graphics @ {
      Circle[
        {0, h},
        angleMarkerRadius,
        {-Pi/2, -Pi/2 + gamma}
      ],
      Text[
        "\[Gamma]" // textStyleLabel
        , {0, h} + XYPolar[angleMarkerRadius, -Pi/2 + gamma/2]
        , {-0.6, 0.7}
      ],
      {}
    },
    (* Tracing contact angle *)
    Graphics @ {
      Circle[
        {d, hTracing},
        angleMarkerRadiusTracing,
        {-Pi/2, -Pi/2 + gammaTracing}
      ],
      Text[
        Subscript["\[Gamma]", Style["\[NegativeVeryThinSpace]\[Bullet]", Magnification -> 1.8]]
          // textStyleLabel
        , {d, hTracing} + XYPolar[angleMarkerRadius, -Pi/2 + gammaTracing/2]
        , {-0.6, 0.1}
      ],
      {}
    },
    (* Guiding lines for tracing contact angle *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {
        {-wallThickness/2, hTracing},
        {d, hTracing}
      },
      Line @ {
        {d, 0},
        {d, hTracing}
      },
      {}
    },
    (* Wall *)
    Graphics @ {
      Line @ {{-wallThickness/2, tMin}, {-wallThickness/2, tMaxAxis}},
      Text[
        Italicise["T"] // textStyleLabel
        , {-wallThickness/2, tMaxAxis}
        , {0, -1.2}
      ],
      {}
    },
    Graphics @ {BoundaryTracingStyle["Wall"],
      Rectangle[{-wallThickness, tMin}, {0, tMaxWall}]
    },
    Graphics @ {
      Text[
        Italicise["h"] // textStyleLabel
        , {-wallThickness, h}
        , {2.25, -0.25}
      ],
      Text[
        Subscript[Italicise["h"], Style["\[NegativeVeryThinSpace]\[Bullet]", Magnification -> 1.8]]
          // textStyleLabel
        , {-wallThickness, hTracing}
        , {1.5, 0}
      ],
      Text[
        0 // textStyleLabel
        , {-wallThickness, 0}
        , {2.5, 0}
      ],
      {}
    },
    (* Base *)
    Graphics @ {
      Line @ {{-wallThickness, tMin}, {xiMax, tMin}},
      Text[
        "\[Xi]" // textStyleLabel
        , {xiMax, tMin}
        , {-3.5, 0}
      ],
      {}
    },
    Graphics @ {
      Text[
        0 // textStyleLabel
        , {-wallThickness/2, 0}
        , {0, 1}
      ],
      Text[
        Row @ {
          Italicise["d"],
          "(" // textStylePointBracket,
          "\[NegativeVeryThinSpace]",
          "\[Gamma]",
          ",\[ThinSpace]",
          Subscript["\[Gamma]", Style["\[NegativeVeryThinSpace]\[Bullet]", Magnification -> 1.8]],
          ")" // textStylePointBracket
        } // textStyleLabel
        , {d, 0}
        , {0, 1}
      ],
      {}
    },
    {}
  ]
] // Ex["wedge_acute-different-angle-offset-side.pdf"]


(* ::Section:: *)
(*Figure: family of candidate boundaries (wedge_acute-candidates)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracingValues,
    xCriticalMin, xCriticalMax,
    xMax, yMax, rMax,
    more, xMaxMore,
    xyTracedAss,
    gammaTracing, xCritical,
    derList, p, q, grad2, f, vi,
    sMax,
    xArrowMin, xArrowMax,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 40;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Tracing contact angles *)
  gpdTracingValues = Range[gpd, 85, 5];
  (* Plot range *)
  xCriticalMin = x0[tNumerical, Min[gpdTracingValues] Degree];
  xCriticalMax = x0[tNumerical, Max[gpdTracingValues] Degree];
  xMax = Ceiling[1.5 xCriticalMax];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.1;
  xMaxMore = xMax + more;
  (* Compute traced boundary candidates *)
  (* (association from gpdTracing to traced-boundary *)
  xyTracedAss =
    Association @ Table[
      (* Tracing contact angle *)
      gammaTracing = gpdTracing * Degree;
      (* Critical terminal point x_0 *)
      xCritical = x0[tNumerical, gammaTracing];
      (* Derivative list for boundary tracing *)
      {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
      (* Traced boundary (upper branch) *)
      sMax = 20;
      gpdTracing ->
        ContactTracedBoundary[derList][
          {xCritical, 0}, 0, {0, sMax}
          , -1, 1
          , -Infinity, Function[{x, y}, x > xMaxMore]
        ]
      , {gpdTracing, gpdTracingValues}
    ];
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMaxMore {1, Tan[alpha]},
        {0, 0},
        xMaxMore {1, -Tan[alpha]}
      }
    },
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        EvaluatePair[xy, Abs[s], Sign[s]]
          // Evaluate
        , {s, -DomainEnd[xy], DomainEnd[xy]}
        , PlotPoints -> 3
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedAss}
    ],
    (* Manually extend gammaTracing == gamma case *)
    Table[
      ParametricPlot[
        Way[
          EvaluatePair[xy, DomainEnd[xy]],
          xMaxMore {1, Tan[alpha]}
          , prop
        ]
          // IncludeYReflection
          // Evaluate
        , {prop, 0, 1}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedAss[[{1}]]}
    ],
    (* Tracing gamma arrow *)
    xArrowMin = Way[xCriticalMin, xCriticalMax, -0.1];
    xArrowMax = Way[xCriticalMin, xCriticalMax, +1.25];
    Graphics @ {
      GeneralStyle["Translucent"], GeneralStyle["Thick"],
      Arrowheads[Medium],
      Arrow @ {{xArrowMin, 0}, {xArrowMax, 0}}
    },
    Graphics @ {
      Text[
        Subscript["\[Gamma]", Style["\[NegativeVeryThinSpace]\[Bullet]", Magnification -> 1.8]]
          // LaTeXStyle
          // Style[#, LabelSize["Label"]] &
        , {xArrowMax, 0}
        , {-1.6, -0.2}
      ]
    },
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["wedge_acute-candidates.pdf"]


(* ::Section:: *)
(*Figure: candidate boundaries grouped by tracing contact angle*)
(*(wedge_acute-candidates-by-tracing-angle)*)


Module[
  {
    apd, gpdTracing,
    alpha, gammaTracing,
    gpdValues,
    tNumericalAss,
    xCriticalMin, xCriticalMax,
    xMax, yMax, rMax,
    more, xMaxMore,
    xyTracedAss,
    gamma, tNumerical, xCritical,
    derList, p, q, grad2, f, vi,
    sMax,
    xArrowMin, xArrowMax,
    dummyForTrailingCommas
  },
  (* Prescribed angles *)
  apd = 40;
  gpdTracing = 70;
  alpha = apd * Degree;
  gammaTracing = gpdTracing * Degree;
  (* Known solution angles *)
  gpdValues = Range[90 - apd, gpdTracing, 5];
  (* Import numerical solutions *)
  (* (association from gpd to numerical solution *)
  tNumericalAss =
    Association @ Table[
      gpd -> (
        Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Uncompress // First
      )
      , {gpd, gpdValues}
    ];
  (* Plot range *)
  xCriticalMin = x0[tNumericalAss @ Max[gpdValues], gammaTracing];
  xCriticalMax = x0[tNumericalAss @ Min[gpdValues], gammaTracing];
  xMax = Ceiling[1.5 xCriticalMax, 0.5];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.1;
  xMaxMore = xMax + more;
  (* Compute traced boundary candidates *)
  (* (association from gpd to traced boundary *)
  xyTracedAss =
    Association @ Table[
      (* Solution *)
      gamma = gpd * Degree;
      tNumerical = tNumericalAss[gpd];
      (* Critical terminal point x_0 *)
      xCritical = x0[tNumerical, gammaTracing];
      (* Derivative list for boundary tracing *)
      {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
      (* Traced boundary (upper branch) *)
      sMax = 20;
      gpd ->
        ContactTracedBoundary[derList][
          {xCritical, 0}, 0, {0, sMax}
          , -1, 1
          , -Infinity, Function[{x, y}, x > xMaxMore]
        ]
      , {gpd, gpdValues}
    ];
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMaxMore {1, Tan[alpha]},
        {0, 0},
        xMaxMore {1, -Tan[alpha]}
      }
    },
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        EvaluatePair[xy, Abs[s], Sign[s]]
          // Evaluate
        , {s, -DomainEnd[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedAss}
    ],
    (* Manually extend gammaTracing == gamma case *)
    Table[
      ParametricPlot[
        Way[
          EvaluatePair[xy, DomainEnd[xy]],
          xMaxMore {1, Tan[alpha]}
          , prop
        ]
          // IncludeYReflection
          // Evaluate
        , {prop, 0, 1}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedAss[gpdTracing] // List}
    ],
    (* Solution gamma arrow *)
    xArrowMin = Way[xCriticalMin, xCriticalMax, -0.45];
    xArrowMax = Way[xCriticalMin, xCriticalMax, +1.25];
    Graphics @ {
      GeneralStyle["Translucent"], GeneralStyle["Thick"],
      Arrowheads[Medium],
      Arrow @ {{xArrowMax, 0}, {xArrowMin, 0}}
    },
    Graphics @ {
      Text[
        "\[Gamma]"
          // LaTeXStyle
          // Style[#, LabelSize["Label"]] &
        , {xArrowMax, 0}
        , {-2.5, -0.2}
      ]
    },
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["wedge_acute-candidates-by-tracing-angle.pdf"]


(* ::Section:: *)
(*Figure: candidate boundaries (grouped by tracing contact angle) with offset applied*)
(*(wedge_acute-candidates-offset-alpha-*_degree)*)


(* (This is slow (~10 sec).) *)
(*
  NOTE: a different approach is taken here compared to the above.
  Here we use symbols with assigned values rather than associations.
*)
Module[
  {
    apdValues, gpdTracing, gammaTracing,
    alpha, gpdValues,
    tNumerical,
    sMax,
    xyTraced,
    gamma, xCritical,
    derList, p, q, grad2, f, vi,
    d, xCriticalOffset,
    plotRangeFactor,
    xMax, yMax, xMaxMore,
    graphicsWedgeWalls, xy, plotTracedBoundaries,
    plot,
    xCriticalOffsetMinGamma, xCriticalOffsetMaxGamma,
    arrowDirection,
    xArrowBase, xArrowTip,
    xGammaLabel, gammaLabelXShift,
    xMinInset, xMaxInset, yMinInset, yMaxInset,
    insetRectangleStyle, insetRectangle, insetPlot,
    dummyForTrailingCommas
  },
  (* Prescribed angles *)
  apdValues = {30, 45, 60};
  gpdTracing = 75;
  gammaTracing = gpdTracing * Degree;
  (* For various wedge half-angles: *)
  Table[
    alpha = apd * Degree;
    (* Known solution angles *)
    gpdValues = Range[90 - apd, gpdTracing, 5];
    Which[
      apd == apdValues[[1]],
        gpdValues = Part[gpdValues, All];
      ,
      apd == apdValues[[2]],
        gpdValues = Part[gpdValues, {1, -1}];
      ,
      apd == apdValues[[3]],
        gpdValues = Part[gpdValues, {1, -4, -2, -1}];
    ];
    (* Import numerical solutions *)
    Table[
      tNumerical[gpd] =
        Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Uncompress // First;
      , {gpd, gpdValues}
    ];
    (* Maximum arc-length for boundary tracing *)
    (* (NOT a proven upper bound, but probably enough) *)
    sMax = 4 x0[tNumerical @ Min[gpdValues], gammaTracing];
    (* Compute traced boundary candidates *)
    (* (association from gpd to traced boundary *)
    Table[
      gamma = gpd * Degree;
      (* Critical terminal point x_0 *)
      xCritical[gpd] = x0[tNumerical[gpd], gammaTracing];
      (* Offset critical terminal point x'_0  *)
      d[gpd] = DHalfPlane[gamma, gammaTracing];
      xCriticalOffset[gpd] = xCritical[gpd] - d[gpd] / Sin[alpha];
      (* Derivative list for boundary tracing *)
      {p, q, grad2, f, vi} = derList =
        ContactDerivativeList[tNumerical[gpd], gammaTracing];
      (* Traced boundary (upper branch) *)
      xyTraced[gpd] =
        ContactTracedBoundary[derList][
          {xCritical[gpd], 0}, 0, {0, sMax}
          , -1, 1
          , -Infinity
        ];
      , {gpd, gpdValues}
    ];
    (* Plot range *)
    plotRangeFactor =
      Which[
        apd == apdValues[[1]], 1.8,
        apd == apdValues[[2]], 2,
        apd == apdValues[[3]], 2.2
      ];
    xMax = plotRangeFactor * Max[xCriticalOffset /@ gpdValues];
    yMax = xMax;
    xMaxMore = 1.1 xMax;
    (* Plots re-used for inset *)
    graphicsWedgeWalls =
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          xMaxMore {1, Tan[alpha]},
          {0, 0},
          xMaxMore {1, -Tan[alpha]}
        }
      };
    plotTracedBoundaries =
      Table[
        xy = xyTraced[gpd];
        ParametricPlot[
          Evaluate @ Subtract[
            EvaluatePair[xy, Abs[s], Sign[s]],
            {d[gpd] / Sin[alpha], 0}
          ]
          , {s, -DomainEnd[xy], DomainEnd[xy]}
          , PlotPoints -> If[apd < apdValues[[3]], 5, 4]
          , PlotStyle -> Directive[
              BoundaryTracingStyle["Traced"],
              GeneralStyle["DefaultThick"]
            ]
        ]
        , {gpd, gpdValues}
      ];
    (* Main plot *)
    plot = Show[
      EmptyFrame[{0, xMax}, {-yMax, yMax}
      ],
      (* Wedge walls *)
      graphicsWedgeWalls,
      (* Traced boundaries *)
      plotTracedBoundaries,
      (* Solution gamma arrow *)
      xCriticalOffsetMinGamma = xCriticalOffset[Min @ gpdValues];
      xCriticalOffsetMaxGamma = xCriticalOffset[Max @ gpdValues];
      arrowDirection = Sign[xCriticalOffsetMaxGamma - xCriticalOffsetMinGamma];
      xArrowBase = xCriticalOffsetMinGamma - 0.1 xMax * arrowDirection;
      xArrowTip = xCriticalOffsetMaxGamma + 0.2 xMax * arrowDirection;
      Graphics @ {
        GeneralStyle["Translucent"], GeneralStyle["Thick"],
        Arrowheads[Medium],
        Arrow @ {{xArrowBase, 0}, {xArrowTip, 0}}
      },
      (* Solution gamma label *)
      xGammaLabel = Max[xArrowBase, xArrowTip];
      gammaLabelXShift = If[xArrowBase < xArrowTip, -1.6, -2.5];
      Graphics @ {
        Text[
          "\[Gamma]"
            // LaTeXStyle
            // Style[#, LabelSize["Label"]] &
          , {xGammaLabel, 0}
          , {gammaLabelXShift, -0.2}
        ]
      },
      {}
      , ImageSize -> 0.31 ImageSizeTextWidth
    ];
    (* Inset to show bulging *)
    If[apd == apdValues[[1]],
      (* Inset bounding box coordinates *)
      xMinInset = xCriticalOffsetMinGamma;
      xMaxInset = Way[xMinInset, xMax, 0.525];
      yMinInset = xMinInset Tan[alpha];
      yMaxInset = Way[yMinInset, yMax, 0.3];
      (* Inset rectangle *)
      insetRectangleStyle = Directive[EdgeForm[Black], FaceForm[None]];
      insetRectangle = Rectangle[{xMinInset, yMinInset}, {xMaxInset, yMaxInset}];
      (* Inset plot *)
      insetPlot =
        Show[
          graphicsWedgeWalls,
          plotTracedBoundaries,
          Graphics @ {insetRectangleStyle, insetRectangle},
          {}
          , ImageSize -> 0.15 ImageSizeTextWidth
          , PlotRange -> {{xMinInset, xMaxInset}, {yMinInset, yMaxInset}}
        ];
      (* Put inset plot into main plot *)
      plot =
        Show[
          plot,
          Graphics @ {insetRectangleStyle, insetRectangle},
          {}
          , Epilog -> {
              Inset[insetPlot, Scaled @ {0.4, 0.852}]
            }
        ];
    ];
    (* Export table of offset critical terminal point coordinates *)
    Table[{gpd, xCriticalOffset[gpd]}, {gpd, gpdValues}]
      // TableForm[#,
           TableHeadings -> {
             None,
             {"gpd", "xCriticalOffset"}
           }
         ] &
      // Ex @ FString["wedge_acute-candidates-offset-alpha-{apd}_degree-table.pdf"]
    ;
    (* Export plot *)
    Show[
      plot
      , FrameLabel -> {
          Superscript[
            Italicise["x"],
            Style["\[NegativeVeryThinSpace]\[Prime]", Magnification -> 1.3]
          ] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ] // Ex @ FString["wedge_acute-candidates-offset-alpha-{apd}_degree.pdf"]
    , {apd, apdValues}
  ]
]


(* ::Section:: *)
(*Figure: well-approximating circular arc (wedge_acute-well-approximating-arc)*)


Module[
  {
    alpha,
    wedgeRMax, wedgeXMax, wedgeYMax,
    clearanceProportion, xMin, xMax, yMax,
    xTip, rad, xCentre, angMin, angMax,
    xTangency, yTangency,
    orthogonalityMarkerLength,
    angleMarkerRadius, textSize, textStyle,
    dummyForTrailingCommas
  },
  (* Wedge *)
  alpha = 30 Degree;
  wedgeRMax = 1;
  {wedgeXMax, wedgeYMax} = XYPolar[wedgeRMax, alpha];
  (* Plot range *)
  clearanceProportion = 1/10;
  xMin = Way[0, wedgeXMax, -clearanceProportion];
  xMax = Way[0, wedgeXMax, 1 + clearanceProportion];
  yMax = Way[0, wedgeYMax, 1 + clearanceProportion];
  (* Tip of circle (offset critical terminal point x'_0) *)
  xTip = Way[xMin, xMax, 0.45];
  (* Radius, centre, and angular range of circular arc *)
  rad = xTip Sin[alpha] / (1 - Sin[alpha]);
  xCentre = rad / Sin[alpha];
  angMin = Pi/2 + alpha;
  angMax = 2 Pi - angMin;
  (* Point of tangency unto upper wedge wall *)
  {xTangency, yTangency} = XYPolar[rad / Tan[alpha], alpha];
  (* Etc. *)
  orthogonalityMarkerLength = 0.125 rad;
  angleMarkerRadius = 0.3 rad;
  textSize = LabelSize["Label"];
  textStyle = Style[#, textSize] & @* LaTeXStyle;
  (* Diagram *)
  Show[
    EmptyAxes[{xMin, xMax}, {-yMax, yMax}
      , AspectRatio -> Automatic
      , AxesLabel -> {
          Superscript[
            Italicise["x"],
            Style["\[Prime]", Magnification -> 1.3]
          ] // Margined @ {{0, 0}, {0, -5}},
          Italicise["y"]
        }
      , LabelStyle -> LatinModernLabelStyle @ textSize
      , Method -> {"AxesInFront" -> False}
      , Ticks -> None
    ],
    (* Radius unto tangent at upper wedge wall *)
    Graphics @ {
      Line @ {{xCentre, 0}, {xTangency, yTangency}}
    },
    Graphics @ {
      Text[
        Italicise["R"] // textStyle
        , Way[{xCentre, 0}, {xTangency, yTangency}]
        , {-1.7, -0.5}
      ]
    },
    Graphics @ {
      Line[orthogonalityMarkerLength {{1, 0}, {1, -1}, {0, -1}}]
        // Rotate[#, alpha, {0, 0}] &
        // Translate[#, {xTangency, yTangency}] &
    },
    (* Wedge half-angle *)
    Graphics @ {
      Circle[{0, 0}, angleMarkerRadius, {0, alpha}]
    },
    Graphics @ {
      Text[
        "\[Alpha]" // textStyle
        , XYPolar[angleMarkerRadius, alpha / 2]
        , {-2.5, -0.3}
      ]
    },
    (* Wedge walls *)
    Graphics @ {
      BoundaryTracingStyle["Wall"],
      Line @ {{wedgeXMax, wedgeYMax}, {0, 0}, {wedgeXMax, -wedgeYMax}}
    },
    (* Circular arc *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Circle[{xCentre, 0}, rad, {angMin, angMax}]
    },
    (* Tip of circle (offset critical terminal point x'_0) *)
    Graphics @ {GeneralStyle["Point"],
      Point @ {xTip, 0}
    },
    Graphics @ {
      Text[
        Subsuperscript[
          Italicise["x"],
          0,
          Style["\[Prime]", Magnification -> 1.3]
        ] // textStyle
        , {xTip, 0}
        , {1.25, 1.1}
      ]
    },
    (* Centre of circle *)
    Graphics @ {GeneralStyle["Point"],
      Point @ {xCentre, 0}
    },
    Graphics @ {
      Text[
        SeparatedRow["/"][
          Italicise["R"],
          SeparatedRow["Thin"]["sin", "\[Alpha]"]
        ] // textStyle
        , {xCentre, 0}
        , {0, 1.1}
      ]
    },
    {}
    , ImageSize -> 0.45 ImageSizeTextWidth
  ]
] // Ex["wedge_acute-well-approximating-arc.pdf"]


(* ::Section:: *)
(*Table: well-approximating radius versus height rise*)


(* (This is slow (~6 sec).) *)
Module[
  {
    apd, gpdTracing,
    alpha, gammaTracing,
    gpdValues,
    tNumerical,
    gamma, xCritical,
    d, xCriticalOffset,
    rad, height,
    dummyForTrailingCommas
  },
  (* Prescribed angles *)
  apd = 60;
  gpdTracing = 75;
  alpha = apd * Degree;
  gammaTracing = gpdTracing * Degree;
  (* For various known solution angles: *)
  gpdValues = Range[90 - apd, gpdTracing, 5];
  Table[
    gamma = gpd * Degree;
    (* Import numerical solutions *)
    tNumerical[gpd] =
      Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    (* Compute critical terminal point x_0 *)
    xCritical[gpd] = x0[tNumerical[gpd], gammaTracing];
    (* Offset critical terminal point x'_0  *)
    d[gpd] = DHalfPlane[gamma, gammaTracing];
    xCriticalOffset[gpd] = xCritical[gpd] - d[gpd] / Sin[alpha];
    (* Well-approximating rounding radius *)
    rad[gpd] = xCriticalOffset[gpd] * Sin[alpha] / (1 - Sin[alpha]);
    (* Height rise *)
    height[gpd] = tNumerical[gpd][xCritical[gpd], 0];
    , {gpd, gpdValues}
  ];
  (* Export table *)
  Table[
    {gpd, rad[gpd], height[gpd]}
    , {gpd, gpdValues}
  ] // TableForm[#,
    TableHeadings -> {
      None,
      {"gpd", "radius", "height"}
    }
  ] &
] // Ex["wedge_acute-well-approximating-table.pdf"]


(* ::Section:: *)
(*Figure: height-rise profiles (wedge_acute-height-rise-profiles)*)


(* (This is slow (~5 sec).) *)
Module[
  {
    apd, gpdTracing,
    alpha, gammaTracing, h,
    gpdValues,
    sMax,
    xySharp,
    tNumerical, xCritical,
    derList, p, q, grad2, f, vi,
    xyTraced,
    dummyForTrailingCommas
  },
  (* Prescribed angles *)
  apd = 60;
  gpdTracing = 75;
  alpha = apd * Degree;
  gammaTracing = gpdTracing * Degree;
  h = HHalfPlane[gammaTracing];
  (* Known solution gamma *)
  gpdValues = {90 - apd, gpdTracing};
  (* Plot range *)
  sMax = 4;
  (* Sharp corner (original walls) *)
  xySharp[s_] := XYPolar[Abs[s], Sign[s] alpha];
  (* For the chosen known solution gamma: *)
  Table[
    (* Import numerical solution *)
    tNumerical[gpd] =
      Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    (* Critical terminal point x_0 *)
    xCritical[gpd] = x0[tNumerical[gpd], gammaTracing];
    (* Derivative list for boundary tracing *)
    {p, q, grad2, f, vi} = derList =
      ContactDerivativeList[tNumerical[gpd], gammaTracing];
    (* Traced boundary (upper branch) *)
    xyTraced[gpd] =
      ContactTracedBoundary[derList][
        {xCritical[gpd], 0}, 0, {0, sMax}
        , -1, 1
        , -Infinity
      ];
    , {gpd, gpdValues}
  ];
  (* Plot *)
  Plot[
    Join[
      (* Half-plane wall height *)
      h // List,
      (* Sharp corner *)
      tNumerical[gpdTracing] @@ xySharp[s] // List,
      (* Rounded corners *)
      Table[
        tNumerical[gpd] @@ EvaluatePair[xyTraced[gpd], Abs[s], Sign[s]]
        , {gpd, gpdValues}
      ],
      {}
    ] // Evaluate
    , {s, -sMax, sMax}
    , AspectRatio -> 0.5
    , AxesLabel -> Italicise /@ {"s", "T"}
    , ImageSize -> 0.6 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotPoints -> 4
    , PlotRange -> {0, All}
    , PlotStyle ->
        Join[
          BoundaryTracingStyle["Contour"] // List,
          BoundaryTracingStyle["Wall"] // List,
          ConstantArray[
            Directive[BoundaryTracingStyle["Traced"], GeneralStyle["DefaultThick"]],
            Length[gpdValues]
          ],
          {}
        ]
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes] // Evaluate
  ]
] // Ex["wedge_acute-height-rise-profiles.pdf"]


(* ::Section:: *)
(*Figure: radius vs height rise comparison (wedge_acute-radius-height-rise-comparison)*)


Module[
  {
    boundaryTracingTable, conventionalTable,
    heightMin, heightMax,
    dummyForTrailingCommas
  },
  (* Import radius vs height rise tables *)
  conventionalTable = Import["wedge_acute-radius-height-rise-conventional.csv"];
  boundaryTracingTable = Import["wedge_acute-radius-height-rise-tracing.csv"] // Reverse;
  (* Nicer vertical plot range *)
  {heightMin, heightMax} = MinMax @ conventionalTable[[All, 2]];
  heightMin = Floor[heightMin, 0.05];
  heightMax = Ceiling[heightMax, 0.05];
  (* Plot *)
  ListPlot[
    {conventionalTable, boundaryTracingTable}
    , AspectRatio -> 0.75
    , AxesLabel -> {"Radius", "Height"}
    , ImageSize -> 0.7 ImageSizeTextWidth
    , Joined -> True
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotLegends -> Placed[
        PointLegend[
          {"Conventional", "Boundary tracing"}
          , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
          , LegendLayout -> "ReversedColumn"
        ]
        , Scaled @ {0.7, 0.7}
      ]
    , PlotMarkers -> {{Automatic, 6}, {"OpenMarkers", 6}}
    , PlotRange -> {heightMin, heightMax}
    , PlotStyle -> Black
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes] // Evaluate
  ]
] // Ex["wedge_acute-radius-height-rise-comparison.pdf"]


(* ::Section:: *)
(*Figure: modified boundary contours*)


Module[
  {
    alpha, gamma,
    xCentreValues, hValues, xyContourList,
    numContours,
    tNumerical, vi,
    meshes, rUpperLowerJoinPairs,
    xMax, yMax,
    xMaxWall, yMaxWall, rMaxWall,
    style,
      xy,
      rUpperJoin, rLowerJoin,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {alpha, gamma} = {apdMod, gpdMod} Degree;
  (* Import contours *)
  {xCentreValues, hValues, xyContourList} =
    Import["modification/wedge_acute-modification-contours.txt"] // Uncompress;
  numContours = Length[xyContourList] - 1;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
      // Uncompress // First;
  vi = Last @ ContactDerivativeList[tNumerical, gamma];
  (* Import meshes *)
  rUpperLowerJoinPairs =
    Import["modification/wedge_acute-modification-meshes.txt"]
      // Uncompress // #[[All, 3 ;; 4]] &;
  (* Make plot *)
  xMax = 1.2 Max[xCentreValues];
  yMax = xMax Tan[alpha];
  {xMaxWall, yMaxWall} = 1.1 {xMax, yMax};
  rMaxWall = RPolar[xMaxWall, yMaxWall];
  style[n_] := If[
    n > numContoursViable,
    BoundaryTracingStyle["ContourPlain"],
    Black
  ];
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        {xMaxWall, +yMaxWall},
        {0, 0},
        {xMaxWall, -yMaxWall}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, 0, xMaxWall}
      , {y, -yMaxWall, yMaxWall}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Contours *)
    Table[
      {
        xy = xyContourList[[n]];
        ParametricPlot[
          EvaluatePair[xy, s] // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> style[n]
        ],
        {rUpperJoin, rLowerJoin} = rUpperLowerJoinPairs[[n]];
        Graphics @ {
          style[n], GeneralStyle["DefaultThick"],
          Line @ {XYPolar[rUpperJoin, +alpha], XYPolar[rMaxWall, +alpha]},
          Line @ {XYPolar[rLowerJoin, -alpha], XYPolar[rMaxWall, -alpha]},
          {}
        },
        {}
      }
      , {n, numContours}
    ],
    {}
    , FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
        Italicise["y"]
      }
    , ImageSize -> 0.45 * 0.7 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , FrameTicksStyle -> LabelSize["Tick"]
  ]
] // Ex["wedge_acute-modification-contours.pdf"]
