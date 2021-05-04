(* ::Package:: *)

(*
  Basically we extract the singularity
  and treat polar coordinates as "formally rectangular".
  See c1 (manuscripts/capillary-1-small.pdf).
*)


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
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "wedge_small"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Finite element mesh*)


(* Wedge half-angles (apd means "alpha per degree") *)
apdMeshValues = {15, 30, 45, 60};
rMaxMesh = 10;


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString["mesh/wedge_small-mesh-apd-{apd}.txt"],
    Module[
      {
        rMax, alpha,
        fineLengthScale, coarseLengthScale,
        boundaryPoints, numBoundaryPoints, boundaryElements,
        boundaryMesh, mesh,
        predicateUpperWall, predicateLowerWall, predicateCorner,
        dummyForTrailingCommas
      },
      (* Geometry *)
      rMax = rMaxMesh;
      alpha = apd * Degree;
      (* Refinement length scales *)
      (* (Note that we are using the same refinement for phi as for r) *)
      fineLengthScale = 0.01;
      coarseLengthScale = Min[0.1, alpha/5];
      (* Boundary points *)
      boundaryPoints =
        Join[
          (* Corner segment (r == 0) *)
          Table[{0, phi}, {phi, UniformRange[alpha, -alpha, -fineLengthScale]}] // Most,
          (* Lower wedge wall (\[Phi] == -\[Alpha]) *)
          Table[{r, -alpha}, {r, UniformRange[0, rMaxMesh, fineLengthScale]}] // Most,
          (* Far arc at infinity (r == INF) *)
          Table[{rMaxMesh, phi}, {phi, UniformRange[-alpha, alpha, coarseLengthScale]}] // Most,
          (* Upper wedge wall (\[Phi] == +\[Alpha]) *)
          Table[{r, +alpha}, {r, UniformRange[rMaxMesh, 0, -fineLengthScale]}] // Most,
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
      mesh = ToElementMesh[boundaryMesh
        , "ImproveBoundaryPosition" -> True
        , MeshRefinementFunction -> MeshRefinementUniform[coarseLengthScale]
      ];
      (* Predicate for wetting boundaries (wedge walls) *)
      predicateUpperWall = Function[{r, phi}, phi == +alpha // Evaluate];
      predicateLowerWall = Function[{r, phi}, phi == -alpha // Evaluate];
      (* Predicate for corner segment (r == 0) *)
      predicateCorner = Function[{r, phi}, r == 0 // Evaluate];
      (* Return *)
      {mesh, predicateUpperWall, predicateLowerWall, predicateCorner} // Compress
    ]
  ]
  , {apd, apdMeshValues}
]


(* ::Subsection:: *)
(*Solve capillary boundary value problem (fixed-point iteration)*)


(* See Page c1-4 of manuscripts/capillary-1-small.pdf *)


(* (This is slow (~3.5 hr), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
$HistoryLength = 0;
With[{r = \[FormalR], phi = \[FormalPhi], h = \[FormalCapitalH]},
Module[
  {
    gpdValues,
    mesh, predicateUpperWall, predicateLowerWall, predicateCorner,
    meshCoordinates,
    alpha, gamma, k,
    grad, inactiveGrad, inactiveDiv,
    hCornerSegment, hPrevious,
      nMax, absoluteChangeTolerance,
      pTilde, qTilde, c,
      kMatrix, vVector,
      hCurrent,
      absoluteChange,
      returnList,
    dummyForTrailingCommas
  },
  (* For each value of alpha: *)
  Table[
    (* Solution contact angles (gpd means "gamma per degree") *)
    gpdValues = Join[{1}, Range[5, 90 - apd, 5]];
    (* If all values of contact angle have been account for: *)
    If[
      Table[
        FString @ "solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"
        , {gpd, gpdValues}
      ] // AllTrue[FileExistsQ]
      ,
      (* Then skip; *)
      Null
      ,
      (* Otherwise import mesh *)
      {mesh, predicateUpperWall, predicateLowerWall, predicateCorner} =
        FString["mesh/wedge_small-mesh-apd-{apd}.txt"] // Import // Uncompress;
      meshCoordinates = mesh["Coordinates"];
      (* For each value of gamma: *)
      Table[
        (* Angular parameters etc. *)
        alpha = apd * Degree;
        gamma = gpd * Degree;
        k = Sin[alpha] / Cos[gamma];
        (* Solve and export *)
        ExportIfNotExists[FString @ "solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt",
          (* Formally rectangular derivative operators (see (c1.27)) *)
          grad = Grad[#, {r, phi}] &;
          inactiveGrad = Inactive[Grad][#, {r, phi}] &;
          inactiveDiv = Inactive[Div][#, {r, phi}] &;
          (* Boundary data along corner segment (c1.32) *)
          hCornerSegment = Function[{r, phi}, (Cos[phi] - Sqrt[k^2 - Sin[phi]^2]) / k];
          (* May as well use that as the initial guess *)
          hPrevious = hCornerSegment;
          (* Fixed-point iteration *)
          nMax = 1000;
          absoluteChangeTolerance = 10^-4;
          Do[
            (* Components of gradient de-singularised (see (c1.23) and (c1.24)) *)
            pTilde = r * D[hPrevious[r, phi], r] - hPrevious[r, phi];
            qTilde = D[hPrevious[r, phi], phi];
            (* Common scalar coefficient (c1.26) *)
            c = 1 / Sqrt[r^4 + pTilde^2 + qTilde^2];
            (* Diffusion coefficient matrix (c1.28) *)
            kMatrix = {{r^2 * c, 0}, {0, c}};
            (* Conservative convection coefficients (c1.29) *)
            vVector = {r * c, 0};
            (* Solve linearised boundary value problem in H *)
            hCurrent =
              NDSolveValue[
                {
                  (* Left hand side of (c1.30), negated *)
                  inactiveDiv[
                    -kMatrix . inactiveGrad[h[r, phi]]
                    + Inactive[Times][vVector, h[r, phi]]
                  ] ==
                    (* Right hand side of (c1.30), negated *)
                    - h[r, phi]
                    (* Neumann boundary condition (c1.31) *)
                    + NeumannValue[Cos[gamma], predicateUpperWall[r, phi]]
                    + NeumannValue[Cos[gamma], predicateLowerWall[r, phi]]
                  ,
                  (* Dirichlet condition along corner segment (c1.32) *)
                  DirichletCondition[
                    h[r, phi] == hCornerSegment[r, phi],
                    predicateCorner[r, phi]
                  ]
                },
                h,
                Element[{r, phi}, mesh]
              ];
            (* Compute absolute change *)
            absoluteChange =
              Max @ Table[
                hCurrent @@ coord - hPrevious @@ coord // Abs
                , {coord, meshCoordinates}
              ];
            (* Store return list *)
            returnList = {hCurrent, n, absoluteChange};
            (* If global convergence reached: *)
            If[absoluteChange < absoluteChangeTolerance,
              (* Stop iteration *)
              Break[];
              ,
              (* Otherwise update solution *)
              hPrevious = hCurrent;
            ];
            , {n, nMax}
          ];
          (* Return *)
          returnList // Compress
        ]
        , {gpd, gpdValues}
      ]
    ]
    , {apd, apdMeshValues}
  ]
]
]
$HistoryLength = Infinity;


(* ::Subsection:: *)
(*Boundary tracing*)


(* ::Subsubsection:: *)
(*List of pure-function derivatives*)


(*
  See (c1.23), (c1.24), (c1.37), and (c1.40).
  Analogue to LaplaceYoung`ContactDerivativeList.
*)
tracingDerivativeList[hNumerical_, gammaTracing_?NumericQ] :=
  Module[
    {
      hRDerivative, hPhiDerivative,
      pTilde, qTilde, fTilde, phiTilde,
      dummyForTrailingCommas
    },
    (* Straight derivatives *)
    hRDerivative = Derivative[1, 0][hNumerical];
    hPhiDerivative = Derivative[0, 1][hNumerical];
    (* De-singularised components of gradient (P-tilde, Q-tilde) *)
    pTilde = Function[{r, phi},
      r * hRDerivative[r, phi] - hNumerical[r, phi]
        // Evaluate
    ];
    qTilde = hPhiDerivative;
    (* De-singularised flux function (F-tilde) *)
    fTilde = Function[{r, phi},
      Cos[gammaTracing]
      Sqrt[r^4 + pTilde[r, phi]^2 + qTilde[r, phi]^2]
        // Evaluate
    ];
    (* De-singularised viability function (\[CapitalPhi]-tilde) *)
    phiTilde = Function[{r, phi},
      Sin[gammaTracing]^2 (pTilde[r, phi]^2 + qTilde[r, phi]^2)
        - r^4 * Cos[gammaTracing]^2
        // Evaluate
    ];
    (* Return *)
    {pTilde, qTilde, fTilde, phiTilde}
  ];


(* ::Subsubsection:: *)
(*Elliptic critical terminal point*)


r0[pTilde_, gammaTracing_] :=
  Quiet[
    SeekRoot[
      pTilde[#, 0] / #^2 + Cot[gammaTracing] &,
      {0, rMaxMesh}
    ]
    , {Power::infy}
  ];


(* ::Subsubsection:: *)
(*Traced boundary (arc-length parametrisation)*)


(* See (c1.42) & (c1.43). *)
tracedBoundary[derList_][{r0_, phi0_}, s0_, {sStart_, sEnd_}
  , sSign_: 1
  , branchSign_: 1
  , terminationPhiTilde_: 0
  , terminationFunction_: (False &)
] :=
  Module[
    {
      pTilde, qTilde, fTilde, phiTilde,
      rDerivative, phiDerivative,
      dummyForTrailingCommas
    },
    (* Derivative list *)
    {pTilde, qTilde, fTilde, phiTilde} = derList;
    (* Right hand sides of boundary tracing system of ODEs *)
    rDerivative[r_, phi_] :=
      Divide[
        -qTilde[r, phi] fTilde[r, phi] + branchSign * pTilde[r, phi] Re @ Sqrt @ phiTilde[r, phi],
        pTilde[r, phi]^2 + qTilde[r, phi]^2
      ];
    phiDerivative[r_, phi_] :=
      Divide[
        +pTilde[r, phi] fTilde[r, phi] + branchSign * qTilde[r, phi] Re @ Sqrt @ phiTilde[r, phi],
        pTilde[r, phi]^2 + qTilde[r, phi]^2
      ] / r;
    (* Solve boundary tracing system of ODEs *)
    (* NOTE:
      Using `phi = \[FormalPhi]` results in a Transpose::nmtx error.
      See <https://mathematica.stackexchange.com/q/236730>.
    *)
    With[{r = \[FormalR], phi = \[FormalF], s = \[FormalS]},
      NDSolveValue[
        {
          r'[s] == sSign * rDerivative[r[s], phi[s]],
          phi'[s] == sSign * phiDerivative[r[s], phi[s]],
          r[s0] == r0,
          phi[s0] == phi0,
          WhenEvent[
            {
              phiTilde[r[s], phi[s]] < terminationPhiTilde,
              terminationFunction[r[s], phi[s]]
            },
            "StopIntegration"
          ]
        }
        , {r, phi}
        , {s, sStart, sEnd}
        , NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*Contour (arc-length parametrisation)*)


contour[derList_][{r0_, phi0_}, s0_, {sStart_, sEnd_}
  , sSign_: 1
  , terminationFunction_: (False &)
] :=
  Module[{pTilde, qTilde, slopeTilde, rDerivative, phiDerivative},
    (* Derivatives *)
    {pTilde, qTilde} = derList[[1 ;; 2]];
    slopeTilde = Function[{r, phi}, Sqrt[pTilde[r, phi]^2 + qTilde[r, phi]^2] // Evaluate];
    (* Right hand sides of contour system of ODEs *)
    rDerivative[r_, phi_] := +qTilde[r, phi] / slopeTilde[r, phi];
    phiDerivative[r_, phi_] := -pTilde[r, phi] / slopeTilde[r, phi] / r;
    (* Solve contour system of ODEs *)
    With[{r = \[FormalR], phi = \[FormalF], s = \[FormalS]},
      NDSolveValue[
        {
          r'[s] == sSign * rDerivative[r[s], phi[s]],
          phi'[s] == sSign * phiDerivative[r[s], phi[s]],
          r[s0] == r0,
          phi[s0] == phi0,
          WhenEvent[
            terminationFunction[r[s], phi[s]],
            "StopIntegration"
          ]
        }
        , {r, phi}
        , {s, sStart, sEnd}
        , NoExtrapolation
      ]
    ]
  ]


(* ::Subsection:: *)
(*Modifications (viable contours)*)


(* ::Subsubsection:: *)
(*Angular parameters*)


{apdMod, gpdMod} = {45, 30};
numContoursViable = 4;


(* ::Subsubsection:: *)
(*Contours*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["modification/wedge_small-modification-contours.txt",
  Module[
    {
      gamma, hNumerical, derList, tNumerical,
      rCritical,
      hWall, rCentreWall,
      numberUpToCritical, rStep, numberUpToWall,
      rCentreValues, heightValues,
      sMax, rphiContourList,
      dummyForTrailingCommas
    },
    (* Numerical solution *)
    gamma = gpdMod * Degree;
    hNumerical =
      Import @ FString["solution/wedge_small-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
        // Uncompress // First;
    derList = tracingDerivativeList[hNumerical, gamma];
    tNumerical[r_, phi_] := hNumerical[r, phi] / r // Evaluate;
    (* Critical terminal point r_0 *)
    rCritical = r0[derList // First, gamma];
    (* Point along centreline which has wall height *)
    hWall = HHalfPlane[gamma];
    rCentreWall =
      SeekRoot[
        tNumerical[#, 0] - hWall &,
        {0, 10}
        , 5
      ];
    (* Centreline r-values and height rises for contours *)
    numberUpToCritical = numContoursViable;
    rStep = rCritical / numberUpToCritical;
    numberUpToWall = Floor[rCentreWall / rStep];
    rCentreValues = rCritical * Range[numberUpToWall] / numberUpToCritical;
    rCentreValues = Take[rCentreValues, UpTo[7]];
    heightValues = Table[tNumerical[r, 0], {r, rCentreValues}];
    (* Compute T-contours *)
    sMax = 4; (* probably enough *)
    rphiContourList =
      Table[
        Quiet[
          contour[derList][{r, 0}, 0, sMax {-1, 1}]
          , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
        ]
        , {r, rCentreValues}
      ];
    (* Return lists to export *)
    {rCentreValues, heightValues, rphiContourList} // Compress
  ]
]


(* ::Subsubsection:: *)
(*Finite element meshes*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["modification/wedge_small-modification-meshes.txt",
  Module[
    {
      rMax, alpha,
      rCentreValues, heightValues, rphiContourList,
      fineLengthScale,
      coarseLengthScale,
        sContourUpper, sContourLower, sContourValues,
        xyContourPoints, rphiContourPoints,
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
    {rCentreValues, heightValues, rphiContourList} =
      Import["modification/wedge_small-modification-contours.txt"] // Uncompress;
    (* Mesh length scales *)
    fineLengthScale = 0.01;
    coarseLengthScale = 0.5;
    (* For each contour *)
    Table[
      (* Raw T-contour points *)
      sContourUpper = DomainEnd[rphi];
      sContourLower = DomainStart[rphi];
      sContourValues = Join[
        UniformRange[sContourUpper, 0, -fineLengthScale],
        UniformRange[0, sContourLower, -fineLengthScale] // Rest,
        {}
      ];
      rphiContourPoints = Table[EvaluatePair[rphi, s], {s, sContourValues}];
      xyContourPoints = Table[XYPolar @@ pair, {pair, rphiContourPoints}];
      (* Interpolate linearly in polar coordinates unto upper wall *)
      {rUpper1, phiUpper1} = rphiContourPoints[[1]];
      {rUpper2, phiUpper2} = rphiContourPoints[[2]];
      rUpperJoin = Rescale[alpha, {phiUpper2, phiUpper1}, {rUpper2, rUpper1}];
      (* Interpolate linearly in polar coordinates unto lower wall *)
      {rLower1, phiLower1} = rphiContourPoints[[-1]];
      {rLower2, phiLower2} = rphiContourPoints[[-2]];
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
      , {rphi, rphiContourList}
    ] // Compress
  ]
]


(* ::Subsubsection:: *)
(*Solve capillary BVP and extract height rise profiles*)


(* (This is slow (~40 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["modification/wedge_small-modification-profiles.txt",
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
    meshStuffList = Import["modification/wedge_small-modification-meshes.txt"] // Uncompress;
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
] // AbsoluteTiming


(* ::Section:: *)
(*Finite element mesh check*)


Table[
  Module[
    {
      mesh,
        predicateUpperWall, predicateLowerWall, predicateCorner,
      numElements,
      boundaryCoordinates,
        upperWallBoundaryCoordinates,
        lowerWallBoundaryCoordinates,
        cornerBoundaryCoordinates,
      plot,
      dummyForTrailingCommas
    },
    (* Import mesh *)
    {mesh, predicateUpperWall, predicateLowerWall, predicateCorner} =
      FString["mesh/wedge_small-mesh-apd-{apd}.txt"] // Import // Uncompress;
    numElements = Length @ mesh[[2, 1, 1]];
    (* Determine predicate-satisfying boundary coordinates *)
    boundaryCoordinates =
      Part[
        mesh["Coordinates"],
        mesh["BoundaryElements"][[1, 1]] // Flatten // Union
      ];
    upperWallBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateUpperWall @@ # &
      ];
    lowerWallBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateLowerWall @@ # &
      ];
    cornerBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateCorner @@ # &
      ];
    (* Make plot *)
    plot =
      Show[
        mesh["Wireframe"],
        ListPlot[
          {
            upperWallBoundaryCoordinates,
            lowerWallBoundaryCoordinates,
            cornerBoundaryCoordinates
          }
          , PlotMarkers -> {Automatic, Small}
        ],
        {}
        , ImageSize -> Large
        , PlotLabel -> FString @ "{numElements} mesh elements"
      ];
    (* Export *)
    {
      (* Full plot *)
      plot
        // Ex @ FString @ "mesh/wedge_small-mesh-apd-{apd}.png"
      ,
      (* Zoom near origin to verify wetting predicate *)
      Show[plot
        , PlotLabel -> None
        , PlotRange -> {{0, 1}, {1/2, 1}} apd * Degree
      ]
        // Ex @ FString @ "mesh/wedge_small-mesh-apd-{apd}-zoom.png"
    }
  ]
  , {apd, apdMeshValues}
]


(* ::Section:: *)
(*Numerical solution check*)


(* ::Subsection:: *)
(*Table*)


(* (This is slow (~15 sec) from all the calls to Import.) *)
Module[
  {
    gpdValues,
    numIterates, absoluteChange,
    dummyForTrailingCommas
  },
  Join @@ Table[
    gpdValues = Join[{1}, Range[5, 90 - apd, 5]];
    Table[
      {numIterates, absoluteChange} =
        FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // Rest;
      {apd * Degree, gpd * Degree, numIterates, absoluteChange}
      , {gpd, gpdValues}
    ]
    , {apd, apdMeshValues}
  ] // TableForm[#
    , TableHeadings -> {None,
        {"alpha", "gamma", "Iterates", "Absolute change"}
      }
  ] &
] // Ex["wedge_small-solution-table.pdf"]


(* ::Subsection:: *)
(*H and T plots*)


(* (This is slow (~30 sec) from all the calls to Import.) *)
Module[
  {
    alpha, gpdValues,
    gamma,
    hSmall, tSmallPolar,
    rMax, phiValues, opts,
    hPlot, tPlot,
    dummyForTrailingCommas
  },
  (* For various alpha *)
  Table[
    alpha = apd * Degree;
    (* For various gamma *)
    gpdValues = Join[{1}, Range[5, 90 - apd, 5]];
    Table[
      gamma = gpd * Degree;
      (* Import wedge_small solution H == H (r, \[Phi]) *)
      hSmall =
        FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // First;
      tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
      (* Plotting constants *)
      rMax = rMaxMesh;
      phiValues = Subdivide[0, alpha, 4];
      opts = {
        ImageSize -> 240,
        PlotLabel -> {"apd", apd, "gpd", gpd},
        PlotOptions[Axes]
      };
      (* Make plots *)
      hPlot =
        Plot[
          Table[hSmall[r, phi], {phi, phiValues}] // Evaluate
          , {r, 0, rMax}
          , AxesLabel -> Italicise /@ {"r", "H"}
          , PlotLegends -> LineLegend[phiValues, LegendLabel -> "\[Phi]"]
          , opts // Evaluate
        ];
      tPlot =
        Plot[
          Table[tSmallPolar[r, phi], {phi, phiValues}]
            // Append @ HHalfPlane[gamma]
            // Evaluate
          , {r, 0, rMax}
          , AxesLabel -> Italicise /@ {"r", "T"}
          , PlotStyle -> (
              ConstantArray[Automatic, Length[phiValues]]
                // Append[Dashed]
            )
          , opts // Evaluate
        ];
      (* Export *)
      Row @ {hPlot, tPlot}
        // Ex @ FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.png"]
      , {gpd, gpdValues}
    ]
    , {apd, apdMeshValues}
  ]
]


(* ::Subsection:: *)
(*Comparison with Cartesian solver, borderline case*)


(* (This is slow (~5 sec) from all the calls to Import.) *)
Module[
  {
    alpha, gpdValues,
    gamma,
    hSmall, tSmallPolar,
    tAcute, tAcutePolar,
    rMax, phiValues, opts,
    plot,
    dummyForTrailingCommas
  },
  (* For various alpha *)
  Table[
    alpha = apd * Degree;
    (* For various gamma *)
    gpdValues = {90 - apd};
    Table[
      gamma = gpd * Degree;
      (* Import wedge_small solution H == H (r, \[Phi]) *)
      hSmall =
        FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // First;
      tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
      (* Import wedge_acute solution T == T (x, y) *)
      tAcute =
        FString["../wedge_acute/solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // First;
      tAcutePolar[r_, phi_] := tAcute @@ XYPolar[r, phi];
      (* Plotting constants *)
      rMax = rMaxMesh;
      phiValues = Subdivide[0, alpha, 4];
      opts = {
        ImageSize -> 270,
        PlotLabel -> {"apd", apd, "gpd", gpd},
        PlotOptions[Axes]
      };
      (* Make plot and export *)
      plot =
        LogPlot[
          Table[
            tSmallPolar[r, phi] / tAcutePolar[r, phi] - 1 // Abs
            , {phi, phiValues}
          ] // Evaluate
          , {r, 0, rMax}
          , AxesLabel -> {Italicise["r"], "Rel. discr."}
          , AxesOrigin -> {0, 0}
          , PlotLegends -> LineLegend[phiValues, LegendLabel -> "\[Phi]"]
          , PlotRange -> Full
          , opts // Evaluate
        ];
      plot
        // Ex @ FString["solution/wedge_small-solution-comparison-apd-{apd}-gpd-{gpd}.png"]
      , {gpd, gpdValues}
    ]
    , {apd, apdMeshValues}
  ]
]


(* ::Section:: *)
(*Modifications (viable contours)*)


(* ::Subsection:: *)
(*Visualise contours*)


Module[
  {
    alpha, gamma,
    rCentreValues, heightValues, rphiContourList,
    hNumerical, derList, tNumerical, phiTilde,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {alpha, gamma} = {apdMod, gpdMod} Degree;
  (* Import contours *)
  {rCentreValues, heightValues, rphiContourList} =
    Import["modification/wedge_small-modification-contours.txt"] // Uncompress;
  (* Import numerical solution *)
  hNumerical =
    Import @ FString["solution/wedge_small-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
      // Uncompress // First;
  derList = tracingDerivativeList[hNumerical, gamma];
  tNumerical[r_, phi_] := hNumerical[r, phi] / r // Evaluate;
  phiTilde = Last @ tracingDerivativeList[hNumerical, gamma];
  (* Make plot *)
  Show[
    (* Wedge walls *)
    Graphics @ {Purple,
      HalfLine[{0, 0}, XYPolar[1, alpha]],
      HalfLine[{0, 0}, XYPolar[1, -alpha]]
    },
    (* Non-viable domain *)
    RegionPlot[
      phiTilde @@ RPhiPolar[x, y] < 0
      , {x, 0, 1.5}
      , {y, -1.5, 1.5}
      , BoundaryStyle -> Pink
      , PlotStyle -> Directive[Opacity[0.7], LightGray]
    ],
    (* Contours *)
    Table[
      ParametricPlot[
        XYPolar @@ EvaluatePair[rphi, s] // Evaluate
        , {s, DomainStart[rphi], DomainEnd[rphi]}
        , PlotStyle -> Black
      ]
      , {rphi, rphiContourList}
    ],
    (* Centreline starting points *)
    ListPlot[
      Table[{r, 0}, {r, rCentreValues}]
    ],
    {}
    , Frame -> True
    , PlotRange -> {{0, 1.4}, {-1.4, 1.4}}
  ]
] // Ex["modification/wedge_small-modification-contours.pdf"]


(* ::Subsection:: *)
(*Check contours are actually contours*)


Module[
  {
    alpha, gamma,
    rCentreValues, heightValues, rphiContourList,
    hNumerical, tNumerical,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {alpha, gamma} = {apdMod, gpdMod} Degree;
  (* Import contours *)
  {rCentreValues, heightValues, rphiContourList} =
    Import["modification/wedge_small-modification-contours.txt"] // Uncompress;
  (* Import numerical solution *)
  hNumerical =
    Import @ FString["solution/wedge_small-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
      // Uncompress // First;
  tNumerical[r_, phi_] := hNumerical[r, phi] / r // Evaluate;
  (* Make plot *)
  Show[
    (* Pre-computed heights *)
    Plot[
      heightValues
      , {s, -2, 2}
      , PlotStyle -> Black
      , PlotOptions[Axes] // Evaluate
    ],
    (* Evaluated heights along curves *)
    Table[
      Plot[
        tNumerical @@ EvaluatePair[rphi, s]
        , {s, DomainStart[rphi], DomainEnd[rphi]}
        , PlotStyle -> Directive[Dotted, Green]
      ]
      , {rphi, rphiContourList}
    ],
    {}
    , AxesLabel -> {Italicise["s"], "Height"}
    , AxesOrigin -> {0, 0}
    , PlotRange -> Full
  ]
] // Ex["modification/wedge_small-modification-contours-check.pdf"]


(* ::Subsection:: *)
(*Visualise finite element meshes*)


Module[
  {
    rCentreValues, heightValues, rphiContourList,
    numContours,
    meshes,
      rphi, mesh,
      rUpperJoin, rLowerJoin,
    dummyForTrailingCommas
  },
  (* Import contours *)
  {rCentreValues, heightValues, rphiContourList} =
    Import["modification/wedge_small-modification-contours.txt"] // Uncompress;
  numContours = Length[rphiContourList];
  (* Import meshes *)
  meshes = Import["modification/wedge_small-modification-meshes.txt"] // Uncompress;
  (* For each contour: *)
  Table[
    (* Make plot *)
    Show[
      (* Mesh *)
      mesh = meshes[[n, 1]];
      mesh["Wireframe"],
      (* Contour *)
      rphi = rphiContourList[[n]];
      ParametricPlot[
        XYPolar @@ EvaluatePair[rphi, s] // Evaluate
        , {s, DomainStart[rphi], DomainEnd[rphi]}
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
      // Ex @ FString["modification/wedge_small-modification-mesh-{n}.png"]
    , {n, numContours}
  ]
]


(* ::Subsection:: *)
(*Height comparison*)


Module[
  {
    heightValues,
    profilesStuff, profiles,
      xyList, heightList, yList,
    dummyForTrailingCommas
  },
  (* Import contours *)
  heightValues =
    Import["modification/wedge_small-modification-contours.txt"]
      // Uncompress // #[[2]] &;
  (* Import profiles *)
  profilesStuff =
    Import["modification/wedge_small-modification-profiles.txt"]
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
    Plot[heightValues, {y, -2, 2}
      , PlotStyle -> Dashed
    ],
    (* Profile heights *)
    ListPlot[profiles
      , Joined -> True
    ],
    {}
    , AxesLabel -> {Italicise["y"], "Height"}
    , AxesOrigin -> {-2, 1}
    , PlotRange -> Full
    , PlotOptions[Axes] // Evaluate
  ]
] // Ex["modification/wedge_small-modification-contours-comparison.pdf"]


(* ::Section:: *)
(*Figure: original wedge & formal rectangle (wedge_small-domain-*)*)


Module[
  {
    rMax, alpha, numSubdivisions,
    textStyleLabel, opts, paddingLabel, paddingNothing,
    plot,
    xArrow, yArrow,
    dummyForTrailingCommas
  },
  (* Geometry *)
  rMax = 2;
  alpha = 25 Degree;
  numSubdivisions = 4;
  (* Options *)
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  opts = {Nothing
    , ImageSize -> 0.3 ImageSizeTextWidth
    , PlotRange -> {{0, rMax}, rMax Sin[alpha] {-1, 1}}
  };
  paddingLabel = Scaled[0.12];
  paddingNothing = Scaled[0.02];
  (* Original wedge diagram *)
  plot["original-wedge"] = Show[
    Graphics @ {
      (* Constant r (circular arcs) *)
      Table[
        Circle[{0, 0}, r, {-alpha, alpha}]
        , {r, Subdivide[0, rMax, numSubdivisions]}
      ],
      (* Constant \[Phi] (radial rays) *)
      Table[
        Line @ {{0, 0}, XYPolar[rMax, phi]}
        , {phi, Subdivide[-alpha, alpha, numSubdivisions]}
      ],
      (* Labels *)
      Text[
        Italicise["r"] // textStyleLabel
        , XYPolar[rMax / 2, -alpha]
        , {1.7, 0.8}
      ],
      Text[
        "\[Phi]" // textStyleLabel
        , {rMax, 0}
        , {-2.7, 0}
      ],
      {}
    },
    {}
    , PlotRangePadding -> {{paddingNothing, paddingLabel}, paddingNothing}
    , opts
  ];
  (* Formal rectangle diagram *)
  plot["formal-rectangle"] = Show[
    Graphics @ {
      (* Constant r (circular arcs) *)
      Table[
        Line @ {{r, -alpha}, {r, alpha}}
        , {r, Subdivide[0, rMax, numSubdivisions]}
      ],
      (* Constant \[Phi] (radial rays) *)
      Table[
        Line @ {{0, phi}, {rMax, phi}}
        , {phi, Subdivide[-alpha, alpha, numSubdivisions]}
      ],
      (* Labels *)
      Text[
        Italicise["r"] // textStyleLabel
        , {rMax / 2, -alpha}
        , {0, 1}
      ],
      Text[
        "\[Phi]" // textStyleLabel
        , {0, 0}
        , {3, 0}
      ],
      {}
    },
    {}
    , PlotRangePadding -> {{paddingLabel, paddingNothing}, paddingNothing}
    , opts
  ];
  (* Double-arrow between *)
  xArrow = 1;
  yArrow = 0.13;
  plot["double-arrow"] = Show[
    Graphics @ {
      GeneralStyle["VeryThick"], Arrowheads[0.2],
      LightGray,
      Arrow @ {{-xArrow, yArrow}, {xArrow, yArrow}},
      Arrow @ {{xArrow, -yArrow}, {-xArrow, -yArrow}}
    },
    {}
    , ImageSize -> 0.25 ImageSizeTextWidth
    , PlotRange -> All
    , PlotRangePadding -> {yArrow, 2 yArrow}
  ];
  (* Export *)
  Table[
    plot[case]
      // Ex @ FString["wedge_small-domain-{case}.pdf"]
    , {case, {"original-wedge", "formal-rectangle", "double-arrow"}}
  ]
]


(* ::Section:: *)
(*Figure: mesh detail (wedge_small-mesh-detail)*)


Module[
  {
    apd, alpha, rMax,
    mesh, meshWireframe,
    dummyForTrailingCommas
  },
  (* Geometry *)
  apd = 30;
  alpha = apd * Degree;
  rMax = rMaxMesh;
  (* Import mesh *)
  mesh =
    Import @ FString["mesh/wedge_small-mesh-apd-{apd}.txt"]
      // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  (* Mesh detail *)
  (* (not bothering with full mesh) *)
  Show[
    meshWireframe
    , Frame -> {{True, False}, {True, True}}
    , FrameLabel -> {
        Italicise["r"] // Margined @ {{0, 0}, {0, -0.03 ImageSizeTextWidth}},
        "\[Phi]" // LaTeXStyle
      }
    , FrameTicks -> {Automatic,
        Table[
          {
            phipd * Degree,
            SeparatedRow["VeryThin"][phipd, Style["\[Degree]", Magnification -> 1.2]]
          }
          , {phipd, Subdivide[-apd, apd, 4]}
        ]
      }
    , FrameTicksStyle -> LabelSize["Tick"]
    , ImageSize -> 0.8 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {{0, 2}, {-alpha, alpha}}
    , PlotRangeClipping -> True
    , PlotOptions[Frame] // Evaluate
  ]
    // Ex["wedge_small-mesh-detail.pdf"]
]


(* ::Section:: *)
(*Figure: numerical solution (wedge_small-solution-*)*)


Module[
  {
    apd, alpha,
    gpdValues, gpdCases,
    rMax, phiValues,
    opts,
    case,
    hSolution, tSolution,
    rArrow,
    hArrowMin, hArrowMax,
    tArrowMin, tArrowMax,
    arrowStyle, textStyle,
    hPlot, tPlot,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 30;
  alpha = apd * Degree;
  (* Contact angles *)
  gpdValues = {45, 60};
  gpdCases = {"small", "borderline"};
  (* Plot range *)
  rMax = 6;
  phiValues = Range[0, apd, 10] Degree;
  opts = {Nothing
    , ImagePadding -> {{0.5, 0.5}, {0.5, 0.7}} ImageSizeCentimetre
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotPoints -> 2
    , PlotStyle -> Black
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes]
  };
  (* For each solution contact angle: *)
  Table[
    case = AssociationThread[gpdValues -> gpdCases][gpd];
    (* Import numerical solution *)
    hSolution =
      FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Import // Uncompress // First;
    tSolution[r_, phi_] := hSolution[r, phi] / r;
    (* Phi arrow *)
    rArrow = 1.7;
    hArrowMin = hSolution[rArrow, 0] - 0.4;
    hArrowMax = hSolution[rArrow, alpha] + 0.7;
    tArrowMin = tSolution[rArrow, 0] - 0.25;
    tArrowMax = tSolution[rArrow, alpha] + 0.45;
    arrowStyle = Directive[
      GeneralStyle["Translucent"],
      GeneralStyle["Thick"],
      Arrowheads[Medium]
    ];
    textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
    (* Make plots *)
    hPlot =
      Plot[
        Table[hSolution[r, phi], {phi, phiValues}] // Evaluate
        , {r, 0, rMax}
        , AxesLabel -> Italicise /@ {"r", "H"}
        , Epilog -> {
            {arrowStyle,
              Arrow @ {{rArrow, hArrowMin}, {rArrow, hArrowMax}}
            },
            {
              Text["\[Phi]" // textStyle
                , {rArrow, hArrowMax}
                , {0, -1}
              ]
            }
          }
        , PlotRange -> {0, 3}
        , opts // Evaluate
      ];
    tPlot =
      Plot[
        Table[tSolution[r, phi], {phi, phiValues}] // Evaluate
        , {r, 0, rMax}
        , AxesLabel -> Italicise /@ {"r", "T"}
        , Epilog -> {
            {arrowStyle,
              Arrow @ {{rArrow, tArrowMin}, {rArrow, tArrowMax}}
            },
            {
              Text["\[Phi]" // textStyle
                , {rArrow, tArrowMax}
                , {0, -1}
              ]
            }
          }
        , PlotRange -> {0, 2.1}
        , opts // Evaluate
      ];
    GraphicsColumn[{hPlot, tPlot}
      , ImageSize -> 0.45 ImageSizeTextWidth
      , ItemAspectRatio -> 0.75
    ]
      // Ex @ FString["wedge_small-solution-{case}.pdf"]
    , {gpd, gpdValues}
  ]
]


(* ::Section:: *)
(*Figure: borderline comparison (wedge_small-borderline-comparison)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    tAcute, tAcutePolar,
    rMax, rValues,
    relativeDiscrepancyTable,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  apd = 30;
  gpd = 60;
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Import wedge_acute solution T == T (x, y) *)
  tAcute =
    FString["../wedge_acute/solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tAcutePolar[r_, phi_] := tAcute @@ XYPolar[r, phi];
  (* Values to be plotted *)
  rMax = 5;
  rValues = Subdivide[0, rMax, 20] // Rest // Prepend[10^-6];
  relativeDiscrepancyTable =
    Table[
      Table[
        {
          r,
          tSmallPolar[r, phi] / tAcutePolar[r, phi] - 1
            // Abs
        }
        , {r, rValues}
      ]
      , {phi, {0, alpha}}
    ];
  (* Make plot and export *)
  ListLogPlot[
    relativeDiscrepancyTable
    , AxesLabel -> {
        Italicise["r"] // Margined @ {{0, 1}, {5, 0}},
        Style["Relative discrepancy", Smaller]
      }
    , ImageSize -> {Automatic, 0.375 ImageSizeTextWidth}
    , Joined -> True
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotLegends -> LineLegend[
        {
          LaTeXStyle["\[Phi]"] == 0,
          LaTeXStyle["\[Phi]"] == LaTeXStyle["\[Alpha]"]
        }
        , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
      ]
    , PlotMarkers -> {{Automatic, 6}, {"OpenMarkers", 5}}
    , PlotStyle -> Black
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes] // Evaluate
  ]
    // Ex["wedge_small-borderline-comparison.pdf"]
]


(* ::Section:: *)
(*Figure: (non-)viable domain  (wedge_small-viable)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    derivativeList,
    pTilde, qTilde, fTilde, phiTilde,
    rCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore, rMaxMore,
    plot,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {30, 45};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Derivative list for boundary tracing *)
  derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
    tracingDerivativeList[hSmall, gamma];
  (* Critical terminal point r_0 *)
  rCritical = r0[pTilde, gamma];
  (* Plot range *)
  xMax = Ceiling[2 rCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  rMaxMore = rMax + more;
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
      Line @ {
        xMaxMore {1, Tan[alpha]},
        {0, 0},
        xMaxMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      phiTilde @@ RPhiPolar[x, y] < 0
      , {x, rCritical, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 11
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
  ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Wall", "Terminal"},
      {"wedge wall", "terminal curve"}
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
      , ImageSize -> 0.45 ImageSizeTextWidth
    ]
      // Ex["wedge_small-viable.pdf"]
    ,
    GraphicsColumn[Join[legendCurves, legendRegions]
      , Alignment -> Left
      , ImageSize -> 0.3 ImageSizeTextWidth
      , ItemAspectRatio -> 0.135
      , Spacings -> {0, 0}
    ]
      // Ex["wedge_small-viable-legend.pdf"]
  }
]


(* ::Section:: *)
(*Figure: generic traced boundaries  (wedge_small-traced-boundaries)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    derivativeList,
    pTilde, qTilde, fTilde, phiTilde,
    rCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore, rMaxMore,
    nMax, nValues, correctionFactor,
    rStartList, rphiStartList,
      rStartTerminal, phiStartTerminal,
    rphiTracedList,
    plot,
    textStyle, textStyleBracket,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {30, 45};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Derivative list for boundary tracing *)
  derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
    tracingDerivativeList[hSmall, gamma];
  (* Critical terminal point r_0 *)
  rCritical = r0[pTilde, gamma];
  (* Plot range *)
  xMax = Ceiling[1.5 rCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  rMaxMore = rMax + more;
  (*
    NOTE: wall integration doesn't work
    because the traced boundary ventures slightly outside the wedge domain
    and Mathematica complains about Indeterminate values.
    So the wall boundaries are drawn manually at the end.
  *)
  (* Starting points for boundary tracing *)
  nMax = 10;
  nValues = Join[
    Range[nMax - 4],
    {nMax - 3.1, nMax - 2.2, nMax - 1/2}
  ];
  correctionFactor = 1.05; (* to make boundary pass through (x_0, 0) *)
  rStartList = Table[correctionFactor * n / nMax * rMax, {n, nValues}];
  rphiStartList = Table[{r, -alpha}, {r, rStartList}];
    (* Append a starting point along the terminal curve *)
    rStartTerminal = Way[rCritical, rMax, 0.56];
    phiStartTerminal = SeekRoot[
      phiTilde[rStartTerminal, #] &,
      {0, alpha},
      5
    ];
    rphiStartList = rphiStartList // Append @ {rStartTerminal, phiStartTerminal};
  (* Traced boundaries (upper branch) *)
  rphiTracedList =
    Table[
      tracedBoundary[derivativeList][
        rphiStart, 0, {0, 2 rMax}
        , -1, 1
        , 0, Function[{r, phi}, r Cos[phi] > xMaxMore]
      ]
      , {rphiStart, rphiStartList}
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
      Line @ {
        xMaxMore {1, Tan[alpha]},
        {0, 0},
        xMaxMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      phiTilde @@ RPhiPolar[x, y] < 0
      , {x, rCritical, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 6
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        XYPolar @@ EvaluatePair[rphi, s]
          // IncludeYReflection
          // Evaluate
        , {s, DomainStart[rphi], DomainEnd[rphi]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {rphi, rphiTracedList}
    ],
    (* Traced boundary walls *)
    ParametricPlot[
      XYPolar[r, alpha]
        // IncludeYReflection
        // Evaluate
      , {r, 0, rMaxMore}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {rCritical, 0}
    },
    textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
    textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
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
        {rCritical, 0},
        {-1.4, -0.1}
      ] // textStyle,
      {}
    },
    {}
  ];
  (* Export *)
  (* (re-use legend from "wedge_acute-traced-boundaries-legend.pdf") *)
  Show[plot
    , ImageSize -> 0.45 ImageSizeTextWidth
  ] // Ex["wedge_small-traced-boundaries.pdf"]
]


(* ::Section:: *)
(*Figure: traced boundaries through (x_0, 0) (wedge_small-traced-boundaries-hyperbolic)*)


(* ::Subsection:: *)
(*Both*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    derivativeList,
    pTilde, qTilde, fTilde, phiTilde,
    rCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore, rMaxMore,
    rphi,
    textStyle, textStyleBracket,
    plot,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {30, 45};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Derivative list for boundary tracing *)
  derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
    tracingDerivativeList[hSmall, gamma];
  (* Critical terminal point r_0 *)
  rCritical = r0[pTilde, gamma];
  (* Plot range *)
  xMax = Ceiling[2 rCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  rMaxMore = rMax + more;
  (* Traced boundary through (x_0, 0) (upper branch) *)
  rphi =
    Quiet[
      tracedBoundary[derivativeList][
        {rCritical, 0}, 0, 2 rMax {-1, 1}
        , -1, 1
        , -Infinity, Function[{r, phi}, r Cos[phi] > xMaxMore]
      ]
      , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Plot *)
  plot = Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMaxMore {1, Tan[alpha]},
        {0, 0},
        xMaxMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)(*
    RegionPlot[
      phiTilde @@ RPhiPolar[x, y] < 0
      , {x, rCritical, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 11
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],*)
    (* Traced boundaries *)
    ParametricPlot[
      XYPolar @@ EvaluatePair[rphi, s]
        // IncludeYReflection
        // Evaluate
      , {s, DomainStart[rphi], DomainEnd[rphi]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {rCritical, 0}
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
        {rCritical, 0},
        {1.45, -0.1}
      ] // textStyle,
      {}
    },
    {}
  ];
  (* Export *)
  Show[plot
    , ImageSize -> 0.45 * 0.75 ImageSizeTextWidth
  ]
    // Ex["wedge_small-traced-boundaries-hyperbolic-both.pdf"]
]


(* ::Subsection:: *)
(*Rounding*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    derivativeList,
    pTilde, qTilde, fTilde, phiTilde,
    rCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore, rMaxMore,
    rphi,
    textStyle, textStyleBracket,
    plot,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {30, 45};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Derivative list for boundary tracing *)
  derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
    tracingDerivativeList[hSmall, gamma];
  (* Critical terminal point r_0 *)
  rCritical = r0[pTilde, gamma];
  (* Plot range *)
  xMax = Ceiling[2 rCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  rMaxMore = rMax + more;
  (* Traced boundary through (x_0, 0) (upper branch) *)
  rphi =
    Quiet[
      tracedBoundary[derivativeList][
        {rCritical, 0}, 0, 2 rMax {0, 1}
        , -1, 1
        , -Infinity, Function[{r, phi}, r Cos[phi] > xMaxMore]
      ]
      , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Plot *)
  plot = Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
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
      XYPolar @@ EvaluatePair[rphi, s]
        // IncludeYReflection
        // Evaluate
      , {s, DomainStart[rphi], DomainEnd[rphi]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {rCritical, 0}
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
        {rCritical, 0},
        {1.45, -0.1}
      ] // textStyle,
      {}
    },
    {}
  ];
  (* Export *)
  Show[plot
    , ImageSize -> 0.45 * 0.75 ImageSizeTextWidth
  ]
    // Ex["wedge_small-traced-boundaries-hyperbolic-rounding.pdf"]
]


(* ::Section:: *)
(*Figure: traced boundaries, different contact angle   (wedge_small-traced-boundaries-different-angle-*)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    gpdTracingValues, caseNameList,
    gammaTracing,
    derivativeList, pTilde, qTilde, fTilde, phiTilde,
    rCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    rStartList, rphiStartList,
    rphiTracedList,
    plot,
    caseName,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {30, 45};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Different tracing contact angle *)
  gpdTracingValues = {40, 55};
  caseNameList = {"less", "more"};
  Table[
    (* BEGIN Table content *)
    gammaTracing = gpdTracing * Degree;
    (* Derivative list for boundary tracing *)
    derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
      tracingDerivativeList[hSmall, gammaTracing];
    (* Critical terminal point r_0 *)
    rCritical = r0[pTilde, gammaTracing];
    (* Plot range *)
    xMax =
      If[gammaTracing < gamma,
        2.5 rCritical,
        2 rCritical
      ];
    yMax = xMax Tan[alpha];
    rMax = RPolar[xMax, yMax];
    (* Plot range but more *)
    more = 0.1;
    xMaxMore = xMax + more;
    yMaxMore = yMax + more;
    (* Starting points *)
    If[
      gammaTracing < gamma,
        rStartList = Subdivide[rCritical, 4] // Rest;
        rphiStartList = Table[{r, -alpha}, {r, rStartList}];
      ,
        rStartList = Subdivide[rMax, 7] // Rest;
        rphiStartList = Join[
          Table[{r, +alpha}, {r, rStartList[[2 ;; ;; 2]]}],
          Table[{r, -alpha}, {r, rStartList}],
          {{10^-6, 0}},
          {}
        ];
    ];
    (* Traced boundaries (upper branch) *)
    rphiTracedList =
      Table[
        Quiet[
          tracedBoundary[derivativeList][
            rphiStart, 0, {0, 2 rMax}
            , -1, 1
            , 0, Function[{r, phi}, r Cos[phi] > xMaxMore]
          ]
          , {InterpolatingFunction::femdmval, NDSolveValue::nlnum, NDSolveValue::nrnum1}
        ]
        , {rphiStart, rphiStartList}
      ];
    (* Make plot *)
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
               phiTilde @@ RPhiPolar[x, y] < 0
               , {x, rCritical, xMaxMore}
               , {y, -yMaxMore, yMaxMore}
               , BoundaryStyle -> None
               , PlotPoints -> 5
               , PlotStyle -> BoundaryTracingStyle["NonViable"]
             ],
            ContourPlot[
               phiTilde @@ RPhiPolar[x, y] == 0
               , {x, rCritical, yMaxMore}
               , {y, -yMaxMore, yMaxMore}
               , ContourLabels -> None
               , ContourStyle -> BoundaryTracingStyle["Terminal"]
               , PlotPoints -> 7
             ]
           },
           {
             RegionPlot[
               phiTilde @@ RPhiPolar[x, y] < 0
               , {x, rCritical, xMaxMore}
               , {y, -yMax, yMax}
               , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
               , PlotPoints -> 10
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
            XYPolar @@ EvaluatePair[rphi, s]
              // IncludeYReflection
              // Evaluate
            , {s, DomainStart[rphi], DomainEnd[rphi]}
            , PlotPoints -> 2
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ]
          , {rphi, rphiTracedList}
        ],
        {}
      ];
    (* Export *)
    caseName = AssociationThread[gpdTracingValues, caseNameList][gpdTracing];
    Show[plot
      , ImageSize -> 0.43 ImageSizeTextWidth
    ] // Ex @ FString["wedge_small-traced-boundaries-different-angle-{caseName}.pdf"]
    (* END Table content *)
    , {gpdTracing, gpdTracingValues}
  ]
]


(* ::Section:: *)
(*Figure: different contact angle offset, top view (wedge_small-different-angle-offset-top)*)


Module[
  {
    apd, gpd, gpdTracing,
    alpha, gamma, gammaTracing,
    hSmall, tSmallPolar,
    derivativeList,
    pTilde, qTilde, fTilde, phiTilde,
    rCritical,
    rEnd, phiEnd, sMax,
    rphi, sEnd,
    xEnd, yEnd,
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
  (* Angular parameters *)
  {apd, gpd, gpdTracing} = {30, 45, 60};
  {alpha, gamma, gammaTracing} = {apd, gpd, gpdTracing} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Derivative list for boundary tracing *)
  derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
    tracingDerivativeList[hSmall, gammaTracing];
  (* Critical terminal point r_0 *)
  rCritical = r0[pTilde, gammaTracing];
  (* Ending r for traced boundary *)
  rEnd = 2.5 rCritical;
  sMax = 2 rEnd;
  (* Traced boundary through (x_0, 0) (upper branch) *)
  rphi =
    Quiet[
      tracedBoundary[derivativeList][
        {rCritical, 0}, 0, {0, sMax}
        , -1, 1
        , -Infinity, Function[{r, phi}, r > rEnd]
      ]
      , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
    ];
  sEnd = DomainEnd[rphi];
  {rEnd, phiEnd} = EvaluatePair[rphi, sEnd];
  {xEnd, yEnd} = XYPolar[rEnd, phiEnd];
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
  arrowLength = 5 d;
  orthogonalityMarkerLength = 0.5 d;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[None]];
  wallThicknessCorrection = -0.1;
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
          , {-2.5, -0.5}
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
          , {-1.6, -0.4}
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
          , {-1.3, 0.3}
          , XYPolar[1, alpha]
        ]
      },
      (* Offset wedge *)
      (* (separate lines to ensure dotting at the corner) *)
      Graphics @ {BoundaryTracingStyle["Contour"],
        Line @ {
          {d / Sin[alpha], 0},
          {xMaxOffset, +yMaxOffset} * 1.01 (* optical correction *)
        },
        Line @ {
          {d / Sin[alpha], 0},
          {xMaxOffset, -yMaxOffset}
        },
        {}
      },
      (* Traced boundaries *)
      ParametricPlot[
        XYPolar @@ EvaluatePair[rphi, s]
          // IncludeYReflection
          // Evaluate
        , {s, 0, 0.98 sEnd}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
      (* Critical terminal point (x_0, 0) *)
      Graphics @ {
        GeneralStyle["Point"],
        Point @ {rCritical, 0}
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
          , {rCritical, 0}
          , {-1.5, -0.2}
        ],
        {}
      },
      {}
      , PlotRange -> All
      , PlotRangePadding -> 0.1
    ];
  (* Export *)
  Show[plot
    , ImageSize -> 0.5 ImageSizeTextWidth
  ]
    // Ex["wedge_small-different-angle-offset-top.pdf"]
]


(* ::Section:: *)
(*Figure: family of candidate boundaries (wedge_small-candidates)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    gpdTracingValues,
    pTildeTemp, xCriticalMin, xCriticalMax,
    xMax, yMax, rMax,
    more, xMaxMore,
    gammaTracing,
    derivativeList, pTilde, qTilde, fTilde, phiTilde,
    rCritical,
    sMax, rphiTraced,
    rphi,
    xArrowMin, xArrowMax,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {30, 45};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Tracing contact angles *)
  gpdTracingValues = Range[gpd, 75, 5];
  (* Plot range *)
  pTildeTemp = tracingDerivativeList[hSmall, gamma] // First;
  xCriticalMin = r0[pTildeTemp, Min[gpdTracingValues] Degree];
  xCriticalMax = r0[pTildeTemp, Max[gpdTracingValues] Degree];
  xMax = 1.5 xCriticalMax;
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.1;
  xMaxMore = xMax + more;
  (* Compute traced boundary candidates *)
  Table[
    (* BEGIN Table content *)
    gammaTracing = gpdTracing * Degree;
    (* Derivative list for boundary tracing *)
    derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
      tracingDerivativeList[hSmall, gammaTracing];
    (* Critical terminal point x_0 *)
    rCritical = r0[pTilde, gammaTracing];
    (* Traced boundary (upper branch) *)
    sMax = 2 rMax;
    rphiTraced[gpdTracing] =
      Quiet[
        tracedBoundary[derivativeList][
          {rCritical, 0}, 0, {0, sMax}
          , -1, 1
          , -Infinity, Function[{r, phi}, r Cos[phi] > xMaxMore]
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ];
    (* END Table content *)
    , {gpdTracing, gpdTracingValues}
  ];
  (* Make plot *)
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
      rphi = rphiTraced[gpdTracing];
      ParametricPlot[
        XYPolar @@ EvaluatePair[rphi, Abs[s], Sign[s]]
          // Evaluate
        , {s, -DomainEnd[rphi], DomainEnd[rphi]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {gpdTracing, gpdTracingValues}
    ],
    (* Manually extend gammaTracing == gamma case *)
    {
      rphi = rphiTraced[gpd];
      ParametricPlot[
        Way[
          XYPolar @@ EvaluatePair[rphi, DomainEnd[rphi]],
          xMaxMore {1, Tan[alpha]}
          , prop
        ]
          // IncludeYReflection
          // Evaluate
        , {prop, 0, 1}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    },
    (* Tracing gamma arrow *)
    xArrowMin = Way[xCriticalMin, xCriticalMax, -0.13];
    xArrowMax = Way[xCriticalMin, xCriticalMax, +1.27];
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
    , ImageSize -> 0.5 ImageSizeTextWidth
  ] // Ex["wedge_small-candidates.pdf"]
]


(* ::Section:: *)
(*Figure: small regime candidate boundaries grouped with offset applied (wedge_small-small-candidates-offset-alpha-*_degree)*)


(* (This is slow (~10 sec).) *)
Module[
  {
    apdValues, gpdTracing, gammaTracing,
    alpha, gpdValues,
    gamma,
    hSmall, tSmallPolar,
    derivativeList, pTilde, qTilde, fTilde, phiTilde,
    rCritical, d, xCriticalOffset,
    sMax, rphiTraced,
    plotRangeFactor, xMax, yMax, xMaxMore,
    graphicsWedgeWalls, rphi, plotTracedBoundaries,
    mainPlot,
    xCriticalOffsetMinGamma, xCriticalOffsetMaxGamma,
    arrowDirection, xArrowBase, xArrowTip,
    xGammaLabel, gammaLabelXShift,
    xyFun, sBulging, xBulging, yBulging,
    deltaXInset, deltaYInset,
    xMinInset, xMaxInset, yMinInset, yMaxInset,
    insetScaledCoordinates,
    insetRectangleStyle, insetRectangle,
    insetPlot, combinedPlot, plot,
    dummyForTrailingCommas
  },
  (* Prescribed angles *)
  apdValues = {30, 45, 60};
  gpdTracing = 25;
  gammaTracing = gpdTracing * Degree;
  (* For various wedge half-angles: *)
  Table[
    alpha = apd * Degree;
    (* Known solution angles *)
    gpdValues = Range[5, gpdTracing, 5];
    Which[
      apd == apdValues[[1]],
        gpdValues = Drop[gpdValues, {2, 3}];
      ,
      apd == apdValues[[2]],
        gpdValues = Drop[gpdValues, {2, 3}];
      ,
      apd == apdValues[[3]],
        gpdValues = Part[gpdValues, {1, -1}];
    ];
    (* For various known solution angles: *)
    Table[
      gamma = gpd * Degree;
      (* Import numerical solutions *)
      hSmall[gpd] =
        FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // First;
      tSmallPolar[gpd][r_, phi_] := hSmall[gpd][r, phi] / r;
      (* Derivative list for boundary tracing *)
      derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
        tracingDerivativeList[hSmall[gpd], gammaTracing];
      (* Critical terminal point (x_0, 0) *)
      rCritical[gpd] = r0[pTilde, gammaTracing];
      (* Offset distance *)
      d[gpd] = DHalfPlane[gamma, gammaTracing];
      (* Offset critical terminal point (x'_0, 0) *)
      xCriticalOffset[gpd] = rCritical[gpd] - d[gpd] / Sin[alpha];
      (* Compute candidate boundary *)
      sMax = 2; (* Probably enough, NOT a proven upper bound *)
      rphiTraced[gpd] =
        Quiet[
          tracedBoundary[derivativeList][
            {rCritical[gpd], 0}, 0, {0, sMax}
            , -1, 1
            , -Infinity
          ]
          , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
        ];
      , {gpd, gpdValues}
    ];
    (* Plot range *)
    plotRangeFactor =
      Which[
        apd == apdValues[[1]], 1.85,
        apd == apdValues[[2]], 1.7,
        apd == apdValues[[3]], 2.5
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
        rphi = rphiTraced[gpd];
        ParametricPlot[
          Subtract[
            XYPolar @@ EvaluatePair[rphi, Abs[s], Sign[s]],
            {d[gpd] / Sin[alpha], 0}
          ] // Evaluate
          , {s, -DomainEnd[rphi], DomainEnd[rphi]}
          , PlotPoints -> 4
          , PlotStyle -> Directive[
              BoundaryTracingStyle["Traced"],
              GeneralStyle["DefaultThick"]
            ]
        ]
        , {gpd, gpdValues}
      ];
    (* Make main plot *)
    mainPlot = Show[
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
    ];
    (* Make inset to show bulging *)
    If[apd == apdValues[[3]], combinedPlot = mainPlot,
      (* Position of bulge *)
      rphi = rphiTraced[Min @ gpdValues];
      xyFun =
        Subtract[
          XYPolar @@ EvaluatePair[rphi, #],
          {d[Min @ gpdValues] / Sin[alpha], 0}
        ] &;
      sBulging =
        SeekFirstRootBisection[
          PhiPolar @@ xyFun[#] - alpha &,
          {0, 0.6},
          10
        ];
      {xBulging, yBulging} = xyFun[sBulging];
      (* Inset bounding box size and scaled positioning coordinates *)
      Which[
        apd == apdValues[[1]],
          deltaXInset = deltaYInset = 0.07;
          insetScaledCoordinates = {0.33, 0.85};
          ,
        apd == apdValues[[2]],
          deltaXInset = deltaYInset = 0.03;
          insetScaledCoordinates = {0.26, 0.87};
      ];
      (* Inset rectangle *)
      xMinInset = xBulging - deltaXInset;
      xMaxInset = xBulging + deltaXInset;
      yMinInset = yBulging - deltaYInset;
      yMaxInset = yBulging + deltaYInset;
      insetRectangleStyle = Directive[EdgeForm[Black], FaceForm[None]];
      insetRectangle = Rectangle[{xMinInset, yMinInset}, {xMaxInset, yMaxInset}];
      (* Inset *)
      insetPlot =
        Show[
          graphicsWedgeWalls,
          plotTracedBoundaries,
          Graphics @ {insetRectangleStyle, insetRectangle},
          {}
          , ImageSize -> 0.105 ImageSizeTextWidth
          , PlotRange -> {{xMinInset, xMaxInset}, {yMinInset, yMaxInset}}
          , PlotRangePadding -> Scaled[0.01]
        ];
      (* Combined *)
      combinedPlot =
        Show[
          mainPlot,
          Graphics @ {insetRectangleStyle, insetRectangle},
          {}
          , Epilog -> {Inset[insetPlot, Scaled[insetScaledCoordinates]]}
        ];
    ];
    (* Final plot (with frame) *)
    plot =
      Show[combinedPlot
        , FrameLabel -> {
            Superscript[
              Italicise["x"],
              Style["\[NegativeVeryThinSpace]\[Prime]", Magnification -> 1.3]
            ] // Margined @ {{0, 0}, {0, -15}},
            Italicise["y"]
          }
        , FrameTicksStyle -> LabelSize["Tick"]
        , ImageSize -> 0.31 ImageSizeTextWidth
        , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      ];
    (* Export plot *)
    {
      {"apd", apd},
      {"gpdValues", gpdValues},
      plot // Ex @ FString[
        "wedge_small-small-candidates-offset-alpha-{apd}_degree.pdf"
      ]
    }
    , {apd, apdValues}
  ] // Column
] // Ex["wedge_small-small-candidates-offset-log.pdf"]


(* ::Section:: *)
(*Figure: moderate regime candidate boundaries grouped with offset applied (wedge_small-moderate-candidates-offset-alpha-*_degree)*)


(* (This is slow (~30 sec).) *)
(*
  Here we plot the 'small-regime-known-solution boundaries'
  on top of the 'moderate-regime-known-solution boundaries'
  in "wedge_acute-candidates-offset-alpha-{apd}_degree.pdf".
*)
Module[
  {
    x0,
    apdValues, gpdTracing, gammaTracing,
    alpha,
      gpdValuesModerate, gpdValuesSmall,
      xMaxModerate,
      xCriticalOffsetMinGamma, xCriticalOffsetMaxGamma,
    moderatePlot, smallPlot, plot,
    arrowDirection, xArrowBase, xArrowTip,
    xGammaLabel, gammaLabelXShift,
    dummyForTrailingCommas
  },
  (* Definition of x0 from wedge_acute.wl *)
  x0[tSol_, gammaTra_] :=
    Module[{p, xMax},
      (* \[PartialD]T/\[PartialD]x *)
      p = Derivative[1, 0][tSol];
      (* Find x == x_0 for which -\[PartialD]T/\[PartialD]x == cot(gamma-tilde) *)
      xMax = 10;
      SeekRoot[p[#, 0] + Cot[gammaTra] &, {0, xMax}]
    ];
  (* Prescribed angles *)
  apdValues = {30, 45, 60};
  gpdTracing = 75;
  gammaTracing = gpdTracing * Degree;
  (* For various wedge half-angles: *)
  Table[
    alpha = apd * Degree;
    (*
      ----------------------------------------------------------------
      PART 1: produce the existing moderate-regime plots (grey)
      ----------------------------------------------------------------
    *)
    Module[
      {
        gpdValues,
        tNumerical,
        sMax,
        xyTraced,
        gamma, xCritical,
        derList, p, q, grad2, f, vi,
        d, xCriticalOffset,
        plotRangeFactor,
        xMax, yMax, xMaxMore,
        xy,
        dummyForTrailingCommasModerate
      },
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
      gpdValuesModerate = gpdValues;
      (* Import numerical solutions *)
      Table[
        tNumerical[gpd] =
          Import @ FString[
            "../wedge_acute/solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"
          ] // Uncompress // First;
        , {gpd, gpdValues}
      ];
      (* Maximum arc-length for boundary tracing *)
      sMax = 3; (* Probably enough, NOT a proven upper bound *)
      (* Compute traced boundary candidates *)
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
      xCriticalOffsetMaxGamma = xCriticalOffset[Max @ gpdValues];
      (* Plot range *)
      plotRangeFactor =
        Which[
          apd == apdValues[[1]], 1.8,
          apd == apdValues[[2]], 2,
          apd == apdValues[[3]], 2.2
        ];
      xMax = plotRangeFactor * Max[xCriticalOffset /@ gpdValues];
      xMaxModerate = xMax;
      yMax = xMax;
      xMaxMore = 1.1 xMax;
      (* Make plot *)
      moderatePlot = Show[
        EmptyFrame[{0, xMax}, {-yMax, yMax}
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
          xy = xyTraced[gpd];
          ParametricPlot[
            Evaluate @ Subtract[
              EvaluatePair[xy, Abs[s], Sign[s]],
              {d[gpd] / Sin[alpha], 0}
            ]
            , {s, -DomainEnd[xy], DomainEnd[xy]}
            , PlotPoints -> 4
            , PlotStyle -> BoundaryTracingStyle["Background"]
          ]
          , {gpd, gpdValues}
        ],
        {}
      ];
    ];
    (*
      ----------------------------------------------------------------
      PART 2: produce the new small-regime boundaries (black)
      ----------------------------------------------------------------
    *)
    Module[
      {
        gpdValues,
        gamma,
        hSmall, tSmallPolar,
        derivativeList, pTilde, qTilde, fTilde, phiTilde,
        rCritical, d, xCriticalOffset,
        sMax, rphiTraced,
        rphi,
        dummyForTrailingCommasSmall
      },
      (* Known solution angles *)
      gpdValues = Range[5, 90 - apd, 5];
      Which[
        apd == apdValues[[1]],
          gpdValues = Part[gpdValues, {1, -5, -2, -1}];
        ,
        apd == apdValues[[2]],
          gpdValues = Part[gpdValues, {1, -1}];
        ,
        apd == apdValues[[3]],
          gpdValues = Part[gpdValues, {-1}];
      ];
      gpdValuesSmall = gpdValues;
      (* For various known solution angles: *)
      Table[
        gamma = gpd * Degree;
        (* Import numerical solutions *)
        hSmall[gpd] =
          FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
            // Import // Uncompress // First;
        tSmallPolar[gpd][r_, phi_] := hSmall[gpd][r, phi] / r;
        (* Derivative list for boundary tracing *)
        derivativeList = {pTilde, qTilde, fTilde, phiTilde} =
          tracingDerivativeList[hSmall[gpd], gammaTracing];
        (* Critical terminal point (x_0, 0) *)
        rCritical[gpd] = r0[pTilde, gammaTracing];
        (* Offset distance *)
        d[gpd] = DHalfPlane[gamma, gammaTracing];
        (* Offset critical terminal point (x'_0, 0) *)
        xCriticalOffset[gpd] = rCritical[gpd] - d[gpd] / Sin[alpha];
        (* Compute candidate boundary *)
        sMax = 3; (* Probably enough, NOT a proven upper bound *)
        rphiTraced[gpd] =
          Quiet[
            tracedBoundary[derivativeList][
              {rCritical[gpd], 0}, 0, {0, sMax}
              , -1, 1
              , -Infinity
            ]
            , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
          ];
        , {gpd, gpdValues}
      ];
      xCriticalOffsetMinGamma = xCriticalOffset[Min @ gpdValues];
      (* Make plot *)
      smallPlot = Show[
        Table[
          rphi = rphiTraced[gpd];
          ParametricPlot[
            Subtract[
              XYPolar @@ EvaluatePair[rphi, Abs[s], Sign[s]],
              {d[gpd] / Sin[alpha], 0}
            ] // Evaluate
            , {s, -DomainEnd[rphi], DomainEnd[rphi]}
            , PlotPoints -> 4
            , PlotStyle -> Directive @@ {
                BoundaryTracingStyle["Traced"],
                GeneralStyle["DefaultThick"],
                If[apd + gpd == 90, AbsoluteDashing @ {5, 8}, Nothing]
              }
          ]
          , {gpd, gpdValues}
        ]
      ];
    ];
    (*
      ----------------------------------------------------------------
      PART 3: combine and export
      ----------------------------------------------------------------
    *)
    (* Make plot *)
    plot = Show[
      moderatePlot,
      smallPlot,
      (* Solution gamma arrow *)
      arrowDirection = Sign[xCriticalOffsetMaxGamma - xCriticalOffsetMinGamma];
      xArrowBase = xCriticalOffsetMinGamma - 0.1 xMaxModerate * arrowDirection;
      xArrowTip = xCriticalOffsetMaxGamma + 0.2 xMaxModerate * arrowDirection;
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
      , FrameLabel -> {
          Superscript[
            Italicise["x"],
            Style["\[NegativeVeryThinSpace]\[Prime]", Magnification -> 1.3]
          ] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , ImageSize -> 0.31 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ];
    (* Export *)
    {
      {"apd", apd},
      {"gpdValuesSmall", gpdValuesSmall},
      {"gpdValuesModerate", gpdValuesModerate},
      plot // Ex @ FString["wedge_small-moderate-candidates-offset-alpha-{apd}_degree.pdf"]
    }
    , {apd, apdValues}
  ] // Column
] // Ex["wedge_small-moderate-candidates-offset-log.pdf"]


(* ::Section:: *)
(*Figure: modified boundary contours*)


Module[
  {
    alpha, gamma,
    rCentreValues, heightValues, rphiContourList,
    numContours,
    hNumerical, derList, tNumerical, phiTilde,
    xMax, yMax,
    xMaxWall, yMaxWall, rMaxWall,
    style,
      rphi,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {alpha, gamma} = {apdMod, gpdMod} Degree;
  (* Import contours *)
  {rCentreValues, heightValues, rphiContourList} =
    Import["modification/wedge_small-modification-contours.txt"] // Uncompress;
  numContours = Length[rphiContourList] - 1;
  (* Import numerical solution *)
  hNumerical =
    Import @ FString["solution/wedge_small-solution-apd-{apdMod}-gpd-{gpdMod}.txt"]
      // Uncompress // First;
  derList = tracingDerivativeList[hNumerical, gamma];
  tNumerical[r_, phi_] := hNumerical[r, phi] / r // Evaluate;
  phiTilde = Last @ tracingDerivativeList[hNumerical, gamma];
  (* Make plot *)
  xMax = 1.4 Max[rCentreValues];
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
      phiTilde @@ RPhiPolar[x, y] < 0
      , {x, 0, xMaxWall}
      , {y, -yMaxWall, yMaxWall}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Contours *)
    Table[
      {
        rphi = rphiContourList[[n]];
        ParametricPlot[
          XYPolar @@ EvaluatePair[rphi, s] // Evaluate
          , {s, DomainStart[rphi], DomainEnd[rphi]}
          , PlotPoints -> 2
          , PlotStyle -> style[n]
        ],
        {}
      }
      , {n, numContours}
    ],
    {}
    , FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -18}},
        Italicise["y"]
      }
    , ImageSize -> 0.45 * 0.7 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , FrameTicksStyle -> LabelSize["Tick"]
  ]
] // Ex["wedge_small-modification-contours.pdf"]


(* ::Section:: *)
(*Figure: modified boundary height rise profiles*)


Module[
  {
    contourHeightValues,
    profilesStuff,
    numContours,
      xyList, heightList,
      yList, yMin, yMax,
      profilePointsList,
    yMaxPlot, heightMin,
    style,
      textStyle,
      yBracket, bracketWidth,
      heightMinBracket, heightMaxBracket,
    dummyForTrailingCommas
  },
  (* Import contours *)
  contourHeightValues =
    Import["modification/wedge_small-modification-contours.txt"]
      // Uncompress // #[[2]] &;
  (* Import profiles *)
  profilesStuff =
    Import["modification/wedge_small-modification-profiles.txt"]
      // Uncompress;
  (* Process contours to be shown *)
  numContours = Length[profilesStuff] - 1;
  Table[
    {xyList, heightList} = profilesStuff[[n]];
    yList = xyList[[All, 2]];
    {yMin[n], yMax[n]} = MinMax[yList];
    profilePointsList[n] = {yList, heightList} // Transpose;
    , {n, numContours}
  ];
  (* Make plot *)
  yMaxPlot = 1;
  heightMin = 1;
  style[n_] := If[
    n > numContoursViable,
    BoundaryTracingStyle["ContourPlain"],
    Black
  ];
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    (* Main *)
    Table[
      {
        (* Contour heights *)
        Plot[contourHeightValues[[n]], {y, yMin[n], yMax[n]}
          , PlotStyle -> Directive[style[n], Dotted]
        ],
        (* Profiles *)
        ListPlot[profilePointsList[n]
          , Joined -> True
          , PlotRange -> All
          , PlotStyle -> style[n]
        ],
        {}
      }
      , {n, numContours}
    ],
    (* Viable bracket label *)
    yBracket = Way[
      yMax[numContoursViable],
      yMax[numContoursViable + 1]
      , 2/3
    ];
    bracketWidth = 2/3 yMax[1];
    heightMinBracket = Way[
      contourHeightValues[[numContoursViable]],
      contourHeightValues[[numContoursViable + 1]]
      , 1/2
    ];
    heightMaxBracket = contourHeightValues[[1]] + 2/3 bracketWidth;
    Graphics @ {
      Line @ {
        {yBracket - bracketWidth, heightMinBracket},
        {yBracket, heightMinBracket},
        {yBracket, heightMaxBracket},
        {yBracket - bracketWidth, heightMaxBracket},
        Nothing
      },
      Text["viable" // textStyle
        , {yBracket, Way[heightMinBracket, heightMaxBracket]}
        , {-1.3, 0}
      ],
      {}
    },
    {}
    , AspectRatio -> Automatic
    , AxesLabel -> {Italicise["y"], "Height"}
    , AxesOrigin -> {-yMaxPlot, heightMin}
    , ImageSize -> 0.5 * 0.9 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {{-yMaxPlot, yMaxPlot}, {heightMin, Full}}
    , PlotRangeClipping -> False
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["wedge_small-modification-height-rise-profiles.pdf"]
