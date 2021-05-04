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
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "wedge_obtuse"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Finite element mesh*)


(* Wedge half-angles (apd means "alpha per degree") *)
apdMeshValues = {120, 135, 150};
rMaxMesh = 10;


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.txt",
    Module[
      {
        rMax, alpha,
        fineLengthScale, coarseLengthScale,
        boundaryPoints, numBoundaryPoints, boundaryElements,
        boundaryMesh, mesh,
        predicateWet,
        dummyForTrailingCommas
      },
      (* Geometry *)
      rMax = rMaxMesh;
      alpha = apd * Degree;
      (* Refinement length scales *)
      fineLengthScale = 0.01;
      coarseLengthScale = 0.2;
      (* Boundary points and elements *)
      boundaryPoints =
        Join[
          (* Lower wedge wall *)
          Table[
            XYPolar[r, -alpha]
            , {r, UniformRange[0, rMax, fineLengthScale]}
          ]
            // Most,
          (* Far arc at infinity *)
          Table[
            XYPolar[rMax, phi]
            , {phi, UniformRange[-alpha, alpha, coarseLengthScale / rMax]}
          ]
            // Most,
          (* Upper wedge wall *)
          Table[
            XYPolar[r, +alpha]
            , {r, UniformRange[rMax, 0, -fineLengthScale]}
          ]
            // Most,
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
      (* Predicate function for wetting boundaries *)
      predicateWet =
        Function[{x, y},
          Abs[y] <= rMax Sin[Pi - alpha] && x <= 0
          // Evaluate
        ];
      (* Return *)
      {mesh, predicateWet} // Compress
    ]
  ]
  , {apd, apdMeshValues}
]


(* ::Subsection:: *)
(*Solve capillary BVP*)


(* Solution contact angles (gpd means "gamma per degree") *)
gpdSolutionValues = Join[{1}, Range[5, 85, 5]];


(* (These are slow (~10 min), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
(* For each value of alpha: *)
Table[
  (* If all values of contact angle have been account for: *)
  If[
    Table[
      FString @ "solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"
      , {gpd, gpdSolutionValues}
    ] // AllTrue[FileExistsQ]
    ,
    (* Then skip: *)
    Null
    ,
    (* Otherwise: *)
    Module[{mesh, predicateWet, gamma, evaluationData},
      (* Import mesh *)
      {mesh, predicateWet} =
        FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.txt"
          // Import // Uncompress;
      (* For each value of gamma: *)
      Table[
        gamma = gpd * Degree;
        ExportIfNotExists[
          FString @ "solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"
          ,
          evaluationData = EvaluationData @ SolveLaplaceYoung[gamma, mesh, predicateWet];
          evaluationData /@ {"Result", "Success", "FailureType", "MessagesText", "AbsoluteTiming"}
            // Compress
        ]
        , {gpd, gpdSolutionValues}
      ]
    ]
  ]
  , {apd, apdMeshValues}
]


(* ::Subsection:: *)
(*Elliptic critical terminal point*)


x0[tNumerical_, gammaTracing_] :=
  Module[{p},
    (* \[PartialD]T/\[PartialD]x *)
    p = Derivative[1, 0][tNumerical];
    (* Find x == x_0 such that \[PartialD]T/\[PartialD]x == cot(gamma-tracing) *)
    SeekRoot[
      p[#, 0] + Cot[gammaTracing] &,
      {0, rMaxMesh}
    ]
  ];


(* ::Subsection:: *)
(*Dip-coating with indentations*)


(* ::Subsubsection:: *)
(*Angular parameters*)


apdIndentations = 135;
gpdIndentations = 20;
gpdTracingIndentations = 55;


{alphaIndentations, gammaIndentations, gammaTracingIndentations} =
  {apdIndentations, gpdIndentations, gpdTracingIndentations} Degree;


(* ::Subsubsection:: *)
(*Numerical solution*)


tNumericalIndentations =
  Import @ FString[
    "solution/wedge_obtuse-solution-apd-{apdIndentations}-gpd-{gpdIndentations}.txt"
  ] // Uncompress // First;


(* ::Subsubsection:: *)
(*Contour to be approximated*)


xStartContourIndentations =
  Module[{hTracing},
    hTracing = HHalfPlane[gammaTracingIndentations];
    SeekRoot[
      tNumericalIndentations[#, 0] - hTracing &,
      {0, rMaxMesh}
      , 20
    ]
  ];


xyContourIndentations =
  Quiet[
    ContourByArcLength[tNumericalIndentations][
      {xStartContourIndentations, 0}
      , 0, {-1, 1} 2 rMaxMesh
    ]
    , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
  ];


(* ::Subsubsection:: *)
(*Roughness profile*)


(* ::Subsubsubsection:: *)
(*Raw version*)


(* rho == rho (x, y) *)
rhoIndentationsRaw =
  Module[{gammaTracing, grad2},
    gammaTracing = gammaTracingIndentations;
    grad2 = ContactDerivativeList[tNumericalIndentations, gammaTracing][[3]];
    Function[{x, y},
      Divide[
        Sqrt @ grad2[x, y],
        Cos[gammaTracing] Sqrt[1 + grad2[x, y]]
      ] // Evaluate
    ]
  ];


(* rho == rho (s) *)
rhoIndentationsRawProfile =
  Function[{s},
    rhoIndentationsRaw @@ EvaluatePair[xyContourIndentations, s]
      // Evaluate
  ];


(* ::Subsubsubsection:: *)
(*Fitted version*)


(* rho == rho (s) model *)
rhoIndentationsFittedProfileModel =
  Module[
    {
      sStart, sEnd,
      sDataCutoff, sDataValues, data,
      rhoModel,
      dummyForTrailingCommas
    },
    (* Domain in arc-length *)
    sStart = DomainStart[xyContourIndentations];
    sEnd = DomainEnd[xyContourIndentations];
    (* Cutoff for data points when roughness is near 1 *)
    sDataCutoff = SeekRoot[
      rhoIndentationsRawProfile[#] - (1 + 10^-4) &,
      {0, sEnd}
    ];
    (* Build table of data for fitting *)
    sDataValues = Subdivide[-sDataCutoff, sDataCutoff, 2 * 50];
    data = Table[{s, rhoIndentationsRawProfile[s]}, {s, sDataValues}];
    (* Determine fit *)
    With[{s = \[FormalS], rho0 = \[FormalRho], sChar = \[FormalA]},
      rhoModel =
        NonlinearModelFit[
          data,
          1 + (rho0 - 1) Exp[-Abs[s] / sChar],
          {{rho0, 1.4}, {sChar, 1}},
          s
        ]
    ];
    (* Return fit *)
    rhoModel
  ];


(* rho == rho (s) *)
rhoIndentationsFittedProfile =
  Function[{s}, rhoIndentationsFittedProfileModel[s] // Evaluate];


(* ::Subsubsection:: *)
(*Effective contact angle profile*)


(* ::Subsubsubsection:: *)
(*Raw version*)


(* gamma_e == gamma_e (x, y) *)
gammaEffectiveRaw =
  Module[{gammaTracing, grad2},
    gammaTracing = gammaTracingIndentations;
    grad2 = ContactDerivativeList[tNumericalIndentations, gammaTracing][[3]];
    Function[{x, y},
      ArcCos @ Divide[
        Sqrt @ grad2[x, y],
        Sqrt[1 + grad2[x, y]]
      ] // Evaluate
    ]
  ];


(* gamma_e == gamma_e (s) *)
gammaEffectiveRawProfile =
  Function[{s},
    gammaEffectiveRaw @@ EvaluatePair[xyContourIndentations, s]
      // Evaluate
  ];


(* ::Subsubsection:: *)
(*Triangular groove width values*)


sigmaValues = {0.2, 0.1, 0.05, 0.02, 0.01};


(* ::Subsubsection:: *)
(*Triangular groove angle value*)


phiIndentations =
  Max[
    ArcCos[1 / rhoIndentationsRawProfile[0]],
    ArcCos[1 / rhoIndentationsFittedProfile[0]]
  ] / Degree // Ceiling[#, 5] & // # Degree &


phiIndentations < gammaTracingIndentations


(* ::Subsubsection:: *)
(*Triangular groove positions (arc-length)*)


(*
  Note that we use the smoothed version of roughness
  to avoid issues with rho < 1 (resulting in lambda < 0).
  We only bother with the top half (y > 0) due to symmetry.
  We compute pairs {sLower, sUpper} for the endpoints
  of the triangular grooves.
  Note that sUpper - sLower is not exactly equal
  to the groove width sigma, because the boundary is curved.
*)


(* Straight-line distance between points along contour *)
distanceContourIndentations[s1_, s2_] :=
  RPolar @@ Subtract[
    EvaluatePair[xyContourIndentations, s2],
    EvaluatePair[xyContourIndentations, s1]
  ] // Evaluate;


(* Compute positions *)
sGroovePairsList[];
Table[
  Module[
    {
      sMax, nMax,
      sLower, sUpper, sPairsList,
      rho, lambda, sLowerNext,
      dummyForTrailingCommas
    },
    (* Maximum arc-length cutoff *)
    sMax = 10;
    (* Crude upper bound for numer of grooves *)
    nMax = Ceiling[sMax / sigma];
    (* Initialise pairs list with central groove *)
    sUpper = SeekRoot[
      distanceContourIndentations[#, -#] - sigma &,
      {0, 2 sigma}, 4
    ];
    sLower = -sUpper;
    sPairsList = {{sLower, sUpper}};
    (* Iteratively compute groove locations *)
    Do[
      (* Current roughness *)
      rho = rhoIndentationsFittedProfile @ Way[sLower, sUpper];
      (* Separation distance *)
      lambda = Divide[1 / Cos[phiIndentations] - rho, rho - 1] sigma;
      (* Next location *)
      sLowerNext = sUpper + lambda;
      (* Check whether finished *)
      If[sLowerNext > sMax - 2 sigma,
        Break[]
      ];
      (* Update location *)
      sLower = sLowerNext;
      sUpper = SeekRoot[
        distanceContourIndentations[sLower, #] - sigma &,
        {sLower, sLower + 2 sigma}, 4
      ];
      (* Append to pairs list *)
      sPairsList = sPairsList // Append @ {sLower, sUpper};
      , {n, 0, nMax}
    ];
    (* Store values *)
    sGroovePairsList[sigma] = sPairsList;
  ];
  , {sigma, sigmaValues}
];


(* ::Subsubsection:: *)
(*Triangular-groove indented boundary points*)


xyIndentationList[];
Table[
  Module[
    {
      fineLengthScale,
      grooveHeight, grooveHypotenuse,
      numGrooveSubdivisions,
      sPairsList, nMax,
      boundaryPointList,
      sLower, sUpper,
      xyLower, xyUpper,
      xyDisplacement, xyHalfway,
      xyVertex,
      sLowerNext,
      dummyForTrailingCommas
    },
    (* Fine length scale *)
    fineLengthScale = sigma / 5;
    (* Groove geometry *)
    grooveHeight = sigma/2 * Tan[phiIndentations];
    grooveHypotenuse = sigma/2 / Cos[phiIndentations];
    (* Number of subdivisions required for groove walls *)
    numGrooveSubdivisions = Ceiling[grooveHypotenuse / fineLengthScale];
    (* Get grooves list *)
    sPairsList = sGroovePairsList[sigma];
    nMax = Length[sPairsList];
    (* Initialise boundary point list *)
    boundaryPointList = {};
    (* For each groove: *)
    Do[
      (* Groove endpoint arc-lengths *)
      {sLower, sUpper} = sPairsList[[n]];
      (* Groove endpoint positions *)
      xyLower = EvaluatePair[xyContourIndentations, sLower];
      xyUpper = EvaluatePair[xyContourIndentations, sUpper];
      (* Displacement and halfway therebetween *)
      xyDisplacement = xyUpper - xyLower;
      xyHalfway = Way[xyLower, xyUpper];
      (* Groove vertex position *)
      xyVertex = xyHalfway + grooveHeight Cross @ Normalize[xyDisplacement];
      (* Append groove points *)
      boundaryPointList = Join[boundaryPointList,
        Subdivide[xyLower, xyVertex, numGrooveSubdivisions] // Most,
        Subdivide[xyVertex, xyUpper, numGrooveSubdivisions] // Most,
        {}
      ];
      (* Next groove position *)
      If[n < nMax,
        sLowerNext = sPairsList[[n + 1]] // First;
        ,
        sLowerNext = DomainEnd[xyContourIndentations];
      ];
      (* Append points up till next groove *)
      boundaryPointList = Join[boundaryPointList,
        Table[
          EvaluatePair[xyContourIndentations, s]
          , {s, UniformRange[sUpper, sLowerNext, fineLengthScale]}
        ]
      ];
      , {n, nMax}
    ];
    (* Symmetrise about y == 0 *)
    boundaryPointList = boundaryPointList[[numGrooveSubdivisions + 1 ;;]];
    boundaryPointList =
      Join[
        boundaryPointList /. {x_, y_} :> {x, -y} // Reverse,
        boundaryPointList // Rest
      ];
    xyIndentationList[sigma] = boundaryPointList;
  ]
  , {sigma, sigmaValues}
];


(* ::Subsubsection:: *)
(*Indented finite element mesh*)


(* (These are not slow, nevertheless compute once and store. *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString @ "mesh/wedge_obtuse-mesh-indented-sigma-{sigma}.txt",
    Module[
      {
        boundaryPointsIndented, xMaxIndented,
        rCorner, phiCorner,
        coarseLengthScale, boundaryPointsFarArc,
        boundaryPoints, numBoundaryPoints, boundaryElements,
        boundaryMesh, mesh,
        predicateWet,
        dummyForTrailingCommas
      },
      (* Indented boundary points *)
      boundaryPointsIndented = xyIndentationList[sigma];
      xMaxIndented = boundaryPointsIndented[[All, 1]] // Max;
      (* Corner where indented boundary meets far arc *)
      {rCorner, phiCorner} = RPhiPolar @@ Last[boundaryPointsIndented];
      (* Far arc boundary points *)
      coarseLengthScale = 0.5;
      boundaryPointsFarArc =
        Table[
          XYPolar[rCorner, ang]
          , {ang, UniformRange[-phiCorner, phiCorner, coarseLengthScale / rCorner]}
        ] // Rest // Most;
      (* Combined boundary points *)
      boundaryPoints = Join[
        boundaryPointsIndented // Reverse,
        boundaryPointsFarArc
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
      (* Predicate function for wetting boundaries *)
      predicateWet =
        Function[{x, y},
          And[
            Abs[y] <= rCorner Sin[Pi - phiCorner],
            x < xMaxIndented + sigma
          ] // Evaluate
        ];
      (* Return *)
      {mesh, predicateWet} // Compress
    ]
  ]
  , {sigma, sigmaValues}
]


(* ::Subsubsection:: *)
(*Indented BVP numerical solutions*)


(* (This is slow (~45 sec), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
Table[
  ExportIfNotExists[FString @ "solution/wedge_obtuse-solution-indented-sigma-{sigma}.txt",
    Module[{mesh, predicateWet, tNumerical},
      (* Import mesh *)
      {mesh, predicateWet} =
        FString @ "mesh/wedge_obtuse-mesh-indented-sigma-{sigma}.txt"
          // Import // Uncompress;
      (* Compute numerical solution to capillary BVP *)
      tNumerical = SolveLaplaceYoung[gammaTracingIndentations, mesh, predicateWet];
      (* Return *)
      tNumerical // Compress
    ]
  ]
  , {sigma, sigmaValues}
]


(* ::Subsubsection:: *)
(*Non-indented finite element mesh*)


(* (This is not slow, nevertheless compute once and store. *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists[FString @ "mesh/wedge_obtuse-mesh-non_indented.txt",
  Module[
    {
      fineLengthScale, sStart, sEnd,
      boundaryPointsWetting, xMaxWetting,
      rCorner, phiCorner,
      coarseLengthScale, boundaryPointsFarArc,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      boundaryMesh, mesh,
      predicateWet,
      dummyForTrailingCommas
    },
    (* Wetting boundary points (which is the T-contour) *)
    fineLengthScale = 0.01;
    sStart = DomainStart[xyContourIndentations];
    sEnd = DomainEnd[xyContourIndentations];
    boundaryPointsWetting =
      Table[
        EvaluatePair[xyContourIndentations, s]
        , {s, UniformRange[sStart, sEnd, fineLengthScale]}
      ];
    xMaxWetting = boundaryPointsWetting[[All, 1]] // Max;
    (* Corner where wetting boundary meets far arc *)
    {rCorner, phiCorner} = RPhiPolar @@ Last[boundaryPointsWetting];
    (* Far arc boundary points *)
    coarseLengthScale = 0.5;
    boundaryPointsFarArc =
      Table[
        XYPolar[rCorner, ang]
        , {ang, UniformRange[-phiCorner, phiCorner, coarseLengthScale / rCorner]}
      ] // Rest // Most;
    (* Combined boundary points *)
    boundaryPoints = Join[
      boundaryPointsWetting // Reverse,
      boundaryPointsFarArc
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
    (* Predicate function for wetting boundaries *)
    predicateWet =
      Function[{x, y},
        And[
          Abs[y] <= rCorner Sin[Pi - phiCorner],
          x <= xMaxWetting
        ] // Evaluate
      ];
    (* Return *)
    {mesh, predicateWet} // Compress
  ]
]


(* ::Subsubsection:: *)
(*Non-indented BVP numerical solution*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["solution/wedge_obtuse-solution-non_indented.txt",
  Module[{mesh, predicateWet, tNumerical},
    (* Import mesh *)
    {mesh, predicateWet} =
      "mesh/wedge_obtuse-mesh-non_indented.txt" // Import // Uncompress;
    (* Compute numerical solution to capillary BVP *)
    tNumerical = SolveLaplaceYoung[gammaTracingIndentations, mesh, predicateWet];
    (* Return *)
    tNumerical // Compress
  ]
]


(* ::Section:: *)
(*Finite element mesh check*)


(* (This is slow (~5 sec).) *)
Table[
  Module[
    {
      mesh, predicateWet, numElements,
      boundaryCoordinates, wetBoundaryCoordinates,
      plot,
      dummyForTrailingCommas
    },
    (* Import mesh *)
    {mesh, predicateWet} =
      FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.txt"
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
    (* Make plot *)
    plot =
      Show[
        mesh["Wireframe"],
        ListPlot[boundaryCoordinates, PlotStyle -> Directive[Red, PointSize[Medium]]],
        ListPlot[wetBoundaryCoordinates, PlotStyle -> Directive[Blue, PointSize[Medium]]],
        {}
        , PlotLabel -> FString @ "{numElements} mesh elements"
      ];
    (* Export *)
    {
      (* Full plot *)
      plot
        // Ex @ FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.png"
      ,
      (* Zoom near origin to verify wetting predicate *)
      Show[plot
        , ImageSize -> 480
        , PlotLabel -> None
        , PlotRange -> {{-0.15, 0.15}, {-0.15, 0.15}}
      ]
        // Ex @ FString @ "mesh/wedge_obtuse-mesh-apd-{apd}-zoom.png"
    }
  ]
  , {apd, apdMeshValues}
]


(* ::Section:: *)
(*Numerical solution*)


(* ::Subsection:: *)
(*Table*)


(* (This is slow (~30 sec) from all the calls to Import.) *)
Join @@ Table[
  Table[
    Module[{tNumerical, messages, timing, heightCorner},
      (* Import numerical solution etc. *)
      {tNumerical, messages, timing} =
        FString @ "solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"
          // Import // Uncompress // Part[#, {1, 4, 5}] &;
      (* Filter out uninformative messages *)
      messages = DeleteCases[messages, string_ /; StringContainsQ[string, "StringForm"]];
      (* Make pretty *)
      messages = StringRiffle[messages, "\n"];
      messages = messages /. {"" -> "\[HappySmiley]"};
      (* Compute corner height *)
      heightCorner = tNumerical[0, 0];
      (* Return row of table *)
      {apd * Degree, gpd * Degree, heightCorner, timing, messages}
    ]
    , {gpd, gpdSolutionValues}
  ]
  , {apd, apdMeshValues}
] // TableForm[#
  , TableAlignments -> {Left, Center}
  , TableDepth -> 2
  , TableHeadings -> {
      None,
      {"alpha", "gamma", "corner height", "absolute time", "messages"}
    }
] & // Ex["wedge_obtuse-solution-table.pdf"]


(* ::Subsection:: *)
(*Table for theoretical vs computed corner slopes*)


(* (This is slow (~10 sec) from all the calls to Import.) *)
Module[
  {
    apd, alpha,
    gamma, k,
    infiniteSlopeQ, theoreticalSlope,
    tNumerical, tXDerivative, tYDerivative,
    computedSlope, slopeRelError,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 135;
  alpha = apd * Degree;
  (* For each contact angle: *)
  Table[
    (* Theoretical corner slope *)
    gamma = gpd Degree;
    k = Sin[alpha] / Cos[gamma];
    infiniteSlopeQ = alpha >= Pi/2 + gamma;
    theoreticalSlope = If[infiniteSlopeQ, Infinity, 1 / Sqrt[k^2 - 1] // N];
    (* Import numerical solution *)
    tNumerical =
      Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    tXDerivative = Derivative[1, 0][tNumerical];
    tYDerivative = Derivative[0, 1][tNumerical];
    computedSlope = Sqrt[tXDerivative[0, 0]^2 + tYDerivative[0, 0]^2];
    slopeRelError = If[infiniteSlopeQ, "N/A", computedSlope / theoreticalSlope - 1];
    (* Build row of table *)
    {gamma, theoreticalSlope, computedSlope, slopeRelError}
    , {gpd, gpdSolutionValues}
  ]
    // TableForm[#
      , TableHeadings -> {
          None,
          {"gamma", "theoretical slope", "computed slope", "rel. error"}
        }
    ] &
] // Ex["wedge_obtuse-slope-comparison.pdf"]


(* ::Subsection:: *)
(*Corner neighbourhood test*)


(* A visual check of local linearity. *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, tXDerivative,
    coordMax, coordSliceValues,
    imageSize,
    dummyForTrailingCommas
  },
  (* Parameters *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  tXDerivative = Derivative[1, 0][tNumerical];
  (* Plot *)
  coordMax = 0.01;
  coordSliceValues = Subdivide[-coordMax, coordMax, 4];
  imageSize = 360;
  {
    (* Height *)
    Plot[
      Table[tNumerical[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "height"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["solution/wedge_obtuse-solution-neighbourhood-height.png"]
    ,
    (* x-derivative *)
    Plot[
      Table[tXDerivative[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "slope"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["solution/wedge_obtuse-solution-neighbourhood-slope.png"]
  }
]


(* ::Section:: *)
(*Slope dodginess test*)


(*
  Numerical slope near a re-entrant corner is dodgy.
  Here, we examine the effect of additional refinement for:
  * alpha == {45, 135} Degree
  * gamma == 60 Degree
  * global coarse length scale 0.2
  * wall fine length scale 0.01
  * near-corner (r < 0.01) fine length scales {0.01, 0.005, 0.002, 0.001}
*)


(* ::Subsection:: *)
(*Compute data and store*)


(* (This is slow (~1.5 min).) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
Module[
  {
    gamma,
    globalCoarseLengthScale, wallFineLengthScale, nearCornerRadius,
    apdValues, nearCornerFineLengthScaleValues,
      rMax, alpha,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      equilateralTriangleArea, refinementFunction,
      boundaryMesh, mesh, numMeshElements,
      predicateWet,
      tNumerical, tXDerivative, tYDerivative, tSlope,
      rValues,
        symmetryHeightData, symmetrySlopeData,
        wallHeightData, wallSlopeData,
    dummyForTrailingCommas
  },
  (* Constants *)
  gamma = 60 Degree;
  globalCoarseLengthScale = 0.2;
  wallFineLengthScale = 0.01;
  nearCornerRadius = 0.01;
  (* Values to test *)
  apdValues = {45, 135};
  nearCornerFineLengthScaleValues = {0.01, 0.001, 0.0005, 0.0002, 0.0001};
  (* For each pair of test values: *)
  Table[
    (* Build mesh *)
    rMax = rMaxMesh;
    alpha = apd * Degree;
    (* Boundary points and elements *)
    boundaryPoints =
      Join[
        (* Lower wedge wall *)
        Table[
          XYPolar[r, -alpha]
          , {r, UniformRange[0, rMax, wallFineLengthScale]}
        ]
          // Most,
        (* Far arc at infinity *)
        Table[
          XYPolar[rMax, phi]
          , {phi, UniformRange[-alpha, alpha, globalCoarseLengthScale / rMax]}
        ]
          // Most,
        (* Upper wedge wall *)
        Table[
          XYPolar[r, +alpha]
          , {r, UniformRange[rMax, 0, -wallFineLengthScale]}
        ]
          // Most,
        {}
      ];
    numBoundaryPoints = Length[boundaryPoints];
    boundaryElements =
      LineElement /@ {
        Table[{n, n + 1}, {n, numBoundaryPoints}]
          // Mod[#, numBoundaryPoints, 1] &
      };
    (* Special refinement *)
    equilateralTriangleArea = Function[{len}, Sqrt[3]/4 * len^2];
    refinementFunction =
      Function[{vertices, area},
        Or[
          area > equilateralTriangleArea[globalCoarseLengthScale],
          Norm @ Mean[vertices] < nearCornerRadius
            && area > equilateralTriangleArea[nearCornerFineLengthScale]
        ]
      ];
    (* Build mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> boundaryPoints,
        "BoundaryElements" -> boundaryElements,
        {}
      ];
    mesh = ToElementMesh[boundaryMesh
      , "ImproveBoundaryPosition" -> True
      , MeshRefinementFunction -> refinementFunction
    ];
    (* Export mesh wireframe for visual check *)
    numMeshElements = Length @ mesh[[2, 1, 1]];
    Show[
      mesh["Wireframe"]
      , Frame -> True
      , ImageSize -> 240
      , PlotLabel -> Column @ {
          FString @ "apd: {apd}",
          FString @ "nearCornerFineLengthScale: {nearCornerFineLengthScale}",
          FString @ "numMeshElements: {numMeshElements}",
          Nothing
        }
      , PlotRange -> 0.1 {{-1, 1}, {-1, 1}}
      , PlotRangeClipping -> True
    ] // Ex @ FString[
      "slope_test/wedge_obtuse-slope_test-apd-{apd}-fls-{nearCornerFineLengthScale}.png"
    ];
    (* Predicate function for wetting boundaries *)
    predicateWet =
      If[alpha < Pi/2,
        Function[{x, y}, x <= rMax Cos[alpha] // Evaluate],
        Function[{x, y}, Abs[y] <= rMax Sin[Pi - alpha] && x <= 0 // Evaluate]
      ];
    (* Compute numerical solution to capillary BVP *)
    tNumerical = SolveLaplaceYoung[gamma, mesh, predicateWet];
    tXDerivative = Derivative[1, 0][tNumerical];
    tYDerivative = Derivative[0, 1][tNumerical];
    tSlope = Function[{x, y}, Sqrt[tXDerivative[x, y]^2 + tYDerivative[x, y]^2]];
    rValues =
      Join[
        Range[0, 0.1, 0.001],
        Range[0.1, 1, 0.01] // Rest,
        Range[1, 5, 0.1] // Rest,
        {}
      ];
    (* Compute data (along line of symmetry) *)
    symmetryHeightData = Table[{r, tNumerical @@ XYPolar[r, 0]}, {r, rValues}];
    symmetrySlopeData = Table[{r, tSlope @@ XYPolar[r, 0]}, {r, rValues}];
    (* Compute data (along wall) *)
    wallHeightData = Table[{r, tNumerical @@ XYPolar[r, alpha]}, {r, rValues}];
    wallSlopeData = Table[{r, tSlope @@ XYPolar[r, alpha]}, {r, rValues}];
    (* Return row *)
    {
      apd, nearCornerFineLengthScale,
      symmetryHeightData, symmetrySlopeData,
      wallHeightData, wallSlopeData,
      Nothing
    }
    , {apd, apdValues}
    , {nearCornerFineLengthScale, nearCornerFineLengthScaleValues}
  ] // Compress
] // Ex["slope_test/wedge_obtuse-slope_test-all-data.txt"]


(* ::Subsection:: *)
(*Visualise*)


Module[
  {
    allData,
    apdValues, nearCornerFineLengthScaleValues,
    symmetryHeightData, symmetrySlopeData,
    cornerHeight, cornerSlope,
    wallHeightData, wallSlopeData,
    dummyForTrailingCommas
  },
  (* Import all data *)
  allData = Import["slope_test/wedge_obtuse-slope_test-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apdValues = allData[[All, 1]] // Union;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd_, ncfls_, shData_, ssData_, whData_, wsData_} :> (
      (* Along line of symmetry *)
      symmetryHeightData[apd, ncfls] = shData;
        cornerHeight[apd, ncfls] = shData[[1, 2]];
      symmetrySlopeData[apd, ncfls] = ssData;
        cornerSlope[apd, ncfls] = ssData[[1, 2]];
      (* Along wall *)
      wallHeightData[apd, ncfls] = whData;
      wallSlopeData[apd, ncfls] = wsData;
    )
  ];
  (* For each wedge half-angle: *)
  Table[
    List[
      (* Corner height plot *)
      ListLogLinearPlot[
        Table[{ncfls, cornerHeight[apd, ncfls]}, {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"ncfls", "corner height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> All
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-corner-height-apd-{apd}.png"]
      ,
      (* Symmetry height plot *)
      ListPlot[
        Table[symmetryHeightData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "symm height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.01}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-symm-height-apd-{apd}.png"]
      ,
      (* Corner slope plot *)
      ListLogLinearPlot[
        Table[{ncfls, -cornerSlope[apd, ncfls]}, {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"ncfls", "corner slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> All
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-corner-slope-apd-{apd}.png"]
      ,
      (* Symmetry slope plot *)
      ListPlot[
        Table[symmetrySlopeData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "symm slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.01}, All}
        , PlotRangeClipping -> False
        , PlotRangePadding -> {None, Scaled[0.1]}
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-symm-slope-apd-{apd}.png"]
      ,
      (* Wall height *)
      ListPlot[
        Table[wallHeightData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "wall height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.02}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-wall-height-apd-{apd}.png"]
      ,
      (* Wall slope *)
      ListPlot[
        Table[wallSlopeData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "wall slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.02}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-wall-slope-{apd}.png"]
    ]
    , {apd, apdValues}
  ]
]


(* ::Subsection:: *)
(*Corner neighbourhood test*)


(* (This is slow (~15 sec).) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
Module[
  {
    apd, gpd,
    alpha, gamma,
      globalCoarseLengthScale, wallFineLengthScale, nearCornerRadius,
      nearCornerFineLengthScale,
      rMax,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      equilateralTriangleArea, refinementFunction,
      boundaryMesh, mesh, numMeshElements,
      predicateWet,
    tNumerical, tXDerivative, tYDerivative,
    coordMax, coordSliceValues,
    imageSize,
    dummyForTrailingCommas
  },
  (* Parameters *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Build mesh and compute numerical solution *)
  (* (copied from above) *)
    (* Parameters *)
    globalCoarseLengthScale = 0.2;
    wallFineLengthScale = 0.01;
    nearCornerRadius = 0.01;
    nearCornerFineLengthScale = 0.0001;
    (* Geometry *)
    rMax = rMaxMesh;
    (* Boundary points and elements *)
    boundaryPoints =
      Join[
        (* Lower wedge wall *)
        Table[
          XYPolar[r, -alpha]
          , {r, UniformRange[0, rMax, wallFineLengthScale]}
        ]
          // Most,
        (* Far arc at infinity *)
        Table[
          XYPolar[rMax, phi]
          , {phi, UniformRange[-alpha, alpha, globalCoarseLengthScale / rMax]}
        ]
          // Most,
        (* Upper wedge wall *)
        Table[
          XYPolar[r, +alpha]
          , {r, UniformRange[rMax, 0, -wallFineLengthScale]}
        ]
          // Most,
        {}
      ];
    numBoundaryPoints = Length[boundaryPoints];
    boundaryElements =
      LineElement /@ {
        Table[{n, n + 1}, {n, numBoundaryPoints}]
          // Mod[#, numBoundaryPoints, 1] &
      };
    (* Special refinement *)
    equilateralTriangleArea = Function[{len}, Sqrt[3]/4 * len^2];
    refinementFunction =
      Function[{vertices, area},
        Or[
          area > equilateralTriangleArea[globalCoarseLengthScale],
          Norm @ Mean[vertices] < nearCornerRadius
            && area > equilateralTriangleArea[nearCornerFineLengthScale]
        ]
      ];
    (* Build mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> boundaryPoints,
        "BoundaryElements" -> boundaryElements,
        {}
      ];
    mesh = ToElementMesh[boundaryMesh
      , "ImproveBoundaryPosition" -> True
      , MeshRefinementFunction -> refinementFunction
    ];
    (* Predicate function for wetting boundaries *)
    predicateWet =
      If[alpha < Pi/2,
        Function[{x, y}, x <= rMax Cos[alpha] // Evaluate],
        Function[{x, y}, Abs[y] <= rMax Sin[Pi - alpha] && x <= 0 // Evaluate]
      ];
    (* Compute numerical solution to capillary BVP *)
    tNumerical = SolveLaplaceYoung[gamma, mesh, predicateWet];
    tXDerivative = Derivative[1, 0][tNumerical];
  (* Plot *)
  coordMax = 0.01;
  coordSliceValues = Subdivide[-coordMax, coordMax, 4];
  imageSize = 360;
  {
    (* Height *)
    Plot[
      Table[tNumerical[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "height"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["slope_test/wedge_obtuse-slope_test-neighbourhood-height.png"]
    ,
    (* x-derivative *)
    Plot[
      Table[tXDerivative[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "slope"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["slope_test/wedge_obtuse-slope_test-neighbourhood-slope.png"]
  }
] // AbsoluteTiming


(* ::Subsection:: *)
(*Table for theoretical vs computed corner slopes*)


Module[
  {
    allData,
    apd, alpha,
    gpd, gamma, k,
    theoreticalSlope,
    apdValues, nearCornerFineLengthScaleValues,
    computedSlope, slopeRelError,
    dummyForTrailingCommas
  },
  (* Import all data *)
  allData = Import["slope_test/wedge_obtuse-slope_test-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Selected value of wedge half-angle *)
  apd = 135;
  alpha = apd * Degree;
  (* Theoretical corner slope *)
  gpd = 60;
  gamma = gpd Degree;
  k = Sin[alpha] / Cos[gamma];
  theoreticalSlope = 1 / Sqrt[k^2 - 1];
  (* Extract computed corner slope data *)
  apdValues = allData[[All, 1]] // Union;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd, ncfls_, _, symmetrySlopeData_, __} :> (
      computedSlope[ncfls] = symmetrySlopeData[[1, 2]];
      slopeRelError[ncfls] = computedSlope[ncfls] / theoreticalSlope - 1;
    )
  ];
  (* Build table *)
  Table[
    {ncfls, computedSlope[ncfls], slopeRelError[ncfls]}
    , {ncfls, nearCornerFineLengthScaleValues // Reverse}
  ]
    // TableForm[#,
      , TableHeadings -> {
          None,
          {"ncfls", "computed slope", "rel. error"}
        }
    ] &
] // Ex["wedge_obtuse-slope-comparison-additional-refinement.pdf"]


(* ::Section:: *)
(*Slope dodginess test, large angle case*)


(*
  Same as "Slope dodginess test" above,
  but for the large angle case.
  We examine the effect of additional refinement for:
  * alpha == {135} Degree
  * gamma == 30 Degree
  * global coarse length scale 0.2
  * wall fine length scale 0.01
  * near-corner (r < 0.01) fine length scales {0.01, 0.005, 0.002, 0.001}
*)


(* ::Subsection:: *)
(*Compute data and store*)


(* (This is slow (~2 min).) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
Module[
  {
    gamma,
    globalCoarseLengthScale, wallFineLengthScale, nearCornerRadius,
    apdValues, nearCornerFineLengthScaleValues,
      rMax, alpha,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      equilateralTriangleArea, refinementFunction,
      boundaryMesh, mesh, numMeshElements,
      predicateWet,
      tNumerical, tXDerivative, tYDerivative, tSlope,
      rValues,
        symmetryHeightData, symmetrySlopeData,
        wallHeightData, wallSlopeData,
    dummyForTrailingCommas
  },
  (* Constants *)
  gamma = 30 Degree;
  globalCoarseLengthScale = 0.2;
  wallFineLengthScale = 0.01;
  nearCornerRadius = 0.01;
  (* Values to test *)
  apdValues = {135};
  nearCornerFineLengthScaleValues = {0.01, 0.001, 0.0005, 0.0002, 0.0001};
  (* For each pair of test values: *)
  Table[
    (* Build mesh *)
    rMax = rMaxMesh;
    alpha = apd * Degree;
    (* Boundary points and elements *)
    boundaryPoints =
      Join[
        (* Lower wedge wall *)
        Table[
          XYPolar[r, -alpha]
          , {r, UniformRange[0, rMax, wallFineLengthScale]}
        ]
          // Most,
        (* Far arc at infinity *)
        Table[
          XYPolar[rMax, phi]
          , {phi, UniformRange[-alpha, alpha, globalCoarseLengthScale / rMax]}
        ]
          // Most,
        (* Upper wedge wall *)
        Table[
          XYPolar[r, +alpha]
          , {r, UniformRange[rMax, 0, -wallFineLengthScale]}
        ]
          // Most,
        {}
      ];
    numBoundaryPoints = Length[boundaryPoints];
    boundaryElements =
      LineElement /@ {
        Table[{n, n + 1}, {n, numBoundaryPoints}]
          // Mod[#, numBoundaryPoints, 1] &
      };
    (* Special refinement *)
    equilateralTriangleArea = Function[{len}, Sqrt[3]/4 * len^2];
    refinementFunction =
      Function[{vertices, area},
        Or[
          area > equilateralTriangleArea[globalCoarseLengthScale],
          Norm @ Mean[vertices] < nearCornerRadius
            && area > equilateralTriangleArea[nearCornerFineLengthScale]
        ]
      ];
    (* Build mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> boundaryPoints,
        "BoundaryElements" -> boundaryElements,
        {}
      ];
    mesh = ToElementMesh[boundaryMesh
      , "ImproveBoundaryPosition" -> True
      , MeshRefinementFunction -> refinementFunction
    ];
    (* Export mesh wireframe for visual check *)
    numMeshElements = Length @ mesh[[2, 1, 1]];
    Show[
      mesh["Wireframe"]
      , Frame -> True
      , ImageSize -> 240
      , PlotLabel -> Column @ {
          FString @ "apd: {apd}",
          FString @ "nearCornerFineLengthScale: {nearCornerFineLengthScale}",
          FString @ "numMeshElements: {numMeshElements}",
          Nothing
        }
      , PlotRange -> 0.1 {{-1, 1}, {-1, 1}}
      , PlotRangeClipping -> True
    ] // Ex @ FString[
      "slope_test_large/wedge_obtuse-slope_test_large-apd-{apd}-fls-{nearCornerFineLengthScale}.png"
    ];
    (* Predicate function for wetting boundaries *)
    predicateWet =
      If[alpha < Pi/2,
        Function[{x, y}, x <= rMax Cos[alpha] // Evaluate],
        Function[{x, y}, Abs[y] <= rMax Sin[Pi - alpha] && x <= 0 // Evaluate]
      ];
    (* Compute numerical solution to capillary BVP *)
    tNumerical = SolveLaplaceYoung[gamma, mesh, predicateWet];
    tXDerivative = Derivative[1, 0][tNumerical];
    tYDerivative = Derivative[0, 1][tNumerical];
    tSlope = Function[{x, y}, Sqrt[tXDerivative[x, y]^2 + tYDerivative[x, y]^2]];
    rValues =
      Join[
        Range[0, 0.1, 0.001],
        Range[0.1, 1, 0.01] // Rest,
        Range[1, 5, 0.1] // Rest,
        {}
      ];
    (* Compute data (along line of symmetry) *)
    symmetryHeightData = Table[{r, tNumerical @@ XYPolar[r, 0]}, {r, rValues}];
    symmetrySlopeData = Table[{r, tSlope @@ XYPolar[r, 0]}, {r, rValues}];
    (* Compute data (along wall) *)
    wallHeightData = Table[{r, tNumerical @@ XYPolar[r, alpha]}, {r, rValues}];
    wallSlopeData = Table[{r, tSlope @@ XYPolar[r, alpha]}, {r, rValues}];
    (* Return row *)
    {
      apd, nearCornerFineLengthScale,
      symmetryHeightData, symmetrySlopeData,
      wallHeightData, wallSlopeData,
      Nothing
    }
    , {apd, apdValues}
    , {nearCornerFineLengthScale, nearCornerFineLengthScaleValues}
  ] // Compress
] // Ex["slope_test_large/wedge_obtuse-slope_test_large-all-data.txt"]


(* ::Subsection:: *)
(*Visualise*)


Module[
  {
    allData,
    apdValues, nearCornerFineLengthScaleValues,
    symmetryHeightData, symmetrySlopeData,
    cornerHeight, cornerSlope,
    wallHeightData, wallSlopeData,
    dummyForTrailingCommas
  },
  (* Import all data *)
  allData = Import["slope_test_large/wedge_obtuse-slope_test_large-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apdValues = allData[[All, 1]] // Union;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd_, ncfls_, shData_, ssData_, whData_, wsData_} :> (
      (* Along line of symmetry *)
      symmetryHeightData[apd, ncfls] = shData;
        cornerHeight[apd, ncfls] = shData[[1, 2]];
      symmetrySlopeData[apd, ncfls] = ssData;
        cornerSlope[apd, ncfls] = ssData[[1, 2]];
      (* Along wall *)
      wallHeightData[apd, ncfls] = whData;
      wallSlopeData[apd, ncfls] = wsData;
    )
  ];
  (* For each wedge half-angle: *)
  Table[
    List[
      (* Corner height plot *)
      ListLogLinearPlot[
        Table[{ncfls, cornerHeight[apd, ncfls]}, {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"ncfls", "corner height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> All
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test_large/wedge_obtuse-slope_test_large-corner-height-apd-{apd}.png"]
      ,
      (* Symmetry height plot *)
      ListPlot[
        Table[symmetryHeightData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "symm height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.01}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test_large/wedge_obtuse-slope_test_large-symm-height-apd-{apd}.png"]
      ,
      (* Corner slope plot *)
      ListLogLinearPlot[
        Table[{ncfls, -cornerSlope[apd, ncfls]}, {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"ncfls", "corner slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> All
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test_large/wedge_obtuse-slope_test_large-corner-slope-apd-{apd}.png"]
      ,
      (* Symmetry slope plot *)
      ListPlot[
        Table[symmetrySlopeData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "symm slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.01}, All}
        , PlotRangeClipping -> False
        , PlotRangePadding -> {None, Scaled[0.1]}
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test_large/wedge_obtuse-slope_test_large-symm-slope-apd-{apd}.png"]
      ,
      (* Wall height *)
      ListPlot[
        Table[wallHeightData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "wall height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.02}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test_large/wedge_obtuse-slope_test_large-wall-height-apd-{apd}.png"]
      ,
      (* Wall slope *)
      ListPlot[
        Table[wallSlopeData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "wall slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.02}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test_large/wedge_obtuse-slope_test_large-wall-slope-{apd}.png"]
    ]
    , {apd, apdValues}
  ]
]


(* ::Section:: *)
(*Roughness profile fit check*)


Module[{sMax, plot},
  sMax = 5;
  plot = Plot[
    {
      rhoIndentationsRawProfile[s],
      rhoIndentationsFittedProfile[s]
    }
    , {s, -sMax, sMax}
    , ImageSize -> Medium
    , PlotRange -> {0, All}
    , PlotLegends -> {"Raw", "Fitted"}
  ];
  Column @ {
    rhoIndentationsFittedProfileModel,
    Table[
      {prop, rhoIndentationsFittedProfileModel[prop]}
      , {prop, {"ParameterTable", "RSquared"}}
    ] // TableForm,
    plot,
    Nothing
  }
] // Ex["wedge_obtuse-roughness-profile-fit-check.pdf"]


(* ::Section:: *)
(*Indented finite element mesh check*)


(* (This is slow (~5 sec).) *)
Table[
  Module[
    {
      mesh, predicateWet, numElements,
      boundaryCoordinates, wetBoundaryCoordinates,
      plot,
      dummyForTrailingCommas
    },
    (* Import mesh *)
    {mesh, predicateWet} =
      FString @ "mesh/wedge_obtuse-mesh-indented-sigma-{sigma}.txt"
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
    (* Make plot *)
    plot =
      Show[
        mesh["Wireframe"],
        ListPlot[boundaryCoordinates, PlotStyle -> Directive[Red, PointSize[Small]]],
        ListPlot[wetBoundaryCoordinates, PlotStyle -> Directive[Blue, PointSize[Small]]],
        {}
        , PlotLabel -> FString @ "{numElements} mesh elements"
      ];
    (* Export *)
    {
      (* Full plot *)
      plot
        // Ex @ FString @ "mesh/wedge_obtuse-mesh-indented-sigma-{sigma}.png"
      ,
      (* Zoom near origin to verify wetting predicate *)
      Show[plot
        , ImageSize -> 480
        , PlotLabel -> None
        , PlotRange -> 8 sigma {{-1, 1}, {-1, 1}}
      ]
        // Ex @ FString @ "mesh/wedge_obtuse-mesh-indented-sigma-{sigma}-zoom.png"
    }
  ]
  , {sigma, sigmaValues}
]


(* ::Section:: *)
(*Indented numerical solution check*)


(* (This is slow (~10 sec).) *)
Table[
  Module[{tNumerical, mesh, x, y},
    (* Import numerical solution *)
    tNumerical =
      FString @ "solution/wedge_obtuse-solution-indented-sigma-{sigma}.txt"
        // Import // Uncompress;
    mesh = tNumerical["ElementMesh"];
    (* Make plot *)
    Plot3D[
      tNumerical[x, y], Element[{x, y}, mesh]
      , PlotRange -> Full
      , RegionFunction -> Function[{x, y}, RPolar[x, y] < 3]
    ] // Ex @ FString["solution/wedge_obtuse-solution-indented-sigma-{sigma}.png"]
  ]
  , {sigma, sigmaValues}
]


(* ::Section:: *)
(*Non-indented finite element mesh check*)


Module[
  {
    mesh, predicateWet, numElements,
    boundaryCoordinates, wetBoundaryCoordinates,
    plot,
    dummyForTrailingCommas
  },
  (* Import mesh *)
  {mesh, predicateWet} = "mesh/wedge_obtuse-mesh-non_indented.txt" // Import // Uncompress;
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
  (* Make plot *)
  plot =
    Show[
      mesh["Wireframe"],
      ListPlot[boundaryCoordinates, PlotStyle -> Directive[Red, PointSize[Medium]]],
      ListPlot[wetBoundaryCoordinates, PlotStyle -> Directive[Blue, PointSize[Medium]]],
      {}
      , PlotLabel -> FString @ "{numElements} mesh elements"
    ];
  (* Export *)
  {
    (* Full plot *)
    plot
      // Ex["mesh/wedge_obtuse-mesh-non_indented.png"]
    ,
    (* Zoom near origin to verify wetting predicate *)
    Show[plot
      , PlotLabel -> None
      , PlotRange -> 0.3 {{-0.2, 1}, {-1, 1}}
    ]
      // Ex["mesh/wedge_obtuse-mesh-non_indented-zoom.png"]
  }
]


(* ::Section:: *)
(*Non-indented numerical solution check*)


Module[{tNumerical, mesh, x, y},
  (* Import numerical solution *)
  tNumerical =
    "solution/wedge_obtuse-solution-non_indented.txt"
      // Import // Uncompress;
  mesh = tNumerical["ElementMesh"];
  (* Make plot *)
  Plot3D[
    tNumerical[x, y], Element[{x, y}, mesh]
    , PlotRange -> Full
  ] // Ex["solution/wedge_obtuse-solution-non_indented.png"]
]


(* ::Section:: *)
(*Figure: numerical wedge domain (wedge_obtuse-numerical-domain)*)


Module[
  {
    rMax, alpha,
    angleMarkerRadius,
    textStyle,
    dummyForTrailingCommas
  },
  (* Geometry *)
  rMax = rMaxMesh;
  alpha = 135 Degree;
  (* Diagram *)
  angleMarkerRadius = 0.15 rMax;
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    (* Domain boundaries *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      (* Wedge walls *)
      Line @ {XYPolar[rMax, -alpha], {0, 0}, XYPolar[rMax, +alpha]},
      (* Far arc at infinity *)
      Circle[{0, 0}, rMax, {-alpha, alpha}],
      {}
    },
    (* Wedge radius *)
    Graphics @ {
      Text[
        rMaxMesh // textStyle
        , XYPolar[rMax / 2, alpha]
        , {2, 0.25}
      ]
    },
    (* Wedge half-angle *)
    Graphics @ {
      Circle[{0, 0}, angleMarkerRadius, {0, alpha}]
    },
    Graphics @ {
      Text[
        "\[Alpha]" // textStyle
        , XYPolar[angleMarkerRadius, alpha/2]
        , {-2.5, -0.6}
      ]
    },
    {}
    , Axes -> True
    , AxesLabel -> Italicise /@ {"x", "y"}
    , ImageSize -> 0.48 * 0.9 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , Method -> {"AxesInFront" -> False}
    , PlotRangePadding -> Scaled[0.067]
    , Ticks -> None
  ]
] // Ex["wedge_obtuse-numerical-domain.pdf"]


(* ::Section:: *)
(*Figure: finite element mesh detail (wedge_obtuse-mesh-detail)*)


Module[
  {
    apd, alpha,
    mesh, meshWireframe,
    rMax, xMin, yMax, xMax,
    rTickValues, tickLength,
    tickTextStyle, axisTextStyle,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 135;
  alpha = apd * Degree;
  (* Import mesh *)
  mesh = Import @ FString["mesh/wedge_obtuse-mesh-apd-{apd}.txt"] // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  rMax = 0.25;
  {xMin, yMax} = XYPolar[rMax, alpha];
  xMax = Abs[xMin];
  (* Make plot *)
  tickTextStyle = Style[#, LabelSize["Tick"]] & @* LaTeXStyle;
  axisTextStyle = Style[#, LabelSize["Axis"]] & @* LaTeXStyle;
  Show[
    meshWireframe,
    (* Ticks *)
    rTickValues = Range[0, rMax, 0.1] // Rest;
    tickLength = 0.03 rMax;
    Graphics @ {
      Table[
        {
          Line @ {
            XYPolar[r, alpha],
            XYPolar[r, alpha] + XYPolar[tickLength, alpha + Pi/2]
          },
          Text[
            r // tickTextStyle
            , XYPolar[r, alpha] + XYPolar[1 tickLength, alpha + Pi/2]
            , {1.2, 0.65}
          ],
          {}
        }
        , {r, rTickValues}
      ]
    },
    {}
    , ImageSize -> 0.48 * 0.85 ImageSizeTextWidth
    , PlotRange -> {{xMin, xMax}, {-yMax, yMax}}
  ]
] // Ex["wedge_obtuse-mesh-detail.pdf"]


(* ::Section:: *)
(*Figure: numerical wedge solution (wedge_obtuse-solution)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    rMax, xArcCorner, yArcCorner,
    tNumerical,
    x, y, mesh,
    verticalEdge, wallHeight,
    dummyForTrailingCommas
  },
  (* Parameters *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Plot range *)
  rMax = 3;
  {xArcCorner, yArcCorner} = XYPolar[rMax, alpha];
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Plot *)
  mesh = tNumerical["ElementMesh"];
  wallHeight = Ceiling[HHalfPlane[gamma], 0.2];
  verticalEdge @ {x_, y_} :=
    Line @ {{x, y, tNumerical[x, y]}, {x, y, wallHeight}};
  Show[
    (* Numerical solution *)
    Plot3D[
      tNumerical[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"],
          Italicise["y"] // Margined @ {{10, 0}, {0, -15}},
          Italicise["T"] // Margined @ {{0, 5}, {0, 0}}
        }
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Boxed -> {Back, Bottom, Left}
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Solution3D"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , MeshStyle -> BoundaryTracingStyle["Edge3D"]
      , PlotPoints -> 20
      , PlotRange -> {0, wallHeight}
      , PlotRangePadding -> {{Scaled[0.125], Automatic}, Automatic, None}
      , PlotStyle -> BoundaryTracingStyle["Solution3D"]
      , TicksStyle -> LabelSize["Tick"]
      , RegionFunction -> Function[{x, y},
          RPolar[x, y] < rMax
        ]
      , ViewPoint -> {3.8, -1.2, 1.7}
    ],
    (* Wedge walls *)
    Plot3D[
      wallHeight, {x, xArcCorner, 0}, {y, -yArcCorner, yArcCorner}
      , BoundaryStyle -> BoundaryTracingStyle["Edge3D"]
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Wall3D"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , Mesh -> None
      , PlotPoints -> 20
      , PlotStyle -> BoundaryTracingStyle["Wall3D"]
      , RegionFunction -> Function[{x, y},
          Abs[y] < x Tan[alpha]
        ]
    ],
    (* Manually drawn solution base edges *)
    ParametricPlot3D[
      rMax {Cos[phi], Sin[phi], 0}
      , {phi, -alpha, alpha}
      , PlotStyle -> Directive[Black, BoundaryTracingStyle["Edge3D"]]
    ],
    (* Manually drawn wedge wall vertical edges *)
    Graphics3D @ {
      BoundaryTracingStyle["Edge3D"],
      verticalEdge /@ {
        {xArcCorner, -yArcCorner},
        {xArcCorner, +yArcCorner},
        {0, 0}
      },
      {}
    },
    {}
    , ImageSize -> 0.6 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Figure: mesh with additional refinement (wedge_obtuse-mesh-additional-refinement)*)


Module[
  {
    globalCoarseLengthScale, wallFineLengthScale, nearCornerRadius,
    apd, nearCornerFineLengthScale,
      rMax, alpha,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      equilateralTriangleArea, refinementFunction,
      boundaryMesh, mesh, numMeshElements,
      predicateWet,
      tNumerical, tXDerivative, tYDerivative, tSlope,
      rValues,
        symmetryHeightData, symmetrySlopeData,
        wallHeightData, wallSlopeData,
    rMaxPlot, xMinPlot, yMaxPlot, xMaxPlot,
    rTickValues, tickLength,
    tickTextStyle, axisTextStyle,
    dummyForTrailingCommas
  },
  (* Constants *)
  globalCoarseLengthScale = 0.2;
  wallFineLengthScale = 0.01;
  nearCornerRadius = 0.01;
  (* Constants (2) *)
  apd = 135;
  nearCornerFineLengthScale = 0.001;
  (* Build mesh *)
    (* Geometry *)
    rMax = rMaxMesh;
    alpha = apd * Degree;
    (* Boundary points and elements *)
    boundaryPoints =
      Join[
        (* Lower wedge wall *)
        Table[
          XYPolar[r, -alpha]
          , {r, UniformRange[0, rMax, wallFineLengthScale]}
        ]
          // Most,
        (* Far arc at infinity *)
        Table[
          XYPolar[rMax, phi]
          , {phi, UniformRange[-alpha, alpha, globalCoarseLengthScale / rMax]}
        ]
          // Most,
        (* Upper wedge wall *)
        Table[
          XYPolar[r, +alpha]
          , {r, UniformRange[rMax, 0, -wallFineLengthScale]}
        ]
          // Most,
        {}
      ];
    numBoundaryPoints = Length[boundaryPoints];
    boundaryElements =
      LineElement /@ {
        Table[{n, n + 1}, {n, numBoundaryPoints}]
          // Mod[#, numBoundaryPoints, 1] &
      };
    (* Special refinement *)
    equilateralTriangleArea = Function[{len}, Sqrt[3]/4 * len^2];
    refinementFunction =
      Function[{vertices, area},
        Or[
          area > equilateralTriangleArea[globalCoarseLengthScale],
          Norm @ Mean[vertices] < nearCornerRadius
            && area > equilateralTriangleArea[nearCornerFineLengthScale]
        ]
      ];
    (* Build mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> boundaryPoints,
        "BoundaryElements" -> boundaryElements,
        {}
      ];
    mesh = ToElementMesh[boundaryMesh
      , "ImproveBoundaryPosition" -> True
      , MeshRefinementFunction -> refinementFunction
    ];
  (* Plot *)
  rMaxPlot = 0.055;
  {xMinPlot, yMaxPlot} = XYPolar[rMaxPlot, alpha];
  xMaxPlot = Abs[xMinPlot];
  (* Make plot *)
  tickTextStyle = Style[#, LabelSize["Tick"]] & @* LaTeXStyle;
  axisTextStyle = Style[#, LabelSize["Axis"]] & @* LaTeXStyle;
  Show[
    mesh["Wireframe"],
    (* Ticks *)
    rTickValues = Range[0, rMaxPlot, 0.02] // Rest;
    tickLength = 0.03 rMaxPlot;
    Graphics @ {
      Table[
        {
          Line @ {
            XYPolar[r, alpha],
            XYPolar[r, alpha] + XYPolar[tickLength, alpha + Pi/2]
          },
          Text[
            r // tickTextStyle
            , XYPolar[r, alpha] + XYPolar[1 tickLength, alpha + Pi/2]
            , {1.2, 0.65}
          ],
          {}
        }
        , {r, rTickValues}
      ]
    },
    {}
    , ImageSize -> 0.45 ImageSizeTextWidth
    , PlotRange -> {{xMinPlot, xMaxPlot}, {-yMaxPlot, yMaxPlot}}
  ]
] // Ex["wedge_obtuse-mesh-detail-additional-refinement.pdf"]


(* ::Section:: *)
(*Figure: near-corner slope with additional refinement (wedge_obtuse-moderate-slope-additional-refinement-*)*)


Module[
  {
    rMax, casesWithinRange,
    allData,
    slopeMin, slopeMax,
    apd, nearCornerFineLengthScaleValues,
    symmetrySlopeData, wallSlopeData,
    opts,
    symmetryPlot, legend, wallPlot,
    dummyForTrailingCommas
  },
  (* Horizontal plot range *)
  rMax = 0.015;
  casesWithinRange = Cases[{r_, _} /; r <= rMax];
  (* Import all data *)
  allData = Import["slope_test/wedge_obtuse-slope_test-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apd = 135;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd, ncfls_, _, ssData_, _, wsData_} :> (
      symmetrySlopeData[ncfls] = ssData // casesWithinRange;
      wallSlopeData[ncfls] = wsData // casesWithinRange;
    )
  ];
  (* Vertical plot range *)
  {slopeMin, slopeMax} = MinMax[
    Table[
      {
        symmetrySlopeData[ncfls][[All, 2]],
        wallSlopeData[ncfls][[All, 2]]
      }
      , {ncfls, nearCornerFineLengthScaleValues}
    ]
  ];
  (* Plot options *)
  opts = {
    AspectRatio -> 0.5,
    AxesLabel -> {Italicise["r"], "Slope"},
    ImageSize -> 0.4 ImageSizeTextWidth,
    Joined -> True,
    LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"],
    TicksStyle -> LabelSize["Tick"],
    PlotMarkers -> {{Automatic, 6}, {"OpenMarkers", 5}},
    PlotRange -> {slopeMin, slopeMax},
    PlotRangeClipping -> False,
    PlotRangePadding -> {{None, 0.001}, {0.05, 0.08}},
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  };
  (* Slope along line of symmetry *)
  (*
    First use PlotLegends to get a properly styled legend,
    then extract it for separate export.
  *)
  symmetryPlot =
    ListPlot[
      Table[symmetrySlopeData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , PlotLegends -> Placed[
          LineLegend[
            nearCornerFineLengthScaleValues
            , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
            (* Don't bother with LegendLabel; just write \ell_\ultra in LaTeX *)
            , LegendLayout -> "Row"
            , LegendMarkerSize -> {4, (* magic *) 0.6} LabelSize["Legend"]
          ]
          , Scaled @ {0.5, 0.5}
        ]
      , opts
    ];
  symmetryPlot /. {
    Legended[p_, Placed[l_, ___]] :> (
      symmetryPlot = p;
      legend = l;
    )
  };
  legend = GraphicsRow[{legend}
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.1
    , Spacings -> -ImageSizeTextWidth (* (no idea why this works) *)
  ];
  (* Slope along wall *)
  wallPlot =
    ListPlot[
      Table[wallSlopeData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Export *)
  {
    symmetryPlot // Ex["wedge_obtuse-moderate-slope-additional-refinement-symmetry.pdf"],
    wallPlot // Ex["wedge_obtuse-moderate-slope-additional-refinement-wall.pdf"],
    legend // Ex["wedge_obtuse-moderate-slope-additional-refinement-legend.pdf"],
    Nothing
  }
]


(* ::Section:: *)
(*Figure: near-corner height with additional refinement (wedge_obtuse-moderate-height-additional-refinement-*)*)


Module[
  {
    rMax, casesWithinRange,
    allData,
    heightMin, heightMax,
    apd, nearCornerFineLengthScaleValues,
    symmetryHeightData, wallHeightData,
    opts,
    symmetryPlot, wallPlot,
    dummyForTrailingCommas
  },
  (* Horizontal plot range *)
  rMax = 0.015;
  casesWithinRange = Cases[{r_, _} /; r <= rMax];
  (* Import all data *)
  allData = Import["slope_test/wedge_obtuse-slope_test-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apd = 135;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd, ncfls_, shData_, _, whData_, _} :> (
      symmetryHeightData[ncfls] = shData // casesWithinRange;
      wallHeightData[ncfls] = whData // casesWithinRange;
    )
  ];
  (* Vertical plot range *)
  heightMin = 0;
  heightMax = 1.2 Max[
    Table[
      {
        symmetryHeightData[ncfls][[All, 2]],
        wallHeightData[ncfls][[All, 2]]
      }
      , {ncfls, nearCornerFineLengthScaleValues}
    ]
  ];
  (* Plot options *)
  opts = {
    AspectRatio -> 0.5,
    AxesLabel -> {Italicise["r"], "Height"},
    ImageSize -> 0.4 ImageSizeTextWidth,
    Joined -> True,
    LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"],
    TicksStyle -> LabelSize["Tick"],
    PlotMarkers -> {{Automatic, 6}, {"OpenMarkers", 5}},
    PlotRange -> {heightMin, heightMax},
    PlotRangeClipping -> False,
    PlotRangePadding -> {{None, 0.001}, None},
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  };
  (* Slope along line of symmetry *)
  symmetryPlot =
    ListPlot[
      Table[symmetryHeightData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Slope along wall *)
  wallPlot =
    ListPlot[
      Table[wallHeightData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Export *)
  {
    symmetryPlot // Ex["wedge_obtuse-moderate-height-additional-refinement-symmetry.pdf"],
    wallPlot // Ex["wedge_obtuse-moderate-height-additional-refinement-wall.pdf"],
    Nothing
  }
]


(* ::Section:: *)
(*Figure: near-corner slope with additional refinement, large angle case (wedge_obtuse-large-slope-additional-refinement-*)*)


Module[
  {
    rMax, casesWithinRange,
    allData,
    slopeMin, slopeMax,
    apd, nearCornerFineLengthScaleValues,
    symmetrySlopeData, wallSlopeData,
    opts,
    symmetryPlot, wallPlot,
    dummyForTrailingCommas
  },
  (* Horizontal plot range *)
  rMax = 0.015;
  casesWithinRange = Cases[{r_, _} /; r <= rMax];
  (* Import all data *)
  allData = Import["slope_test_large/wedge_obtuse-slope_test_large-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apd = 135;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd, ncfls_, _, ssData_, _, wsData_} :> (
      symmetrySlopeData[ncfls] = ssData // casesWithinRange;
      wallSlopeData[ncfls] = wsData // casesWithinRange;
    )
  ];
  (* Vertical plot range *)
  slopeMin = 0;
  slopeMax = Max[
    Table[
      {
        symmetrySlopeData[ncfls][[All, 2]],
        wallSlopeData[ncfls][[All, 2]]
      }
      , {ncfls, nearCornerFineLengthScaleValues}
    ]
  ];
  (* Plot options *)
  opts = {
    AspectRatio -> 0.5,
    AxesLabel -> {Italicise["r"], "Slope"},
    ImageSize -> 0.4 ImageSizeTextWidth,
    Joined -> True,
    LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"],
    TicksStyle -> LabelSize["Tick"],
    PlotMarkers -> {{Automatic, 6}, {"OpenMarkers", 5}},
    PlotRange -> {slopeMin, slopeMax},
    PlotRangeClipping -> False,
    PlotRangePadding -> {{None, 0.001}, {None, Scaled[0.15]}},
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  };
  (* Slope along line of symmetry *)
  (*
    Don't bother with legend;
    just re-use "wedge_obtuse-moderate-slope-additional-refinement-legend.pdf".
  *)
  symmetryPlot =
    ListPlot[
      Table[symmetrySlopeData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Slope along wall *)
  wallPlot =
    ListPlot[
      Table[wallSlopeData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Export *)
  {
    symmetryPlot // Ex["wedge_obtuse-large-slope-additional-refinement-symmetry.pdf"],
    wallPlot // Ex["wedge_obtuse-large-slope-additional-refinement-wall.pdf"],
    Nothing
  }
]


(* ::Section:: *)
(*Figure: near-corner height with additional refinement, large angle case (wedge_obtuse-large-height-additional-refinement-*)*)


Module[
  {
    rMax, casesWithinRange,
    allData,
    heightMin, heightMax,
    apd, nearCornerFineLengthScaleValues,
    symmetryHeightData, wallHeightData,
    opts,
    symmetryPlot, wallPlot,
    dummyForTrailingCommas
  },
  (* Horizontal plot range *)
  rMax = 0.015;
  casesWithinRange = Cases[{r_, _} /; r <= rMax];
  (* Import all data *)
  allData = Import["slope_test_large/wedge_obtuse-slope_test_large-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apd = 135;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd, ncfls_, shData_, _, whData_, _} :> (
      symmetryHeightData[ncfls] = shData // casesWithinRange;
      wallHeightData[ncfls] = whData // casesWithinRange;
    )
  ];
  (* Vertical plot range *)
  heightMin = 0;
  heightMax = 1.25 Max[
    Table[
      {
        symmetryHeightData[ncfls][[All, 2]],
        wallHeightData[ncfls][[All, 2]]
      }
      , {ncfls, nearCornerFineLengthScaleValues}
    ]
  ];
  (* Plot options *)
  opts = {
    AspectRatio -> 0.5,
    AxesLabel -> {Italicise["r"], "Height"},
    ImageSize -> 0.4 ImageSizeTextWidth,
    Joined -> True,
    LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"],
    TicksStyle -> LabelSize["Tick"],
    PlotMarkers -> {{Automatic, 6}, {"OpenMarkers", 5}},
    PlotRange -> {heightMin, heightMax},
    PlotRangeClipping -> False,
    PlotRangePadding -> {{None, 0.001}, None},
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  };
  (* Slope along line of symmetry *)
  symmetryPlot =
    ListPlot[
      Table[symmetryHeightData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Slope along wall *)
  wallPlot =
    ListPlot[
      Table[wallHeightData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Export *)
  {
    symmetryPlot // Ex["wedge_obtuse-large-height-additional-refinement-symmetry.pdf"],
    wallPlot // Ex["wedge_obtuse-large-height-additional-refinement-wall.pdf"],
    Nothing
  }
]


(* ::Section:: *)
(*Figure: viable domain  (wedge_obtuse-viable)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.1 xCritical, 0.05];
  xMin = Round[-3 xMax, 0.05];
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
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
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 3
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
    , ImageSize -> 0.5 * 0.7 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-viable.pdf"]


(* ::Section:: *)
(*Figure: generic traced boundaries  (wedge_obtuse-traced-boundaries)*)


(* (This is a little slow (~5 sec).) *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    xyStartListWall, xyStartListSymmetry, xyStartList,
    sMax, xyTracedList,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.1 xCritical, 0.05];
  xMin = Round[-3 xMax, 0.05];
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary starting points (upper branch) *)
  xyStartListWall = Table[XYPolar[r, alpha], {r, Subdivide[rMaxWall, 12]}];
  xyStartListSymmetry = Table[{x, 0}, {x, {0.27, 0.52, 0.78} xCritical}];
  xyStartList = Join[xyStartListWall, xyStartListSymmetry];
  (* Traced boundaries (upper branch) *)
  sMax = RPolar[xMax - xMin, 2 yMax];
  xyTracedList =
    Table[
      Quiet[
        ContactTracedBoundary[derList][xyStart, 0, {-1, 1} sMax
          , -1, 1
          , 0
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ]
      , {xyStart, xyStartList}
    ];
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
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
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
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
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedList}
    ],
    {}
    , ImageSize -> 0.5 * 0.7 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: viable domain, larger scale  (wedge_obtuse-viable-larger-scale)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = 2 xCritical;
  xMin = -7 xMax;
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
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
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 7
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
    , ImageSize -> 0.36 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-viable-larger-scale.pdf"]


(* ::Section:: *)
(*Figure: terminal points  (wedge_obtuse-terminal-points)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    xOrdinary, yOrdinary,
    sMax, xyContourList,
    textStyle, textStyleBracket,
    plot,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = 3.1 xCritical;
  xMin = -0.65 xMax;
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Ordinary terminal point *)
  xOrdinary = Way[0, xCritical, 0.05];
  yOrdinary = SeekRoot[vi[xOrdinary, #] &, {0, yMax}, 5];
  (* T-contours *)
  sMax = RPolar[xMax - xMin, 2 yMax];
  xyContourList =
    Table[
      ContourByArcLength[tNumerical][xyInitial, 0, sMax {-1, 1}, 1]
      , {xyInitial, {{xCritical, 0}, {xOrdinary, yOrdinary}}}
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Make plot *)
  plot = Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , Frame -> None
    ],
    (* Wedge walls *)
    (*Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },*)
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> None
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
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
      , {x, xMinMore, xMaxMore}
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
        {-1.7, -0.82}
      ] // textStyle,
      Text[
        "critical",
        {xCritical, 0},
        {-1.55, 0.82}
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
        {-1.4, -0.25}
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
    , ItemAspectRatio -> 2 yMax / (xMax - xMin)
    , ImageSize -> 0.63 ImageSizeTextWidth
    , Spacings -> {0, 0}
  ]
] // Ex["wedge_obtuse-terminal-points.pdf"]


(* ::Section:: *)
(*Figure: generic traced boundary branches  (wedge_obtuse-traced-boundaries-*-branch)*)


(* (This is a little slow (~5 sec).) *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    rStartListWall, xyStartListWall, xyStartListSymmetry, xyStartList,
    sMax, xyTracedList,
    branchCases, branchYSigns, ySign,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = 1.2 xCritical;
  xMin = -5 xCritical;
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary starting points (upper branch) *)
  rStartListWall = (1/12 + Subdivide[9]) rMaxWall // Prepend[0];
  xyStartListWall = Table[XYPolar[r, alpha], {r, rStartListWall}];
  xyStartListSymmetry = Table[{x, 0}, {x, {0.2, 0.35, 0.5, 0.7, 0.9} xCritical}];
  xyStartList = Join[xyStartListWall, xyStartListSymmetry];
  (* Traced boundaries (upper branch) *)
  sMax = RPolar[xMax - xMin, 2 yMax];
  xyTracedList =
    Table[
      Quiet[
        ContactTracedBoundary[derList][xyStart, 0, {-1, 1} sMax
          , -1, 1
          , 0, Function[{x, y}, x < xMinMore]
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ]
      , {xyStart, xyStartList}
    ];
  (* Make plot *)
  branchCases = {"upper", "lower"};
  branchYSigns = {1, -1};
  Table[
    ySign = AssociationThread[branchCases -> branchYSigns][case];
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax}
        , Frame -> None
      ],
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          xMinMore {1, ySign * Tan[alpha]},
          {0, 0}
        }
      },
      (* Non-viable domain *)
      RegionPlot[
        vi[x, y] < 0
        , {x, xMinMore, xMaxMore}
        , {y, -yMaxMore, yMaxMore}
        , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
        , PlotPoints -> 7
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
          /. line_Line :> {Arrowheads @ {{-0.1, 0.33}}, Arrow[line]}
        , {xy, xyTracedList}
      ],
      {}
      , ImageSize -> 0.45 * 0.6 ImageSizeTextWidth
    ] // Ex @ FString["wedge_obtuse-traced-boundaries-{case}-branch.pdf"]
    , {case, branchCases}
  ]
]


(* ::Section:: *)
(*Figure: traced boundaries, different contact angle (wedge_obtuse-traced-boundaries-different-angle-*)*)


(* (This is a little slow (~6 sec).) *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracingValues, caseNameList,
    gammaTracing,
    d,
    xCritical,
    xTerminal, yTerminal, rTerminal,
    xMin, xMax, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    clearance, xyStartList,
    sMax, xyTracedList,
    plot,
    caseName,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Different tracing contact angle *)
  gpdTracingValues = {58, 63};
  caseNameList = {"less", "more"};
  Table[
    gammaTracing = gpdTracing * Degree;
    (* Half-plane offset distance *)
    d = DHalfPlane[gamma, gammaTracing];
    (* Critical terminal point x_0 *)
    xCritical = x0[tNumerical, gammaTracing];
    (* Derivative list for boundary tracing *)
    {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
    (* Plot range *)
    If[
      gammaTracing < gamma,
        (* Intersection of terminal curve and wall *)
        xTerminal = SeekRoot[vi[#, # Tan[alpha]] &, {0, -10 xCritical}];
        yTerminal = xTerminal Tan[alpha];
        rTerminal = RPolar[xTerminal, yTerminal];
        (* Plot range *)
        xMin = 1.15 xTerminal;
        xMax = 1.5 xCritical;
        yMax = xMin Tan[alpha];
      ,
        xMax = 1.7 xCritical;
        xMin = -4.5 xCritical;
        yMax = xMin Tan[alpha]
        (*yMax = 0.65 SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20]*);
    ];
    rMaxWall = xMin Sec[alpha];
    (* Plot range but more *)
    more = 0.15;
    xMinMore = xMin (1 + more);
    xMaxMore = xMax (1 + more);
    yMaxMore = yMax (1 + more);
    (* Starting points *)
    clearance = 10^-6;
    If[
      gammaTracing < gamma,
        xyStartList = Join[
          Table[XYPolar[r, +alpha], {r, {0, 0.09, 0.21, 0.34, 0.48, 0.63} rTerminal}],
          Table[XYPolar[r, -alpha + clearance], {r, {0.31, 0.65} rTerminal}],
          {{0.5 xCritical, 0}},
          {}
        ];
      ,
        xyStartList = Join[
          {{0, 0}},
          Table[XYPolar[r, +alpha], {r, (Range[6] 0.16 - 0.02) rMaxWall}],
          Table[XYPolar[r, -alpha + clearance], {r, (Range[3, 9, 3] 0.16 - 0.02) rMaxWall}],
          {XYPolar[2 rMaxWall, -alpha] + XYPolar[d, -alpha + Pi/2]},
          Table[
            {
              xStart,
              SeekRoot[vi[xStart, #] &, {-xStart Tan[alpha], -yMax}, 10] + clearance
            }
            , {xStart, {-0.2, -1.25} xCritical}
          ],
          {}
        ];
    ];
    (* Traced boundaries (upper branch) *)
    sMax = rMaxMesh;
    xyTracedList =
      Table[
        Quiet[
          ContactTracedBoundary[derList][xyStart, 0, {-1, 1} sMax
            , -1, 1
            , 0, Function[{x, y}, x > xMaxMore]
          ]
          , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
        ]
        , {xyStart, xyStartList}
      ];
    (* Plot *)
    plot =
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax}
          , FrameLabel -> {
              Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
              Italicise["y"]
            }
          , FrameTicksStyle -> LabelSize["Tick"]
          , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
        ],
        (* Non-viable domain *)
        RegionPlot[
          vi[x, y] < 0
          , {x, xMinMore, xMaxMore}
          , {y, -yMaxMore, yMaxMore}
          , BoundaryStyle -> None
          , PlotPoints -> 5
          , PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        ContourPlot[
          vi[x, y] == 0
          , {x, xMinMore, xMaxMore}
          , {y, -yMaxMore, yMaxMore}
          , ContourLabels -> None
          , ContourStyle -> BoundaryTracingStyle["Terminal"]
          , PlotPoints -> 5
        ],
        (* Wedge walls *)
        Graphics @ {BoundaryTracingStyle["Wall"],
          Line @ {
            xMinMore {1, Tan[alpha]},
            {0, 0},
            xMinMore {1, -Tan[alpha]}
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
      , ImageSize -> {Automatic, 1.4 * 0.45 ImageSizeTextWidth}
    ] // Ex @ FString["wedge_obtuse-traced-boundaries-different-angle-{caseName}.pdf"]
    , {gpdTracing, gpdTracingValues}
  ]
]


(* ::Section:: *)
(*Figure: traced boundary branches, different contact angle (wedge_obtuse-traced-boundaries-different-angle-*-branch)*)


(* (This is a little slow (~5 sec).) *)
Module[
  {
    apd, gpd, gpdTracing,
    alpha, gamma, gammaTracing,
    d,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    clearance, xyStartList,
    sMax, xyTracedList,
    branchCases, branchYSigns, ySign,
    nonManifoldStyle,
    arrowify,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd, gpdTracing} = {135, 60, 63};
  {alpha, gamma, gammaTracing} = {apd, gpd, gpdTracing} Degree;
  (* Half-plane offset distance *)
  d = DHalfPlane[gamma, gammaTracing];
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gammaTracing];
  (* Plot range *)
  xMax = 1.2 xCritical;
  xMin = -5 xCritical;
  yMax = xMin Tan[alpha];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary starting points (upper branch) *)
  clearance = 10^-6;
  xyStartList["general"] = Join[
    (* collides with vertex *)
    {{0, 0}},
    (* collides with upper wall *)
    Table[XYPolar[r, +alpha], {r, (Range[5] 0.16 - 0.02) rMaxWall}],
    {}
  ];
  xyStartList["lower"] =
    (* deflects into lower wall *)
    Table[XYPolar[r, -alpha + clearance], {r, (Range[2, 8, 2] 0.16 - 0.02) rMaxWall}];
  xyStartList["terminal"] =
    (* deflects into lower half of terminal curve *)
    Table[
      {
        xStart,
        SeekRoot[vi[xStart, #] &, {-xStart Tan[alpha], -yMax}, 10] + clearance
      }
      , {xStart, {-0.2, -1.25} xCritical}
    ];
  xyStartList["manifold"] =
    (* close enough to unstable manifold *)
    {XYPolar[2 rMaxWall, -alpha] + XYPolar[d, -alpha + Pi/2]};
  (* Traced boundaries (upper branch) *)
  sMax = rMaxMesh;
  Table[
    xyTracedList[type] =
      Table[
        Quiet[
          ContactTracedBoundary[derList][xyStart, 0, {-1, 1} sMax
            , -1, 1
            , 0, Function[{x, y}, x < xMinMore]
          ]
          , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
        ]
        , {xyStart, xyStartList[type]}
      ]
    , {type, {"general", "lower", "terminal", "manifold"}}
  ];
  (* Make plot *)
  branchCases = {"upper", "lower"};
  branchYSigns = {1, -1};
  arrowify[proportion_] :=
    # /. {line_Line :> {Arrowheads @ {{-0.09, proportion}}, Arrow[line]}} &;
  nonManifoldStyle = Directive[
    BoundaryTracingStyle["Traced"],
    BoundaryTracingStyle["Background"] // Last
  ];
  Table[
    ySign = AssociationThread[branchCases -> branchYSigns][case];
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax}
        , Frame -> None
        , PlotRangePadding -> None
      ],
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          xMinMore {1, Tan[alpha]},
          {0, 0},
          xMinMore {1, -Tan[alpha]}
        }
      },
      (* Non-viable domain *)
      RegionPlot[
        vi[x, ySign * y] < 0
        , {x, xMinMore, xMaxMore}
        , {y, -yMaxMore, yMaxMore}
        , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
        , PlotPoints -> 7
        , PlotStyle -> BoundaryTracingStyle["NonViable"]
      ],
      (* Traced boundaries (general) *)
      Table[
        ParametricPlot[
          EvaluatePair[xy, s, ySign] // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> nonManifoldStyle
        ]
          // arrowify[0.3]
        , {xy, xyTracedList["general"]}
      ],
      (* Traced boundaries (deflects into lower wall) *)
      Table[
        ParametricPlot[
          EvaluatePair[xy, s, ySign] // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> nonManifoldStyle
        ]
          // arrowify[0.45]
        , {xy, xyTracedList["lower"]}
      ],
      (* Traced boundaries (deflects into terminal) *)
      Table[
        ParametricPlot[
          EvaluatePair[xy, s, ySign] // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> nonManifoldStyle
        ]
          // arrowify[0.15]
        , {xy, xyTracedList["terminal"]}
      ],
      (* Traced boundaries (unstable manifold) *)
      Table[
        ParametricPlot[
          EvaluatePair[xy, s, ySign] // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
          // arrowify[0.56]
        , {xy, xyTracedList["manifold"]}
      ],
      {}
      , ImageSize -> 0.45 * 0.75 ImageSizeTextWidth
    ] // Ex @ FString["wedge_obtuse-traced-boundaries-different-angle-{case}-branch.pdf"]
    , {case, branchCases}
  ]
]


(* ::Section:: *)
(*Figure: pseudo-rounding construction (wedge_obtuse-pseudo-rounding-construction-*)*)


Module[
  {
    apd, gpd, gpdTracing,
    alpha, gamma, gammaTracing,
    d,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xArrowBase, yArrowBase,
    xInitial, yInitial,
    xMinWall, yMaxWall,
    xMinOffset, yMaxOffset,
    arrowLength,
    orthogonalityMarkerLength, orthogonalityMarkerStyle,
    wallThicknessCorrection,
    textStyleLabel, textStylePoint, textStylePointBracket,
    xyTraced, sCorner, xyCorner,
    plot, opts,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd, gpdTracing} = {135, 60, 63};
  {alpha, gamma, gammaTracing} = {apd, gpd, gpdTracing} Degree;
  (* Half-plane offset distance *)
  d = DHalfPlane[gamma, gammaTracing];
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gammaTracing];
  (* Base of wall coordinate (xi) arrow *)
  xArrowBase = -7 xCritical;
  yArrowBase = xArrowBase Tan[alpha];
  (* End point of traced boundary (used for initial condition) *)
  {xInitial, yInitial} = {xArrowBase, yArrowBase} + XYPolar[d, alpha - Pi/2];
  (* End point of wall to be shown *)
  {xMinWall, yMaxWall} = 1.2 {xArrowBase, yArrowBase};
  (* End point of offset wedge to be shown *)
  {xMinOffset, yMaxOffset} = {xMinWall, yMaxWall} + XYPolar[d, alpha - Pi/2];
  (* Etc. *)
  arrowLength = 2.7 d;
  orthogonalityMarkerLength = 0.35 d;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[]];
  wallThicknessCorrection = -0.1;
  (* Text style *)
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStylePoint = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStylePointBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Compute traced boundary (LOWER branch, which is physically higher) *)
  xyTraced =
    Quiet[
      ContactTracedBoundary[derList][
        {xInitial, yInitial}, 0, {0, 2 yMaxWall}
        , 1, -1
      ]
      , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
    ];
  (* Corner point (y == 0) *)
  sCorner = SeekRoot[xyTraced[[2]], {DomainStart[xyTraced], DomainEnd[xyTraced]}];
  xyCorner = EvaluatePair[xyTraced, sCorner];
  (* Wedge walls *)
  plot["walls"] = Show[
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        {xMinWall, +yMaxWall},
        {0, 0},
        {xMinWall, -yMaxWall}
      }
    },
    {}
  ];
  (* Plot showing direction *)
  plot["direction"] = Show[
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
        , {-2.25, -0.3}
      ],
      {}
    },
    Graphics @ {orthogonalityMarkerStyle,
      Rectangle[
        {xArrowBase, yArrowBase},
        {xArrowBase, yArrowBase}
          + orthogonalityMarkerLength {1 + wallThicknessCorrection, -1}
      ] // Rotate[#, alpha - Pi/2, {xArrowBase, yArrowBase}] &
    },
    (* Wedge walls *)
    plot["walls"],
    Graphics @ {
      Text[
        "\[Xi]" == SeparatedRow[]["\[VeryThinSpace]", 0] // textStyleLabel
        , {xMinWall, +yMaxWall}
        , {1.4, -0.35}
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
        , {xMinOffset, +yMaxOffset}
        , {0.5, -1.2}
      ]
    },
    (* Offset wedge *)
    (* (separate lines to ensure dotting at the corner) *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {
        {d / Sin[alpha], 0},
        {xMinOffset, +yMaxOffset}
      },
      Line @ {
        {d / Sin[alpha], 0},
        {xMinOffset, -yMaxOffset} * 1.02 (* (optical correction) *)
      },
      {}
    },
    (* Traced boundaries *)
    ParametricPlot[
      EvaluatePair[xyTraced, s] // IncludeYReflection // Evaluate
      , {s, DomainStart[xyTraced], DomainEnd[xyTraced]}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ] /. {line_Line :> {Arrowheads @ {{0.1, 0.85}}, Arrow[line]}},
    {}
  ];
  (* Plot showing psuedo-rounding corner *)
  plot["corner"] = Show[
    (* Wedge walls *)
    plot["walls"],
    (* Traced boundaries *)
    ParametricPlot[
      EvaluatePair[xyTraced,
        sCorner - Abs[s - sCorner]
        , -Sign[s - sCorner]
      ] // Evaluate
      , {s, 0, 2 sCorner}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    (* Corner point (x_c, 0) *)
    Graphics @ {
      Text[
        Row @ {
          "(" // textStylePointBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], "\[VeryThinSpace]\[VeryThinSpace]c"],
          ",\[ThinSpace]",
          0,
          ")" // textStylePointBracket
        } // textStylePoint
        , xyCorner
        , {-1.25, -0.15}
      ],
    },
    {}
  ];
  (* Options (fix vertical sizing) *)
  opts = {
    ImageSize -> {Automatic, 0.45 * 1 ImageSizeTextWidth},
    PlotRange -> {{All, All}, {-1.2 yMaxWall, 1.35 yMaxWall}},
    {}
  };
  (* Export *)
  Table[
    Show[plot[case], opts]
      // Ex @ FString["wedge_obtuse-pseudo-rounding-construction-{case}.pdf"]
    , {case, {"direction", "corner"}}
  ]
]


(* ::Section:: *)
(*Table: pseudo-rounding corner angles (wedge_obtuse-pseudo-rounding-corner-angle-table)*)


(* (This is slow (~6 sec).) *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracingValues,
    xCriticalMin, xCriticalMax,
    xMin, xMax, yMax,
    more, xMinMore,
    gammaTracing, d,
    xInitial, yInitial,
    p, q, grad2, f, vi, derList,
    xyTraced, sCorner,
    xDerivative, yDerivative,
    cornerAngle, cornerAngleDeviation,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Tracing contact angles *)
  gpdTracingValues = Range[gpd + 5, 85, 5];
  (* Plot range *)
  xCriticalMin = x0[tNumerical, Min[gpdTracingValues] Degree];
  xCriticalMax = x0[tNumerical, Max[gpdTracingValues] Degree];
  xMin = -2 xCriticalMax;
  xMax = 1.5 xCriticalMax;
  yMax = 1.2 xMin Tan[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  (* Compute pseudo-roundings *)
  Table[
    (* Tracing contact angle *)
    gammaTracing = gpdTracing * Degree;
    (* Half-plane offset distance *)
    d = DHalfPlane[gamma, gammaTracing];
    (* Initial point *)
    {xInitial, yInitial} = 1.5 xMin {1, Tan[alpha]} + XYPolar[d, alpha - Pi/2];
    (* Derivative list for boundary tracing *)
    {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
    (* Compute traced boundary (LOWER branch, which be physically higher) *)
    xyTraced[gpdTracing] =
      Quiet[
        ContactTracedBoundary[derList][
          {xInitial, yInitial}, 0, {0, 2 yInitial}
          , 1, -1
          , -Infinity
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ];
    (* Corner point (y == 0) *)
    sCorner[gpdTracing] =
      SeekRoot[
        xyTraced[gpdTracing][[2]],
        {DomainStart @ xyTraced[gpdTracing], DomainEnd @ xyTraced[gpdTracing]}
        , 20
      ];
    (* Corner angle *)
    xDerivative = xyTraced[gpdTracing][[1]]' @ sCorner[gpdTracing];
    yDerivative = xyTraced[gpdTracing][[2]]' @ sCorner[gpdTracing];
    cornerAngle = 2 ArcTan[yDerivative / xDerivative // Abs];
    (* Deviation from a straight angle *)
    cornerAngleDeviation = cornerAngle - Pi;
    (* Return table row *)
    {
      gpdTracing,
      cornerAngle / Degree // SignificantFiguresForm[4],
      cornerAngleDeviation / Degree // DecimalPlacesForm[1],
      Nothing
    }
    , {gpdTracing, gpdTracingValues}
  ]
    // TableForm[#
      , TableHeadings -> {
          None,
          {"tracing angle / \[Degree]", "corner angle / \[Degree]", "deviation"}
        }
    ] &
    // Column @ {FString["contact angle: {gpd}\[Degree]"], #} &
] // Ex["wedge_obtuse-pseudo-rounding-corner-angle-table.pdf"]


(* ::Section:: *)
(*Figure: family of pseudo-roundings (wedge_obtuse-pseudo-roundings)*)


(* (This is slow (~6 sec).) *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    gpdTracingValues,
    xCriticalMin, xCriticalMax,
    xMin, xMax, yMax,
    more, xMinMore,
    gammaTracing, d,
    xInitial, yInitial,
    p, q, grad2, f, vi, derList,
    xyTraced, sCorner,
    xy, sMax,
    xArrowMin, xArrowMax,
    dummyForTrailingCommas
  },
  (* Known solution angles *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Tracing contact angles *)
  gpdTracingValues = Range[gpd + 5, 85, 5];
  (* Plot range *)
  xCriticalMin = x0[tNumerical, Min[gpdTracingValues] Degree];
  xCriticalMax = x0[tNumerical, Max[gpdTracingValues] Degree];
  xMin = -2 xCriticalMax;
  xMax = 1.5 xCriticalMax;
  yMax = 1.1 xMin Tan[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  (* Compute pseudo-roundings *)
  Table[
    (* Tracing contact angle *)
    gammaTracing = gpdTracing * Degree;
    (* Half-plane offset distance *)
    d = DHalfPlane[gamma, gammaTracing];
    (* Initial point *)
    {xInitial, yInitial} = 1.5 xMin {1, Tan[alpha]} + XYPolar[d, alpha - Pi/2];
    (* Derivative list for boundary tracing *)
    {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
    (* Compute traced boundary (LOWER branch, which be physically higher) *)
    xyTraced[gpdTracing] =
      Quiet[
        ContactTracedBoundary[derList][
          {xInitial, yInitial}, 0, {0, 2 yInitial}
          , 1, -1
          , -Infinity
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ];
    (* Corner point (y == 0) *)
    sCorner[gpdTracing] =
      SeekRoot[
        xyTraced[gpdTracing][[2]],
        {DomainStart @ xyTraced[gpdTracing], DomainEnd @ xyTraced[gpdTracing]}
        , 20
      ];
    , {gpdTracing, gpdTracingValues}
  ];
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -24}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Traced boundaries *)
    Table[
      xy = xyTraced[gpdTracing];
      sMax = sCorner[gpdTracing];
      ParametricPlot[
        EvaluatePair[xy, sMax - Abs[s - sMax], -Sign[s - sMax]]
          // Evaluate
        , {s, 0, 2 sMax}
        , PlotPoints -> 3
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {gpdTracing, gpdTracingValues}
    ],
    (* Tracing gamma arrow *)
    xArrowMin = Way[xCriticalMin, xCriticalMax, -0.06];
    xArrowMax = Way[xCriticalMin, xCriticalMax, +1.5];
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
        , {1, -1.2}
      ]
    },
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-pseudo-roundings.pdf"]


(* ::Section:: *)
(*Figure: pseudo-roundings (grouped by contact angle) with offset applied (wedge_obtuse-pseudo-roundings-offset)*)


(* (This is slow (~6 sec).) *)
Module[
  {
    apd, gpdTracing,
    alpha, gammaTracing,
    gpdValues,
    gamma, d,
    tNumerical,
    xCriticalOffset,
    xMin, xMax, yMax,
    more, xMinMore,
    p, q, grad2, f, vi, derList,
    xInitial, yInitial, xyTraced, sCorner,
    xy, sMax,
    xCriticalOffsetMinGamma, xCriticalOffsetMaxGamma,
    xArrowMin, xArrowMax,
    dummyForTrailingCommas
  },
  (* Prescribed angles *)
  apd = 135;
  gpdTracing = 60;
  {alpha, gammaTracing} = {apd, gpdTracing} Degree;
  (* Solution contact angles *)
  gpdValues = Range[10, 50, 10];
  (* For various solution contact angles (1): *)
  Table[
    (* Solution contact angle *)
    gamma = gpd * Degree;
    (* Half-plane offset distance *)
    d[gpd] = DHalfPlane[gamma, gammaTracing];
    (* Import numerical solution *)
    tNumerical[gpd] =
      Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    (* Offset critical terminal point (for plot range) *)
    xCriticalOffset[gpd] = x0[tNumerical[gpd], gammaTracing] - d[gpd] / Sin[alpha];
    , {gpd, gpdValues}
  ];
  (* Plot range *)
  xMin = 1.5 Min[xCriticalOffset /@ gpdValues];
  xMax = 0;
  yMax = xMin Tan[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  (* For various solution contact angles (2): *)
  Table[
    (* Initial point *)
    {xInitial, yInitial} = 1.5 xMin {1, Tan[alpha]} + XYPolar[d[gpd], alpha - Pi/2];
    (* Derivative list for boundary tracing *)
    {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical[gpd], gammaTracing];
    (* Compute traced boundary (LOWER branch, which be physically higher) *)
    xyTraced[gpd] =
      Quiet[
        ContactTracedBoundary[derList][
          {xInitial, yInitial}, 0, {0, 2 yInitial}
          , 1, -1
          , -Infinity
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ];
    (* Corner point (y == 0) *)
    sCorner[gpd] =
      SeekRoot[
        xyTraced[gpd][[2]],
        {DomainStart @ xyTraced[gpd], DomainEnd @ xyTraced[gpd]}
        , 20
      ];
    , {gpd, gpdValues}
  ];
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Superscript[
            Italicise["x"],
            Style["\[NegativeVeryThinSpace]\[Prime]", Magnification -> 1.3]
          ] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Traced boundaries *)
    Table[
      xy = xyTraced[gpd];
      sMax = sCorner[gpd];
      ParametricPlot[
        Subtract[
          EvaluatePair[xy, sMax - Abs[s - sMax], -Sign[s - sMax]],
          {d[gpd] / Sin[alpha], 0}
        ]
          // Evaluate
        , {s, 0, 2 sMax}
        , PlotPoints -> 3
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {gpd, gpdValues}
    ],
    (* Solution gamma arrow *)
    xCriticalOffsetMinGamma = xCriticalOffset[Min @ gpdValues];
    xCriticalOffsetMaxGamma = xCriticalOffset[Max @ gpdValues];
    xArrowMin = Way[xCriticalOffsetMinGamma, xCriticalOffsetMaxGamma, -0.3];
    xArrowMax = Way[xCriticalOffsetMinGamma, xCriticalOffsetMaxGamma, +1.4];
    Graphics @ {
      GeneralStyle["Translucent"], GeneralStyle["Thick"],
      Arrowheads[Medium],
      Arrow @ {{xArrowMin, 0}, {xArrowMax, 0}}
    },
    Graphics @ {
      Text[
        "\[Gamma]"
          // LaTeXStyle
          // Style[#, LabelSize["Label"]] &
        , {xArrowMin, 0}
        , {2.5, -0.2}
      ]
    },
    {}
    , ImageSize -> 0.35 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-pseudo-roundings-offset.pdf"]


(* ::Section:: *)
(*Figure: height-rise profiles (wedge_obtuse-height-rise-profiles)*)


(* (This is slow (~6 sec).) *)
Module[
  {
    apd, gpdTracing,
    alpha, gammaTracing, h,
    gpdValues,
    gamma, d,
    tNumerical,
    xCriticalOffset,
    xMin, xMax, yMax,
    sMaxPlot,
    more, xMinMore,
    p, q, grad2, f, vi, derList,
    xInitial, yInitial, xyTraced, sCorner,
    xy, sMax,
    xySharp, tCorner,
    dummyForTrailingCommas
  },
  (* Prescribed angles *)
  apd = 135;
  gpdTracing = 60;
  {alpha, gammaTracing} = {apd, gpdTracing} Degree;
  (* Half-plane wall height *)
  h = HHalfPlane[gammaTracing];
  (* Solution contact angles *)
  gpdValues = {10, 40};
  (* For various solution contact angles (1): *)
  Table[
    (* Solution contact angle *)
    gamma = gpd * Degree;
    (* Half-plane offset distance *)
    d[gpd] = DHalfPlane[gamma, gammaTracing];
    (* Import numerical solution *)
    tNumerical[gpd] =
      Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    (* Offset critical terminal point (for plot range) *)
    xCriticalOffset[gpd] = x0[tNumerical[gpd], gammaTracing] - d[gpd] / Sin[alpha];
    , {gpd, gpdValues // Append[gpdTracing]}
  ];
  (* Plot range *)
  xMin = 1.5 Min[xCriticalOffset /@ gpdValues];
  xMax = 0;
  yMax = xMin Tan[alpha];
  sMaxPlot = 4;
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  (* For various solution contact angles (2): *)
  Table[
    (* Initial point *)
    {xInitial, yInitial} = XYPolar[0.7 rMaxMesh, alpha] + XYPolar[d[gpd], alpha - Pi/2];
    (* Derivative list for boundary tracing *)
    {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical[gpd], gammaTracing];
    (* Compute traced boundary (LOWER branch, which be physically higher) *)
    xyTraced[gpd] =
      Quiet[
        ContactTracedBoundary[derList][
          {xInitial, yInitial}, 0, {0, rMaxMesh}
          , 1, -1
          , -Infinity
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ];
    (* Corner point (y == 0) *)
    sCorner[gpd] =
      SeekRoot[
        xyTraced[gpd][[2]],
        {DomainStart @ xyTraced[gpd], DomainEnd @ xyTraced[gpd]}
        , 20
      ];
    , {gpd, gpdValues}
  ];
  (* Sharp corner (original walls) *)
  xySharp[s_] := XYPolar[Abs[s], Sign[s] alpha];
  (* Sharp corner (dip) height *)
  tCorner = tNumerical[gpdTracing][0, 0];
  (* Make plot *)
  Plot[
    Join[
      (* Half-plane wall height *)
      h // List,
      (* Sharp corner *)
      tNumerical[gpdTracing] @@ xySharp[s] // List,
      (* Rounded corners *)
      Table[
        xy = xyTraced[gpd];
        sMax = sCorner[gpd];
        tNumerical[gpd] @@
          EvaluatePair[xy, sMax - Abs[s], -Sign[s]]
        , {gpd, gpdValues}
      ],
      {}
    ] // Evaluate
    , {s, -sMaxPlot, sMaxPlot}
    , AspectRatio -> 0.5
    , AxesLabel -> Italicise /@ {"s", "T"}
    (*, AxesOrigin -> {-sMaxPlot, Automatic}*)
    , ImageSize -> 0.6 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotPoints -> 4
    , PlotRange -> {Floor[tCorner, 0.1], All}
    (*, PlotRangePadding -> {None, Scaled /@ {0.01, 0.05}}*)
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
] // Ex["wedge_obtuse-height-rise-profiles.pdf"]


(* ::Section:: *)
(*Figure: approximating a contour with traced boundaries (wedge_obtuse-traced-boundaries-approximation-*)*)


(* (This is rather slow (~10 sec).) *)
Module[
  {
    apd, gpd, gpdTracing,
    alpha, gamma, gammaTracing,
    d, hTracing,
    tNumerical,
    xCritical,
    p, q, grad2, f, vi, derList,
    xMinWall, yMaxWall, rMaxWall,
    xMax, xMin, yMax,
    sMax,
    xStartContour, xyContour,
    xyStartListUpper, xyTracedListUpper,
    xyStartListLower, xyTracedListLower,
    intersectionCurveAss, sStart, sEnd,
    yMaxTraced, regionFunctionTraced, plot,
    yContourLabel, sContourLabel,
    dummyForTrailingCommas
  },
  (* Angular paramters *)
  {apd, gpd, gpdTracing} = {135, 20, 55};
  {alpha, gamma, gammaTracing} = {apd, gpd, gpdTracing} Degree;
  (* Half-plane offset distance *)
  d = DHalfPlane[gamma, gammaTracing];
  (* Half-plane wall height for tracing contact angle *)
  hTracing = HHalfPlane[gammaTracing];
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gammaTracing];
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gammaTracing];
  (* Wall extremities *)
  xMinWall = - xCritical;
  yMaxWall = xMinWall Tan[alpha];
  rMaxWall = RPolar[xMinWall, yMaxWall];
  (* Plot range *)
  xMax = 1.1 xCritical;
  xMin = xMinWall;
  yMax = 1.65 yMaxWall;
  (* Maxmimum arc length *)
  sMax = 4 rMaxWall;
  (* Contour T == h_tracing *)
  xStartContour = SeekRoot[tNumerical[#, 0] - hTracing &, {0, xCritical}, 20];
  xyContour = ContourByArcLength[tNumerical][{xStartContour, 0}, 0, {-1, 1} sMax];
  (* Traced boundaries (upper branch) *)
  xyStartListUpper = Join[
    Table[XYPolar[r, +alpha], {r, Range[0, 8] rMaxWall / 6}],
    Table[XYPolar[r, -alpha], {r, Range[0, 10] rMaxWall / 3 // Rest}],
    {}
  ] // SortBy[#, Last, Greater] &;
  xyTracedListUpper =
    Table[
      Quiet[
        ContactTracedBoundary[derList][xyStart, 0, {0, sMax}, -1, 1]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum, NDSolveValue::nrnum1}
      ]
      , {xyStart, xyStartListUpper}
    ];
  (* Traced boundaries (lower branch) *)
  xyStartListLower = {#[[1]], -#[[2]]} & /@ xyStartListUpper;
  xyTracedListLower =
    Table[
      Quiet[
        ContactTracedBoundary[derList][xyStart, 0, {0, sMax}, 1, -1]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum, NDSolveValue::nrnum1}
      ]
      , {xyStart, xyStartListLower}
    ];
  (* Determine intersections (hard-coded) *)
  (* (Too much work to try and automate this) *)
  intersectionCurveAss = Association[
    5 -> {15, 16},
    6 -> {14, 15},
    7 -> {12, 14},
    8 -> {11, 12},
    9 -> {10, 11},
    10 -> {9, 10},
    9 -> {10, 11},
    10 -> {9, 10},
    11 -> {8, 9},
    12 -> {7, 8},
    14 -> {6, 7},
    15 -> {5, 6},
    Nothing
  ];
  Table[
    sStart[n] = First @ SeekParametricIntersection[
      xyTracedListUpper[[n]],
      xyTracedListLower[[intersectionCurveAss[n] // First]]
    ];
    sEnd[n] = First @ SeekParametricIntersection[
      xyTracedListUpper[[n]],
      xyTracedListLower[[intersectionCurveAss[n] // Last]]
    ];
    , {n, Keys[intersectionCurveAss]}
  ];
  (* Make plots *)
  yMaxTraced = Way[yMaxWall, yMax, 0.42];
  regionFunctionTraced = Function[{x, y}, Abs[y] < yMaxTraced // Evaluate];
  plot["common"] = Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -17}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line[2 {{xMinWall, +yMaxWall}, {0, 0}, {xMinWall, -yMaxWall}}]
    },
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        EvaluatePair[xy, s]
          // IncludeYReflection
          // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Background"]
        , RegionFunction -> regionFunctionTraced
      ]
      , {xy, xyTracedListUpper}
    ],
    (* Contour *)
    Table[
      ParametricPlot[
        EvaluatePair[xy, s] // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Contour"]
      ]
      , {xy, {xyContour}}
    ],
    {}
    , ImageSize -> 0.45 * 0.9 ImageSizeTextWidth
  ];
  plot["contour"] = Show[
    plot["common"],
    (* Label for contour *)
    yContourLabel = Way[yMaxTraced, yMax, 0.78];
    sContourLabel = SeekRoot[
      xyContour[[2]][#] - yContourLabel &,
      {0, DomainEnd[xyContour]}
      , 10
    ];
    Graphics @ {
      Text[
        (*Italicise["T"]
          == Subscript[Italicise["h"], Style["\[NegativeVeryThinSpace]\[Bullet]", Magnification -> 1.8]]*)
        Row @ {Italicise["T"], "\[Hyphen]contour"}
          // LaTeXStyle
          // Style[#, LabelSize["Straight"]] &
        , EvaluatePair[xyContour, sContourLabel]
        , {-1.2, 0}
      ]
    },
    {}
  ];
  plot["serrated"] = Show[
    plot["common"],
    (* Traced boundaries (indentations) *)
    Table[
      ParametricPlot[
        EvaluatePair[xyTracedListUpper[[n]], s]
          // IncludeYReflection
          // Evaluate
        , {s, sStart[n], sEnd[n]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {n, Keys[intersectionCurveAss]}
    ],
    {}
  ];
  (* Export *)
  Table[
    plot[case]
      // Ex @ FString["wedge_obtuse-traced-boundaries-approximation-{case}.pdf"]
    , {case, {"contour", "serrated"}}
  ]
]


(* ::Section:: *)
(*Figure: roughness profile along contour (wedge_obtuse-roughness-profile-contour)*)


Module[
  {
    sMax,
    dummyForTrailingCommas
  },
  (* Plot range *)
  sMax = 5;
  (* Make plot *)
  Plot[
    {1, rhoIndentationsRawProfile[s]}
    , {s, -sMax, sMax}
    , AspectRatio -> 0.7
    , AxesLabel -> {Italicise["s"], "\[Rho]" // LaTeXStyle}
    , ImageSize -> 0.5 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotPoints -> 2
    , PlotRange -> {0, All}
    , PlotStyle -> {BoundaryTracingStyle["Contour"], Black}
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["wedge_obtuse-roughness-profile-contour.pdf"]


(* ::Section:: *)
(*Figure: triangular groove indentations domain and mesh (wedge_obtuse-with-grooves-*)*)


Module[
  {
    sigma,
    xMin, xMax, yMax,
    boundaryPoints, mesh,
    plot,
    dummyForTrailingCommas
  },
  (* Value of sigma *)
  sigma = 0.1;
  (* Plot range *)
  xMin = -0.65;
  xMax = 0.4;
  yMax = 1.05;
  (* Boundary points for domain plot *)
  boundaryPoints = xyIndentationList[sigma];
  (* Mesh for mesh plot *)
  mesh =
    FString["mesh/wedge_obtuse-mesh-indented-sigma-{sigma}.txt"]
      // Import // Uncompress // First;
  (* Common plot *)
  plot["common"] =
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , ImageSize -> 0.85 * 0.45 ImageSizeTextWidth
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ];
  (* Make domain plot *)
  plot["domain"] =
    Show[
      plot["common"],
      Graphics @ {
        Directive[
          FaceForm[BoundaryTracingStyle["Wall"] // Last],
          EdgeForm[GeneralStyle["DefaultThick"]]
        ],
        Polygon[boundaryPoints]
      },
      {}
    ];
  (* Make mesh plot *)
  plot["mesh"] =
    Show[
      plot["common"],
      mesh["Wireframe"],
      {}
    ];
  (* Export *)
  Table[
    plot[case] // Ex @ FString["wedge_obtuse-with-grooves-{case}.pdf"]
    , {case, {"domain", "mesh"}}
  ]
]


(* ::Section:: *)
(*Figure: triangular groove indentations height rise profiles (wedge_obtuse-height-rise-profiles-grooves)*)


Module[
  {
    yMax, tMax,
    makeData,
    tNumerical,
      fineLengthScale, sStart, sEnd,
    boundaryPoints, data,
    plotList, plot,
    dummyForTrailingCommas
  },
  (* Plot range *)
  yMax = 3;
  tMax = 0.85;
  (* Function converting <{x, y}-list, T> to <{y, T}-list> *)
  makeData[xyList_, tFun_] :=
    {
      (* y-coordinates *)
      xyList[[All, 2]],
      (* T-values *)
      tFun @@@ xyList
    } // Transpose // Select[Abs @ #[[1]] <= yMax &];
  (* Import non-indented numerical solution *)
  tNumerical["non"] =
    "solution/wedge_obtuse-solution-non_indented.txt"
      // Import // Uncompress;
  (* Non-indented boundary (i.e. T-contour) plotting points *)
  fineLengthScale = 0.01;
  sStart = DomainStart[xyContourIndentations];
  sEnd = DomainEnd[xyContourIndentations];
  boundaryPoints["non"] =
    Table[
      EvaluatePair[xyContourIndentations, s]
      , {s, UniformRange[sStart, sEnd, fineLengthScale][[;; ;; 10]]}
    ];
  (* Non-indented data to be plotted *)
  data["non"] = makeData[boundaryPoints["non"], tNumerical["non"]];
  (* Make list of plots for various sigma: *)
  plotList = Table[
    (* Import indented numerical solution *)
    tNumerical[sigma] =
      FString @ "solution/wedge_obtuse-solution-indented-sigma-{sigma}.txt"
        // Import // Uncompress;
    (* Indented data to be plotted *)
    data[sigma] = makeData[xyIndentationList[sigma], tNumerical[sigma]];
    (* Make plot *)
    plot = ListPlot[
      {data["non"], data[sigma]}
      , AspectRatio -> 0.2
      , AxesLabel -> Italicise /@ {"y", "T"}
      , AxesOrigin -> {-yMax, Automatic}
      , Epilog -> {
          (* PlotLabel gets screwed up by GraphicsColumn *)
          Text[
            LaTeXStyle["\[Sigma]"] == sigma
              // LaTeXStyle
              // Style[#, LabelSize["Label"]] &
            , {0, tMax}
            , {0, -0.5}
          ]
        }
      , Joined -> True
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , PlotRange -> {{-yMax, yMax}, {0, tMax}}
      , PlotRangeClipping -> None
      , PlotStyle -> {Gray, Directive[Black, AbsoluteThickness[1]]}
      , TicksStyle -> LabelSize["Tick"]
    ]
    , {sigma, sigmaValues // Rest // Most}
  ];
  (* Export *)
  GraphicsColumn[plotList
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.3
    , Spacings -> 0.01 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-height-rise-profiles-grooves.pdf"]


(* ::Section:: *)
(*Figure: effective contact angle profile along contour (wedge_obtuse-effective-contact-angle-profile-contour)*)


Module[
  {
    sMax,
    dummyForTrailingCommas
  },
  (* Plot range *)
  sMax = 5;
  (* Make plot *)
  Plot[
    {gammaTracingIndentations, gammaEffectiveRawProfile[s]}
    , {s, -sMax, sMax}
    , AxesLabel -> {
        Italicise["s"],
        Subscript["\[Gamma]" // LaTeXStyle, "e"]
      }
    , ImageSize -> 0.5 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotPoints -> 2
    , PlotRange -> {0, (gpdTracingIndentations + 5) Degree}
    , PlotStyle -> {BoundaryTracingStyle["Contour"], Black}
    , Ticks -> {Automatic,
        Table[
          {
            gpd * Degree,
            SeparatedRow["VeryThin"][gpd, Style["\[Degree]", Magnification -> 1.2]]
          }
          , {gpd, Range[0, gpdTracingIndentations + 5, 10]}
        ]
      }
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["wedge_obtuse-effective-contact-angle-profile-contour.pdf"]
