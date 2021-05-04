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
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "dip_coating"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Contact angle*)


gamma = 10 Degree;


(* ::Subsection:: *)
(*Domain*)


(* ::Subsubsection:: *)
(*Dipped object (rectangular prism)*)


{prismXMax, prismYMax} = {3, 2} / 2;


(* ::Subsubsection:: *)
(*Vat (cylinder)*)


vatRMax = 5;


(* ::Subsubsection:: *)
(*Mesh for region occupied by liquid*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["dip_coating-mesh.txt",
  Module[
    {
      fineLengthScale, coarseLengthScale,
      innerBoundaryPoints, innerBoundaryNumber,
      outerBoundaryPoints, outerBoundaryNumber,
      allBoundaryPoints, allBoundaryElements,
      boundaryMesh, mesh,
      dummyForTrailingCommas
    },
    (* Length scales *)
    fineLengthScale = 0.01;
    coarseLengthScale = 0.5;
    (* Inner boundary (prism) *)
    innerBoundaryPoints =
      Join[
        Table[{prismXMax, y}, {y, UniformRange[-prismYMax, prismYMax, +fineLengthScale]}] // Rest,
        Table[{x, prismYMax}, {x, UniformRange[prismXMax, -prismXMax, -fineLengthScale]}] // Rest,
        Table[{-prismXMax, y}, {y, UniformRange[prismYMax, -prismYMax, -fineLengthScale]}] // Rest,
        Table[{x, -prismYMax}, {x, UniformRange[-prismXMax, prismXMax, +fineLengthScale]}] // Rest,
        {}
      ];
    innerBoundaryNumber = Length[innerBoundaryPoints];
    (* Outer boundary (rim of vat) *)
    outerBoundaryPoints =
      Table[
        XYPolar[vatRMax, phi]
        , {phi, UniformRange[0, 2 Pi, coarseLengthScale / vatRMax]}
      ] // Rest;
    outerBoundaryNumber = Length[outerBoundaryPoints];
    (* Full boundary *)
    allBoundaryPoints = Join[innerBoundaryPoints, outerBoundaryPoints];
    allBoundaryElements =
      LineElement /@ {
        Table[{n, n + 1}, {n, innerBoundaryNumber}]
          // Mod[#, innerBoundaryNumber, 1] &
          ,
        Table[{n, n + 1}, {n, outerBoundaryNumber}]
          // Mod[#, outerBoundaryNumber, 1] &
          // # + innerBoundaryNumber &
          ,
        Nothing
      };
    (* Build mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> allBoundaryPoints,
        "BoundaryElements" -> allBoundaryElements,
        {}
      ];
    mesh =
      ToElementMesh[
        boundaryMesh
        , "ImproveBoundaryPosition" -> True
        , MeshRefinementFunction -> MeshRefinementUniform[coarseLengthScale]
        , "RegionHoles" -> {0, 0}
      ];
    mesh // Compress
  ]
]


(* ::Subsubsection:: *)
(*Predicate for wetting boundary*)


predicateWet =
  Function[{x, y},
    RPolar[x, y] < Way[RPolar[prismXMax, prismYMax], vatRMax] // Evaluate
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Solve PDE*)


(* (This is slow (~5 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["dip_coating-solution.txt",
  Module[{mesh, tSol, x, y},
    mesh = Import["dip_coating-mesh.txt"] // Uncompress;
    tSol = SolveLaplaceYoung[10 Degree, mesh, predicateWet];
    tSol // Compress
  ]
]


(* ::Section:: *)
(*Mesh*)


Module[{mesh},
  mesh = Import["dip_coating-mesh.txt"] // Uncompress;
  Show[
    mesh["Wireframe"]
    , PlotLabel -> StringTemplate["`` elements"] @ Length @ mesh[[2, 1, 1]]
  ]
] // Ex["dip_coating-mesh.pdf"]


(* ::Section:: *)
(*Solution*)


Module[{tSol, mesh, x, y},
  tSol = Import["dip_coating-solution.txt"] // Uncompress;
  mesh = tSol["ElementMesh"];
  Plot3D[tSol[x, y], Element[{x, y}, mesh]
    , PlotRange -> Full
  ]
] // Ex["dip_coating-solution.png"]


(* ::Section:: *)
(*Figure: coating profile ideal vs actual (dip_coating-*)*)


Module[
  {
    tSol, mesh,
    tSolMin, tSolMax,
    zMin, zMax,
    tIdeal, tActual,
    globalLighting,
    coatingOffset, coatingOptions,
    tCoating,
    dummyForTrailingCommas
  },
  (* Import numerical solution *)
  tSol = Import["dip_coating-solution.txt"] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Extremal values of solution *)
  tSolMin = tSol[prismXMax, prismYMax];
  tSolMax = Max[tSol[prismXMax, 0], tSol[0, prismYMax]];
  (* Vertical plot range *)
  zMin = Way[tSolMin, tSolMax, -0.1];
  zMax = Way[tSolMin, tSolMax, +1.3];
  (* Options *)
  globalLighting = {{"Ambient"}, White};
  coatingOffset = 0.001;
  coatingOptions = {
    BoundaryStyle -> BoundaryTracingStyle["Edge3D"],
    Lighting -> globalLighting,
    Mesh -> None,
    PlotPoints -> 50,
    PlotStyle -> BoundaryTracingStyle["Solution3D"],
    Nothing
  };
  (* Ideal vs actual *)
  tCoating["ideal"] = Evaluate @ Way[tSolMin, tSolMax, 0.6] &;
  tCoating["actual"] = tSol;
  (* Plot and export *)
  Table[
    Show[
      (* Prism *)
      Graphics3D @ {
        BoundaryTracingStyle["Wall3D"],
        EdgeForm @ BoundaryTracingStyle["Edge3D"],
        Cuboid[
          {-prismXMax, -prismYMax, zMin},
          {+prismXMax, +prismYMax, zMax}
        ]
      },
      (* Coating along x == x_max *)
      ParametricPlot3D[
        {
          prismXMax + coatingOffset, y,
          Way[zMin, tCoating[type][prismXMax, y], p]
        }
        , {y, -prismYMax, prismYMax}
        , {p, 0, 1}
        , coatingOptions // Evaluate
      ],
      (* Coating along y == -y_max *)
      ParametricPlot3D[
        {
          x, -prismYMax - coatingOffset,
          Way[zMin, tCoating[type][x, -prismYMax], p]
        }
        , {x, -prismXMax, prismXMax}
        , {p, 0, 1}
        , coatingOptions // Evaluate
      ],
      {}
      , Lighting -> globalLighting
      , Boxed -> False
      , BoxRatios -> {Automatic, Automatic, 1.6}
      , ImageSize -> 0.35 ImageSizeTextWidth
      , ViewPoint -> {2.4, -3.5, 1.6}
    ]
      //
        Ex[
          StringTemplate["dip_coating-``.png"][type]
          , Background -> None
          , ImageResolution -> 4 BasicImageResolution
        ]
    , {type, {"ideal", "actual"}}
  ]
]
