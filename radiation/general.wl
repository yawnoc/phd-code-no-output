(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "general"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Section:: *)
(*Figure: radiation--conduction BVP (radiation-conduction-bvp)*)


(* RUN ON WINDOWS FOR ITALICISED OMEGA *)
Module[
  {
    conductionA, conductionB, conductionPhi,
    conductionXY, conductionNormalPhi,
    bathX, bathY, bathA, bathB,
    textStyle, textStyleGreek,
    radiationArrowClearanceFactor,
    dummyForTrailingCommas
  },
  (* Conduction ellipse *)
  {conductionA, conductionB} = {8, 5};
  conductionPhi = 20 Degree;
  (* Geometry *)
  (* (see <https://math.stackexchange.com/a/990013> for angle of normal) *)
  conductionXY[ang_] := {conductionA Cos[ang], conductionB Sin[ang]};
  conductionNormalPhi[ang_] := ArcTan[conductionB Cos[ang], conductionA Sin[ang]];
  (* Heat bath ellipse *)
  {bathX, bathY} = {-3, -2};
  {bathA, bathB} = {2.5, 1.5};
  (* Diagram *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleGreek = Style[#, LabelSize["LabelOmega"]] & @* LaTeXStyle;
  Show[
    (* Conduction ellipse *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Circle[{0, 0}, {conductionA, conductionB}]
        // Rotate[#, conductionPhi] &
    },
    Graphics @ {
      Text[
        SeparatedRow["\[VeryThinSpace]" // textStyle] @@ {
          "conduction" // textStyle,
          Italicise["\[CapitalOmega]"] // textStyleGreek
        }
        , -{0.3 bathX, 0.8bathY}
      ]
    },
    (* Radiation arrows *)
    radiationArrowClearanceFactor = 1.1;
    Graphics @ {Arrowheads[0.02],
      Table[
        SquigglyArrow[
          radiationArrowClearanceFactor * conductionXY[ang]
          , conductionNormalPhi[ang]
          , 2
        ]
        , {ang, Subdivide[0, 2 Pi, 8]}
      ]
        // Rotate[#, conductionPhi] &
    },
    Graphics @ {
      Text[
        Column[
          {
            "radiation" // textStyle,
            SeparatedRow[
              If[$OperatingSystem == "Windows",
                StringRepeat["\[NegativeThickSpace]", 9] <> "\[NegativeVeryThinSpace]" // textStyle,
                ""
              ]
            ] @@ {
              "\[PartialD]" // textStyleGreek,
              Italicise["\[CapitalOmega]"] // textStyleGreek
            }
          }
          , Alignment -> Right
          , Spacings -> If[$OperatingSystem == "Windows", 0, 0.3]
        ]
        , conductionXY[5/8 Pi] // RotationTransform[conductionPhi]
        , {1.1, If[$OperatingSystem == "Windows", -1, -1.15]}
      ]
    },
    (* Heat bath ellipse *)
    Graphics @ {Directive[FaceForm @ GrayLevel[0.9], EdgeForm[Black]],
      Disk[{bathX, bathY}, {bathA, bathB}]
    },
    Graphics @ {
      Text[
        "heat" // textStyle
        , {bathX, bathY}
        , {0, If[$OperatingSystem == "Windows", -0.15, -0.1]}
      ]
    },
    {}
    , ImageSize -> 0.45 ImageSizeTextWidth
  ]
] // Ex["radiation-conduction-bvp.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    conductionA, conductionB, conductionPhi,
    conductionXY, conductionNormalPhi,
    bathX, bathY, bathA, bathB,
    textStyle, textStyleGreek,
    radiationArrowClearanceFactor,
    dummyForTrailingCommas
  },
  (* Conduction ellipse *)
  {conductionA, conductionB} = {8, 5};
  conductionPhi = 20 Degree;
  (* Geometry *)
  (* (see <https://math.stackexchange.com/a/990013> for angle of normal) *)
  conductionXY[ang_] := {conductionA Cos[ang], conductionB Sin[ang]};
  conductionNormalPhi[ang_] := ArcTan[conductionB Cos[ang], conductionA Sin[ang]];
  (* Heat bath ellipse *)
  {bathX, bathY} = {-3, -2};
  {bathA, bathB} = {2.5, 1.5};
  (* Diagram *)
  textStyle = Style[#, 10] &;
  textStyleGreek = Style[#, 10] &;
  Show[
    (* Conduction ellipse *)
    Graphics @ {
      EdgeForm @ Directive[BoundaryTracingStyle["Traced"], SlidesStyle["Boundary"]],
      FaceForm @ SlidesStyle["InteriorRegion"],
      Disk[{0, 0}, {conductionA, conductionB}]
        // Rotate[#, conductionPhi] &
    },
    Graphics @ {
      Text[
        SeparatedRow["\[VeryThinSpace]" // textStyle] @@ {
          "conduction" // textStyle,
          Italicise["\[CapitalOmega]"] // textStyleGreek
        }
        , -{0.3 bathX, 0.8bathY}
      ]
    },
    (* Radiation arrows *)
    radiationArrowClearanceFactor = 1.1;
    Graphics @ {Arrowheads[0.03],
      Table[
        SquigglyArrow[
          radiationArrowClearanceFactor * conductionXY[ang]
          , conductionNormalPhi[ang]
          , 2.5
        ]
        , {ang, Subdivide[0, 2 Pi, 8]}
      ]
        // Rotate[#, conductionPhi] &
    },
    Graphics @ {
      Text[
        Column[
          {
            "radiation" // textStyle,
            SeparatedRow[
              If[$OperatingSystem == "Windows",
                StringRepeat["\[NegativeThickSpace]", 9] <> "\[NegativeVeryThinSpace]" // textStyle,
                ""
              ]
            ] @@ {
              "\[PartialD]" // textStyleGreek,
              Italicise["\[CapitalOmega]"] // textStyleGreek
            }
          }
          , Alignment -> Right
          , Spacings -> 0
        ]
        , conductionXY[4.8/8 Pi] // RotationTransform[conductionPhi]
        , {1.1, -1}
      ]
    },
    (* Heat bath ellipse *)
    Graphics @ {
      FaceForm @ SlidesStyle["SourceRegion"],
      EdgeForm @ SlidesStyle["Source"],
      Disk[{bathX, bathY}, {bathA, bathB}]
    },
    Graphics @ {
      Text[
        "heat" // textStyle
        , {bathX, bathY}
        , {0, 0}
      ]
    },
    {}
    , ImageSize -> 0.5 ImageSizeTextWidthBeamer
  ]
] // Ex["radiation-conduction-bvp-slides.pdf"];


(* ::Section:: *)
(*Figure: self viewing radiation elements (self-viewing-radiation-elements)*)


Module[
  {
    r, n, rStar, nStar, dStar, theta, thetaStar,
    normalVectorLength, elementWidth,
    angleMarkerRadius, angleMarkerRadiusStar,
    nTip, nStarTip,
    zVector, elementGeneric, element, elementStar,
    sweep,
    centre, eye, canvasNormal, canvasPoint,
    project,
    elementStyle, positionStyle,
    normalVectorStyle, displacementVectorStyle, angleMarkerStyle,
    starred, textStyle,
    dummyForTrailingCommas
  },
  (* Coordinates *)
  r = {-1, 0, 0};
  n = {0, 0, 1} // Normalize;
  rStar = {0.7, 0, 0.3};
  nStar = {-0.4, -0.2, 0.7} // Normalize;
  dStar = r - rStar;
  theta = ArcCos[n . (-dStar) / Norm[dStar]];
  thetaStar = ArcCos[nStar . dStar / Norm[dStar]];
  (* Plotting constants *)
  normalVectorLength = 0.8;
  elementWidth = 0.35;
  angleMarkerRadius = 0.2 normalVectorLength;
  angleMarkerRadiusStar = 0.25 normalVectorLength;
  nTip = r + normalVectorLength * n;
  nStarTip = rStar + normalVectorLength * nStar;
  (* Element coordinates *)
  zVector = {0, 0, 1};
  elementGeneric = elementWidth/2 {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
  element = (elementGeneric
    // RotationTransform[{zVector, n}]
    // TranslationTransform[r]
  );
  elementStar = (elementGeneric
    // RotationTransform[{zVector, nStar}]
    // TranslationTransform[rStar]
  );
  (* Sweeping transformation (for angle markers) *)
  sweep[origin_, angle_][startVector_, endVector_][vector_] := (
    vector
      // RotationTransform[angle, {startVector, endVector}, {0, 0, 0}]
      // TranslationTransform[origin]
  );
  (* Perspective geometry *)
  centre = Way[r, rStar];
  canvasNormal = {0.15, -1, 0.4};
  canvasPoint = centre + canvasNormal;
  eye = centre + 10 canvasNormal;
  project[point_] :=
    OnePointPerspective[2][eye, canvasNormal, canvasPoint, zVector][point];
  (* Plot styles *)
  elementStyle = Directive[FaceForm[None], EdgeForm[Black]];
  positionStyle = PointSize[Medium];
  normalVectorStyle = Directive[GeneralStyle["Thick"], Arrowheads[0.06]];
  displacementVectorStyle = Arrowheads @ {{0.05, 0.55}};
  angleMarkerStyle = Directive[Black, Thickness[Small]];
  (* Text functions *)
  starred[expr_] := Superscript[expr, "\[NegativeVeryThinSpace]\[FivePointedStar]"];
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Make scene *)
  Show[
    (* Element dA *)
    Graphics @ {
      positionStyle, Point[project /@ {r}],
        Text[Embolden["r"] // textStyle
          , project[r]
          , {2.35, -0.27}
        ],
      elementStyle, Polygon[project /@ element],
        Text[SeparatedRow[""] @ {"d", Italicise["A"]} // textStyle
          , project[r]
          , {0.5, 1.6}
        ],
      normalVectorStyle, Arrow[project /@ {r, nTip}],
        Text[Embolden["n"] // textStyle
          , project[nTip]
          , {-2.5, 0}
        ],
      {}
    },
    (* Angle theta *)
    ParametricPlot[
      sweep[r, ang][n, -dStar][angleMarkerRadius * n]
        // project
        // Evaluate
      , {ang, 0, theta}
      , PlotPoints -> 2
      , PlotStyle -> angleMarkerStyle
    ],
    Graphics @ {
      Text["\[Theta]" // textStyle
        , sweep[r, theta/2][n, -dStar][angleMarkerRadius * n] // project
        , {-1.5, -0.6}
      ]
    },
    (* Element dA-star *)
    Graphics @ {
      positionStyle, Point[project /@ {rStar}],
        Text[
          Embolden["r"] // starred // textStyle
          , project[rStar]
          , {-1.5, -0.65}
        ],
      elementStyle, Polygon[project /@ elementStar],
        Text[SeparatedRow[""] @ {"d", Italicise["A"] // starred} // textStyle
          , project[rStar]
          , {-0.5, 1.8}
        ],
      normalVectorStyle, Arrow[project /@ {rStar, nStarTip}],
        Text[Embolden["n"] // starred // textStyle
          , project[nStarTip]
          , {-1, -0.9}
        ],
      {}
    },
    (* Angle theta-star *)
    ParametricPlot[
      sweep[rStar, ang][nStar, dStar][angleMarkerRadiusStar * nStar]
        // project
        // Evaluate
      , {ang, 0, thetaStar}
      , PlotPoints -> 2
      , PlotStyle -> angleMarkerStyle
    ],
    Graphics @ {
      Text["\[Theta]" // starred // textStyle
        , sweep[rStar, theta/2][nStar, dStar][angleMarkerRadiusStar * nStar] // project
        , {0.9, -0.3}
      ]
    },
    (* Displacement d-star *)
    Graphics @ {
      displacementVectorStyle, Arrow[project /@ {rStar, r}],
        Text[
          Embolden["d"] // starred // textStyle
          , project @ Way[rStar, r]
          , {-0.3, 1.3}
        ],
      {}
    },
    {}
    , ImageSize -> 0.5 ImageSizeTextWidth
  ]
] // Ex["self-viewing-radiation-elements.pdf"]


(* ::Section:: *)
(*Figure: self viewing radiation elements for a fin (self-viewing-radiation-elements-fin)*)


Module[
  {
    xMaxFin, yMaxFin, zMinFin,
    yFinScaled, yFin,
    xMinAxis,
    xMaxAxis, yMaxAxis, zMinAxis, zMaxAxis,
    xElement, yElement, zElement,
      dx, dy, dz,
    xStar, yStar, zStar,
      dxStar, dyStar, dzStar,
    normalVectorLength, normalVector,
      normalBase, normalTip,
      normalStarBase, normalStarTip,
    centre, canvasNormal, canvasPoint, eye, verticalVector, project,
    axisStyle, finStyle,
    elementStyle, positionStyle, normalVectorStyle,
    differentialStyle, guideStyle,
    starred, textStyle,
    dummyForTrailingCommas
  },
  (* Fin *)
  {xMaxFin, yMaxFin, zMinFin} = {1.6, 0.48, -2.3};
  yFinScaled[m_][x_] := m x + (1-m) x^2;
  yFin[x_] := yMaxFin * yFinScaled[0.13][x / xMaxFin];
  (* Axes *)
  xMinAxis = -0.12 xMaxFin;
  xMaxAxis = 1.25 xMaxFin;
  yMaxAxis = 1.3 yMaxFin;
  zMinAxis = Way[zMinFin, 0, 0];
  zMaxAxis = Way[zMinFin, 0, +1.15];
  (* Non-starred element *)
  {xElement, zElement} = {0.13 xMaxFin, 0.13 zMinFin};
  yElement = yFin[xElement];
  dx = 0.16 xMaxFin;
  dz = -1.3 dx;
  dy = yFin[xElement + dx] - yFin[xElement];
  (* Starred element *)
  {xStar, zStar} = {0.64 xMaxFin, 0.45 zMinFin};
  yStar = yFin[xStar];
  dxStar = 1.15 dx;
  dzStar = dz;
  dyStar = yFin[xStar + dxStar] - yFin[xStar];
  (* Normal vectors *)
  normalVectorLength = 1.3 yMaxFin;
  normalVector[{x_, y_, z_}, dx_] :=
    normalVectorLength * Normalize @ {-yFin'[x + dx/2], 1, 0};
  normalBase = {xElement, yElement, zElement} + 1/2 {dx, dy, dz};
  normalTip = normalBase + normalVector[normalBase, dx];
  normalStarBase = {xStar, yStar, zStar} + 1/2 {dxStar, dyStar, dzStar};
  normalStarTip = normalStarBase + normalVector[normalStarBase, dxStar];
  (* Perspective geometry *)
  centre = 1/2 {xMaxFin, yMaxFin, zMinFin};
  canvasNormal = {0.36, 0.5, 1};
  canvasPoint = centre + canvasNormal;
  eye = centre + 20 canvasNormal;
  verticalVector = {0, 1, 0};
  project[point_] :=
    OnePointPerspective[2][eye, canvasNormal, canvasPoint, verticalVector][point];
  (* Plotting styles *)
  axisStyle = Darker[Gray];
  finStyle = Directive[Black, Thickness[Medium]];
  elementStyle = Directive[FaceForm[None], EdgeForm[Black]];
  positionStyle = PointSize[Medium];
  normalVectorStyle = Directive[GeneralStyle["Thick"], Arrowheads[0.06]];
  differentialStyle = Directive[CapForm["Round"], Thickness[0.015]];
  guideStyle = Directive[Black, Thickness[Medium], Dotted];
  (* Text functions *)
  starred[expr_] := Superscript[expr, "\[NegativeVeryThinSpace]\[FivePointedStar]"];
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Make scene *)
  Show[
    (* Axes *)
    Graphics @ {axisStyle,
      Line[project /@ {{xMinAxis, 0, 0}, {xMaxAxis, 0, 0}}],
        Text[
          Italicise["x"] // textStyle
          , project @ {xMaxAxis, 0, 0}
          , {-2, -0.1}
        ],
      Line[project /@ {{xMinAxis, 0, 0}, {xMinAxis, yMaxAxis, 0}}],
        Text[
          Italicise["y"] // textStyle
          , project @ {xMinAxis, yMaxAxis, 0}
          , {0, -1.25}
        ],
      Line[project /@ {{xMinAxis, 0, zMinAxis}, {xMinAxis, 0, zMaxAxis}}],
        Text[
          Italicise["z"] // textStyle
          , project @ {xMinAxis, 0, zMaxAxis}
          , {2, 0.3}
        ],
      {}
    },
    (* x-axis marks *)
    Graphics @ {axisStyle,
      Text[
        Subscript[Italicise["x"], 1] // textStyle
        , project @ {0, 0, 0}
        , {0.3, 0.65}
      ],
      Text[
        Subscript[Italicise["x"], 2] // textStyle
        , project @ {xMaxFin, 0, 0}
        , {1.4, 0.4}
      ],
      {}
    },
    (* Fin *)
    {
      (* Curved edges *)
      ParametricPlot[
        project /@ {
          {x, yFin[x], 0}, (* near upper edge *)
          {x, -yFin[x], 0}, (* near lower edge *)
          {x, yFin[x], zMinFin}, (* near far edge *)
          Nothing
        }
        , {x, 0, xMaxFin}
        , PlotPoints -> 2
        , PlotStyle -> finStyle
      ],
      (* Straight edges *)
      Graphics @ {finStyle,
        Line[project /@ {
          {xMaxFin, yMaxFin, 0},
          {xMaxFin, -yMaxFin, 0},
          {xMaxFin, -yMaxFin, zMinFin},
          {xMaxFin, yMaxFin, zMinFin},
          {xMaxFin, yMaxFin, 0}
        }],
        Line[project /@ {{0, 0, 0}, {0, 0, zMinFin}}]
      },
      {}
    },
    (* Non-starred element *)
    Graphics @ {elementStyle,
      Polygon[project /@ {
        {xElement, yElement, zElement},
        {xElement + dx, yElement + dy, zElement},
        {xElement + dx, yElement + dy, zElement + dz},
        {xElement, yElement, zElement + dz},
        Nothing
      }],
      Text[
        SeparatedRow[""] @ {"d", Italicise["A"]} // textStyle
        , project[normalBase]
        , {-3.25, -0.5}
      ],
      {}
    },
    (* Starred element *)
    Graphics @ {elementStyle,
      Polygon[project /@ {
        {xStar, yStar, zStar},
        {xStar + dxStar, yStar + dyStar, zStar},
        {xStar + dxStar, yStar + dyStar, zStar + dzStar},
        {xStar, yStar, zStar + dzStar},
        Nothing
      }],
      Text[
        SeparatedRow[""] @ {"d", Italicise["A"] // starred} // textStyle
        , project[normalStarBase]
        , {1.8, -0.4}
      ],
      {}
    },
    (* Normal vectors *)
    Graphics @ {
      positionStyle, Point[project /@ {normalBase}],
      normalVectorStyle, Arrow[project /@ {normalBase, normalTip}],
      Text[
        Embolden["n"] // textStyle
        , project[normalTip]
        , {1.8, -0.1}
      ],
      {}
    },
    Graphics @ {
      positionStyle, Point[project /@ {normalStarBase}],
      normalVectorStyle, Arrow[project /@ {normalStarBase, normalStarTip}],
      Text[
        Embolden["n"] // starred // textStyle
        , project[normalStarTip]
        , {-2.1, 0}
      ],
      {}
    },
    (* Starred differentials (not z) *)
    Graphics @ {
      differentialStyle,
      Line[project /@ {
        {xStar, yStar, 0},
        {xStar + dxStar, yStar, 0},
        {xStar + dxStar, yStar + dyStar, 0},
        {xStar, yStar, 0}
      }],
      guideStyle,
        Line[project /@ {
          {xStar, yStar, zStar},
          {xStar, yStar, 0}
        }],
        Line[project /@ {
          {xStar + dxStar, yStar + dyStar, zStar},
          {xStar + dxStar, yStar + dyStar, 0}
        }],
      Text[
        SeparatedRow[""] @ {"d", Italicise["x"] // starred} // textStyle
        , project[{xStar + dxStar/2, yStar, 0}]
        , {0.04, 0.85}
      ],
      Text[
        SeparatedRow[""] @ {"d", Italicise["y"] // starred} // textStyle
        , project[{xStar + dxStar, yStar + dyStar/2, 0}]
        , {-1.35, -0.08}
      ],
      Text[
        SeparatedRow[""] @ {"d", Italicise["s"] // starred} // textStyle
        , project @ {xStar + dxStar/2, yStar + dyStar/2, 0}
        , {-1.05, -1.15}
      ],
      {}
    },
    (* Starred differentials (z) *)
    Graphics @ {
      differentialStyle,
      Line[project /@ {
        {xMaxFin, yMaxFin, zStar},
        {xMaxFin, yMaxFin, zStar + dzStar}
      }],
      Text[
        SeparatedRow[""] @ {"d", Italicise["z"] // starred} // textStyle
        , project[{xMaxFin, yMaxFin, zStar + dzStar/2}]
        , {-1.6, 0}
      ],
      {}
    },
    {
      ParametricPlot[
        project /@ {
          {x, yFin[x], zStar},
          {x, yFin[x], zStar + dzStar},
          Nothing
        }
        , {x, xStar + dxStar, xMaxFin}
        , PlotStyle -> guideStyle
      ]
    },
    {}
    , ImageSize -> 0.55 ImageSizeTextWidth
  ]
] // Ex["self-viewing-radiation-elements-fin.pdf"]


(* ::Subsection:: *)
(*Version for slides*)


Module[
  {
    xMaxFin, yMaxFin, zMinFin,
    yFinScaled, yFin,
    xMinAxis,
    xMaxAxis, yMaxAxis, zMinAxis, zMaxAxis,
    xElement, yElement, zElement,
      dx, dy, dz,
    xStar, yStar, zStar,
      dxStar, dyStar, dzStar,
    normalVectorLength, normalVector,
      normalBase, normalTip,
      normalStarBase, normalStarTip,
    centre, canvasNormal, canvasPoint, eye, verticalVector, project,
    axisStyle, radiationStyle, sourceStyle,
    elementStyle, positionStyle, normalVectorStyle,
    differentialStyle, guideStyle,
    starred, textStyle,
    dummyForTrailingCommas
  },
  (* Fin *)
  {xMaxFin, yMaxFin, zMinFin} = {1.6, 0.48, -2.3};
  yFinScaled[m_][x_] := m x + (1-m) x^2;
  yFin[x_] := yMaxFin * yFinScaled[0.13][x / xMaxFin];
  (* Axes *)
  xMinAxis = -0.12 xMaxFin;
  xMaxAxis = 1.25 xMaxFin;
  yMaxAxis = 1.3 yMaxFin;
  zMinAxis = Way[zMinFin, 0, 0];
  zMaxAxis = Way[zMinFin, 0, +1.15];
  (* Non-starred element *)
  {xElement, zElement} = {0.13 xMaxFin, 0.13 zMinFin};
  yElement = yFin[xElement];
  dx = 0.16 xMaxFin;
  dz = -1.3 dx;
  dy = yFin[xElement + dx] - yFin[xElement];
  (* Starred element *)
  {xStar, zStar} = {0.64 xMaxFin, 0.45 zMinFin};
  yStar = yFin[xStar];
  dxStar = 1.15 dx;
  dzStar = dz;
  dyStar = yFin[xStar + dxStar] - yFin[xStar];
  (* Normal vectors *)
  normalVectorLength = 1.5 yMaxFin;
  normalVector[{x_, y_, z_}, dx_] :=
    normalVectorLength * Normalize @ {-yFin'[x + dx/2], 1, 0};
  normalBase = {xElement, yElement, zElement} + 1/2 {dx, dy, dz};
  normalTip = normalBase + normalVector[normalBase, dx];
  normalStarBase = {xStar, yStar, zStar} + 1/2 {dxStar, dyStar, dzStar};
  normalStarTip = normalStarBase + normalVector[normalStarBase, dxStar];
  (* Perspective geometry *)
  centre = 1/2 {xMaxFin, yMaxFin, zMinFin};
  canvasNormal = {0.36, 0.5, 1};
  canvasPoint = centre + canvasNormal;
  eye = centre + 20 canvasNormal;
  verticalVector = {0, 1, 0};
  project[point_] :=
    OnePointPerspective[2][eye, canvasNormal, canvasPoint, verticalVector][point];
  (* Plotting styles *)
  axisStyle = Darker[Gray];
  radiationStyle = Directive[SlidesStyle["Boundary"], Thickness[Medium]];
  sourceStyle = Directive[SlidesStyle["Source"], Thickness[Medium]];
  elementStyle = Directive[
    FaceForm[None],
    EdgeForm @ SlidesStyle["Boundary"],
    SlidesStyle["Boundary"]
  ];
  normalVectorStyle = Directive[Thick, Arrowheads[0.08]];
  differentialStyle = Directive[CapForm["Round"], Thickness[0.015]];
  guideStyle = Directive[Black, Thickness[Medium], Dotted];
  (* Text functions *)
  starred[expr_] := Superscript[expr, "\[NegativeVeryThinSpace]\[FivePointedStar]"];
  textStyle = Style[#, 8] &;
  (* Make scene *)
  Show[
    (* Axes *)
    Graphics @ {axisStyle,
      Line[project /@ {{xMinAxis, 0, 0}, {xMaxAxis, 0, 0}}],
        Text[
          Italicise["x"] // textStyle
          , project @ {xMaxAxis, 0, 0}
          , {-2, -0.1}
        ],
      Line[project /@ {{xMinAxis, 0, 0}, {xMinAxis, yMaxAxis, 0}}],
        Text[
          Italicise["y"] // textStyle
          , project @ {xMinAxis, yMaxAxis, 0}
          , {0, -1.25}
        ],
      Line[project /@ {{xMinAxis, 0, zMinAxis}, {xMinAxis, 0, zMaxAxis}}],
        Text[
          Italicise["z"] // textStyle
          , project @ {xMinAxis, 0, zMaxAxis}
          , {2, 0.3}
        ],
      {}
    },
    (* x-axis marks *)
    Graphics @ {axisStyle,
      Text[
        Subscript[Italicise["x"], 1] // textStyle
        , project @ {0, 0, 0}
        , {0.3, 0.65}
      ],
      Text[
        Subscript[Italicise["x"], 2] // textStyle
        , project @ {xMaxFin, 0, 0}
        , {1.4, 0.4}
      ],
      {}
    },
    (* Fin *)
    {
      (* Radiation edges *)
      ParametricPlot[
        project /@ {
          {x, yFin[x], 0}, (* near upper edge *)
          {x, -yFin[x], 0}, (* near lower edge *)
          {x, yFin[x], zMinFin}, (* near far edge *)
          Nothing
        }
        , {x, 0, xMaxFin}
        , PlotPoints -> 2
        , PlotStyle -> radiationStyle
      ],
      Graphics @ {radiationStyle,
        Line[project /@ {{0, 0, 0}, {0, 0, zMinFin}}]
      },
      (* Source edges *)
      Graphics @ {sourceStyle,
        Line[project /@ {
          {xMaxFin, yMaxFin, 0},
          {xMaxFin, -yMaxFin, 0},
          {xMaxFin, -yMaxFin, zMinFin},
          {xMaxFin, yMaxFin, zMinFin},
          {xMaxFin, yMaxFin, 0}
        }],
      },
      {}
    },
    (* Non-starred element *)
    Graphics @ {elementStyle,
      Polygon[project /@ {
        {xElement, yElement, zElement},
        {xElement + dx, yElement + dy, zElement},
        {xElement + dx, yElement + dy, zElement + dz},
        {xElement, yElement, zElement + dz},
        Nothing
      }],
      Text[
        Row @ {"d", Italicise["A"]} // textStyle
        , project[normalBase]
        , {-3.25, -0.5}
      ],
      {}
    },
    (* Starred element *)
    Graphics @ {elementStyle,
      Polygon[project /@ {
        {xStar, yStar, zStar},
        {xStar + dxStar, yStar + dyStar, zStar},
        {xStar + dxStar, yStar + dyStar, zStar + dzStar},
        {xStar, yStar, zStar + dzStar},
        Nothing
      }],
      Text[
        Row @ {"d", Italicise["A"] // starred} // textStyle
        , project[normalStarBase]
        , {1.8, -0.4}
      ],
      {}
    },
    (* Normal vectors *)
    Graphics @ {
      normalVectorStyle, Arrow[project /@ {normalBase, normalTip}],
      Text[
        Embolden["n"] // textStyle
        , project[normalTip]
        , {1.8, -0.1}
      ],
      {}
    },
    Graphics @ {
      normalVectorStyle, Arrow[project /@ {normalStarBase, normalStarTip}],
      Text[
        Embolden["n"] // starred // textStyle
        , project[normalStarTip]
        , {-2.1, 0}
      ],
      {}
    },
    {}
    , ImageSize -> 0.5 ImageSizeTextWidthBeamer
  ]
] // Ex["self-viewing-radiation-elements-fin-slides.pdf"];
