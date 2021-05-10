(* ::Package:: *)

(* ::Text:: *)
(*FigureStyles.wl (for Wolfram Mathematica)*)


(* ::Text:: *)
(*Semantic styling for figures.*)
(*ABSOLUTELY NO WARRANTY, i.e. "GOD SAVE YOU"*)


(* ::Section:: *)
(*Package definitions*)


(* ::Subsection:: *)
(*Start of package*)


BeginPackage["FigureStyles`", {"Conway`"}];


(* ::Subsection:: *)
(*Clear existing definitions if any*)


Unprotect["FigureStyles`*"];
ClearAll["FigureStyles`*"];
ClearAll["FigureStyles`*`*"];


(* ::Subsection:: *)
(*Mention non-private symbols*)


{
  GeneralStyle,
  BoundaryTracingStyle,
  ImageSizeTextWidth,
  ImageSizeTextWidthBeamer,
  LabelSize,
  SquigglyArrow,
  SlidesStyle,
  {}
};


(* ::Subsection:: *)
(*Private scope*)


(* ::Subsubsection:: *)
(*Start of private scope*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*GeneralStyle*)


GeneralStyle[type_String : Automatic] := Association[
  "AmbientLighting" -> {{"Ambient"}, White},
  "Dashed" -> AbsoluteDashing @ {5, 4},
  "DefaultThick" -> AbsoluteThickness[1.5],
  "Dotted" -> AbsoluteDashing @ {0.3, 4},
  "Point" -> PointSize[Large],
  "Thick" -> AbsoluteThickness[2.5],
  "Translucent" -> Opacity[0.7],
  "VeryThick" -> AbsoluteThickness[5],
  Automatic -> Automatic
][type];


(* ::Subsubsection:: *)
(*BoundaryTracingStyle*)


BoundaryTracingStyle[type : (_String | Automatic) : Automatic] := Association[
  "Background" -> Directive[
    GeneralStyle["DefaultThick"],
    GrayLevel[0.7]
  ],
  "Contour" -> Directive[
    GeneralStyle["DefaultThick"],
    GeneralStyle["Dotted"],
    Black
  ],
  "ContourImportant" -> Directive[
    GeneralStyle["Thick"],
    Gray
  ],
  "ContourPlain" -> Directive[
    GeneralStyle["DefaultThick"],
    GrayLevel[0.6]
  ],
  "Edge3D" -> Thick,
  "NonViable" -> Directive[
    GeneralStyle["Translucent"],
    GrayLevel[0.8]
  ],
  "Solution3D" -> White,
  "Terminal" -> Directive[
    GeneralStyle["DefaultThick"],
    GeneralStyle["Dashed"],
    Black
  ],
  "Traced" -> Directive[
    GeneralStyle["Thick"],
    Black
  ],
  "Unphysical" -> Black,
  "Viable" -> White,
  "Wall" -> Directive[
    GeneralStyle["DefaultThick"],
    Gray
  ],
  "Wall3D" -> GrayLevel[0.6],
  Automatic -> Black
][type];


BoundaryTracingStyle[typeSeq__String] :=
  Directive @@ BoundaryTracingStyle /@ {typeSeq};


(* ::Subsubsection:: *)
(*ImageSizeTextWidth*)


ImageSizeTextWidth =
  Module[{a4width, bindingoffset, inner, outer},
    a4width = 21 Conway`ImageSizeCentimetre;
    bindingoffset = 2 Conway`ImageSizeCentimetre;
    inner = 2 Conway`ImageSizeCentimetre;
    outer = 2.5 Conway`ImageSizeCentimetre;
    a4width - (bindingoffset + inner + outer)
  ];


(* ::Subsubsection:: *)
(*ImageSizeTextWidthBeamer*)


ImageSizeTextWidthBeamer = 10.8 Conway`ImageSizeCentimetre;


(* ::Subsubsection:: *)
(*LabelSize*)


LabelSize[type_String : Automatic] := Association[
  "Axis" -> 12,
  "Label" -> 12,
  "LabelOmega" -> If[$OperatingSystem == "Windows", 13, 15],
  "Legend" -> 10,
  "Point" -> 11,
  "PointBracket" -> 15,
  "Straight" -> 10,
  "Tick" -> 9,
  Automatic -> Automatic
][type];


(* ::Subsubsection:: *)
(*SquigglyArrow*)


SquigglyArrow[{xBase_, yBase_}, phi_: 0, size_: 1] :=
  Module[{yList, listLength, xList, pointList},
    (* List of y-coordinates (for x-coordinates 0, ..., 1) *)
    (*
      NOTE:
      - Nonzero y-coordinate in the ConstantArrays is an optical correction
        to make the straight parts of the arrow appear in line.
      - Greater number of elements in the second ConstantArray is a correction
        to make the arrow appear balanced when the arrowhead is included.
    *)
    yList =
      {
        ConstantArray[-0.1, 5],
        1, 0, -1, 0, 1, 0, -1,
        ConstantArray[+0.1, 7],
        {}
      }
        // 0.15 # &
        // Flatten;
    (* List of x-coordinates *)
    listLength = Length[yList];
    xList = Subdivide[0, 1, listLength - 1];
    (* Make B-spline arrow *)
    pointList = {xList, yList} // Transpose;
    Arrow @ BSplineCurve[
      pointList
        // ScalingTransform @ {size, size}
        // RotationTransform[phi, {0, 0}]
        // TranslationTransform @ {xBase, yBase}
    ]
  ];


(* ::Subsubsection:: *)
(*SlidesStyle*)


SlidesStyle[type_String : Automatic] := Association[
  "Boundary" -> RGBColor["#9400D3"],
  "InteriorRegion" -> RGBColor["#DEF0FF"],
  "Source" -> Red,
  "SourceRegion" -> RGBColor["#FFD9D9"],
  Automatic -> Automatic
][type];


(* ::Subsubsection:: *)
(*End of private scope*)


End[];


(* ::Subsection:: *)
(*Protect definitions*)


Protect["FigureStyles`*"];


(* ::Subsection:: *)
(*End of package*)


EndPackage[];
