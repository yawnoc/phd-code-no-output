(* ::Package:: *)

(* ::Text:: *)
(*Curvilinear.wl (for Wolfram Mathematica)*)


(* ::Text:: *)
(*A Wolfram Language Package (.wl) for curvilinear coordinates.*)
(*ABSOLUTELY NO WARRANTY, i.e. "GOD SAVE YOU"*)


(* ::Text:: *)
(*See Section 3 (Page 5 of manuscripts/boundary-tracing.pdf).*)


(* ::Section:: *)
(*Package definitions*)


(* ::Subsection:: *)
(*Start of package*)


BeginPackage["Curvilinear`"];


(* ::Subsection:: *)
(*Clear existing definitions if any*)


Unprotect["Curvilinear`*"];
ClearAll["Curvilinear`*"];
ClearAll["Curvilinear`*`*"];


(* ::Subsection:: *)
(*Mention non-private symbols*)


(* ::Subsubsection:: *)
(*Bipolar coordinates*)


{
  (* Coordinate transformations *)
  XBipolar,
  YBipolar,
  XYBipolar,
  (* Inverse coordinate transformations *)
  UBipolar,
  VBipolar,
  UVBipolar,
  (* Scale factors *)
  HBipolar,
  HUBipolar,
  HVBipolar,
  (* Abbreviations *)
  CBipolar,
  SBipolar,
  (* Cartesian components of local orthonormal basis *)
  AUBipolar,
  AVBipolar
};


(* ::Subsubsection:: *)
(*Polar coordinates*)


{
  (* Coordinate transformations *)
  XPolar,
  YPolar,
  XYPolar,
  (* Inverse coordinate transformations *)
  RPolar,
  PhiPolar,
  RPhiPolar,
  (* Scale factors *)
  HRPolar,
  HPhiPolar,
  (* Cartesian components of local orthonormal basis *)
  ARPolar,
  APhiPolar
};


(* ::Subsection:: *)
(*Private scope*)


(* ::Subsubsection:: *)
(*Start of private scope*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*Bipolar coordinates (scaled)*)


(* ::Text:: *)
(*See r3: Section 3 (Page r3-5 of manuscripts/radiation-3-bipolar.pdf).*)


(* Coordinate transformations *)
XBipolar[u_, v_] := Sinh[v] / (Cosh[v] - Cos[u]);
YBipolar[u_, v_] := Sin[u] / (Cosh[v] - Cos[u]);
XYBipolar[u_, v_] := {XBipolar, YBipolar} @@ {u, v} // Through // Evaluate;


(* Inverse coordinate transformations *)
UBipolar[x_, y_] := ArcTan[x^2 + y^2 - 1, 2 y];
VBipolar[x_, y_] := ArcTanh[2 x / (x^2 + y^2 + 1)];
UVBipolar[x_, y_] := {UBipolar, VBipolar} @@ {x, y} // Through // Evaluate;


(* Scale factors (both are the same) *)
HBipolar[u_, v_] := 1 / (Cosh[v] - Cos[u]);
HUBipolar[u_, v_] := HBipolar[u, v] // Evaluate;
HVBipolar[u_, v_] := HBipolar[u, v] // Evaluate;


(* Abbreviations *)
CBipolar[u_, v_] := Cos[u] Cosh[v] - 1;
SBipolar[u_, v_] := Sin[u] Sinh[v];


(* Cartesian components of local orthonormal basis *)
AUBipolar[u_, v_] :=
  HBipolar[u, v] {-SBipolar[u, v], CBipolar[u, v]} // Evaluate;
AVBipolar[u_, v_] :=
  HBipolar[u, v] {-CBipolar[u, v], -SBipolar[u, v]} // Evaluate;


(* ::Subsubsection:: *)
(*Polar coordinates*)


(* ::Text:: *)
(*See r2: Section 1 (Page r2-1 of manuscripts/radiation-2-line.pdf).*)


(* Coordinate transformations *)
XPolar[r_, phi_] := r Cos[phi];
YPolar[r_, phi_] := r Sin[phi];
XYPolar[r_, phi_] := {XPolar, YPolar} @@ {r, phi} // Through // Evaluate;


(* Inverse coordinate transformations *)
RPolar[x_, y_] := Sqrt[x^2 + y^2];
PhiPolar[x_, y_] := ArcTan[x, y];
RPhiPolar[x_, y_] := {RPolar, PhiPolar} @@ {x, y} // Through // Evaluate;


(* Scale factors *)
HRPolar[r_, phi_] := 1;
HPhiPolar[r_, phi_] := r;


(* Cartesian components of local orthonormal basis *)
ARPolar[r_, phi_] := {Cos[phi], Sin[phi]};
APhiPolar[r_, phi_] := {-Sin[phi], Cos[phi]};


(* ::Subsubsection:: *)
(*End of private scope*)


End[];


(* ::Subsection:: *)
(*Protect definitions*)


Protect["Curvilinear`*"];


(* ::Subsection:: *)
(*End of package*)


EndPackage[];
