(* ::Package:: *)

(* ::Section:: *)
(*List of installed fonts*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/18945*)
(*Output needs to be tweaked on Mac.*)


ReplaceAll[
  FE`Evaluate @ FEPrivate`GetPopupList["MenuListFonts"],
  (string_ -> family_) :> Style[string, FontFamily -> family]
]


(* ::Section:: *)
(*LaTeX fonts*)


(* ::Text:: *)
(*See https://tex.stackexchange.com/q/524996*)
(*and http://mirrors.ctan.org/macros/latex/contrib/unicode-math/unimath-symbols.pdf*)
(*On my current machine (Debian GNU/Linux 9.12 (stretch)) are installed:*)
(*  "Latin Modern Math" which is to be used with special unicode symbols*)
(*  "Latin Modern Roman" which seems to handle everything except lowercase italic Greek*)


(* ::Text:: *)
(*$ fc-match "Latin Modern Math" file*)
(*:file=/usr/share/texmf/fonts/opentype/public/lm-math/latinmodern-math.otf*)
(*This is identical to latinmodern-math-1959/otf/latinmodern-math.otf from*)
(*<http://www.gust.org.pl/projects/e-foundry/lm-math/download>.*)


(* ::Text:: *)
(*$ fc-match "Latin Modern Roman" file*)
(*:file=/usr/share/texmf/fonts/opentype/public/lm/lmroman10-regular.otf*)
(*This is identical to lmroman10-regular.otf from*)
(*<http://www.gust.org.pl/projects/e-foundry/latin-modern/download>.*)
(*On the Windows desktop machine supplied by my university,*)
(*this is only recognised under the name "LM Roman 10".*)


(* ::Text:: *)
(*Note that "Latin Modern Math" can't handle math italic "h" (see below).*)
(*The best approach appears to be using "Latin Modern Roman" by default,*)
(*and only using "Latin Modern Math" for Greek.*)


(* Fonts *)
lmRoman10[expr_] := Style[expr, FontFamily -> "LM Roman 10"];
latinModernMath[expr_] := Style[expr, FontFamily -> "Latin Modern Math"];
SetAttributes[{lmRoman10, latinModernMath}, Listable];


(* Formatting *)
formatItalic[expr_] := Style[expr, Italic];
formatBold[expr_] := Style[expr, Bold];
SetAttributes[{formatItalic, formatBold}, Listable];


(* Table for comparing the two approaches *)
comparisonTableForm[romanString_, mathString_] :=
  TableForm[
    {
      lmRoman10[romanString],
      latinModernMath[mathString]
    },
    TableHeadings -> {{"LM Roman 10", "Latin Modern Math"}, None},
    TableSpacing -> {2, 1}
  ] // Style[#, FontSize -> 24] &;


(* ::Subsection:: *)
(*Digits*)


comparisonTableForm[
  CharacterRange["0", "9"],
  CharacterRange["0", "9"]
]


(* ::Subsection:: *)
(*Latin*)


(* ::Subsubsection:: *)
(*Uppercase*)


comparisonTableForm[
  CharacterRange["A", "Z"],
  CharacterRange["A", "Z"]
]


(* ::Subsubsection:: *)
(*Lowercase*)


comparisonTableForm[
  CharacterRange["a", "z"],
  CharacterRange["a", "z"]
]


(* ::Subsection:: *)
(*Latin, italic*)


(* ::Subsubsection:: *)
(*Uppercase*)


comparisonTableForm[
  CharacterRange["A", "Z"] // formatItalic,
  CharacterRange[16^^1D434, 16^^1D44D]
]


(* ::Subsubsection:: *)
(*Lowercase*)


(*
  https://tex.stackexchange.com/questions/524996/#comment1328261_525087
    For the missing "h", blame the Unicode Consortium;
    they already defined the math italic "h" as the Planck constant U+210E
    and didn't fill the hole in the Mathematical Italic block,
    which is a silly decision.
    \[Dash] egreg
 *)
comparisonTableForm[
  CharacterRange["a", "z"] // formatItalic,
  CharacterRange[16^^1D44E, 16^^1D467] // Append @ FromCharacterCode[16^^210E]
]


(* ::Subsection:: *)
(*Latin, bold*)


(* ::Subsubsection:: *)
(*Uppercase*)


comparisonTableForm[
  CharacterRange["A", "Z"] // formatBold,
  CharacterRange[16^^1D400, 16^^1D419]
]


(* ::Subsubsection:: *)
(*Lowercase*)


comparisonTableForm[
  CharacterRange["a", "z"] // formatBold,
  CharacterRange[16^^1D41A, 16^^1D433]
]


(* ::Subsection:: *)
(*Greek*)


(* ::Subsubsection:: *)
(*Uppercase*)


comparisonTableForm[
  CharacterRange["\[CapitalAlpha]", "\[CapitalOmega]"],
  CharacterRange[16^^00391, 16^^003A9]
]


(* ::Subsubsection:: *)
(*Lowercase*)


comparisonTableForm[
  CharacterRange["\[Alpha]", "\[Omega]"],
  CharacterRange[16^^003B1, 16^^003C9]
]


comparisonTableForm[
  {"\[CurlyTheta]", "\[Phi]", "\[CurlyPi]", "\[CurlyKappa]", "\[CurlyRho]", "", "\[Epsilon]"},
  FromCharacterCode /@ {16^^003D1, 16^^003D5, 16^^003D6, 16^^003F0, 16^^003F1, 16^^003F4, 16^^003F5}
]


(* ::Subsection:: *)
(*Greek, italic*)


(* ::Subsubsection:: *)
(*Uppercase*)


comparisonTableForm[
  CharacterRange["\[CapitalAlpha]", "\[CapitalOmega]"] // formatItalic,
  CharacterRange[16^^1D6E2, 16^^1D6FA]
]


(* ::Subsubsection:: *)
(*Lowercase*)


comparisonTableForm[
  CharacterRange["\[Alpha]", "\[Omega]"] // formatItalic,
  CharacterRange[16^^1D6FC, 16^^1D714]
]


comparisonTableForm[
  {},
  CharacterRange[16^^1D715, 16^^1D71B]
]


(* ::Subsection:: *)
(*Accidentals*)


comparisonTableForm[
  {"\[Flat]", "\[Natural]", "\[Sharp]"},
  CharacterRange[16^^0266D, 16^^0266F]
]
