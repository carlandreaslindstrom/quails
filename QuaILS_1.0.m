(* ::Package:: *)

(************************************************************************)
(***    QuaILS  -  Quadrupole-based Intelligent Lattice Solver        ***)
(************************************************************************)
(***                                                                  ***)
(***	AUTHOR                                                        ***)
(***	Carl A. Lindstr\[OSlash]m                                             ***)
(***	University of Oslo / SLAC National Accelerator Laboratory     ***)
(***	Mar 2015 - Sep 2016                                           ***)
(***                                                                  ***)
(***	DESCRIPTION:                                                  ***)
(***	Framework to find apochromatic (to your required order)       ***)
(***	quadrupole-only lattices.                                     ***)
(***                                                                  ***)
(***	DISCLAIMER / HOW TO CITE:                                     ***)
(***	If you use this framework for your academic work, please      ***)
(***	cite our article on the method implemented:                   ***)
(***                                                                  ***)
(***	"Design of general apochromatic drift-quadrupole beam lines"  ***)
(***	C. A. Lindstrom and E. Adli, PRAB 19, 071002 (2016)           ***)
(***    DOI: http://dx.doi.org/10.1103/PhysRevAccelBeams.19.071002    ***)
(***                                                                  ***)
(***    REQUIREMENTS / DEPENDENCIES                                   ***)
(***    \[Dash] MATLAB must be installed: not included (duh)                ***)
(***    \[Dash] MATLink (http://matlink.org): not included                  ***)
(***    \[Dash] ELEGANT (by M. Borland): is included                        ***)
(***    \[Dash] ToMATLAB.m library: is included                             ***)
(***                                                                  ***)
(************************************************************************)



BeginPackage["QuaILS`"];

(*** USER CONFIGURATION ***)

(* Please change "HOME" to the absolute path for this file [e.g. /Users/macuser/QuaILS] *)
HOME = "/Users/calind/PhD/QuaILS"; (* no / at the end *)



(*** ALL AVAILABLE PUBLIC FUNCTION ***)

(* Lattice definitions *)
Lattice::usage = "Lattice given alternating drifts and quadrupoles.";
ThickLattice::usage = "Define thick lens lattice.";

(* Extract info from the defined lattice *)
DS::usage = "Extract drift spaces from a lattice.";
KS::usage = "Extract quadrupole strengths from a lattice.";
LS::usage = "Lens thicknesses.";
IsQuads::usage = "Extracts whether lenses are quadrupoles (x-y asym) or capillary discharge tubes (x-y sym).";
IsThick::usage = "Returns whether lattice is thick or thin.";
LatticeLength::usage = "Length of  a given lattice.";

(* Mirror or repeat the lattice *)
SymmetricLattice::usage = "Appends a mirror image of the lattice with identical quad strengths.";
AntiSymmetricLattice::usage = "Appends a mirror image of the lattice with inversed quad strengths.";
RepLattice::usage = "Repeats a lattice periodically.";

(* Transfer R-matrix of a lattice *)
R::usage = "R-matrix for a given lattice.";

(* Useful lattice functions *)
\[Beta]x::usage = "Beta function in x-direction given a lattice and an initial beta.";
\[Beta]y::usage = "Beta function in y-direction given a lattice and an initial beta.";
Sqrt\[Beta]x::usage = "Square root of beta function in x-direction.";
Sqrt\[Beta]y::usage = "Square root of beta function in y-direction.";
\[Alpha]x::usage = "Alpha function in x-direction given a lattice and an initial beta.";
\[Alpha]y::usage = "Alpha function in y-direction given a lattice and an initial beta.";
\[Gamma]x::usage = "Gamma function in x-direction given a lattice and an initial beta.";
\[Gamma]y::usage = "Gamma function in y-direction given a lattice and an initial beta.";
\[Mu]x::usage = "Phase advance in x-direction given a lattice and an initial beta.";
\[Mu]y::usage = "Phase advance in y-direction given a lattice and an initial beta.";
Qx::usage = "Tune in x-direction.";
Qy::usage = "Tune in y-direction.";

(* Twiss parameter input/output stucture *)
TwissParams::usage = "Input/output format for twiss parameters.";
TW\[Beta]x::usage = "Extract Twiss beta in x.";
TW\[Beta]y::usage = "Extract Twiss beta in y.";
TW\[Alpha]x::usage = "Extract Twiss alpha in x.";
TW\[Alpha]y::usage = "Extract Twiss alpha in y.";
TW\[Gamma]x::usage = "Extract Twiss gamma in x.";
TW\[Gamma]y::usage = "Extract Twiss gamma in y.";

(* Chromatic Taylor expansion on the R-matrix *)
Rexpand::usage = "Chromatic expansion of R-matrix for a given lattice.";

(* Value of lattice functions vs. s *)
LatticeFunction::usage = "Piecewise lattice function versus s.";
ChromaticLatticeFunction::usage = "Piecewise chromatic lattice function versus s.";

(* Plot of lattice functions vs. s *)
PlotLatticeFunctions::usage = "Plot lattice functions along beamline distance s.";
PlotChromaticLatticeFunctions::usage = "Plot chromatic lattice function.";
PlotBeamSize3D::usage = "3D plot of actual beam size.";
PlotWvector3D::usage = "3D parametric plot of W-function in chromatic A/B space.";
PlotWvector2D::usage = "2D parametric plot of W-function in chromatic A/B space.";

(* Useful chromatic lattice functions *)
Wx::usage = "Chromatic amplitude in x-direction.";
Wy::usage = "Chromatic amplitude in y-direction.";
Ax::usage = "Chromatic A-vector in x";
Bx::usage = "Chromatic B-vector in x";
Ay::usage = "Chromatic A-vector in y";
By::usage = "Chromatic B-vector in y";
\[CapitalPhi]x::usage = "Chromatic phase in x-direction.";
\[CapitalPhi]y::usage = "Chromatic phase in y-direction.";
\[Xi]x::usage = "Phase chromaticity in x-direction.";
\[Xi]y::usage = "Phase chromaticity in x-direction.";
\[CapitalDelta]\[Epsilon]x\[Epsilon]0::usage = "Projected emittance growth in x, given some energy spread.";
\[CapitalDelta]\[Epsilon]y\[Epsilon]0::usage = "Projected emittance growth in y, given some energy spread.";

(* Emittance growths *)
EmittanceGrowth::usage = "Value of emittance growth";
PrintEmittanceGrowths::usage = "Prints emittance growths for several energy spreads.";
PlotEmittanceGrowth::usage = "Plot emittance growth vs. energy spread.";

(* Luminosity for final focus systems *)
PlotRelativeLuminosity::usage = "Plots relative luminosity as a function of energy spread.";

(* Visualizing the chromatic focusing dependence of the lattice *) 
PlotChromaticDependence::usage = "Plot and print alpha and beta chromatic expansions.";
PrintChromaticExpansion::usage = "Print the chromatic Taylor series of the Twiss parameters";

(* Looking at the lattice with different energies *)
TriBeamSize::usage = "Function for illustration in a paper";
LatticeFigure::usage = "Shows a MADX-style diagram of the lattice.";

(* Visualizing phase space ellipses *)
PhaseSpaceEllipse::usage = "Phase space ellipse.";

(* MATLAB interface used for numerical solving *)
WriteMatlabConstraints::usage = "Write matching constraints to Matlab file.";
WriteMatlabMeritFunction::usage = "Write matching merit function to Matlab file.";
WriteMatlabMatchingParams::usage = "Write matching parameters.";
WriteMatlabInitialGuess::usage = "Write initial guess.";
ReadMatlabSolutions::usage = "Read matching solutions from Matlab file.";
RunMatlabMinimizer::usage = "Run Matlab minimizer for matching with constraints and merit function (interfaces using MatLink).";

(* ELEGANT interface for particle tracking [for verification] *)
TrackElegantLattice::usage = "Tracks particles in Elegant via Matlab.";
TrackSingleParticle::usage = "Track single particle in Elegant.";



Begin["QuaILS`Private`"];


(*** LATTICES UTILITIES ***)

(* Lattice definition: thin and thick lenses *)
Lattice[ds_,kls_,quads_:True] := {ds,kls,quads}; (* thin-lens lattice *)
ThickLattice[ds_, ks_, ls_, quads_:True] := {ds, ks, ls, quads}; (* thick-lens lattice *)

(* Extracting lattice parameters *)
DS[lat_] := lat[[1]]; (* extract drift space lengths from a lattice *)
KS[lat_] := lat[[2]]; (* extract quad strengths from a lattice *)
LS[lat_] := If[IsThick[lat],lat[[3]],{}]; (* extract quad lengths from a lattice *)
IsQuads[lat_] := lat[[Length[lat]]]; (* is the lattice quads or radial lenses? *)
IsThick[lat_] := (Length[lat]==4); (* are the lenses used thick or thin? *)

(* Total length of the lattice *)
LatticeLength[lat_] := Total[DS[lat]] + If[IsThick[lat],Total[LS[lat]],0];

(* Mirror symmetric lattice *)
SymmetricLattice[lat_] := Module[{ds=DS[lat], ks=KS[lat], ls=LS[lat], isQ=IsQuads[lat], dsym, ksym, lsym},
	If[Length[ds]==Length[ks],{
		dsym = Join[ds,Reverse[ds]];
		If[IsThick[lat],{
			ksym = Join[Drop[ks,-1],{Last[ks]},Reverse[Drop[ks,-1]]],
			lsym = Join[Drop[ls,-1],{2*Last[ls]},Reverse[Drop[ls,-1]]];
		},{
			ksym = Join[Drop[ks,-1],{2*Last[ks]},Reverse[Drop[ks,-1]]]
		}];
	},{
		dsym = Join[Drop[ds,-1],{2*Last[ds]},Reverse[Drop[ds,-1]]];
		ksym = Join[ks,Reverse[ks]];
		lsym = Join[ls,Reverse[ls]];
	}];
	If[IsThick[lat],
		ThickLattice[dsym,ksym,lsym,isQ],
		Lattice[dsym,ksym,isQ]
	]
];

(* Mirror anti-symmetric lattice [polarity switched, only use for quads] *)
AntiSymmetricLattice[lat_] := Module[{ds=DS[lat], ks=KS[lat], ls=LS[lat], isQ=IsQuads[lat], dsym, ksym, lsym},
	If[Length[ds]==Length[ks],{
		dsym = Join[ds,{0},Reverse[ds]];
	},{
		dsym = Join[Drop[ds,-1],{2*Last[ds]},Reverse[Drop[ds,-1]]];
	}];
	ksym = Join[ks,-Reverse[ks]];
	lsym = Join[ls,Reverse[ls]];
	If[IsThick[lat],
		ThickLattice[dsym,ksym,lsym,isQ],
		Lattice[dsym,ksym,isQ]
	]
];

(* Periodic repetition of identical lattices *)
RepLattice[lat_, n_] := Module[{ds=DS[lat],ks,ls},
	ds = Flatten[{ds[[1]],Riffle[ConstantArray[Drop[Drop[ds,1],-1],n],ds[[1]]+ds[[-1]]],ds[[-1]]}];
	ks = Flatten[ConstantArray[KS[lat],n]];
	ls = Flatten[ConstantArray[LS[lat],n]];
	If[IsThick[lat],
		Lattice[ds,ks,ls,IsQuads[lat]],
		Lattice[ds,ks,IsQuads[lat]]
	]
];
	


(*** TRANSFER R-MATRICES ***)
           
(* quad lattice R-matrix *)
R[lat_, \[Delta]_:0]:= Apply[Dot,Reverse[Riffle[
	Rd/@DS[lat], 
	If[IsThick[lat],
		If[IsQuads[lat],RQf,RCDf]@@@Transpose[{KS[lat]/(1+\[Delta]),LS[lat]}],
		If[IsQuads[lat],Rq,Rcd]/@(KS[lat]/(1+\[Delta]))]
	]
]];

(* R-matrix element [i,j] *)
Rij[M_, i_, j_] := M[[i,j]];

(* INTERNAL \[Dash] Drift space R-matrix *)
Rd[L_]:={{1,L,0,0},{0,1,0,0},{0,0,1,L},{0,0,0,1}};

(* INTERNAL \[Dash] thin quad R-matrix *)
Rq[K_] := {{1,0,0,0},{-K,1,0,0},{0,0,1,0},{0,0,K,1}};

(* INTERNAL \[Dash] Thin capillary R-matrix *)
Rcd[K_]:={{1,0,0,0},{-K,1,0,0},{0,0,1,0},{0,0,-K,1}};

(* INTERNAL\[NonBreakingSpace]\[Dash]\[NonBreakingSpace]Thick quadrupole R-matrix *)
RQf[K_, L_] := {{Cos[Sqrt[K]L], Sin[Sqrt[K]L]/Sqrt[K], 0, 0},
	           {-Sqrt[K] Sin[Sqrt[K]L], Cos[Sqrt[K]L], 0, 0},
	           {0, 0, Cosh[Sqrt[K]L], Sinh[Sqrt[K]L]/Sqrt[K]},
	           {0, 0, Sqrt[K] Sinh[Sqrt[K]L], Cosh[Sqrt[K]L]}};

(* INTERNAL\[NonBreakingSpace]\[Dash]\[NonBreakingSpace]Thick capillary R-matrix *)
RCDf[K_, L_] := {{Cos[Sqrt[K]L], Sin[Sqrt[K]L]/Sqrt[K], 0, 0},
	           {-Sqrt[K] Sin[Sqrt[K]L], Cos[Sqrt[K]L], 0, 0},
	           {0, 0, Cos[Sqrt[K]L], Sin[Sqrt[K]L]/Sqrt[K]},
	           {0, 0, -Sqrt[K] Sin[Sqrt[K]L], Cos[Sqrt[K]L]}};



(*** TWISS PARAMETERS ***)

(* Twiss parameter \[Beta] as well as Sqrt[\[Beta]] for beam size *)
\[Beta]x[M_, input_] := Rij[M,1,1]^2*TW\[Beta]x[input] - 2*Rij[M,1,1]*Rij[M,1,2]*TW\[Alpha]x[input] + Rij[M,1,2]^2*TW\[Gamma]x[input];
\[Beta]y[M_, input_] := Rij[M,3,3]^2*TW\[Beta]y[input] - 2*Rij[M,3,3]*Rij[M,3,4]*TW\[Alpha]y[input] + Rij[M,3,4]^2*TW\[Gamma]y[input];
Sqrt\[Beta]x[M_, input_] := Sqrt[\[Beta]x[M, input]];
Sqrt\[Beta]y[M_, input_] := Sqrt[\[Beta]y[M, input]];

(* Twiss parameter \[Alpha] *)
\[Alpha]x[M_, input_] := - Rij[M,1,1]*Rij[M,2,1]*TW\[Beta]x[input] + (Rij[M,1,1]*Rij[M,2,2]+R12[R]*Rij[M,2,1])*TW\[Alpha]x[input] - Rij[M,1,2]Rij[M,2,2]*TW\[Gamma]x[input];
\[Alpha]y[M_, input_] := - Rij[M,3,3]*Rij[M,4,3]*TW\[Beta]y[input] + (Rij[M,3,3]*Rij[M,4,4]+Rij[M,3,4]*Rij[M,4,3])*TW\[Alpha]y[input] - Rij[M,3,4]Rij[M,4,4]*TW\[Gamma]y[input];

(* Twiss parameter \[Gamma] = (1-\[Alpha]^2)/\[Beta] *)
\[Gamma]x[M_, input_] := Rij[M,2,1]^2*TW\[Beta]x[input] - 2*Rij[M,2,1]*Rij[M,2,2]*TW\[Alpha]x[input] + Rij[M,2,2]^2*TW\[Gamma]x[input];
\[Gamma]y[M_, input_] := Rij[M,4,3]^2*TW\[Beta]y[input] - 2*Rij[M,4,3]*Rij[M,4,4]*TW\[Alpha]y[input] + Rij[M,4,4]^2*TW\[Gamma]y[input];

(* Phase advance *)
\[Mu]x[M_, input_] := Mod[ArcTan[(Rij[M,1,1]*TW\[Beta]x[input] - Rij[M,1,2]*TW\[Alpha]x[input]),Rij[M,1,2]], 2Pi];
\[Mu]y[M_, input_] := Mod[ArcTan[(Rij[M,3,3]*TW\[Beta]y[input] - Rij[M,3,4]*TW\[Alpha]y[input]),Rij[M,3,4]], 2Pi];

(* Tune [normalized by 2\[Pi]] *)
Qx[M_, input_] := ArcTan[(Rij[M,1,1]*TW\[Beta]x[input] - Rij[M,1,2]*TW\[Alpha]x[input]),Rij[M,1,2]]/(2Pi);
Qy[M_, input_] := ArcTan[(Rij[M,3,3]*TW\[Beta]y[input] - Rij[M,3,4]*TW\[Alpha]y[input]),Rij[M,3,4]]/(2Pi);



(*** TWISS FORMAT ***)

(* Twiss parameter definition [input/outputs] *)
TwissParams[\[Beta]x_, \[Beta]y_, \[Alpha]x_:0, \[Alpha]y_:0] := {{\[Beta]x,\[Beta]y},{\[Alpha]x,\[Alpha]y}};

(* Extracting specific Twiss parameters *)
TW\[Beta]x[twiss_] := twiss[[1,1]];
TW\[Beta]y[twiss_] := twiss[[1,2]];
TW\[Alpha]x[twiss_] := twiss[[2,1]];
TW\[Alpha]y[twiss_] := twiss[[2,2]];
TW\[Gamma]x[twiss_] := (1+TW\[Alpha]x[twiss]^2)/TW\[Beta]x[twiss];
TW\[Gamma]y[twiss_] := (1+TW\[Alpha]y[twiss]^2)/TW\[Beta]y[twiss];



(*** CHROMATIC EXPANSION ***)

(* Series expansion of the R-matrix to desired order *)
Rexpand[lat_, \[Delta]_, order_] := Series[R[lat,\[Delta]],{\[Delta],0,order}];



(*** EMITTANCE GROWTH ***)

(* Emittance growth value in [x,y] *)
EmittanceGrowth[lat_, input_, \[Sigma]E_, order_:10] := Module[{Rlat, \[Delta], \[Mu]},
	Rlat = Rexpand[lat, \[Delta], order];
	\[Mu][g_] := Sum[SeriesCoefficient[g[Rlat, input[\[Delta]]],{\[Delta],0,n}]*\[Sigma]E^n*(n-1)!!,{n,0,order,2}];
	Return[{Sqrt[\[Mu][\[Beta]x]*\[Mu][\[Gamma]x]-\[Mu][\[Alpha]x]^2]-1, Sqrt[\[Mu][\[Beta]y]*\[Mu][\[Gamma]y]-\[Mu][\[Alpha]y]^2]-1}];
];

(* Printout of relative emittance growth for many energy spreads *)
PrintEmittanceGrowths[lat_, input_, \[Sigma]Es_] := Module[{\[CapitalDelta]\[Epsilon]xs, \[CapitalDelta]\[Epsilon]ys},
	{\[CapitalDelta]\[Epsilon]xs, \[CapitalDelta]\[Epsilon]ys} = EmittanceGrowth[lat, input, \[Sigma]Es, 20];
	If[IsQuads[lat],
		Do[Print["Energy spread ", N[\[Sigma]Es[[i]]*100,2], "% rms \[Rule] \[CapitalDelta]\[Epsilon]/\[Epsilon] = (", Re[N[\[CapitalDelta]\[Epsilon]xs[[i]], 3]], ", ", Re[N[\[CapitalDelta]\[Epsilon]ys[[i]], 3]], ")"], {i,1,Length[\[Sigma]Es]}],
		Do[Print["Energy spread ", N[\[Sigma]Es[[i]]*100,2], "% rms \[Rule] \[CapitalDelta]\[Epsilon]/\[Epsilon] = ", Re[N[\[CapitalDelta]\[Epsilon]xs[[i]], 3]]], {i,1,Length[\[Sigma]Es]}]
	];
];

(* Plot of relative emittance growth vs. energy spread *)
PlotEmittanceGrowth[lat_, input_, \[Sigma]Emin_:1*^-5, \[Sigma]Emax_:0.1, order_:8] := Module[{\[Delta], Rlat, \[Mu], \[CapitalDelta]\[Epsilon]x, \[CapitalDelta]\[Epsilon]y, opts},
	Rlat = Rexpand[lat, \[Delta], order];
	\[Mu][g_,\[Sigma]E_] := Sum[SeriesCoefficient[g[Rlat, input[\[Delta]]],{\[Delta],0,n}]*\[Sigma]E^n*(n-1)!!,{n,0,order,2}];
	opts = Sequence[PlotLabel->"Relative emittance growth vs. energy spread", FrameLabel->{"Energy spread (rms), \!\(\*SubscriptBox[\(\[Sigma]\), \(E\)]\)", "Relative emittance growth, \[CapitalDelta]\[Epsilon]/\[Epsilon]"}, ImageSize->Medium, Frame->True, GridLinesStyle->LightGray, GridLines->Full];
	If[IsQuads[lat],
		LogLogPlot[{Sqrt[\[Mu][\[Beta]x,\[Sigma]E]*\[Mu][\[Gamma]x,\[Sigma]E]-\[Mu][\[Alpha]x,\[Sigma]E]^2]-1, Sqrt[\[Mu][\[Beta]y,\[Sigma]E]*\[Mu][\[Gamma]y,\[Sigma]E]-\[Mu][\[Alpha]y,\[Sigma]E]^2]-1}, {\[Sigma]E,\[Sigma]Emin,\[Sigma]Emax}, Evaluate[opts], PlotLegends->{"x", "y"}],
		LogLogPlot[{Sqrt[\[Mu][\[Beta]x,\[Sigma]E]*\[Mu][\[Gamma]x,\[Sigma]E]-\[Mu][\[Alpha]x,\[Sigma]E]^2]-1}, {\[Sigma]E,\[Sigma]Emin,\[Sigma]Emax}, Evaluate[opts]]
	]
];



(*** LUMINOSITY for FINAL FOCUS STUDIES ***)

(* Plot the relative luminosity vs energy spread *)
PlotRelativeLuminosity[lat_, input_, output_, \[Delta]max_:0.05, order_:10] := Module[{\[Delta],Rlat,\[Sigma]Emax=\[Delta]max/1.5,\[Delta]s,\[Sigma]Es,\[CapitalDelta]L\[Delta]={},\[CapitalDelta]L\[Sigma]E={},\[Mu],\[Mu]rms,opts,\[Mu]\[Beta]xs,\[Mu]\[Beta]ys,\[Mu]\[Beta]xrms,\[Mu]\[Beta]yrms,acceptance},
	Rlat = Normal[Rexpand[lat, \[Delta], order]];
	
	\[Delta]s = Range[-\[Delta]max(1-1*^-6),\[Delta]max,\[Delta]max/100];	
	\[Mu][\[Beta]_]:=Sum[SeriesCoefficient[\[Beta][Rlat, input],{\[Delta],0,n}]*\[Delta]s^n,{n,0,order}];
	{\[Mu]\[Beta]xs, \[Mu]\[Beta]ys} = {\[Mu][\[Beta]x], \[Mu][\[Beta]y]};
	Do[AppendTo[\[CapitalDelta]L\[Delta],{\[Delta]s[[i]],Sqrt[TW\[Beta]x[output]TW\[Beta]y[output]/(\[Mu]\[Beta]xs[[i]]\[Mu]\[Beta]ys[[i]])]}],{i,1,Length[\[Delta]s]}];

	\[Sigma]Es = Range[\[Sigma]Emax*(1*^-6),\[Sigma]Emax,\[Sigma]Emax/100];
	\[Mu]rms[\[Beta]_] := Sum[SeriesCoefficient[\[Beta][Rlat, input],{\[Delta],0,n}]*\[Sigma]Es^n*(n-1)!!,{n,0,order,2}];
	{\[Mu]\[Beta]xrms,\[Mu]\[Beta]yrms} = {\[Mu]rms[\[Beta]x], \[Mu]rms[\[Beta]y]};
	Do[AppendTo[\[CapitalDelta]L\[Sigma]E,{\[Sigma]Es[[i]],Sqrt[TW\[Beta]x[output]TW\[Beta]y[output]/(\[Mu]\[Beta]xrms[[i]]\[Mu]\[Beta]yrms[[i]])]}],{i,1,Length[\[Sigma]Es]}];
	
	opts = Sequence[PlotRange->{0,All}, Frame->True, PlotRangePadding->0, ImageSize->Medium, AspectRatio->0.6];
	Grid[{{
		ListLinePlot[\[CapitalDelta]L\[Sigma]E,FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(E\)]\)","\[GothicCapitalL]/\!\(\*SubscriptBox[\(\[GothicCapitalL]\), \(0\)]\)"}, PlotLabel->"Relative luminosity vs. rms energy spread", Evaluate[opts]],
		ListLinePlot[\[CapitalDelta]L\[Delta],FrameLabel->{"\[Delta]","\[GothicCapitalL]/\!\(\*SubscriptBox[\(\[GothicCapitalL]\), \(0\)]\)"}, PlotLabel->"Relative luminosity vs. energy offset", Evaluate[opts]]
	}}]
];



(*** LATTICE FUNCTIONS ALONG LATTICE: f(s) ***)

(* VALUE of non-chromatic lattice functions [e.g. Twiss params] along s *)
LatticeFunction[tlat_, input_, f_, \[Delta]0_:0] := Module[{ctlats, g, Ds, Rs},
	ctlats = CumulLattices[tlat];
	Ds = LatticeLength /@ ctlats;
	Rs = R[#,\[Delta]0]& /@ Drop[ctlats, -1];
	g[s_] := If[IsThick[tlat],
		Piecewise[{If[Length[DS[ctlats[[#]]]]==Length[KS[ctlats[[#]]]],
			f[Rd[s-Ds[[#]]].Rs[[#]], input[]],
			f[If[IsQuads[tlat],RQf,RCDf][Last[KS[ctlats[[#+1]]]]/(1+\[Delta]0),s-Ds[[#]]].Rs[[#]], input[]]
		],(Ds[[#]]<=s<=Ds[[#+1]])}& /@ Range[Length[Rs]]],
		Piecewise[{f[Rd[s-Ds[[#]]].Rs[[#]], input[]],(Ds[[#]]<=s<=Ds[[#+1]])}& /@ Range[Length[Rs]]]
	];
	Return[g];
];

(* VALUE of chromatic lattice functions [e.g. W-function] along s *)
ChromaticLatticeFunction[tlat_, input_, f_, order_] := Module[{\[Delta], ctlats, g, Ds, Rs},
	ctlats = CumulLattices[tlat];
	Ds = LatticeLength /@ ctlats;
	Rs = Rexpand[#,\[Delta],order]& /@ Drop[ctlats,-1];
	g[s_] := If[IsThick[tlat],
		Piecewise[{If[Length[DS[ctlats[[#]]]]==Length[KS[ctlats[[#]]]],
			f[Rd[s-Ds[[#]]].Rs[[#]], input, \[Delta]],
			f[If[IsQuads[tlat],RQf,RCDf][Series[Last[KS[ctlats[[#+1]]]]/(1+\[Delta]),{\[Delta],0,order}],s-Ds[[#]]].Rs[[#]], input, \[Delta]]
		],(Ds[[#]]<=s<=Ds[[#+1]])}& /@ Range[Length[Rs]]],
		Piecewise[{f[Rd[s-Ds[[#]]].Rs[[#]], input, \[Delta]],(Ds[[#]]<=s<=Ds[[#+1]])}& /@ Range[Length[Rs]]]
	];
	Return[g];
];

(* INTERNAL/PRIVATE\[NonBreakingSpace]\[Dash]\[NonBreakingSpace]Cumulative sub-lattices *)
CumulLattices[lat_] := If[IsThick[lat], 
	Riffle[
		Table[ThickLattice[Take[DS[lat],n],Take[KS[lat],n],Take[LS[lat],n],IsQuads[lat]],{n,0,Length[KS[lat]]}],
		Table[ThickLattice[Take[DS[lat],n+1],Take[KS[lat],n],Take[LS[lat],n],IsQuads[lat]],{n,0,Length[KS[lat]]}]
	],
	Append[Table[Lattice[Take[DS[lat],n],Take[KS[lat],n],IsQuads[lat]],{n,0,Length[KS[lat]]}],lat]
];



(*** USEFUL LATTICE FUNCTIONS ***)

(* Chromatic amplitude in x/y *)
Wx[R_, input_, \[Delta]_] := Sqrt[Ax[R, input, \[Delta]]^2+Bx[R, input, \[Delta]]^2]; 
Wy[R_, input_, \[Delta]_] := Sqrt[Ay[R, input, \[Delta]]^2+By[R, input, \[Delta]]^2];

(* Chromatic phase in x/y [normalized by 2\[Pi] such that \[CapitalDelta]\[CapitalPhi]=1 is one oscillation] *)
\[CapitalPhi]x[R_, input_, \[Delta]_] := ArcTan[N[Bx[R, input, \[Delta]]], Ax[R, input, \[Delta]]]/(2Pi);
\[CapitalPhi]y[R_, input_, \[Delta]_] := ArcTan[N[By[R, input, \[Delta]]], Ay[R, input, \[Delta]]]/(2Pi);

(* Chromatic components A and B [where W^2 = A^2 + B^2] *)
Bx[R_, input_, \[Delta]_] := Module[{f\[Beta]x = \[Beta]x[R, input[\[Delta]]]}, SeriesCoefficient[f\[Beta]x,{\[Delta],0,1}]/SeriesCoefficient[f\[Beta]x,{\[Delta],0,0}]];
By[R_, input_, \[Delta]_] := Module[{f\[Beta]y = \[Beta]y[R, input[\[Delta]]]}, SeriesCoefficient[f\[Beta]y,{\[Delta],0,1}]/SeriesCoefficient[f\[Beta]y,{\[Delta],0,0}]];
Ax[R_, input_, \[Delta]_] := Module[{f\[Beta]x = \[Beta]x[R, input[\[Delta]]], f\[Alpha]x = \[Alpha]x[R, input[\[Delta]]], bx},
	bx = SeriesCoefficient[f\[Beta]x,{\[Delta],0,1}]/SeriesCoefficient[f\[Beta]x,{\[Delta],0,0}];
	Return[SeriesCoefficient[f\[Alpha]x,{\[Delta],0,1}] - bx*SeriesCoefficient[f\[Alpha]x,{\[Delta],0,0}]];
];
Ay[R_, input_, \[Delta]_] := Module[{f\[Beta]y = \[Beta]y[R, input[\[Delta]]], f\[Alpha]y = \[Alpha]y[R, input[\[Delta]]], by},
	by = SeriesCoefficient[f\[Beta]y,{\[Delta],0,1}]/SeriesCoefficient[f\[Beta]y,{\[Delta],0,0}];
	Return[SeriesCoefficient[f\[Alpha]y,{\[Delta],0,1}] - by*SeriesCoefficient[f\[Alpha]y,{\[Delta],0,0}]];
];

(* Chromaticity in x/y *)
\[Xi]x[R_, input_, \[Delta]_] := SeriesCoefficient[Qx[R, input[\[Delta]]],{\[Delta],0,1}];
\[Xi]y[R_, input_, \[Delta]_] := SeriesCoefficient[Qy[R, input[\[Delta]]],{\[Delta],0,1}];

(* Relative projected emittance growth in x/y *)
\[CapitalDelta]\[Epsilon]x\[Epsilon]0[R_, input_, \[Delta]_, \[Sigma]E_, order_:6] := Module[{\[Mu]},
	\[Mu][g_] := Sum[SeriesCoefficient[g[R, input],{\[Delta],0,n}]*\[Sigma]E^n*(n-1)!!,{n,0,order,2}];
	Return[Sqrt[\[Mu][\[Beta]x]*\[Mu][\[Gamma]x]-\[Mu][\[Alpha]x]^2]-1];
];
\[CapitalDelta]\[Epsilon]y\[Epsilon]0[R_, input_, \[Delta]_, \[Sigma]E_, order_:6] := Module[{\[Mu]},
	\[Mu][g_] := Sum[SeriesCoefficient[g[R, input],{\[Delta],0,n}]*\[Sigma]E^n*(n-1)!!,{n,0,order,2}];
	Return[Sqrt[\[Mu][\[Beta]y]*\[Mu][\[Gamma]y]-\[Mu][\[Alpha]y]^2]-1];
];



(* PLOTTING LATTICE FUNCTIONS *)

(* Plot of non-chromatic lattice functions [e.g. Twiss params] along s *)
PlotLatticeFunctions[lat_, input_, fs_, \[Delta]0_:0, range_:All] := Module[{Lfs, labels, opts},
	Lfs = LatticeFunction[lat, input, #, \[Delta]0]& /@ fs;
	labels = StringSplit[ToString/@fs,"["][[All,1]];
	opts = Sequence[Exclusions->None, PlotRange->range, AxesLabel->{"s / m",labels}, PlotLegends->{labels}, ImageSize->{600,200}, AspectRatio->1/3];
	Return[Plot[Evaluate[#[s]& /@ Lfs],{s,0,LatticeLength[lat]}, Evaluate[opts]]];
];

(* Plot of chromatic lattice functions [e.g. W-function] along s *)
PlotChromaticLatticeFunctions[lat_, input_, fs_, order_:3] := Module[{Lfs, labels, opts},
	Lfs = ChromaticLatticeFunction[lat, input, #, order]& /@ fs;
	labels = StringSplit[ToString/@fs,"["][[All,1]];
	opts = Sequence[Exclusions->None, PlotRange->All, AxesLabel->{"s",labels}, PlotLegends->{labels}, ImageSize->{600,200}, AspectRatio->1/3];
	Return[Plot[Evaluate[#[s]& /@ Lfs],{s,0,LatticeLength[lat]}, Evaluate[opts]]];
];

(* 3D plot of beam size Sqrt[\[Beta]] in x/y along s *)
PlotBeamSize3D[lat_, input_] := Module[{sq\[Beta]x, sq\[Beta]y, opts},
	{sq\[Beta]x, sq\[Beta]y} = LatticeFunction[lat, input, #]& /@ {Sqrt\[Beta]x, Sqrt\[Beta]y};
	opts = Sequence[Exclusions->None, BoxRatios->{3, 1, 1}, ImageSize->Medium, PlotLabel->"Beam size evolution", AxesLabel->{"s","\!\(\*SqrtBox[SubscriptBox[\(\[Beta]\), \(x\)]]\)","\!\(\*SqrtBox[SubscriptBox[\(\[Beta]\), \(y\)]]\)"}, PlotStyle->ColorData[97][[2]],PerformanceGoal->"Quality"];
	ParametricPlot3D[{s, Cos[\[Theta]]*sq\[Beta]x[s], Sin[\[Theta]]*sq\[Beta]y[s]},{s,0,LatticeLength[lat]},{\[Theta],0,2Pi}, Evaluate[opts]]
];

(* 3D plot of the chromatic vector [W and \[CapitalPhi]] *)
PlotWvector3D[lat_, input_] := Module[{\[Delta], ax, bx, ay, by, opts},
	{ax,bx,ay,by} = ChromaticLatticeFunction[lat,input,#,1]& /@ {Ax, Bx, Ay, By};
	opts = Sequence[Exclusions->None, BoxRatios->{3, 1, 1}, ImageSize->Medium, AxesLabel->{"s","A","B"}, PlotLabel->"W-vector (3D)"];
	If[IsQuads[lat],
		ParametricPlot3D[{{s, bx[s], ax[s]},{s,by[s], ay[s]}}, {s,0,LatticeLength[lat]}, Evaluate[opts]],
		ParametricPlot3D[{s, bx[s], ax[s]}, {s,0,LatticeLength[lat]}, Evaluate[opts], PlotLegends->{"x","y"}]
	]
];

(* 2D plot of the chromatic vector [W and \[CapitalPhi]] *)
PlotWvector2D[lat_, input_] := Module[{\[Delta], ax, bx, ay, by, opts},
	{ax,bx,ay,by} = ChromaticLatticeFunction[lat,input,#,1]& /@ {Ax, Bx, Ay, By};
	opts = Sequence[Exclusions->None, AspectRatio->1, ImageSize->Medium, AxesLabel->{"B","A"}, PlotLabel->"W-vector (2D)"];
	If[IsQuads[lat],
		ParametricPlot[{{bx[s], ax[s]},{by[s], ay[s]}}, {s,0,LatticeLength[lat]}, Evaluate[opts], PlotLegends->{"x","y"}],
		ParametricPlot[{bx[s], ax[s]}, {s,0,LatticeLength[lat]}, Evaluate[opts]]
	]
];



(* PLOTTING CHROMATIC EXPANSIONS *)

(* Plot the chromatic dependence of Twiss parameters *)
PlotChromaticDependence[lat_, input_, output_, \[Delta]max_:0.03, order_:4] := Module[{Rlat, \[Delta], f\[Beta]x, f\[Beta]y, f\[Alpha]x, f\[Alpha]y, opts},
	Rlat = Normal[Rexpand[lat, \[Delta], order]];
	
	opts = Sequence[FrameLabel->{"\[Delta]","\[Beta]/\[Beta]*, \[Alpha]"}, PlotLabel->"Twiss chromatic dependence", PlotRange->All, AxesOrigin->{0,0}, Frame->True, PlotRangePadding->0, ImageSize->Medium, Exclusions->None, PlotStyle->{ColorData[97,1],{ColorData[97,1],Dotted},ColorData[97,2],{ColorData[97,2],Dotted}}];	
	
	f\[Beta]x = \[Beta]x[Rlat, input[\[Delta]]]/TW\[Beta]x[output[]];
	f\[Alpha]x = \[Alpha]x[Rlat, input[\[Delta]]];
	If[IsQuads[lat],
		f\[Beta]y = \[Beta]y[Rlat, input[\[Delta]]] / TW\[Beta]y[output[]];
		f\[Alpha]y = \[Alpha]y[Rlat, input[\[Delta]]];
		Plot[{f\[Beta]x,f\[Alpha]x,f\[Beta]y,f\[Alpha]y},{\[Delta],-\[Delta]max,\[Delta]max}, Evaluate[opts], PlotLegends->{"\!\(\*SubscriptBox[\(\[Beta]\), \(x\)]\)/\[Beta]*","\!\(\*SubscriptBox[\(\[Alpha]\), \(y\)]\)","\!\(\*SubscriptBox[\(\[Beta]\), \(y\)]\)/\[Beta]*","\!\(\*SubscriptBox[\(\[Alpha]\), \(x\)]\)"}]
	,
		Plot[{f\[Beta]x,f\[Alpha]x},{\[Delta],-\[Delta]max,\[Delta]max}, Evaluate[opts], PlotLegends->{"\[Beta]/\[Beta]*","\[Alpha]"}]
	]
];

(* Print the chromatic Taylor series of the Twiss parameters *)
PrintChromaticExpansion[lat_, input_, output_, \[Delta]max_:0.03, order_:4] := Module[{Rlat, \[Delta], f\[Beta]x, f\[Beta]y, f\[Alpha]x, f\[Alpha]y},
	Rlat = Normal[Rexpand[lat, \[Delta], order]];
	
	f\[Beta]x = \[Beta]x[Rlat, input[\[Delta]]]/TW\[Beta]x[output[]];
	f\[Alpha]x = \[Alpha]x[Rlat, input[\[Delta]]];
	Print["\[Beta]x/\[Beta]x* = ", Quiet[Series[f\[Beta]x,{\[Delta],0,order}]/.\[Delta]->"\[Delta]"]];
	Print["\[Alpha]x = ", Quiet[Series[f\[Alpha]x,{\[Delta],0,order}]/.\[Delta]->"\[Delta]"]];
	If[IsQuads[lat],{
		f\[Beta]y = \[Beta]y[Rlat, input[\[Delta]]] / TW\[Beta]y[output[]];
		f\[Alpha]y = \[Alpha]y[Rlat, input[\[Delta]]];
		Print["\[Beta]y/\[Beta]y* = ", Quiet[Series[f\[Beta]y,{\[Delta],0,order}]/.\[Delta]->"\[Delta]"]];
		Print["\[Alpha]y = ", Quiet[Series[f\[Alpha]y,{\[Delta],0,order}]/.\[Delta]->"\[Delta]"]];}
	];
];

(* Beam size Sqrt[\[Beta]] for three different energies *)
TriBeamSize[lat_, input_, \[Delta]offset_, order_:4] := Module[{Rlat,\[Delta],Lrx,Lgx,Lbx,Lry,Lgy,Lby,opts1,opts1y,opts1xy,opts2,col,colstyle,off},
	Rlat = Normal[Rexpand[lat, \[Delta], order]];
	
	off=ToString[Round[\[Delta]offset*100]];
	col = {ColorData[67,8],ColorData[3,4],ColorData[3,6]}; (* red, green, blue *) 
	colstyle = {col[[1]],col[[2]],col[[3]],{Dotted,col[[1]]},{Dotted,col[[2]]},{Dotted,col[[3]]}};
	opts1 = Sequence[PlotLegends->{"\[Delta] = -"<>off<>"%","\[Delta] = 0","\[Delta] = +"<>off<>"%"}];
	opts1xy = Sequence[PlotLegends->{"x, \[Delta] = -"<>off<>"%","x, \[Delta] = 0","x, \[Delta] = +"<>off<>"%","y, \[Delta] = -"<>off<>"%","y, \[Delta] = 0","y, \[Delta] = +"<>off<>"%"}];
	opts2 = Sequence[AspectRatio->1/2.3, ImageSize->Large, PlotRange->{0,All}, FrameLabel->{"s","\!\(\*SqrtBox[\(\[Beta]\)]\)"}, PlotLabel->"Beam size \!\(\*SqrtBox[\(\[Beta]\)]\) vs s for 3 energies", PlotStyle->colstyle, Exclusions->None, Frame->True, PlotRangePadding->0];
	
	Lrx = LatticeFunction[lat, input, #, -\[Delta]offset]& /@ {Sqrt\[Beta]x}; (* red *)
	Lgx = LatticeFunction[lat, input, #, 0]& /@ {Sqrt\[Beta]x}; (* green *)
	Lbx = LatticeFunction[lat, input, #, \[Delta]offset]& /@ {Sqrt\[Beta]x}; (* blue *)
	If[IsQuads[lat],
		Lry = LatticeFunction[lat, input, #, -\[Delta]offset]& /@ {Sqrt\[Beta]y}; (* red *)
		Lgy = LatticeFunction[lat, input, #, 0]& /@ {Sqrt\[Beta]y}; (* green *)
		Lby = LatticeFunction[lat, input, #, \[Delta]offset]& /@ {Sqrt\[Beta]y}; (* blue *)
		Plot[{Evaluate[#[s]& /@ Lrx],Evaluate[#[s]& /@ Lgx],Evaluate[#[s]& /@ Lbx],Evaluate[#[s]& /@ Lry],Evaluate[#[s]& /@ Lgy],Evaluate[#[s]& /@ Lby]},{s,0,LatticeLength[lat]}, Evaluate[opts1xy], Evaluate[opts2]]
	,
		Plot[{Evaluate[#[s]& /@ Lrx],Evaluate[#[s]& /@ Lgx],Evaluate[#[s]& /@ Lbx]},{s,0,LatticeLength[lat]}, Evaluate[opts1], Evaluate[opts2]]
	]
];

(* Generates an image representing the lattice, similar to MADX *)
LatticeFigure[lat_] := Module[{beamline,fig,height=55*2,length=435*2,l,h,isFocusing,s0,l0,box,line},
	l = LatticeLength[lat];	
	h = l*height/(3*length);
	beamline = Line[{{0,h},{l,h}}];
	fig = {Black,Thick,beamline,Black,Thick,beamline};
	If[IsThick[lat],
		Do[
			isFocusing = -Sign[KS[lat][[i]]];
			s0 = Total[DS[lat][[1;;i]]]+Total[LS[lat][[1;;(i-1)]]];
			l0 = LS[lat][[i]];
			box = Rectangle[{s0,h},{s0+l0, h*(1-isFocusing)}];
			fig = Prepend[fig,{EdgeForm[{Thin,Black}],FaceForm[White],box}];
			,{i,1,Length[KS[lat]]}
		],
		Do[
			isFocusing = Sign[KS[lat][[i]]];
			s0 = Total[DS[lat][[1;;i]]];
				line = Line[{{s0,h},{s0, h*(1-isFocusing)}}];
			fig = Prepend[fig,{line,Black,Thick}];
			,{i,1,Length[KS[lat]]}
		]
	];
	Graphics[fig, ImageSize->{length,height}]
];

(* Plots the phase space ellipse(s) at a given location in the lattice *)
(* TODO: something is wrong here, giving non-elliptical results *)
PhaseSpaceEllipse[lat_, input_, \[Delta]_:0, s_:Infinity] := Module[{opts, L = LatticeLength[lat], f\[Beta]x, f\[Alpha]x, f\[Gamma]x, f\[Beta]y, f\[Alpha]y, f\[Gamma]y, bb,a, cc, f, size = 1.3},
	L = If[s>L, L, s];
	f\[Beta]x = Re[LatticeFunction[lat, input, \[Beta]x, \[Delta]][L]];
	f\[Alpha]x = Re[LatticeFunction[lat, input, \[Alpha]x, \[Delta]][L]];
	f\[Gamma]x = Re[LatticeFunction[lat, input, \[Gamma]x, \[Delta]][L]];
	f\[Beta]y = Re[LatticeFunction[lat, input, \[Beta]y, \[Delta]][L]];
	f\[Alpha]y = Re[LatticeFunction[lat, input, \[Alpha]y, \[Delta]][L]];
	f\[Gamma]y = Re[LatticeFunction[lat, input, \[Gamma]y, \[Delta]][L]];
	opts = Sequence[FrameLabel->{"x / \!\(\*SqrtBox[\(\[Beta]\[Epsilon]\)]\)","x' / \!\(\*SqrtBox[\(\[Gamma]\[Epsilon]\)]\)"}, PlotLabel-> "\[Delta] = "<>ToString[NumberForm[\[Delta],2]], PlotLegends->LineLegend[{"x (beam ellipse)","y (beam ellipse)"},LabelStyle->{FontSize->15}], ImagePadding -> {{75, 0}, {60,5}}, BaseStyle->{FontSize->16},AxesOrigin->{0,0},ContourStyle->{Black,ColorData[97, 2]}];
	f[\[Delta]\[Delta]_]:=Module[{b,c,\[Epsilon]=1*^-6*0.511*^-3,it},
		it = TrackSingleParticle[lat,(1+\[Delta]\[Delta]),-Sqrt[TW\[Beta]x[input[]]*\[Epsilon]],Sqrt[TW\[Gamma]x[input[]]*\[Epsilon]],Sqrt[TW\[Beta]y[input[]]*\[Epsilon]],-Sqrt[TW\[Gamma]x[input[]]*\[Epsilon]],\[Delta]\[Delta]];
		b = ListPlot[{{it[[1]][[1]]/Sqrt[4*TW\[Beta]x[input[]]*\[Epsilon]],it[[1]][[2]]/Sqrt[4*TW\[Gamma]x[input[]]*\[Epsilon]]}}, PlotStyle->Black, PlotLegends->PointLegend[{"x (tracked particle)"},LabelStyle->{FontSize->15}]];
		c = ListPlot[{{it[[2]][[1]]/Sqrt[4*TW\[Beta]y[input[]]*\[Epsilon]],it[[2]][[2]]/Sqrt[4*TW\[Gamma]y[input[]]*\[Epsilon]]}}, PlotStyle->ColorData[97, 2], PlotLegends->PointLegend[{"y (tracked particle)"},LabelStyle->{FontSize->15}]]; 
		Return[{b,c}]
	];
	{bb,cc}=f[\[Delta]];
	a = ContourPlot[{
		(1+f\[Alpha]x^2)/f\[Beta]x*x^2+2*f\[Alpha]x*xp*x + f\[Beta]x*xp^2==1,
		(1+f\[Alpha]y^2)/f\[Beta]y*x^2+2*f\[Alpha]y*xp*x + f\[Beta]y*xp^2==1},
		{x,-size Sqrt[f\[Beta]x], size Sqrt[f\[Beta]x]},
		{xp,-size Sqrt[f\[Gamma]x], size Sqrt[f\[Gamma]x]},
		Evaluate[opts], PerformanceGoal->"Quality",PlotRange->{{-size,size},{-size,size}},AxesOrigin->{0,0}];
	Show[a,bb,cc]
];



(*** MATLAB INTERFACE for faster numerical solving ***)

(* NOTE: QuaILS imports "ToMatlab.m" for translating Mathematica to MATLAB *)
Import[FileNameJoin[{HOME,"lib/ToMatlab.m"}]];

(* Writes the constraints for the minimization routine to file *)
WriteMatlabConstraints[constraints_, variables_, fileID_:""] := Module[{str},
	Needs["MatlabUtils`ToMatlab`"];
	str = "% CONSTRAINTS AUTO-GENERATED BY MATHEMATICA \n" <> "function [c, ceq] = constraints_"<>fileID<>"(x)\n"
	   <> StringJoin[("\t" <> ToString[variables[[#]]] <> " = x(" <> ToString[#] <> ");\n")&/@Range[Length[variables]]] <> "\n"
	   <> StringJoin[("\t" <> "ceq" <> ToString[#] <> " = " <> ToString[constraints[[#]]//ToMatlab] <> "\n")& /@ Range[Length[constraints]]]
	   <> "\t" <> "c = [];\n\t" <> "ceq = [" <> StringJoin[("ceq" <> ToString[#] <> " ")& /@ Range[Length[constraints]]] <> "];\n" <> "end\n";
	Export[HOME<>"/matlab/constraints/constraints_"<>fileID<>".m",str,"Text"]
];

(* Writes the merit function for the minimization routine to file *)
WriteMatlabMeritFunction[meritfunction_, variables_, fileID_:""] := Module[{str},
	Needs["MatlabUtils`ToMatlab`"];
	str = "% MERIT FUNCTION AUTO-GENERATED BY MATHEMATICA \n" <> "function f = meritfunction_"<>fileID<>"(x)\n"
	   <> StringJoin[("\t" <> ToString[variables[[#]]] <> " = x(" <> ToString[#] <> ");\n")& /@ Range[Length[variables]]] <> "\n"
	   <> "\t" <> "f = " <> ToString[meritfunction//ToMatlab] <> "\n" <> "end\n";
	Export[HOME<>"/matlab/meritfunctions/meritfunction_"<>fileID<>".m",str,"Text"]
];

(* Writes the minimization routine parameters for the minimization routine to file *)
WriteMatlabMatchingParams[lowerBounds_, upperBounds_, maxStepShots_:100, maxShots_:500, funEvals_:2000, funTol_:1*^-10, xTol_:1*^-10] := Module[{str},
	Needs["MatlabUtils`ToMatlab`"];
	str = "% MATCHING PARAMETERS AUTO-GENERATED BY MATHEMATICA \n"
	   <> "function [lb, ub, maxStepShots, maxShots, funEvals, funTol, xTol] = matchingParams()\n"
	   <> "\t" <> "lb = " <> ToString[lowerBounds//ToMatlab]
	   <> "\t" <> "ub = " <> ToString[upperBounds//ToMatlab]
	   <> "\t" <> "maxStepShots = " <> ToString[maxStepShots//ToMatlab]
	   <> "\t" <> "maxShots = " <> ToString[maxShots//ToMatlab]
	   <> "\t" <> "funEvals = " <> ToString[funEvals//ToMatlab]	
	   <> "\t" <> "funTol = " <> ToString[funTol//ToMatlab]
	   <> "\t" <> "xTol = " <> ToString[xTol//ToMatlab] <> "end\n";
	Export[HOME<>"/matlab/matchingParams.m",str,"Text"]
];

(* Writes the initial guess for the minimization routine to file *)
WriteMatlabInitialGuess[guess_, fileID_:""] := Module[{str, guessStr},
	guessStr = If[Length[guess]==0, "[];\n", ToString[guess//ToMatlab]];
	str = "% MATCHING PARAMETERS AUTO-GENERATED BY MATHEMATICA \n"
	   <> "function [ x0 ] = initialGuess_"<>fileID<>"() \n"
       <> "    x0 = " <> guessStr
       <> "end";
	Export[HOME<>"/matlab/guesses/initialGuess_"<>fileID<>".m",str,"Text"]
];

(* Run the MATLAB minimization routine [see code in /matlab/quailsMinimizer.m] *)
RunMatlabMinimizer[vars_, guess_:{}, fileID_:""] := (
	Needs["MATLink`"];
	Quiet[MATLink`OpenMATLAB[]];
	MATLink`MEvaluate["cd '"<>HOME<>"/matlab'"];
	WriteMatlabInitialGuess[guess, fileID];
	Print[MATLink`MEvaluate["quailsMinimizer('"<>fileID<>"')"]];
	Return[ReadMatlabSolutions[vars, fileID]];
);

(* Extract the solutions that MATLAB found: run after RunMatlabMinimizer[] *)
ReadMatlabSolutions[variables_, fileID_:""] := Module[{nums},
	nums = Import[HOME<>"/matlab/solutions/solutions_"<>fileID<>".dat","List"];
	If[Length[variables]==Length[nums], (variables[[#]]->nums[[#]])&/@Range[Length[nums]], {}]
];



(*** ELEGANT INTERFACE FOR PARTICLE TRACKING ***)

(* NOTE: QuaILS requires an ELEGANT binary in the "elegant/bin"-directory *)

(* Tracks a beam of N particles using ELEGANT *)
TrackElegantLattice[lat_, E_, N_, \[Sigma]E_, \[Epsilon]x_, \[Epsilon]y_, input_, \[Delta]_:0] := Module[{track},
	Needs["MATLink`"];
	Quiet[MATLink`OpenMATLAB[]];
	WriteElegantLattice[lat];
	WriteElegantParams[E, N, \[Sigma]E, \[Epsilon]x, \[Epsilon]y, input, \[Delta]];
	MATLink`MEvaluate["cd '"<>HOME<>"/elegant'"];
	track = MATLink`MFunction["trackLattice"];
	Return[track[HOME]];
];

(* Tracks a single particle using ELEGANT *)
TrackSingleParticle[lat_, E_, x_, xp_, y_, yp_, \[Delta]_:0] := Module[{track,x1,xp1,y1,yp1,x0,xp0,y0,yp0},
	Needs["MATLink`"];
	Quiet[MATLink`OpenMATLAB[]];
	WriteElegantLattice[lat];
	WriteElegantParams[E, 1, 0, 0, 0, TwissParams[1,1], \[Delta], x, xp, y, yp];
	MATLink`MEvaluate["cd '"<>HOME<>"/elegant'"];
	track = MATLink`MFunction["trackParticle"];
	{x1,xp1,y1,yp1,x0,xp0,y0,yp0} = track[HOME];
	Return[{{x1,-xp1},{y1,-yp1}}];
];

(* INTERNAL - Writes ELEGANT lattice [directly from the QuaILS lattice] to file *)
WriteElegantLattice[lat_] := Module[{str, dlen, qLen = 1*^-6},
	
	str = "! QUAILS AUTO-GENERATED TRACKING LATTICE \n";
	
	str = str <> "\n! DRIFT SPACES \n";
	Do[(
	    If[IsThick[lat],
			dlen = DS[lat][[i]],
			dlen = If[Or[i==1,i==Length[DS[lat]]],(DS[lat][[i]]-qLen/2),(DS[lat][[i]]-qLen)]];
		str = str <> "D" <> ToString[i] <> " : DRIF, L = " <> ToString[dlen//CForm] <> "\n";);
	,{i,1,Length[DS[lat]]}];
	
	str = str <> "\n! QUADRUPOLES \n";
	
	Do[
		str = If[IsThick[lat],
			str <> "Q" <> ToString[i] <> " : KQUAD, L = " <> ToString[LS[lat][[i]]//CForm] <> ", "
		        <> "INTEGRATION_ORDER = 4, ISR = 0, SYNCH_RAD = 0, N_KICKS = 100, "
		        <> "K1 = " <> ToString[KS[lat][[i]]//CForm] <> "\n",
			str <> "Q" <> ToString[i] <> " : KQUAD, L = " <> ToString[qLen//CForm] <> ", "
		        <> "INTEGRATION_ORDER = 4, ISR = 0, SYNCH_RAD = 0, N_KICKS = 1, "
		        <> "K1 = " <> ToString[(KS[lat][[i]]/qLen)//CForm] <> "\n"]
	,{i,1,Length[KS[lat]]}];
	
	str = str <> "\n! LINE \n"
	          <> "LATTICE : LINE = (";
	Do[
		str = str <> "D" <> ToString[i] <> ", Q" <> ToString[i] <> ", ";
	,{i,1,Length[KS[lat]]}];
	str = str <> "D" <> ToString[Length[DS[lat]]] <> ")\n\n";

	Export[HOME<>"/elegant/lattice.lte", str, "Text"]
];

(* INTERNAL - Writes ELEGANT parameters to file *)
WriteElegantParams[E_, N_, \[Sigma]E_, \[Epsilon]x_, \[Epsilon]y_, input_, \[Delta]_:0, x_:0, xp_:0, y_:0, yp_:0] := Module[{str},
	str = "! QUAILS AUTO-GENERATED TRACKING PARAMETERS \n\n";
    
	str = str <> "&run_setup \n"
              <> "    lattice       = \"lattice.lte\",\n"
              <> "    use_beamline  = LATTICE,\n"
              <> "    output        = \"final.bun\",\n"
              <> "    p_central_mev = " <> ToString[E*1000//CForm] <> ",\n"
              <> "&end \n\n";

    str = str <> "&run_control \n"
              <> "&end \n\n";
              
    str = str <> "&bunched_beam \n"
              <> "    bunch                  = \"initial.bun\",\n"
              <> "    n_particles_per_bunch  = " <> ToString[N] <> ",\n"
              <> "    distribution_cutoff[0] = 5,5,5,\n"
              <> "    Po                     = " <> ToString[((E/(1+\[Delta]))/0.000510998903076601)//CForm] <> ",\n"
              <> "    sigma_dp               = " <> ToString[\[Sigma]E//CForm] <> ",\n"
              <> "    emit_nx                = " <> ToString[\[Epsilon]x//CForm] <> ",\n"
              <> "    emit_ny                = " <> ToString[\[Epsilon]y//CForm] <> ",\n"
              <> "    beta_x                 = " <> ToString[TW\[Beta]x[input]//CForm] <> ",\n"
              <> "    beta_y                 = " <> ToString[TW\[Beta]y[input]//CForm] <> ",\n"
              <> "    alpha_x                = " <> ToString[TW\[Alpha]x[input]//CForm] <> ",\n"
              <> "    alpha_y                = " <> ToString[TW\[Alpha]y[input]//CForm] <> ",\n"
              <> "    centroid[0]            = " <> ToString[x//CForm] <> "," 
                                                 <> ToString[xp//CForm] <> "," 
                                                 <> ToString[y//CForm] <> "," 
                                                 <> ToString[yp//CForm] <> ",0,0,\n"
              <> "&end \n\n";
    
    str = str <> "&track \n"
              <> "&end \n\n";
             
	Export[HOME<>"/elegant/params.ele", str, "Text"];
];




End[];
EndPackage[];
