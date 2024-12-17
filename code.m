(* ::Package:: *)

(* ::Subsection:: *)
(*Info*)


CellPrint[TextCell["A Mathematica code to generate Hilbert Series.
Authors:  Rodrigo Alonso & Shakeel Ur Rahaman  
IPPP, Durham University                 
arxiv [2412.09463]",CellFrame->True(*,CellMargins->{{10, 800}, {10, 10}}*),Background->Purple,FontSize->20,FontColor->White]];


(* ::Input::Initialization:: *)
Quiet[Remove["Global`*"]]
$Assumptions=\[DoubleStruckCapitalD]>0 && q >0;


CellPrint[TextCell[TextGrid[{{"HEFT:"," ","   ","SMEFT:"," "},{"h","Physical Higgs.","   ","H(d)","Higgs Doublet (its conjugate)."},
{"Vp","(+) goldstone","   ","LL(d)","Left handed lepton doublet (its conjugate)."},{"Vm","(-) goldstone.","   ","QL(d):","Left handed quark doublet (its conjugate)."},
{"Vz","Neutral goldstone boson.","   ","",""},
{"eL(d)","Left handed electron (its conjugate).","   ","",""},
{"\[Nu]L(d)","Left handed neutrino (its conjugate).","   ","",""},
{"eR(d)","Right handed electron (its conjugate).","   ","eR(d)","Right handed electron (its conjugate)."},
{"uL(d)","Left handed up quark (its conjugate).","   ","",""},
{"uR(d)","Right handed up quark (its conjugate).","   ","uR(d)","Right handed up quark (its conjugate)."},
{"dL(d)","Left handed down quark (its conjugate).","   ","",""},
{"dR(d)","Right handed down quark (its conjugate).","   ","dR(d)","Right handed down quark (its conjugate)."},
{"GL","Left handed SU(3)field strength tensor.","   ","GL","Left handed SU(3)field strength tensor."},
{"GR","Right handed SU(3)field strength tensor.","   ","GR","Right handed SU(3)field strength tensor."},
{"WpL","Left handed (+) W boson field strength tensor.","   ","WL","Left handed W boson field strength tensor."},
{"WpR","Right handed (+) W boson field strength tensor.","   ","WR","Right handed W boson field strength tensor."},
{"WmL","Left handed (-) W boson field strength tensor.","   ","BL","Left handed B field strength tensor."},
{"WmR","Right handed (-) W boson field strength tensor.","   ","BR","Right handed B field strength tensor."},
{"ZL","Left handed Z boson field strength tensor.","   ","",""},
{"ZR","Right handed Z boson field strength tensor.","   ",""," *** Generic Field names to define operator classes ***"},
{"AL","Left handed photon field strength tensor.","   ","S","Scalar field"},
{"AR","Right handed photon field strength tensor.","   ","F","Fermion Field"},
{"U(d)","Spurion (its conjugate).","   ","X","Field Strenght tensor"},
{"V","Goldstone in Linear representation.","   ","\[DoubleStruckCapitalD]","Covariant Derivative"}},Alignment->{Center,Baseline},FrameStyle->Black,Background->LightBrown,Frame->True](*,CellMargins->{{10, 800}, {10, 10}}*),FontSize->10]]


(* ::Subsection::Closed:: *)
(*Haar Measures*)


(* ::Input::Initialization:: *)
HaarSO4[\[Alpha]_,\[Beta]_]:=1/(4\[Alpha] \[Beta]) (1-\[Alpha]^2)(1-\[Alpha]^-2)(1-\[Beta]^2)(1-\[Beta]^-2);
HaarU1[y_]:=1/y;
HaarSU2[y_]:=1/(2y) (1-y^2)(1-y^-2);
HaarSU3[y1_,y2_]:=1/(6y1 y2) (1-y1 y2)(1-y1^2/y2)(1-y2^2/y1)(1-1/(y1 y2))(1-y1/y2^2)(1-y2/y1^2);


(* ::Subsection::Closed:: *)
(*Momentum Generating function*)


(* ::Input::Initialization:: *)
P[q_,\[Alpha]_,\[Beta]_]:=1/((1-q \[Alpha] \[Beta])(1-q \[Alpha]^-1 \[Beta]^-1)(1-q \[Alpha] \[Beta]^-1)(1-q \[Alpha]^-1 \[Beta]));


(* ::Subsection::Closed:: *)
(*Gauge Characters*)


(* ::Input::Initialization:: *)
CharU1[x_,y_]:=y^(6*x);
CharSU2fund[z_]:=z+1/z;
CharSU2adj[z_]:=CharSU2fund[z]*CharSU2fund[z]-1;
CharSU3fund[y1_,y2_]:=y1+1/y2+y2/y1;
CharSU3antifund[y1_,y2_]:=y2+1/y1+y1/y2;
CharSU3adj[y1_,y2_]:=CharSU3fund[y1,y2]*CharSU3antifund[y1,y2]-1;
CharSO4Scalar[\[Alpha]_,\[Beta]_]:= 1;
CharSO4LFerm[\[Alpha]_,\[Beta]_]:=\[Alpha]+\[Alpha]^-1;
CharSO4RFerm[\[Alpha]_,\[Beta]_]:=\[Beta]+\[Beta]^-1;
CharSO4LFS[\[Alpha]_,\[Beta]_]:=\[Alpha]^2+1+\[Alpha]^-2;
CharSO4RFS[\[Alpha]_,\[Beta]_]:=\[Beta]^2+1+\[Beta]^-2;
CharSO4V[\[Alpha]_,\[Beta]_]:= \[Alpha] \[Beta]+\[Alpha]^-1 \[Beta]^-1+\[Alpha] \[Beta]^-1+\[Alpha]^-1 \[Beta];
CharSO4F[\[Alpha]_,\[Beta]_]:=\[Alpha]^2+1+\[Alpha]^-2+\[Beta]^2+1+\[Beta]^-2;


(* ::Subsection::Closed:: *)
(*Conformal Characters*)


(* ::Input::Initialization:: *)
CharCScalar[q_,\[Alpha]_,\[Beta]_]:=q(1-q^2)P[q,\[Alpha],\[Beta]];
CharCGoldstone[q_,\[Alpha]_,\[Beta]_]:=(1-q^2)P[q,\[Alpha],\[Beta]]-1;
CharCLFerm[q_,\[Alpha]_,\[Beta]_]:=q^(3/2) (CharSO4LFerm[\[Alpha],\[Beta]]-q CharSO4RFerm[\[Alpha],\[Beta]])P[q,\[Alpha],\[Beta]];
CharCRFerm[q_,\[Alpha]_,\[Beta]_]:=q^(3/2) (CharSO4RFerm[\[Alpha],\[Beta]]-q CharSO4LFerm[\[Alpha],\[Beta]])P[q,\[Alpha],\[Beta]];
CharCLFS[q_,\[Alpha]_,\[Beta]_]:=q^2 (CharSO4LFS[\[Alpha],\[Beta]]-q CharSO4V[\[Alpha],\[Beta]]+q^2)P[q,\[Alpha],\[Beta]];
CharCRFS[q_,\[Alpha]_,\[Beta]_]:=q^2 (CharSO4RFS[\[Alpha],\[Beta]]-q CharSO4V[\[Alpha],\[Beta]]+q^2)P[q,\[Alpha],\[Beta]];


(* ::Subsection::Closed:: *)
(*Goldstones*)


(* ::Text:: *)
(*Singlets*)


(* ::Input::Initialization:: *)
CharVp[Vp_,q_,\[Alpha]_,\[Beta]_,w_]:=Vp CharCGoldstone[q,\[Alpha],\[Beta]]CharU1[1,w];
CharVm[Vm_,q_,\[Alpha]_,\[Beta]_,w_]:=Vm CharCGoldstone[q,\[Alpha],\[Beta]]CharU1[-1,w];
CharVz[Vz_,q_,\[Alpha]_,\[Beta]_]:=Vz CharCGoldstone[q,\[Alpha],\[Beta]];


(* ::Text:: *)
(*Triplet*)


(* ::Input::Initialization:: *)
CharV[V_,q_,\[Alpha]_,\[Beta]_,s_]:=V CharCGoldstone[q,\[Alpha],\[Beta]]CharSU2adj[s];


(* ::Subsection::Closed:: *)
(*Scalar*)


(* ::Text:: *)
(*Physical Higgs*)


(* ::Input::Initialization:: *)
Charh[h_,q_,\[Alpha]_,\[Beta]_]:=h CharCScalar[q,\[Alpha],\[Beta]];


(* ::Text:: *)
(*Higgs doublet*)


(* ::Input::Initialization:: *)
CharH[H_,q_,\[Alpha]_,\[Beta]_,s_,w_]:=H CharCScalar[q,\[Alpha],\[Beta]] CharSU2fund[s] CharU1[1/2,w];
CharHd[Hd_,q_,\[Alpha]_,\[Beta]_,s_,w_]:=Hd CharCScalar[q,\[Alpha],\[Beta]] CharSU2fund[s] CharU1[-1/2,w];


(* ::Text:: *)
(*Physical Higgs as Goldstone*)


(* ::Input::Initialization:: *)
CharhGB[GB_,q_,\[Alpha]_,\[Beta]_]:=GB CharCGoldstone[q,\[Alpha],\[Beta]];


(* ::Subsection:: *)
(*Spurions*)


(* ::Text:: *)
(*Triplet*)


(* ::Input::Initialization:: *)
CharS[S_,s_]:= S CharSU2adj[s];


(* ::Text:: *)
(*Doublet with U(1)*)


(* ::Input::Initialization:: *)
CharT[T_,s_,w_]:= T CharSU2fund[s] CharU1[1/2,w];
CharTd[Td_,s_,w_]:= Td CharSU2fund[s]CharU1[-1/2,w];


(* ::Subsection::Closed:: *)
(*Fermions*)


(* ::Text:: *)
(*Quarks*)


(* ::Input::Initialization:: *)
CharuL[uL_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nu uL CharCLFerm[q,\[Alpha],\[Beta]]CharSU3fund[z1,z2]CharU1[2/3,w];
CharuLd[uLd_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nu uLd CharCRFerm[q,\[Alpha],\[Beta]]CharSU3antifund[z1,z2]CharU1[-2/3,w];

ChardL[dL_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nd dL CharCLFerm[q,\[Alpha],\[Beta]]CharSU3fund[z1,z2]CharU1[-1/3,w];
ChardLd[dLd_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nd dLd CharCRFerm[q,\[Alpha],\[Beta]]CharSU3antifund[z1,z2]CharU1[1/3,w];

CharuR[uR_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nu uR CharCRFerm[q,\[Alpha],\[Beta]]CharSU3fund[z1,z2]CharU1[2/3,w];
CharuRd[uRd_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nu uRd CharCLFerm[q,\[Alpha],\[Beta]]CharSU3antifund[z1,z2]CharU1[-2/3,w];

ChardR[dR_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nd dR CharCRFerm[q,\[Alpha],\[Beta]]CharSU3fund[z1,z2]CharU1[-1/3,w];
ChardRd[dRd_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,w_]:=nd dRd CharCLFerm[q,\[Alpha],\[Beta]]CharSU3antifund[z1,z2]CharU1[1/3,w];


(* ::Text:: *)
(*Leptons*)


(* ::Input::Initialization:: *)
Char\[Nu]L[\[Nu]L_,q_,\[Alpha]_,\[Beta]_]:=n\[Nu] \[Nu]L CharCLFerm[q,\[Alpha],\[Beta]];
Char\[Nu]Ld[\[Nu]Ld_,q_,\[Alpha]_,\[Beta]_]:=n\[Nu] \[Nu]Ld CharCRFerm[q,\[Alpha],\[Beta]];

ChareL[eL_,q_,\[Alpha]_,\[Beta]_,w_]:=ne eL CharCLFerm[q,\[Alpha],\[Beta]]CharU1[-1,w];
ChareLd[eLd_,q_,\[Alpha]_,\[Beta]_,w_]:=ne eLd CharCRFerm[q,\[Alpha],\[Beta]]CharU1[1,w];

ChareR[eR_,q_,\[Alpha]_,\[Beta]_,w_]:=ne eR CharCRFerm[q,\[Alpha],\[Beta]]CharU1[-1,w];
ChareRd[eRd_,q_,\[Alpha]_,\[Beta]_,w_]:=ne eRd CharCLFerm[q,\[Alpha],\[Beta]]CharU1[1,w];

Char\[Nu]R[\[Nu]R_,q_,\[Alpha]_,\[Beta]_]:=n\[Nu] \[Nu]R CharCRFerm[q,\[Alpha],\[Beta]];
Char\[Nu]Rd[\[Nu]Rd_,q_,\[Alpha]_,\[Beta]_]:=n\[Nu] \[Nu]Rd CharCLFerm[q,\[Alpha],\[Beta]];


(* ::Text:: *)
(*Quark Doublets*)


(* ::Input::Initialization:: *)
CharQL[QL_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,s_,w_]:=nQ QL CharCLFerm[q,\[Alpha],\[Beta]]CharSU3fund[z1,z2]CharSU2fund[s]CharU1[1/6,w];
CharQLd[QLd_,q_,\[Alpha]_,\[Beta]_,z1_,z2_,s_,w_]:=nQ QLd CharCRFerm[q,\[Alpha],\[Beta]]CharSU3antifund[z1,z2]CharSU2fund[s]CharU1[-1/6,w];


(* ::Text:: *)
(*Lepton Doublet*)


(* ::Input::Initialization:: *)
CharLL[LL_,q_,\[Alpha]_,\[Beta]_,s_,w_]:=nL LL CharCLFerm[q,\[Alpha],\[Beta]]CharSU2fund[s] CharU1[-1/2,w];
CharLLd[LLd_,q_,\[Alpha]_,\[Beta]_,s_,w_]:=nL LLd CharCRFerm[q,\[Alpha],\[Beta]]  CharSU2fund[s]CharU1[1/2,w];


(* ::Subsection::Closed:: *)
(*Field Strength Tensors*)


(* ::Input::Initialization:: *)
CharGL[GL_,q_,\[Alpha]_,\[Beta]_,z1_,z2_]:=GL CharCLFS[q,\[Alpha],\[Beta]]CharSU3adj[z1,z2];
CharGR[GR_,q_,\[Alpha]_,\[Beta]_,z1_,z2_]:=GR CharCRFS[q,\[Alpha],\[Beta]]CharSU3adj[z1,z2];

CharWpL[WpL_,q_,\[Alpha]_,\[Beta]_,w_]:=WpL CharCLFS[q,\[Alpha],\[Beta]]CharU1[1,w];
CharWpR[WpR_,q_,\[Alpha]_,\[Beta]_,w_]:=WpR CharCRFS[q,\[Alpha],\[Beta]]CharU1[1,w];

CharWmL[WmL_,q_,\[Alpha]_,\[Beta]_,w_]:=WmL CharCLFS[q,\[Alpha],\[Beta]]CharU1[-1,w];
CharWmR[WmR_,q_,\[Alpha]_,\[Beta]_,w_]:=WmR CharCRFS[q,\[Alpha],\[Beta]]CharU1[-1,w];

CharZL[ZL_,q_,\[Alpha]_,\[Beta]_]:=ZL CharCLFS[q,\[Alpha],\[Beta]];
CharZR[ZR_,q_,\[Alpha]_,\[Beta]_]:=ZR CharCRFS[q,\[Alpha],\[Beta]];

CharAL[AL_,q_,\[Alpha]_,\[Beta]_]:=AL CharCLFS[q,\[Alpha],\[Beta]];
CharAR[AR_,q_,\[Alpha]_,\[Beta]_]:=AR CharCRFS[q,\[Alpha],\[Beta]];

CharWL[WL_,q_,\[Alpha]_,\[Beta]_,s_]:=WL CharCLFS[q,\[Alpha],\[Beta]]CharSU2adj[s];
CharWR[WR_,q_,\[Alpha]_,\[Beta]_,s_]:=WR CharCRFS[q,\[Alpha],\[Beta]]CharSU2adj[s];

CharBL[BL_,q_,\[Alpha]_,\[Beta]_]:=BL CharCLFS[q,\[Alpha],\[Beta]];
CharBR[BR_,q_,\[Alpha]_,\[Beta]_]:=BR CharCRFS[q,\[Alpha],\[Beta]];


(* ::Subsection:: *)
(*Plethystic Arguments*)


(* ::Subsubsection:: *)
(*Fermionic Plethystic Argument*)


(* ::Text:: *)
(*Fermion Left Singlets*)


(* ::Input::Initialization:: *)
argPEFermionSingletLeft[\[CapitalDelta]_]:=Refine[Sum[(-1)^(r+1)/r CharuL[(uL*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]+(-1)^(r+1)/r CharuLd[(uLd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]+
(-1)^(r+1)/r ChardL[(dL*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]+(-1)^(r+1)/r ChardLd[(dLd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]+
(-1)^(r+1)/r Char\[Nu]L[(\[Nu]L*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]+(-1)^(r+1)/r Char\[Nu]Ld[(\[Nu]Ld*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r] +
(-1)^(r+1)/r ChareL[(eL*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+(-1)^(r+1)/r ChareLd[(eLd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]
,{r,1,Floor[2\[CapitalDelta]/3]}],Assumptions->q>0];


(* ::Text:: *)
(*Fermion Right Singlets*)


(* ::Input::Initialization:: *)
argPEFermionSingletRight[\[CapitalDelta]_]:=Refine[Sum[(-1)^(r+1)/r CharuR[(uR*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]+(-1)^(r+1)/r CharuRd[(uRd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]
+(-1)^(r+1)/r ChardR[(dR*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]+(-1)^(r+1)/r ChardRd[(dRd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,w^r]
+(-1)^(r+1)/r ChareR[(eR*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+(-1)^(r+1)/r ChareRd[(eRd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]
,{r,1,Floor[2\[CapitalDelta]/3]}],Assumptions->q>0];


(* ::Input::Initialization:: *)
argPERHN[\[CapitalDelta]_]:=Refine[Sum[
(-1)^(r+1)/r Char\[Nu]R[(\[Nu]R*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]+(-1)^(r+1)/r Char\[Nu]Rd[(\[Nu]Rd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]
,{r,1,Floor[2\[CapitalDelta]/3]}],Assumptions->q>0];


(* ::Text:: *)
(*Fermion Left Doublets*)


(* ::Input::Initialization:: *)
argPEFermionDoublets[\[CapitalDelta]_]:=Refine[Sum[(-1)^(r+1)/r CharQL[(QL*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,s^r,w^r]+(-1)^(r+1)/r CharQLd[(QLd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r,s^r,w^r]+(-1)^(r+1)/r CharLL[(LL*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,s^r,w^r]+(-1)^(r+1)/r CharLLd[(LLd*\[DoubleStruckCapitalD]^(-3/2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,s^r,w^r],{r,1,Floor[2\[CapitalDelta]/3]}],Assumptions->q>0];


(* ::Subsubsection::Closed:: *)
(*Scalar Plethystic Argument*)


(* ::Text:: *)
(*Scalar Singlet*)


(* ::Input::Initialization:: *)
argPEHiggsPhysical[\[CapitalDelta]_]:=Sum[1/r Charh[(h/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r],{r,1,\[CapitalDelta]}];


(* ::Input::Initialization:: *)
argPEHiggsPhysicalGB[\[CapitalDelta]_]:=Sum[1/r CharhGB[(GB/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r],{r,1,\[CapitalDelta]}];


(* ::Text:: *)
(*Scalar Doublet*)


(* ::Input::Initialization:: *)
argPEHiggsDoublet[\[CapitalDelta]_]:=Sum[(1/r) CharH[(H/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,s^r,w^r]+(1/r) CharHd[(Hd/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,s^r,w^r],{r,1,\[CapitalDelta]}];


(* ::Subsubsection::Closed:: *)
(*Goldstone Plethystic Argument*)


(* ::Text:: *)
(*Goldstones Singlets*)


(* ::Input::Initialization:: *)
argPEGoldstoneSinglet[\[CapitalDelta]_]:=Sum[1/r CharVp[(Vp/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+1/r CharVm[(Vm/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+1/r CharVz[(Vz/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r],{r,1,\[CapitalDelta]}];


(* ::Text:: *)
(*Goldstone Triplet*)


(* ::Input::Initialization:: *)
argPEGoldstoneTriplet[\[CapitalDelta]_]:=Sum[1/r CharV[(V/\[DoubleStruckCapitalD])^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,s^r],{r,1,\[CapitalDelta]}];


(* ::Subsubsection:: *)
(*Spurion Plethystic Argument*)


(* ::Text:: *)
(*Spurion*)


(* ::Input::Initialization:: *)
argPESpurionTriplet[C_]:=Sum[1/r CharS[S^r,s^r],{r,1,C}];
argPESpurion[C_]:=Sum[1/r CharT[U^r,s^r,w^r]+1/r CharTd[Ud^r,s^r,w^r],{r,1,C}];


(* ::Text:: *)
(*Spurion Generating Function*)


(* ::Input::Initialization:: *)
SpurionGenFunc2[dim_]:=Normal[Series[Exp[argPESpurionTriplet[dim+1]],{S,0,dim}]];
SpurionGenFunc[dim_]:=Normal[Series[Normal[Series[Exp[(argPESpurion[dim])],{U,0,dim}]],{Ud,0,dim}]];
(*SpurionGenFunc[dim_]:=Normal[Series[Normal[Series[1/(((1-T w)(1-T s w)(1 - T w /s))((1-Td/ w)(1-Td s/ w)(1 - Td /(w s)))),{T,0,dim}]],{Td,0,dim}]];*)


(* ::Subsubsection::Closed:: *)
(*FS Plethystic Argument*)


(* ::Text:: *)
(*Field Strength*)


(* ::Input::Initialization:: *)
argPEFSHEFT[\[CapitalDelta]_]:=Sum[1/r CharGL[(GL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r]+1/r CharGR[(GR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r]+1/r CharWpL[(WpL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+1/r CharWpR[(WpR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+
1/r CharWmL[(WmL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+1/r CharWmR[(WmR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,w^r]+
	1/r CharZL[(ZL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]+1/r CharZR[(ZR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]+
	1/r CharAL[(AL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]+1/r CharAR[(AR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]
,{r,1,Floor[\[CapitalDelta]/2]}];


(* ::Input::Initialization:: *)
argPEFSSMEFT[\[CapitalDelta]_]:=Sum[1/r CharGL[(GL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r]+1/r CharGR[(GR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,z1^r,z2^r]+
1/r CharWL[(WL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,s^r]+1/r CharWR[(WR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r,s^r]+
	1/r CharBL[(BL*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]+1/r CharBR[(BR*\[DoubleStruckCapitalD]^(-2))^r,(q*\[DoubleStruckCapitalD])^r,\[Alpha]^r,\[Beta]^r]
,{r,1,Floor[\[CapitalDelta]/2]}];


(* ::Subsubsection:: *)
(*Sum of Arguments*)


(* ::Input::Initialization:: *)
argHEFT[\[CapitalDelta]_]:=argPEFermionSingletLeft[\[CapitalDelta]]+argPEFermionSingletRight[\[CapitalDelta]]+argPEHiggsPhysical[\[CapitalDelta]]+argPEGoldstoneSinglet[\[CapitalDelta]]+argPEFSHEFT[\[CapitalDelta]];


(* ::Input::Initialization:: *)
(*argHEFTXcheck[\[CapitalDelta]_]:=argPEFermionSingletLeft[\[CapitalDelta]]+argPEFermionSingletRight[\[CapitalDelta]]+argPEHiggsPhysical[\[CapitalDelta]]+argPEGoldstoneSinglet[\[CapitalDelta]]+argPEFSHEFT[\[CapitalDelta]]+argPERHN[\[CapitalDelta]];*)


(* ::Input::Initialization:: *)
argSMEFT[\[CapitalDelta]_]:=argPEFermionDoublets[\[CapitalDelta]]+argPEFermionSingletRight[\[CapitalDelta]]+argPEHiggsDoublet[\[CapitalDelta]]+argPEFSSMEFT[\[CapitalDelta]](*+argPERHN[\[CapitalDelta]]*);


(* ::Input::Initialization:: *)
argHEFTSpurion[\[CapitalDelta]_]:=argPEFermionDoublets[\[CapitalDelta]]+argPEFermionSingletRight[\[CapitalDelta]]+argPEHiggsPhysical[\[CapitalDelta]]+argPEGoldstoneTriplet[\[CapitalDelta]]+argPEFSSMEFT[\[CapitalDelta]];


(* ::Input::Initialization:: *)
argHEFTMod[\[CapitalDelta]_]:=argPEFermionSingletLeft[\[CapitalDelta]]+argPEFermionSingletRight[\[CapitalDelta]]+argPEHiggsPhysicalGB[\[CapitalDelta]]+argPEGoldstoneSinglet[\[CapitalDelta]]+argPEFSHEFT[\[CapitalDelta]](*+argPERHN[\[CapitalDelta]]*);


(* ::Subsection:: *)
(*Integrand*)


(* ::Input::Initialization:: *)
\[CapitalDelta]H[1]=0;
\[CapitalDelta]H[2]=\[DoubleStruckCapitalD] Vz;
\[CapitalDelta]H[3]=\[DoubleStruckCapitalD]^2 h;
\[CapitalDelta]H[4] =-AL \[DoubleStruckCapitalD]^2-AR \[DoubleStruckCapitalD]^2-\[DoubleStruckCapitalD]^4+\[DoubleStruckCapitalD] dL dLd nd^2+\[DoubleStruckCapitalD] dR dRd nd^2+\[DoubleStruckCapitalD] eL eLd ne^2+\[DoubleStruckCapitalD] eR eRd ne^2+\[DoubleStruckCapitalD] nu^2 uL uLd+\[DoubleStruckCapitalD] nu^2 uR uRd-\[DoubleStruckCapitalD]^2 Vm Vp+AL \[DoubleStruckCapitalD] Vz+AR \[DoubleStruckCapitalD] Vz+\[DoubleStruckCapitalD]^3 Vz+\[DoubleStruckCapitalD] Vm Vp Vz+\[DoubleStruckCapitalD] Vp WmL+\[DoubleStruckCapitalD] Vp WmR+\[DoubleStruckCapitalD] Vm WpL+\[DoubleStruckCapitalD] Vm WpR-\[DoubleStruckCapitalD]^2 ZL+\[DoubleStruckCapitalD] Vz ZL-\[DoubleStruckCapitalD]^2 ZR+\[DoubleStruckCapitalD] Vz ZR+\[DoubleStruckCapitalD] n\[Nu]^2 \[Nu]L \[Nu]Ld;


(* ::Input::Initialization:: *)
KineticTermsHEFT:=dL dLd nd^2 \[DoubleStruckCapitalD]+dR dRd nd^2 \[DoubleStruckCapitalD]+eL eLd ne^2 \[DoubleStruckCapitalD]+eR eRd ne^2 \[DoubleStruckCapitalD]+nu^2 uL uLd \[DoubleStruckCapitalD]+nu^2 uR uRd \[DoubleStruckCapitalD]+n\[Nu]^2 \[DoubleStruckCapitalD] \[Nu]L \[Nu]Ld + \[DoubleStruckCapitalD]^2 h^2;
\[CapitalDelta]HSMEFT:=2dR dRd nd^2 \[DoubleStruckCapitalD]+2eR eRd ne^2 \[DoubleStruckCapitalD]+2LL LLd nL^2 \[DoubleStruckCapitalD]+2nQ^2 QL QLd \[DoubleStruckCapitalD]+2nu^2 uR uRd \[DoubleStruckCapitalD]-BL \[DoubleStruckCapitalD]^2-BR \[DoubleStruckCapitalD]^2+2H Hd \[DoubleStruckCapitalD]^2-\[DoubleStruckCapitalD]^4


(* ::Input::Initialization:: *)
\[CapitalDelta]HSpurion[1]:=0;
\[CapitalDelta]HSpurion[2]:=U Ud V \[DoubleStruckCapitalD];
\[CapitalDelta]HSpurion[3]:=h \[DoubleStruckCapitalD]^2;
\[CapitalDelta]HSpurion[4]:=dR dRd nd^2 \[DoubleStruckCapitalD]+eR eRd ne^2 \[DoubleStruckCapitalD]+LL LLd nL^2 \[DoubleStruckCapitalD]+nQ^2 QL QLd \[DoubleStruckCapitalD]+LL LLd nL^2 U Ud \[DoubleStruckCapitalD]+nQ^2 QL QLd U Ud \[DoubleStruckCapitalD]+nu^2 uR uRd \[DoubleStruckCapitalD]+BL U Ud V \[DoubleStruckCapitalD]+BR U Ud V \[DoubleStruckCapitalD]+V^3 \[DoubleStruckCapitalD]+V WL \[DoubleStruckCapitalD]+U Ud V WL \[DoubleStruckCapitalD]+U^2 Ud^2 V WL \[DoubleStruckCapitalD]+V WR \[DoubleStruckCapitalD]+U Ud V WR \[DoubleStruckCapitalD]+U^2 Ud^2 V WR \[DoubleStruckCapitalD]+U Ud V \[DoubleStruckCapitalD]^3(**)-BL \[DoubleStruckCapitalD]^2-BR \[DoubleStruckCapitalD]^2-U Ud V^2 \[DoubleStruckCapitalD]^2-U Ud WL \[DoubleStruckCapitalD]^2-U Ud WR \[DoubleStruckCapitalD]^2-\[DoubleStruckCapitalD]^4;
KineticTermsHEFTSpurion:= dR dRd nd^2 \[DoubleStruckCapitalD]+ eR eRd ne^2 \[DoubleStruckCapitalD]+ LL LLd nL^2 \[DoubleStruckCapitalD]+nQ^2 QL QLd \[DoubleStruckCapitalD]+ nu^2 uR uRd \[DoubleStruckCapitalD]+ \[DoubleStruckCapitalD]^2 h^2;


(* ::Input::Initialization:: *)
arg\[CapitalDelta]H:=Expand[Sum[\[CapitalDelta]H[i]*q^i,{i,1,4}]];


(* ::Input::Initialization:: *)
arg\[CapitalDelta]HMod:=Expand[((1-q h)Sum[\[CapitalDelta]H[i]*q^i,{i,1,4}])/.h->GB/(q \[DoubleStruckCapitalD])];


(* ::Input::Initialization:: *)
HEFTIntegrand[\[CapitalDelta]_]:=Module[{Integrand},
Integrand=SeriesCoefficient[1/P[q*\[DoubleStruckCapitalD],\[Alpha],\[Beta]] Exp[argHEFT[\[CapitalDelta]]],{q,0,\[CapitalDelta]}];
Integrand//Expand
];


(* ::Input::Initialization:: *)
HEFTformfactorIntegrand[\[CapitalDelta]_]:=Module[{Integrand},
Integrand=SeriesCoefficient[(1-q*h)(1/P[q*\[DoubleStruckCapitalD],\[Alpha],\[Beta]] Exp[argHEFT[\[CapitalDelta]]]+arg\[CapitalDelta]H),{q,0,\[CapitalDelta]}];
Integrand//Expand
];


(* ::Input::Initialization:: *)
HEFTIntegrandModded[\[CapitalDelta]_]:=Module[{Integrand,expr},
Integrand=SeriesCoefficient[1/P[q*\[DoubleStruckCapitalD],\[Alpha],\[Beta]] Exp[argHEFTMod[\[CapitalDelta]]],{q,0,\[CapitalDelta]}];
Integrand//Expand
];


(* ::Input::Initialization:: *)
SMEFTIntegrand[\[CapitalDelta]_]:=Module[{Integrand},
Integrand=SeriesCoefficient[1/P[q*\[DoubleStruckCapitalD],\[Alpha],\[Beta]] Exp[argSMEFT[\[CapitalDelta]]],{q,0,\[CapitalDelta]}];
Integrand//Expand
];


(* ::Input::Initialization:: *)
HEFTSpurionIntegrand[\[CapitalDelta]_]:=Module[{Integrand},
Integrand=SeriesCoefficient[1/P[q*\[DoubleStruckCapitalD],\[Alpha],\[Beta]] Exp[argHEFTSpurion[\[CapitalDelta]]],{q,0,\[CapitalDelta]}];
Integrand//Expand
];


(* ::Subsection:: *)
(*List for Simplification*)


(* ::Input::Initialization:: *)
mergeList={Vp->V,Vm->V,Vz->V,nd->nf,nu->nf,nL->nf,nQ->nf,uL->Q,QL->Q,dL->Q,uR->Q,dR->Q,uLd->Qbar,QLd->Qbar,uRd->Qbar,dRd->Qbar,dLd->Qbar,n\[Nu]->nf,ne->nf,eRd->Lbar,eR->L,eLd->Lbar,eL->L,LL->L,LLd->Lbar,\[Nu]L->L,\[Nu]Ld->Lbar,\[Nu]R->L,\[Nu]Rd->Lbar};


(* ::Input::Initialization:: *)
mergefurther={L->F,Lbar->F,Q->F,Qbar->F,GL->X,GR->X,WL->X,WR->X,BL->X,BR->X,WmL->X,WmR->X,WpL->X,WpR->X,ZL->X,ZR->X,AL->X,AR->X,H->S,Hd->S,h->S};
numlist={S->1,F->1,X->1,\[DoubleStruckCapitalD]->1,V->1};


(* ::Input::Initialization:: *)
XcheckList={Vp->V,Vm->V,Vz->V,nd->nf,nu->nf,nL->nf,nQ->nf,uL->Q,QL->Q,dL->Q,uR->Q,dR->Q,uLd->Qbar,QLd->Qbar,uRd->Qbar,dRd->Qbar,dLd->Qbar,n\[Nu]->nf,ne->nf,eRd->Lbar,eR->L,eLd->Lbar,eL->L,LL->L,LLd->Lbar,\[Nu]L->L,\[Nu]Ld->Lbar,\[Nu]R->L,\[Nu]Rd->Lbar,GL->X,GR->X,WL->X,WR->X,BL->X,BR->X,WmL->X,WmR->X,WpL->X,WpR->X,ZL->X,ZR->X,AL->X,AR->X};


(* ::Input::Initialization:: *)
(*all2one={uL->1,uLd->1,dL->1,dLd->1,uR->1,uRd->1,dR->1,dRd->1,\[Nu]L->1,\[Nu]Ld->1,eL->1,eLd->1,eR->1,eRd->1,\[Nu]R->1,\[Nu]Rd->1,GL->1,GR->1,WpL->1,WpR->1,WmL->1,WmR->1,ZL->1,ZR->1,AL->1,AR->1,Vp->1,Vm->1,Vz->1,h->1,H->1,Hd->1,V->1,WL->1,WR->1,BL->1,BR->1,LL->1,LLd->1,QL->1,QLd->1,\[DoubleStruckCapitalD]->1};*)


(* ::Input::Initialization:: *)
(*onlyh={uL->0,uLd->0,dL->0,dLd->0,uR->0,uRd->0,dR->0,dRd->0,\[Nu]L->0,\[Nu]Ld->0,eL->0,eLd->0,eR->0,eRd->0,\[Nu]R->0,\[Nu]Rd->0,GL->0,GR->0,WpL->0,WpR->0,WmL->0,WmR->0,ZL->0,ZR->0,AL->0,AR->0,Vp->0,Vm->0,Vz->0};*)


(* ::Subsection:: *)
(*Haar Integrals*)


(* ::Input::Initialization:: *)
HaarIntegralHEFT[term_]:=Module[{x,t1,t2,t3},
If[term===0, Return[0],
t1=Fold[SeriesCoefficient, term*HaarSO4[\[Alpha],\[Beta]],{{\[Alpha],0,-1},{\[Beta],0,-1}}];];
If[t1===0,Return[0],
t2=Fold[SeriesCoefficient, t1*HaarSU3[z1,z2],{{z1,0,-1},{z2,0,-1}}];];
If[t2===0,Return[0],
t3=SeriesCoefficient[ t2*HaarU1[w],{w,0,-1}];];
t3];


(* ::Input::Initialization:: *)
HaarIntegralHEFTSpurion[term_,trunc_]:=Module[{x,t1,t2,t3,t4},
If[term===0, Return[0],
t1=Fold[SeriesCoefficient, term*HaarSO4[\[Alpha],\[Beta]],{{\[Alpha],0,-1},{\[Beta],0,-1}}];];
If[t1===0,Return[0],
t2=Fold[SeriesCoefficient, t1*HaarSU3[z1,z2],{{z1,0,-1},{z2,0,-1}}];];
If[t2===0,Return[0],
t3=SeriesCoefficient[ t2*HaarSU2[s]*SpurionGenFunc[trunc],{s,0,-1}];];
If[t3===0,Return[0],
t4=SeriesCoefficient[ t3*HaarU1[w],{w,0,-1}];];
t4];


(* ::Input::Initialization:: *)
HaarIntegralSMEFT[term_]:=Module[{x,t1,t2,t3,t4},
If[term===0, Return[0],
t1=Fold[SeriesCoefficient, term*HaarSO4[\[Alpha],\[Beta]],{{\[Alpha],0,-1},{\[Beta],0,-1}}];];
If[t1===0,Return[0],
t2=Fold[SeriesCoefficient, t1*HaarSU3[z1,z2],{{z1,0,-1},{z2,0,-1}}];];
If[t2===0,Return[0],
t3=SeriesCoefficient[ t2*HaarSU2[s],{s,0,-1}];];
If[t3===0,Return[0],
t4=SeriesCoefficient[ t3*HaarU1[w],{w,0,-1}];];
t4];


(* ::Subsection:: *)
(*Output Functions*)


(* ::Input::Initialization:: *)
HEFT[massdim_]:=Module[{list,result,rl},
If[massdim>=1 && IntegerQ[massdim],
PrintTemporary["Calculating Hilbert Series..."],
Abort[]];
list=Flatten[HEFTIntegrand[massdim]/.{Plus->List}];
result=Timing[Plus@@ParallelMap[HaarIntegralHEFT,list]];
Print["Time taken: ", result[[1]]];
Which[massdim>4,rl=result[[2]],
massdim < 4, rl=result[[2]]+\[CapitalDelta]H[massdim],
massdim==4,rl=result[[2]]+\[CapitalDelta]H[4]+KineticTermsHEFT];
rl
];


(* ::Input::Initialization:: *)
SMEFT[massdim_]:=Module[{list,result,rl},
If[massdim>=1 && IntegerQ[massdim],
PrintTemporary["Calculating Hilbert Series..."],
Abort[]];
list=Flatten[SMEFTIntegrand[massdim]/.{Plus->List}];
result=Timing[Plus@@ParallelMap[HaarIntegralSMEFT,list]];
Print["Time taken: ", result[[1]]];
If[massdim==4,rl=result[[2]]+\[CapitalDelta]HSMEFT,rl=result[[2]]];
rl
];


(* ::Input::Initialization:: *)
HEFTLinear[massdim_]:=Module[{list,result,trunc=massdim+1,rl,out},
If[massdim>=1 && IntegerQ[massdim],
PrintTemporary["Calculating Hilbert Series..."],
Abort[]];
list=Flatten[HEFTSpurionIntegrand[massdim]/.{Plus->List}];
result=Timing[Plus@@ParallelMap[HaarIntegralHEFTSpurion[#,trunc]&,list]];
Print["Time taken: ", result[[1]]];
rl=Normal[Series[Normal[Series[result[[2]]*(1-U*Ud),{U,0,trunc}]],{Ud,0,trunc}]];
Which[massdim<4,out=Expand[rl]+\[CapitalDelta]HSpurion[massdim],
massdim==4, out=Expand[rl]+\[CapitalDelta]HSpurion[massdim]+KineticTermsHEFTSpurion,
massdim>4, out = rl];
out
];


(* ::Input::Initialization:: *)
HEFTModded[massdim_]:=Module[{list,result,rl},
If[massdim>=1 && IntegerQ[massdim],
PrintTemporary["Calculating Hilbert Series..."],
Abort[]];
list=Flatten[HEFTformfactorIntegrand[massdim]/.{Plus->List}];
result=Timing[Plus@@ParallelMap[HaarIntegralHEFT,list]];
Print["Time taken: ", result[[1]]];
If[massdim==4,rl=result[[2]]+KineticTermsHEFT,rl=result[[2]]];
rl//Expand
];


(* ::Input::Initialization:: *)
HEFTngb[massdim_]:=Module[{list,result,rl},
If[massdim>=1 && IntegerQ[massdim],
PrintTemporary["Calculating Hilbert Series..."],
Abort[]];
list=Flatten[HEFTIntegrandModded[massdim]/.{Plus->List}];
result=Timing[Plus@@ParallelMap[HaarIntegralHEFT,list]];
Print["Time taken: ", result[[1]]];
If[massdim == 4, 
  rl = Expand[(result[[2]] + 
        Coefficient[arg\[CapitalDelta]HMod, q, massdim]) /. {GB -> 
        h*\[DoubleStruckCapitalD]}] + KineticTermsHEFT, 
  rl = (result[[2]] + 
      Coefficient[arg\[CapitalDelta]HMod, q, massdim]) /. {GB -> 
      h*\[DoubleStruckCapitalD]}];
rl // Expand
];


(* ::Input::Initialization:: *)
Counting[expr_]:=Module[{ex2,l1,l2,counting,classlist},
ex2=expr/.mergeList/.mergefurther;
classlist={S,F,X,\[DoubleStruckCapitalD],V};
l1=CoefficientRules[ex2,classlist];
counting=Table[Times@@Table[classlist[[i]]^l1[[len,1,i]],{i,5}]->l1[[len,2]],{len,Length@l1}];
Print["Total operators: ", Plus@@l1[[All,2]]];
counting//MatrixForm
];
