(* ::Package:: *)

(* ::Input::Initialization:: *)
Clear["Global`*"];
mc=1.836;

cc=-4/3;
(*ss=1/4;*)
ss=-3/4;
ccss=cc*ss;
aphscl=0.380175;
alphshyp=1.39568;
a=1.6553;
bexp=0.2204;
b=0.1653;
vc=0.624075;

mu=mc/2;
tau=(2*mu)^bexp/a;


gauss=Exp[-\[Nu]*r^2];
fn=Integrate[gauss^2*4Pi*r^2,{r,0,Infinity},Assumptions->\[Nu]>0]//Sqrt;
fwv[\[Nu]_,r_]:=Evaluate[gauss/fn];
ii=Integrate[fwv[\[Nu]1,r]*fwv[\[Nu]2,r]*4Pi*r^2,{r,0,Infinity},Assumptions->{\[Nu]1>0&&\[Nu]2>0}]
vclmb=Integrate[fwv[\[Nu]1,r]*fwv[\[Nu]2,r]*4Pi*r^2/r,{r,0,Infinity},Assumptions->{\[Nu]1>0&&\[Nu]2>0}];
vln=Integrate[fwv[\[Nu]1,r]*fwv[\[Nu]2,r]*4Pi*r^2*r,{r,0,Infinity},Assumptions->{\[Nu]1>0&&\[Nu]2>0}];
vhyper=Integrate[fwv[\[Nu]1,r]*fwv[\[Nu]2,r]*4Pi*r^2*Exp[-tau^2*r^2],{r,0,Infinity},Assumptions->{\[Nu]1>0&&\[Nu]2>0&&tau>0}];
kntc=Integrate[fwv[\[Nu]1,r]*Laplacian[fwv[\[Nu]2,r],{r,\[Theta],\[Phi]},"Spherical"]*4Pi*r^2,{r,0,Infinity},Assumptions->{\[Nu]1>0&&\[Nu]2>0}];
fii[\[Nu]1_,\[Nu]2_]:=Evaluate[ii];

fvclmb[\[Nu]1_,\[Nu]2_]:=Evaluate[vclmb];
fvvln[\[Nu]1_,\[Nu]2_]:=Evaluate[vln];
fvhyper[\[Nu]1_,\[Nu]2_]:=Evaluate[vhyper];

fv[\[Nu]1_,\[Nu]2_]:=aphscl*fvclmb[\[Nu]1,\[Nu]2]*cc+(-0.75*b*fvvln[\[Nu]1,\[Nu]2]+vc*fii[\[Nu]1,\[Nu]2])*cc-8*Pi*alphshyp/3/mc^2*tau^3/Pi^(3/2)*fvhyper[\[Nu]1,\[Nu]2]*ccss;

ft[\[Nu]1_,\[Nu]2_]:=-1/(2mu)*Evaluate[kntc]+fii[\[Nu]1,\[Nu]2]*2mc;


(* ::Input::Initialization:: *)
nr=20;
rmin=0.0001/0.1973;
rmax=3/0.1973;

lr=Table[rmin*(rmax/rmin)^(i/(nr-1)),{i,0,nr-1}];
l\[Nu]=1/lr^2;

mii=Table[fii[\[Nu]1,\[Nu]2],{\[Nu]1,l\[Nu]},{\[Nu]2,l\[Nu]}];
mhh=Table[fv[\[Nu]1,\[Nu]2]+ft[\[Nu]1,\[Nu]2],{\[Nu]1,l\[Nu]},{\[Nu]2,l\[Nu]}];

Eigenvalues[{mhh,mii}]//Sort
