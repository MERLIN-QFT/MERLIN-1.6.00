(* ::Package:: *)

BeginPackage["MERLIN`"];

 (*Directories*)
SetMERLINDirectory::usage
 (*Main Program*)
EVALUATE::usage
PREPAREDIFFERENTIAL::usage = "DIFFERENTIAL, diferential equations.";
SERIESEXPANSION::usage = "Series expansion of the matrices";
DIFFSERIESEXPANSION::usage = "Series expansion of the differential equations";
DIAGRAMEVALUATION::usage
IMPLICITSYM::usage
FINDLIMITDIRECTION::usage
BUILDMATRICESRED::usage
INITIALIZE::usage
MASSCONFIG::usage
 (*Auxiliary functions*)
SYMMETRYANALYSIS::usage = "Symmetry[incidence, masses, nu1, nu2] checks if two integrals are equivalent.";
SYMMETRY::usage = "Symmetry[incidence,massconfig] return the new Master list with symmetries applied.";
FINDREDUCEDPOSITION::usage = "Internal usage only! Set the master integral basis of the chosen diagram.";
CLEARCACHE::usage
SAVEDIAGRAM::usage
LOADDIAGRAM::usage
RELOAD::usage
SAVEINCIDENCEMATRIX::usage
EVALUATION::usage
DIAGRAM::usage
LOADRESULTS::usage
INFO::usage
MAKERULES::usage
SIMPLIFY::usage
STARTPARALLEL::usage
 (*Total Time function*)
TIMING::usage


(* Declare variables for external access *)
MASTERLIST = {};
MASTERINT = {};
MASS = {};
MATRICES = {};
MASTERSYM = {};
MASTERRED = {};
DIFFCOEFFICIENTS = {};
INTEXPANSION = {};
MATRIXEXPANSION = {};
DIFFEXPANSION = {};
DIAGRAMRESULT = {};
IMPLICITRULES = {};
MATRICESRED = {};
VECDIAGRAM = {};
RESULTS = {};

(* Declare variables for internal access *)
INCIDENCEMATRIX = {};
DIAGRAMVar = {};
FILENAMECACHE = {}
LIMITDIRECTION = {};

(* Default directories *)
If[!ValueQ[MERLINDirectory], MERLINDirectory = NotebookDirectory[]];


Begin["`Private`"];

VERSION = Style["1.6.00",Bold];
Print["MERLIN, ver. ",VERSION," loaded successfully."];

(*Directories*)
	SetMERLINDirectory[vacdiffdir_]:= (
		MERLINDirectory = vacdiffdir;
		Print["MERLIN files directory updated."]
	);

(*Build Matrices*)
	(*Auxiliary cache values*)
		If[!ValueQ[CachedMaster0], CachedMaster0 = {}];
		If[!ValueQ[CachedMaster], CachedMaster = {}];
		If[!ValueQ[CachedMassU], CachedMassU = {}];
		If[!ValueQ[CachedMatricesA], CachedMatricesA = {}];
	(*Auxiliary cache values*)
		If[!ValueQ[loaded], loaded = False];
INITIALIZE[filename_]:=Module[{filenamemaster,filenamemasterlist,filenamematrix,Master0,Master,MatricesA,mass,filenameim,temp,filedata,data},
	MERLIN`STARTPARALLEL;
	If[!loaded,
		Print["Loading initial data..."];
		SetDirectory[MERLINDirectory];
		filedata=StringJoin["data/",filename,".dat"];
		filenamemaster=StringJoin["master-integrals/",filename,"-master",".dat"];
		filenamemasterlist=StringJoin["master-integrals/",filename,"-list",".dat"];
		filenamematrix=StringJoin["matrices/",filename,"-matrices",".dat"];
		filenameim=StringJoin["matrices/",filename,"-incidence-matrix",".dat"];
		temp = FileExistsQ[filedata];
		If[temp===True, data=Get[filedata]; Print[data]];
		Master0 = Get[filenamemasterlist];
		Master = Get[filenamemaster];
		Print["Master Integrals loaded."];
		MatricesA = Get[filenamematrix];
		Print["Matrices loaded."];
		mass = Table[Symbol["u" <> ToString[i]], {i, Length[Master0[[1]]]}];
		Print["Default mass created: ", mass,"."];
		temp = FileExistsQ[filenameim];
		If[!temp, Print[Style["Incidence matrix file not found. Define the variable INCIDENCEMATRIX.",Red]],
			INCIDENCEMATRIX = Get[filenameim];
			Print["Incidence matrix loaded:", MatrixForm[INCIDENCEMATRIX],"."]
		];
		CachedMaster0 = Master0;
		CachedMaster = Master;
		CachedMassU = mass;
		CachedMatricesA = MatricesA;
		MASTERLIST = CachedMaster0;
		MASTERINT = CachedMaster;
		MASS = CachedMassU;
		MATRICES = Transpose[CachedMatricesA,2<->3];
		FILENAMECACHE = filename;
		
		Print["Initial data succesfully loaded."];
		loaded = True,
		Print["Initial data already loaded."];
     ];
];

RELOAD := (
	Print["Setting reload state..."];
	loaded=False;
	FILENAMECACHE ={};
	CachedMassConfig = {};
	CachedMassRule = {};
	OriginalMassConfig = {};
	ReverseMassRule = {};
	CachedMastersym = {};
	CachedMasterred = {};
	CachedIMPLICITRULE = {};
	MERLIN`CLEARCACHE;
	Print["Reload state set."];
);

SYMMETRYANALYSIS[incidence_,masses_,\[Nu]1_,\[Nu]2_] := Module[{NV,NP0,NP,NE,CandidatePermutationsVertices,CandidatePermutationsEdges,n1,n2,ContractedEdges1,ContractedEdges2,m1,m2,SV,SP,list,B1,B2,contract},
	NP0=Count[Total[incidence],2]; 
	NP=NP0-Count[\[Nu]1,0]; 
	NV=Length[incidence]-Count[\[Nu]1,0]; 
	NE=Count[Total[incidence],1] ;
	If[(NP0==Length[\[Nu]1] &&
		NP0==Length[\[Nu]2]&&
		NP0==Length[masses]&&
		NE+NP0==Length[Transpose[incidence]]&& 
		OrderedQ[Reverse@Total[incidence]])==False, Print["incorrect input"];Return[False]
	] ;
	If [Count[\[Nu]1,0]!=Count[\[Nu]2,0],Return[False]];
	CandidatePermutationsEdges=Permutations[Range[NP]];
	ContractedEdges1=Position[\[Nu]1,0];
	ContractedEdges2=Position[\[Nu]2,0];
	n1=Delete[\[Nu]1,ContractedEdges1];
	n2=Delete[\[Nu]2,ContractedEdges2];
	CandidatePermutationsEdges=Select[CandidatePermutationsEdges,n1[[#]]==n2&];
	If[CandidatePermutationsEdges=={},Return[False]];
	m1=Delete[masses,ContractedEdges1];
	m2=Delete[masses,ContractedEdges2];
	CandidatePermutationsEdges=Select[CandidatePermutationsEdges,m1[[#]]==m2&];
	If[CandidatePermutationsEdges=={},Return[False]];
	CandidatePermutationsEdges=Join[#,Range[NP+1,NP+NE]]&/@CandidatePermutationsEdges;
	CandidatePermutationsVertices=Permutations[Range[NV]];
	
	B1=incidence;
	For[i=1,i<=Count[\[Nu]1,0],i++,
		contract=B1[[All,ContractedEdges1[[i,1]]]];
		If[MemberQ[contract,2],B1={};Break[]];
			B1=Join[{contract},Table[UnitVector[Length[contract],x],{x,Flatten[Position[contract,0]]}]] . B1;
		];
	If[B1!={},B1=Mod[Transpose[Delete[Transpose[B1],Position[\[Nu]1,0]]],2]];
	
	B2=incidence;
	For[i=1,i<=Count[\[Nu]2,0],i++,
		contract=B2[[All,ContractedEdges2[[i,1]]]];
		If[MemberQ[contract,2],B2={};Break[]];
			B2=Join[{contract},Table[UnitVector[Length[contract],x],{x,Flatten[Position[contract,0]]}]] . B2;
		];
	If[B2!={},B2=Mod[Transpose[Delete[Transpose[B2],Position[\[Nu]2,0]]],2]];
	If [B1=={} && B2=={}, Return[True]];
	If[Dimensions[B1]!=Dimensions[B2],Return[False]];
	MemberQ[Table[Transpose[B1[[x]]][[y]]==Transpose[B2],{x,CandidatePermutationsVertices},{y,CandidatePermutationsEdges}],True,2]	
];

	(*Auxiliary cache values*)
		If[!ValueQ[CachedMasterred], CachedMasterred = {}];
		If[!ValueQ[CachedMastersym], CachedMastersym = {}];
(*Symmetry application function*)
SYMMETRY:= Module[{equivalences,Mastersym,Masterred},
	If[CachedMaster0==={}, Print["Master Integral list not found. Execute LOAD[]."];Abort[]];
	If[INCIDENCEMATRIX==={} || CachedMassConfig==={}, Print["Mass Configuration or Incident Matrix not found.\n", "No symmetry applied."],
		equivalences=Gather[CachedMaster0,MERLIN`SYMMETRYANALYSIS[INCIDENCEMATRIX,CachedMassConfig,#1,#2]&]//Sort;
		Mreplacements=Flatten[Table[x->y[[1]],{y,equivalences},{x,y}]];
		Mastersym = CachedMaster/. Mreplacements;
		Masterred = DeleteDuplicates[Mastersym];
		Print["Symmetries applied for Mass configuration: ",OriginalMassConfig,"."];
		CachedMastersym = Mastersym;
		CachedMasterred = Masterred;
		MASTERSYM = CachedMastersym;
		MASTERRED = CachedMasterred;
	];
];

(*Differential function*)
	(*Auxiliary cache values*)
		If[!ValueQ[CachedCoefficientsDIFF], CachedCoefficientsDIFF = {}];
(*Main function*)

PREPAREDIFFERENTIAL:=Module[{e, DN, index, Eqn, result,M, TempMass, tempIndex, TempMatricesA,currentDerivative, totalDerivativeIndex,currentIndex, vecNumber, vecTotal},
   
	If[CachedMatricesA === {}, Print["Matrices not found. Execute LOAD[]."];Abort[]];
	If[DIAGRAMVar === {}, Print["Diagram to be evaluated not found. Define DIAGRAM."];Abort[]];
	
	e[n_]=UnitVector[Length[CachedMaster0],n];
	MERLIN`FINDREDUCEDPOSITION[CachedMaster0,DIAGRAMVar];
	index[y_]:=If[#-1<0,0,#-1]&/@y;
	tempIndex = index[DIAGRAMVar];
	M = Table[Symbol["m" <> ToString[i]], {i, Length[CachedMassConfig]}];
	Do[If [tempIndex[[i]]>= 1, M[[i]]=M[[i]],M[[i]]=0], {i,1,Length[M]}];

	TempMatricesA[i_,M_] := Module[{rule,matricesB},
	rule = Thread[CachedMassU -> CachedMassConfig+Global`t*LIMITDIRECTION + M];
	matricesB = CachedMatricesA[[i]] /.rule;
	matricesB
	];

	DN[i_, x_, n_] := With[{var = M[[i]], mat = TempMatricesA[i, M]},
		Module[{r = x, new, k, j, l},
			vecTotal = Length[x];
			For[k = 1, k <= n, k++,
				currentIndex = i;
				currentDerivative = k;
				totalDerivativeIndex = n;
				new = Table[0, {vecTotal}];
				For[j = 1, j <= vecTotal, j++,
					If[Mod[Floor[100*j/vecTotal],10] == 0 ||j == vecTotal,vecNumber = Floor[100*j/vecTotal];];
					new[[j]] = D[r[[j]], var] - Sum[mat[[j,l]]*r[[l]], {l,1,vecTotal}];];
			r = new;];(-1)^n/n! r]];
	Eqn[q_] := Module[{expr},
		expr = e[masterpos];
			currentDerivative = 0;
			totalDerivativeIndex = 0;
			currentIndex = 0;
			vecNumber = 0;
			vecTotal = Length[CachedMaster0];
		Monitor[Do[If[q[[i]] =!= 0,expr = DN[i, expr, q[[i]]] /. M[[i]] -> 0;M = M /. M[[i]] -> 0;],{i, 1, Length[CachedMaster0[[1]]]}],
			Column[{Row[{currentDerivative, "/", totalDerivativeIndex," Derivative of index ", currentIndex}],ProgressIndicator[vecNumber,{0,100}],Row[{vecNumber,"%"}]}]];
		expr];
		Print["Starting the derivatives for diagram: I [",DIAGRAMVar,"]..."];
	result=AbsoluteTiming[Eqn[tempIndex]];
	If[Length[result[[2]]] =!= Length[CachedMaster],Print[Style["ERROR: Coefficients have been incorrectly built.",Red]];Abort[],
		Print["Derivatives completed, time: ",result[[1]]," seconds."]];
	CachedCoefficientsDIFF = result[[2]];
	MERLIN`TIMING[1,result[[1]]][False];
	DIFFCOEFFICIENTS = CachedCoefficientsDIFF;
];

	(*Auxiliary cache values*)
		If[!ValueQ[masterpos], masterpos = {}];
(*Auxiliary function for Differential function*)
FINDREDUCEDPOSITION[list_,b_]:=Module[{pos,reduceVector},
	pos=Position[list,b/. {x_/;x>0->1},1,1]//Flatten;
	If[pos==={},Print[Style["ERROR: The chosen diagram do not exists!",Red]]; Abort[]  ,masterpos = First[pos]];
];
	(*Auxiliary cache values*)
		If[!ValueQ[CachedMATRIXEXPANSION], CachedMATRIXEXPANSION = {}];
		If[!ValueQ[CachedINTEXPANSION], CachedINTEXPANSION = {}];
(*Series expansion*)

SERIESEXPANSION[k_] := Module[{An1, AS, J, n, time1, time2, time3, progress, stage},

  If[CachedMatricesA === {}, Print["Matrices not found. Execute LOAD[]."]; Abort[]];
  If[LIMITDIRECTION === {}, Print["Velocity vector not found. Execute FINDLIMITDIRECTION."]; Abort[]];
  If[CachedMassConfig === {}, Print["Mass configuration not found. Define MASSCONFIG."]; Abort[]];

  Print["Starting Matrices and Master Integral series expansion, up to leading order: ", k, "..."];
  progress = 0;
  stage = "Starting";
  Monitor[
    J = Table[0, k + 1];
    AS = Table[0, k + 1];
    CachedMATRIXEXPANSION = ConstantArray[0, k + 2];
    stage = "Computing A[-1]";
    progress = 0;
    time1 = AbsoluteTiming[
      An1 =
        SeriesCoefficient[
          Transpose[LIMITDIRECTION . CachedMatricesA] /.
            Thread[CachedMassU -> (CachedMassConfig + Global`t*LIMITDIRECTION)],
          {Global`t, 0, -1}
        ];
      progress = 10;
    ];
    stage = "Expanding matrices";
    time2 = AbsoluteTiming[Do[AS[[n]] =SeriesCoefficient[Transpose[LIMITDIRECTION . CachedMatricesA] /.Thread[CachedMassU -> (CachedMassConfig + Global`t*LIMITDIRECTION)],{Global`t, 0, n - 1}];
		If[Mod[Floor[100*n/(k + 1)], 10] == 0 || n == k + 1,progress = Floor[100*n/(k + 1)];],{n, 1, k + 1}]];
    J[[1]] = CachedMaster;
    stage = "Expanding master integrals";
    progress = 0;
    time3 = AbsoluteTiming[Do[J[[n]] =-Inverse[(n - 1)*IdentityMatrix[Length[CachedMaster]] + An1] . Sum[AS[[i]] . J[[n - i]], {i, 1, n - 1}];
    If[Mod[Floor[100*(n - 1)/k], 10] == 0 || n == k + 1,progress = Floor[100*(n - 1)/k];],{n, 2, k + 1}]];
    Print["Series expansion complete, time: ",time1[[1]] + time2[[1]] + time3[[1]]," seconds."];
    CachedINTEXPANSION = J;
    INTEXPANSION = CachedINTEXPANSION;
    stage = "Saving matrix expansion";
    progress = 0;
    Do[If[i === 1,
       CachedMATRIXEXPANSION[[i]] = An1,
       CachedMATRIXEXPANSION[[i]] = AS[[i - 1]]
      ];
      If[Mod[Floor[100*i/(k + 2)], 10] == 0 || i == k + 2,progress = Floor[100*i/(k + 2)];],{i, 1, k + 2}];
    MATRIXEXPANSION = CachedMATRIXEXPANSION;
    MERLIN`TIMING[3, time1[[1]] + time2[[1]] + time3[[1]]][False],
		Column[{Row[{stage, ": ", progress, "%"}],ProgressIndicator[progress, {0, 100}]}]]
];
	(*Auxiliary cache values*)
		If[!ValueQ[CachedDIFFEXPANSION], CachedDIFFEXPANSION = {}];
DIFFSERIESEXPANSION:= Module[{temp1,temp2,res,time,rule,maxlength,coeffs,needSimplify,count,total,progress},
	If[CachedMatricesA === {}, Print["Matrices not found. Execute LOAD[]."];Abort[]];
	If[CachedCoefficientsDIFF === {}, Print["Differential Coefficients not found. Execute PREPAREDIFFERENTIAL."];Abort[]];
	If[LIMITDIRECTION === {}, Print["Velocity vector not found. Execute FINDLIMITDIRECTION."];Abort[]];
	If[CachedMassConfig === {}, Print["Mass configuration not found. Define MASSCONFIG."];Abort[]];
	Print["Preparing the diagram I[",DIAGRAMVar,"] for series expansion"];
	needSimplify = Length[DeleteDuplicates[CachedMassConfig]] > 1;
	temp1 = AbsoluteTiming[coeffs = If[needSimplify,MERLIN`SIMPLIFY[CachedCoefficientsDIFF],CachedCoefficientsDIFF]];
	Print["Preparation complete, time: ",temp1[[1]]," seconds."];
	Print["Starting series expansion for the diagram: I[",DIAGRAMVar,"], with mass configuration: ",OriginalMassConfig,"..."];
	count = 0;
	progress = 0;
	total = Length[coeffs];
	temp2 = AbsoluteTiming[Monitor[Map[(count++;
        If[Mod[Floor[100*count/total],10] == 0 ||count == total,progress = Floor[100*count/total];];
        CoefficientList[Normal[Series[#, {Global`t,0,0}]],1/Global`t]) &,coeffs],
        Column[{Row[{"Expanding: ", progress, "%"}],ProgressIndicator[progress,{0,100}]}]]];
	maxlength = Max[Map[Length[#] &, temp2[[2]]]];
	time = AbsoluteTiming[Transpose[    Map[PadRight[#, maxlength] &, temp2[[2]]]]];
	CachedDIFFEXPANSION = time[[2]];
		Print["Leading order: ",-(Length[CachedDIFFEXPANSION]-1),"."];
	DIFFEXPANSION = CachedDIFFEXPANSION;
		Print["Series expansion complete, time: ",time[[1]]+temp2[[1]]," seconds."];
	MERLIN`TIMING[2,time[[1]]+temp1[[1]]+temp2[[1]]][False];
];

(*Auxiliary function Total Time*)
	(*Auxiliary cache values*)
		If[!ValueQ[TIME], TIME = Table[0,10]];
(*Total Time*)
TIMING[k_,time_][show_]:= Module[{i},
	If[!show,
		TIME[[k]] = time,
		TIME[[k]] = time;
		Print["Total CPU time: ", Sum[TIME[[i]],{i,1,10}], " seconds."];
	]
];
(*Diagram result*)
	(*Auxiliary cache values*)
		If[!ValueQ[CachedDIAGRAM], CachedDIAGRAM = {}];
DIAGRAMEVALUATION:=Module[{temp,time},
	MERLIN`DIFFSERIESEXPANSION;
	MERLIN`SERIESEXPANSION[Length[CachedDIFFEXPANSION]-1];
	If[CachedMastersym ==={},Print["Continuing evaluation without symmetry..."];CachedMastersym = CachedMaster;CachedMasterred = CachedMaster];
		Print["Preparing final result..."];
	temp = AbsoluteTiming[Sum[CachedDIFFEXPANSION[[i]] . CachedINTEXPANSION[[i]],{i,1,Length[CachedDIFFEXPANSION]}]/.Thread[CachedMaster->CachedMastersym]];
	
	time = AbsoluteTiming[CachedDIAGRAM = temp[[2]] /. CachedIMPLICITRULE
			];
	DIAGRAMRESULT = CachedDIAGRAM;
		Print["Final result complete, time: ",temp[[1]]+time[[1]]," seconds."];
	MERLIN`TIMING[4,temp[[1]]][True]
];

	(*Auxiliary cache values*)
		If[!ValueQ[CachedIMPLICITRULE], CachedIMPLICITRULE = {}];

IMPLICITSYM:= Module[{temp,tempvec},
	If[CachedMATRIXEXPANSION==={},Print["Matrices series expansion and Master series expansion not found, execute SERIESEXPANSION[k]."];Abort[]];
	temp = Simplify[CachedMATRIXEXPANSION[[1]] . CachedINTEXPANSION[[1]]]/.Thread[CachedMaster->CachedMastersym];
	If[temp=!=ConstantArray[0, Length[temp]],Print["Implicit symmetry found:"];
		pos = DeleteCases[DeleteDuplicates[First/@Position[temp,_?(Not[PossibleZeroQ[#]]&)] ],0];
		tempvec = Table[0,Length[pos]];
		Do[tempvec[[i]]=Simplify[Solve[temp[[pos[[i]]]]==0,CachedMastersym[[pos[[i]]]]]],{i,Length[pos]}];
		CachedIMPLICITRULE = Flatten[DeleteDuplicates[tempvec]];Print[CachedIMPLICITRULE/. ReversedMassRule],
			Print["No Implicit symmetry found."];CachedIMPLICITRULE = {x_->x};
	];
	IMPLICITRULES = CachedIMPLICITRULE;
];

CLEARCACHE:= (
	Clear[CachedDIAGRAM];
	Clear[CachedDIFFEXPANSION];
	Clear[CachedCoefficientsDIFF];
	Clear[CachedMATRIXEXPANSION];
	Clear[CachedINTEXPANSION];
	Clear[CachedMatricesB];
	Clear[TIME];
  
	CachedDIAGRAM = {};
	CachedDIFFEXPANSION = {};
	CachedCoefficientsDIFF = {};
	CachedINTEXPANSION = {};
	CachedMATRIXEXPANSION = {};
	CachedMatricesB = {};
	TIME = Table[0,10];
);

(*Find LIMITDIRECTION*)
FINDLIMITDIRECTION := Module[{n, TESTVECTOR, t, indeterminateFound, positions, 
   attemptsLeft, possibleVectors, usedVectors, selectedVector, minNZ}, 
	If[CachedMatricesA === {}, Print["Matrices not found. Execute LOAD[]."];Abort[]];
	If[CachedMassConfig === {}, Print["Mass configuration not found. Define MASSCONFIG."];Abort[]];
	
	n = Length[CachedMassConfig];
	TESTVECTOR = Table[0, {n}];
	TESTVECTOR[[1]] = 1;
	attemptsLeft = 100;
	indeterminateFound = True;
	usedVectors = {TESTVECTOR};

  While[indeterminateFound && attemptsLeft > 0,
    attemptsLeft--;
    positions = Quiet[
      Position[Simplify[CachedMatricesA /. Thread[CachedMassU -> CachedMassConfig + Global`t*TESTVECTOR]], Indeterminate]
    ];
    If[positions === {},
      LIMITDIRECTION = TESTVECTOR;
      Print["Limit direction created: ", LIMITDIRECTION, "."];
      Return[Null], 
      possibleVectors = DeleteDuplicates[
        Join[
          Table[ReplacePart[TESTVECTOR, i -> #] & /@ {-1, 1}, {i, 1, n}] // Flatten[#, 1] &,
          Table[ReplacePart[TESTVECTOR, {i -> -1, j -> 1}], {i, 1, n}, {j, i + 1, n}] // Flatten[#, 1] &
        ]
      ];
      possibleVectors = Complement[possibleVectors, usedVectors];
      possibleVectors = SortBy[possibleVectors, Count[#,_?(# != 0 &)] &];
      If[possibleVectors =!= {},
        minNZ = Count[First[possibleVectors], _?(# != 0 &)];
        possibleVectors = Select[possibleVectors, Count[#,_?(# != 0 &)] == minNZ &];
      ];
      If[possibleVectors =!= {},
        selectedVector = First[possibleVectors];
        AppendTo[usedVectors, selectedVector];
        TESTVECTOR = selectedVector,
        Return[Null] 
      ];
    ];
  ];
  Return[Null]; 
];

SAVEDIAGRAM:=Module[{diagramStr,massConfigStr,filename,diagram,temp,tempname,sanitizeFilename},
	SetDirectory[MERLINDirectory];
	sanitizeFilename[str_String] := StringReplace[str,{"*" -> " ","\\" -> "","/"  -> "",":"  -> "","?"  -> "","\"" -> "","<"  -> "",">"  -> "","|"  -> ""}];
	diagramStr=sanitizeFilename @ StringReplace[ToString[DIAGRAMVar,InputForm]," "->""];
	massConfigStr=sanitizeFilename @ StringReplace[ToString[CachedMassConfig,InputForm]," "->""];
	filename=StringJoin["results/",FILENAMECACHE,diagramStr,massConfigStr,".dat"];
	tempname = FileExistsQ[filename];
	If[!tempname,
		If[CachedDIAGRAM==={},"No result found. Aborting save function.";Abort[]];
		Print["Preparing the result to be saved..."];
		temp = AbsoluteTiming[Expand[Flatten[CachedDIAGRAM]]];
		diagram = AbsoluteTiming[Collect[temp[[2]],CachedMasterred]];
		Print["Complete, time: ",diagram[[1]]+temp[[1]]," seconds."];
		Print["Saving Diagram result..."];
		Put[diagram[[2]],filename];
		Print["Save complete, filename: ",filename,"."],
		Print["File alerady exist."]
	];
];

LOADDIAGRAM[choosediagram_,massconfig_]:=Module[{diagram,diagramStr,massConfigStr,filename,sanitizeFilename},
	SetDirectory[MERLINDirectory];
	sanitizeFilename[str_String] := StringReplace[str,{"*" -> " ","\\" -> "","/"  -> "",":"  -> "","?"  -> "","\"" -> "","<"  -> "",">"  -> "","|"  -> ""}];
	Print["Loading Diagram result: I", choosediagram," with mass configuration: ",OriginalMassConfig,"..."];
		diagramStr=sanitizeFilename @ StringReplace[ToString[choosediagram,InputForm]," "->""];
		massConfigStr=sanitizeFilename @ StringReplace[ToString[massconfig,InputForm]," "->""];
		filename=StringJoin["results/",FILENAMECACHE,diagramStr,massConfigStr,".dat"];
		diagram = Get[filename];
	Print["Result sucessfully loaded."];
		Return[diagram];
];
SAVEINCIDENCEMATRIX[filename_] := Module[{filenameim},
	SetDirectory[MERLINDirectory];
	filenameim=StringJoin["matrices/",filename,"-incidence-matrix",".dat"];
	Print["Saving Incidence Matrix..."];
	Put[INCIDENCEMATRIX,filenameim];
	Print["Save complete, filename: ",filenameim,"."];
];
	If[!ValueQ[LOADFUNC], LOADFUNC = {}];
	If[!ValueQ[EVALUATED], EVALUATED = True];

EVALUATION[save_:False] := Module[{diagramStr,massConfigStr,filename,temp,sanitizeFilename},
		If[!EVALUATED, MERLIN`CLEARCACHE];
	sanitizeFilename[str_String] := StringReplace[str,{"*" -> " ","\\" -> "","/"  -> "",":"  -> "","?"  -> "","\"" -> "","<"  -> "",">"  -> "","|"  -> ""}];
	diagramStr=sanitizeFilename @ StringReplace[ToString[DIAGRAMVar,InputForm]," "->""];
	massConfigStr=sanitizeFilename @ StringReplace[ToString[CachedMassConfig,InputForm]," "->""];
	filename=StringJoin["results/",FILENAMECACHE,diagramStr,massConfigStr,".dat"];
	temp = FileExistsQ[filename];
	If[!temp,	
		Print[Style[Row[{  "Starting the evaluation of the diagram: I[", DIAGRAMVar, "], with mass configuration: ", OriginalMassConfig, "..."}], Bold]];
		MERLIN`PREPAREDIFFERENTIAL;
		MERLIN`DIAGRAMEVALUATION;
		If[save===True, MERLIN`SAVEDIAGRAM];
		LOADFUNC = False,
		Print["The result for the diagram: I[",DIAGRAMVar,"], with mass configuration: ",OriginalMassConfig," already exists."];
		LOADFUNC = True;
	];
	EVALUATED = False;
];
EVALUATE := Module[{temp,mult},
		i=1;
		Do[DIAGRAMVar = VECDIAGRAM[[i]];
		If[Length[CachedMassU]=!=Length[VECDIAGRAM[[i]]], Print["The diagram I[",VECDIAGRAM[[i]],"] does not correspond to \"",FILENAMECACHE,"\" topology."];Abort[]];
		If[Length[CachedMassConfig]=!=Length[VECDIAGRAM[[i]]], Print["The diagram I[",VECDIAGRAM[[i]],"] does not correspond to mass configuration: ",OriginalMassConfig,"."];Abort[]];
		MERLIN`EVALUATION[True],{i,1,Length[VECDIAGRAM]}];
		MERLIN`LOADRESULTS;
];
	
DIAGRAM[diagrams__List] := Module[{temp},
	temp = {diagrams};
	VECDIAGRAM = temp;
];

LOADRESULTS:= Module[{diagram,temp},
	temp = Dimensions[VECDIAGRAM];
	If[Length[temp] === 2, diagram=ConstantArray[0, Length[VECDIAGRAM]];
	Do[diagram[[i]]=MERLIN`LOADDIAGRAM[VECDIAGRAM[[i]],CachedMassConfig];Print["The result for the diagram can be found in variable ", Style["RESULTS",Bold],Style["[[",Bold],i,Style["]]",Bold],"."],{i,1,Length[VECDIAGRAM]}],
	diagram=MERLIN`LOADDIAGRAM[VECDIAGRAM,CachedMassConfig];Print["The result for the diagram can be found in variable ", Style["RESULTS",Bold],"."];
	];
	RESULTS = diagram/.ReversedMassRule;
];
INFO:= Module[{k},
	k=1;
	MERLIN`FINDLIMITDIRECTION;
	Block[{Print = (Null &)},
	MERLIN`SERIESEXPANSION[k]];
	MERLIN`SYMMETRY;
	MERLIN`IMPLICITSYM;
];
	If[!ValueQ[CachedMassConfig], CachedMassConfig = {}];
	If[!ValueQ[CachedMassRule], CachedMassRule = {}];
	If[!ValueQ[OriginalMassConfig], OriginalMassConfig = {}];
	If[!ValueQ[ReversedMassRule], ReverseMassRule = {}];
MASSCONFIG[massconfig__] := Module[{tempmass},
	tempmass={massconfig};
	CachedMassRule = MERLIN`MAKERULES[tempmass,CachedMassU];
	CachedMassConfig = tempmass /. CachedMassRule;
	ReversedMassRule = Reverse/@CachedMassRule ;
	OriginalMassConfig = tempmass;
	MERLIN`INFO
];
MAKERULES[lhs_List, rhs_List] := Module[{bases, rhsTrim},
	bases = DeleteDuplicates @ Cases[lhs, s_Symbol :> s, {1, Infinity}];
	rhsTrim = Take[rhs, Length@bases];
	Thread[bases -> rhsTrim]
];

SIMPLIFY[expr_List] := ParallelMap[
	Function[z,
		If[z === 0 || LeafCount[z] < 5000, z, Cancel[Together[z]]]], 
		expr,
		Method -> "CoarsestGrained"
	];
STARTPARALLEL := Module[{},
	If[$KernelCount == 0, LaunchKernels[]];Print["Parallel kernels active: ", $KernelCount];
	DistributeDefinitions[MERLIN`SIMPLIFY];
	ParallelEvaluate[$HistoryLength = 0];
];


End[];
EndPackage[];
