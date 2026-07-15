(* ::Package:: *)

BeginPackage["GRIMOIRE`"];

(* Declare exported functions *)
COMMANDSFIRE::usage
 (*Directories*)
SetFIREFilesDirectories::usage = "SetDirectories[dirPROG_, dirStart_, dirTables_] updates the working directories.";
ResetFIREFilesDirectories::usage = "ResetDirectories resets the working directories to their default values.";
SetFIREDirectory::usage = "SetFIREDirectory[dirFIRE] updates FIRE Path.";
CurrentFIREFilesDirectories::usage
SetGRIMOIREDirectory::usage
 (*Main Program*)
BUILDMATRICES::usage = "BUILDMATRICES[] constructs matrices based on FIRE6 master integrals and check the matrices.";
START::usage
FIRSTDATA::usage
INCIDENCEMATRIX::usage
 (*Auxiliary functions*)
FILENAME::usage
CLEARFIRE::usage
MAKERULES::usage
CONFIGTEMPLATE::usage
PREPARE::usage
 (*Total Time function*)


(* Declare variables for external access *)
INTERNAL::usage
EXTERNAL::usage
PROPAGATORS::usage
REPLACEMENT::usage
(* Declare variables for internal access *)


(* Default directories *)
If[!ValueQ[dirPROG], dirPROG = NotebookDirectory[]];
If[!ValueQ[dirStart], dirStart = {}];
If[!ValueQ[dirTables], dirTables = {}];
If[!ValueQ[dirStartFolder], dirStartFolder = {}];
If[!ValueQ[dirDerivative], dirDerivative = {}];
If[!ValueQ[dirConfig], dirConfig = {}];
If[!ValueQ[dirMatrix], dirMatrix = {}];
If[!ValueQ[dirMaster], dirMaster = {}];
If[!ValueQ[dirMasterList], dirMasterList = {}];
If[!ValueQ[dirIncidence], dirIncidence = {}];
If[!ValueQ[FIREDirectory], FIREDirectory = {}];
If[!ValueQ[GRIMOIREDirectory], GRIMOIREDirectory = NotebookDirectory[]];


Begin["`Private`"];

VERSION = Style["0.0.00",Bold];
Print["GRIMOIRE, ver. ",VERSION," loaded successfully."];
Print["Use ",Style["COMMANDS",Bold], ", for GRIMOIRE commands."];
Print["---------------------------------------------"];


COMMANDSFIRE := (
	Print[Style["FIRE relate commands:",Bold,Italic,Blue]];
		Print[Style["     BUILDMATRICES[FileName][CHECK]",Bold], " - Entries: FileName: Name of the data to be save; CHECK: True -> Check the matrices, False -> Skip check.\n             This function is required to built the matrices and master integrals files.\n             This function will call FIRE6 and build the matrices and master integral list.\n",Style["             FIRE6 ",Italic,Bold],"(https://gitlab.com/feynmanintegrals/fire) is required!"];
		Print[Style["     SetFIREDirectory[FIRE$Path]",Bold], " - Set FIRE directory. Must be called before BUILDMATRICES[][]."];
		Print[Style["     SetFIREFilesDirectories[PROG$Path,START$Path,TABLES$Path]",Bold], " - Entries:\n             1st: The mains directory of start and tables files; \n             2nd: The subdirectory and name of start file and problem number (see default path);\n             3rd: The subdirectory and name of tables file"];
		Print[Style["     CurrentFIREFilesDirectories",Bold], " - Shows current directories."];
		Print[Style["     ResetFIREFilesDirectories",Bold], " - resets directories to default:\n             dirPROG = NotebookDirectory[];\n             dirStart = {\"FIRE-files/start/master\", 1};\n             dirTables = \"FIRE-files/tables/master.tables\";"];
);

(*Directories*)
	SetGRIMOIREDirectory[vacdiffdir_]:= (
		GRIMOIREDirectory = vacdiffdir;
		Print["GRIMOIRE files directory updated."]
	);
	SetFIREDirectory[FireDirectory_]:= (
		FIREDirectory = FireDirectory;
		Print["FIRE6 directory updated."]
	);
	If[!ValueQ[filename], filename = {}];
	If[!ValueQ[TEMP], TEMP = {}];
	FILENAME[tempfilename_] := (
		dirStartFolder="FIRE-files/start/";
		dirStart={StringJoin["FIRE-files/start/",tempfilename],1};
		dirTables=StringJoin["FIRE-files/tables/",tempfilename,TEMP,".tables"];
		dirDerivative=StringJoin["FIRE-files/list-master/",tempfilename,TEMP,".m"];
		dirConfig=StringJoin["FIRE-files/config/",tempfilename,TEMP,".config"];
		dirMaster=StringJoin["master-integrals/",tempfilename,"-master",".dat"];
		dirMasterList=StringJoin["master-integrals/",tempfilename,"-list",".dat"];
		dirMatrix=StringJoin["matrices/",tempfilename,"-matrices",".dat"];
		dirIncidence=StringJoin["matrices/",tempfilename,"-incidence-matrix",".dat"];
		filename = tempfilename;
	);
 
(*Build Matrices*)
	(*Auxiliary cache values*)
		If[!ValueQ[FIRE6Loaded], FIRE6Loaded = False];
		If[!ValueQ[CachedMaster0], CachedMaster0 = {}];
		If[!ValueQ[CachedMaster], CachedMaster = {}];
		If[!ValueQ[CachedMassU], CachedMassU = {}];
		If[!ValueQ[CachedMatricesA], CachedMatricesA = {}];
(*Main function*)		
BUILDMATRICES := Module[{Master0, Master, ListMaster, mass, MatricesA, AA, results, Crosscheck},
	
	If[FIREDirectory === {}, Print["FIRE6 directory error!\n", "Set the directory of FIRE6 with the following command, after loading the GRIMOIRE Package: SetFIREDirectory[$Path]."];Abort[]];
	Get["FIRE6.m",Path->FIREDirectory];
		Print["Loading FIRE6 tables..."];
		FIRE`LoadStart[dirStart[[1]], dirStart[[2]]];
		FIRE`Burn[];
		FIRE`LoadTables[dirTables];
		Print["Tables loaded successfully."];
		FIREMaster0 = FIRE`MasterIntegrals[];
		FIREMaster = Table[Global`G[x[[1]], x[[2]]], {x, FIREMaster0}];
		Master0 = Table[x[[2]], {x, FIREMaster0}];
		Master = Table[Global`G[x[[2]]], {x, FIREMaster0}];
		mass = Table[Symbol["u" <> ToString[i]], {i, Length[Master0[[1]]]}];
		
			Print["Starting to build the matrices..."];
		AA[i_] := Transpose[Table[Coefficient[x, FIREMaster], {x,Table[y[[2, i]] FIRE`F[FIREMaster0[[1]][[1]], y[[2]] + UnitVector[Length[Master0[[1]]], i]], {y, List @@@ FIREMaster}]}]];
		MatricesA = Table[AA[i], {i, Length[Master0[[1]]]}];

        CachedMaster0 = Master0;
        CachedMaster = Master;
        CachedMassU = mass;
        CachedMatricesA = Simplify[MatricesA];
		results = Table[Count[D[MatricesA[[i1]], mass[[i2]]] - D[MatricesA[[i2]], mass[[i1]]] + (MatricesA[[i1]] . MatricesA[[i2]] - MatricesA[[i2]] . MatricesA[[i1]]) // Simplify, 0, 2] == Length[Master0]^2, {i1, Length[Master0[[1]]]}, {i2, Length[Master0[[1]]]}];
		If[AllTrue[Flatten[results], TrueQ], If[!FIRE6Loaded,Print["Crosscheck Passed.\n", "The matrices have been successfully built."],Print["Crosscheck passed.\n","The matrices have already been successfully built."]], Print[Style["Crosscheck Error.\n","Matrices have been built incorrectly.",Red]];Abort[]];
		Print["Creating files..."];
		Put[CachedMaster0,dirMasterList];
		Print["Master Integrals list file created, path: ", Style[dirMasterList,Bold]," ."];
		Put[CachedMaster,dirMaster];
		Print["Master Integrals file created, path: ",Style[dirMaster,Bold]," ."];
		Put[CachedMatricesA,dirMatrix];
		Print["Matrices file created, path: ",Style[dirMatrix,Bold]," ."];
		Quiet[Remove["FIRE`*"]];
		Print["FIRE6 closed."];
];
CLEARFIRE:= (
	Clear[CachedMaster0];
	Clear[CachedMaster];
	Clear[CachedMassU];
	Clear[CachedMatricesA];
	);
	
START := (	
	If[FIREDirectory==={},Print["FIRE directory not found. Execute SetFIREDirectory[$path]."];Abort[]];
	SetDirectory[dirPROG];
		Get["FIRE6.m",Path->FIREDirectory];
			
		FIRE`Internal = CachedInternMomenta;
		FIRE`External = OriginalExt;
		FIRE`Propagators = CachedProp;
		FIRE`Replacements = CachedReplacement;
		FIRE`PrepareIBP[];
		FIRE`Prepare[];
		FIRE`SaveStart[dirStart[[1]], dirStart[[2]]];
		Print["Start file saved to: ", Style[ dirStart[[1]]<>".start", Bold]];
		GRIMOIRE`CONFIGTEMPLATE[True][{}][Length[CachedProp]];
				        
    Quiet[Remove["FIRE`*"]];
    Print["FIRE6 closed."];
);

MAKERULES[lhs_List, rhs_List] := Module[{bases, rhsTrim},
	bases = DeleteDuplicates @ Cases[lhs, s_Symbol :> s, {1, Infinity}];
	rhsTrim = Take[rhs, Length@bases];
	Thread[bases -> rhsTrim]
];
	(*Auxiliary cache values*)
		If[!ValueQ[CachedIntern], CachedIntern = {}];
		If[!ValueQ[CachedInternRule], CachedInternRule = {}];
		If[!ValueQ[CachedInternMomenta], CachedInternMomenta = {}];
		If[!ValueQ[ReversedInternRule], ReversedInternRule = {}];
		If[!ValueQ[OriginalIntern], OriginalIntern = {}];
INTERNAL[intern__] := Module[{tempintern},
	tempintern={intern};
    CachedIntern = Table[Symbol["k" <> ToString[i]], {i, Length[tempintern]}];
	CachedInternRule = GRIMOIRE`MAKERULES[tempintern,CachedIntern];
	CachedInternMomenta = tempintern /. CachedInternRule;
	ReversedInternRule = Reverse/@CachedInternRule ;
	OriginalIntern= tempintern;
];
		If[!ValueQ[CachedExt], CachedExt = {}];
		If[!ValueQ[CachedExtRule], CachedExtRule = {}];
		If[!ValueQ[CachedExtMomenta], CachedExtMomenta = {}];
		If[!ValueQ[ReversedExtRule], ReversedExtRule = {}];
		If[!ValueQ[OriginalExt], OriginalExt = {}];
EXTERNAL[exter__] := Module[{temp},
	temp=exter;
    CachedExt = Table[Symbol["p" <> ToString[i]], {i, Length[temp]}];
	CachedExtRule = GRIMOIRE`MAKERULES[temp,CachedExt];
	CachedExtMomenta = temp /. CachedExtRule;
	ReversedExtRule = Reverse/@CachedExtRule ;
	OriginalExt= Table[Symbol["p" <> ToString[i]], {i, Length[CachedExtRule]}]
];
	If[!ValueQ[CachedProp], CachedProp = {}];
PROPAGATORS[prop__]:= Module[{temp,temp2,ext,temp3},
	temp = {prop} /. CachedInternRule;
	ext=(temp/. Thread[CachedInternMomenta->0])=!=ConstantArray[0,Length[temp]];
	If[!ext, 
		CachedExtMomenta = {};
		OriginalExt = {},
		Print["External Momenta Detected!"];
		temp2 = DeleteDuplicates[temp/. Thread[CachedInternMomenta->0]];
		GRIMOIRE`EXTERNAL[DeleteCases[temp2,0]];
	];
	temp3 = temp /. CachedExtRule;
	CachedProp = -temp3^2 + Table[Symbol["u" <> ToString[i]], {i, Length[temp]}];
	If[!ext, Print["Internal Momenta:",CachedInternMomenta,"\n","Propagators:", CachedProp],
	Print["Internal Momenta:",CachedInternMomenta,"\n","External Momenta:",OriginalExt ,"\n","Propagators:", CachedProp,"\n","Replacements:",CachedReplacement];];
	If[FileExistsQ[dirIncidence],Print["The incidence matrix for ",Style[filename, Bold], " already exists.\n","File path: ", Style[dirIncidence, Bold]],GRIMOIRE`INCIDENCEMATRIX[CachedInternMomenta,temp3]];
];
	If[!ValueQ[CachedReplacement], CachedReplacement = {}];
	If[!ValueQ[CachedReplacement], CachedReplacement = {}];
REPLACEMENT[replace__] := (CachedReplacement = {replace};
	CachedReplace = DeleteDuplicates@Cases[CachedReplacement[[All,2]],_Symbol,Infinity];);

CONFIGTEMPLATE[temp_:False][master_][length_] := Module[{content, vars, varLine, dir, configline, problemline, integralsline, outputline, folderPath, outputPath, integralsPath},
	If[!temp, TEMP = {};GRIMOIRE`FILENAME[filename], TEMP = "-temp"; GRIMOIRE`FILENAME[filename]];
	n = length;
	vars = Join[{"d"}, Table["u" <> ToString[i], {i, n}],CachedReplace];
	varLine = "#variables         " <> StringRiffle[vars, ","];
	dir = GRIMOIREDirectory;
	folderPath = FileNameJoin[{dir, dirStartFolder}];
	outputPath = FileNameJoin[{dir, dirTables}];
	integralsPath = FileNameJoin[{dir, dirDerivative}];
	configline = "#folder            " <> ToWSLPath[folderPath] <> "/";
	problemline = "#problem           " <>	StringRiffle[{dirStart[[2]], filename <> ".start"}, " "];
	integralsline = "#integrals         " <> ToWSLPath[integralsPath];
	outputline = "#output            " <> ToWSLPath[outputPath];
	content = StringRiffle[{
		"#threads           4",
		"#fthreads          4",
		"#clean_databases",
		varLine,
		"#start",
		configline,
		problemline,
		integralsline,
		outputline
		},
		"\n"
	];
	Export[dirConfig, content, "Text"];
	If[!temp,Put[master, dirDerivative],Put[{{dirStart[[2]],ReplacePart[ConstantArray[2, Length[CachedProp]], -1 -> 1]}},dirDerivative]];
	Print["Config file saved to: ", Style[ dirConfig, Bold]];
];

ToWSLPath[path_] := Module[{p = ToString[path], drive, rest},
  If[StringMatchQ[p, LetterCharacter ~~ ":" ~~ ___],
    drive = ToLowerCase @ StringTake[p, 1];
    rest  = StringDrop[p, 2];
    "/mnt/" <> drive <> StringReplace[rest, "\\" -> "/"],
    p
  ]
];

PREPARE:=(
	If[FIREDirectory==={},Print["FIRE directory not found. Execute SetFIREDirectory[$path]."];Abort[]];
	Get["FIRE6.m",Path->FIREDirectory];
		TEMP = "-temp";
		GRIMOIRE`FILENAME[filename];
		FIRE`LoadStart[dirStart[[1]], dirStart[[2]]];
		FIRE`Burn[];
		FIRE`LoadTables[dirTables];
		MASTER0 = Table[FIRE`MasterIntegrals[][[i]][[2]],{i,1,Length[FIRE`MasterIntegrals[]]}];
		MASTER = Flatten[Table[If[y[[j]]>0,{dirStart[[2]],ReplacePart[y,j->y[[j]]+1]},Nothing],{y,MASTER0},{j,Length[y]}],1];
		TEMP={};
		GRIMOIRE`CONFIGTEMPLATE[False][MASTER][Length[MASTER[[1]][[2]]]];
			
	Quiet[Remove["FIRE`*"]];
	Print["FIRE6 closed."];
);

INCIDENCEMATRIX[loopmomenta_,propagatormomenta_] :=Module[{NP,NL,NV,possiblevertices,positions,vertices,validpartialcolumn,validpartialincidence,validcolumn,validincidence,func,list,valid,mat,mat1,k,listpositions,i},
	k=1;
	NP=Length[propagatormomenta];
	NL=Length[loopmomenta];
	NV=NP-NL+1;
	possiblevertices=DeleteCases[Tuples[{-1,0,1},NP],ConstantArray[0,NP]];

	positions=Position[Product[ Map[Boole[FreeQ[#,k]]&,possiblevertices . propagatormomenta],{k,loopmomenta}],1]//Flatten;
	vertices=possiblevertices[[positions]];
		validpartialcolumn[col_]:=Count[col,-1]<=1&&Count[col,1]<=1;
	validpartialincidence[inc_]:=AllTrue[Transpose[inc],validpartialcolumn[#]&];
	validcolumn[col_]:=((Count[col,1]==1&&Count[col,-1]==1)||Count[col,0]==Length[col]);
	validincidence[inc_]:=AllTrue[Transpose[inc],validcolumn[#]&]&&MatrixRank[inc]==NV-1;

	func[partialincidence_,candidatevertices_]:=Module[{i,x,indices},
	indices=Select[candidatevertices,(validpartialincidence[vertices[[Join[partialincidence,{#}]]]])&];
	{Table[Join[partialincidence,{x}],{x,indices}],Table[Drop[indices,i],{i,1,Length[indices]}]}//Transpose];
	list=func[{},Range[Length[vertices]]];
	For[k=1,k<=NV-1,k++,
	list=Flatten[Table[func[x[[1]],x[[2]]],{x,list}],1]];
	list;
	valid=Select[Table[vertices[[x]],{x,list[[All,1]]}],validincidence[#]&];

	listpositions=Table[Complement[Range[NV],Position[x . propagatormomenta,0]//Flatten],{x,valid}];
i=Flatten[Position[Length/@listpositions,Min[Length/@listpositions]]][[1]];
mat=valid[[i]];
positions=listpositions[[i]];

	mat1=Join[Transpose[mat],Table[UnitVector[NV,i],{i,positions}]]//Transpose;
	mat1=Abs[mat1];
	Put[mat1,dirIncidence];
	Print["Incidence matrix created."];
	Print[mat1 // MatrixForm];
	Print["Incidence matrix saved to:",dirIncidence]
	
];


End[];
EndPackage[];
