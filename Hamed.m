
(*************************************************************************)
(*************************************************************************)
(**** : HamedPartitions.m:                            ********************)
(**** : Author:  Maxie D. Schmidt (maxieds@gmail.com) ********************)
(**** : Created: 07.11.2018                           ********************) 
(*************************************************************************)
(*************************************************************************)

LocalPartitionsPackageName = "PartitionsIntoParts`";
BeginPackage[LocalPartitionsPackageName, {"Combinatorica`"}]

(**** : Clear out any old definitions and stored data on reload: ****)
Clear[HPByPartitions, HPByGF, HP, SetComputationMethod, HamedPartitions, PartitionRank];
Clear[GetPartitionT, IdentifyHPInequality, CheckHPDiff, TestHPDiff, PrintHP]; 

(***************************************************************************)
(**** : Package configuration for the core functionality offered here : ****) 
(***************************************************************************)

COMPMETHOD_SETLENGTH::usage = "Named constant to denote that the partition function HP should be computed by set length " <> 
                              "using Partitions and the Combinatorica package. This method may not scale well " <> 
                              "when n is large due to the amount of memory required to store a large array of lists.\n" <>
                              "See also COMPMETHOD_GENFUNC, COMPMETHOD_ADAPTIVE, and SetComputationMethod[...].";
COMPMETHOD_SETLENGTH = 0;

COMPMETHOD_GENFUNC::usage = "Named constant to denote that the partition function HP should be computed by taking " <> 
                            "coefficients of the truncated (to order n) reciprocal product generating function for " <> 
                            "HP_{p,a}(n). This method should in principle be more efficient than the related " <> 
                            "COMPMETHOD_SETLENGTH method in the package..\n" <>
                            "See also COMPMETHOD_SETLENGTH, COMPMETHOD_ADAPTIVE, and SetComputationMethod[...].";
COMPMETHOD_GENFUNC = 1;

COMPMETHOD_ADAPTIVE::usage = "A hybrid of the two methods COMPMETHOD_SETLENGTH and COMPMETHOD_GENFUNC which " <> 
                             "prefers the latter method when n is large-ish, say (by informal heuristic only) " <> 
                             "n \[RightTriangleEqual] 25 to be safe. This parameter can be fine tuned later.\n" <> 
                             "See also SetComputationMethod[...] and the other COMPMETHOD_* named constants above. " <> 
                             "To set the threshold of n values between which we distinguish the particular " <> 
                             "computation method to use, see the local package variable AdaptiveMethodThreshold.";
COMPMETHOD_ADAPTIVE = 4;

Begin["`Private`"];
     PartitionFunctionComputationMethod = COMPMETHOD_ADAPTIVE;
     AdaptiveMethodThreshold = 25;
End[];

SetComputationMethod::usage = "Default settings are SetComputationMethod[COMPMETHOD_ADAPTIVE, 25]. See the " <> 
                              "COMPMETHOD_* named constants usage for additional possibilities.";
SetComputationMethod[compMethod_, adaptiveThreshold_:25] := Block[{}, 
     If[!MemberQ[{COMPMETHOD_SETLENGTH, COMPMETHOD_GENFUNC, COMPMETHOD_ADAPTIVE}, compMethod] || adaptiveThreshold <= 0, 
          PrintError["Invalid parameters. See ?SetComputationMethod for usage instructions."];
          Return[-1];
     ];
     Private`PartitionFunctionComputationMethod = compMethod;
     Private`AdaptiveMethodThreshold = adaptiveThreshold;
     Return[compMethod];
];


(***************************************************************************)
(**** : The crux of the key new functionality offered by the package  : ****) 
(**** : Namely, the code to generate the partitions of our desired    : ****)
(**** : forms and count them appropriately. Partition statistics and  : ****)
(**** : specialized functions for exploratory computations with       : ****)
(**** : Turan-type inequalities and sums of the rank statistic are    : ****)
(**** : implelemented below.                                          : ****)
(***************************************************************************)
HP::usage = "Hamed's special partition function: \n" <> "
            "Computes the number of partitions of n of the form pt+a for fixed primes p and " <> 
            "0 \[LeftTriangleEqual] a < p. The computation method of these integer counts can be changed " <> 
            "to affect the performance speed of the calculations using the function " <> 
            "SetComputationMethod[COMPMETHOD_SETLENGTH|COMPMETHOD_GENFUNC|COMPMETHOD_ADAPTIVE].\n" <> 
            "See also HPByPartitions and HPByGF.";
HP[p_, a_, n_] := HP[p, a, n] = Module[{hpValue}, 
     hpValue = 
     Which[Private`PartitionFunctionComputationMethod == COMPMETHOD_SETLENGTH, HPByPartitions[p, a, n], 
           Private`PartitionFunctionComputationMethod == COMPMETHOD_GENFUNC, HPByGF[p, a, n], 
           Private`PartitionFunctionComputationMethod == COMPMETHOD_ADAPTIVE 
           && n >= AdaptiveMethodThreshold, HPByGF[p, a, n], 
           True, HPByPartitions[p, a, n]];
     Return[hpValue];
];

HPByPartitions::usage = "Same as HP[p, a, n] computed with Combinatorica's Partitions function. In other words, " <> 
                        "the partition count defining this special partition function is generated by " <> 
                        "first computing the explicit partition sets and then taking the length of this set of sets. " <> 
                        "Why is this useful is it doesn't scale well for large n? Well, it implicitly stores the " <> 
                        "computed partitions associated with HP[p, a, n] by Mathematica's built-in dynamic " <> 
                        "programming method, which then means that you are effectively computing both HP[p, a, n] and " <> 
                        "HamedPartitions[p, a, n] at the same time. This is obviously useful if you intend on " <> 
                        "recovering not only the partition function numbers but additional information about the " <> 
                        "associated partitions later, such as, for example, their rank.\n" <> 
                        "For obvious speed and memory consumption reasons use of this function is " <> 
                        "*deprecated* in place of HPByGF[p, a, n] for large-ish n unless the user truely needs to " <> 
                        "recover said more detailed structure of the underlying partitions for later analysis.";
HPByPartitions[p_, a_, n_] := HPByPartitions[p, a, n] = 
     Return[Length[HamedPartitions[p, a, n]]];
];

HPByGF::usage = "Same as HP[p, a, n] computed by means of a q-Pochhammer product generating function.";
HPByGF[p_, a_, n_] := HPByGF[p, a, n] = 
     SeriesCoefficient[Power[q, a] * QPochhammer[Power[q, a], Power[q, p]], {q, 0, n}];

HPByGF2::usage = "Same as HP[p, a, n] computed by means of a truncated reciprocal product generating function. " <> 
                 "It is not clear whether this is more efficient than HPByGF[p, a, n], but I suspect that the " <> 
                 "former function is probably less error prone since it uses the built-in Mathematica routine for the " <> 
                 "full infinite q-Pochhamer symbol. At any rate, this generating function method is also implemented " <> 
                 "here for testing, verification, and comparison.";
HPByGF2[p_, a_, n_] := HPByGF2[p, a, n] = 
     SeriesCoefficient[1 / Product[1 - Power[q, p * t + a], {t, 0, n + 1}], {q, 0, n}]

HamedPartitions::usage = "Returns the actual components of the partitions of said form as a list of lists. " <> 
                         "This function is essentially the analog to Partitions[n] from Combinatorica, " <> 
                         "except that we distinguish that the key components must be of the form pt+a.\n" <>
                         "See also the related function PrintPartitionStats[p, a, n] for a visual summary of these " <> 
                         "associated partitions and their summary statistics.\n" <> 
                         "HP[p, a, n] somewhat in contrast only requires the length of the list returned by this " <> 
                         "function call. As such, if/when HPByPartitions[p, a, n] is called these partitions are " <> 
                         "automatically precomputed in the calculation of the partition number so that the " <>
                         "more detailed structure of the partitions can be analyzed later.\n" <>
                         "To make sure this function is called every time HP[p, a, n] is called, run the " <> 
                         "following sequence of package configuration commands:\n" <> 
                         "(1) SetComputationMethod[COMPMETHOD_SETLENGTH]\n" <>
                         "(** OR ALTERNATELY: **)\n" <>
                         "(2) SetComputationMethod[COMPMETHOD_ADAPTIVE, Infinity].";
HamedPartitions[p_, a_, n_] := HamedPartitions[p, a, n] = 
     Module[{otf, ofTheFormFunc, partsOfTheForm}, 
     otf = IntegerQ[(# - a) / p]&; 
     ofTheFormFunc = Union[Map[otf, #]]&; 
     partsOfTheForm = Map[ofTheFormFunc, Partitions[n]];
     partsOfTheForm = Select[partsOfTheForm, #=={True}&];
     Return[partsOfTheForm];

];

Begin["`Testing`"];

UnitTestingFeaturesUsage = "Used for internal testing, verification, and what we will loosely call our " <> 
                           "first-order approximation to Unit Tests in Mathematica. " <> 
                           "Probably wise not to call this function directly unless you know what you are " <> 
                           "doing and why.\n" <> 
                           "Consider running the public wrapper function RunUnitTests[...] for a sanity check " <> 
                           "that this package is minimally computing what we expect it to compute instead.";

CheckHPDiff::usage = UnitTestingFeaturesUsage;
CheckHPDiff[p_, a_, n_] := HPByPartitions[p, a, n] == HPByGF[p, a, n] && HPByPartitions[p, a, n] == HPByGF2[p, a, n]

CheckHPDiff2::usage = UnitTestingFeaturesUsage;
CheckHPDiff2[p_, n_] := Sum[HP[p, a, n], {a, 0, p - 1}] == PartitionsP[n]

CheckHPDiff3::usage = UnitTestingFeaturesUsage;
CheckHPDiff3[p_, n_] := Sum[HP[p, a, n], {a, 1, p - 1}, 2] == PartitionsQ[n]

CheckTableValid[tableInput_] := SameQ[Union[Flatten[tableInput]], {True}];

TestHPDiff::usage = UnitTestingFeaturesUsage;
TestHPDiff[numPrimes_:15, nupper_:50] := 
     CheckTableValid[Table[CheckHPDiff[p, a, n], {p, Table[Prime[m], {m, 1, numPrimes}, {a, 0, p - 1}, {n, 1, nupper}]];

End[]; (* End Testing Section *)

RunUnitTests::usage = "A basic internal package sanity check on the correctness of the core functions " <> 
                      "which we have coded / implemented above. This is the pretty-print notebook wrapper " <> 
                      "around all of the numbered functions in the `Testing` subpackage. This is the " <> 
                      "high-level function you want to call to verify that the package is configured " <> 
                      "properly and generating at least minimally correct results.";
RunUnitTests[numPrimes_:15, nupper_:50] := Module[{}, 
     test1Result = "[Test 1] Verifying HP[p, a, n] values match across all computation methods: " <> 
                   getUnitTestPassString[Testing`TestHPDiff[numPrimes, nupper]];
     test2Result = "[Test 2] Verifying HP[p, a, n] values against known functions (sum over all a yields p(n)): " <> 
                   getUnitTestPassString[Testing`TestHPDiff2[numPrimes, nupper]];
     test3Result = "[Test 2] Verifying HP[p, a, n] values against known functions (sum over all odd a yields q(n)): " <> 
                   getUnitTestPassString[Testing`TestHPDiff3[numPrimes, nupper]];
     testResultsInit = {test1Result, test2Result, test3Result};
     bulletPoints = {"\[ClubSuit]", "\[SpadeSuit]", "\[BlackQueen]"};
     testResultsDisplayString = MapIndexed[StringJoin[bulletPoints[[First[#1]]], #1]&, testResultsInit];
     PrintNotebookNotification[Blend[Red, Yellow]] @@ {testResultsDisplayString};
];

(***************************************************************************)
(**** : Partition function summary statistics and descriptive         : ****) 
(**** : (verbose) notebook printing routines for in depth analysis    : ****)
(**** : of the component partitions counted to construct HP[p, a, n]. : ****)
(***************************************************************************)

PartitionRank::usage = "Standard summary statistic for an individual partition: MaxPart - NumberOfParts";
PartitionRank[part_] := Max[part] - Length[part];

PartitionOnes::usage = "Summary statistic: The number of parts in the partition which are equal to one.";
PartitionOnes[part_] := Length[Select[part, # == 1&]];

PartitionMu::usage = "Summary statistic used in the computation of the *crank* of a partition (see below).";
PartitionMu[part_] := Length[Select[part, # > PartitionOnes[part]&]];

PartitionCrank::usage = "Standard summary statistic for an individual partition introduced by Freeman Dyson.";
PartitionCrank[part_] := With[{ellLength = Max[part], omegaOnes = PartitionOnes[part], muCount = PartitionMu[part]}, 
     If[omega == 0, ellLength, muCount - omegaOnes]
]

GetPartitionIndexParameter::usage = "The parameter t if the component represents pt+a.";
GetPartitionIndexParameter[p_, a_, comp_] := (comp - a) / p;

GetSinglePartitionStats[p_, a_, part_, includeFormatting_:True] := 
Module[{partitionStats, partSummarySpec, outputDesc}, 
     partitionStats = {p, GetPartitionIndexParameter[#], a, Plus @@ #, 
                       PartitionRank[#], PartitionCrank[#], Min[#], Max[#], 
                       PartitionOnes[#], PartitionMu[#]}&[part];
     partSummarySpec = "\[Lambda] = `1`\[FilledVerySmallCircle]`2`+`3` \[LeftGuillmet] `4`; " <> 
                       "rank=`5`, crank=`6`, spart=`7`, lpart=`8`, ones=`9`, \[Mu](\[Lambda])=`10`";
     outputDesc = ToString[StringForm[partSummarySpec, ##]]& @@ partitionStats;
     If[includeFormatting, 
          outputDesc = " \[FivePointedStar] " <> outputDesc <> "\n";
     ];
     Return[outputDesc];
];

PrintPartitionStats[p_, a_, n_] := Module[{headerString, summaryString}, 
     headerString = ToString[StringForm["Partitions of n = `1` into parts of the form `2`t+`3`:\n\n", n, p, a]];
     headerString = headerString <> 
                    ToString[StringForm["HP(p=`1`, a=`2`; n=`3`) = `4`", p, a, n, HPByPartitions[p, a, n]]];
     headerString = headerString <> 
                    ToString[StringForm["For reference: p(`1`)=`2`, q(`1`)=`3`.\n", n, PartitionsP[n], PartitionsQ[n]]];
     headerString = headerString <> "\n==========================================================================\n\n";
     summaryString = StringJoin @@ Map[GetSinglePartitionStats, HamedPartitions[p, a, n]];
     fullDescString = headerString <> summaryString;
     PrintNotebookNotification[RGBColor[0.27, 0.71, 1.0]] @@ {fullDescString};
]; 


(*************************************************************************)
(**** : Exploratory procedures used to experiment with sums of the  : ****) 
(**** : ranks of the partitions of these forms and congruences      : ****)
(**** : satisfied thereof. Used by Hamed and Maxie post conference. : ****)
(*************************************************************************)
GetPartitionRankSum::usage = "The sum of the ranks of all component partitions of n into parts of " <> 
                             "the form pt+a. Used experimentally for exploratory analysis by Hamed and Maxie.";
GetPartitionRankSum[p_, a_, n_] := GetPartitionRankSum[p, a, n] = 
     Plus @@ Map[PartitionRank, HamedPartitions[p, a, n]];
     
GetPartitionRankSumModuloBases::usage = "Computes all relevant primes p for which HP[p, a, n] === 0 (mod p), i.e., " <> 
                                        "all such primes in the range 2 <= p <= HP[p, a, n].";
GetPartitionRankSumModuloBases[p_, a_, n_] := With[{hp = HP[p, a, n]}, 
     Map[First, Select[Table[{p, Mod[hp, p, 0] == 0}, {p, Table[Prime[m], {m, 1, hp}]}], #[[1]] <= hp && #[[2]]&]];
];

(*************************************************************************)
(**** : Exploratory procedures used to experiment the convexity of  : ****)
(**** : Turan-type inequalities. In analog to Ono, Rolen, and       : ****)
(**** : Zagier's semi-recent papers proving the convexity of the    : ****)
(**** : ordinary partition function p(n) for sufficiently large n.  : ****)
(**** : This is for Hamed, since Maxie personally gets bored and    : ****)
(**** : incomprensibly stiffled by working with these inequalities  : ****)
(**** : in the context of partitions :)                             : ****)
(*************************************************************************)

IdentifyHPInequality::usage = "Determines if the convexity condition on the Turan-type inequality is " <> 
                              "satisfied at this n: HP(p, a; n) \[GreaterEqual] HP(p, a; n-1) HP(p, a; n+1)??? " <> 
                              "Also generates a usable summary string description of the direction of the " <> 
                              "inequality. The returned value is of the form " <> 
                              "{IsConvex -> True|False, InequalityDirection -> directionString}.";
IdentifyHPInequality[p_, a_, n_] := Module[{p2, p1, p0, ttIneqLHS, ttIneqRHS, isConvex, areEqual, ineqDirStr}, 
     p0 = HP[p, a, n - 1];
     p1 = HP[p, a, n];
     p2 = HP[p, a, n + 1];
     ttIneqLHS = p2; 
     ttIneqRHS = p0 * p1;
     isConvex = ttIneqLHS >= ttIneqRHS;
     isEqual = ttIneqLHS == ttIneqRHS;
     ineqDirStr = ToString[
                       StringForm["`1` `2` `3`", ttIneqLHS, 
                       Which[isEqual, "\[Congruent]", 
                             isConvex, "\[RightTriangleEqual]", 
                             True, "\[NotRightTriangleEqual]"], ttIneqRHS]
                  ];
     Return[{IsConvex -> isConvex, InequalityDirection -> ineqDirStr}];
];
                  


(*************************************************************************)
(**** : Text processing, display and general utility functions      : ****) 
(*************************************************************************)

PrintNotebookNotification[textColor_, backgroundColor_, borderColor_] := 
If[$Notebooks,
     CellPrint[Cell[packageAnnouncement, "Print", 
                    FontColor -> textColor, 
                    FontSize -> 14, 
                    CellFrame -> 2.5, 
                    CellFrameColor -> borderColor, 
                    Background -> backgroundColor, 
                    CellFrameMargins -> 14, 
                    ShowCellBracket -> False
                   ]
              ],
     Print[#1]; 
];

PrintNotebookNotification[baseColor_] := 
     PrintNotebookNotification[Darker[baseColor], Lighter[baseColor], Glow[baseColor]];

PrintError[errorText_] := 
     PrintNotebookNotification @@ StringJoin["\[ScriptX]\n", " !!! ERROR : >>> \n", errorText, "\n\[ScriptX]"]

(*************************************************************************)
(**** : Package help and information printing routines              : ****) 
(*************************************************************************)



(*************************************************************************)
(**** : Perform pretty printing of the package details and          : ****) 
(**** : revsision information if the package is loaded in a notebook: ****)
(*************************************************************************)

(** Local package loaded announcement: **)
packageAnnouncement = "Hamed.m Package (2018.07.11-v2) Loaded ... \n\n";
packageUsage = "Examples:\n" <> Map[StringJoin["\[RightTriangle] ", #1, "\n"]&, GetPackageUsage[]];
PrintNotebookNotification @@ {packageAnnouncementStrs <> packageUsage}

EndPackage[]

(*************************************************************************)
(*************************************************************************)
(*************************************************************************)
(*************************************************************************)
