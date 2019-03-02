
(*************************************************************************)
(*************************************************************************)
(**** : Hamed.m:                                      ********************)
(**** : Author:  Maxie D. Schmidt (maxieds@gmail.com) ********************)
(**** : Created: 07.11.2018                           ********************) 
(*************************************************************************)
(*************************************************************************)

(**** : Clear out any old definitions and stored data on reload: ****) 
ClearAll["PartitionsIntoParts`*"];

LocalPartitionsPackageName = "PartitionsIntoParts`";
BeginPackage[LocalPartitionsPackageName, {"Combinatorica`"}]

StartMemoryInUse = Once[MemoryInUse[]];

(***************************************************************************)
(**** : Package configuration for the core functionality offered here : ****) 
(***************************************************************************)

COMPMETHODSETLENGTHUsage = "Named constant to denote that the partition function HP should be computed by set length " <> 
                              "using Partitions and the Combinatorica package. This method may not scale well " <> 
                              "when n is large due to the amount of memory required to store a large array of lists.\n" <>
                              "See also COMPMETHODGENFUNC, COMPMETHODADAPTIVE, and SetComputationMethod[...].";
COMPMETHODSETLENGTH = 0;

COMPMETHODGENFUNCUsage = "Named constant to denote that the partition function HP should be computed by taking " <> 
                            "coefficients of the truncated (to order n) reciprocal product generating function for " <> 
                            "HP_{p,a}(n). This method should in principle be more efficient than the related " <> 
                            "COMPMETHOD_SETLENGTH method in the package..\n" <>
                            "See also COMPMETHODSETLENGTH, COMPMETHODADAPTIVE, and SetComputationMethod[...].";
COMPMETHODGENFUNC = 1;

COMPMETHODADAPTIVEUsage = "A hybrid of the two methods COMPMETHOD_SETLENGTH and COMPMETHOD_GENFUNC which " <> 
                             "prefers the latter method when n is large-ish, say (by informal heuristic only) " <> 
                             "n \[RightTriangleEqual] 25 to be safe. This parameter can be fine tuned later.\n" <> 
                             "See also SetComputationMethod[...] and the other COMPMETHOD_* named constants above. " <> 
                             "To set the threshold of n values between which we distinguish the particular " <> 
                             "computation method to use, see the local package variable AdaptiveMethodThreshold.";
COMPMETHODADAPTIVE = 4;

PartitionFunctionComputationMethod = COMPMETHODADAPTIVE;
AdaptiveMethodThreshold = 25;

SetComputationMethod::usage = "Default settings are SetComputationMethod[COMPMETHOD_ADAPTIVE, 25]. See the " <> 
                              "COMPMETHOD* named constants usage below for additional possibilities.\n\n" <> 
                              "COMPMETHODSETLENGTH: " <> ToString[COMPMETHODSETLENGTHUsage] <> "\n\n" <> 
                              "COMPMETHODGENFUNC: " <> ToString[COMPMETHODGENFUNCUsage] <> "\n\n" <> 
                              "COMPMETHODADAPTIVE: " <> ToString[COMPMETHODADAPTIVEUsage] <> "\n\n";
SetComputationMethod[compMethod_, adaptiveThreshold_:25] := Block[{}, 
     If[!MemberQ[{COMPMETHODSETLENGTH, COMPMETHODGENFUNC, COMPMETHODADAPTIVE}, compMethod] || adaptiveThreshold <= 0, 
          PrintError["Invalid parameters. See ?SetComputationMethod for usage instructions."];
          Return[-1];
     ];
     PartitionFunctionComputationMethod = compMethod;
     AdaptiveMethodThreshold = adaptiveThreshold;
     Update[HP];
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
HP::usage = "Hamed's special partition function: \n" <> 
            "Computes the number of partitions of n of the form pt+a for fixed primes p and " <> 
            "0 \[LeftTriangleEqual] a < p. The computation method of these integer counts can be changed " <> 
            "to affect the performance speed of the calculations using the function " <> 
            "SetComputationMethod[COMPMETHODSETLENGTH|COMPMETHODGENFUNC|COMPMETHODADAPTIVE].\n" <> 
            "See also HPByPartitions and HPByGF.";
HP[p_, a_, n_] := HP[p, a, n] = Module[{hpValue}, 
     hpValue = 
     Which[PartitionFunctionComputationMethod == COMPMETHODSETLENGTH, HPByPartitions[p, a, n], 
           PartitionFunctionComputationMethod == COMPMETHODGENFUNC, HPByGF[p, a, n], 
           PartitionFunctionComputationMethod == COMPMETHODADAPTIVE 
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
     Length[HamedPartitions[p, a, n]];

HPByGF::usage = "Same as HP[p, a, n] computed by means of a truncated reciprocal product generating function. ";
HPByGF[p_, a_, n_] := HPByGF2[p, a, n] = 
     SeriesCoefficient[1 / Product[1 - Power[q, p * t + a], {t, If[a != 0, 0, 1], n + 1}], {q, 0, n}]

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
                         "(1) SetComputationMethod[COMPMETHODSETLENGTH]\n" <>
                         "(** OR ALTERNATELY: **)\n" <>
                         "(2) SetComputationMethod[COMPMETHODADAPTIVE, Infinity].";
HamedPartitions[p_, a_, n_] := HamedPartitions[p, a, n] = 
     Module[{otf, ofTheFormFunc, partsOfTheForm, partitionIsGoodQ}, 
     otf = {#, IntegerQ[(# - a) / p]}&; 
     ofTheFormFunc = Map[otf, #]&; 
     partsOfTheForm = Map[ofTheFormFunc, Partitions[n]];
     partitionIsGoodQ[part_] := Union[Map[#1[[2]]&, part]] === {True};
     goodPartsList = Select[Map[{Map[First, #1], partitionIsGoodQ[#1]}&, partsOfTheForm], #[[2]] === True&];
     Return[Map[First, goodPartsList]];
];

HPSymmetric::usage = "Hamed's special partition function: \n" <> 
            "Computes the number of partitions of n of the form pt+/-a for fixed primes p and " <> 
            "0 \[LeftTriangleEqual] a < p. The computation method of these integer counts can be changed " <> 
            "to affect the performance speed of the calculations using the function " <> 
            "SetComputationMethod[COMPMETHODSETLENGTH|COMPMETHODGENFUNC|COMPMETHODADAPTIVE].\n" <> 
            "See also HPByPartitions and HPByGF.";
HPSymmetric[p_, a_, n_] := HPSymmetric[p, a, n] = Module[{hpValue}, 
     hpValue = 
     Which[PartitionFunctionComputationMethod == COMPMETHODSETLENGTH, HPSymmetricByPartitions[p, a, n], 
           PartitionFunctionComputationMethod == COMPMETHODGENFUNC, HPSymmetricByGF[p, a, n], 
           PartitionFunctionComputationMethod == COMPMETHODADAPTIVE 
           && n >= AdaptiveMethodThreshold, HPSymmetricByGF[p, a, n], 
           True, HPSymmetricByPartitions[p, a, n]];
     Return[hpValue];
];

HPSymmetricByPartitions[p_, a_, n_] := HPSymmetricByPartitions[p, a, n] = 
     Length[HamedPartitionsSymmetric[p, a, n]];

HPSymmetricByGF[p_, a_, n_] := HPSymmetricByGF2[p, a, n] = 
     Module[{gfIndices, gfProduct, absA},
     absA = Abs[a]; 
     gfIndices = Union[Flatten[Table[{p * t + a, p * t - a}, {t, 0, n + 1}]]];
     gfIndices = Select[gfIndices, (!SameQ[#, 0] && # > 0)&];
     gfProduct = Times @@ Map[(1 - Power[q, #1])&, gfIndices];
     SeriesCoefficient[1 / gfProduct, {q, 0, n}]
];

HamedPartitionsSymmetric::usage = "Symmetric version of HamedPartitions[...].";
HamedPartitionsSymmetric[p_, a_, n_] := HamedPartitionsSymmetric[p, a, n] = 
     Module[{plusAParts, minusAParts, absA},
     absA = Abs[a]; 
     plusAParts = HamedPartitions[p, absA, n];
     minusAParts = HamedPartitions[p, -1 * absA, n];
     Return[Union[plusAParts, minusAParts]];
];

UnitTestingFeaturesUsage = "Used for internal testing, verification, and what we will loosely call our " <> 
                           "first-order approximation to Unit Tests in Mathematica. " <> 
                           "Probably wise not to call this function directly unless you know what you are " <> 
                           "doing and why.\n" <> 
                           "Consider running the public wrapper function RunUnitTests[...] for a sanity check " <> 
                           "that this package is minimally computing what we expect it to compute instead.";

CheckHPDiff::usage = UnitTestingFeaturesUsage;
CheckHPDiff[p_, a_, n_] := HPByPartitions[p, a, n] == HPByGF[p, a, n]

CheckHPDiff2::usage = UnitTestingFeaturesUsage;
CheckHPDiff2[p_, n_] := HP[1, 0, n] == PartitionsP[n]

CheckHPDiff3::usage = UnitTestingFeaturesUsage;
CheckHPDiff3[p_, n_] := HP[2, 1, n] == PartitionsQ[n]

CheckTableValid[tableInput_] := SameQ[Union[Flatten[tableInput]], {True}];

TestHPDiff::usage = UnitTestingFeaturesUsage;
TestHPDiff[numPrimes_:15, nupper_:50] := 
     CheckTableValid[Table[CheckHPDiff[p, a, n], {p, Table[Prime[m], {m, 1, numPrimes}]}, {a, 0, p - 1}, {n, 1, nupper}]];
     
TestHPDiff2::usage = UnitTestingFeaturesUsage;
TestHPDiff2[numPrimes_:15, nupper_:50] := 
     CheckTableValid[Table[CheckHPDiff2[p, n], {p, Table[Prime[m], {m, 1, numPrimes}]}, {n, 1, nupper}]];
     
TestHPDiff3::usage = UnitTestingFeaturesUsage;
TestHPDiff3[numPrimes_:15, nupper_:50] := 
     CheckTableValid[Table[CheckHPDiff3[p, n], {p, Table[Prime[m], {m, 1, numPrimes}]}, {n, 1, nupper}]];

RunUnitTests::usage = "A basic internal package sanity check on the correctness of the core functions " <> 
                      "which we have coded / implemented above. This is the pretty-print notebook wrapper " <> 
                      "around all of the numbered functions in the `Testing` subpackage. This is the " <> 
                      "high-level function you want to call to verify that the package is configured " <> 
                      "properly and generating at least minimally correct results.";
RunUnitTests[numPrimes_:15, nupper_:50] := 
Module[{checkOrXFunc, getUnitTestPassString, test1Result, test2Result, test3Result, testResultsInit, bulletPoints}, 
     checkOrXFunc[boolean_] := Which[boolean, "\[Checkmark] (PASSED)", !boolean, "\[ScriptX] (FAILED)"];
     getUnitTestPassString[timingData_] := ToString[StringForm["`1` and run in `2` seconds\n", 
                                                               checkOrXFunc[Last[timingData]], First[timingData]]];
     test1Result = "[Test 1] Verifying HP[p, a, n] values match across all computation methods: " <> 
                   getUnitTestPassString[Timing[TestHPDiff[numPrimes, nupper]]];
     test2Result = "[Test 2] Verifying HP[p, a, n] values against known functions p(n): " <> 
                   getUnitTestPassString[Timing[TestHPDiff2[numPrimes, nupper]]];
     test3Result = "[Test 3] Verifying HP[p, a, n] values against known functions q(n): " <> 
                   getUnitTestPassString[Timing[TestHPDiff3[numPrimes, nupper]]];
     testResultsInit = {test1Result, test2Result, test3Result};
     bulletPoints = {"\[ClubSuit]", "\[SpadeSuit]", "\[BlackQueen]"};
     testResultsDisplayString = StringJoin @@ MapIndexed[StringJoin[bulletPoints[[ First[#2] ]], #1]&, testResultsInit];
     PrintNotebookNotification[Orange][testResultsDisplayString];
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
     If[omegaOnes == 0, ellLength, muCount - omegaOnes]
]

GetPartitionIndexParameter::usage = "The parameter t if the component represents pt+a.";
GetPartitionIndexParameter[p_, a_, comp_] := (comp - a) / p;

GetSinglePartitionStats[p_, a_, part_, includeFormatting_:True] := 
Module[{partitionStats, partSummarySpec, outputDesc, getPartitionIndexParameters},
     getPartitionIndexParameters := StringJoin @@ Map[ToString[StringForm["`1`\[FilledSmallCircle]`2`+`3`; ", 
                                                                                 p, GetPartitionIndexParameter[p, a, #1], a]]&, part];
     partitionStats = {getPartitionIndexParameters, Plus @@ #, 
                       PartitionRank[#], PartitionCrank[#], Min[#], Max[#], 
                       PartitionOnes[#], PartitionMu[#]}&[part];
     partSummarySpec = "\[Lambda] = `1` \[RightGuillemet] `2` ;;;" <> 
                       "rank=`3`, crank=`4`, spart=`5`, lpart=`6`, ones=`7`, \[Mu](\[Lambda])=`8`";
     outputDesc = ToString[StringForm[partSummarySpec, ##]]& @@ partitionStats;
     If[includeFormatting, 
          outputDesc = " \[FivePointedStar] " <> outputDesc <> "\n";
     ];
     Return[ToString[outputDesc]];
];

PrintPartitionStats[p_, a_, n_] := Module[{headerString, summaryString}, 
     headerString = ToString[StringForm["Partitions of n = `1` into parts of the form `2`t+`3`:\n\n", n, p, a]];
     headerString = headerString <> 
                    ToString[StringForm["HP(p=`1`, a=`2`; n=`3`) = `4`\n", p, a, n, HPByPartitions[p, a, n]]];
     headerString = headerString <> 
                    ToString[StringForm["For reference: p(`1`)=`2`, q(`1`)=`3`.\n", n, PartitionsP[n], PartitionsQ[n]]];
     headerString = headerString <> "\n==========================================================================\n\n";
     summaryString = StringJoin @@ Map[GetSinglePartitionStats[p, a, #1]&, HamedPartitions[p, a, n]];
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
GetPartitionRankSumModuloBases[p_, a_, n_] := With[{hp = GetPartitionRankSum[p, a, n]}, 
     With[{congList = Map[First, Select[Table[{pv, Mod[hp, pv, 0] == 0}, {pv, Table[Prime[m], {m, 1, hp}]}], #[[1]] <= Abs[hp] && #[[2]]&]]}, 
          If[congList === {}, "None", congList]
     ]
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
     ttIneqLHS = Power[p1, 2]; 
     ttIneqRHS = p0 * p2;
     isConvex = ttIneqLHS >= ttIneqRHS;
     isEqual = ttIneqLHS == ttIneqRHS;
     ineqDirStr = ToString[
                       StringForm["`1` `2` `3`", ttIneqLHS, 
                       Which[isEqual, "\[Congruent]", 
                             isConvex, "\[RightTriangleEqual]", 
                             True, "\[NotRightTriangleEqual]"], ttIneqRHS]
                  ];
     Return[{isConvex, ineqDirStr}];
];

(*************************************************************************)
(**** : A possible generalization of Stanley's theorem?             : ****)
(**** : For exploring Maxie's new suggestion that it may be         : ****)
(**** : reasonable or directly possible to form a generalized       : ****)
(**** : analog to Stanley's theorem which interprets the count      : ****)
(**** : statistic of the number of ones in all partitions of n      : ****)
(**** : into our specified special forms.                           : ****)
(*************************************************************************)

StanleysTheoremComponents::usage = "Computes (a) the number of all ones in all the partitions returned by " <> 
                                   "HamedPartitions[p, a, n]; and (b) Attempts to correlate this to the " <> 
                                   "second count in Stanley's theorem: namely, the sums of numbers of parts in " <> 
                                   "all such partitions. We seek to determine how (if at all) these components " <> 
                                   "are related in this new setting of (pt+a)-formed partitions of n.";
StanleysTheoremComponents[p_, a_, n_] := StanleysTheoremComponents[p, a, n] = 
Module[{hpParts, numberOfOnes, countsOfParts, secondCountStatistic}, 
     hpParts = HamedPartitions[p, a, n];
     numberOfOnes = Plus @@ Map[PartitionOnes, hpParts];
     countsOfParts = Map[Length, hpParts];
     secondCountStatistic = Plus @@ countsOfParts;
     Return[{numberOfOnes, secondCountStatistic, countsOfParts}];
];

StanleysTheoremQ::usage = "Boolean-valued indicator of whether Stanley's theorem (applied directly) is also " <> 
                          "true for these parameter values of (p, a; n).\n" <> 
                          "See also StanleysTheoremComponents[p, a, n].";
StanleysTheoremQ[p_, a_, n_] := With[{stcomps = StanleysTheoremComponents[p, a, n]}, 
     stcomps[[1]] == stcomps[[2]]
];

(*************************************************************************)
(**** : Text processing, display and general utility functions      : ****) 
(*************************************************************************)

GetDatestamp[] := ToString[StringForm["`1`/`2`/`3` : `4`:`5`:`6`", ##]] @@ 
     Map[DateValue, {"Month", "Day", "Year", "Hour", "Minute", "Second"}]

PrintNotebookNotificationDetailed[textColor_, backgroundColor_, borderColor_] := 
Function[cellText, 
     If[$Notebooks,
     CellPrint[Cell[cellText, "Print", 
                    FontColor -> textColor, 
                    FontSize -> 14, 
                    CellFrame -> 2.5, 
                    CellFrameColor -> borderColor, 
                    Background -> backgroundColor, 
                    CellFrameMargins -> 14, 
                    ShowCellBracket -> False
                   ]
              ],
     Print[cellText]; 
    ]
];

PrintNotebookNotification[baseColor_:Green] := 
     PrintNotebookNotificationDetailed[Darker[baseColor, 0.8], Lighter[baseColor, 0.5], Black];

PrintError[errorText_] := 
     PrintNotebookNotification[Red][StringJoin["\[ScriptX]\n", " !!! ERROR : >>> @ ", 
                                               GetDatestamp[], "\n", errorText, "\n\[ScriptX]"]]

(*************************************************************************)
(**** : Package help and information printing routines              : ****) 
(*************************************************************************) 

PackageSingleStringFromList[headerStr_, listComps_, bulletMarker_:"\[RightTriangle]"] := 
Module[{constructListElementFunc}, 
     constructListElementFunc = (" " <> bulletMarker <> " " <> ToString[#] <> "\n")&;
     Return[headerStr <> "\n" <> StringJoin @@ Map[constructListElementFunc, listComps]];
];

Options[GetHPPackageUsageString] = {BulletPointMarker -> "\[RightTriangle]"};
GetHPPackageUsageString[OptionsPattern[]] := Module[{lineFunc, usageDescList}, 
     lineFunc = ToString[ToString[" "] <> ToString[OptionValue[BulletPointMarker]] <> " " <> ToString[#]]&;
     usageDescList = {"Loading the Package: \n", 
                      lineFunc["SetDirectory[\"/path/to/saved/HamedDotMFile\"]\n"], 
                      lineFunc["<< Hamed.m\n\n"], 
                      "Core Functions Provided by the Package: \n", 
                      lineFunc["HP[p, a, n], HamedPartitions[p, a, n];\n"], 
                      lineFunc["HPSymmetric[p, a, n], HamedPartitionsSymmetric[p, a, n];\n"], 
                      lineFunc["PrintPartitionStats[p, a, n]; RunUnitTests[] (if you so please); etc.;\n\n"], 
                      "Helpful References: \n", 
                      lineFunc["See HPPackageHelp[], HPPackageExamples[], and the package maintainer's website at " <> 
                      "https://github.com/maxieds/PartitionsIntoParts for more documentation."]
     };
     Return[StringJoin @@ usageDescList];
];

Options[HPPackageUsage] = {NotebookDisplayColor -> RGBColor[1.0, 0.0, 0.54]};
HPPackageUsage[OptionsPattern[]] := 
     PrintNotebookNotification[OptionValue[NotebookDisplayColor]][GetHPPackageUsageString[]];

Options[HPPackageHelp] = {NotebookDisplayColor -> Yellow};
HPPackageHelp[OptionsPattern[]] := Module[{helpHeader, helpStr}, 
     helpHeader = "Helful Information About the Package: ";
     helpStr = "You are fortunate in using this package in so much as its authors have decorated the " <> 
               "core functions you will be using with Fn::usage package strings. This means that " <> 
               "helpful documentation for a particular constant or package feature can be accessed by " <> 
               "running ?ConstantOrFunctionName in your working notebook. If you are unsure where to " <> 
               "start, or just want a detailed list of what exactly we have implemented here, you can " <> 
               "evaluate the following in your notebook with the package loaded " <> 
               "?" <> LocalPartitionsPackageName <> "`* ... " <> 
               "Alternatively, we have added the separate helper functions " <> 
               "HPPackageUsage[] and HPPackageExamples[] for you to use and explore. Enjoy!";
    PrintNotebookNotification[OptionValue[NotebookDisplayColor]][helpHeader <> "\n\n" <> helpStr];
];             

Options[HPPackageExamples] = {NotebookDisplayColor -> Cyan};
HPPackageExamples[OptionsPattern[]] := 
     PrintNotebookNotification[OptionValue[NotebookDisplayColor]][
     PackageSingleStringFromList["Examples of the Package in Use: ", {
     "Please see the sample notebook at " <> 
     "https://github.com/maxieds/PartitionsIntoParts/blob/master/hameds-partition-function.nb."
     }]];

(*************************************************************************)
(**** : Perform pretty printing of the package details and          : ****) 
(**** : revsision information if the package is loaded in a notebook: ****)
(*************************************************************************)

(** Local package loaded announcement: **)
packageAnnouncement = ToString[LocalPartitionsPackageName] <> " Package (2018.07.11-v2) Loaded!\n\n";
PrintNotebookNotification[][packageAnnouncement <> GetHPPackageUsageString[]]

EndPackage[]

(*************************************************************************)
(*************************************************************************)
(*************************************************************************)
(*************************************************************************)
