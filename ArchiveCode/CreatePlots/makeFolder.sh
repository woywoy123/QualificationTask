#!/bin/bash 

#rm -r Figures
mkdir Figures 
mkdir Figures/Appendix
mkdir Figures/Debug 
mkdir Figures/Estimation 
mkdir Figures/FLost 
mkdir Figures/Introduction 
mkdir Figures/Truth

mkdir Figures/Appendix/Debugging 
mkdir Figures/Appendix/Debugging/Ratio
mkdir Figures/Appendix/Debugging/Ratio/Nominal
mkdir Figures/Appendix/Debugging/Ratio/Nominal/IBL
mkdir Figures/Appendix/Debugging/Ratio/Nominal/layer1
mkdir Figures/Appendix/Debugging/Ratio/Nominal/layer2

mkdir Figures/Appendix/Debugging/Ratio/Smooth
mkdir Figures/Appendix/Debugging/Ratio/Smooth/IBL
mkdir Figures/Appendix/Debugging/Ratio/Smooth/layer1
mkdir Figures/Appendix/Debugging/Ratio/Smooth/layer2

mkdir Figures/Appendix/Debugging/Ratio/SmoothSubtract
mkdir Figures/Appendix/Debugging/Ratio/SmoothSubtract/IBL
mkdir Figures/Appendix/Debugging/Ratio/SmoothSubtract/layer1
mkdir Figures/Appendix/Debugging/Ratio/SmoothSubtract/layer2

mkdir Figures/Appendix/Debugging/Ratio/Subtract
mkdir Figures/Appendix/Debugging/Ratio/Subtract/IBL
mkdir Figures/Appendix/Debugging/Ratio/Subtract/layer1
mkdir Figures/Appendix/Debugging/Ratio/Subtract/layer2

mkdir Figures/Debug/Ratio
mkdir Figures/Debug/Ratio/Blayer
mkdir Figures/Debug/Ratio/Blayer_Smooth
mkdir Figures/Debug/Ratio/Blayer_SmoothSubtract
mkdir Figures/Debug/Ratio/Blayer_Subtract

mkdir Figures/Debug/IllustrationOfCases

DBG_AP=Figures/Appendix/Debugging/Ratio
DBG_MAIN=Figures/Debug/Ratio
DBG_NOM=Output_Debug/_Debug
DBG_SMO=Output_Debug/_Debug_Smooth
DBG_SUB=Output_Debug/_Debug_Subtract
DBG_SMSB=Output_Debug/_Debug_Subtract_Smooth


#Appendix Debug Ratio 
## Nominal 
cp -r $DBG_NOM/RatioMatrix/IBL/* $DBG_AP/Nominal/IBL
cp -r $DBG_NOM/RatioMatrix/layer1/* $DBG_AP/Nominal/layer1
cp -r $DBG_NOM/RatioMatrix/layer2/* $DBG_AP/Nominal/layer2

## Smooth
cp -r $DBG_SMO/RatioMatrix/IBL/* $DBG_AP/Smooth/IBL
cp -r $DBG_SMO/RatioMatrix/layer1/* $DBG_AP/Smooth/layer1
cp -r $DBG_SMO/RatioMatrix/layer2/* $DBG_AP/Smooth/layer2

## Subtract 
cp -r $DBG_SUB/RatioMatrix/IBL/* $DBG_AP/Subtract/IBL
cp -r $DBG_SUB/RatioMatrix/layer1/* $DBG_AP/Subtract/layer1
cp -r $DBG_SUB/RatioMatrix/layer2/* $DBG_AP/Subtract/layer2

## Smooth + Subtract
cp -r $DBG_SMSB/RatioMatrix/IBL/* $DBG_AP/SmoothSubtract/IBL
cp -r $DBG_SMSB/RatioMatrix/layer1/* $DBG_AP/SmoothSubtract/layer1
cp -r $DBG_SMSB/RatioMatrix/layer2/* $DBG_AP/SmoothSubtract/layer2

#Main Debug Ratio
cp -r $DBG_NOM/RatioMatrix/Blayer/* $DBG_MAIN/Blayer
cp -r $DBG_SMO/RatioMatrix/Blayer/* $DBG_MAIN/Blayer_Smooth
cp -r $DBG_SUB/RatioMatrix/Blayer/* $DBG_MAIN/Blayer_Subtract
cp -r $DBG_SMSB/RatioMatrix/Blayer/* $DBG_MAIN/Blayer_SmoothSubtract

cp $DBG_NOM/IllustrationOfCases/Case1/NormalCase1.png Figures/Debug/IllustrationOfCases/NormalCase1.png
cp $DBG_NOM/IllustrationOfCases/Case2/ShiftNormalCase2.png Figures/Debug/IllustrationOfCases/ShiftNormalCase2.png
cp $DBG_NOM/IllustrationOfCases/Case3/ShiftNormalWidthFFTCase3.png Figures/Debug/IllustrationOfCases/ShiftNormalWidthFFTCase3.png

# Estimation 
EST_FitTo=Output/_FitTo/ShapePerformance
EST_MIN=Output/_Minimizer/ShapePerformance


EST_ROOT=Figures/Estimation
mkdir $EST_ROOT/Nominal_Minimizer_Data
mkdir $EST_ROOT/Nominal_FitTo_DirectTruth
mkdir $EST_ROOT/Nominal_FitTo_Data

cp $EST_MIN/trk1/Template-trk1_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_Minimizer_Data
cp $EST_MIN/trk2/Template-trk2_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_Minimizer_Data
cp $EST_MIN/trk3/Template-trk3_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_Minimizer_Data
cp $EST_MIN/trk4/Template-trk4_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_Minimizer_Data

cp $EST_FitTo/trk1/Template-trk1_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_Data
cp $EST_FitTo/trk2/Template-trk2_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_Data
cp $EST_FitTo/trk3/Template-trk3_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_Data
cp $EST_FitTo/trk4/Template-trk4_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_Data


mkdir $EST_ROOT/Nominal_FitTo_DirectTruth/Track-1
mkdir $EST_ROOT/Nominal_FitTo_DirectTruth/Track-2
mkdir $EST_ROOT/Nominal_FitTo_DirectTruth/Track-3
mkdir $EST_ROOT/Nominal_FitTo_DirectTruth/Track-4

cp $EST_FitTo/trk1/Test-trk1_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-1
cp $EST_FitTo/trk2/Test-trk2_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-2
cp $EST_FitTo/trk3/Test-trk3_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-3
cp $EST_FitTo/trk4/Test-trk4_Blayer_Shape_ERROR.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-4

cp Output/_FitTo/Histograms/Test/trk1/Normal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-1/Blayer_800_1000_GeV_Normal.png
cp Output/_FitTo/Histograms/Test/trk2/Normal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-2/Blayer_800_1000_GeV_Normal.png
cp Output/_FitTo/Histograms/Test/trk3/Normal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-3/Blayer_800_1000_GeV_Normal.png
cp Output/_FitTo/Histograms/Test/trk4/Normal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-4/Blayer_800_1000_GeV_Normal.png

cp Output/_FitTo/Histograms/Test/trk1/ShiftNormal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-1/Blayer_800_1000_GeV_ShiftNormal.png
cp Output/_FitTo/Histograms/Test/trk2/ShiftNormal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-2/Blayer_800_1000_GeV_ShiftNormal.png
cp Output/_FitTo/Histograms/Test/trk3/ShiftNormal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-3/Blayer_800_1000_GeV_ShiftNormal.png
cp Output/_FitTo/Histograms/Test/trk4/ShiftNormal/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-4/Blayer_800_1000_GeV_ShiftNormal.png

cp Output/_FitTo/Histograms/Test/trk1/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-1/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png
cp Output/_FitTo/Histograms/Test/trk2/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-2/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png
cp Output/_FitTo/Histograms/Test/trk3/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-3/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png
cp Output/_FitTo/Histograms/Test/trk4/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $EST_ROOT/Nominal_FitTo_DirectTruth/Track-4/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png

mkdir Figures/Appendix/ErrorDirectTruthFits_FitTo
mkdir Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-1
mkdir Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-2
mkdir Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-3
mkdir Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-4

cp $EST_FitTo/trk1/Test-trk1_IBL_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-1
cp $EST_FitTo/trk2/Test-trk2_IBL_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-2
cp $EST_FitTo/trk3/Test-trk3_IBL_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-3
cp $EST_FitTo/trk4/Test-trk4_IBL_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-4

cp $EST_FitTo/trk1/Test-trk1_layer1_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-1
cp $EST_FitTo/trk2/Test-trk2_layer1_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-2
cp $EST_FitTo/trk3/Test-trk3_layer1_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-3
cp $EST_FitTo/trk4/Test-trk4_layer1_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-4

cp $EST_FitTo/trk1/Test-trk1_layer2_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-1
cp $EST_FitTo/trk2/Test-trk2_layer2_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-2
cp $EST_FitTo/trk3/Test-trk3_layer2_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-3
cp $EST_FitTo/trk4/Test-trk4_layer2_Shape_ERROR.png Figures/Appendix/ErrorDirectTruthFits_FitTo/Track-4


mkdir Figures/Appendix/ErrorOfDataTemplateFits
mkdir Figures/Appendix/ErrorOfDataTemplateFits/FitTo
mkdir Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-1
mkdir Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-2
mkdir Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-3
mkdir Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-4

cp $EST_FitTo/trk1/Template-trk1_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-1
cp $EST_FitTo/trk2/Template-trk2_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-2
cp $EST_FitTo/trk3/Template-trk3_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-3
cp $EST_FitTo/trk4/Template-trk4_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-4

cp $EST_FitTo/trk1/Template-trk1_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-1
cp $EST_FitTo/trk2/Template-trk2_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-2
cp $EST_FitTo/trk3/Template-trk3_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-3
cp $EST_FitTo/trk4/Template-trk4_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-4

cp $EST_FitTo/trk1/Template-trk1_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-1
cp $EST_FitTo/trk2/Template-trk2_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-2
cp $EST_FitTo/trk3/Template-trk3_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-3
cp $EST_FitTo/trk4/Template-trk4_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/FitTo/Track-4


mkdir Figures/Appendix/ErrorOfDataTemplateFits/Minimizer
mkdir Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-1
mkdir Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-2
mkdir Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-3
mkdir Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-4

cp $EST_MIN/trk1/Template-trk1_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-1
cp $EST_MIN/trk2/Template-trk2_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-2
cp $EST_MIN/trk3/Template-trk3_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-3
cp $EST_MIN/trk4/Template-trk4_IBL_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-4

cp $EST_MIN/trk1/Template-trk1_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-1
cp $EST_MIN/trk2/Template-trk2_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-2
cp $EST_MIN/trk3/Template-trk3_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-3
cp $EST_MIN/trk4/Template-trk4_layer1_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-4

cp $EST_MIN/trk1/Template-trk1_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-1
cp $EST_MIN/trk2/Template-trk2_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-2
cp $EST_MIN/trk3/Template-trk3_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-3
cp $EST_MIN/trk4/Template-trk4_layer2_Shape_ERROR.png Figures/Appendix/ErrorOfDataTemplateFits/Minimizer/Track-4

mkdir Figures/Appendix/ExampleOfExperimental
cp Output/_Minimizer/Histograms/Template/trk1/Experimental/Blayer_800_1000_GeV.png Figures/Appendix/ExampleOfExperimental/Track-1_Blayer_800_1000_GeV.png
cp Output/_Minimizer/Histograms/Template/trk2/Experimental/Blayer_800_1000_GeV.png Figures/Appendix/ExampleOfExperimental/Track-2_Blayer_800_1000_GeV.png
cp Output/_Minimizer/Histograms/Template/trk3/Experimental/Blayer_800_1000_GeV.png Figures/Appendix/ExampleOfExperimental/Track-3_Blayer_800_1000_GeV.png
cp Output/_Minimizer/Histograms/Template/trk4/Experimental/Blayer_800_1000_GeV.png Figures/Appendix/ExampleOfExperimental/Track-4_Blayer_800_1000_GeV.png

mkdir Figures/Appendix/FLost
mkdir Figures/Appendix/FLost/AlgoCombination_FLost2
mkdir Figures/Appendix/FLost/AlgoCombination_FLost2/IBL
mkdir Figures/Appendix/FLost/AlgoCombination_FLost2/layer1
mkdir Figures/Appendix/FLost/AlgoCombination_FLost2/layer2
mkdir Figures/FLost/AlgoCombination_FLost2

cp -r Final/BestAlgos/Minimizer/FLost2_IBL_* Figures/Appendix/FLost/AlgoCombination_FLost2/IBL
cp -r Final/BestAlgos/Minimizer/FLost2_layer1_* Figures/Appendix/FLost/AlgoCombination_FLost2/layer1
cp -r Final/BestAlgos/Minimizer/FLost2_layer2_* Figures/Appendix/FLost/AlgoCombination_FLost2/layer2
cp -r Final/BestAlgos/Minimizer/FLost2_Blayer_* Figures/FLost/AlgoCombination_FLost2

mkdir Figures/Appendix/FLost/AlgoCombination_FLost3
mkdir Figures/Appendix/FLost/AlgoCombination_FLost3/IBL
mkdir Figures/Appendix/FLost/AlgoCombination_FLost3/layer1
mkdir Figures/Appendix/FLost/AlgoCombination_FLost3/layer2
mkdir Figures/FLost/AlgoCombination_FLost3

cp -r Final/BestAlgos/Minimizer/FLost3_IBL_* Figures/Appendix/FLost/AlgoCombination_FLost3/IBL
cp -r Final/BestAlgos/Minimizer/FLost3_layer1_* Figures/Appendix/FLost/AlgoCombination_FLost3/layer1
cp -r Final/BestAlgos/Minimizer/FLost3_layer2_* Figures/Appendix/FLost/AlgoCombination_FLost3/layer2
cp -r Final/BestAlgos/Minimizer/FLost3_Blayer_* Figures/FLost/AlgoCombination_FLost3

mkdir Figures/Appendix/FLost/AlgoWithBestMini_FLost2
mkdir Figures/FLost/AlgoWithBestMini_FLost2

cp Final/BestAlgos/OptimalMiniWithAlgo/FLost2_IBL_AlgoWithBestMin.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost2_layer1_AlgoWithBestMin.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost2_layer2_AlgoWithBestMin.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost2_Blayer_AlgoWithBestMin.png Figures/FLost/AlgoWithBestMini_FLost2
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost2_Verdict.txt Figures/FLost/AlgoWithBestMini_FLost2

mkdir Figures/Appendix/FLost/AlgoWithBestMini_FLost3
mkdir Figures/FLost/AlgoWithBestMini_FLost3

cp Final/BestAlgos/OptimalMiniWithAlgo/FLost3_IBL_AlgoWithBestMin.png Figures/Appendix/FLost/AlgoWithBestMini_FLost3
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost3_layer1_AlgoWithBestMin.png Figures/Appendix/FLost/AlgoWithBestMini_FLost3
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost3_layer2_AlgoWithBestMin.png Figures/Appendix/FLost/AlgoWithBestMini_FLost3
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost3_Blayer_AlgoWithBestMin.png Figures/FLost/AlgoWithBestMini_FLost3
cp Final/BestAlgos/OptimalMiniWithAlgo/FLost3_Verdict.txt Figures/FLost/AlgoWithBestMini_FLost3

mkdir Figures/Appendix/FLost/BestCombinations_FLost2
mkdir Figures/FLost/BestCombinations_FLost2
cp Final/BestAlgos/FLost2_Error/FLost2_Blayer_Best.png Figures/FLost/BestCombinations_FLost2
cp Final/BestAlgos/FLost2_Error/FLost2_ScoreMatrix.txt Figures/FLost/BestCombinations_FLost2

cp Final/BestAlgos/FLost2_Error/FLost2_layer2_Best.png Figures/Appendix/FLost/BestCombinations_FLost2
cp Final/BestAlgos/FLost2_Error/FLost2_layer1_Best.png Figures/Appendix/FLost/BestCombinations_FLost2
cp Final/BestAlgos/FLost2_Error/FLost2_IBL_Best.png Figures/Appendix/FLost/BestCombinations_FLost2


mkdir Figures/Appendix/FLost/BestCombinations_FLost3
mkdir Figures/FLost/BestCombinations_FLost3
cp Final/BestAlgos/FLost3_Error/FLost3_Blayer_Best.png Figures/FLost/BestCombinations_FLost3
cp Final/BestAlgos/FLost3_Error/FLost3_ScoreMatrix.txt Figures/FLost/BestCombinations_FLost3

cp Final/BestAlgos/FLost3_Error/FLost3_layer2_Best.png Figures/Appendix/FLost/BestCombinations_FLost3
cp Final/BestAlgos/FLost3_Error/FLost3_layer1_Best.png Figures/Appendix/FLost/BestCombinations_FLost3
cp Final/BestAlgos/FLost3_Error/FLost3_IBL_Best.png Figures/Appendix/FLost/BestCombinations_FLost3

mkdir Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets
mkdir Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk1
mkdir Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk2
mkdir Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk3

cp Output/_Minimizer/Histograms/Template/trk1/Experimental/Blayer_800_1000_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk1
cp Output/_Minimizer/Histograms/Template/trk1/Experimental/Blayer_1000_1200_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk1
cp Output/_Minimizer/Histograms/Template/trk1/Experimental/Blayer_1200_1400_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk1
cp Output/_Minimizer/Histograms/Template/trk1/Experimental/Blayer_1400_1600_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk1

cp Output/_Minimizer/Histograms/Template/trk2/Experimental/Blayer_800_1000_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk2
cp Output/_Minimizer/Histograms/Template/trk2/Experimental/Blayer_1000_1200_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk2
cp Output/_Minimizer/Histograms/Template/trk2/Experimental/Blayer_1200_1400_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk2
cp Output/_Minimizer/Histograms/Template/trk2/Experimental/Blayer_1400_1600_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk2

cp Output/_Minimizer/Histograms/Template/trk3/Experimental/Blayer_800_1000_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk3
cp Output/_Minimizer/Histograms/Template/trk3/Experimental/Blayer_1000_1200_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk3
cp Output/_Minimizer/Histograms/Template/trk3/Experimental/Blayer_1200_1400_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk3
cp Output/_Minimizer/Histograms/Template/trk3/Experimental/Blayer_1400_1600_GeV.png Figures/Appendix/FLost/AlgoWithBestMini_FLost2/ExperimentalJets/trk3
