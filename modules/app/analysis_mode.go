package app

import (
	"atsp_aco_msa/modules/analysis/msaHeuristicTours"
	"atsp_aco_msa/modules/analysis/structuralComparison"
	"atsp_aco_msa/modules/analysis/tuningSummary"
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"image/color"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

func runAnalysisMode(atspsData []AtspData, analysisScope string, tuningHeuristics []string, workers int) error {
	if analysisScope == analysisScopeTuning {
		if err := saveTuningSummary(resultsDirectoryName, atspsData, tuningHeuristics); err != nil {
			return err
		}
		fmt.Printf("Tuning summary saved to %s\n", filepath.Join(resultsDirectoryName, tuningSummary.ReportFileName))
		return nil
	}

	if analysisScope == analysisScopeAll {
		if err := saveTuningSummary(resultsDirectoryName, atspsData, tuningHeuristics); err != nil {
			return err
		}
		fmt.Printf("Tuning summary saved to %s\n", filepath.Join(resultsDirectoryName, tuningSummary.ReportFileName))
	}

	gksDeviationReportPath := filepath.Join(finalResultsDirectoryName, "gks_deviation.md")
	gksDeviationAtspData, err := loadSelectedAtspData(instanceSetAllKnown)
	if err != nil {
		return err
	}

	if analysisScope == analysisScopeGksDeviation {
		if err := ensureMsaHeuristicCache(gksDeviationAtspData, workers, false); err != nil {
			return err
		}
		if err := saveGksDeviationReport(gksDeviationReportPath, gksDeviationAtspData, gksDeviationMsaPatchBiases); err != nil {
			return err
		}
		fmt.Printf("GKS deviation report saved to %s using %d all-known instance(s)\n", gksDeviationReportPath, len(gksDeviationAtspData))
		return nil
	}

	if err := ensureMsaHeuristicArtifacts(atspsData, workers, false); err != nil {
		return err
	}

	msaHeuristicTourConfigs := make([]msaHeuristicTours.InstanceConfig, 0, len(atspsData))
	structuralConfigs := make([]structuralComparison.InstanceConfig, 0, len(atspsData))
	for _, atspData := range atspsData {
		msaHeuristicTourConfigs = append(msaHeuristicTourConfigs, msaHeuristicTours.InstanceConfig{
			Name:                                atspData.name,
			Dimension:                           len(atspData.matrix),
			MsaHeuristicDirectoryPath:           atspData.msaHeuristicDirectoryPath,
			OptimalToursCsvPath:                 atspData.optimalUniqueToursCsvPath,
			ToursHeatmapPath:                    atspData.toursHeatmapPlotPath,
			ToursHistogramPath:                  atspData.toursHistogramPlotPath,
			MsaHeuristicToursOverlapHeatmapPath: atspData.msaHeuristicToursOverlapHeatmapPlotPath,
		})

		structuralConfigs = append(structuralConfigs, structuralComparison.InstanceConfig{
			Name:                      atspData.name,
			Dimension:                 len(atspData.matrix),
			Matrix:                    atspData.matrix,
			MsaHeuristicDirectoryPath: atspData.msaHeuristicDirectoryPath,
			OptimalToursCsvPath:       atspData.optimalUniqueToursCsvPath,
		})
	}

	if err := msaHeuristicTours.AnalyzeInstances(msaHeuristicTours.Config{Instances: msaHeuristicTourConfigs}); err != nil {
		return err
	}

	structuralAnalyses, err := structuralComparison.AnalyzeInstances(structuralComparison.Config{
		Instances:     structuralConfigs,
		HighThreshold: 1.0,
		MsaPatchBias:  finalCycleCoverMsaPatchingMsaPatchBias,
	})
	if err != nil {
		return err
	}

	structuralSimilarityReportPath := filepath.Join(finalResultsDirectoryName, "structural_similarity.md")
	if err := saveStructuralSimilarityReport(structuralSimilarityReportPath, structuralAnalyses); err != nil {
		return err
	}

	heuristicOverlapReportPath := filepath.Join(finalResultsDirectoryName, "msa_cycle_cover_overlap.md")
	if err := saveMsaHeuristicCycleCoverOverlapReport(heuristicOverlapReportPath, structuralAnalyses); err != nil {
		return err
	}

	msaCountScalingReportPath := filepath.Join(finalResultsDirectoryName, "msa_count_scaling.md")
	if err := saveMsaCountScalingReport(msaCountScalingReportPath, atspsData); err != nil {
		return err
	}

	finalControlsRootPath := finalControlsResultsRootPath(finalResultsDirectoryName)
	randomSparseControlReportPath := filepath.Join(finalControlsRootPath, "random_sparse_control.md")
	randomSparseControlReportSaved, err := saveRandomSparseControlReport(randomSparseControlReportPath, atspsData, finalResultsDirectoryName, finalControlsRootPath)
	if err != nil {
		return err
	}
	distanceRankedSparseControlReportPath := filepath.Join(finalControlsRootPath, "distance_ranked_sparse_control.md")
	distanceRankedSparseControlReportSaved, err := saveDistanceRankedSparseControlReport(distanceRankedSparseControlReportPath, atspsData, finalResultsDirectoryName, finalControlsRootPath)
	if err != nil {
		return err
	}
	shuffledMsaControlReportPath := filepath.Join(finalControlsRootPath, "shuffled_msa_control.md")
	shuffledMsaControlReportSaved, err := saveShuffledMsaControlReport(shuffledMsaControlReportPath, atspsData, finalResultsDirectoryName, finalControlsRootPath)
	if err != nil {
		return err
	}
	if err := removeFileIfExists(filepath.Join(finalResultsDirectoryName, "random_sparse_control.md")); err != nil {
		return err
	}
	if err := removeFileIfExists(filepath.Join(finalResultsDirectoryName, "distance_ranked_sparse_control.md")); err != nil {
		return err
	}
	if err := removeFileIfExists(filepath.Join(finalResultsDirectoryName, "shuffled_msa_control.md")); err != nil {
		return err
	}

	if err := ensureMsaHeuristicCache(gksDeviationAtspData, workers, false); err != nil {
		return err
	}

	if err := saveGksDeviationReport(gksDeviationReportPath, gksDeviationAtspData, gksDeviationMsaPatchBiases); err != nil {
		return err
	}

	finalResultsSummaryPath, finalRows, finalSummarySaved, err := runFinalResultsAnalysis(atspsData, structuralAnalyses, finalResultsDirectoryName)
	if err != nil {
		return err
	}

	finalThreeOptResultsSummaryPath, finalThreeOptRows, finalThreeOptSummarySaved, err := runFinalResultsAnalysis(atspsData, structuralAnalyses, finalThreeOptResultsDirectoryName)
	if err != nil {
		return err
	}

	fmt.Printf("Structural similarity report saved to %s\n", structuralSimilarityReportPath)
	fmt.Printf("MSA heuristic/cycle-cover overlap report saved to %s\n", heuristicOverlapReportPath)
	fmt.Printf("MSA count scaling report saved to %s\n", msaCountScalingReportPath)
	if randomSparseControlReportSaved {
		fmt.Printf("Random sparse control report saved to %s\n", randomSparseControlReportPath)
	}
	if distanceRankedSparseControlReportSaved {
		fmt.Printf("Distance-ranked sparse control report saved to %s\n", distanceRankedSparseControlReportPath)
	}
	if shuffledMsaControlReportSaved {
		fmt.Printf("Shuffled MSA control report saved to %s\n", shuffledMsaControlReportPath)
	}
	fmt.Printf("GKS deviation report saved to %s using %d all-known instance(s)\n", gksDeviationReportPath, len(gksDeviationAtspData))
	if finalSummarySaved {
		fmt.Printf("Final results summary saved to %s\n", finalResultsSummaryPath)
		fmt.Printf("Pairwise performance report saved to %s\n", filepath.Join(finalResultsDirectoryName, "pairwise_performance.md"))
		fmt.Printf("Convergence summary report saved to %s\n", filepath.Join(finalResultsDirectoryName, "convergence_summary.md"))
		fmt.Printf("Structural/performance link report saved to %s\n", filepath.Join(finalResultsDirectoryName, "structural_performance_link.md"))
	}
	if finalThreeOptSummarySaved {
		fmt.Printf("Final+3opt results summary saved to %s\n", finalThreeOptResultsSummaryPath)
		fmt.Printf("Final+3opt pairwise performance report saved to %s\n", filepath.Join(finalThreeOptResultsDirectoryName, "pairwise_performance.md"))
		fmt.Printf("Final+3opt convergence summary report saved to %s\n", filepath.Join(finalThreeOptResultsDirectoryName, "convergence_summary.md"))
		fmt.Printf("Final+3opt structural/performance link report saved to %s\n", filepath.Join(finalThreeOptResultsDirectoryName, "structural_performance_link.md"))
	}
	if finalSummarySaved && finalThreeOptSummarySaved {
		threeOptComparisonPath := filepath.Join(finalThreeOptResultsDirectoryName, "comparison_to_final.md")
		if err := saveFinalThreeOptComparisonReport(threeOptComparisonPath, finalRows, finalThreeOptRows); err != nil {
			return err
		}
		fmt.Printf("Final+3opt comparison report saved to %s\n", threeOptComparisonPath)
	}
	return nil
}

func runFinalResultsAnalysis(atspsData []AtspData, structuralAnalyses []structuralComparison.InstanceAnalysis, resultsRootPath string) (string, []finalResultsSummaryRow, bool, error) {
	finalAtspsData := make([]AtspData, 0, len(atspsData))
	missingInstances := make([]string, 0)

	for _, atspData := range atspsData {
		finalAtspData := withExperimentOutputRoot(atspData, resultsRootPath)
		if _, err := os.Stat(finalAtspData.resultFilePath); err != nil {
			if os.IsNotExist(err) {
				missingInstances = append(missingInstances, atspData.name)
				continue
			}

			return "", nil, false, err
		}

		finalAtspsData = append(finalAtspsData, finalAtspData)
	}

	if len(finalAtspsData) == 0 {
		fmt.Printf("Final results summary skipped: no final result files found in %s\n", resultsRootPath)
		return "", nil, false, nil
	}

	if len(missingInstances) != 0 {
		return "", nil, false, fmt.Errorf("cannot create final results summary; missing final result.csv for: %s", strings.Join(missingInstances, ", "))
	}

	finalRows, err := readFinalResultsSummaryRows(finalAtspsData)
	if err != nil {
		return "", nil, false, err
	}

	finalResultsSummaryPath := filepath.Join(resultsRootPath, "summary.md")
	if err := saveFinalResultsSummaryRows(finalRows, finalResultsSummaryPath); err != nil {
		return "", nil, false, err
	}
	if err := saveFinalPairwisePerformanceReport(filepath.Join(resultsRootPath, "pairwise_performance.md"), finalRows); err != nil {
		return "", nil, false, err
	}
	if err := saveFinalConvergenceSummaryReport(filepath.Join(resultsRootPath, "convergence_summary.md"), finalRows); err != nil {
		return "", nil, false, err
	}
	if len(structuralAnalyses) != 0 {
		if err := saveStructuralPerformanceLinkReport(filepath.Join(resultsRootPath, "structural_performance_link.md"), finalRows, structuralAnalyses); err != nil {
			return "", nil, false, err
		}
	}
	if err := removeFileIfExists(filepath.Join(resultsRootPath, "summary.csv")); err != nil {
		return "", nil, false, err
	}

	return finalResultsSummaryPath, finalRows, true, nil
}

func saveExperimentPlots(statistics []ExperimentsDataStatistics, plotTitle, plotPathPrefix string) {
	bestStatistic := statistics[0]
	includeMsaPatchBias := shouldIncludeMsaPatchBiasInPlots(statistics)

	for _, statistic := range statistics {
		if statistic.alpha != bestStatistic.alpha ||
			statistic.beta != bestStatistic.beta ||
			statistic.rho != bestStatistic.rho {
			continue
		}

		minDeviationPlotData := utilities.LinePlotData{Name: "min deviation", Color: color.RGBA{G: 255, A: 255}, Values: statistic.minDeviationPerIteration}
		avgDeviationPlotData := utilities.LinePlotData{Name: "avg deviation", Color: color.RGBA{B: 255, A: 255}, Values: statistic.averageDeviationPerIteration}
		maxDeviationPlotData := utilities.LinePlotData{Name: "max deviation", Color: color.RGBA{R: 255, A: 255}, Values: statistic.maxDeviationPerIteration}
		lines := []utilities.LinePlotData{minDeviationPlotData, avgDeviationPlotData, maxDeviationPlotData}

		titleParameters := fmt.Sprintf("alpha=%.2f, beta=%.2f, rho=%.2f, heuristicWeight=%.2f",
			statistic.alpha, statistic.beta, statistic.rho, statistic.heuristicWeight)

		heuristicWeightPlotSuffix := "_heuristicWeight=" + strconv.Itoa(int(100*statistic.heuristicWeight)) + "%"
		plotPath := plotPathPrefix + heuristicWeightPlotSuffix
		if includeMsaPatchBias {
			titleParameters += fmt.Sprintf(", msaPatchBias=%.2f", statistic.msaPatchBias)
			plotPath += "_msaPatchBias=" + strconv.Itoa(int(100*statistic.msaPatchBias)) + "%"
		}
		plotPath += ".png"

		utilities.SaveLinePlotFromData(lines, plotTitle+" ("+titleParameters+")", plotPath)
	}
}

func shouldIncludeMsaPatchBiasInPlots(statistics []ExperimentsDataStatistics) bool {
	seen := make(map[float64]struct{})
	for _, statistic := range statistics {
		seen[statistic.msaPatchBias] = struct{}{}
		if statistic.msaPatchBias != 0 {
			return true
		}
	}

	return len(seen) > 1
}

func removeExperimentPlotsForHeuristic(atspData AtspData, heuristic string) error {
	pattern := resultPlotFilePrefixForHeuristic(atspData, heuristic) + "_heuristicWeight=*.png"
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return err
	}

	for _, match := range matches {
		if err := os.Remove(match); err != nil {
			return err
		}
	}

	return nil
}
