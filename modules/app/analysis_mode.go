package app

import (
	"atsp_aco_msa/modules/analysis/msaHeuristicTours"
	"atsp_aco_msa/modules/analysis/structuralComparison"
	"atsp_aco_msa/modules/analysis/tuningSummary"
	"atsp_aco_msa/modules/project"
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
		if err := saveTuningSummary(project.ResultsDirectoryName, atspsData, tuningHeuristics); err != nil {
			return err
		}
		fmt.Printf("Tuning summary saved to %s\n", filepath.Join(project.ResultsDirectoryName, tuningSummary.ReportFileName))
		return nil
	}

	if analysisScope == analysisScopeAll {
		if err := saveTuningSummary(project.ResultsDirectoryName, atspsData, tuningHeuristics); err != nil {
			return err
		}
		fmt.Printf("Tuning summary saved to %s\n", filepath.Join(project.ResultsDirectoryName, tuningSummary.ReportFileName))
	}

	gksDeviationReportPath := filepath.Join(project.FinalResultsDirectoryName, "gks_deviation.md")
	gksDeviationAtspData, err := project.LoadSelectedAtspData(instanceSetAllKnown)
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
	if err := ensureCycleCoverCache(atspsData, workers); err != nil {
		return err
	}

	msaHeuristicTourConfigs := make([]msaHeuristicTours.InstanceConfig, 0, len(atspsData))
	structuralConfigs := make([]structuralComparison.InstanceConfig, 0, len(atspsData))
	for _, atspData := range atspsData {
		msaHeuristicTourConfigs = append(msaHeuristicTourConfigs, msaHeuristicTours.InstanceConfig{
			Name:                                atspData.Name,
			Dimension:                           len(atspData.Matrix),
			MsaHeuristicDirectoryPath:           atspData.MsaHeuristicDirectoryPath,
			OptimalToursCsvPath:                 atspData.OptimalUniqueToursCsvPath,
			ToursHeatmapPath:                    atspData.ToursHeatmapPlotPath,
			ToursHistogramPath:                  atspData.ToursHistogramPlotPath,
			MsaHeuristicToursOverlapHeatmapPath: atspData.MsaHeuristicToursOverlapHeatmapPlotPath,
		})

		structuralConfigs = append(structuralConfigs, structuralComparison.InstanceConfig{
			Name:                      atspData.Name,
			Dimension:                 len(atspData.Matrix),
			Matrix:                    atspData.Matrix,
			CycleCoverDirectoryPath:   atspData.CycleCoverDirectoryPath,
			MsaHeuristicDirectoryPath: atspData.MsaHeuristicDirectoryPath,
			OptimalToursCsvPath:       atspData.OptimalUniqueToursCsvPath,
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

	structuralSimilarityReportPath := filepath.Join(project.FinalResultsDirectoryName, "structural_similarity.md")
	if err := saveStructuralSimilarityReport(structuralSimilarityReportPath, structuralAnalyses); err != nil {
		return err
	}

	heuristicOverlapReportPath := filepath.Join(project.FinalResultsDirectoryName, "msa_cycle_cover_overlap.md")
	if err := saveMsaHeuristicCycleCoverOverlapReport(heuristicOverlapReportPath, structuralAnalyses); err != nil {
		return err
	}

	msaCountScalingReportPath := filepath.Join(project.FinalResultsDirectoryName, "msa_count_scaling.md")
	if err := saveMsaCountScalingReport(msaCountScalingReportPath, atspsData); err != nil {
		return err
	}

	finalControlsRootPath := finalControlsResultsRootPath(project.FinalResultsDirectoryName)
	randomSparseControlReportPath := filepath.Join(finalControlsRootPath, "random_sparse_control.md")
	randomSparseControlReportSaved, err := saveRandomSparseControlReport(randomSparseControlReportPath, atspsData, project.FinalResultsDirectoryName, finalControlsRootPath)
	if err != nil {
		return err
	}
	distanceRankedSparseControlReportPath := filepath.Join(finalControlsRootPath, "distance_ranked_sparse_control.md")
	distanceRankedSparseControlReportSaved, err := saveDistanceRankedSparseControlReport(distanceRankedSparseControlReportPath, atspsData, project.FinalResultsDirectoryName, finalControlsRootPath)
	if err != nil {
		return err
	}
	shuffledMsaControlReportPath := filepath.Join(finalControlsRootPath, "shuffled_msa_control.md")
	shuffledMsaControlReportSaved, err := saveShuffledMsaControlReport(shuffledMsaControlReportPath, atspsData, project.FinalResultsDirectoryName, finalControlsRootPath)
	if err != nil {
		return err
	}
	if err := removeFileIfExists(filepath.Join(project.FinalResultsDirectoryName, "random_sparse_control.md")); err != nil {
		return err
	}
	if err := removeFileIfExists(filepath.Join(project.FinalResultsDirectoryName, "distance_ranked_sparse_control.md")); err != nil {
		return err
	}
	if err := removeFileIfExists(filepath.Join(project.FinalResultsDirectoryName, "shuffled_msa_control.md")); err != nil {
		return err
	}

	if err := ensureMsaHeuristicCache(gksDeviationAtspData, workers, false); err != nil {
		return err
	}

	if err := saveGksDeviationReport(gksDeviationReportPath, gksDeviationAtspData, gksDeviationMsaPatchBiases); err != nil {
		return err
	}

	finalResultsSummaryPath, finalRows, finalSummarySaved, err := runFinalResultsAnalysis(atspsData, structuralAnalyses, project.FinalResultsDirectoryName)
	if err != nil {
		return err
	}

	finalThreeOptResultsSummaryPath, finalThreeOptRows, finalThreeOptSummarySaved, err := runFinalResultsAnalysis(atspsData, structuralAnalyses, project.FinalThreeOptResultsDirectoryName)
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
		fmt.Printf("Pairwise performance report saved to %s\n", filepath.Join(project.FinalResultsDirectoryName, "pairwise_performance.md"))
		fmt.Printf("Convergence summary report saved to %s\n", filepath.Join(project.FinalResultsDirectoryName, "convergence_summary.md"))
		fmt.Printf("Structural/performance link report saved to %s\n", filepath.Join(project.FinalResultsDirectoryName, "structural_performance_link.md"))
	}
	if finalThreeOptSummarySaved {
		fmt.Printf("Final+3opt results summary saved to %s\n", finalThreeOptResultsSummaryPath)
		fmt.Printf("Final+3opt pairwise performance report saved to %s\n", filepath.Join(project.FinalThreeOptResultsDirectoryName, "pairwise_performance.md"))
		fmt.Printf("Final+3opt convergence summary report saved to %s\n", filepath.Join(project.FinalThreeOptResultsDirectoryName, "convergence_summary.md"))
		fmt.Printf("Final+3opt structural/performance link report saved to %s\n", filepath.Join(project.FinalThreeOptResultsDirectoryName, "structural_performance_link.md"))
	}
	if finalSummarySaved && finalThreeOptSummarySaved {
		threeOptComparisonPath := filepath.Join(project.FinalThreeOptResultsDirectoryName, "comparison_to_final.md")
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
		finalAtspData := project.WithExperimentOutputRoot(atspData, resultsRootPath)
		if _, err := os.Stat(finalAtspData.ResultFilePath); err != nil {
			if os.IsNotExist(err) {
				missingInstances = append(missingInstances, atspData.Name)
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
		if statistic.Alpha != bestStatistic.Alpha ||
			statistic.Beta != bestStatistic.Beta ||
			statistic.Rho != bestStatistic.Rho {
			continue
		}

		minDeviationPlotData := utilities.LinePlotData{Name: "min deviation", Color: color.RGBA{G: 255, A: 255}, Values: statistic.MinDeviationPerIteration}
		avgDeviationPlotData := utilities.LinePlotData{Name: "avg deviation", Color: color.RGBA{B: 255, A: 255}, Values: statistic.AverageDeviationPerIteration}
		maxDeviationPlotData := utilities.LinePlotData{Name: "max deviation", Color: color.RGBA{R: 255, A: 255}, Values: statistic.MaxDeviationPerIteration}
		lines := []utilities.LinePlotData{minDeviationPlotData, avgDeviationPlotData, maxDeviationPlotData}

		titleParameters := fmt.Sprintf("alpha=%.2f, beta=%.2f, rho=%.2f, heuristicWeight=%.2f",
			statistic.Alpha, statistic.Beta, statistic.Rho, statistic.HeuristicWeight)

		heuristicWeightPlotSuffix := "_heuristicWeight=" + strconv.Itoa(int(100*statistic.HeuristicWeight)) + "%"
		plotPath := plotPathPrefix + heuristicWeightPlotSuffix
		if includeMsaPatchBias {
			titleParameters += fmt.Sprintf(", msaPatchBias=%.2f", statistic.MsaPatchBias)
			plotPath += "_msaPatchBias=" + strconv.Itoa(int(100*statistic.MsaPatchBias)) + "%"
		}
		plotPath += ".png"

		utilities.SaveLinePlotFromData(lines, plotTitle+" ("+titleParameters+")", plotPath)
	}
}

func shouldIncludeMsaPatchBiasInPlots(statistics []ExperimentsDataStatistics) bool {
	seen := make(map[float64]struct{})
	for _, statistic := range statistics {
		seen[statistic.MsaPatchBias] = struct{}{}
		if statistic.MsaPatchBias != 0 {
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
