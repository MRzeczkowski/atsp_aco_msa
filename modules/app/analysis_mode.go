package app

import (
	"atsp_aco_msa/modules/analysis/reports"
	"atsp_aco_msa/modules/analysis/structure"
	"atsp_aco_msa/modules/analysis/tours"
	"atsp_aco_msa/modules/analysis/tuning"
	"atsp_aco_msa/modules/experiments"
	"atsp_aco_msa/modules/project"
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"image/color"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

func evaluationReportsConfig() reports.EvaluationReportsConfig {
	return reports.EvaluationReportsConfig{
		Heuristics:                     evaluationResultsSummaryHeuristics,
		BaselineHeuristic:              heuristicBaseline,
		StrictMsaHeuristic:             heuristicStrictMsa,
		CycleCoverHeuristic:            heuristicCycleCover,
		CycleCoverMsaPatchingHeuristic: heuristicCycleCoverMsaPatching,
		DisplayName:                    heuristicDisplayName,
	}
}

func controlReportsConfig() reports.ControlReportsConfig {
	return reports.ControlReportsConfig{
		StrictMsaHeuristic:                 heuristicStrictMsa,
		RandomSparseHeuristic:              heuristicRandomSparse,
		DistanceRankedSparseHeuristic:      heuristicDistanceRankedSparse,
		ShuffledMsaHeuristic:               heuristicShuffledMsa,
		EvaluationStrictMsaHeuristicWeight: evaluationStrictMsaHeuristicWeight,
		ReadStatistics:                     experiments.ReadStatistics,
		ReadHeuristicStatistics:            experiments.ReadHeuristicStatistics,
		ResultFilePathForHeuristic:         resultFilePathForHeuristic,
	}
}

func runAnalysisMode(atspsData []AtspData, analysisScope string, tuningHeuristics []string, workers int) error {
	if analysisScope == analysisScopeTuning {
		if err := saveTuningSummary(project.ResultsDirectoryName, atspsData, tuningHeuristics); err != nil {
			return err
		}
		fmt.Printf("Tuning summary saved to %s\n", filepath.Join(project.ResultsDirectoryName, tuning.ReportFileName))
		return nil
	}

	if analysisScope == analysisScopeAll {
		if err := saveTuningSummary(project.ResultsDirectoryName, atspsData, tuningHeuristics); err != nil {
			return err
		}
		fmt.Printf("Tuning summary saved to %s\n", filepath.Join(project.ResultsDirectoryName, tuning.ReportFileName))
	}

	gksDeviationReportPath := filepath.Join(project.EvaluationResultsDirectoryName, "gks_deviation.md")
	gksDeviationAtspData, err := project.LoadSelectedAtspData(instanceSetAllKnown)
	if err != nil {
		return err
	}

	if analysisScope == analysisScopeGksDeviation {
		if err := ensureMsaHeuristicCache(gksDeviationAtspData, workers, false); err != nil {
			return err
		}
		if err := ensureCycleCoverCache(gksDeviationAtspData, workers); err != nil {
			return err
		}
		if err := reports.SaveGksDeviation(gksDeviationReportPath, gksDeviationAtspData, gksDeviationMsaPatchBiases); err != nil {
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

	solutionTourConfigs := make([]tours.InstanceConfig, 0, len(atspsData))
	structuralConfigs := make([]structure.InstanceConfig, 0, len(atspsData))
	for _, atspData := range atspsData {
		solutionTourConfigs = append(solutionTourConfigs, tours.InstanceConfig{
			Name:                                atspData.Name,
			Dimension:                           len(atspData.Matrix),
			MsaHeuristicDirectoryPath:           atspData.MsaHeuristicDirectoryPath,
			CycleCoverDirectoryPath:             atspData.CycleCoverDirectoryPath,
			OptimalToursCsvPath:                 atspData.OptimalUniqueToursCsvPath,
			ToursHeatmapPath:                    atspData.ToursHeatmapPlotPath,
			ToursHistogramPath:                  atspData.ToursHistogramPlotPath,
			MsaHeuristicToursOverlapHeatmapPath: atspData.MsaHeuristicToursOverlapHeatmapPlotPath,
			CycleCoverToursOverlapHeatmapPath:   atspData.CycleCoverToursOverlapHeatmapPlotPath,
		})

		structuralConfigs = append(structuralConfigs, structure.InstanceConfig{
			Name:                      atspData.Name,
			Dimension:                 len(atspData.Matrix),
			Matrix:                    atspData.Matrix,
			CycleCoverDirectoryPath:   atspData.CycleCoverDirectoryPath,
			MsaHeuristicDirectoryPath: atspData.MsaHeuristicDirectoryPath,
			OptimalToursCsvPath:       atspData.OptimalUniqueToursCsvPath,
		})
	}

	if err := tours.AnalyzeInstances(tours.Config{Instances: solutionTourConfigs}); err != nil {
		return err
	}

	structuralAnalyses, err := structure.AnalyzeInstances(structure.Config{
		Instances:     structuralConfigs,
		HighThreshold: 1.0,
		MsaPatchBias:  evaluationCycleCoverMsaPatchingMsaPatchBias,
	})
	if err != nil {
		return err
	}

	structuralSimilarityReportPath := filepath.Join(project.EvaluationResultsDirectoryName, "structural_similarity.md")
	if err := reports.SaveStructuralSimilarity(structuralSimilarityReportPath, structuralAnalyses); err != nil {
		return err
	}

	heuristicOverlapReportPath := filepath.Join(project.EvaluationResultsDirectoryName, "msa_cycle_cover_overlap.md")
	if err := reports.SaveMsaHeuristicCycleCoverOverlap(heuristicOverlapReportPath, structuralAnalyses); err != nil {
		return err
	}

	msaCountScalingReportPath := filepath.Join(project.EvaluationResultsDirectoryName, "msa_count_scaling.md")
	if err := reports.SaveMsaCountScaling(msaCountScalingReportPath, atspsData, msaCountScalingCounts); err != nil {
		return err
	}

	evaluationControlsRootPath := evaluationControlsResultsRootPath(project.EvaluationResultsDirectoryName)
	randomSparseControlReportPath := filepath.Join(evaluationControlsRootPath, "random_sparse_control.md")
	randomSparseControlReportSaved, err := reports.SaveRandomSparseControlReport(randomSparseControlReportPath, atspsData, project.EvaluationResultsDirectoryName, evaluationControlsRootPath, controlReportsConfig())
	if err != nil {
		return err
	}
	distanceRankedSparseControlReportPath := filepath.Join(evaluationControlsRootPath, "distance_ranked_sparse_control.md")
	distanceRankedSparseControlReportSaved, err := reports.SaveDistanceRankedSparseControlReport(distanceRankedSparseControlReportPath, atspsData, project.EvaluationResultsDirectoryName, evaluationControlsRootPath, controlReportsConfig())
	if err != nil {
		return err
	}
	shuffledMsaControlReportPath := filepath.Join(evaluationControlsRootPath, "shuffled_msa_control.md")
	shuffledMsaControlReportSaved, err := reports.SaveShuffledMsaControlReport(shuffledMsaControlReportPath, atspsData, project.EvaluationResultsDirectoryName, evaluationControlsRootPath, controlReportsConfig())
	if err != nil {
		return err
	}

	if err := ensureMsaHeuristicCache(gksDeviationAtspData, workers, false); err != nil {
		return err
	}
	if err := ensureCycleCoverCache(gksDeviationAtspData, workers); err != nil {
		return err
	}

	if err := reports.SaveGksDeviation(gksDeviationReportPath, gksDeviationAtspData, gksDeviationMsaPatchBiases); err != nil {
		return err
	}

	evaluationResultsSummaryPath, evaluationRows, evaluationSummarySaved, err := runEvaluationResultsAnalysis(atspsData, structuralAnalyses, project.EvaluationResultsDirectoryName)
	if err != nil {
		return err
	}

	evaluationThreeOptResultsSummaryPath, evaluationThreeOptRows, evaluationThreeOptSummarySaved, err := runEvaluationResultsAnalysis(atspsData, structuralAnalyses, project.EvaluationThreeOptResultsDirectoryName)
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
	if evaluationSummarySaved {
		fmt.Printf("Evaluation results summary saved to %s\n", evaluationResultsSummaryPath)
		fmt.Printf("Pairwise performance report saved to %s\n", filepath.Join(project.EvaluationResultsDirectoryName, "pairwise_performance.md"))
		fmt.Printf("Convergence summary report saved to %s\n", filepath.Join(project.EvaluationResultsDirectoryName, "convergence_summary.md"))
		fmt.Printf("Structural/performance link report saved to %s\n", filepath.Join(project.EvaluationResultsDirectoryName, "structural_performance_link.md"))
	}
	if evaluationThreeOptSummarySaved {
		fmt.Printf("Evaluation+3opt results summary saved to %s\n", evaluationThreeOptResultsSummaryPath)
		fmt.Printf("Evaluation+3opt pairwise performance report saved to %s\n", filepath.Join(project.EvaluationThreeOptResultsDirectoryName, "pairwise_performance.md"))
		fmt.Printf("Evaluation+3opt convergence summary report saved to %s\n", filepath.Join(project.EvaluationThreeOptResultsDirectoryName, "convergence_summary.md"))
		fmt.Printf("Evaluation+3opt structural/performance link report saved to %s\n", filepath.Join(project.EvaluationThreeOptResultsDirectoryName, "structural_performance_link.md"))
	}
	if evaluationSummarySaved && evaluationThreeOptSummarySaved {
		threeOptComparisonPath := filepath.Join(project.EvaluationThreeOptResultsDirectoryName, "comparison_to_evaluation.md")
		if err := reports.SaveEvaluationThreeOptComparisonReport(threeOptComparisonPath, evaluationRows, evaluationThreeOptRows, evaluationReportsConfig()); err != nil {
			return err
		}
		fmt.Printf("Evaluation+3opt comparison report saved to %s\n", threeOptComparisonPath)
	}
	return nil
}

func runEvaluationResultsAnalysis(atspsData []AtspData, structuralAnalyses []structure.InstanceAnalysis, resultsRootPath string) (string, []reports.EvaluationResultsSummaryRow, bool, error) {
	evaluationAtspsData := make([]AtspData, 0, len(atspsData))
	missingInstances := make([]string, 0)

	for _, atspData := range atspsData {
		evaluationAtspData := project.WithExperimentOutputRoot(atspData, resultsRootPath)
		if _, err := os.Stat(evaluationAtspData.ResultFilePath); err != nil {
			if os.IsNotExist(err) {
				missingInstances = append(missingInstances, atspData.Name)
				continue
			}

			return "", nil, false, err
		}

		evaluationAtspsData = append(evaluationAtspsData, evaluationAtspData)
	}

	if len(evaluationAtspsData) == 0 {
		fmt.Printf("Evaluation results summary skipped: no evaluation result files found in %s\n", resultsRootPath)
		return "", nil, false, nil
	}

	if len(missingInstances) != 0 {
		return "", nil, false, fmt.Errorf("cannot create evaluation results summary; missing evaluation result.csv for: %s", strings.Join(missingInstances, ", "))
	}

	evaluationRows, err := reports.ReadEvaluationResultsSummaryRows(evaluationAtspsData)
	if err != nil {
		return "", nil, false, err
	}

	evaluationResultsSummaryPath := filepath.Join(resultsRootPath, "summary.md")
	if err := reports.SaveEvaluationResultsSummaryRows(evaluationRows, evaluationResultsSummaryPath, evaluationReportsConfig()); err != nil {
		return "", nil, false, err
	}
	if err := reports.SaveEvaluationPairwisePerformanceReport(filepath.Join(resultsRootPath, "pairwise_performance.md"), evaluationRows, evaluationReportsConfig()); err != nil {
		return "", nil, false, err
	}
	if err := reports.SaveEvaluationConvergenceSummaryReport(filepath.Join(resultsRootPath, "convergence_summary.md"), evaluationRows, evaluationReportsConfig()); err != nil {
		return "", nil, false, err
	}
	if len(structuralAnalyses) != 0 {
		if err := reports.SaveStructuralPerformanceLinkReport(filepath.Join(resultsRootPath, "structural_performance_link.md"), evaluationRows, structuralAnalyses, evaluationReportsConfig()); err != nil {
			return "", nil, false, err
		}
	}

	return evaluationResultsSummaryPath, evaluationRows, true, nil
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
