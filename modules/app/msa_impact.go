package app

import (
	"atsp_aco_msa/modules/algorithms/heuristics"
	"atsp_aco_msa/modules/analysis/msaHeuristicTours"
	"atsp_aco_msa/modules/models"
	"errors"
	"fmt"
	"html"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// This file contains the temporary msa-impact prototype pipeline. It is kept
// working for now, but it is intentionally isolated because these experiments
// are expected to be removed in a later cleanup pass.

// Temporary prototype pipeline for quick MSA heuristic experiments and controls.
// This is deliberately separate from tuning/final runs and should be removed
// together with artifacts/msa_impact after the MSA integration experiments end.
func runMsaImpactMode(atspsData []AtspData, workers int) error {
	workers = maxIntValue(1, workers)
	configurations := msaImpactExperimentConfigurations()

	if finalConfigurationsUseMsaHeuristic(configurations) {
		if err := ensureMsaHeuristicCache(atspsData, workers, finalConfigurationsUseRootedMsa(configurations)); err != nil {
			return err
		}
	}
	if err := removeLegacyFinalReports(msaImpactResultsDirectoryName); err != nil {
		return err
	}
	if err := removeMsaImpactResultFiles(atspsData); err != nil {
		return err
	}

	fmt.Printf("Running MSA impact pipeline with %d run(s) per configuration and up to %d parameter worker(s)\n",
		msaImpactNumberOfExperiments,
		min(workers, finalExperimentParameterSetCount(configurations)))
	if err := runMsaImpactExperimentInstances(atspsData, msaImpactResultsDirectoryName, configurations, workers, msaImpactNumberOfExperiments); err != nil {
		return err
	}
	if err := runMsaImpactControls(atspsData, workers); err != nil {
		return err
	}

	return saveMsaImpactReports(atspsData)
}

func runMsaImpactExperimentInstances(atspsData []AtspData, resultsRootPath string, configurations []finalExperimentConfiguration, workers, numberOfExperiments int) error {
	return runBoundedInstanceJobs(atspsData, 1, func(atspData AtspData) error {
		finalAtspData := withExperimentOutputRoot(atspData, resultsRootPath)
		return runFinalExperimentForInstanceWithParameterWorkers(finalAtspData, resultsRootPath, false, configurations, numberOfExperiments, workers)
	})
}

func runMsaImpactControls(atspsData []AtspData, workers int) error {
	workers = maxIntValue(1, workers)
	if err := removeMsaImpactControlFiles(atspsData); err != nil {
		return err
	}

	return runBoundedInstanceJobs(atspsData, 1, func(atspData AtspData) error {
		strictMetric, err := readBestMsaImpactMetric(atspData, heuristicStrictMsa)
		if err != nil {
			return err
		}
		rootedMetric, err := readBestMsaImpactMetric(atspData, heuristicRootedMsa)
		if err != nil {
			return err
		}

		controlAtspData := withExperimentOutputRoot(atspData, msaImpactControlsDirectoryName)
		configurations := append(
			msaImpactControlExperimentConfigurations(heuristicStrictMsa, strictMetric.heuristicWeight),
			msaImpactControlExperimentConfigurations(heuristicRootedMsa, rootedMetric.heuristicWeight)...)
		return runFinalExperimentForInstanceWithParameterWorkers(controlAtspData, msaImpactControlsDirectoryName, false, configurations, msaImpactNumberOfExperiments, workers)
	})
}

func readBestMsaImpactMetric(atspData AtspData, heuristic string) (finalResultsSummaryMetric, error) {
	impactAtspData := withExperimentOutputRoot(atspData, msaImpactResultsDirectoryName)
	metrics, err := readFinalResultSummaryMetrics(impactAtspData.resultFilePath)
	if err != nil {
		return finalResultsSummaryMetric{}, fmt.Errorf("%s: read MSA impact results: %w", atspData.name, err)
	}

	metric, ok := metrics[heuristic]
	if !ok {
		return finalResultsSummaryMetric{}, fmt.Errorf("%s: missing %s impact heuristic metric", atspData.name, heuristic)
	}

	return metric, nil
}

func removeMsaImpactResultFiles(atspsData []AtspData) error {
	configurations := msaImpactExperimentConfigurations()
	for _, atspData := range atspsData {
		impactAtspData := withExperimentOutputRoot(atspData, msaImpactResultsDirectoryName)
		if err := removeFileIfExists(impactAtspData.resultFilePath); err != nil {
			return err
		}
		for _, configuration := range configurations {
			if err := removeFileIfExists(resultFilePathForHeuristic(impactAtspData, configuration.heuristic)); err != nil {
				return err
			}
		}
	}
	return nil
}

func removeMsaImpactControlFiles(atspsData []AtspData) error {
	controlHeuristics := []string{
		heuristicRandomSparse,
		heuristicDistanceRankedSparse,
		heuristicShuffledMsa,
		heuristicStrictDistanceRanked,
		heuristicStrictShuffledMsa,
		heuristicRootedDistanceRanked,
		heuristicRootedShuffledMsa,
	}
	for _, atspData := range atspsData {
		controlAtspData := withExperimentOutputRoot(atspData, msaImpactControlsDirectoryName)
		for _, heuristic := range controlHeuristics {
			if err := removeFileIfExists(resultFilePathForHeuristic(controlAtspData, heuristic)); err != nil {
				return err
			}
		}
	}

	reportPaths := []string{
		filepath.Join(msaImpactControlsDirectoryName, "random_sparse_control.md"),
		filepath.Join(msaImpactControlsDirectoryName, "distance_ranked_sparse_control.md"),
		filepath.Join(msaImpactControlsDirectoryName, "shuffled_msa_control.md"),
		filepath.Join(msaImpactControlsDirectoryName, "strict_distance_ranked_sparse_control.md"),
		filepath.Join(msaImpactControlsDirectoryName, "strict_shuffled_msa_control.md"),
		filepath.Join(msaImpactControlsDirectoryName, "rooted_distance_ranked_sparse_control.md"),
		filepath.Join(msaImpactControlsDirectoryName, "rooted_shuffled_msa_control.md"),
		filepath.Join(msaImpactControlsDirectoryName, "msa_weight_control_summary.md"),
		filepath.Join(msaImpactControlsDirectoryName, "msa_distance_ranked_edge_categories.md"),
	}
	for _, path := range reportPaths {
		if err := removeFileIfExists(path); err != nil {
			return err
		}
	}

	return nil
}

func saveMsaImpactReports(atspsData []AtspData) error {
	summaryPath := filepath.Join(msaImpactResultsDirectoryName, "msa_impact_summary.md")
	if err := saveMsaImpactSummary(summaryPath, atspsData); err != nil {
		return err
	}
	fmt.Printf("MSA impact summary saved to %s\n", summaryPath)

	structureReportPath := filepath.Join(msaImpactResultsDirectoryName, "msa_impact_structure.md")
	if err := saveMsaImpactStructureReport(structureReportPath, atspsData); err != nil {
		return err
	}
	fmt.Printf("MSA impact structure report saved to %s\n", structureReportPath)

	strictDistanceRankedControlReportPath := filepath.Join(msaImpactControlsDirectoryName, "strict_distance_ranked_sparse_control.md")
	strictDistanceRankedControlReportSaved, err := saveDistanceRankedSparseControlReportForReference(strictDistanceRankedControlReportPath, atspsData, msaImpactResultsDirectoryName, msaImpactControlsDirectoryName, heuristicStrictMsa, "Strict MSA", heuristicStrictDistanceRanked, "strict distance-ranked sparse")
	if err != nil {
		return err
	}
	strictShuffledControlReportPath := filepath.Join(msaImpactControlsDirectoryName, "strict_shuffled_msa_control.md")
	strictShuffledControlReportSaved, err := saveSeededSparseControlReport(strictShuffledControlReportPath, atspsData, msaImpactResultsDirectoryName, msaImpactControlsDirectoryName, seededSparseControlReportConfig{
		title:              "Strict MSA Shuffled Control",
		description:        fmt.Sprintf("This sanity check compares Strict MSA against deterministic shuffles of the strict MSA mask. Each shuffle preserves the number of strict-MSA boosted directed edges, but assigns them to shuffled directed edges. The control %s.", controlWeightDescription(msaImpactResultsDirectoryName, "Strict MSA")),
		referenceHeuristic: heuristicStrictMsa,
		referenceName:      "Strict MSA",
		controlName:        "strict shuffled MSA",
		controlHeuristic:   heuristicStrictShuffledMsa,
		missingResultLabel: "strict-shuffled-MSA result CSV",
	})
	if err != nil {
		return err
	}
	rootedDistanceRankedControlReportPath := filepath.Join(msaImpactControlsDirectoryName, "rooted_distance_ranked_sparse_control.md")
	rootedDistanceRankedControlReportSaved, err := saveDistanceRankedSparseControlReportForReference(rootedDistanceRankedControlReportPath, atspsData, msaImpactResultsDirectoryName, msaImpactControlsDirectoryName, heuristicRootedMsa, "Rooted MSA", heuristicRootedDistanceRanked, "rooted distance-ranked sparse")
	if err != nil {
		return err
	}
	rootedShuffledControlReportPath := filepath.Join(msaImpactControlsDirectoryName, "rooted_shuffled_msa_control.md")
	rootedShuffledControlReportSaved, err := saveSeededSparseControlReport(rootedShuffledControlReportPath, atspsData, msaImpactResultsDirectoryName, msaImpactControlsDirectoryName, seededSparseControlReportConfig{
		title:              "Rooted MSA Shuffled Control",
		description:        fmt.Sprintf("This sanity check compares Rooted MSA against deterministic per-root shuffles of the rooted MSA masks. Each root preserves the number of boosted directed edges from that root's MSA, but assigns them to shuffled directed edges. The control %s.", controlWeightDescription(msaImpactResultsDirectoryName, "Rooted MSA")),
		referenceHeuristic: heuristicRootedMsa,
		referenceName:      "Rooted MSA",
		controlName:        "rooted shuffled MSA",
		controlHeuristic:   heuristicRootedShuffledMsa,
		missingResultLabel: "rooted-shuffled-MSA result CSV",
	})
	if err != nil {
		return err
	}
	weightControlReportPath := filepath.Join(msaImpactControlsDirectoryName, "msa_weight_control_summary.md")
	weightControlReportSaved, err := saveMsaImpactControlWeightSummary(weightControlReportPath, atspsData)
	if err != nil {
		return err
	}
	distanceRankedCategoryReportPath := filepath.Join(msaImpactControlsDirectoryName, "msa_distance_ranked_edge_categories.md")
	if err := saveMsaDistanceRankedCategoryReport(distanceRankedCategoryReportPath, atspsData); err != nil {
		return err
	}

	if strictDistanceRankedControlReportSaved {
		fmt.Printf("Strict distance-ranked sparse control report saved to %s\n", strictDistanceRankedControlReportPath)
	}
	if strictShuffledControlReportSaved {
		fmt.Printf("Strict shuffled MSA control report saved to %s\n", strictShuffledControlReportPath)
	}
	if rootedDistanceRankedControlReportSaved {
		fmt.Printf("Rooted distance-ranked sparse control report saved to %s\n", rootedDistanceRankedControlReportPath)
	}
	if rootedShuffledControlReportSaved {
		fmt.Printf("Rooted shuffled MSA control report saved to %s\n", rootedShuffledControlReportPath)
	}
	if weightControlReportSaved {
		fmt.Printf("MSA weight control summary saved to %s\n", weightControlReportPath)
	}
	fmt.Printf("MSA/distance-ranked edge category report saved to %s\n", distanceRankedCategoryReportPath)
	return nil
}

type msaImpactControlWeightRow struct {
	heuristicWeight             float64
	instance                    string
	msaAverageBestDeviation     float64
	controlAverageBestDeviation float64
	averageBestDeviationDelta   float64
	msaSuccessRate              float64
	controlSuccessRate          float64
	successRateDelta            float64
}

type msaImpactControlWeightSummaryRow struct {
	heuristicWeight                 float64
	instanceCount                   int
	msaWins                         int
	controlWins                     int
	ties                            int
	meanMsaAverageBestDeviation     float64
	meanControlAverageBestDeviation float64
	meanAverageBestDeviationDelta   float64
	meanMsaSuccessRate              float64
	meanControlSuccessRate          float64
	meanSuccessRateDelta            float64
	signTestPValue                  float64
}

func saveMsaImpactControlWeightSummary(path string, atspsData []AtspData) (bool, error) {
	strictDistanceRankedRows, err := buildMsaImpactSingleControlWeightRows(atspsData, heuristicStrictMsa, heuristicStrictDistanceRanked)
	if err != nil {
		return false, err
	}
	strictShuffledRows, err := buildMsaImpactSeededControlWeightRows(atspsData, heuristicStrictMsa, heuristicStrictShuffledMsa)
	if err != nil {
		return false, err
	}
	rootedDistanceRankedRows, err := buildMsaImpactSingleControlWeightRows(atspsData, heuristicRootedMsa, heuristicRootedDistanceRanked)
	if err != nil {
		return false, err
	}
	rootedShuffledRows, err := buildMsaImpactSeededControlWeightRows(atspsData, heuristicRootedMsa, heuristicRootedShuffledMsa)
	if err != nil {
		return false, err
	}
	if len(strictDistanceRankedRows) == 0 && len(strictShuffledRows) == 0 && len(rootedDistanceRankedRows) == 0 && len(rootedShuffledRows) == 0 {
		return false, nil
	}
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return false, err
	}

	var builder strings.Builder
	builder.WriteString("# MSA Weight Control Summary\n\n")
	builder.WriteString("This report checks whether each MSA-impact variant keeps beating its matched control heuristics at the best variant weight selected for each instance. Rows are grouped by weight because different instances may select different best weights. The p-value is a two-sided sign test over per-instance average-best-deviation wins and losses; ties are ignored.\n\n")
	writeMsaImpactControlWeightTable(&builder, "Strict MSA vs Strict Distance-ranked Sparse", "Strict MSA", "Strict distance-ranked sparse", strictDistanceRankedRows)
	writeMsaImpactControlWeightTable(&builder, "Strict MSA vs Strict Shuffled MSA", "Strict MSA", "Strict shuffled MSA", strictShuffledRows)
	writeMsaImpactControlWeightTable(&builder, "Rooted MSA vs Rooted Distance-ranked Sparse", "Rooted MSA", "Rooted distance-ranked sparse", rootedDistanceRankedRows)
	writeMsaImpactControlWeightTable(&builder, "Rooted MSA vs Rooted Shuffled MSA", "Rooted MSA", "Rooted shuffled MSA", rootedShuffledRows)

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildMsaImpactSeededControlWeightRows(atspsData []AtspData, referenceHeuristic, controlHeuristic string) ([]msaImpactControlWeightRow, error) {
	rows := make([]msaImpactControlWeightRow, 0, len(atspsData))
	for _, atspData := range atspsData {
		msaMetric, ok, err := readMsaControlMetric(atspData, msaImpactResultsDirectoryName, referenceHeuristic)
		if err != nil {
			return nil, err
		}
		if !ok {
			continue
		}

		controlAtspData := withExperimentOutputRoot(atspData, msaImpactControlsDirectoryName)
		controlStatistics, err := readStatistics(resultFilePathForHeuristic(controlAtspData, controlHeuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				continue
			}
			return nil, err
		}

		controlStatisticsForWeight := statisticsForHeuristicWeightAll(controlStatistics, msaMetric.heuristicWeight)
		if len(controlStatisticsForWeight) == 0 {
			continue
		}

		controlAverageBestDeviation, controlSuccessRate := averageControlStatistics(controlStatisticsForWeight)
		rows = append(rows, msaImpactControlWeightRow{
			heuristicWeight:             msaMetric.heuristicWeight,
			instance:                    atspData.name,
			msaAverageBestDeviation:     msaMetric.averageMinDeviation,
			controlAverageBestDeviation: controlAverageBestDeviation,
			averageBestDeviationDelta:   msaMetric.averageMinDeviation - controlAverageBestDeviation,
			msaSuccessRate:              msaMetric.successRate,
			controlSuccessRate:          controlSuccessRate,
			successRateDelta:            msaMetric.successRate - controlSuccessRate,
		})
	}

	return rows, nil
}

func buildMsaImpactSingleControlWeightRows(atspsData []AtspData, referenceHeuristic, controlHeuristic string) ([]msaImpactControlWeightRow, error) {
	rows := make([]msaImpactControlWeightRow, 0, len(atspsData))
	for _, atspData := range atspsData {
		msaMetric, ok, err := readMsaControlMetric(atspData, msaImpactResultsDirectoryName, referenceHeuristic)
		if err != nil {
			return nil, err
		}
		if !ok {
			continue
		}

		controlAtspData := withExperimentOutputRoot(atspData, msaImpactControlsDirectoryName)
		controlStatistics, err := readStatistics(resultFilePathForHeuristic(controlAtspData, controlHeuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				continue
			}
			return nil, err
		}

		controlStatistic, ok := statisticsForHeuristicWeight(controlStatistics, msaMetric.heuristicWeight)
		if !ok {
			continue
		}

		rows = append(rows, msaImpactControlWeightRow{
			heuristicWeight:             msaMetric.heuristicWeight,
			instance:                    atspData.name,
			msaAverageBestDeviation:     msaMetric.averageMinDeviation,
			controlAverageBestDeviation: controlStatistic.averageBestDeviation,
			averageBestDeviationDelta:   msaMetric.averageMinDeviation - controlStatistic.averageBestDeviation,
			msaSuccessRate:              msaMetric.successRate,
			controlSuccessRate:          controlStatistic.successRate,
			successRateDelta:            msaMetric.successRate - controlStatistic.successRate,
		})
	}

	return rows, nil
}

func readMsaImpactMsaMetrics(atspData AtspData) ([]finalResultsSummaryMetric, error) {
	impactAtspData := withExperimentOutputRoot(atspData, msaImpactResultsDirectoryName)
	statistics, err := readHeuristicStatistics(impactAtspData.resultFilePath)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return nil, nil
		}
		return nil, err
	}

	metrics := make([]finalResultsSummaryMetric, 0, len(statistics))
	for _, statistic := range statistics {
		if statistic.heuristic != heuristicRootedMsa {
			continue
		}

		metric := finalResultsSummaryMetric{
			averageMinDeviation:  statistic.statistics.averageBestDeviation,
			successRate:          statistic.statistics.successRate,
			averageBestIteration: statistic.statistics.averageBestAtIteration,
			heuristicWeight:      statistic.statistics.heuristicWeight,
			iterations:           statistic.statistics.iterations,
		}
		current, ok := metricForHeuristicWeight(metrics, metric.heuristicWeight)
		if !ok {
			metrics = append(metrics, metric)
			continue
		}
		if finalResultsSummaryMetricIsBetter(metric, current) {
			for i := range metrics {
				if math.Abs(metrics[i].heuristicWeight-metric.heuristicWeight) < 1e-9 {
					metrics[i] = metric
					break
				}
			}
		}
	}

	return metrics, nil
}

func metricForHeuristicWeight(metrics []finalResultsSummaryMetric, heuristicWeight float64) (finalResultsSummaryMetric, bool) {
	for _, metric := range metrics {
		if math.Abs(metric.heuristicWeight-heuristicWeight) < 1e-9 {
			return metric, true
		}
	}

	return finalResultsSummaryMetric{}, false
}

func averageControlStatistics(statistics []ExperimentsDataStatistics) (averageBestDeviation, successRate float64) {
	for _, statistic := range statistics {
		averageBestDeviation += statistic.averageBestDeviation
		successRate += statistic.successRate
	}
	count := float64(len(statistics))
	return averageBestDeviation / count, successRate / count
}

func writeMsaImpactControlWeightTable(builder *strings.Builder, title, referenceName, controlName string, rows []msaImpactControlWeightRow) {
	builder.WriteString("## " + title + "\n\n")
	if len(rows) == 0 {
		builder.WriteString("No comparable rows were found.\n\n")
		return
	}

	summaryRows := summarizeMsaImpactControlWeightRows(rows)
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Weight</th><th>Instances</th><th>")
	builder.WriteString(html.EscapeString(referenceName))
	builder.WriteString(" wins</th><th>Control wins</th><th>Ties</th><th>p-value</th><th>")
	builder.WriteString(html.EscapeString(referenceName))
	builder.WriteString(" avg dev. [%]</th><th>")
	builder.WriteString(html.EscapeString(controlName))
	builder.WriteString(" avg dev. [%]</th><th>Delta [pp]</th><th>")
	builder.WriteString(html.EscapeString(referenceName))
	builder.WriteString(" success [%]</th><th>")
	builder.WriteString(html.EscapeString(controlName))
	builder.WriteString(" success [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range summaryRows {
		fmt.Fprintf(builder,
			"<tr><td align=\"right\">%.2f</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%.6f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			row.heuristicWeight,
			row.instanceCount,
			row.msaWins,
			row.controlWins,
			row.ties,
			row.signTestPValue,
			row.meanMsaAverageBestDeviation,
			row.meanControlAverageBestDeviation,
			formatSignedFloat(row.meanAverageBestDeviationDelta),
			row.meanMsaSuccessRate,
			row.meanControlSuccessRate)
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n\n")
}

func summarizeMsaImpactControlWeightRows(rows []msaImpactControlWeightRow) []msaImpactControlWeightSummaryRow {
	summaryRows := make([]msaImpactControlWeightSummaryRow, 0, len(msaImpactHeuristicWeights))
	for _, heuristicWeight := range msaImpactHeuristicWeights {
		var summary msaImpactControlWeightSummaryRow
		summary.heuristicWeight = heuristicWeight
		for _, row := range rows {
			if math.Abs(row.heuristicWeight-heuristicWeight) >= 1e-9 {
				continue
			}

			summary.instanceCount++
			if row.averageBestDeviationDelta < -1e-9 {
				summary.msaWins++
			} else if row.averageBestDeviationDelta > 1e-9 {
				summary.controlWins++
			} else {
				summary.ties++
			}

			summary.meanMsaAverageBestDeviation += row.msaAverageBestDeviation
			summary.meanControlAverageBestDeviation += row.controlAverageBestDeviation
			summary.meanAverageBestDeviationDelta += row.averageBestDeviationDelta
			summary.meanMsaSuccessRate += row.msaSuccessRate
			summary.meanControlSuccessRate += row.controlSuccessRate
			summary.meanSuccessRateDelta += row.successRateDelta
		}
		if summary.instanceCount == 0 {
			continue
		}

		count := float64(summary.instanceCount)
		summary.meanMsaAverageBestDeviation /= count
		summary.meanControlAverageBestDeviation /= count
		summary.meanAverageBestDeviationDelta /= count
		summary.meanMsaSuccessRate /= count
		summary.meanControlSuccessRate /= count
		summary.meanSuccessRateDelta /= count
		summary.signTestPValue = twoSidedSignTestPValue(summary.msaWins, summary.controlWins)
		summaryRows = append(summaryRows, summary)
	}

	return summaryRows
}

type msaDistanceRankedCategoryRow struct {
	instance              string
	dimension             int
	foundOptimalTourCount int
	optimalEdges          int

	bothEdges           int
	msaOnlyEdges        int
	distanceOnlyEdges   int
	neitherEdges        int
	bothOptimal         int
	msaOnlyOptimal      int
	distanceOnlyOptimal int
	neitherOptimal      int
}

func saveMsaDistanceRankedCategoryReport(path string, atspsData []AtspData) error {
	rows, err := buildMsaDistanceRankedCategoryRows(atspsData)
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	var builder strings.Builder
	builder.WriteString("# MSA vs Distance-ranked Edge Categories\n\n")
	builder.WriteString("This report splits directed edges into four categories: edges boosted by at least one rooted MSA-impact matrix and at least one rooted distance-ranked sparse control matrix, edges boosted only by rooted MSA, edges boosted only by rooted distance-ranked sparse control, and edges boosted by neither. The distance-ranked set uses the same per-root boosted-edge count as the rooted MSA heuristic and the same deterministic ordering as the control heuristic.\n\n")
	builder.WriteString("Instances without found optimal tours in `solutions.csv` are omitted, because precision and recall need a reference edge set.\n\n")
	if len(rows) == 0 {
		builder.WriteString("No instances with found optimal tours were available.\n")
		return os.WriteFile(path, []byte(builder.String()), 0644)
	}

	total := totalMsaDistanceRankedCategoryRow(rows)
	writeMsaDistanceRankedCategoryFindings(&builder, total, len(rows))
	builder.WriteString("\n")
	writeMsaDistanceRankedCategoryPooledTable(&builder, total)
	builder.WriteString("\n")
	writeMsaDistanceRankedCategoryInstanceTable(&builder, rows)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildMsaDistanceRankedCategoryRows(atspsData []AtspData) ([]msaDistanceRankedCategoryRow, error) {
	rows := make([]msaDistanceRankedCategoryRow, 0, len(atspsData))
	for _, atspData := range atspsData {
		tours, err := msaHeuristicTours.ReadOptimalTours(atspData.optimalUniqueToursCsvPath)
		if err != nil {
			return nil, fmt.Errorf("%s: read found optimal tours: %w", atspData.name, err)
		}
		optimalEdges := buildAnalysisTourEdgeSet(tours)
		if len(optimalEdges) == 0 {
			continue
		}

		rootedMsas, err := readMsaImpactRootedMsaHeuristics(atspData)
		if err != nil {
			return nil, fmt.Errorf("%s: read rooted MSA impact heuristics: %w", atspData.name, err)
		}

		msaEdges := rootedBoostedModifierEdgeSet(heuristics.BuildRootedMsaHeuristicModifiers(rootedMsas))
		distanceRankedEdges := rootedBoostedModifierEdgeSet(heuristics.BuildRootedDistanceRankedSparseModifiers(atspData.matrix, rootedMsas))
		rows = append(rows, calculateMsaDistanceRankedCategoryRow(atspData.name, len(atspData.matrix), len(tours), optimalEdges, msaEdges, distanceRankedEdges))
	}

	sort.SliceStable(rows, func(i, j int) bool {
		return rows[i].instance < rows[j].instance
	})

	return rows, nil
}

func calculateMsaDistanceRankedCategoryRow(instance string, dimension, foundOptimalTourCount int, optimalEdges, msaEdges, distanceRankedEdges map[models.Edge]struct{}) msaDistanceRankedCategoryRow {
	row := msaDistanceRankedCategoryRow{
		instance:              instance,
		dimension:             dimension,
		foundOptimalTourCount: foundOptimalTourCount,
		optimalEdges:          len(optimalEdges),
	}

	for edge := range msaEdges {
		if _, ok := distanceRankedEdges[edge]; ok {
			row.bothEdges++
		} else {
			row.msaOnlyEdges++
		}
	}
	for edge := range distanceRankedEdges {
		if _, ok := msaEdges[edge]; !ok {
			row.distanceOnlyEdges++
		}
	}

	totalDirectedEdges := dimension * maxIntValue(dimension-1, 0)
	row.neitherEdges = totalDirectedEdges - row.bothEdges - row.msaOnlyEdges - row.distanceOnlyEdges

	for edge := range optimalEdges {
		_, inMsa := msaEdges[edge]
		_, inDistanceRanked := distanceRankedEdges[edge]
		switch {
		case inMsa && inDistanceRanked:
			row.bothOptimal++
		case inMsa:
			row.msaOnlyOptimal++
		case inDistanceRanked:
			row.distanceOnlyOptimal++
		default:
			row.neitherOptimal++
		}
	}

	return row
}

func boostedModifierEdgeSet(modifiers [][]float64) map[models.Edge]struct{} {
	edges := make(map[models.Edge]struct{})
	for from := 0; from < len(modifiers); from++ {
		for to, value := range modifiers[from] {
			if from != to && value > 0.0 {
				edges[models.Edge{From: from, To: to}] = struct{}{}
			}
		}
	}

	return edges
}

func rootedBoostedModifierEdgeSet(rootedModifiers [][][]float64) map[models.Edge]struct{} {
	edges := make(map[models.Edge]struct{})
	for _, modifiers := range rootedModifiers {
		for edge := range boostedModifierEdgeSet(modifiers) {
			edges[edge] = struct{}{}
		}
	}

	return edges
}

func totalMsaDistanceRankedCategoryRow(rows []msaDistanceRankedCategoryRow) msaDistanceRankedCategoryRow {
	total := msaDistanceRankedCategoryRow{instance: "Total"}
	for _, row := range rows {
		total.foundOptimalTourCount += row.foundOptimalTourCount
		total.optimalEdges += row.optimalEdges
		total.bothEdges += row.bothEdges
		total.msaOnlyEdges += row.msaOnlyEdges
		total.distanceOnlyEdges += row.distanceOnlyEdges
		total.neitherEdges += row.neitherEdges
		total.bothOptimal += row.bothOptimal
		total.msaOnlyOptimal += row.msaOnlyOptimal
		total.distanceOnlyOptimal += row.distanceOnlyOptimal
		total.neitherOptimal += row.neitherOptimal
	}

	return total
}

func writeMsaDistanceRankedCategoryFindings(builder *strings.Builder, total msaDistanceRankedCategoryRow, instanceCount int) {
	msaBoostedEdges := total.bothEdges + total.msaOnlyEdges
	sharedMsaEdges := ratio(total.bothEdges, msaBoostedEdges)
	msaOnlyPrecision := ratio(total.msaOnlyOptimal, total.msaOnlyEdges)
	distanceOnlyPrecision := ratio(total.distanceOnlyOptimal, total.distanceOnlyEdges)
	msaOnlyRecall := ratio(total.msaOnlyOptimal, total.optimalEdges)
	distanceOnlyRecall := ratio(total.distanceOnlyOptimal, total.optimalEdges)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **Analyzed %d instances with found optimal tours.**\n", instanceCount)
	fmt.Fprintf(builder, "- **%.2f%% of MSA-impact boosted edges are also distance-ranked sparse edges.**\n", 100*sharedMsaEdges)
	fmt.Fprintf(builder, "- **MSA-only precision %.2f%% vs distance-only precision %.2f%%.**\n", 100*msaOnlyPrecision, 100*distanceOnlyPrecision)
	fmt.Fprintf(builder, "- **MSA-only recall %.2f%% vs distance-only recall %.2f%%.**\n", 100*msaOnlyRecall, 100*distanceOnlyRecall)
	fmt.Fprintf(builder, "- **Found-optimal edge distribution: both %d, MSA-only %d, distance-only %d, neither %d.**\n",
		total.bothOptimal,
		total.msaOnlyOptimal,
		total.distanceOnlyOptimal,
		total.neitherOptimal)
}

func writeMsaDistanceRankedCategoryPooledTable(builder *strings.Builder, total msaDistanceRankedCategoryRow) {
	builder.WriteString("## Pooled Categories\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Category</th><th>Edges</th><th>Found-optimal edges</th><th>Precision [%]</th><th>Recall [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	writeMsaDistanceRankedCategoryPooledRow(builder, "Both", total.bothEdges, total.bothOptimal, total.optimalEdges)
	writeMsaDistanceRankedCategoryPooledRow(builder, "MSA only", total.msaOnlyEdges, total.msaOnlyOptimal, total.optimalEdges)
	writeMsaDistanceRankedCategoryPooledRow(builder, "Distance-ranked only", total.distanceOnlyEdges, total.distanceOnlyOptimal, total.optimalEdges)
	writeMsaDistanceRankedCategoryPooledRow(builder, "Neither", total.neitherEdges, total.neitherOptimal, total.optimalEdges)
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeMsaDistanceRankedCategoryPooledRow(builder *strings.Builder, category string, edgeCount, optimalEdgeCount, totalOptimalEdges int) {
	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
		html.EscapeString(category),
		edgeCount,
		optimalEdgeCount,
		100*ratio(optimalEdgeCount, edgeCount),
		100*ratio(optimalEdgeCount, totalOptimalEdges))
}

func writeMsaDistanceRankedCategoryInstanceTable(builder *strings.Builder, rows []msaDistanceRankedCategoryRow) {
	builder.WriteString("## Per-instance Counts\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>n</th><th>Optimal edges</th><th>Both edges</th><th>MSA-only edges</th><th>Distance-only edges</th><th>Optimal both</th><th>Optimal MSA-only</th><th>Optimal distance-only</th><th>Optimal neither</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td></tr>\n",
			html.EscapeString(row.instance),
			row.dimension,
			row.optimalEdges,
			row.bothEdges,
			row.msaOnlyEdges,
			row.distanceOnlyEdges,
			row.bothOptimal,
			row.msaOnlyOptimal,
			row.distanceOnlyOptimal,
			row.neitherOptimal)
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func saveMsaImpactSummary(path string, atspsData []AtspData) error {
	impactAtspsData := make([]AtspData, 0, len(atspsData))
	for _, atspData := range atspsData {
		impactAtspsData = append(impactAtspsData, withExperimentOutputRoot(atspData, msaImpactResultsDirectoryName))
	}

	rows, err := readFinalResultsSummaryRows(impactAtspsData)
	if err != nil {
		return err
	}
	if len(rows) == 0 {
		return fmt.Errorf("no MSA impact rows to summarize")
	}

	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	var baselineDeviationTotal, baselineSuccessTotal float64
	variantSummaries := map[string]*msaImpactVariantSummary{
		heuristicStrictMsa: newMsaImpactVariantSummary(),
		heuristicRootedMsa: newMsaImpactVariantSummary(),
	}
	variants := []string{heuristicStrictMsa, heuristicRootedMsa}
	const epsilon = 1e-9

	for _, row := range rows {
		baselineMetric, baselineOk := row.metrics[heuristicBaseline]
		if !baselineOk {
			return fmt.Errorf("%s: missing baseline metric", row.instance)
		}

		baselineDeviationTotal += baselineMetric.averageMinDeviation
		baselineSuccessTotal += baselineMetric.successRate

		for _, heuristic := range variants {
			metric, ok := row.metrics[heuristic]
			if !ok {
				return fmt.Errorf("%s: missing %s metric", row.instance, heuristic)
			}

			summary := variantSummaries[heuristic]
			summary.deviationTotal += metric.averageMinDeviation
			summary.successTotal += metric.successRate
			summary.bestWeightCounts[metric.heuristicWeight]++

			delta := metric.averageMinDeviation - baselineMetric.averageMinDeviation
			if math.Abs(delta) < epsilon {
				summary.ties++
			} else if delta < 0 {
				summary.wins++
			} else {
				summary.losses++
			}
		}
	}

	instanceCount := len(rows)
	averageBaselineDeviation := baselineDeviationTotal / float64(instanceCount)
	averageBaselineSuccess := baselineSuccessTotal / float64(instanceCount)

	var builder strings.Builder
	builder.WriteString("# MSA Impact Summary\n\n")
	builder.WriteString("Negative deviation delta means the MSA variant had lower average best deviation than the baseline. Strict MSA uses only edges present in every rooted MSA. Rooted MSA gives each ant the MSA rooted at its start vertex. Both variants use the best heuristic weight from the 0.10-1.00 sweep.\n\n")
	builder.WriteString("## Findings\n\n")
	for _, heuristic := range variants {
		summary := variantSummaries[heuristic]
		averageDeviation := summary.deviationTotal / float64(instanceCount)
		averageSuccess := summary.successTotal / float64(instanceCount)
		fmt.Fprintf(&builder, "- **%s improved average best deviation in %d/%d instances, tied in %d, and was worse in %d.**\n",
			heuristicDisplayName(heuristic), summary.wins, instanceCount, summary.ties, summary.losses)
		fmt.Fprintf(&builder, "- **Mean average best deviation: baseline %.2f%%, %s %.2f%%, delta %s pp; success-rate delta %s pp.**\n",
			averageBaselineDeviation,
			heuristicDisplayName(heuristic),
			averageDeviation,
			formatSignedFloat(averageDeviation-averageBaselineDeviation),
			formatSignedFloat(averageSuccess-averageBaselineSuccess))
		fmt.Fprintf(&builder, "- **Best %s weights: %s.**\n",
			heuristicDisplayName(heuristic),
			formatMsaImpactBestWeightCounts(summary.bestWeightCounts))
	}
	builder.WriteString("\n")

	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>Baseline avg best dev. [%]</th><th>Strict MSA avg best dev. [%]</th><th>Strict weight</th><th>Strict delta [pp]</th><th>Rooted MSA avg best dev. [%]</th><th>Rooted weight</th><th>Rooted delta [pp]</th><th>Baseline success [%]</th><th>Strict success [%]</th><th>Rooted success [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		baselineMetric := row.metrics[heuristicBaseline]
		strictMetric := row.metrics[heuristicStrictMsa]
		rootedMetric := row.metrics[heuristicRootedMsa]
		strictDelta := strictMetric.averageMinDeviation - baselineMetric.averageMinDeviation
		rootedDelta := rootedMetric.averageMinDeviation - baselineMetric.averageMinDeviation
		fmt.Fprintf(&builder,
			"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(row.instance),
			baselineMetric.averageMinDeviation,
			strictMetric.averageMinDeviation,
			strictMetric.heuristicWeight,
			formatSignedFloat(strictDelta),
			rootedMetric.averageMinDeviation,
			rootedMetric.heuristicWeight,
			formatSignedFloat(rootedDelta),
			baselineMetric.successRate,
			strictMetric.successRate,
			rootedMetric.successRate)
	}
	strictSummary := variantSummaries[heuristicStrictMsa]
	rootedSummary := variantSummaries[heuristicRootedMsa]
	averageStrictDeviation := strictSummary.deviationTotal / float64(instanceCount)
	averageRootedDeviation := rootedSummary.deviationTotal / float64(instanceCount)
	averageStrictSuccess := strictSummary.successTotal / float64(instanceCount)
	averageRootedSuccess := rootedSummary.successTotal / float64(instanceCount)
	fmt.Fprintf(&builder,
		"<tr><td><strong>Average</strong></td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td></td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td></td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
		averageBaselineDeviation,
		averageStrictDeviation,
		formatSignedFloat(averageStrictDeviation-averageBaselineDeviation),
		averageRootedDeviation,
		formatSignedFloat(averageRootedDeviation-averageBaselineDeviation),
		averageBaselineSuccess,
		averageStrictSuccess,
		averageRootedSuccess)
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

type msaImpactVariantSummary struct {
	deviationTotal     float64
	successTotal       float64
	wins, ties, losses int
	bestWeightCounts   map[float64]int
}

func newMsaImpactVariantSummary() *msaImpactVariantSummary {
	return &msaImpactVariantSummary{
		bestWeightCounts: make(map[float64]int),
	}
}

func formatMsaImpactBestWeightCounts(counts map[float64]int) string {
	weights := make([]float64, 0, len(counts))
	for weight := range counts {
		weights = append(weights, weight)
	}
	sort.Float64s(weights)

	parts := make([]string, 0, len(weights))
	for _, weight := range weights {
		parts = append(parts, fmt.Sprintf("%.2f (%d)", weight, counts[weight]))
	}

	return strings.Join(parts, ", ")
}

type msaImpactStructureRow struct {
	instance                  string
	dimension                 int
	boostedEdges              float64
	boostedEdgesPerTreeEdge   float64
	missingOutgoingVertices   float64
	missingIncomingVertices   float64
	averageBestDeviationDelta float64
	successRateDelta          float64
}

func saveMsaImpactStructureReport(path string, atspsData []AtspData) error {
	rows, err := readMsaImpactStructureRows(atspsData)
	if err != nil {
		return err
	}
	if len(rows) == 0 {
		return fmt.Errorf("no MSA impact structure rows to summarize")
	}

	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	var boostedRatioTotal, missingOutgoingTotal, missingIncomingTotal float64
	var improvedMissingTotal, nonImprovedMissingTotal float64
	var improvedCount, nonImprovedCount int
	const epsilon = 1e-9

	for _, row := range rows {
		boostedRatioTotal += row.boostedEdgesPerTreeEdge
		missingOutgoingTotal += row.missingOutgoingVertices
		missingIncomingTotal += row.missingIncomingVertices
		missingEndpoints := row.missingOutgoingVertices + row.missingIncomingVertices

		if row.averageBestDeviationDelta < -epsilon {
			improvedMissingTotal += missingEndpoints
			improvedCount++
		} else {
			nonImprovedMissingTotal += missingEndpoints
			nonImprovedCount++
		}
	}

	count := float64(len(rows))
	averageBoostedRatio := boostedRatioTotal / count
	averageMissingOutgoing := missingOutgoingTotal / count
	averageMissingIncoming := missingIncomingTotal / count
	averageImprovedMissing := safeAverage(improvedMissingTotal, improvedCount)
	averageNonImprovedMissing := safeAverage(nonImprovedMissingTotal, nonImprovedCount)

	var builder strings.Builder
	builder.WriteString("# MSA Impact Structure\n\n")
	builder.WriteString("This report describes the average rooted MSA modifier matrix used by the MSA impact pipeline. Each ant uses the MSA rooted at its start vertex. Negative deviation delta means rooted MSA had lower average best deviation than the baseline.\n\n")
	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(&builder, "- **Average boosted-edge ratio against n-1: %.2f.**\n", averageBoostedRatio)
	fmt.Fprintf(&builder, "- **Average missing outgoing vertices: %.2f; average missing incoming vertices: %.2f.**\n", averageMissingOutgoing, averageMissingIncoming)
	if improvedCount != 0 || nonImprovedCount != 0 {
		fmt.Fprintf(&builder, "- **Average missing endpoints: improved cases %.2f, tied/worse cases %.2f.**\n", averageImprovedMissing, averageNonImprovedMissing)
	}
	builder.WriteString("\n")

	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>n</th><th>Avg boosted edges</th><th>Boosted/(n-1)</th><th>Avg missing outgoing</th><th>Avg missing incoming</th><th>Dev. delta [pp]</th><th>Success delta [pp]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		fmt.Fprintf(&builder,
			"<tr><td>%s</td><td align=\"right\">%d</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
			html.EscapeString(row.instance),
			row.dimension,
			row.boostedEdges,
			row.boostedEdgesPerTreeEdge,
			row.missingOutgoingVertices,
			row.missingIncomingVertices,
			formatSignedFloat(row.averageBestDeviationDelta),
			formatSignedFloat(row.successRateDelta))
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func readMsaImpactStructureRows(atspsData []AtspData) ([]msaImpactStructureRow, error) {
	rows := make([]msaImpactStructureRow, 0, len(atspsData))
	for _, atspData := range atspsData {
		impactAtspData := withExperimentOutputRoot(atspData, msaImpactResultsDirectoryName)
		metrics, err := readFinalResultSummaryMetrics(impactAtspData.resultFilePath)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read MSA impact result metrics: %w", atspData.name, err)
		}

		baselineMetric, baselineOk := metrics[heuristicBaseline]
		if !baselineOk {
			return nil, fmt.Errorf("%s: missing baseline metric", atspData.name)
		}

		rootedMsas, err := readMsaImpactRootedMsaHeuristics(atspData)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read rooted MSA impact heuristics: %w", atspData.name, err)
		}

		rootedModifiers := heuristics.BuildRootedMsaHeuristicModifiers(rootedMsas)
		boostedEdges, missingOutgoingVertices, missingIncomingVertices := rootedMsaModifierStructure(rootedModifiers)
		msaMetric, msaOk := metrics[heuristicRootedMsa]
		if !msaOk {
			return nil, fmt.Errorf("%s: missing rooted MSA heuristic metric", atspData.name)
		}

		rows = append(rows, msaImpactStructureRow{
			instance:                  atspData.name,
			dimension:                 len(atspData.matrix),
			boostedEdges:              boostedEdges,
			boostedEdgesPerTreeEdge:   boostedEdgesPerTreeEdge(boostedEdges, len(atspData.matrix)),
			missingOutgoingVertices:   missingOutgoingVertices,
			missingIncomingVertices:   missingIncomingVertices,
			averageBestDeviationDelta: msaMetric.averageMinDeviation - baselineMetric.averageMinDeviation,
			successRateDelta:          msaMetric.successRate - baselineMetric.successRate,
		})
	}

	return rows, nil
}

func rootedMsaModifierStructure(rootedModifiers [][][]float64) (boostedEdges, missingOutgoingVertices, missingIncomingVertices float64) {
	if len(rootedModifiers) == 0 {
		return 0, 0, 0
	}

	for _, modifiers := range rootedModifiers {
		rootBoostedEdges, rootMissingOutgoingVertices, rootMissingIncomingVertices := msaModifierStructure(modifiers)
		boostedEdges += float64(rootBoostedEdges)
		missingOutgoingVertices += float64(rootMissingOutgoingVertices)
		missingIncomingVertices += float64(rootMissingIncomingVertices)
	}

	count := float64(len(rootedModifiers))
	return boostedEdges / count, missingOutgoingVertices / count, missingIncomingVertices / count
}

func msaModifierStructure(modifiers [][]float64) (boostedEdges, missingOutgoingVertices, missingIncomingVertices int) {
	dimension := len(modifiers)
	hasOutgoing := make([]bool, dimension)
	hasIncoming := make([]bool, dimension)

	for from := 0; from < dimension; from++ {
		for to := 0; to < len(modifiers[from]); to++ {
			if from == to || modifiers[from][to] <= 0.0 {
				continue
			}

			boostedEdges++
			hasOutgoing[from] = true
			if to < dimension {
				hasIncoming[to] = true
			}
		}
	}

	for vertex := 0; vertex < dimension; vertex++ {
		if !hasOutgoing[vertex] {
			missingOutgoingVertices++
		}
		if !hasIncoming[vertex] {
			missingIncomingVertices++
		}
	}

	return boostedEdges, missingOutgoingVertices, missingIncomingVertices
}

func boostedEdgesPerTreeEdge(boostedEdges float64, dimension int) float64 {
	if dimension <= 1 {
		return 0
	}
	return boostedEdges / float64(dimension-1)
}

func safeAverage(total float64, count int) float64 {
	if count == 0 {
		return 0
	}
	return total / float64(count)
}
