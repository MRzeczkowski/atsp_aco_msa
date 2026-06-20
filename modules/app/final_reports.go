package app

import (
	"atsp_aco_msa/modules/analysis/reports"
	"atsp_aco_msa/modules/analysis/structure"
)

func finalReportsConfig() reports.FinalReportsConfig {
	return reports.FinalReportsConfig{
		Heuristics:                     finalResultsSummaryHeuristics,
		BaselineHeuristic:              heuristicBaseline,
		StrictMsaHeuristic:             heuristicStrictMsa,
		CycleCoverHeuristic:            heuristicCycleCover,
		CycleCoverMsaPatchingHeuristic: heuristicCycleCoverMsaPatching,
		DisplayName:                    heuristicDisplayName,
	}
}

func saveFinalResultsSummary(atspsData []AtspData, summaryPath string) error {
	return reports.SaveFinalResultsSummary(atspsData, summaryPath, finalReportsConfig())
}

func readFinalResultsSummaryRows(atspsData []AtspData) ([]finalResultsSummaryRow, error) {
	return reports.ReadFinalResultsSummaryRows(atspsData)
}

func saveFinalResultsSummaryRows(rows []finalResultsSummaryRow, summaryPath string) error {
	return reports.SaveFinalResultsSummaryRows(rows, summaryPath, finalReportsConfig())
}

func saveFinalPairwisePerformanceReport(path string, rows []finalResultsSummaryRow) error {
	return reports.SaveFinalPairwisePerformanceReport(path, rows, finalReportsConfig())
}

func saveFinalConvergenceSummaryReport(path string, rows []finalResultsSummaryRow) error {
	return reports.SaveFinalConvergenceSummaryReport(path, rows, finalReportsConfig())
}

func saveFinalThreeOptComparisonReport(path string, finalRows, finalThreeOptRows []finalResultsSummaryRow) error {
	return reports.SaveFinalThreeOptComparisonReport(path, finalRows, finalThreeOptRows, finalReportsConfig())
}

func saveStructuralPerformanceLinkReport(path string, rows []finalResultsSummaryRow, analyses []structure.InstanceAnalysis) error {
	return reports.SaveStructuralPerformanceLinkReport(path, rows, analyses, finalReportsConfig())
}

func finalResultsSummaryMetricCell(value float64, bold bool) string {
	return reports.MetricCell(value, bold)
}

func formatSignedFloat(value float64) string {
	return reports.FormatSignedFloat(value)
}
