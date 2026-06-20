package reports

import (
	"atsp_aco_msa/modules/analysis/structure"
	"atsp_aco_msa/modules/project"
	"fmt"
	"html"
	"math"
	"os"
	"path/filepath"
	"strings"
)

type FinalReportsConfig struct {
	Heuristics                     []string
	BaselineHeuristic              string
	StrictMsaHeuristic             string
	CycleCoverHeuristic            string
	CycleCoverMsaPatchingHeuristic string
	DisplayName                    func(string) string
}

type FinalResultsSummaryRow struct {
	Instance string
	Metrics  map[string]FinalResultSummaryMetric
}

func (config FinalReportsConfig) displayName(heuristic string) string {
	if config.DisplayName == nil {
		return heuristic
	}

	return config.DisplayName(heuristic)
}

func (config FinalReportsConfig) comparisonHeuristics() []string {
	if len(config.Heuristics) <= 1 {
		return nil
	}

	return config.Heuristics[1:]
}

func SaveFinalResultsSummary(atspsData []project.AtspData, summaryPath string, config FinalReportsConfig) error {
	if err := os.MkdirAll(filepath.Dir(summaryPath), 0700); err != nil {
		return err
	}

	rows, err := ReadFinalResultsSummaryRows(atspsData)
	if err != nil {
		return err
	}

	return SaveFinalResultsSummaryRows(rows, summaryPath, config)
}

func ReadFinalResultsSummaryRows(atspsData []project.AtspData) ([]FinalResultsSummaryRow, error) {
	rows := make([]FinalResultsSummaryRow, 0, len(atspsData))
	for _, atspData := range atspsData {
		metrics, err := ReadFinalResultSummaryMetrics(atspData.ResultFilePath)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read final result metrics: %w", atspData.Name, err)
		}

		rows = append(rows, FinalResultsSummaryRow{Instance: atspData.Name, Metrics: metrics})
	}

	return rows, nil
}

func SaveFinalResultsSummaryRows(rows []FinalResultsSummaryRow, summaryPath string, config FinalReportsConfig) error {
	if err := os.MkdirAll(filepath.Dir(summaryPath), 0700); err != nil {
		return err
	}

	averageMetrics := averageFinalResultsSummaryMetrics(rows, config)
	var builder strings.Builder
	builder.WriteString("# Final Results Summary\n\n")
	writeFinalResultsSummaryFindings(&builder, rows, averageMetrics, config)
	builder.WriteString("\n")
	writeFinalResultsSummaryTable(&builder, rows, averageMetrics, config)

	return os.WriteFile(summaryPath, []byte(builder.String()), 0644)
}

func averageFinalResultsSummaryMetrics(rows []FinalResultsSummaryRow, config FinalReportsConfig) map[string]FinalResultSummaryMetric {
	totals := make(map[string]FinalResultSummaryMetric, len(config.Heuristics))
	counts := make(map[string]int, len(config.Heuristics))
	for _, row := range rows {
		for _, heuristic := range config.Heuristics {
			metric, ok := row.Metrics[heuristic]
			if !ok {
				continue
			}

			total := totals[heuristic]
			total.AverageMinDeviation += metric.AverageMinDeviation
			total.SuccessRate += metric.SuccessRate
			total.AverageBestIteration += metric.AverageBestIteration
			total.Iterations += metric.Iterations
			totals[heuristic] = total
			counts[heuristic]++
		}
	}

	averageMetrics := make(map[string]FinalResultSummaryMetric, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		count := counts[heuristic]
		if count == 0 {
			continue
		}

		total := totals[heuristic]
		averageMetrics[heuristic] = FinalResultSummaryMetric{
			AverageMinDeviation:  total.AverageMinDeviation / float64(count),
			SuccessRate:          total.SuccessRate / float64(count),
			AverageBestIteration: total.AverageBestIteration / float64(count),
			Iterations:           int(math.Round(float64(total.Iterations) / float64(count))),
		}
	}

	return averageMetrics
}

func writeFinalResultsSummaryFindings(builder *strings.Builder, rows []FinalResultsSummaryRow, averageMetrics map[string]FinalResultSummaryMetric, config FinalReportsConfig) {
	averageBestDeviationHighlights, averageBestSuccessHighlights := finalResultsSummaryHighlights(averageMetrics, config)
	bestDeviationHeuristics := highlightedHeuristicDisplayNames(averageBestDeviationHighlights, config)
	bestSuccessHeuristics := highlightedHeuristicDisplayNames(averageBestSuccessHighlights, config)
	deviationWinCounts := make(map[string]int, len(config.Heuristics))

	for _, row := range rows {
		bestDeviationHighlights, _ := finalResultsSummaryHighlights(row.Metrics, config)
		for _, heuristic := range config.Heuristics {
			if bestDeviationHighlights[heuristic] {
				deviationWinCounts[heuristic]++
			}
		}
	}

	builder.WriteString("## Findings\n\n")
	if len(bestDeviationHeuristics) != 0 {
		bestDeviation := averageMetrics[bestDeviationHeuristics[0].heuristic].AverageMinDeviation
		fmt.Fprintf(builder, "- **%s has the lowest average best deviation overall: %.2f%%.**\n",
			joinHeuristicDisplayNames(bestDeviationHeuristics),
			bestDeviation)
	}
	if len(bestSuccessHeuristics) != 0 {
		bestSuccess := averageMetrics[bestSuccessHeuristics[0].heuristic].SuccessRate
		fmt.Fprintf(builder, "- **%s has the highest average success rate overall: %.2f%%.**\n",
			joinHeuristicDisplayNames(bestSuccessHeuristics),
			bestSuccess)
	}

	fmt.Fprintf(builder, "- **Best-or-tied average best deviation counts: %s.**\n", heuristicCountList(deviationWinCounts, len(rows), config))
}

func writeFinalResultsSummaryTable(builder *strings.Builder, rows []FinalResultsSummaryRow, averageMetrics map[string]FinalResultSummaryMetric, config FinalReportsConfig) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th rowspan=\"2\">Instance</th>")
	for _, heuristic := range config.Heuristics {
		fmt.Fprintf(builder, "<th colspan=\"2\">%s</th>", html.EscapeString(config.displayName(heuristic)))
	}
	builder.WriteString("</tr>\n")

	builder.WriteString("<tr>")
	for range config.Heuristics {
		builder.WriteString("<th>Avg best dev. [%]</th><th>Success [%]</th>")
	}
	builder.WriteString("</tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, row := range rows {
		writeFinalResultsSummaryTableRow(builder, row.Instance, row.Metrics, false, config)
	}
	writeFinalResultsSummaryTableRow(builder, "Average", averageMetrics, true, config)
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeFinalResultsSummaryTableRow(builder *strings.Builder, instance string, metrics map[string]FinalResultSummaryMetric, boldInstance bool, config FinalReportsConfig) {
	deviationHighlights, successHighlights := finalResultsSummaryHighlights(metrics, config)
	instanceCell := html.EscapeString(instance)
	if boldInstance {
		instanceCell = "<strong>" + instanceCell + "</strong>"
	}

	fmt.Fprintf(builder, "<tr><td>%s</td>", instanceCell)
	for _, heuristic := range config.Heuristics {
		metric, ok := metrics[heuristic]
		if !ok {
			builder.WriteString("<td></td><td></td>")
			continue
		}

		fmt.Fprintf(builder, "<td align=\"right\">%s</td><td align=\"right\">%s</td>",
			finalResultsSummaryMetricCell(metric.AverageMinDeviation, deviationHighlights[heuristic]),
			finalResultsSummaryMetricCell(metric.SuccessRate, successHighlights[heuristic]))
	}
	builder.WriteString("</tr>\n")
}

func finalResultsSummaryHighlights(metrics map[string]FinalResultSummaryMetric, config FinalReportsConfig) (map[string]bool, map[string]bool) {
	deviationHighlights := make(map[string]bool, len(config.Heuristics))
	successHighlights := make(map[string]bool, len(config.Heuristics))
	minDeviation := math.Inf(1)
	maxSuccess := math.Inf(-1)

	for _, heuristic := range config.Heuristics {
		metric, ok := metrics[heuristic]
		if !ok {
			continue
		}

		if metric.AverageMinDeviation < minDeviation {
			minDeviation = metric.AverageMinDeviation
		}
		if metric.SuccessRate > maxSuccess {
			maxSuccess = metric.SuccessRate
		}
	}

	for _, heuristic := range config.Heuristics {
		metric, ok := metrics[heuristic]
		if !ok {
			continue
		}

		if math.Abs(metric.AverageMinDeviation-minDeviation) < 1e-9 {
			deviationHighlights[heuristic] = true
		}
		if maxSuccess > 0 && math.Abs(metric.SuccessRate-maxSuccess) < 1e-9 {
			successHighlights[heuristic] = true
		}
	}

	return deviationHighlights, successHighlights
}

type heuristicDisplay struct {
	heuristic string
	display   string
}

func highlightedHeuristicDisplayNames(highlights map[string]bool, config FinalReportsConfig) []heuristicDisplay {
	names := make([]heuristicDisplay, 0, len(highlights))
	for _, heuristic := range config.Heuristics {
		if highlights[heuristic] {
			names = append(names, heuristicDisplay{
				heuristic: heuristic,
				display:   config.displayName(heuristic),
			})
		}
	}

	return names
}

func joinHeuristicDisplayNames(heuristics []heuristicDisplay) string {
	names := make([]string, len(heuristics))
	for i, heuristic := range heuristics {
		names[i] = heuristic.display
	}

	return strings.Join(names, ", ")
}

func heuristicCountList(counts map[string]int, total int, config FinalReportsConfig) string {
	parts := make([]string, 0, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		parts = append(parts, fmt.Sprintf("%s %d/%d", config.displayName(heuristic), counts[heuristic], total))
	}
	return strings.Join(parts, ", ")
}

func heuristicFloatList(values map[string]float64, format string, config FinalReportsConfig) string {
	parts := make([]string, 0, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		value, ok := values[heuristic]
		if !ok {
			continue
		}
		parts = append(parts, fmt.Sprintf("%s "+format, config.displayName(heuristic), value))
	}
	return strings.Join(parts, ", ")
}

func finalResultsSummaryMetricCell(value float64, bold bool) string {
	return metricCell(value, bold)
}

func SaveFinalPairwisePerformanceReport(path string, rows []FinalResultsSummaryRow, config FinalReportsConfig) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	comparisons := finalPairwisePerformanceComparisons(rows, config)

	var builder strings.Builder
	builder.WriteString("# Pairwise Performance Summary\n\n")
	builder.WriteString("Negative average-best-deviation delta means the first heuristic in the comparison had lower deviation.\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Comparison</th><th>Avg best dev. delta [pp]</th><th>Wins</th><th>Ties</th><th>Losses</th><th>Success delta [pp]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, comparison := range comparisons {
		fmt.Fprintf(&builder,
			"<tr><td>%s vs %s</td><td align=\"right\">%s</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%s</td></tr>\n",
			html.EscapeString(config.displayName(comparison.left)),
			html.EscapeString(config.displayName(comparison.right)),
			formatSignedFloat(comparison.averageBestDeviationDelta),
			comparison.wins,
			comparison.ties,
			comparison.losses,
			formatSignedFloat(comparison.successRateDelta))
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func finalPairwisePerformanceComparisons(rows []FinalResultsSummaryRow, config FinalReportsConfig) []finalPairwisePerformanceComparison {
	comparisons := make([]finalPairwisePerformanceComparison, 0)
	for _, heuristic := range config.comparisonHeuristics() {
		comparison := calculateFinalPairwisePerformanceComparison(rows, heuristic, config.BaselineHeuristic)
		if comparison.count > 0 {
			comparisons = append(comparisons, comparison)
		}
	}
	heuristics := config.comparisonHeuristics()
	for i := 0; i < len(heuristics); i++ {
		for j := i + 1; j < len(heuristics); j++ {
			comparison := calculateFinalPairwisePerformanceComparison(rows, heuristics[i], heuristics[j])
			if comparison.count > 0 {
				comparisons = append(comparisons, comparison)
			}
		}
	}
	return comparisons
}

type finalPairwisePerformanceComparison struct {
	left                      string
	right                     string
	count                     int
	wins                      int
	ties                      int
	losses                    int
	averageBestDeviationDelta float64
	successRateDelta          float64
}

func calculateFinalPairwisePerformanceComparison(rows []FinalResultsSummaryRow, left, right string) finalPairwisePerformanceComparison {
	const epsilon = 1e-9
	comparison := finalPairwisePerformanceComparison{left: left, right: right}

	for _, row := range rows {
		leftMetric, leftOk := row.Metrics[left]
		rightMetric, rightOk := row.Metrics[right]
		if !leftOk || !rightOk {
			continue
		}

		comparison.count++
		deviationDelta := leftMetric.AverageMinDeviation - rightMetric.AverageMinDeviation
		comparison.averageBestDeviationDelta += deviationDelta
		comparison.successRateDelta += leftMetric.SuccessRate - rightMetric.SuccessRate

		if math.Abs(deviationDelta) < epsilon {
			comparison.ties++
		} else if deviationDelta < 0 {
			comparison.wins++
		} else {
			comparison.losses++
		}
	}

	if comparison.count > 0 {
		comparison.averageBestDeviationDelta /= float64(comparison.count)
		comparison.successRateDelta /= float64(comparison.count)
	}

	return comparison
}

func SaveFinalConvergenceSummaryReport(path string, rows []FinalResultsSummaryRow, config FinalReportsConfig) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	averages := averageConvergenceByHeuristic(rows, nil, config)

	var builder strings.Builder
	builder.WriteString("# Convergence Summary\n\n")
	builder.WriteString("Each value is the average iteration where the best run solution was found, expressed as a percentage of the configured iteration budget. Lower values mean earlier convergence.\n\n")
	writeFinalConvergenceFindings(&builder, rows, averages, config)
	builder.WriteString("\n")
	writeFinalConvergenceTable(&builder, rows, averages, config)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeFinalConvergenceFindings(builder *strings.Builder, rows []FinalResultsSummaryRow, averages map[string]float64, config FinalReportsConfig) {
	winCounts := make(map[string]int, len(config.Heuristics))
	for _, row := range rows {
		highlights := lowestConvergenceHighlights(row.Metrics, config)
		for _, heuristic := range config.Heuristics {
			if highlights[heuristic] {
				winCounts[heuristic]++
			}
		}
	}

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder,
		"- **Average best-iteration position: %s.**\n",
		heuristicFloatList(averages, "%.2f%%", config))
	fmt.Fprintf(builder,
		"- **Earliest-or-tied convergence counts: %s.**\n",
		heuristicCountList(winCounts, len(rows), config))
}

func writeFinalConvergenceTable(builder *strings.Builder, rows []FinalResultsSummaryRow, averages map[string]float64, config FinalReportsConfig) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th>")
	for _, heuristic := range config.Heuristics {
		fmt.Fprintf(builder, "<th>%s best iter [%%]</th>", html.EscapeString(config.displayName(heuristic)))
	}
	builder.WriteString("</tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		highlights := lowestConvergenceHighlights(row.Metrics, config)
		fmt.Fprintf(builder, "<tr><td>%s</td>", html.EscapeString(row.Instance))
		for _, heuristic := range config.Heuristics {
			metric, ok := row.Metrics[heuristic]
			if !ok {
				builder.WriteString("<td></td>")
				continue
			}
			fmt.Fprintf(builder, "<td align=\"right\">%s</td>",
				finalResultsSummaryMetricCell(convergencePercent(metric), highlights[heuristic]))
		}
		builder.WriteString("</tr>\n")
	}

	averageHighlights := lowestFloatHighlights(averages, false)
	builder.WriteString("<tr><td><strong>Average</strong></td>")
	for _, heuristic := range config.Heuristics {
		average, ok := averages[heuristic]
		if !ok {
			builder.WriteString("<td></td>")
			continue
		}
		fmt.Fprintf(builder, "<td align=\"right\">%s</td>",
			finalResultsSummaryMetricCell(average, averageHighlights[heuristic]))
	}
	builder.WriteString("</tr>\n")
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func averageConvergenceByHeuristic(rows []FinalResultsSummaryRow, allowedInstances map[string]struct{}, config FinalReportsConfig) map[string]float64 {
	sums := make(map[string]float64, len(config.Heuristics))
	counts := make(map[string]int, len(config.Heuristics))
	for _, row := range rows {
		if allowedInstances != nil {
			if _, ok := allowedInstances[row.Instance]; !ok {
				continue
			}
		}

		for _, heuristic := range config.Heuristics {
			metric, ok := row.Metrics[heuristic]
			if !ok {
				continue
			}
			if metric.Iterations <= 0 {
				continue
			}

			sums[heuristic] += convergencePercent(metric)
			counts[heuristic]++
		}
	}

	averages := make(map[string]float64, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		if counts[heuristic] == 0 {
			continue
		}
		averages[heuristic] = sums[heuristic] / float64(counts[heuristic])
	}

	return averages
}

func lowestConvergenceHighlights(metrics map[string]FinalResultSummaryMetric, config FinalReportsConfig) map[string]bool {
	values := make(map[string]float64, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		metric, ok := metrics[heuristic]
		if !ok || metric.Iterations <= 0 {
			continue
		}
		values[heuristic] = convergencePercent(metric)
	}

	return lowestFloatHighlights(values, false)
}

func convergencePercent(metric FinalResultSummaryMetric) float64 {
	if metric.Iterations <= 0 {
		return 0
	}

	return 100.0 * metric.AverageBestIteration / float64(metric.Iterations)
}

func SaveFinalThreeOptComparisonReport(path string, finalRows, finalThreeOptRows []FinalResultsSummaryRow, config FinalReportsConfig) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	finalAverages := averageFinalResultsSummaryMetrics(finalRows, config)
	finalThreeOptAverages := averageFinalResultsSummaryMetrics(finalThreeOptRows, config)

	var builder strings.Builder
	builder.WriteString("# Reduced 3-Opt Impact\n\n")
	builder.WriteString("This report compares the final MMAS experiments without local search against the final MMAS experiments with reduced 3-opt enabled.\n\n")
	writeFinalThreeOptComparisonFindings(&builder, finalAverages, finalThreeOptAverages, config)
	builder.WriteString("\n")
	writeFinalThreeOptImpactTable(&builder, finalAverages, finalThreeOptAverages, config)
	builder.WriteString("\n")
	writeFinalThreeOptSignalTable(&builder, finalAverages, finalThreeOptAverages, config)
	builder.WriteString("\n")
	builder.WriteString("Deviation gain vs baseline is `baseline average best deviation - heuristic average best deviation`, so positive values mean that the heuristic improved over the baseline. Signal remaining is the share of this deviation gain still visible after enabling reduced 3-opt.\n")

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeFinalThreeOptComparisonFindings(builder *strings.Builder, finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, config FinalReportsConfig) {
	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder,
		"- **Reduced 3-opt average-best-deviation deltas: %s.**\n",
		finalThreeOptDeviationDeltaList(finalAverages, finalThreeOptAverages, config))
	fmt.Fprintf(builder,
		"- **Reduced 3-opt success-rate deltas: %s.**\n",
		finalThreeOptSuccessDeltaList(finalAverages, finalThreeOptAverages, config))
	fmt.Fprintf(builder,
		"- **Deviation gain over baseline with and without 3-opt: %s.**\n",
		finalThreeOptGainList(finalAverages, finalThreeOptAverages, config))
	fmt.Fprintf(builder, "- **Signal remaining after enabling 3-opt: %s.**\n", finalThreeOptSignalRemainingList(finalAverages, finalThreeOptAverages, config))
	builder.WriteString("- **This supports treating reduced 3-opt as a strong local-search layer that partially hides the construction heuristic effect.**\n")
}

func writeFinalThreeOptImpactTable(builder *strings.Builder, finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, config FinalReportsConfig) {
	builder.WriteString("## Overall Effect\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, heuristic := range config.Heuristics {
		without, withoutOk := finalAverages[heuristic]
		with, withOk := finalThreeOptAverages[heuristic]
		if !withoutOk || !withOk {
			continue
		}
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(config.displayName(heuristic)),
			without.AverageMinDeviation,
			with.AverageMinDeviation,
			without.SuccessRate,
			with.SuccessRate)
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeFinalThreeOptSignalTable(builder *strings.Builder, finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, config FinalReportsConfig) {
	builder.WriteString("## Heuristic Signal\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Heuristic</th><th>Dev. gain vs baseline without 3-opt [pp]</th><th>Dev. gain vs baseline with 3-opt [pp]</th><th>Signal remaining [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, heuristic := range config.comparisonHeuristics() {
		if !hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages, heuristic, config) {
			continue
		}

		gainWithout := deviationGainVsBaseline(finalAverages, heuristic, config)
		gainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristic, config)

		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(config.displayName(heuristic)),
			formatSignedFloat(gainWithout),
			formatSignedFloat(gainWith),
			signalRemainingPercent(gainWithout, gainWith))
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func finalThreeOptDeviationDeltaList(finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, config FinalReportsConfig) string {
	parts := make([]string, 0, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		if !hasPairedMetrics(finalAverages, finalThreeOptAverages, heuristic) {
			continue
		}
		delta := finalAverages[heuristic].AverageMinDeviation - finalThreeOptAverages[heuristic].AverageMinDeviation
		parts = append(parts, fmt.Sprintf("%s %s pp", config.displayName(heuristic), formatSignedFloat(delta)))
	}
	return strings.Join(parts, ", ")
}

func finalThreeOptSuccessDeltaList(finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, config FinalReportsConfig) string {
	parts := make([]string, 0, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		if !hasPairedMetrics(finalAverages, finalThreeOptAverages, heuristic) {
			continue
		}
		delta := finalThreeOptAverages[heuristic].SuccessRate - finalAverages[heuristic].SuccessRate
		parts = append(parts, fmt.Sprintf("%s %s pp", config.displayName(heuristic), formatSignedFloat(delta)))
	}
	return strings.Join(parts, ", ")
}

func finalThreeOptGainList(finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, config FinalReportsConfig) string {
	parts := make([]string, 0, len(config.comparisonHeuristics()))
	for _, heuristic := range config.comparisonHeuristics() {
		if !hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages, heuristic, config) {
			continue
		}
		gainWithout := deviationGainVsBaseline(finalAverages, heuristic, config)
		gainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristic, config)
		parts = append(parts, fmt.Sprintf("%s %s -> %s pp", config.displayName(heuristic), formatSignedFloat(gainWithout), formatSignedFloat(gainWith)))
	}
	return strings.Join(parts, ", ")
}

func finalThreeOptSignalRemainingList(finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, config FinalReportsConfig) string {
	parts := make([]string, 0, len(config.comparisonHeuristics()))
	for _, heuristic := range config.comparisonHeuristics() {
		if !hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages, heuristic, config) {
			continue
		}
		gainWithout := deviationGainVsBaseline(finalAverages, heuristic, config)
		gainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristic, config)
		parts = append(parts, fmt.Sprintf("%s %.2f%%", config.displayName(heuristic), signalRemainingPercent(gainWithout, gainWith)))
	}
	return strings.Join(parts, ", ")
}

func hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages map[string]FinalResultSummaryMetric, heuristic string, config FinalReportsConfig) bool {
	return hasPairedMetrics(finalAverages, finalThreeOptAverages, heuristic) &&
		hasPairedMetrics(finalAverages, finalThreeOptAverages, config.BaselineHeuristic)
}

func hasPairedMetrics(left, right map[string]FinalResultSummaryMetric, heuristic string) bool {
	_, leftOk := left[heuristic]
	_, rightOk := right[heuristic]
	return leftOk && rightOk
}

func deviationGainVsBaseline(averages map[string]FinalResultSummaryMetric, heuristic string, config FinalReportsConfig) float64 {
	return averages[config.BaselineHeuristic].AverageMinDeviation - averages[heuristic].AverageMinDeviation
}

func signalRemainingPercent(without, with float64) float64 {
	if math.Abs(without) < 1e-9 {
		return 0
	}

	return 100.0 * math.Abs(with) / math.Abs(without)
}

func SaveStructuralPerformanceLinkReport(path string, rows []FinalResultsSummaryRow, analyses []structure.InstanceAnalysis, config FinalReportsConfig) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	structuralRows := filterAnalysesWithFoundOptimalEdges(sortedStructuralAnalyses(analyses))
	structuralTotals := structuralSimilarityTotals(structuralRows)
	allowedInstances := make(map[string]struct{}, len(structuralRows))
	for _, row := range structuralRows {
		allowedInstances[row.Instance] = struct{}{}
	}

	performance := averagePerformanceByHeuristic(rows, allowedInstances, config)
	msaPrecision := ratio(structuralTotals.msaOptimalEdges, structuralTotals.msaEdges)
	cycleCoverPrecision := ratio(structuralTotals.cycleCoverOptimalEdges, structuralTotals.cycleCoverEdges)
	patchingPrecision := ratio(structuralTotals.patchingOptimalEdges, structuralTotals.patchingEdges)
	msaRecall := ratio(structuralTotals.msaOptimalEdges, structuralTotals.foundOptimalEdges)
	cycleCoverRecall := ratio(structuralTotals.cycleCoverOptimalEdges, structuralTotals.foundOptimalEdges)
	patchingRecall := ratio(structuralTotals.patchingOptimalEdges, structuralTotals.foundOptimalEdges)

	structuralPrecision := map[string]float64{
		config.StrictMsaHeuristic:             msaPrecision,
		config.CycleCoverHeuristic:            cycleCoverPrecision,
		config.CycleCoverMsaPatchingHeuristic: patchingPrecision,
	}
	structuralRecall := map[string]float64{
		config.StrictMsaHeuristic:             msaRecall,
		config.CycleCoverHeuristic:            cycleCoverRecall,
		config.CycleCoverMsaPatchingHeuristic: patchingRecall,
	}

	includedHeuristics := make([]string, 0, len(config.comparisonHeuristics()))
	precisionValues := make(map[string]float64)
	recallValues := make(map[string]float64)
	deviationValues := make(map[string]float64)
	successValues := make(map[string]float64)
	for _, heuristic := range config.comparisonHeuristics() {
		metric, ok := performance[heuristic]
		if !ok {
			continue
		}
		precision, hasStructuralMetrics := structuralPrecision[heuristic]
		if !hasStructuralMetrics {
			continue
		}

		includedHeuristics = append(includedHeuristics, heuristic)
		precisionValues[heuristic] = precision
		recallValues[heuristic] = structuralRecall[heuristic]
		deviationValues[heuristic] = metric.AverageMinDeviation
		successValues[heuristic] = metric.SuccessRate
	}

	precisionHighlights := highestFloatHighlights(precisionValues)
	recallHighlights := highestFloatHighlights(recallValues)
	deviationHighlights := lowestFloatHighlights(deviationValues, true)
	successHighlights := highestFloatHighlights(successValues)

	var builder strings.Builder
	builder.WriteString("# Structural Similarity And Performance\n\n")
	builder.WriteString("This table links structural similarity to found optimal tours with final MMAS performance. Both structural and performance values are computed only for instances with at least one found optimal tour.\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Heuristic</th><th>Structural precision [%]</th><th>Structural recall [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, heuristic := range includedHeuristics {
		writeStructuralPerformanceLinkRow(&builder, heuristic, structuralPrecision[heuristic], structuralRecall[heuristic], performance[heuristic], precisionHighlights, recallHighlights, deviationHighlights, successHighlights, config)
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeStructuralPerformanceLinkRow(builder *strings.Builder, heuristic string, precision, recall float64, performance FinalResultSummaryMetric, precisionHighlights, recallHighlights, deviationHighlights, successHighlights map[string]bool, config FinalReportsConfig) {
	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		html.EscapeString(config.displayName(heuristic)),
		finalResultsSummaryMetricCell(100*precision, precisionHighlights[heuristic]),
		finalResultsSummaryMetricCell(100*recall, recallHighlights[heuristic]),
		finalResultsSummaryMetricCell(performance.AverageMinDeviation, deviationHighlights[heuristic]),
		finalResultsSummaryMetricCell(performance.SuccessRate, successHighlights[heuristic]))
}

func averagePerformanceByHeuristic(rows []FinalResultsSummaryRow, allowedInstances map[string]struct{}, config FinalReportsConfig) map[string]FinalResultSummaryMetric {
	totals := make(map[string]FinalResultSummaryMetric, len(config.Heuristics))
	counts := make(map[string]int, len(config.Heuristics))
	for _, row := range rows {
		if allowedInstances != nil {
			if _, ok := allowedInstances[row.Instance]; !ok {
				continue
			}
		}

		for _, heuristic := range config.Heuristics {
			metric, ok := row.Metrics[heuristic]
			if !ok {
				continue
			}

			total := totals[heuristic]
			total.AverageMinDeviation += metric.AverageMinDeviation
			total.SuccessRate += metric.SuccessRate
			totals[heuristic] = total
			counts[heuristic]++
		}
	}

	averages := make(map[string]FinalResultSummaryMetric, len(config.Heuristics))
	for _, heuristic := range config.Heuristics {
		count := counts[heuristic]
		if count == 0 {
			continue
		}

		total := totals[heuristic]
		averages[heuristic] = FinalResultSummaryMetric{
			AverageMinDeviation: total.AverageMinDeviation / float64(count),
			SuccessRate:         total.SuccessRate / float64(count),
		}
	}

	return averages
}

func highestFloatHighlights(values map[string]float64) map[string]bool {
	return floatHighlights(values, true, false)
}

func lowestFloatHighlights(values map[string]float64, allowZero bool) map[string]bool {
	return floatHighlights(values, false, allowZero)
}

func floatHighlights(values map[string]float64, higherIsBetter, allowZero bool) map[string]bool {
	const epsilon = 1e-9
	highlights := make(map[string]bool, len(values))
	if len(values) == 0 {
		return highlights
	}

	best := math.Inf(1)
	if higherIsBetter {
		best = math.Inf(-1)
	}

	for _, value := range values {
		if !allowZero && value <= 0 {
			continue
		}
		if higherIsBetter {
			if value > best {
				best = value
			}
		} else if value < best {
			best = value
		}
	}

	if math.IsInf(best, 0) {
		return highlights
	}

	for key, value := range values {
		if !allowZero && value <= 0 {
			continue
		}
		if math.Abs(value-best) < epsilon {
			highlights[key] = true
		}
	}

	return highlights
}

func formatSignedFloat(value float64) string {
	return FormatSignedFloat(value)
}
