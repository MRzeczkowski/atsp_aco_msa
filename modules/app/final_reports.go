package app

import (
	"atsp_aco_msa/modules/analysis/structure"
	"fmt"
	"html"
	"math"
	"os"
	"path/filepath"
	"strings"
)

func saveFinalResultsSummary(atspsData []AtspData, summaryPath string) error {
	if err := os.MkdirAll(filepath.Dir(summaryPath), 0700); err != nil {
		return err
	}

	rows, err := readFinalResultsSummaryRows(atspsData)
	if err != nil {
		return err
	}

	return saveFinalResultsSummaryRows(rows, summaryPath)
}

func readFinalResultsSummaryRows(atspsData []AtspData) ([]finalResultsSummaryRow, error) {
	rows := make([]finalResultsSummaryRow, 0, len(atspsData))
	for _, atspData := range atspsData {
		metrics, err := readFinalResultSummaryMetrics(atspData.ResultFilePath)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read final result metrics: %w", atspData.Name, err)
		}

		rows = append(rows, finalResultsSummaryRow{instance: atspData.Name, metrics: metrics})
	}

	return rows, nil
}

func saveFinalResultsSummaryRows(rows []finalResultsSummaryRow, summaryPath string) error {
	if err := os.MkdirAll(filepath.Dir(summaryPath), 0700); err != nil {
		return err
	}

	averageMetrics := averageFinalResultsSummaryMetrics(rows)
	var builder strings.Builder
	builder.WriteString("# Final Results Summary\n\n")
	writeFinalResultsSummaryFindings(&builder, rows, averageMetrics)
	builder.WriteString("\n")
	writeFinalResultsSummaryTable(&builder, rows, averageMetrics)

	return os.WriteFile(summaryPath, []byte(builder.String()), 0644)
}

func averageFinalResultsSummaryMetrics(rows []finalResultsSummaryRow) map[string]finalResultsSummaryMetric {
	totals := make(map[string]finalResultsSummaryMetric, len(finalResultsSummaryHeuristics))
	counts := make(map[string]int, len(finalResultsSummaryHeuristics))
	for _, row := range rows {
		for _, heuristic := range finalResultsSummaryHeuristics {
			metric, ok := row.metrics[heuristic]
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

	averageMetrics := make(map[string]finalResultsSummaryMetric, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		count := counts[heuristic]
		if count == 0 {
			continue
		}

		total := totals[heuristic]
		averageMetrics[heuristic] = finalResultsSummaryMetric{
			AverageMinDeviation:  total.AverageMinDeviation / float64(count),
			SuccessRate:          total.SuccessRate / float64(count),
			AverageBestIteration: total.AverageBestIteration / float64(count),
			Iterations:           int(math.Round(float64(total.Iterations) / float64(count))),
		}
	}

	return averageMetrics
}

func writeFinalResultsSummaryFindings(builder *strings.Builder, rows []finalResultsSummaryRow, averageMetrics map[string]finalResultsSummaryMetric) {
	averageBestDeviationHighlights, averageBestSuccessHighlights := finalResultsSummaryHighlights(averageMetrics)
	bestDeviationHeuristics := highlightedHeuristicDisplayNames(averageBestDeviationHighlights)
	bestSuccessHeuristics := highlightedHeuristicDisplayNames(averageBestSuccessHighlights)
	deviationWinCounts := make(map[string]int, len(finalResultsSummaryHeuristics))

	for _, row := range rows {
		bestDeviationHighlights, _ := finalResultsSummaryHighlights(row.metrics)
		for _, heuristic := range finalResultsSummaryHeuristics {
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

	fmt.Fprintf(builder, "- **Best-or-tied average best deviation counts: %s.**\n", heuristicCountList(deviationWinCounts, len(rows)))
}

func writeFinalResultsSummaryTable(builder *strings.Builder, rows []finalResultsSummaryRow, averageMetrics map[string]finalResultsSummaryMetric) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th rowspan=\"2\">Instance</th>")
	for _, heuristic := range finalResultsSummaryHeuristics {
		fmt.Fprintf(builder, "<th colspan=\"2\">%s</th>", html.EscapeString(heuristicDisplayName(heuristic)))
	}
	builder.WriteString("</tr>\n")

	builder.WriteString("<tr>")
	for range finalResultsSummaryHeuristics {
		builder.WriteString("<th>Avg best dev. [%]</th><th>Success [%]</th>")
	}
	builder.WriteString("</tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, row := range rows {
		writeFinalResultsSummaryTableRow(builder, row.instance, row.metrics, false)
	}
	writeFinalResultsSummaryTableRow(builder, "Average", averageMetrics, true)
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeFinalResultsSummaryTableRow(builder *strings.Builder, instance string, metrics map[string]finalResultsSummaryMetric, boldInstance bool) {
	deviationHighlights, successHighlights := finalResultsSummaryHighlights(metrics)
	instanceCell := html.EscapeString(instance)
	if boldInstance {
		instanceCell = "<strong>" + instanceCell + "</strong>"
	}

	fmt.Fprintf(builder, "<tr><td>%s</td>", instanceCell)
	for _, heuristic := range finalResultsSummaryHeuristics {
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

func finalResultsSummaryHighlights(metrics map[string]finalResultsSummaryMetric) (map[string]bool, map[string]bool) {
	deviationHighlights := make(map[string]bool, len(finalResultsSummaryHeuristics))
	successHighlights := make(map[string]bool, len(finalResultsSummaryHeuristics))
	minDeviation := math.Inf(1)
	maxSuccess := math.Inf(-1)

	for _, heuristic := range finalResultsSummaryHeuristics {
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

	for _, heuristic := range finalResultsSummaryHeuristics {
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

func highlightedHeuristicDisplayNames(highlights map[string]bool) []heuristicDisplay {
	names := make([]heuristicDisplay, 0, len(highlights))
	for _, heuristic := range finalResultsSummaryHeuristics {
		if highlights[heuristic] {
			names = append(names, heuristicDisplay{
				heuristic: heuristic,
				display:   heuristicDisplayName(heuristic),
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

func heuristicCountList(counts map[string]int, total int) string {
	parts := make([]string, 0, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		parts = append(parts, fmt.Sprintf("%s %d/%d", heuristicDisplayName(heuristic), counts[heuristic], total))
	}
	return strings.Join(parts, ", ")
}

func heuristicFloatList(values map[string]float64, format string) string {
	parts := make([]string, 0, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		value, ok := values[heuristic]
		if !ok {
			continue
		}
		parts = append(parts, fmt.Sprintf("%s "+format, heuristicDisplayName(heuristic), value))
	}
	return strings.Join(parts, ", ")
}

func comparisonHeuristics() []string {
	return finalResultsSummaryHeuristics[1:]
}

func finalResultsSummaryMetricCell(value float64, bold bool) string {
	valueText := fmt.Sprintf("%.2f", value)
	if bold {
		return "<strong>" + valueText + "</strong>"
	}

	return valueText
}

func saveFinalPairwisePerformanceReport(path string, rows []finalResultsSummaryRow) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	comparisons := finalPairwisePerformanceComparisons(rows)

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
			html.EscapeString(heuristicDisplayName(comparison.left)),
			html.EscapeString(heuristicDisplayName(comparison.right)),
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

func finalPairwisePerformanceComparisons(rows []finalResultsSummaryRow) []finalPairwisePerformanceComparison {
	comparisons := make([]finalPairwisePerformanceComparison, 0)
	for _, heuristic := range comparisonHeuristics() {
		comparison := calculateFinalPairwisePerformanceComparison(rows, heuristic, heuristicBaseline)
		if comparison.count > 0 {
			comparisons = append(comparisons, comparison)
		}
	}
	heuristics := comparisonHeuristics()
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

func calculateFinalPairwisePerformanceComparison(rows []finalResultsSummaryRow, left, right string) finalPairwisePerformanceComparison {
	const epsilon = 1e-9
	comparison := finalPairwisePerformanceComparison{left: left, right: right}

	for _, row := range rows {
		leftMetric, leftOk := row.metrics[left]
		rightMetric, rightOk := row.metrics[right]
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

func saveFinalConvergenceSummaryReport(path string, rows []finalResultsSummaryRow) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	averages := averageConvergenceByHeuristic(rows, nil)

	var builder strings.Builder
	builder.WriteString("# Convergence Summary\n\n")
	builder.WriteString("Each value is the average iteration where the best run solution was found, expressed as a percentage of the configured iteration budget. Lower values mean earlier convergence.\n\n")
	writeFinalConvergenceFindings(&builder, rows, averages)
	builder.WriteString("\n")
	writeFinalConvergenceTable(&builder, rows, averages)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeFinalConvergenceFindings(builder *strings.Builder, rows []finalResultsSummaryRow, averages map[string]float64) {
	winCounts := make(map[string]int, len(finalResultsSummaryHeuristics))
	for _, row := range rows {
		highlights := lowestConvergenceHighlights(row.metrics)
		for _, heuristic := range finalResultsSummaryHeuristics {
			if highlights[heuristic] {
				winCounts[heuristic]++
			}
		}
	}

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder,
		"- **Average best-iteration position: %s.**\n",
		heuristicFloatList(averages, "%.2f%%"))
	fmt.Fprintf(builder,
		"- **Earliest-or-tied convergence counts: %s.**\n",
		heuristicCountList(winCounts, len(rows)))
}

func writeFinalConvergenceTable(builder *strings.Builder, rows []finalResultsSummaryRow, averages map[string]float64) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th>")
	for _, heuristic := range finalResultsSummaryHeuristics {
		fmt.Fprintf(builder, "<th>%s best iter [%%]</th>", html.EscapeString(heuristicDisplayName(heuristic)))
	}
	builder.WriteString("</tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		highlights := lowestConvergenceHighlights(row.metrics)
		fmt.Fprintf(builder, "<tr><td>%s</td>", html.EscapeString(row.instance))
		for _, heuristic := range finalResultsSummaryHeuristics {
			metric, ok := row.metrics[heuristic]
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
	for _, heuristic := range finalResultsSummaryHeuristics {
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

func averageConvergenceByHeuristic(rows []finalResultsSummaryRow, allowedInstances map[string]struct{}) map[string]float64 {
	sums := make(map[string]float64, len(finalResultsSummaryHeuristics))
	counts := make(map[string]int, len(finalResultsSummaryHeuristics))
	for _, row := range rows {
		if allowedInstances != nil {
			if _, ok := allowedInstances[row.instance]; !ok {
				continue
			}
		}

		for _, heuristic := range finalResultsSummaryHeuristics {
			metric, ok := row.metrics[heuristic]
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

	averages := make(map[string]float64, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		if counts[heuristic] == 0 {
			continue
		}
		averages[heuristic] = sums[heuristic] / float64(counts[heuristic])
	}

	return averages
}

func lowestConvergenceHighlights(metrics map[string]finalResultsSummaryMetric) map[string]bool {
	values := make(map[string]float64, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		metric, ok := metrics[heuristic]
		if !ok || metric.Iterations <= 0 {
			continue
		}
		values[heuristic] = convergencePercent(metric)
	}

	return lowestFloatHighlights(values, false)
}

func convergencePercent(metric finalResultsSummaryMetric) float64 {
	if metric.Iterations <= 0 {
		return 0
	}

	return 100.0 * metric.AverageBestIteration / float64(metric.Iterations)
}

func saveFinalThreeOptComparisonReport(path string, finalRows, finalThreeOptRows []finalResultsSummaryRow) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	finalAverages := averageFinalResultsSummaryMetrics(finalRows)
	finalThreeOptAverages := averageFinalResultsSummaryMetrics(finalThreeOptRows)

	var builder strings.Builder
	builder.WriteString("# Reduced 3-Opt Impact\n\n")
	builder.WriteString("This report compares the final MMAS experiments without local search against the final MMAS experiments with reduced 3-opt enabled.\n\n")
	writeFinalThreeOptComparisonFindings(&builder, finalAverages, finalThreeOptAverages)
	builder.WriteString("\n")
	writeFinalThreeOptImpactTable(&builder, finalAverages, finalThreeOptAverages)
	builder.WriteString("\n")
	writeFinalThreeOptSignalTable(&builder, finalAverages, finalThreeOptAverages)
	builder.WriteString("\n")
	builder.WriteString("Deviation gain vs baseline is `baseline average best deviation - heuristic average best deviation`, so positive values mean that the heuristic improved over the baseline. Signal remaining is the share of this deviation gain still visible after enabling reduced 3-opt.\n")

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeFinalThreeOptComparisonFindings(builder *strings.Builder, finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric) {
	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder,
		"- **Reduced 3-opt average-best-deviation deltas: %s.**\n",
		finalThreeOptDeviationDeltaList(finalAverages, finalThreeOptAverages))
	fmt.Fprintf(builder,
		"- **Reduced 3-opt success-rate deltas: %s.**\n",
		finalThreeOptSuccessDeltaList(finalAverages, finalThreeOptAverages))
	fmt.Fprintf(builder,
		"- **Deviation gain over baseline with and without 3-opt: %s.**\n",
		finalThreeOptGainList(finalAverages, finalThreeOptAverages))
	fmt.Fprintf(builder, "- **Signal remaining after enabling 3-opt: %s.**\n", finalThreeOptSignalRemainingList(finalAverages, finalThreeOptAverages))
	builder.WriteString("- **This supports treating reduced 3-opt as a strong local-search layer that partially hides the construction heuristic effect.**\n")
}

func writeFinalThreeOptImpactTable(builder *strings.Builder, finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric) {
	builder.WriteString("## Overall Effect\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, heuristic := range finalResultsSummaryHeuristics {
		without, withoutOk := finalAverages[heuristic]
		with, withOk := finalThreeOptAverages[heuristic]
		if !withoutOk || !withOk {
			continue
		}
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(heuristicDisplayName(heuristic)),
			without.AverageMinDeviation,
			with.AverageMinDeviation,
			without.SuccessRate,
			with.SuccessRate)
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeFinalThreeOptSignalTable(builder *strings.Builder, finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric) {
	builder.WriteString("## Heuristic Signal\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Heuristic</th><th>Dev. gain vs baseline without 3-opt [pp]</th><th>Dev. gain vs baseline with 3-opt [pp]</th><th>Signal remaining [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, heuristic := range comparisonHeuristics() {
		if !hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages, heuristic) {
			continue
		}

		gainWithout := deviationGainVsBaseline(finalAverages, heuristic)
		gainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristic)

		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(heuristicDisplayName(heuristic)),
			formatSignedFloat(gainWithout),
			formatSignedFloat(gainWith),
			signalRemainingPercent(gainWithout, gainWith))
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func finalThreeOptDeviationDeltaList(finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric) string {
	parts := make([]string, 0, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		if !hasPairedMetrics(finalAverages, finalThreeOptAverages, heuristic) {
			continue
		}
		delta := finalAverages[heuristic].AverageMinDeviation - finalThreeOptAverages[heuristic].AverageMinDeviation
		parts = append(parts, fmt.Sprintf("%s %s pp", heuristicDisplayName(heuristic), formatSignedFloat(delta)))
	}
	return strings.Join(parts, ", ")
}

func finalThreeOptSuccessDeltaList(finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric) string {
	parts := make([]string, 0, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		if !hasPairedMetrics(finalAverages, finalThreeOptAverages, heuristic) {
			continue
		}
		delta := finalThreeOptAverages[heuristic].SuccessRate - finalAverages[heuristic].SuccessRate
		parts = append(parts, fmt.Sprintf("%s %s pp", heuristicDisplayName(heuristic), formatSignedFloat(delta)))
	}
	return strings.Join(parts, ", ")
}

func finalThreeOptGainList(finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric) string {
	parts := make([]string, 0, len(comparisonHeuristics()))
	for _, heuristic := range comparisonHeuristics() {
		if !hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages, heuristic) {
			continue
		}
		gainWithout := deviationGainVsBaseline(finalAverages, heuristic)
		gainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristic)
		parts = append(parts, fmt.Sprintf("%s %s -> %s pp", heuristicDisplayName(heuristic), formatSignedFloat(gainWithout), formatSignedFloat(gainWith)))
	}
	return strings.Join(parts, ", ")
}

func finalThreeOptSignalRemainingList(finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric) string {
	parts := make([]string, 0, len(comparisonHeuristics()))
	for _, heuristic := range comparisonHeuristics() {
		if !hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages, heuristic) {
			continue
		}
		gainWithout := deviationGainVsBaseline(finalAverages, heuristic)
		gainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristic)
		parts = append(parts, fmt.Sprintf("%s %.2f%%", heuristicDisplayName(heuristic), signalRemainingPercent(gainWithout, gainWith)))
	}
	return strings.Join(parts, ", ")
}

func hasThreeOptComparisonMetrics(finalAverages, finalThreeOptAverages map[string]finalResultsSummaryMetric, heuristic string) bool {
	return hasPairedMetrics(finalAverages, finalThreeOptAverages, heuristic) &&
		hasPairedMetrics(finalAverages, finalThreeOptAverages, heuristicBaseline)
}

func hasPairedMetrics(left, right map[string]finalResultsSummaryMetric, heuristic string) bool {
	_, leftOk := left[heuristic]
	_, rightOk := right[heuristic]
	return leftOk && rightOk
}

func deviationGainVsBaseline(averages map[string]finalResultsSummaryMetric, heuristic string) float64 {
	return averages[heuristicBaseline].AverageMinDeviation - averages[heuristic].AverageMinDeviation
}

func signalRemainingPercent(without, with float64) float64 {
	if math.Abs(without) < 1e-9 {
		return 0
	}

	return 100.0 * math.Abs(with) / math.Abs(without)
}

func saveStructuralPerformanceLinkReport(path string, rows []finalResultsSummaryRow, analyses []structure.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	structuralRows := filterAnalysesWithFoundOptimalEdges(sortedStructuralAnalyses(analyses))
	structuralTotals := structuralSimilarityTotals(structuralRows)
	allowedInstances := make(map[string]struct{}, len(structuralRows))
	for _, row := range structuralRows {
		allowedInstances[row.Instance] = struct{}{}
	}

	performance := averagePerformanceByHeuristic(rows, allowedInstances)
	msaPrecision := ratio(structuralTotals.msaOptimalEdges, structuralTotals.msaEdges)
	cycleCoverPrecision := ratio(structuralTotals.cycleCoverOptimalEdges, structuralTotals.cycleCoverEdges)
	patchingPrecision := ratio(structuralTotals.patchingOptimalEdges, structuralTotals.patchingEdges)
	msaRecall := ratio(structuralTotals.msaOptimalEdges, structuralTotals.foundOptimalEdges)
	cycleCoverRecall := ratio(structuralTotals.cycleCoverOptimalEdges, structuralTotals.foundOptimalEdges)
	patchingRecall := ratio(structuralTotals.patchingOptimalEdges, structuralTotals.foundOptimalEdges)

	structuralPrecision := map[string]float64{
		heuristicStrictMsa:             msaPrecision,
		heuristicCycleCover:            cycleCoverPrecision,
		heuristicCycleCoverMsaPatching: patchingPrecision,
	}
	structuralRecall := map[string]float64{
		heuristicStrictMsa:             msaRecall,
		heuristicCycleCover:            cycleCoverRecall,
		heuristicCycleCoverMsaPatching: patchingRecall,
	}

	includedHeuristics := make([]string, 0, len(comparisonHeuristics()))
	precisionValues := make(map[string]float64)
	recallValues := make(map[string]float64)
	deviationValues := make(map[string]float64)
	successValues := make(map[string]float64)
	for _, heuristic := range comparisonHeuristics() {
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
		writeStructuralPerformanceLinkRow(&builder, heuristic, structuralPrecision[heuristic], structuralRecall[heuristic], performance[heuristic], precisionHighlights, recallHighlights, deviationHighlights, successHighlights)
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeStructuralPerformanceLinkRow(builder *strings.Builder, heuristic string, precision, recall float64, performance finalResultsSummaryMetric, precisionHighlights, recallHighlights, deviationHighlights, successHighlights map[string]bool) {
	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		html.EscapeString(heuristicDisplayName(heuristic)),
		finalResultsSummaryMetricCell(100*precision, precisionHighlights[heuristic]),
		finalResultsSummaryMetricCell(100*recall, recallHighlights[heuristic]),
		finalResultsSummaryMetricCell(performance.AverageMinDeviation, deviationHighlights[heuristic]),
		finalResultsSummaryMetricCell(performance.SuccessRate, successHighlights[heuristic]))
}

func averagePerformanceByHeuristic(rows []finalResultsSummaryRow, allowedInstances map[string]struct{}) map[string]finalResultsSummaryMetric {
	totals := make(map[string]finalResultsSummaryMetric, len(finalResultsSummaryHeuristics))
	counts := make(map[string]int, len(finalResultsSummaryHeuristics))
	for _, row := range rows {
		if allowedInstances != nil {
			if _, ok := allowedInstances[row.instance]; !ok {
				continue
			}
		}

		for _, heuristic := range finalResultsSummaryHeuristics {
			metric, ok := row.metrics[heuristic]
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

	averages := make(map[string]finalResultsSummaryMetric, len(finalResultsSummaryHeuristics))
	for _, heuristic := range finalResultsSummaryHeuristics {
		count := counts[heuristic]
		if count == 0 {
			continue
		}

		total := totals[heuristic]
		averages[heuristic] = finalResultsSummaryMetric{
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
	return fmt.Sprintf("%+.2f", value)
}
