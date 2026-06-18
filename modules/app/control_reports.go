package app

import (
	"atsp_aco_msa/modules/project"
	"errors"
	"fmt"
	"html"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

type randomSparseControlRow struct {
	instance                           string
	msaAverageBestDeviation            float64
	randomAverageBestDeviation         float64
	randomBestAverageBestDeviation     float64
	averageBestDeviationDelta          float64
	msaSuccessRate                     float64
	randomSuccessRate                  float64
	successRateDelta                   float64
	randomSeedCount                    int
	randomBestAverageBestDeviationSeed int64
}

type randomSparseControlMissingData struct {
	instance string
	reason   string
}

func saveRandomSparseControlReport(path string, atspsData []AtspData, finalResultsRootPath, controlResultsRootPath string) (bool, error) {
	rows, missingData, err := buildRandomSparseControlRows(atspsData, finalResultsRootPath, controlResultsRootPath)
	if err != nil {
		return false, err
	}
	if len(rows) == 0 {
		return false, nil
	}
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return false, err
	}

	var builder strings.Builder
	builder.WriteString("# Random Sparse Control\n\n")
	builder.WriteString("This sanity check compares Strict MSA against deterministic random sparse masks. Each random mask boosts the same number of directed edges as Strict MSA and uses the same heuristic weight. The comparison reads Strict MSA from the final results and averages the available final-control random seeds for each instance.\n\n")
	writeRandomSparseControlFindings(&builder, rows)
	builder.WriteString("\n")
	writeRandomSparseControlTable(&builder, rows)
	if len(missingData) != 0 {
		writeRandomSparseControlMissingData(&builder, missingData)
	}

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildRandomSparseControlRows(atspsData []AtspData, finalResultsRootPath, controlResultsRootPath string) ([]randomSparseControlRow, []randomSparseControlMissingData, error) {
	rows := make([]randomSparseControlRow, 0, len(atspsData))
	missingData := make([]randomSparseControlMissingData, 0)

	for _, atspData := range atspsData {
		msaMetric, ok, err := readFinalMsaHeuristicControlMetric(atspData, finalResultsRootPath)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			missingData = append(missingData, randomSparseControlMissingData{instance: atspData.Name, reason: "missing final MSA result"})
			continue
		}

		controlAtspData := project.WithExperimentOutputRoot(atspData, controlResultsRootPath)
		randomStatistics, err := readStatistics(resultFilePathForHeuristic(controlAtspData, heuristicRandomSparse))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, randomSparseControlMissingData{instance: atspData.Name, reason: "missing final-control random-sparse result CSV"})
				continue
			}
			return nil, nil, err
		}
		randomStatisticsForWeight := statisticsForHeuristicWeightAll(randomStatistics, msaMetric.heuristicWeight)
		if len(randomStatisticsForWeight) == 0 {
			missingData = append(missingData, randomSparseControlMissingData{instance: atspData.Name, reason: fmt.Sprintf("missing random-sparse rows for heuristic weight %.2f", msaMetric.heuristicWeight)})
			continue
		}

		randomAverageBestDeviation := 0.0
		randomSuccessRate := 0.0
		randomBestAverageBestDeviation := math.Inf(1)
		randomBestAverageBestDeviationSeed := int64(0)
		for _, statistic := range randomStatisticsForWeight {
			randomAverageBestDeviation += statistic.averageBestDeviation
			randomSuccessRate += statistic.successRate
			if statistic.averageBestDeviation < randomBestAverageBestDeviation {
				randomBestAverageBestDeviation = statistic.averageBestDeviation
				randomBestAverageBestDeviationSeed = statistic.randomSeed
			}
		}
		randomAverageBestDeviation /= float64(len(randomStatisticsForWeight))
		randomSuccessRate /= float64(len(randomStatisticsForWeight))

		rows = append(rows, randomSparseControlRow{
			instance:                           atspData.Name,
			msaAverageBestDeviation:            msaMetric.averageMinDeviation,
			randomAverageBestDeviation:         randomAverageBestDeviation,
			randomBestAverageBestDeviation:     randomBestAverageBestDeviation,
			averageBestDeviationDelta:          msaMetric.averageMinDeviation - randomAverageBestDeviation,
			msaSuccessRate:                     msaMetric.successRate,
			randomSuccessRate:                  randomSuccessRate,
			successRateDelta:                   msaMetric.successRate - randomSuccessRate,
			randomSeedCount:                    len(randomStatisticsForWeight),
			randomBestAverageBestDeviationSeed: randomBestAverageBestDeviationSeed,
		})
	}

	sort.SliceStable(rows, func(i, j int) bool {
		return rows[i].instance < rows[j].instance
	})
	sort.SliceStable(missingData, func(i, j int) bool {
		return missingData[i].instance < missingData[j].instance
	})

	return rows, missingData, nil
}

func statisticsForHeuristicWeight(statistics []ExperimentsDataStatistics, heuristicWeight float64) (ExperimentsDataStatistics, bool) {
	for _, statistic := range statistics {
		if math.Abs(statistic.heuristicWeight-heuristicWeight) < 1e-9 {
			return statistic, true
		}
	}

	return ExperimentsDataStatistics{}, false
}

func statisticsForHeuristicWeightAll(statistics []ExperimentsDataStatistics, heuristicWeight float64) []ExperimentsDataStatistics {
	rows := make([]ExperimentsDataStatistics, 0)
	for _, statistic := range statistics {
		if math.Abs(statistic.heuristicWeight-heuristicWeight) < 1e-9 {
			rows = append(rows, statistic)
		}
	}

	return rows
}

func readFinalMsaHeuristicControlMetric(atspData AtspData, finalResultsRootPath string) (finalResultsSummaryMetric, bool, error) {
	return readMsaControlMetric(atspData, finalResultsRootPath, heuristicStrictMsa)
}

func readMsaControlMetric(atspData AtspData, finalResultsRootPath, referenceHeuristic string) (finalResultsSummaryMetric, bool, error) {
	finalAtspData := project.WithExperimentOutputRoot(atspData, finalResultsRootPath)
	statistics, err := readHeuristicStatistics(finalAtspData.ResultFilePath)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return finalResultsSummaryMetric{}, false, nil
		}
		return finalResultsSummaryMetric{}, false, err
	}

	return msaHeuristicControlMetricForWeight(statistics, referenceHeuristic, finalStrictMsaHeuristicWeight)
}

func msaHeuristicControlMetricForWeight(statistics []HeuristicExperimentStatistics, referenceHeuristic string, heuristicWeight float64) (finalResultsSummaryMetric, bool, error) {
	for _, statistic := range statistics {
		if statistic.heuristic != referenceHeuristic || math.Abs(statistic.statistics.heuristicWeight-heuristicWeight) >= 1e-9 {
			continue
		}

		return finalResultsSummaryMetric{
			averageMinDeviation:  statistic.statistics.averageBestDeviation,
			successRate:          statistic.statistics.successRate,
			averageBestIteration: statistic.statistics.averageBestAtIteration,
			heuristicWeight:      statistic.statistics.heuristicWeight,
			iterations:           statistic.statistics.iterations,
		}, true, nil
	}

	return finalResultsSummaryMetric{}, false, nil
}

func writeRandomSparseControlFindings(builder *strings.Builder, rows []randomSparseControlRow) {
	summary := randomSparseControlSummary(rows)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **Strict MSA had lower average best deviation than the random-sparse mean in %d/%d instances.**\n", summary.msaWins, len(rows))
	fmt.Fprintf(builder, "- **Mean average best deviation: Strict MSA %.2f%%, random sparse %.2f%%, delta %s pp.**\n",
		summary.meanMsaAverageBestDeviation,
		summary.meanRandomAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta))
	fmt.Fprintf(builder, "- **Mean success rate: Strict MSA %.2f%%, random sparse %.2f%%, delta %s pp.**\n",
		summary.meanMsaSuccessRate,
		summary.meanRandomSuccessRate,
		formatSignedFloat(summary.meanSuccessRateDelta))
	fmt.Fprintf(builder, "- Strict MSA also beat the best random seed in %d/%d instances.\n", summary.msaWinsAgainstBestRandomSeed, len(rows))
	fmt.Fprintf(builder, "- Two-sided sign-test p-value for average-best-deviation wins/losses: %.6f.\n", summary.signTestPValue)
}

type randomSparseControlSummaryData struct {
	msaWins                        int
	randomWins                     int
	ties                           int
	msaWinsAgainstBestRandomSeed   int
	meanMsaAverageBestDeviation    float64
	meanRandomAverageBestDeviation float64
	meanAverageBestDeviationDelta  float64
	meanMsaSuccessRate             float64
	meanRandomSuccessRate          float64
	meanSuccessRateDelta           float64
	signTestPValue                 float64
}

func randomSparseControlSummary(rows []randomSparseControlRow) randomSparseControlSummaryData {
	var summary randomSparseControlSummaryData
	for _, row := range rows {
		if row.averageBestDeviationDelta < -1e-9 {
			summary.msaWins++
		} else if row.averageBestDeviationDelta > 1e-9 {
			summary.randomWins++
		} else {
			summary.ties++
		}
		if row.msaAverageBestDeviation < row.randomBestAverageBestDeviation {
			summary.msaWinsAgainstBestRandomSeed++
		}

		summary.meanMsaAverageBestDeviation += row.msaAverageBestDeviation
		summary.meanRandomAverageBestDeviation += row.randomAverageBestDeviation
		summary.meanAverageBestDeviationDelta += row.averageBestDeviationDelta
		summary.meanMsaSuccessRate += row.msaSuccessRate
		summary.meanRandomSuccessRate += row.randomSuccessRate
		summary.meanSuccessRateDelta += row.successRateDelta
	}

	count := float64(len(rows))
	summary.meanMsaAverageBestDeviation /= count
	summary.meanRandomAverageBestDeviation /= count
	summary.meanAverageBestDeviationDelta /= count
	summary.meanMsaSuccessRate /= count
	summary.meanRandomSuccessRate /= count
	summary.meanSuccessRateDelta /= count
	summary.signTestPValue = twoSidedSignTestPValue(summary.msaWins, summary.randomWins)

	return summary
}

func twoSidedSignTestPValue(wins, losses int) float64 {
	n := wins + losses
	if n == 0 {
		return 1.0
	}

	observed := maxIntValue(wins, losses)
	probability := 0.0
	for k := observed; k <= n; k++ {
		probability += binomialCoefficient(n, k) / math.Pow(2, float64(n))
	}

	probability *= 2
	if probability > 1 {
		return 1
	}

	return probability
}

func binomialCoefficient(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}

	coefficient := 1.0
	for i := 1; i <= k; i++ {
		coefficient *= float64(n-k+i) / float64(i)
	}

	return coefficient
}

func writeRandomSparseControlTable(builder *strings.Builder, rows []randomSparseControlRow) {
	summary := randomSparseControlSummary(rows)

	builder.WriteString("## Per-instance comparison\n\n")
	builder.WriteString("Negative delta means Strict MSA had lower average best deviation than the random-sparse mean.\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>Strict MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>Strict MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		msaWins := row.averageBestDeviationDelta < 0
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%.2f (seed %d)</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td></tr>\n",
			html.EscapeString(row.instance),
			finalResultsSummaryMetricCell(row.msaAverageBestDeviation, msaWins),
			finalResultsSummaryMetricCell(row.randomAverageBestDeviation, !msaWins),
			row.randomBestAverageBestDeviation,
			row.randomBestAverageBestDeviationSeed,
			formatSignedFloat(row.averageBestDeviationDelta),
			row.msaSuccessRate,
			row.randomSuccessRate,
			row.randomSeedCount)
	}
	fmt.Fprintf(builder,
		"<tr><td><strong>Average</strong></td><td align=\"right\"><strong>%.2f</strong></td><td align=\"right\"><strong>%.2f</strong></td><td></td><td align=\"right\"><strong>%s</strong></td><td align=\"right\"><strong>%.2f</strong></td><td align=\"right\"><strong>%.2f</strong></td><td></td></tr>\n",
		summary.meanMsaAverageBestDeviation,
		summary.meanRandomAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta),
		summary.meanMsaSuccessRate,
		summary.meanRandomSuccessRate)
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeRandomSparseControlMissingData(builder *strings.Builder, missingData []randomSparseControlMissingData) {
	builder.WriteString("\n## Missing data\n\n")
	builder.WriteString("| Instance | Reason |\n")
	builder.WriteString("|---|---|\n")
	for _, missing := range missingData {
		fmt.Fprintf(builder, "| %s | %s |\n", missing.instance, missing.reason)
	}
}

type distanceRankedSparseControlRow struct {
	instance                    string
	msaAverageBestDeviation     float64
	controlAverageBestDeviation float64
	averageBestDeviationDelta   float64
	msaSuccessRate              float64
	controlSuccessRate          float64
	successRateDelta            float64
}

type distanceRankedSparseControlMissingData struct {
	instance string
	reason   string
}

func saveDistanceRankedSparseControlReport(path string, atspsData []AtspData, finalResultsRootPath, controlResultsRootPath string) (bool, error) {
	return saveDistanceRankedSparseControlReportForReference(path, atspsData, finalResultsRootPath, controlResultsRootPath, heuristicStrictMsa, "Strict MSA", heuristicDistanceRankedSparse, "distance-ranked sparse")
}

func saveDistanceRankedSparseControlReportForReference(path string, atspsData []AtspData, finalResultsRootPath, controlResultsRootPath, referenceHeuristic, referenceName, controlHeuristic, controlName string) (bool, error) {
	rows, missingData, err := buildDistanceRankedSparseControlRows(atspsData, finalResultsRootPath, controlResultsRootPath, referenceHeuristic, controlHeuristic)
	if err != nil {
		return false, err
	}
	if len(rows) == 0 {
		return false, nil
	}
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return false, err
	}

	var builder strings.Builder
	fmt.Fprintf(&builder, "# %s Distance-ranked Sparse Control\n\n", referenceName)
	fmt.Fprintf(&builder, "This sanity check compares %s against a deterministic sparse mask built from the cheapest directed edges. The control preserves the matching edge-count structure and %s.\n\n", referenceName, controlWeightDescription())
	writeDistanceRankedSparseControlFindings(&builder, rows, referenceName, controlName)
	builder.WriteString("\n")
	writeDistanceRankedSparseControlTable(&builder, rows, referenceName, controlName)
	if len(missingData) != 0 {
		writeDistanceRankedSparseControlMissingData(&builder, missingData)
	}

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildDistanceRankedSparseControlRows(atspsData []AtspData, finalResultsRootPath, controlResultsRootPath, referenceHeuristic, controlHeuristic string) ([]distanceRankedSparseControlRow, []distanceRankedSparseControlMissingData, error) {
	rows := make([]distanceRankedSparseControlRow, 0, len(atspsData))
	missingData := make([]distanceRankedSparseControlMissingData, 0)

	for _, atspData := range atspsData {
		msaMetric, ok, err := readMsaControlMetric(atspData, finalResultsRootPath, referenceHeuristic)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.Name, reason: "missing final " + referenceHeuristic + " result"})
			continue
		}

		controlAtspData := project.WithExperimentOutputRoot(atspData, controlResultsRootPath)
		controlStatistics, err := readStatistics(resultFilePathForHeuristic(controlAtspData, controlHeuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.Name, reason: "missing final-control " + controlHeuristic + " result CSV"})
				continue
			}
			return nil, nil, err
		}

		controlStatistic, ok := statisticsForHeuristicWeight(controlStatistics, msaMetric.heuristicWeight)
		if !ok {
			missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.Name, reason: fmt.Sprintf("missing %s row for heuristic weight %.2f", controlHeuristic, msaMetric.heuristicWeight)})
			continue
		}

		rows = append(rows, distanceRankedSparseControlRow{
			instance:                    atspData.Name,
			msaAverageBestDeviation:     msaMetric.averageMinDeviation,
			controlAverageBestDeviation: controlStatistic.averageBestDeviation,
			averageBestDeviationDelta:   msaMetric.averageMinDeviation - controlStatistic.averageBestDeviation,
			msaSuccessRate:              msaMetric.successRate,
			controlSuccessRate:          controlStatistic.successRate,
			successRateDelta:            msaMetric.successRate - controlStatistic.successRate,
		})
	}

	sort.SliceStable(rows, func(i, j int) bool {
		return rows[i].instance < rows[j].instance
	})
	sort.SliceStable(missingData, func(i, j int) bool {
		return missingData[i].instance < missingData[j].instance
	})

	return rows, missingData, nil
}

func writeDistanceRankedSparseControlFindings(builder *strings.Builder, rows []distanceRankedSparseControlRow, referenceName, controlName string) {
	summary := distanceRankedSparseControlSummary(rows)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **%s had lower average best deviation than the %s control in %d/%d instances.**\n", referenceName, controlName, summary.msaWins, len(rows))
	fmt.Fprintf(builder, "- **Mean average best deviation: %s %.2f%%, %s %.2f%%, delta %s pp.**\n",
		referenceName,
		summary.meanMsaAverageBestDeviation,
		controlName,
		summary.meanControlAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta))
	fmt.Fprintf(builder, "- **Mean success rate: %s %.2f%%, %s %.2f%%, delta %s pp.**\n",
		referenceName,
		summary.meanMsaSuccessRate,
		controlName,
		summary.meanControlSuccessRate,
		formatSignedFloat(summary.meanSuccessRateDelta))
	fmt.Fprintf(builder, "- Two-sided sign-test p-value for average-best-deviation wins/losses: %.6f.\n", summary.signTestPValue)
}

type distanceRankedSparseControlSummaryData struct {
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

func distanceRankedSparseControlSummary(rows []distanceRankedSparseControlRow) distanceRankedSparseControlSummaryData {
	var summary distanceRankedSparseControlSummaryData
	for _, row := range rows {
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

	count := float64(len(rows))
	summary.meanMsaAverageBestDeviation /= count
	summary.meanControlAverageBestDeviation /= count
	summary.meanAverageBestDeviationDelta /= count
	summary.meanMsaSuccessRate /= count
	summary.meanControlSuccessRate /= count
	summary.meanSuccessRateDelta /= count
	summary.signTestPValue = twoSidedSignTestPValue(summary.msaWins, summary.controlWins)

	return summary
}

func writeDistanceRankedSparseControlTable(builder *strings.Builder, rows []distanceRankedSparseControlRow, referenceName, controlName string) {
	summary := distanceRankedSparseControlSummary(rows)

	builder.WriteString("## Per-instance comparison\n\n")
	fmt.Fprintf(builder, "Negative delta means %s had lower average best deviation than the %s control.\n\n", referenceName, controlName)
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	fmt.Fprintf(builder, "<tr><th>Instance</th><th>%s avg best dev. [%%]</th><th>%s avg best dev. [%%]</th><th>Delta [pp]</th><th>%s success [%%]</th><th>%s success [%%]</th></tr>\n",
		html.EscapeString(referenceName),
		html.EscapeString(controlName),
		html.EscapeString(referenceName),
		html.EscapeString(controlName))
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		msaWins := row.averageBestDeviationDelta < 0
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(row.instance),
			finalResultsSummaryMetricCell(row.msaAverageBestDeviation, msaWins),
			finalResultsSummaryMetricCell(row.controlAverageBestDeviation, !msaWins),
			formatSignedFloat(row.averageBestDeviationDelta),
			row.msaSuccessRate,
			row.controlSuccessRate)
	}
	fmt.Fprintf(builder,
		"<tr><td><strong>Average</strong></td><td align=\"right\"><strong>%.2f</strong></td><td align=\"right\"><strong>%.2f</strong></td><td align=\"right\"><strong>%s</strong></td><td align=\"right\"><strong>%.2f</strong></td><td align=\"right\"><strong>%.2f</strong></td></tr>\n",
		summary.meanMsaAverageBestDeviation,
		summary.meanControlAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta),
		summary.meanMsaSuccessRate,
		summary.meanControlSuccessRate)
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeDistanceRankedSparseControlMissingData(builder *strings.Builder, missingData []distanceRankedSparseControlMissingData) {
	builder.WriteString("\n## Missing data\n\n")
	builder.WriteString("| Instance | Reason |\n")
	builder.WriteString("|---|---|\n")
	for _, missing := range missingData {
		fmt.Fprintf(builder, "| %s | %s |\n", missing.instance, missing.reason)
	}
}

type seededSparseControlReportConfig struct {
	title              string
	description        string
	referenceHeuristic string
	referenceName      string
	controlName        string
	controlHeuristic   string
	missingResultLabel string
}

type seededSparseControlRow struct {
	instance                            string
	msaAverageBestDeviation             float64
	controlAverageBestDeviation         float64
	controlBestAverageBestDeviation     float64
	averageBestDeviationDelta           float64
	msaSuccessRate                      float64
	controlSuccessRate                  float64
	successRateDelta                    float64
	seedCount                           int
	controlBestAverageBestDeviationSeed int64
}

type seededSparseControlMissingData struct {
	instance string
	reason   string
}

func saveSeededSparseControlReport(path string, atspsData []AtspData, finalResultsRootPath, controlResultsRootPath string, config seededSparseControlReportConfig) (bool, error) {
	rows, missingData, err := buildSeededSparseControlRows(atspsData, finalResultsRootPath, controlResultsRootPath, config)
	if err != nil {
		return false, err
	}
	if len(rows) == 0 {
		return false, nil
	}
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return false, err
	}

	var builder strings.Builder
	fmt.Fprintf(&builder, "# %s\n\n", config.title)
	builder.WriteString(config.description)
	builder.WriteString("\n\n")
	writeSeededSparseControlFindings(&builder, rows, config)
	builder.WriteString("\n")
	writeSeededSparseControlTable(&builder, rows, config)
	if len(missingData) != 0 {
		writeSeededSparseControlMissingData(&builder, missingData)
	}

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildSeededSparseControlRows(atspsData []AtspData, finalResultsRootPath, controlResultsRootPath string, config seededSparseControlReportConfig) ([]seededSparseControlRow, []seededSparseControlMissingData, error) {
	rows := make([]seededSparseControlRow, 0, len(atspsData))
	missingData := make([]seededSparseControlMissingData, 0)

	for _, atspData := range atspsData {
		msaMetric, ok, err := readMsaControlMetric(atspData, finalResultsRootPath, config.referenceHeuristic)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			missingData = append(missingData, seededSparseControlMissingData{instance: atspData.Name, reason: "missing final " + config.referenceHeuristic + " result"})
			continue
		}

		controlAtspData := project.WithExperimentOutputRoot(atspData, controlResultsRootPath)
		controlStatistics, err := readStatistics(resultFilePathForHeuristic(controlAtspData, config.controlHeuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, seededSparseControlMissingData{instance: atspData.Name, reason: "missing final-control " + config.missingResultLabel})
				continue
			}
			return nil, nil, err
		}
		controlStatisticsForWeight := statisticsForHeuristicWeightAll(controlStatistics, msaMetric.heuristicWeight)
		if len(controlStatisticsForWeight) == 0 {
			missingData = append(missingData, seededSparseControlMissingData{instance: atspData.Name, reason: fmt.Sprintf("missing %s rows for heuristic weight %.2f", config.missingResultLabel, msaMetric.heuristicWeight)})
			continue
		}

		controlAverageBestDeviation := 0.0
		controlSuccessRate := 0.0
		controlBestAverageBestDeviation := math.Inf(1)
		controlBestAverageBestDeviationSeed := int64(0)
		for _, statistic := range controlStatisticsForWeight {
			controlAverageBestDeviation += statistic.averageBestDeviation
			controlSuccessRate += statistic.successRate
			if statistic.averageBestDeviation < controlBestAverageBestDeviation {
				controlBestAverageBestDeviation = statistic.averageBestDeviation
				controlBestAverageBestDeviationSeed = statistic.randomSeed
			}
		}
		controlAverageBestDeviation /= float64(len(controlStatisticsForWeight))
		controlSuccessRate /= float64(len(controlStatisticsForWeight))

		rows = append(rows, seededSparseControlRow{
			instance:                            atspData.Name,
			msaAverageBestDeviation:             msaMetric.averageMinDeviation,
			controlAverageBestDeviation:         controlAverageBestDeviation,
			controlBestAverageBestDeviation:     controlBestAverageBestDeviation,
			averageBestDeviationDelta:           msaMetric.averageMinDeviation - controlAverageBestDeviation,
			msaSuccessRate:                      msaMetric.successRate,
			controlSuccessRate:                  controlSuccessRate,
			successRateDelta:                    msaMetric.successRate - controlSuccessRate,
			seedCount:                           len(controlStatisticsForWeight),
			controlBestAverageBestDeviationSeed: controlBestAverageBestDeviationSeed,
		})
	}

	sort.SliceStable(rows, func(i, j int) bool {
		return rows[i].instance < rows[j].instance
	})
	sort.SliceStable(missingData, func(i, j int) bool {
		return missingData[i].instance < missingData[j].instance
	})

	return rows, missingData, nil
}

type seededSparseControlSummaryData struct {
	msaWins                         int
	controlWins                     int
	ties                            int
	msaWinsAgainstBestControlSeed   int
	meanMsaAverageBestDeviation     float64
	meanControlAverageBestDeviation float64
	meanAverageBestDeviationDelta   float64
	meanMsaSuccessRate              float64
	meanControlSuccessRate          float64
	meanSuccessRateDelta            float64
	signTestPValue                  float64
}

func seededSparseControlSummary(rows []seededSparseControlRow) seededSparseControlSummaryData {
	var summary seededSparseControlSummaryData
	for _, row := range rows {
		if row.averageBestDeviationDelta < -1e-9 {
			summary.msaWins++
		} else if row.averageBestDeviationDelta > 1e-9 {
			summary.controlWins++
		} else {
			summary.ties++
		}
		if row.msaAverageBestDeviation < row.controlBestAverageBestDeviation {
			summary.msaWinsAgainstBestControlSeed++
		}

		summary.meanMsaAverageBestDeviation += row.msaAverageBestDeviation
		summary.meanControlAverageBestDeviation += row.controlAverageBestDeviation
		summary.meanAverageBestDeviationDelta += row.averageBestDeviationDelta
		summary.meanMsaSuccessRate += row.msaSuccessRate
		summary.meanControlSuccessRate += row.controlSuccessRate
		summary.meanSuccessRateDelta += row.successRateDelta
	}

	count := float64(len(rows))
	summary.meanMsaAverageBestDeviation /= count
	summary.meanControlAverageBestDeviation /= count
	summary.meanAverageBestDeviationDelta /= count
	summary.meanMsaSuccessRate /= count
	summary.meanControlSuccessRate /= count
	summary.meanSuccessRateDelta /= count
	summary.signTestPValue = twoSidedSignTestPValue(summary.msaWins, summary.controlWins)

	return summary
}

func writeSeededSparseControlFindings(builder *strings.Builder, rows []seededSparseControlRow, config seededSparseControlReportConfig) {
	summary := seededSparseControlSummary(rows)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **%s had lower average best deviation than the %s mean in %d/%d instances.**\n", config.referenceName, config.controlName, summary.msaWins, len(rows))
	fmt.Fprintf(builder, "- **Mean average best deviation: %s %.2f%%, %s %.2f%%, delta %s pp.**\n",
		config.referenceName,
		summary.meanMsaAverageBestDeviation,
		config.controlName,
		summary.meanControlAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta))
	fmt.Fprintf(builder, "- **Mean success rate: %s %.2f%%, %s %.2f%%, delta %s pp.**\n",
		config.referenceName,
		summary.meanMsaSuccessRate,
		config.controlName,
		summary.meanControlSuccessRate,
		formatSignedFloat(summary.meanSuccessRateDelta))
	fmt.Fprintf(builder, "- %s also beat the best %s seed in %d/%d instances.\n", config.referenceName, config.controlName, summary.msaWinsAgainstBestControlSeed, len(rows))
	fmt.Fprintf(builder, "- Two-sided sign-test p-value for average-best-deviation wins/losses: %.6f.\n", summary.signTestPValue)
}

func writeSeededSparseControlTable(builder *strings.Builder, rows []seededSparseControlRow, config seededSparseControlReportConfig) {
	summary := seededSparseControlSummary(rows)

	builder.WriteString("## Per-instance comparison\n\n")
	fmt.Fprintf(builder, "Negative delta means %s had lower average best deviation than the %s mean.\n\n", config.referenceName, config.controlName)
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	fmt.Fprintf(builder, "<tr><th>Instance</th><th>%s avg best dev. [%%]</th><th>%s mean avg best dev. [%%]</th><th>Best %s avg best dev. [%%]</th><th>Delta [pp]</th><th>%s success [%%]</th><th>%s success [%%]</th><th>Seeds</th></tr>\n",
		html.EscapeString(config.referenceName),
		html.EscapeString(config.controlName),
		html.EscapeString(config.controlName),
		html.EscapeString(config.referenceName),
		html.EscapeString(config.controlName))
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		msaWins := row.averageBestDeviationDelta < 0
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%.2f (seed %d)</td><td align=\"right\">%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td></tr>\n",
			html.EscapeString(row.instance),
			finalResultsSummaryMetricCell(row.msaAverageBestDeviation, msaWins),
			finalResultsSummaryMetricCell(row.controlAverageBestDeviation, !msaWins),
			row.controlBestAverageBestDeviation,
			row.controlBestAverageBestDeviationSeed,
			formatSignedFloat(row.averageBestDeviationDelta),
			row.msaSuccessRate,
			row.controlSuccessRate,
			row.seedCount)
	}
	fmt.Fprintf(builder,
		"<tr><td><strong>Average</strong></td><td align=\"right\"><strong>%.2f</strong></td><td align=\"right\"><strong>%.2f</strong></td><td></td><td align=\"right\"><strong>%s</strong></td><td align=\"right\"><strong>%.2f</strong></td><td align=\"right\"><strong>%.2f</strong></td><td></td></tr>\n",
		summary.meanMsaAverageBestDeviation,
		summary.meanControlAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta),
		summary.meanMsaSuccessRate,
		summary.meanControlSuccessRate)
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeSeededSparseControlMissingData(builder *strings.Builder, missingData []seededSparseControlMissingData) {
	builder.WriteString("\n## Missing data\n\n")
	builder.WriteString("| Instance | Reason |\n")
	builder.WriteString("|---|---|\n")
	for _, missing := range missingData {
		fmt.Fprintf(builder, "| %s | %s |\n", missing.instance, missing.reason)
	}
}

func saveShuffledMsaControlReport(path string, atspsData []AtspData, finalResultsRootPath, controlResultsRootPath string) (bool, error) {
	return saveSeededSparseControlReport(path, atspsData, finalResultsRootPath, controlResultsRootPath, seededSparseControlReportConfig{
		title:              "Shuffled MSA Control",
		description:        fmt.Sprintf("This sanity check compares Strict MSA against deterministic shuffles of the strict MSA mask. Each shuffle preserves the number and boost values of strict-MSA boosted directed edges, but assigns them to shuffled directed edges. The control %s.", controlWeightDescription()),
		referenceHeuristic: heuristicStrictMsa,
		referenceName:      "Strict MSA",
		controlName:        "shuffled MSA",
		controlHeuristic:   heuristicShuffledMsa,
		missingResultLabel: "shuffled-MSA result CSV",
	})
}

func controlWeightDescription() string {
	return fmt.Sprintf("uses the same `heuristicWeight=%.2f`", finalStrictMsaHeuristicWeight)
}
