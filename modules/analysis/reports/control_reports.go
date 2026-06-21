package reports

import (
	"atsp_aco_msa/modules/experiments"
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

type ControlReportsConfig struct {
	StrictMsaHeuristic                 string
	RandomSparseHeuristic              string
	DistanceRankedSparseHeuristic      string
	ShuffledMsaHeuristic               string
	EvaluationStrictMsaHeuristicWeight float64
	ReadStatistics                     func(string) ([]experiments.ExperimentsDataStatistics, error)
	ReadHeuristicStatistics            func(string) ([]experiments.HeuristicExperimentStatistics, error)
	ResultFilePathForHeuristic         func(project.AtspData, string) string
}

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

func SaveRandomSparseControlReport(path string, atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath string, config ControlReportsConfig) (bool, error) {
	rows, missingData, err := buildRandomSparseControlRows(atspsData, evaluationResultsRootPath, controlResultsRootPath, config)
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
	builder.WriteString("This sanity check compares Strict MSA against deterministic random sparse masks. Each random mask boosts the same number of directed edges as Strict MSA and uses the same heuristic weight. The comparison reads Strict MSA from the evaluation results and averages the available evaluation-control random seeds for each instance.\n\n")
	writeRandomSparseControlFindings(&builder, rows)
	builder.WriteString("\n")
	writeRandomSparseControlTable(&builder, rows)
	if len(missingData) != 0 {
		writeRandomSparseControlMissingData(&builder, missingData)
	}

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildRandomSparseControlRows(atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath string, config ControlReportsConfig) ([]randomSparseControlRow, []randomSparseControlMissingData, error) {
	rows := make([]randomSparseControlRow, 0, len(atspsData))
	missingData := make([]randomSparseControlMissingData, 0)

	for _, atspData := range atspsData {
		msaMetric, ok, err := readEvaluationMsaHeuristicControlMetric(atspData, evaluationResultsRootPath, config)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			missingData = append(missingData, randomSparseControlMissingData{instance: atspData.Name, reason: "missing evaluation MSA result"})
			continue
		}

		controlAtspData := project.WithExperimentOutputRoot(atspData, controlResultsRootPath)
		randomStatistics, err := config.ReadStatistics(config.ResultFilePathForHeuristic(controlAtspData, config.RandomSparseHeuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, randomSparseControlMissingData{instance: atspData.Name, reason: "missing evaluation-control random-sparse result CSV"})
				continue
			}
			return nil, nil, err
		}
		randomStatisticsForWeight := statisticsForHeuristicWeightAll(randomStatistics, msaMetric.HeuristicWeight)
		if len(randomStatisticsForWeight) == 0 {
			missingData = append(missingData, randomSparseControlMissingData{instance: atspData.Name, reason: fmt.Sprintf("missing random-sparse rows for heuristic weight %.2f", msaMetric.HeuristicWeight)})
			continue
		}

		randomAverageBestDeviation := 0.0
		randomSuccessRate := 0.0
		randomBestAverageBestDeviation := math.Inf(1)
		randomBestAverageBestDeviationSeed := int64(0)
		for _, statistic := range randomStatisticsForWeight {
			randomAverageBestDeviation += statistic.AverageBestDeviation
			randomSuccessRate += statistic.SuccessRate
			if statistic.AverageBestDeviation < randomBestAverageBestDeviation {
				randomBestAverageBestDeviation = statistic.AverageBestDeviation
				randomBestAverageBestDeviationSeed = statistic.RandomSeed
			}
		}
		randomAverageBestDeviation /= float64(len(randomStatisticsForWeight))
		randomSuccessRate /= float64(len(randomStatisticsForWeight))

		rows = append(rows, randomSparseControlRow{
			instance:                           atspData.Name,
			msaAverageBestDeviation:            msaMetric.AverageMinDeviation,
			randomAverageBestDeviation:         randomAverageBestDeviation,
			randomBestAverageBestDeviation:     randomBestAverageBestDeviation,
			averageBestDeviationDelta:          msaMetric.AverageMinDeviation - randomAverageBestDeviation,
			msaSuccessRate:                     msaMetric.SuccessRate,
			randomSuccessRate:                  randomSuccessRate,
			successRateDelta:                   msaMetric.SuccessRate - randomSuccessRate,
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

func statisticsForHeuristicWeight(statistics []experiments.ExperimentsDataStatistics, heuristicWeight float64) (experiments.ExperimentsDataStatistics, bool) {
	for _, statistic := range statistics {
		if math.Abs(statistic.HeuristicWeight-heuristicWeight) < 1e-9 {
			return statistic, true
		}
	}

	return experiments.ExperimentsDataStatistics{}, false
}

func statisticsForHeuristicWeightAll(statistics []experiments.ExperimentsDataStatistics, heuristicWeight float64) []experiments.ExperimentsDataStatistics {
	rows := make([]experiments.ExperimentsDataStatistics, 0)
	for _, statistic := range statistics {
		if math.Abs(statistic.HeuristicWeight-heuristicWeight) < 1e-9 {
			rows = append(rows, statistic)
		}
	}

	return rows
}

func readEvaluationMsaHeuristicControlMetric(atspData project.AtspData, evaluationResultsRootPath string, config ControlReportsConfig) (EvaluationResultSummaryMetric, bool, error) {
	return readMsaControlMetric(atspData, evaluationResultsRootPath, config.StrictMsaHeuristic, config)
}

func readMsaControlMetric(atspData project.AtspData, evaluationResultsRootPath, referenceHeuristic string, config ControlReportsConfig) (EvaluationResultSummaryMetric, bool, error) {
	evaluationAtspData := project.WithExperimentOutputRoot(atspData, evaluationResultsRootPath)
	statistics, err := config.ReadHeuristicStatistics(evaluationAtspData.ResultFilePath)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return EvaluationResultSummaryMetric{}, false, nil
		}
		return EvaluationResultSummaryMetric{}, false, err
	}

	return msaHeuristicControlMetricForWeight(statistics, referenceHeuristic, config.EvaluationStrictMsaHeuristicWeight)
}

func msaHeuristicControlMetricForWeight(statistics []experiments.HeuristicExperimentStatistics, referenceHeuristic string, heuristicWeight float64) (EvaluationResultSummaryMetric, bool, error) {
	for _, statistic := range statistics {
		if statistic.Heuristic != referenceHeuristic || math.Abs(statistic.Statistics.HeuristicWeight-heuristicWeight) >= 1e-9 {
			continue
		}

		return EvaluationResultSummaryMetric{
			AverageMinDeviation:  statistic.Statistics.AverageBestDeviation,
			SuccessRate:          statistic.Statistics.SuccessRate,
			AverageBestIteration: statistic.Statistics.AverageBestAtIteration,
			HeuristicWeight:      statistic.Statistics.HeuristicWeight,
			Iterations:           statistic.Statistics.Iterations,
		}, true, nil
	}

	return EvaluationResultSummaryMetric{}, false, nil
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

	observed := max(wins, losses)
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
			evaluationResultsSummaryMetricCell(row.msaAverageBestDeviation, msaWins),
			evaluationResultsSummaryMetricCell(row.randomAverageBestDeviation, !msaWins),
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

func SaveDistanceRankedSparseControlReport(path string, atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath string, config ControlReportsConfig) (bool, error) {
	return saveDistanceRankedSparseControlReportForReference(path, atspsData, evaluationResultsRootPath, controlResultsRootPath, config.StrictMsaHeuristic, "Strict MSA", config.DistanceRankedSparseHeuristic, "distance-ranked sparse", config)
}

func saveDistanceRankedSparseControlReportForReference(path string, atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath, referenceHeuristic, referenceName, controlHeuristic, controlName string, config ControlReportsConfig) (bool, error) {
	rows, missingData, err := buildDistanceRankedSparseControlRows(atspsData, evaluationResultsRootPath, controlResultsRootPath, referenceHeuristic, controlHeuristic, config)
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
	fmt.Fprintf(&builder, "This sanity check compares %s against a deterministic sparse mask built from the cheapest directed edges. The control preserves the matching edge-count structure and %s.\n\n", referenceName, controlWeightDescription(config))
	writeDistanceRankedSparseControlFindings(&builder, rows, referenceName, controlName)
	builder.WriteString("\n")
	writeDistanceRankedSparseControlTable(&builder, rows, referenceName, controlName)
	if len(missingData) != 0 {
		writeDistanceRankedSparseControlMissingData(&builder, missingData)
	}

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildDistanceRankedSparseControlRows(atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath, referenceHeuristic, controlHeuristic string, config ControlReportsConfig) ([]distanceRankedSparseControlRow, []distanceRankedSparseControlMissingData, error) {
	rows := make([]distanceRankedSparseControlRow, 0, len(atspsData))
	missingData := make([]distanceRankedSparseControlMissingData, 0)

	for _, atspData := range atspsData {
		msaMetric, ok, err := readMsaControlMetric(atspData, evaluationResultsRootPath, referenceHeuristic, config)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.Name, reason: "missing evaluation " + referenceHeuristic + " result"})
			continue
		}

		controlAtspData := project.WithExperimentOutputRoot(atspData, controlResultsRootPath)
		controlStatistics, err := config.ReadStatistics(config.ResultFilePathForHeuristic(controlAtspData, controlHeuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.Name, reason: "missing evaluation-control " + controlHeuristic + " result CSV"})
				continue
			}
			return nil, nil, err
		}

		controlStatistic, ok := statisticsForHeuristicWeight(controlStatistics, msaMetric.HeuristicWeight)
		if !ok {
			missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.Name, reason: fmt.Sprintf("missing %s row for heuristic weight %.2f", controlHeuristic, msaMetric.HeuristicWeight)})
			continue
		}

		rows = append(rows, distanceRankedSparseControlRow{
			instance:                    atspData.Name,
			msaAverageBestDeviation:     msaMetric.AverageMinDeviation,
			controlAverageBestDeviation: controlStatistic.AverageBestDeviation,
			averageBestDeviationDelta:   msaMetric.AverageMinDeviation - controlStatistic.AverageBestDeviation,
			msaSuccessRate:              msaMetric.SuccessRate,
			controlSuccessRate:          controlStatistic.SuccessRate,
			successRateDelta:            msaMetric.SuccessRate - controlStatistic.SuccessRate,
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
			evaluationResultsSummaryMetricCell(row.msaAverageBestDeviation, msaWins),
			evaluationResultsSummaryMetricCell(row.controlAverageBestDeviation, !msaWins),
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

func saveSeededSparseControlReport(path string, atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath string, reportConfig seededSparseControlReportConfig, config ControlReportsConfig) (bool, error) {
	rows, missingData, err := buildSeededSparseControlRows(atspsData, evaluationResultsRootPath, controlResultsRootPath, reportConfig, config)
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
	fmt.Fprintf(&builder, "# %s\n\n", reportConfig.title)
	builder.WriteString(reportConfig.description)
	builder.WriteString("\n\n")
	writeSeededSparseControlFindings(&builder, rows, reportConfig)
	builder.WriteString("\n")
	writeSeededSparseControlTable(&builder, rows, reportConfig)
	if len(missingData) != 0 {
		writeSeededSparseControlMissingData(&builder, missingData)
	}

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildSeededSparseControlRows(atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath string, reportConfig seededSparseControlReportConfig, config ControlReportsConfig) ([]seededSparseControlRow, []seededSparseControlMissingData, error) {
	rows := make([]seededSparseControlRow, 0, len(atspsData))
	missingData := make([]seededSparseControlMissingData, 0)

	for _, atspData := range atspsData {
		msaMetric, ok, err := readMsaControlMetric(atspData, evaluationResultsRootPath, reportConfig.referenceHeuristic, config)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			missingData = append(missingData, seededSparseControlMissingData{instance: atspData.Name, reason: "missing evaluation " + reportConfig.referenceHeuristic + " result"})
			continue
		}

		controlAtspData := project.WithExperimentOutputRoot(atspData, controlResultsRootPath)
		controlStatistics, err := config.ReadStatistics(config.ResultFilePathForHeuristic(controlAtspData, reportConfig.controlHeuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, seededSparseControlMissingData{instance: atspData.Name, reason: "missing evaluation-control " + reportConfig.missingResultLabel})
				continue
			}
			return nil, nil, err
		}
		controlStatisticsForWeight := statisticsForHeuristicWeightAll(controlStatistics, msaMetric.HeuristicWeight)
		if len(controlStatisticsForWeight) == 0 {
			missingData = append(missingData, seededSparseControlMissingData{instance: atspData.Name, reason: fmt.Sprintf("missing %s rows for heuristic weight %.2f", reportConfig.missingResultLabel, msaMetric.HeuristicWeight)})
			continue
		}

		controlAverageBestDeviation := 0.0
		controlSuccessRate := 0.0
		controlBestAverageBestDeviation := math.Inf(1)
		controlBestAverageBestDeviationSeed := int64(0)
		for _, statistic := range controlStatisticsForWeight {
			controlAverageBestDeviation += statistic.AverageBestDeviation
			controlSuccessRate += statistic.SuccessRate
			if statistic.AverageBestDeviation < controlBestAverageBestDeviation {
				controlBestAverageBestDeviation = statistic.AverageBestDeviation
				controlBestAverageBestDeviationSeed = statistic.RandomSeed
			}
		}
		controlAverageBestDeviation /= float64(len(controlStatisticsForWeight))
		controlSuccessRate /= float64(len(controlStatisticsForWeight))

		rows = append(rows, seededSparseControlRow{
			instance:                            atspData.Name,
			msaAverageBestDeviation:             msaMetric.AverageMinDeviation,
			controlAverageBestDeviation:         controlAverageBestDeviation,
			controlBestAverageBestDeviation:     controlBestAverageBestDeviation,
			averageBestDeviationDelta:           msaMetric.AverageMinDeviation - controlAverageBestDeviation,
			msaSuccessRate:                      msaMetric.SuccessRate,
			controlSuccessRate:                  controlSuccessRate,
			successRateDelta:                    msaMetric.SuccessRate - controlSuccessRate,
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
			evaluationResultsSummaryMetricCell(row.msaAverageBestDeviation, msaWins),
			evaluationResultsSummaryMetricCell(row.controlAverageBestDeviation, !msaWins),
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

func SaveShuffledMsaControlReport(path string, atspsData []project.AtspData, evaluationResultsRootPath, controlResultsRootPath string, config ControlReportsConfig) (bool, error) {
	return saveSeededSparseControlReport(path, atspsData, evaluationResultsRootPath, controlResultsRootPath, seededSparseControlReportConfig{
		title:              "Shuffled MSA Control",
		description:        fmt.Sprintf("This sanity check compares Strict MSA against deterministic shuffles of the strict MSA mask. Each shuffle preserves the number and boost values of strict-MSA boosted directed edges, but assigns them to shuffled directed edges. The control %s.", controlWeightDescription(config)),
		referenceHeuristic: config.StrictMsaHeuristic,
		referenceName:      "Strict MSA",
		controlName:        "shuffled MSA",
		controlHeuristic:   config.ShuffledMsaHeuristic,
		missingResultLabel: "shuffled-MSA result CSV",
	}, config)
}

func controlWeightDescription(config ControlReportsConfig) string {
	return fmt.Sprintf("uses the same `heuristicWeight=%.2f`", config.EvaluationStrictMsaHeuristicWeight)
}
