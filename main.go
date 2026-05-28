package main

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"atsp_aco_msa/modules/algorithms/heuristics"
	"atsp_aco_msa/modules/algorithms/hungarian"
	"atsp_aco_msa/modules/algorithms/msaSupport"
	"atsp_aco_msa/modules/analysis/cycleCover"
	"atsp_aco_msa/modules/analysis/msaSupportTours"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/parsing"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
	"flag"
	"fmt"
	"html"
	"image/color"
	"math"
	"os"
	"path/filepath"
	"runtime/pprof"
	"slices"
	"sort"
	"strconv"
	"strings"
	"time"
)

type ExperimentsData struct {
	ExperimentParameters
	results []ExperimentResult
}

type ExperimentParameters struct {
	alpha, beta, rho, heuristicWeight float64
	iterations                        int
}

type ExperimentResult struct {
	bestAtIteration, threeOptImprovementsCount int
	bestTour                                   []int
	deviationPerIteration                      []float64
}

type ExperimentsDataStatistics struct {
	ExperimentParameters
	minBestAtIteration                                                               int
	averageBestAtIteration                                                           float64
	maxBestAtIteration                                                               int
	minThreeOptImprovementsCount                                                     int
	averageThreeOptImprovementsCount                                                 float64
	maxThreeOptImprovementsCount                                                     int
	minBestDeviation, averageBestDeviation, maxBestDeviation, successRate            float64
	minDeviationPerIteration, averageDeviationPerIteration, maxDeviationPerIteration []float64
}

type HeuristicExperimentStatistics struct {
	heuristic  string
	statistics ExperimentsDataStatistics
}

type finalResultsSummaryMetric struct {
	averageMinDeviation  float64
	successRate          float64
	averageBestIteration float64
	iterations           int
}

type finalResultsSummaryRow struct {
	instance string
	metrics  map[string]finalResultsSummaryMetric
}

const (
	instanceSetSmoke    = "smoke"
	instanceSetTiny     = "tiny"
	instanceSetBalanced = "balanced"
	instanceSetLarge    = "large"
	instanceSetAllKnown = "all-known"
)

const (
	runModeExperiment = "experiment"
	runModeAnalyze    = "analyze"
	runModeAll        = "all"
	runModeFinal      = "final"
	runModeFinal3Opt  = "final+3opt"
)

const (
	finalNumberOfExperiments       = 50
	finalMsaSupportWeight          = 0.9
	finalCycleCoverWeight          = 0.8
	defaultExperimentAlpha         = 1.0
	defaultExperimentBeta          = 2.0
	defaultExperimentRho           = 0.8
	defaultExperimentRunCount      = 30
	defaultBaselineHeuristicWeight = 0.0
)

const (
	heuristicBaseline   = "baseline"
	heuristicMsaSupport = "msa-support"
	heuristicCycleCover = "cycle-cover"
)

const finalHeuristicAll = "all"

var finalResultsSummaryHeuristics = []string{
	heuristicBaseline,
	heuristicMsaSupport,
	heuristicCycleCover,
}

var msaCountScalingCounts = []int{1, 2, 4, 8, 16, 32, 64, 0}

var smokeInstanceFiles = []string{
	"ftv64.atsp",
	"crane66_1.atsp",
	"crane66_2.atsp",
	"atex5.atsp",
	"ftv90.atsp",
}

var tinyInstanceFiles = []string{
	"atex1.atsp",
	"br17.atsp",
	"atex3.atsp",
	"ftv33.atsp",
	"ftv35.atsp",
	"ftv38.atsp",
	"p43.atsp",
	"ftv44.atsp",
	"atex4.atsp",
	"ftv47.atsp",
	"ry48p.atsp",
}

var balancedInstanceFiles = []string{
	"ft53.atsp",
	"ftv55.atsp",
	"ftv64.atsp",
	"crane66_0.atsp",
	"crane66_1.atsp",
	"crane66_2.atsp",
	"ft70.atsp",
	"ftv70.atsp",
	"atex5.atsp",
	"ftv90.atsp",
	"crane100_0.atsp",
	"crane100_1.atsp",
	"crane100_2.atsp",
	"ftv100.atsp",
	"td100_1.atsp",
	"ftv110.atsp",
	"dc112.atsp",
	"ftv120.atsp",
	"dc126.atsp",
	"ftv130.atsp",
	"dc134.atsp",
	"ftv140.atsp",
	"ftv150.atsp",
	"ftv160.atsp",
	"ftv170.atsp",
	"dc176.atsp",
	"dc188.atsp",
	"code198.atsp",
}

var largeInstanceFiles = []string{
	"rbg323.atsp",
	// "rbg358.atsp",
	// "rbg403.atsp",
	// "rbg443.atsp",
}

var statisticsCsvHeader = []string{
	"Alpha",
	"Beta",
	"Rho",
	"Heuristic weight",
	"Iterations",
	"Min best at iteration",
	"Avg best at iteration",
	"Max best at iteration",
	"Min local search improvements",
	"Avg local search improvements",
	"Max local search improvements",
	"Min best deviation",
	"Avg best deviation",
	"Max best deviation",
	"Success rate [%]",
}

func generateMarkdownCounts(paramName string, counts map[float64]int) string {
	markdown := fmt.Sprintf("### %s\n\n", paramName)
	markdown += "| Value | Count |\n"
	markdown += "|-------|-------|\n"

	// Sort the keys
	keys := make([]float64, 0, len(counts))
	for key := range counts {
		keys = append(keys, key)
	}
	sort.Float64s(keys)

	// Add rows
	for _, key := range keys {
		markdown += fmt.Sprintf("| %.2f | %d |\n", key, counts[key])
	}

	markdown += "\n"
	return markdown
}

func saveBestParametersInfo(resultsRootPath, fileName string, bestStatistics []ExperimentsDataStatistics) {
	sort.SliceStable(bestStatistics, func(i, j int) bool {
		if bestStatistics[i].averageBestDeviation != bestStatistics[j].averageBestDeviation {
			return bestStatistics[i].averageBestDeviation < bestStatistics[j].averageBestDeviation
		}

		if bestStatistics[i].successRate != bestStatistics[j].successRate {
			return bestStatistics[i].successRate > bestStatistics[j].successRate
		}

		return bestStatistics[i].averageBestAtIteration < bestStatistics[j].averageBestAtIteration
	})

	uniqueParameters := map[ExperimentParameters]int{}

	for _, statistic := range bestStatistics {
		parameters := statistic.ExperimentParameters
		parameters.iterations = 0
		uniqueParameters[parameters]++
	}

	alphaCounts := make(map[float64]int)
	betaCounts := make(map[float64]int)
	rhoCounts := make(map[float64]int)
	heuristicWeightCounts := make(map[float64]int)

	// Count the occurrences
	for params := range uniqueParameters {
		alphaCounts[params.alpha]++
		betaCounts[params.beta]++
		rhoCounts[params.rho]++
		heuristicWeightCounts[params.heuristicWeight]++
	}

	// Markdown content
	markdown := "# Best Parameters Report\n\n"

	// Unique combinations count
	markdown += fmt.Sprintf("Found **%d** best unique parameter combinations.\n\n", len(uniqueParameters))

	// Best parameters
	markdown += "## Best Parameters\n\n"
	markdown += "| Alpha | Beta | Rho | Heuristic weight | Times used |\n"
	markdown += "|-------|------|-----|-------|------------|\n"

	sortedParameters := make([]ExperimentParameters, 0, len(uniqueParameters))
	for parameters := range uniqueParameters {
		sortedParameters = append(sortedParameters, parameters)
	}
	sort.Slice(sortedParameters, func(i, j int) bool {
		left, right := sortedParameters[i], sortedParameters[j]
		if left.alpha != right.alpha {
			return left.alpha < right.alpha
		}
		if left.beta != right.beta {
			return left.beta < right.beta
		}
		if left.rho != right.rho {
			return left.rho < right.rho
		}
		return left.heuristicWeight < right.heuristicWeight
	})

	for _, parameters := range sortedParameters {
		timesUsed := uniqueParameters[parameters]
		markdown += fmt.Sprintf("| %.2f | %.2f | %.2f | %.2f | %d |\n",
			parameters.alpha, parameters.beta, parameters.rho, parameters.heuristicWeight, timesUsed)
	}
	markdown += "\n"

	// Parameter occurrences
	markdown += "## Parameter Values Occurrences\n\n"
	markdown += generateMarkdownCounts("Alpha", alphaCounts)
	markdown += generateMarkdownCounts("Beta", betaCounts)
	markdown += generateMarkdownCounts("Rho", rhoCounts)
	markdown += generateMarkdownCounts("Heuristic weight", heuristicWeightCounts)

	// Parameter ranges
	markdown += "## Parameter Ranges\n\n"
	minAlpha, maxAlpha := findMinMax(alphaCounts)
	minBeta, maxBeta := findMinMax(betaCounts)
	minRho, maxRho := findMinMax(rhoCounts)
	minHeuristicWeight, maxHeuristicWeight := findMinMax(heuristicWeightCounts)

	markdown += fmt.Sprintf("- **Alpha**: %.2f - %.2f\n", minAlpha, maxAlpha)
	markdown += fmt.Sprintf("- **Beta**: %.2f - %.2f\n", minBeta, maxBeta)
	markdown += fmt.Sprintf("- **Rho**: %.2f - %.2f\n", minRho, maxRho)
	markdown += fmt.Sprintf("- **Heuristic weight**: %.2f - %.2f\n", minHeuristicWeight, maxHeuristicWeight)

	// Save to a file
	reportPath := filepath.Join(resultsRootPath, fileName)
	if err := os.MkdirAll(filepath.Dir(reportPath), 0700); err != nil {
		fmt.Printf("Failed to create report directory: %v\n", err)
		return
	}

	err := os.WriteFile(reportPath, []byte(markdown), 0644)
	if err != nil {
		fmt.Printf("Failed to save report: %v\n", err)
	}
}

func findMinMax(counts map[float64]int) (float64, float64) {
	min := math.MaxFloat64
	max := -math.MaxFloat64

	for value := range counts {
		if value < min {
			min = value
		}
		if value > max {
			max = value
		}
	}

	return min, max
}

func readStatistics(csvFilePath string) ([]ExperimentsDataStatistics, error) {
	file, err := os.Open(csvFilePath)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)

	// Read the header and skip it
	header, err := reader.Read() // Read the first line (header)
	if err != nil {
		return nil, fmt.Errorf("failed to read header: %w", err)
	}
	if len(header) != len(statisticsCsvHeader) {
		return nil, fmt.Errorf("invalid header length: %d", len(header))
	}

	var statistics []ExperimentsDataStatistics

	// Parse the remaining rows
	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("failed to read records: %w", err)
	}

	for _, record := range records {
		if len(record) != len(statisticsCsvHeader) {
			return nil, fmt.Errorf("invalid record length: %d", len(record))
		}

		statistic, err := parseStatisticsRecord(record)
		if err != nil {
			return nil, err
		}

		statistics = append(statistics, statistic)
	}

	return statistics, nil
}

func parseStatisticsRecord(record []string) (ExperimentsDataStatistics, error) {
	if len(record) != len(statisticsCsvHeader) {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid record length: %d", len(record))
	}

	alpha, err := strconv.ParseFloat(record[0], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid alpha: %w", err)
	}
	beta, err := strconv.ParseFloat(record[1], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid beta: %w", err)
	}
	rho, err := strconv.ParseFloat(record[2], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid rho: %w", err)
	}
	heuristicWeight, err := strconv.ParseFloat(record[3], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid heuristic weight: %w", err)
	}
	iterations, err := strconv.Atoi(record[4])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid iterations: %w", err)
	}
	minBestAtIteration, err := strconv.Atoi(record[5])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min best at iteration: %w", err)
	}
	averageBestAtIteration, err := strconv.ParseFloat(record[6], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg best at iteration: %w", err)
	}
	maxBestAtIteration, err := strconv.Atoi(record[7])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max best at iteration: %w", err)
	}
	minThreeOptImprovementsCount, err := strconv.Atoi(record[8])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min local search improvements: %w", err)
	}
	averageThreeOptImprovementsCount, err := strconv.ParseFloat(record[9], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg local search improvements: %w", err)
	}
	maxThreeOptImprovementsCount, err := strconv.Atoi(record[10])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max local search improvements: %w", err)
	}
	minBestDeviation, err := strconv.ParseFloat(record[11], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min best deviation: %w", err)
	}
	averageBestDeviation, err := strconv.ParseFloat(record[12], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg best deviation: %w", err)
	}
	maxBestDeviation, err := strconv.ParseFloat(record[13], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max best deviation: %w", err)
	}
	successRate, err := strconv.ParseFloat(record[14], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid success rate: %w", err)
	}

	return ExperimentsDataStatistics{
		ExperimentParameters: ExperimentParameters{
			alpha:           alpha,
			beta:            beta,
			rho:             rho,
			heuristicWeight: heuristicWeight,
			iterations:      iterations,
		},
		minBestAtIteration:               minBestAtIteration,
		averageBestAtIteration:           averageBestAtIteration,
		maxBestAtIteration:               maxBestAtIteration,
		minThreeOptImprovementsCount:     minThreeOptImprovementsCount,
		averageThreeOptImprovementsCount: averageThreeOptImprovementsCount,
		maxThreeOptImprovementsCount:     maxThreeOptImprovementsCount,
		minBestDeviation:                 minBestDeviation,
		averageBestDeviation:             averageBestDeviation,
		maxBestDeviation:                 maxBestDeviation,
		successRate:                      successRate,
	}, nil
}

func saveStatistics(resultCsvPath string, statistics []ExperimentsDataStatistics) {
	if err := os.MkdirAll(filepath.Dir(resultCsvPath), 0700); err != nil {
		fmt.Printf("Failed to create statistics directory: %v\n", err)
		return
	}

	file, err := os.Create(resultCsvPath)
	if err != nil {
		fmt.Printf("Failed to save statistics: %v\n", err)
		return
	}
	defer file.Close()

	writer := csv.NewWriter(file)

	_ = writer.Write(statisticsCsvHeader)

	for _, statistic := range statistics {
		writer.Write(statisticsCsvRecord(statistic))
	}

	writer.Flush()
}

func saveHeuristicStatistics(resultCsvPath string, statistics []HeuristicExperimentStatistics) error {
	if err := os.MkdirAll(filepath.Dir(resultCsvPath), 0700); err != nil {
		return err
	}

	file, err := os.Create(resultCsvPath)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	header := append([]string{"Heuristic"}, statisticsCsvHeader...)
	if err := writer.Write(header); err != nil {
		return err
	}

	sortedStatistics := append([]HeuristicExperimentStatistics(nil), statistics...)
	sort.SliceStable(sortedStatistics, func(i, j int) bool {
		left, right := sortedStatistics[i].statistics, sortedStatistics[j].statistics
		if left.averageBestDeviation != right.averageBestDeviation {
			return left.averageBestDeviation < right.averageBestDeviation
		}
		if left.successRate != right.successRate {
			return left.successRate > right.successRate
		}
		return sortedStatistics[i].heuristic < sortedStatistics[j].heuristic
	})

	for _, statistic := range sortedStatistics {
		record := append([]string{statistic.heuristic}, statisticsCsvRecord(statistic.statistics)...)
		if err := writer.Write(record); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func saveFinalHeuristicStatistics(resultCsvPath string, statistics []HeuristicExperimentStatistics, configurations []finalExperimentConfiguration) error {
	if len(configurations) == len(finalExperimentConfigurations()) {
		return saveHeuristicStatistics(resultCsvPath, statistics)
	}

	existingStatistics, err := readHeuristicStatistics(resultCsvPath)
	if err != nil {
		if os.IsNotExist(err) {
			existingStatistics = nil
		} else {
			return err
		}
	}

	mergedStatistics := mergeHeuristicStatistics(existingStatistics, statistics, configurations)
	return saveHeuristicStatistics(resultCsvPath, mergedStatistics)
}

func readHeuristicStatistics(resultCsvPath string) ([]HeuristicExperimentStatistics, error) {
	file, err := os.Open(resultCsvPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		return nil, err
	}

	expectedHeader := append([]string{"Heuristic"}, statisticsCsvHeader...)
	if !slices.Equal(header, expectedHeader) {
		return nil, fmt.Errorf("invalid heuristic statistics header")
	}

	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	statistics := make([]HeuristicExperimentStatistics, 0, len(records))
	for _, record := range records {
		if len(record) != len(expectedHeader) {
			return nil, fmt.Errorf("invalid heuristic statistics record length: got %d want %d", len(record), len(expectedHeader))
		}

		parsed, err := parseStatisticsRecord(record[1:])
		if err != nil {
			return nil, fmt.Errorf("%s: %w", record[0], err)
		}

		statistics = append(statistics, HeuristicExperimentStatistics{
			heuristic:  record[0],
			statistics: parsed,
		})
	}

	return statistics, nil
}

func mergeHeuristicStatistics(existingStatistics, newStatistics []HeuristicExperimentStatistics, configurations []finalExperimentConfiguration) []HeuristicExperimentStatistics {
	replaceSet := make(map[string]struct{}, len(configurations))
	for _, configuration := range configurations {
		replaceSet[configuration.heuristic] = struct{}{}
	}

	knownFinalHeuristics := make(map[string]struct{}, len(finalExperimentConfigurations()))
	for _, configuration := range finalExperimentConfigurations() {
		knownFinalHeuristics[configuration.heuristic] = struct{}{}
	}

	merged := make([]HeuristicExperimentStatistics, 0, len(existingStatistics)+len(newStatistics))
	for _, statistic := range existingStatistics {
		if _, known := knownFinalHeuristics[statistic.heuristic]; !known {
			continue
		}
		if _, replace := replaceSet[statistic.heuristic]; !replace {
			merged = append(merged, statistic)
		}
	}

	merged = append(merged, newStatistics...)
	return merged
}

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
		metrics, err := readFinalResultSummaryMetrics(atspData.resultFilePath)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read final result metrics: %w", atspData.name, err)
		}

		rows = append(rows, finalResultsSummaryRow{instance: atspData.name, metrics: metrics})
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
			total.averageMinDeviation += metric.averageMinDeviation
			total.successRate += metric.successRate
			total.averageBestIteration += metric.averageBestIteration
			total.iterations += metric.iterations
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
			averageMinDeviation:  total.averageMinDeviation / float64(count),
			successRate:          total.successRate / float64(count),
			averageBestIteration: total.averageBestIteration / float64(count),
			iterations:           int(math.Round(float64(total.iterations) / float64(count))),
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
		bestDeviation := averageMetrics[bestDeviationHeuristics[0].heuristic].averageMinDeviation
		fmt.Fprintf(builder, "- **%s has the lowest average best deviation overall: %.2f%%.**\n",
			joinHeuristicDisplayNames(bestDeviationHeuristics),
			bestDeviation)
	}
	if len(bestSuccessHeuristics) != 0 {
		bestSuccess := averageMetrics[bestSuccessHeuristics[0].heuristic].successRate
		fmt.Fprintf(builder, "- **%s has the highest average success rate overall: %.2f%%.**\n",
			joinHeuristicDisplayNames(bestSuccessHeuristics),
			bestSuccess)
	}

	fmt.Fprintf(builder, "- **Best-or-tied average best deviation counts: Baseline %d/%d, MSA support %d/%d, Cycle cover %d/%d.**\n",
		deviationWinCounts[heuristicBaseline], len(rows),
		deviationWinCounts[heuristicMsaSupport], len(rows),
		deviationWinCounts[heuristicCycleCover], len(rows))
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
			finalResultsSummaryMetricCell(metric.averageMinDeviation, deviationHighlights[heuristic]),
			finalResultsSummaryMetricCell(metric.successRate, successHighlights[heuristic]))
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

		if metric.averageMinDeviation < minDeviation {
			minDeviation = metric.averageMinDeviation
		}
		if metric.successRate > maxSuccess {
			maxSuccess = metric.successRate
		}
	}

	for _, heuristic := range finalResultsSummaryHeuristics {
		metric, ok := metrics[heuristic]
		if !ok {
			continue
		}

		if math.Abs(metric.averageMinDeviation-minDeviation) < 1e-9 {
			deviationHighlights[heuristic] = true
		}
		if maxSuccess > 0 && math.Abs(metric.successRate-maxSuccess) < 1e-9 {
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

	comparisons := []finalPairwisePerformanceComparison{
		calculateFinalPairwisePerformanceComparison(rows, heuristicMsaSupport, heuristicBaseline),
		calculateFinalPairwisePerformanceComparison(rows, heuristicCycleCover, heuristicBaseline),
		calculateFinalPairwisePerformanceComparison(rows, heuristicMsaSupport, heuristicCycleCover),
	}

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
		deviationDelta := leftMetric.averageMinDeviation - rightMetric.averageMinDeviation
		comparison.averageBestDeviationDelta += deviationDelta
		comparison.successRateDelta += leftMetric.successRate - rightMetric.successRate

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
		"- **Average best-iteration position: Baseline %.2f%%, MSA support %.2f%%, cycle cover %.2f%%.**\n",
		averages[heuristicBaseline],
		averages[heuristicMsaSupport],
		averages[heuristicCycleCover])
	fmt.Fprintf(builder,
		"- **Earliest-or-tied convergence counts: Baseline %d/%d, MSA support %d/%d, cycle cover %d/%d.**\n",
		winCounts[heuristicBaseline], len(rows),
		winCounts[heuristicMsaSupport], len(rows),
		winCounts[heuristicCycleCover], len(rows))
}

func writeFinalConvergenceTable(builder *strings.Builder, rows []finalResultsSummaryRow, averages map[string]float64) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>Baseline best iter [%]</th><th>MSA support best iter [%]</th><th>Cycle cover best iter [%]</th></tr>\n")
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
		fmt.Fprintf(builder, "<td align=\"right\">%s</td>",
			finalResultsSummaryMetricCell(averages[heuristic], averageHighlights[heuristic]))
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
			if metric.iterations <= 0 {
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
		if !ok || metric.iterations <= 0 {
			continue
		}
		values[heuristic] = convergencePercent(metric)
	}

	return lowestFloatHighlights(values, false)
}

func convergencePercent(metric finalResultsSummaryMetric) float64 {
	if metric.iterations <= 0 {
		return 0
	}

	return 100.0 * metric.averageBestIteration / float64(metric.iterations)
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
		"- **Reduced 3-opt lowers average best deviation for every variant: Baseline %s pp, MSA support %s pp, cycle cover %s pp.**\n",
		formatSignedFloat(finalAverages[heuristicBaseline].averageMinDeviation-finalThreeOptAverages[heuristicBaseline].averageMinDeviation),
		formatSignedFloat(finalAverages[heuristicMsaSupport].averageMinDeviation-finalThreeOptAverages[heuristicMsaSupport].averageMinDeviation),
		formatSignedFloat(finalAverages[heuristicCycleCover].averageMinDeviation-finalThreeOptAverages[heuristicCycleCover].averageMinDeviation))
	fmt.Fprintf(builder,
		"- **Reduced 3-opt increases success rate for every variant: Baseline %s pp, MSA support %s pp, cycle cover %s pp.**\n",
		formatSignedFloat(finalThreeOptAverages[heuristicBaseline].successRate-finalAverages[heuristicBaseline].successRate),
		formatSignedFloat(finalThreeOptAverages[heuristicMsaSupport].successRate-finalAverages[heuristicMsaSupport].successRate),
		formatSignedFloat(finalThreeOptAverages[heuristicCycleCover].successRate-finalAverages[heuristicCycleCover].successRate))

	msaGainWithout := deviationGainVsBaseline(finalAverages, heuristicMsaSupport)
	msaGainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristicMsaSupport)
	cycleGainWithout := deviationGainVsBaseline(finalAverages, heuristicCycleCover)
	cycleGainWith := deviationGainVsBaseline(finalThreeOptAverages, heuristicCycleCover)
	fmt.Fprintf(builder,
		"- **Deviation gain over baseline shrinks with 3-opt: MSA support %s -> %s pp, cycle cover %s -> %s pp.**\n",
		formatSignedFloat(msaGainWithout),
		formatSignedFloat(msaGainWith),
		formatSignedFloat(cycleGainWithout),
		formatSignedFloat(cycleGainWith))
	fmt.Fprintf(builder,
		"- **Only %.2f%% of the MSA-support deviation gain and %.2f%% of the cycle-cover deviation gain remains visible after enabling 3-opt.**\n",
		signalRemainingPercent(msaGainWithout, msaGainWith),
		signalRemainingPercent(cycleGainWithout, cycleGainWith))
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
		without := finalAverages[heuristic]
		with := finalThreeOptAverages[heuristic]
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(heuristicDisplayName(heuristic)),
			without.averageMinDeviation,
			with.averageMinDeviation,
			without.successRate,
			with.successRate)
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
	for _, heuristic := range []string{heuristicMsaSupport, heuristicCycleCover} {
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

func deviationGainVsBaseline(averages map[string]finalResultsSummaryMetric, heuristic string) float64 {
	return averages[heuristicBaseline].averageMinDeviation - averages[heuristic].averageMinDeviation
}

func signalRemainingPercent(without, with float64) float64 {
	if math.Abs(without) < 1e-9 {
		return 0
	}

	return 100.0 * math.Abs(with) / math.Abs(without)
}

func saveStructuralPerformanceLinkReport(path string, rows []finalResultsSummaryRow, analyses []cycleCover.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	structuralRows := filterAnalysesWithFoundOptimalEdges(sortedCycleCoverAnalyses(analyses))
	structuralTotals := structuralSimilarityTotals(structuralRows)
	allowedInstances := make(map[string]struct{}, len(structuralRows))
	for _, row := range structuralRows {
		allowedInstances[row.Instance] = struct{}{}
	}

	performance := averagePerformanceByHeuristic(rows, allowedInstances)
	msaPrecision := ratio(structuralTotals.msaOptimalEdges, structuralTotals.msaEdges)
	cycleCoverPrecision := ratio(structuralTotals.cycleCoverOptimalEdges, structuralTotals.cycleCoverEdges)
	msaRecall := ratio(structuralTotals.msaOptimalEdges, structuralTotals.foundOptimalEdges)
	cycleCoverRecall := ratio(structuralTotals.cycleCoverOptimalEdges, structuralTotals.foundOptimalEdges)

	precisionHighlights := highestFloatHighlights(map[string]float64{
		heuristicMsaSupport: msaPrecision,
		heuristicCycleCover: cycleCoverPrecision,
	})
	recallHighlights := highestFloatHighlights(map[string]float64{
		heuristicMsaSupport: msaRecall,
		heuristicCycleCover: cycleCoverRecall,
	})
	deviationHighlights := lowestFloatHighlights(map[string]float64{
		heuristicMsaSupport: performance[heuristicMsaSupport].averageMinDeviation,
		heuristicCycleCover: performance[heuristicCycleCover].averageMinDeviation,
	}, true)
	successHighlights := highestFloatHighlights(map[string]float64{
		heuristicMsaSupport: performance[heuristicMsaSupport].successRate,
		heuristicCycleCover: performance[heuristicCycleCover].successRate,
	})

	var builder strings.Builder
	builder.WriteString("# Structural Similarity And Performance\n\n")
	builder.WriteString("This table links structural similarity to found optimal tours with final MMAS performance. Both structural and performance values are computed only for instances with at least one found optimal tour.\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Heuristic</th><th>Structural precision [%]</th><th>Structural recall [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	writeStructuralPerformanceLinkRow(&builder, heuristicMsaSupport, msaPrecision, msaRecall, performance[heuristicMsaSupport], precisionHighlights, recallHighlights, deviationHighlights, successHighlights)
	writeStructuralPerformanceLinkRow(&builder, heuristicCycleCover, cycleCoverPrecision, cycleCoverRecall, performance[heuristicCycleCover], precisionHighlights, recallHighlights, deviationHighlights, successHighlights)
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
		finalResultsSummaryMetricCell(performance.averageMinDeviation, deviationHighlights[heuristic]),
		finalResultsSummaryMetricCell(performance.successRate, successHighlights[heuristic]))
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
			total.averageMinDeviation += metric.averageMinDeviation
			total.successRate += metric.successRate
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
			averageMinDeviation: total.averageMinDeviation / float64(count),
			successRate:         total.successRate / float64(count),
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

func saveStructuralSimilarityReport(path string, analyses []cycleCover.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows := filterAnalysesWithFoundOptimalEdges(sortedCycleCoverAnalyses(analyses))
	totals := structuralSimilarityTotals(rows)

	var builder strings.Builder
	builder.WriteString("# Structural Similarity To Found Optimal Tours\n\n")
	builder.WriteString("This table compares the current MSA support edge set and the minimum cycle-cover edge set against the found optimal tours saved in `solutions.csv`.\n\n")
	builder.WriteString("Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.\n\n")
	writeStructuralSimilarityFindings(&builder, totals)
	builder.WriteString("\n")
	writeStructuralSimilarityTable(&builder, rows, totals)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeStructuralSimilarityFindings(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **Precision vs found-optimal tours: MSA support %.2f%%, cycle cover %.2f%%.**\n", 100*msaPrecision, 100*cycleCoverPrecision)
	fmt.Fprintf(builder, "- **Recall vs found-optimal tours: MSA support %.2f%%, cycle cover %.2f%%.**\n", 100*msaRecall, 100*cycleCoverRecall)
	fmt.Fprintf(builder, "- **Best-or-tied precision counts: MSA support %d/%d, cycle cover %d/%d.**\n",
		totals.msaPrecisionWins,
		totals.instanceCount,
		totals.cycleCoverPrecisionWins,
		totals.instanceCount)
	fmt.Fprintf(builder, "- **Best-or-tied recall counts: MSA support %d/%d, cycle cover %d/%d.**\n",
		totals.msaRecallWins,
		totals.instanceCount,
		totals.cycleCoverRecallWins,
		totals.instanceCount)
}

func writeStructuralSimilarityTable(builder *strings.Builder, rows []cycleCover.InstanceAnalysis, totals structuralSimilaritySummary) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">MSA support</th><th colspan=\"2\">Cycle cover</th></tr>\n")
	builder.WriteString("<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, analysis := range rows {
		writeStructuralSimilarityRow(builder, analysis)
	}
	writeStructuralSimilarityTotalRow(builder, totals)

	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeStructuralSimilarityRow(builder *strings.Builder, analysis cycleCover.InstanceAnalysis) {
	metrics := analysis.Metrics
	msaMetrics := metrics.HighMsaSupportMetrics
	cycleCoverMetrics := metrics.CycleCoverMetrics
	msaPrecisionBold, cycleCoverPrecisionBold := bestPositivePair(msaMetrics.Precision, cycleCoverMetrics.Precision)
	msaRecallBold, cycleCoverRecallBold := bestPositivePair(msaMetrics.Recall, cycleCoverMetrics.Recall)

	writeStructuralSimilarityTableRow(
		builder,
		html.EscapeString(analysis.Instance),
		msaMetrics.Precision,
		msaMetrics.Recall,
		cycleCoverMetrics.Precision,
		cycleCoverMetrics.Recall,
		msaPrecisionBold,
		msaRecallBold,
		cycleCoverPrecisionBold,
		cycleCoverRecallBold)
}

func writeStructuralSimilarityTotalRow(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)
	msaPrecisionBold, cycleCoverPrecisionBold := bestPositivePair(msaPrecision, cycleCoverPrecision)
	msaRecallBold, cycleCoverRecallBold := bestPositivePair(msaRecall, cycleCoverRecall)

	writeStructuralSimilarityTableRow(
		builder,
		"<strong>Total</strong>",
		msaPrecision,
		msaRecall,
		cycleCoverPrecision,
		cycleCoverRecall,
		msaPrecisionBold,
		msaRecallBold,
		cycleCoverPrecisionBold,
		cycleCoverRecallBold)
}

func writeStructuralSimilarityTableRow(builder *strings.Builder, instanceCell string, msaPrecision, msaRecall, cycleCoverPrecision, cycleCoverRecall float64, msaPrecisionBold, msaRecallBold, cycleCoverPrecisionBold, cycleCoverRecallBold bool) {
	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		instanceCell,
		finalResultsSummaryMetricCell(100*msaPrecision, msaPrecisionBold),
		finalResultsSummaryMetricCell(100*msaRecall, msaRecallBold),
		finalResultsSummaryMetricCell(100*cycleCoverPrecision, cycleCoverPrecisionBold),
		finalResultsSummaryMetricCell(100*cycleCoverRecall, cycleCoverRecallBold))
}

type structuralSimilaritySummary struct {
	instanceCount           int
	foundOptimalTours       int
	foundOptimalEdges       int
	tourEdges               int
	msaEdges                int
	msaOptimalEdges         int
	cycleCoverEdges         int
	cycleCoverOptimalEdges  int
	msaPrecisionWins        int
	cycleCoverPrecisionWins int
	msaRecallWins           int
	cycleCoverRecallWins    int
}

func structuralSimilarityTotals(rows []cycleCover.InstanceAnalysis) structuralSimilaritySummary {
	var totals structuralSimilaritySummary
	totals.instanceCount = len(rows)

	for _, analysis := range rows {
		metrics := analysis.Metrics
		msaMetrics := metrics.HighMsaSupportMetrics
		cycleCoverMetrics := metrics.CycleCoverMetrics
		msaPrecisionWins, cycleCoverPrecisionWins := bestPositivePair(msaMetrics.Precision, cycleCoverMetrics.Precision)
		msaRecallWins, cycleCoverRecallWins := bestPositivePair(msaMetrics.Recall, cycleCoverMetrics.Recall)

		totals.foundOptimalTours += metrics.FoundOptimalTourCount
		totals.foundOptimalEdges += metrics.UniqueFoundOptimalEdgeCount
		totals.tourEdges += analysis.Dimension
		totals.msaEdges += msaMetrics.EdgeCount
		totals.msaOptimalEdges += msaMetrics.OptimalEdgeCount
		totals.cycleCoverEdges += cycleCoverMetrics.EdgeCount
		totals.cycleCoverOptimalEdges += cycleCoverMetrics.OptimalEdgeCount
		if msaPrecisionWins {
			totals.msaPrecisionWins++
		}
		if cycleCoverPrecisionWins {
			totals.cycleCoverPrecisionWins++
		}
		if msaRecallWins {
			totals.msaRecallWins++
		}
		if cycleCoverRecallWins {
			totals.cycleCoverRecallWins++
		}
	}

	return totals
}

func saveMsaSupportCycleCoverOverlapReport(path string, analyses []cycleCover.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows := sortedCycleCoverAnalyses(analyses)
	totals := msaSupportCycleCoverOverlapTotals(rows)

	var builder strings.Builder
	builder.WriteString("# MSA Support And Cycle-Cover Overlap\n\n")
	builder.WriteString("This table compares the current MSA support edge set directly with the minimum cycle-cover edge set.\n\n")
	writeMsaSupportCycleCoverOverlapFindings(&builder, totals)
	builder.WriteString("\n")
	writeMsaSupportCycleCoverOverlapTable(&builder, rows, totals)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeMsaSupportCycleCoverOverlapFindings(builder *strings.Builder, totals msaSupportCycleCoverOverlapSummary) {
	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **%.2f%% of MSA support edges are also cycle-cover edges.**\n", 100*ratio(totals.sharedEdges, totals.msaEdges))
	fmt.Fprintf(builder, "- **%.2f%% of cycle-cover edges are also MSA support edges.**\n", 100*ratio(totals.sharedEdges, totals.cycleCoverEdges))
	fmt.Fprintf(builder,
		"- **Found-optimal edge partition: both %d, only MSA support %d, only cycle cover %d, neither %d.**\n",
		totals.optimalBoth,
		totals.optimalOnlyMsa,
		totals.optimalOnlyCycleCover,
		totals.optimalNeither)
}

func writeMsaSupportCycleCoverOverlapTable(builder *strings.Builder, rows []cycleCover.InstanceAnalysis, totals msaSupportCycleCoverOverlapSummary) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, analysis := range rows {
		writeMsaSupportCycleCoverOverlapRow(builder, analysis)
	}
	writeMsaSupportCycleCoverOverlapTotalRow(builder, totals)

	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeMsaSupportCycleCoverOverlapRow(builder *strings.Builder, analysis cycleCover.InstanceAnalysis) {
	metrics := analysis.Metrics
	msaEdges := metrics.HighMsaSupportMetrics.EdgeCount
	cycleCoverEdges := metrics.CycleCoverMetrics.EdgeCount
	sharedEdges := metrics.CycleCoverEdgesWithHighMsaSupport

	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td></tr>\n",
		html.EscapeString(analysis.Instance),
		100*ratio(sharedEdges, msaEdges),
		100*ratio(sharedEdges, cycleCoverEdges),
		metrics.OptimalEdgesInCycleCoverAndHighMsaSupport,
		metrics.OptimalEdgesInHighMsaSupportNotCycleCover,
		metrics.OptimalEdgesInCycleCoverNotHighMsaSupport)
}

func writeMsaSupportCycleCoverOverlapTotalRow(builder *strings.Builder, totals msaSupportCycleCoverOverlapSummary) {
	fmt.Fprintf(builder,
		"<tr><td><strong>Total</strong></td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td></tr>\n",
		100*ratio(totals.sharedEdges, totals.msaEdges),
		100*ratio(totals.sharedEdges, totals.cycleCoverEdges),
		totals.optimalBoth,
		totals.optimalOnlyMsa,
		totals.optimalOnlyCycleCover)
}

type msaSupportCycleCoverOverlapSummary struct {
	msaEdges              int
	cycleCoverEdges       int
	sharedEdges           int
	optimalBoth           int
	optimalOnlyMsa        int
	optimalOnlyCycleCover int
	optimalNeither        int
}

func msaSupportCycleCoverOverlapTotals(rows []cycleCover.InstanceAnalysis) msaSupportCycleCoverOverlapSummary {
	var totals msaSupportCycleCoverOverlapSummary
	for _, analysis := range rows {
		metrics := analysis.Metrics
		totals.msaEdges += metrics.HighMsaSupportMetrics.EdgeCount
		totals.cycleCoverEdges += metrics.CycleCoverMetrics.EdgeCount
		totals.sharedEdges += metrics.CycleCoverEdgesWithHighMsaSupport
		totals.optimalBoth += metrics.OptimalEdgesInCycleCoverAndHighMsaSupport
		totals.optimalOnlyMsa += metrics.OptimalEdgesInHighMsaSupportNotCycleCover
		totals.optimalOnlyCycleCover += metrics.OptimalEdgesInCycleCoverNotHighMsaSupport
		totals.optimalNeither += metrics.OptimalEdgesInNeitherCycleCoverNorHigh
	}

	return totals
}

func saveMsaCountScalingReport(path string, atspsData []AtspData) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows, err := buildMsaCountScalingRows(atspsData, msaCountScalingCounts)
	if err != nil {
		return err
	}

	var builder strings.Builder
	builder.WriteString("# MSA Count Scaling\n\n")
	builder.WriteString("This table shows how the MSA-support signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.\n\n")
	if len(rows) > 0 {
		fmt.Fprintf(&builder, "The table uses pooled totals over %d selected instances. Precision and recall use the %d instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.\n\n", rows[0].instanceCount, rows[0].referenceInstanceCount)
	}
	writeMsaCountScalingFindings(&builder, rows)
	builder.WriteString("\n")
	writeMsaCountScalingTable(&builder, rows)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildMsaCountScalingRows(atspsData []AtspData, requestedCounts []int) ([]msaCountScalingRow, error) {
	rows := make([]msaCountScalingRow, 0, len(requestedCounts))
	for _, requestedCount := range requestedCounts {
		row := msaCountScalingRow{requestedCount: requestedCount}

		for _, atspData := range atspsData {
			msas, err := readOrCreateIndividualMsas(atspData)
			if err != nil {
				return nil, err
			}
			if len(msas) == 0 {
				continue
			}

			selectedRootCount := requestedCount
			if selectedRootCount == 0 || selectedRootCount > len(msas) {
				selectedRootCount = len(msas)
			}
			selectedRoots := selectEvenlySpacedRootIndexes(len(msas), selectedRootCount)
			boostedEdges := buildPartialMsaSupportEdgeSet(msas, selectedRoots)

			row.instanceCount++
			row.boostedEdges += len(boostedEdges)
			row.boostedTargetEdges += maxIntValue(len(msas)-1, 0)

			tours, err := msaSupportTours.ReadOptimalTours(atspData.optimalUniqueToursCsvPath)
			if err != nil {
				return nil, fmt.Errorf("%s: read found optimal tours: %w", atspData.name, err)
			}
			optimalEdges := buildAnalysisTourEdgeSet(tours)
			if len(optimalEdges) == 0 {
				continue
			}

			overlap := countEdgeSetIntersection(boostedEdges, optimalEdges)
			row.referenceInstanceCount++
			row.referenceBoostedEdges += len(boostedEdges)
			row.optimalEdges += len(optimalEdges)
			row.overlapEdges += overlap
		}

		rows = append(rows, row)
	}

	return rows, nil
}

func readOrCreateIndividualMsas(atspData AtspData) ([][][]float64, error) {
	msas, err := msaSupport.ReadMsas(atspData.msaSupportDirectoryPath)
	if err == nil && len(msas) == len(atspData.matrix) {
		return msas, nil
	}

	if _, createErr := msaSupport.Create(atspData.matrix, atspData.msaSupportDirectoryPath); createErr != nil {
		if err != nil {
			return nil, fmt.Errorf("%s: read individual MSAs: %w; create MSA support: %w", atspData.name, err, createErr)
		}
		return nil, fmt.Errorf("%s: create MSA support: %w", atspData.name, createErr)
	}

	msas, err = msaSupport.ReadMsas(atspData.msaSupportDirectoryPath)
	if err != nil {
		return nil, fmt.Errorf("%s: read individual MSAs: %w", atspData.name, err)
	}

	return msas, nil
}

func writeMsaCountScalingFindings(builder *strings.Builder, rows []msaCountScalingRow) {
	if len(rows) == 0 {
		return
	}

	first := rows[0]
	last := rows[len(rows)-1]
	fmt.Fprintf(builder, "## Findings\n\n")
	fmt.Fprintf(builder, "- **Precision changes from %.2f%% with %s MSA to %.2f%% with %s MSAs.**\n",
		100*first.precision(),
		first.label(),
		100*last.precision(),
		last.label())
	fmt.Fprintf(builder, "- **Recall changes from %.2f%% with %s MSA to %.2f%% with %s MSAs.**\n",
		100*first.recall(),
		first.label(),
		100*last.recall(),
		last.label())
}

func writeMsaCountScalingTable(builder *strings.Builder, rows []msaCountScalingRow) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, row := range rows {
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(row.label()),
			100*row.boostedToTarget(),
			100*row.precision(),
			100*row.recall())
	}

	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

type msaCountScalingRow struct {
	requestedCount         int
	instanceCount          int
	referenceInstanceCount int
	boostedEdges           int
	boostedTargetEdges     int
	referenceBoostedEdges  int
	optimalEdges           int
	overlapEdges           int
}

func (row msaCountScalingRow) label() string {
	if row.requestedCount == 0 {
		return "all"
	}

	return strconv.Itoa(row.requestedCount)
}

func (row msaCountScalingRow) boostedToTarget() float64 {
	return ratio(row.boostedEdges, row.boostedTargetEdges)
}

func (row msaCountScalingRow) precision() float64 {
	return ratio(row.overlapEdges, row.referenceBoostedEdges)
}

func (row msaCountScalingRow) recall() float64 {
	return ratio(row.overlapEdges, row.optimalEdges)
}

func buildPartialMsaSupportEdgeSet(msas [][][]float64, selectedRoots []int) map[models.Edge]struct{} {
	selectedRootSet := make(map[int]struct{}, len(selectedRoots))
	support := make(map[models.Edge]int)
	dimension := 0
	if len(msas) != 0 {
		dimension = len(msas[0])
	}

	for _, root := range selectedRoots {
		selectedRootSet[root] = struct{}{}
		for from := 0; from < len(msas[root]); from++ {
			for to, value := range msas[root][from] {
				if from == to || value <= 0 {
					continue
				}

				support[models.Edge{From: from, To: to}]++
			}
		}
	}

	edges := make(map[models.Edge]struct{})
	for edge, count := range support {
		eligibleRoots := len(selectedRoots)
		if _, ok := selectedRootSet[edge.To]; ok {
			eligibleRoots--
		}
		if eligibleRoots > 0 && count == eligibleRoots && edge.From >= 0 && edge.From < dimension && edge.To >= 0 && edge.To < dimension {
			edges[edge] = struct{}{}
		}
	}

	return edges
}

func selectEvenlySpacedRootIndexes(total, count int) []int {
	if count <= 0 || count >= total {
		roots := make([]int, total)
		for i := range roots {
			roots[i] = i
		}
		return roots
	}
	if count == 1 {
		return []int{0}
	}

	roots := make([]int, count)
	for i := 0; i < count; i++ {
		roots[i] = i * (total - 1) / (count - 1)
	}

	return roots
}

func buildAnalysisTourEdgeSet(tours map[string][]int) map[models.Edge]struct{} {
	edges := make(map[models.Edge]struct{})
	tourIds := make([]string, 0, len(tours))
	for tourId := range tours {
		tourIds = append(tourIds, tourId)
	}
	sort.Strings(tourIds)

	for _, tourId := range tourIds {
		for _, edge := range models.ConvertTourToEdges(tours[tourId]) {
			edges[edge] = struct{}{}
		}
	}

	return edges
}

func countEdgeSetIntersection(left, right map[models.Edge]struct{}) int {
	count := 0
	for edge := range left {
		if _, ok := right[edge]; ok {
			count++
		}
	}

	return count
}

func maxIntValue(left, right int) int {
	if left > right {
		return left
	}

	return right
}

func sortedCycleCoverAnalyses(analyses []cycleCover.InstanceAnalysis) []cycleCover.InstanceAnalysis {
	rows := append([]cycleCover.InstanceAnalysis(nil), analyses...)
	sort.SliceStable(rows, func(i, j int) bool {
		return rows[i].Instance < rows[j].Instance
	})

	return rows
}

func filterAnalysesWithFoundOptimalEdges(analyses []cycleCover.InstanceAnalysis) []cycleCover.InstanceAnalysis {
	rows := make([]cycleCover.InstanceAnalysis, 0, len(analyses))
	for _, analysis := range analyses {
		if analysis.Metrics.UniqueFoundOptimalEdgeCount == 0 {
			continue
		}

		rows = append(rows, analysis)
	}

	return rows
}

func bestPositivePair(left, right float64) (bool, bool) {
	const epsilon = 1e-9
	best := math.Max(left, right)
	if best <= 0 {
		return false, false
	}

	return math.Abs(left-best) < epsilon, math.Abs(right-best) < epsilon
}

func ratio(numerator, denominator int) float64 {
	if denominator == 0 {
		return 0
	}

	return float64(numerator) / float64(denominator)
}

func readFinalResultSummaryMetrics(resultCsvPath string) (map[string]finalResultsSummaryMetric, error) {
	file, err := os.Open(resultCsvPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		return nil, err
	}

	heuristicIndex := indexOf(header, "Heuristic")
	iterationsIndex := indexOf(header, "Iterations")
	averageBestIterationIndex := indexOf(header, "Avg best at iteration")
	averageDeviationIndex := indexOf(header, "Avg best deviation")
	successRateIndex := indexOf(header, "Success rate [%]")
	if heuristicIndex == -1 || iterationsIndex == -1 || averageBestIterationIndex == -1 || averageDeviationIndex == -1 || successRateIndex == -1 {
		return nil, fmt.Errorf("missing required summary columns")
	}

	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	metrics := make(map[string]finalResultsSummaryMetric, len(records))
	for _, record := range records {
		if len(record) != len(header) {
			return nil, fmt.Errorf("invalid record length: got %d want %d", len(record), len(header))
		}

		iterations, err := strconv.Atoi(record[iterationsIndex])
		if err != nil {
			return nil, fmt.Errorf("invalid iterations for %s: %w", record[heuristicIndex], err)
		}
		averageBestIteration, err := strconv.ParseFloat(record[averageBestIterationIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid average best iteration for %s: %w", record[heuristicIndex], err)
		}
		averageDeviation, err := strconv.ParseFloat(record[averageDeviationIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid average deviation for %s: %w", record[heuristicIndex], err)
		}
		successRate, err := strconv.ParseFloat(record[successRateIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid success rate for %s: %w", record[heuristicIndex], err)
		}

		metrics[record[heuristicIndex]] = finalResultsSummaryMetric{
			averageMinDeviation:  averageDeviation,
			successRate:          successRate,
			averageBestIteration: averageBestIteration,
			iterations:           iterations,
		}
	}

	return metrics, nil
}

func indexOf(values []string, value string) int {
	for i, candidate := range values {
		if candidate == value {
			return i
		}
	}

	return -1
}

func heuristicDisplayName(heuristic string) string {
	switch heuristic {
	case heuristicBaseline:
		return "Baseline"
	case heuristicMsaSupport:
		return "MSA support"
	case heuristicCycleCover:
		return "Cycle cover"
	default:
		return heuristic
	}
}

func buildHeuristicModifiers(heuristic string, matrix, msaSupport, cycleCover [][]float64, strength float64) [][]float64 {
	switch heuristic {
	case heuristicBaseline:
		return heuristics.BuildNeutralModifiers(heuristicMatrixDimension(matrix, msaSupport, cycleCover))
	case heuristicMsaSupport:
		return heuristics.BuildMsaSupportModifiers(msaSupport, strength)
	case heuristicCycleCover:
		return heuristics.BuildCycleCoverModifiers(cycleCover, strength)
	default:
		return heuristics.BuildNeutralModifiers(heuristicMatrixDimension(matrix, msaSupport, cycleCover))
	}
}

func heuristicMatrixDimension(matrices ...[][]float64) int {
	for _, matrix := range matrices {
		if len(matrix) != 0 {
			return len(matrix)
		}
	}

	return 0
}

func statisticsCsvRecord(statistic ExperimentsDataStatistics) []string {
	floatFormat := "%.2f"
	return []string{
		fmt.Sprintf(floatFormat, statistic.alpha),
		fmt.Sprintf(floatFormat, statistic.beta),
		fmt.Sprintf(floatFormat, statistic.rho),
		fmt.Sprintf(floatFormat, statistic.heuristicWeight),
		strconv.Itoa(statistic.iterations),
		strconv.Itoa(statistic.minBestAtIteration),
		fmt.Sprintf(floatFormat, statistic.averageBestAtIteration),
		strconv.Itoa(statistic.maxBestAtIteration),
		strconv.Itoa(statistic.minThreeOptImprovementsCount),
		fmt.Sprintf(floatFormat, statistic.averageThreeOptImprovementsCount),
		strconv.Itoa(statistic.maxThreeOptImprovementsCount),
		fmt.Sprintf(floatFormat, statistic.minBestDeviation),
		fmt.Sprintf(floatFormat, statistic.averageBestDeviation),
		fmt.Sprintf(floatFormat, statistic.maxBestDeviation),
		fmt.Sprintf(floatFormat, statistic.successRate),
	}
}

func calculateStatistics(experimentsData []ExperimentsData) []ExperimentsDataStatistics {
	statistics := make([]ExperimentsDataStatistics, len(experimentsData))

	for i, data := range experimentsData {
		minBestAtIteration := math.MaxInt
		averageBestAtIteration := 0.0
		maxBestAtIteration := -math.MaxInt

		minThreeOptImprovementsCount := math.MaxInt
		averageThreeOptImprovementsCount := 0.0
		maxThreeOptImprovementsCount := -math.MaxInt

		minBestDeviation := math.MaxFloat64
		averageBestDeviation := 0.0
		maxBestDeviation := -math.MaxFloat64

		successCounter := 0.0
		minDeviationPerIteration := make([]float64, data.iterations)
		averageDeviationPerIteration := make([]float64, data.iterations)
		maxDeviationPerIteration := make([]float64, data.iterations)

		resultsLen := float64(len(data.results))
		for _, result := range data.results {

			if result.bestAtIteration < minBestAtIteration {
				minBestAtIteration = result.bestAtIteration
			}

			averageBestAtIteration += float64(result.bestAtIteration)

			if result.bestAtIteration > maxBestAtIteration {
				maxBestAtIteration = result.bestAtIteration
			}

			if result.threeOptImprovementsCount < minThreeOptImprovementsCount {
				minThreeOptImprovementsCount = result.threeOptImprovementsCount
			}

			averageThreeOptImprovementsCount += float64(result.threeOptImprovementsCount)

			if result.threeOptImprovementsCount > maxThreeOptImprovementsCount {
				maxThreeOptImprovementsCount = result.threeOptImprovementsCount
			}

			bestDeviation := result.deviationPerIteration[result.bestAtIteration]

			if bestDeviation < minBestDeviation {
				minBestDeviation = bestDeviation
				copy(minDeviationPerIteration, result.deviationPerIteration)
			}

			averageBestDeviation += bestDeviation
			for i, deviation := range result.deviationPerIteration {
				averageDeviationPerIteration[i] += deviation / resultsLen
			}

			if bestDeviation > maxBestDeviation {
				maxBestDeviation = bestDeviation
				copy(maxDeviationPerIteration, result.deviationPerIteration)
			}

			if bestDeviation == 0 {
				successCounter++
			}
		}

		averageBestAtIteration /= resultsLen
		averageBestDeviation /= resultsLen
		successRate := 100.0 * successCounter / resultsLen

		statistics[i] = ExperimentsDataStatistics{
			ExperimentParameters:             data.ExperimentParameters,
			minBestAtIteration:               minBestAtIteration,
			averageBestAtIteration:           averageBestAtIteration,
			maxBestAtIteration:               maxBestAtIteration,
			minThreeOptImprovementsCount:     minThreeOptImprovementsCount,
			averageThreeOptImprovementsCount: averageThreeOptImprovementsCount,
			maxThreeOptImprovementsCount:     maxThreeOptImprovementsCount,
			minBestDeviation:                 minBestDeviation,
			averageBestDeviation:             averageBestDeviation,
			maxBestDeviation:                 maxBestDeviation,
			successRate:                      successRate,
			minDeviationPerIteration:         minDeviationPerIteration,
			averageDeviationPerIteration:     averageDeviationPerIteration,
			maxDeviationPerIteration:         maxDeviationPerIteration,
		}
	}

	sort.SliceStable(statistics, func(i, j int) bool {
		if statistics[i].averageBestDeviation != statistics[j].averageBestDeviation {
			return statistics[i].averageBestDeviation < statistics[j].averageBestDeviation
		}

		if statistics[i].successRate != statistics[j].successRate {
			return statistics[i].successRate > statistics[j].successRate
		}

		return statistics[i].averageBestAtIteration < statistics[j].averageBestAtIteration
	})

	return statistics
}

func runExperiments(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, heuristicModifiers [][]float64, useThreeOpt bool) []ExperimentResult {
	results := make([]ExperimentResult, numberOfRuns)

	aco := aco.NewACO(
		parameters.alpha,
		parameters.beta,
		parameters.rho,
		parameters.iterations,
		knownOptimal,
		matrix,
		heuristicModifiers)
	aco.SetUseThreeOpt(useThreeOpt)

	for i := 0; i < numberOfRuns; i++ {

		aco.Run()

		results[i] = ExperimentResult{
			bestAtIteration:           aco.BestAtIteration,
			threeOptImprovementsCount: aco.ThreeOptImprovementsCount,
			bestTour:                  aco.BestTour,
			deviationPerIteration:     aco.DeviationPerIteration,
		}
	}

	return results
}

func setDimensionDependantParameters(dimension int, parameters *ExperimentParameters) {
	var iterations = 100

	// https://sci-hub.se/10.1109/ICICTA.2010.731
	// "If the number of cities is less than 50, t_max=100; if it is between 50 and 100, t_max=500; and if the problem has more than 100 cities, t_max is set to 5000."
	if dimension < 50 {
		iterations = 100
	}

	if 50 <= dimension && dimension < 100 {
		iterations = 500
	}

	if dimension >= 100 {
		iterations = 5000
	}

	parameters.iterations = iterations
}

func generateParameters() []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0)

	for _, alpha := range utilities.GenerateRange(defaultExperimentAlpha, defaultExperimentAlpha, 0.25) {
		for _, beta := range utilities.GenerateRange(defaultExperimentBeta, defaultExperimentBeta, 1.0) {
			for _, rho := range utilities.GenerateRange(defaultExperimentRho, defaultExperimentRho, 0.1) {
				for _, heuristicWeight := range []float64{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0} {

					parameters = append(parameters,
						ExperimentParameters{
							alpha, beta, rho, heuristicWeight, 0,
						})
				}
			}
		}
	}

	return parameters
}

type finalExperimentConfiguration struct {
	heuristic  string
	parameters []ExperimentParameters
}

func finalExperimentConfigurations() []finalExperimentConfiguration {
	return []finalExperimentConfiguration{
		{
			heuristic: heuristicBaseline,
			parameters: []ExperimentParameters{
				newDefaultExperimentParameters(defaultBaselineHeuristicWeight),
			},
		},
		{
			heuristic: heuristicMsaSupport,
			parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalMsaSupportWeight),
			},
		},
		{
			heuristic: heuristicCycleCover,
			parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalCycleCoverWeight),
			},
		},
	}
}

func selectFinalExperimentConfigurations(finalHeuristic string) ([]finalExperimentConfiguration, error) {
	configurations := finalExperimentConfigurations()
	if finalHeuristic == finalHeuristicAll {
		return configurations, nil
	}

	for _, configuration := range configurations {
		if configuration.heuristic == finalHeuristic {
			return []finalExperimentConfiguration{configuration}, nil
		}
	}

	return nil, fmt.Errorf("unsupported final heuristic %q", finalHeuristic)
}

func newDefaultExperimentParameters(heuristicWeight float64) ExperimentParameters {
	return ExperimentParameters{
		alpha:           defaultExperimentAlpha,
		beta:            defaultExperimentBeta,
		rho:             defaultExperimentRho,
		heuristicWeight: heuristicWeight,
	}
}

func buildMinimumCycleCoverMatrix(matrix [][]float64) ([][]float64, float64, error) {
	edges, cost, err := hungarian.MinimumCycleCover(matrix)
	if err != nil {
		return nil, 0, err
	}

	cycleCover := make([][]float64, len(matrix))
	for i := range cycleCover {
		cycleCover[i] = make([]float64, len(matrix))
	}

	for _, edge := range edges {
		cycleCover[edge.From][edge.To] = 1.0
	}

	return cycleCover, cost, nil
}

var resultsDirectoryName = "results"
var finalResultsDirectoryName = filepath.Join(resultsDirectoryName, "final")
var finalThreeOptResultsDirectoryName = filepath.Join(resultsDirectoryName, "final_3opt")
var resultFileName = "result.csv"

type AtspData struct {
	name         string
	matrix       [][]float64
	knownOptimal float64

	msaSupportDirectoryPath,

	msaSupportHeatmapPlotPath, msaSupportHistogramPlotPath,

	resultFilePath,
	resultPlotFilePrefix,

	optimalUniqueToursCsvPath,
	toursHeatmapPlotPath,
	toursHistogramPlotPath,
	msaSupportToursOverlapHeatmapPlotPath,
	msaSupportSolutionAnalysisCsvPath,
	msaSupportSolutionThresholdsCsvPath,
	cycleCoverEdgesCsvPath,
	cycleCoverAnalysisCsvPath,
	cycleCoverThresholdsCsvPath,
	cycleCoverMsaSupportOverlapCsvPath string
}

func makeAtspData(name string, matrix [][]float64, knownOptimal float64) AtspData {
	return makeAtspDataInResultsDirectory(name, matrix, knownOptimal, resultsDirectoryName)
}

func makeAtspDataInResultsDirectory(name string, matrix [][]float64, knownOptimal float64, resultsRootPath string) AtspData {
	name = strings.TrimSuffix(name, ".atsp")
	resultsDirectoryPath := filepath.Join(resultsRootPath, name)
	msaSupportDirectoryPath := filepath.Join(resultsDirectoryPath, "msa_support")
	plotsDirectoryPath := filepath.Join(resultsDirectoryPath, "plots")

	msaSupportHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "msa_support_heatmap.png")
	msaSupportHistogramPlotPath := filepath.Join(plotsDirectoryPath, "msa_support_histogram.png")

	resultFilePath := filepath.Join(resultsDirectoryPath, resultFileName)
	resultPlotFilePrefix := filepath.Join(plotsDirectoryPath, "best_result")

	optimalUniqueToursCsvPath := filepath.Join(resultsDirectoryPath, "solutions.csv")
	toursHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "tours_heatmap.png")
	toursHistogramPlotPath := filepath.Join(plotsDirectoryPath, "tours_histogram.png")
	msaSupportToursOverlapHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "msa_support_tours_overlap_heatmap.png")
	msaSupportSolutionAnalysisCsvPath := filepath.Join(resultsDirectoryPath, "msa_support_solution_analysis.csv")
	msaSupportSolutionThresholdsCsvPath := filepath.Join(resultsDirectoryPath, "msa_support_solution_thresholds.csv")
	cycleCoverEdgesCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_edges.csv")
	cycleCoverAnalysisCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_analysis.csv")
	cycleCoverThresholdsCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_thresholds.csv")
	cycleCoverMsaSupportOverlapCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_msa_support_overlap.csv")

	return AtspData{
		name,
		matrix,
		knownOptimal,

		msaSupportDirectoryPath,

		msaSupportHeatmapPlotPath, msaSupportHistogramPlotPath,

		resultFilePath,
		resultPlotFilePrefix,

		optimalUniqueToursCsvPath,
		toursHeatmapPlotPath,
		toursHistogramPlotPath,
		msaSupportToursOverlapHeatmapPlotPath,
		msaSupportSolutionAnalysisCsvPath,
		msaSupportSolutionThresholdsCsvPath,
		cycleCoverEdgesCsvPath,
		cycleCoverAnalysisCsvPath,
		cycleCoverThresholdsCsvPath,
		cycleCoverMsaSupportOverlapCsvPath,
	}
}

func withExperimentOutputRoot(atspData AtspData, resultsRootPath string) AtspData {
	output := makeAtspDataInResultsDirectory(atspData.name, atspData.matrix, atspData.knownOptimal, resultsRootPath)
	output.msaSupportDirectoryPath = atspData.msaSupportDirectoryPath
	output.optimalUniqueToursCsvPath = atspData.optimalUniqueToursCsvPath
	return output
}

func selectAtspFiles(atspFilePaths []string, instanceSet string) ([]string, error) {
	switch instanceSet {
	case instanceSetSmoke:
		return selectConfiguredAtspFiles(atspFilePaths, smokeInstanceFiles)
	case instanceSetTiny:
		return selectConfiguredAtspFiles(atspFilePaths, tinyInstanceFiles)
	case instanceSetBalanced:
		return selectConfiguredAtspFiles(atspFilePaths, balancedInstanceFiles)
	case instanceSetLarge:
		return selectConfiguredAtspFiles(atspFilePaths, largeInstanceFiles)
	case instanceSetAllKnown:
		selected := make([]string, 0, len(atspFilePaths))
		for _, atspFilePath := range atspFilePaths {
			if parsing.HasKnownOptimalSolution(filepath.Base(atspFilePath)) {
				selected = append(selected, atspFilePath)
			}
		}

		sort.Strings(selected)
		if len(selected) == 0 {
			return nil, fmt.Errorf("no ATSP files with known optima were found")
		}

		return selected, nil
	default:
		return nil, fmt.Errorf("unsupported -instances value %q; use %q, %q, %q, %q, or %q", instanceSet, instanceSetSmoke, instanceSetTiny, instanceSetBalanced, instanceSetLarge, instanceSetAllKnown)
	}
}

func selectConfiguredAtspFiles(atspFilePaths, configuredFiles []string) ([]string, error) {
	pathByFileName := make(map[string]string, len(atspFilePaths))
	for _, atspFilePath := range atspFilePaths {
		pathByFileName[filepath.Base(atspFilePath)] = atspFilePath
	}

	selected := make([]string, 0, len(configuredFiles))
	for _, fileName := range configuredFiles {
		atspFilePath, ok := pathByFileName[fileName]
		if !ok {
			return nil, fmt.Errorf("configured ATSP instance %q was not found", fileName)
		}

		selected = append(selected, atspFilePath)
	}

	return selected, nil
}

func isValidRunMode(mode string) bool {
	return mode == runModeExperiment || mode == runModeAnalyze || mode == runModeAll || mode == runModeFinal || mode == runModeFinal3Opt
}

func isValidHeuristic(heuristic string) bool {
	return heuristic == heuristicBaseline ||
		heuristic == heuristicMsaSupport ||
		heuristic == heuristicCycleCover
}

func heuristicUsesCycleCover(heuristic string) bool {
	return heuristic == heuristicCycleCover
}

func heuristicUsesMsaSupport(heuristic string) bool {
	return heuristic == heuristicMsaSupport
}

func finalConfigurationsUseMsaSupport(configurations []finalExperimentConfiguration) bool {
	for _, configuration := range configurations {
		if heuristicUsesMsaSupport(configuration.heuristic) {
			return true
		}
	}

	return false
}

func heuristicFileSuffix(heuristic string) string {
	switch heuristic {
	case heuristicBaseline:
		return "_baseline"
	case heuristicMsaSupport:
		return ""
	case heuristicCycleCover:
		return "_cycle_cover"
	default:
		return "_" + strings.ReplaceAll(heuristic, "-", "_")
	}
}

func resultFilePathForHeuristic(atspData AtspData, heuristic string) string {
	suffix := heuristicFileSuffix(heuristic)
	if suffix == "" {
		return atspData.resultFilePath
	}

	return strings.TrimSuffix(atspData.resultFilePath, ".csv") + suffix + ".csv"
}

func resultPlotFilePrefixForHeuristic(atspData AtspData, heuristic string) string {
	return atspData.resultPlotFilePrefix + heuristicFileSuffix(heuristic)
}

func bestParametersReportPathForHeuristic(heuristic string) string {
	return "best_parameters_report" + heuristicFileSuffix(heuristic) + ".md"
}

func resultFilePathsForHeuristic(atspsData []AtspData, heuristic string) []string {
	paths := make([]string, 0, len(atspsData))
	for _, atspData := range atspsData {
		paths = append(paths, resultFilePathForHeuristic(atspData, heuristic))
	}
	return paths
}

func shouldRunExperiments(mode string) bool {
	return mode == runModeExperiment || mode == runModeAll
}

func shouldRunAnalysis(mode string) bool {
	return mode == runModeAnalyze || mode == runModeAll
}

func shouldRunFinalExperiments(mode string) bool {
	return mode == runModeFinal || mode == runModeFinal3Opt
}

func finalExperimentOutputRoot(mode string) string {
	if mode == runModeFinal3Opt {
		return finalThreeOptResultsDirectoryName
	}

	return finalResultsDirectoryName
}

func finalExperimentUsesThreeOpt(mode string) bool {
	return mode == runModeFinal3Opt
}

func main() {
	instances := flag.String("instances", instanceSetSmoke, "ATSP instance set to run: smoke, tiny, balanced, large, or all-known")
	mode := flag.String("mode", runModeExperiment, "Run mode: experiment, analyze, all, final, or final+3opt")
	heuristic := flag.String("heuristic", heuristicMsaSupport, "ACO heuristic modifier to use in experiment mode: baseline, msa-support, or cycle-cover")
	finalHeuristic := flag.String("final-heuristic", finalHeuristicAll, "Final-mode heuristic to run: all, baseline, msa-support, or cycle-cover")
	flag.Parse()

	if !isValidRunMode(*mode) {
		fmt.Printf("Unsupported -mode value %q; use %q, %q, %q, %q, or %q\n", *mode, runModeExperiment, runModeAnalyze, runModeAll, runModeFinal, runModeFinal3Opt)
		return
	}

	if !isValidHeuristic(*heuristic) {
		fmt.Printf("Unsupported -heuristic value %q; use %q, %q, or %q\n", *heuristic, heuristicBaseline, heuristicMsaSupport, heuristicCycleCover)
		return
	}

	selectedInstances := *instances
	var finalConfigurations []finalExperimentConfiguration
	if shouldRunFinalExperiments(*mode) {
		selectedInstances = instanceSetBalanced
		configurations, err := selectFinalExperimentConfigurations(*finalHeuristic)
		if err != nil {
			fmt.Println(err)
			return
		}
		finalConfigurations = configurations
	}

	atspsData, err := loadSelectedAtspData(selectedInstances)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Printf("Selected %d ATSP instance(s) with -instances=%s\n", len(atspsData), selectedInstances)

	if shouldRunExperiments(*mode) {
		stopProfiling, err := startCPUProfile()
		if err != nil {
			fmt.Println(err)
			return
		}

		err = runExperimentMode(atspsData, *heuristic)
		stopProfiling()
		if err != nil {
			fmt.Println(err)
			return
		}
	}

	if shouldRunFinalExperiments(*mode) {
		stopProfiling, err := startCPUProfile()
		if err != nil {
			fmt.Println(err)
			return
		}

		err = runFinalExperimentMode(atspsData, finalExperimentOutputRoot(*mode), finalExperimentUsesThreeOpt(*mode), finalConfigurations)
		stopProfiling()
		if err != nil {
			fmt.Println(err)
			return
		}
	}

	if shouldRunAnalysis(*mode) {
		if err := runAnalysisMode(atspsData); err != nil {
			fmt.Println(err)
			return
		}
	}
}

func loadSelectedAtspData(instances string) ([]AtspData, error) {
	tsplibDir := "tsplib_files"
	atspFilesPaths, err := filepath.Glob(filepath.Join(tsplibDir, "*.atsp"))
	if err != nil {
		return nil, err
	}

	atspFilesPaths, err = selectAtspFiles(atspFilesPaths, instances)
	if err != nil {
		return nil, err
	}

	atspsData := make([]AtspData, len(atspFilesPaths))
	for i, atspFilePath := range atspFilesPaths {
		name, matrix, knownOptimal, err := parsing.ParseTSPLIBFile(atspFilePath)
		if err != nil {
			return nil, err
		}
		if knownOptimal == 0 {
			return nil, fmt.Errorf("missing known optimal solution for %s", atspFilePath)
		}

		atspsData[i] = makeAtspData(name, matrix, knownOptimal)
	}

	return atspsData, nil
}

func startCPUProfile() (func(), error) {
	cf, err := os.Create("cpu.prof")
	if err != nil {
		return nil, err
	}

	if err := pprof.StartCPUProfile(cf); err != nil {
		cf.Close()
		return nil, err
	}

	return func() {
		pprof.StopCPUProfile()
		cf.Close()
	}, nil
}

func runExperimentMode(atspsData []AtspData, heuristic string) error {
	if err := ensureMsaSupportArtifacts(atspsData); err != nil {
		return err
	}

	return runExperimentSet(atspsData, resultsDirectoryName, heuristic, generateParameters(), defaultExperimentRunCount)
}

func runFinalExperimentMode(atspsData []AtspData, resultsRootPath string, useThreeOpt bool, configurations []finalExperimentConfiguration) error {
	needsMsaSupport := finalConfigurationsUseMsaSupport(configurations)
	if needsMsaSupport {
		if err := ensureMsaSupportCache(atspsData); err != nil {
			return err
		}
	}

	if err := removeLegacyFinalReports(resultsRootPath); err != nil {
		return err
	}

	for _, atspData := range atspsData {
		finalAtspData := withExperimentOutputRoot(atspData, resultsRootPath)
		if err := runFinalExperimentForInstance(finalAtspData, useThreeOpt, configurations, needsMsaSupport); err != nil {
			return err
		}
	}

	return nil
}

func runFinalExperimentForInstance(atspData AtspData, useThreeOpt bool, configurations []finalExperimentConfiguration, needsMsaSupport bool) error {
	matrix := atspData.matrix
	knownOptimal := atspData.knownOptimal
	dimension := len(matrix)
	instanceStart := time.Now()
	finalRunName := "final"
	if useThreeOpt {
		finalRunName = "final+3opt"
	}

	if len(configurations) == 0 {
		return fmt.Errorf("no final experiment configurations selected")
	}

	if err := removeLegacyFinalResultFiles(atspData); err != nil {
		return err
	}

	var heuristicMatrix [][]float64
	if needsMsaSupport {
		var err error
		heuristicMatrix, err = readMsaSupportMatrixForHeuristic(atspData, heuristicMsaSupport)
		if err != nil {
			return err
		}
	}

	fmt.Printf("Starting %s %s (dimension=%d, heuristics=%d, runs/heuristic=%d)\n",
		finalRunName,
		atspData.name,
		dimension,
		len(configurations),
		finalNumberOfExperiments)

	var cycleCover [][]float64
	var cycleCoverErr error
	cycleCoverReady := false
	finalStatistics := make([]HeuristicExperimentStatistics, 0, len(configurations))

	for _, config := range configurations {
		if heuristicUsesCycleCover(config.heuristic) && !cycleCoverReady {
			var cycleCoverCost float64
			cycleCover, cycleCoverCost, cycleCoverErr = buildMinimumCycleCoverMatrix(matrix)
			if cycleCoverErr != nil {
				return cycleCoverErr
			}
			cycleCoverReady = true
			fmt.Printf("\tMinimum cycle cover cost=%.2f gap=%.2f%%\n",
				cycleCoverCost,
				100*(knownOptimal-cycleCoverCost)/knownOptimal)
		}

		experimentData := make([]ExperimentsData, 0, len(config.parameters))
		for _, parameters := range config.parameters {
			setDimensionDependantParameters(dimension, &parameters)
			parameterStart := time.Now()
			heuristicModifiers := buildHeuristicModifiers(config.heuristic, matrix, heuristicMatrix, cycleCover, parameters.heuristicWeight)
			results := runExperiments(finalNumberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, useThreeOpt)
			data := ExperimentsData{parameters, results}
			experimentData = append(experimentData, data)

			parameterStatistics := calculateStatistics([]ExperimentsData{data})
			if len(parameterStatistics) != 0 {
				statistic := parameterStatistics[0]
				fmt.Printf("\t%s heuristicWeight=%.2f iterations=%d runs=%d elapsed=%s min deviation=%.2f avg deviation=%.2f\n",
					config.heuristic,
					parameters.heuristicWeight,
					parameters.iterations,
					finalNumberOfExperiments,
					time.Since(parameterStart).Round(time.Millisecond),
					statistic.minBestDeviation,
					statistic.averageBestDeviation)
			}
		}

		statistics := calculateStatistics(experimentData)
		if len(statistics) == 0 {
			continue
		}

		finalStatistics = append(finalStatistics, HeuristicExperimentStatistics{
			heuristic:  config.heuristic,
			statistics: statistics[0],
		})

		if err := removeExperimentPlotsForHeuristic(atspData, config.heuristic); err != nil {
			return err
		}
		saveExperimentPlots(statistics, "MMAS deviation per iteration", resultPlotFilePrefixForHeuristic(atspData, config.heuristic))
	}

	if err := saveFinalHeuristicStatistics(atspData.resultFilePath, finalStatistics, configurations); err != nil {
		return err
	}

	fmt.Printf("Finished %s %s in %s\n", finalRunName, atspData.name, time.Since(instanceStart).Round(time.Millisecond))
	return nil
}

func removeLegacyFinalReports(resultsRootPath string) error {
	legacyReports := []string{
		filepath.Join(resultsRootPath, bestParametersReportPathForHeuristic(heuristicBaseline)),
		filepath.Join(resultsRootPath, bestParametersReportPathForHeuristic(heuristicMsaSupport)),
		filepath.Join(resultsRootPath, bestParametersReportPathForHeuristic(heuristicCycleCover)),
	}

	for _, path := range legacyReports {
		if err := removeFileIfExists(path); err != nil {
			return err
		}
	}

	return nil
}

func removeLegacyFinalResultFiles(atspData AtspData) error {
	legacyFiles := []string{
		resultFilePathForHeuristic(atspData, heuristicBaseline),
		resultFilePathForHeuristic(atspData, heuristicCycleCover),
		filepath.Join(filepath.Dir(atspData.resultFilePath), "solutions.csv"),
	}

	for _, path := range legacyFiles {
		if err := removeFileIfExists(path); err != nil {
			return err
		}
	}

	return nil
}

func removeFileIfExists(path string) error {
	if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
		return err
	}

	return nil
}

func runExperimentSet(atspsData []AtspData, resultsRootPath, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments int) error {
	for _, atspData := range atspsData {
		matrix := atspData.matrix
		knownOptimal := atspData.knownOptimal
		dimension := len(matrix)
		instanceStart := time.Now()

		heuristicMatrix, err := readMsaSupportMatrixForHeuristic(atspData, heuristic)
		if err != nil {
			return err
		}

		fmt.Printf("Starting %s (dimension=%d, heuristic=%s, parameters=%d, runs/parameter=%d)\n",
			atspData.name,
			dimension,
			heuristic,
			len(experimentParameters),
			numberOfExperiments)

		var cycleCover [][]float64
		if heuristicUsesCycleCover(heuristic) {
			var cycleCoverCost float64
			cycleCover, cycleCoverCost, err = buildMinimumCycleCoverMatrix(matrix)
			if err != nil {
				return err
			}
			fmt.Printf("\tMinimum cycle cover cost=%.2f gap=%.2f%%\n",
				cycleCoverCost,
				100*(knownOptimal-cycleCoverCost)/knownOptimal)
		}

		experimentData := make([]ExperimentsData, 0)

		for _, parameters := range experimentParameters {
			setDimensionDependantParameters(dimension, &parameters)
			parameterStart := time.Now()
			heuristicModifiers := buildHeuristicModifiers(heuristic, matrix, heuristicMatrix, cycleCover, parameters.heuristicWeight)
			results := runExperiments(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, false)
			data := ExperimentsData{parameters, results}

			experimentData = append(experimentData, data)

			parameterStatistics := calculateStatistics([]ExperimentsData{data})
			if len(parameterStatistics) != 0 {
				statistic := parameterStatistics[0]
				fmt.Printf("\theuristicWeight=%.2f iterations=%d runs=%d elapsed=%s min deviation=%.2f avg deviation=%.2f\n",
					parameters.heuristicWeight,
					parameters.iterations,
					numberOfExperiments,
					time.Since(parameterStart).Round(time.Millisecond),
					statistic.minBestDeviation,
					statistic.averageBestDeviation)
			}
		}

		statistics := calculateStatistics(experimentData)
		if len(statistics) != 0 {
			saveStatistics(resultFilePathForHeuristic(atspData, heuristic), statistics)
			if err := removeExperimentPlotsForHeuristic(atspData, heuristic); err != nil {
				return err
			}
			saveExperimentPlots(statistics, "MMAS deviation per iteration", resultPlotFilePrefixForHeuristic(atspData, heuristic))
		}

		uniqueOptimalTours, err := msaSupportTours.ReadOptimalTours(atspData.optimalUniqueToursCsvPath)
		if err != nil {
			return err
		}

		for _, data := range experimentData {
			for _, result := range data.results {
				if result.deviationPerIteration[result.bestAtIteration] == 0.0 {
					msaSupportTours.AddUniqueTour(uniqueOptimalTours, result.bestTour)
				}
			}
		}

		if err := msaSupportTours.SaveOptimalToursStatistics(atspData.optimalUniqueToursCsvPath, atspData.msaSupportDirectoryPath, uniqueOptimalTours); err != nil {
			return err
		}

		fmt.Printf("Finished %s in %s\n", atspData.name, time.Since(instanceStart).Round(time.Millisecond))
	}

	bestStatistics := getBestStatisticsFromFiles(resultFilePathsForHeuristic(atspsData, heuristic))
	saveBestParametersInfo(resultsRootPath, bestParametersReportPathForHeuristic(heuristic), bestStatistics)

	return nil
}

func readMsaSupportMatrixForHeuristic(atspData AtspData, heuristic string) ([][]float64, error) {
	return msaSupport.Read(atspData.msaSupportDirectoryPath)
}

func ensureMsaSupportArtifacts(atspsData []AtspData) error {
	for _, atspData := range atspsData {
		name := atspData.name
		matrix := atspData.matrix
		msaSupportDirectoryPath := atspData.msaSupportDirectoryPath

		msaSupportMatrix, err := msaSupport.Read(atspData.msaSupportDirectoryPath)

		if err != nil {
			start := time.Now()
			msaSupportMatrix, err = msaSupport.Create(matrix, msaSupportDirectoryPath)
			elapsed := time.Since(start)

			fmt.Printf("\tCreating %s took: %d ms\n", msaSupportDirectoryPath, elapsed.Milliseconds())

			if err != nil {
				return fmt.Errorf("error saving MSA support: %w", err)
			}
		}

		msaSupportHeatmapPlotTitle := name + " MSA support heatmap"

		err = utilities.SaveHeatmapFromMatrix(msaSupportMatrix, msaSupportHeatmapPlotTitle, atspData.msaSupportHeatmapPlotPath)
		if err != nil {
			return err
		}

		dataForHistogram := filterZeroes(flattenMatrix(msaSupportMatrix))
		msaSupportHistogramPlotTitle := name + " MSA support histogram"

		dimension := len(matrix)
		err = utilities.SaveHistogramFromData(dataForHistogram, dimension-1, msaSupportHistogramPlotTitle, atspData.msaSupportHistogramPlotPath)
		if err != nil {
			return err
		}
	}

	return nil
}

func ensureMsaSupportCache(atspsData []AtspData) error {
	for _, atspData := range atspsData {
		if _, err := msaSupport.Read(atspData.msaSupportDirectoryPath); err == nil {
			continue
		}

		start := time.Now()
		if _, err := msaSupport.Create(atspData.matrix, atspData.msaSupportDirectoryPath); err != nil {
			return fmt.Errorf("error saving MSA support: %w", err)
		}

		fmt.Printf("\tCreating %s took: %d ms\n", atspData.msaSupportDirectoryPath, time.Since(start).Milliseconds())
	}

	return nil
}

func runAnalysisMode(atspsData []AtspData) error {
	if err := ensureMsaSupportArtifacts(atspsData); err != nil {
		return err
	}

	msaSupportTourConfigs := make([]msaSupportTours.InstanceConfig, 0, len(atspsData))
	cycleCoverConfigs := make([]cycleCover.InstanceConfig, 0, len(atspsData))
	for _, atspData := range atspsData {
		msaSupportTourConfigs = append(msaSupportTourConfigs, msaSupportTours.InstanceConfig{
			Name:                              atspData.name,
			Dimension:                         len(atspData.matrix),
			MsaSupportDirectoryPath:           atspData.msaSupportDirectoryPath,
			OptimalToursCsvPath:               atspData.optimalUniqueToursCsvPath,
			ToursHeatmapPath:                  atspData.toursHeatmapPlotPath,
			ToursHistogramPath:                atspData.toursHistogramPlotPath,
			MsaSupportToursOverlapHeatmapPath: atspData.msaSupportToursOverlapHeatmapPlotPath,
			AnalysisCsvPath:                   atspData.msaSupportSolutionAnalysisCsvPath,
			ThresholdsCsvPath:                 atspData.msaSupportSolutionThresholdsCsvPath,
		})

		cycleCoverConfigs = append(cycleCoverConfigs, cycleCover.InstanceConfig{
			Name:                        atspData.name,
			Dimension:                   len(atspData.matrix),
			Matrix:                      atspData.matrix,
			KnownOptimal:                atspData.knownOptimal,
			MsaSupportDirectoryPath:     atspData.msaSupportDirectoryPath,
			OptimalToursCsvPath:         atspData.optimalUniqueToursCsvPath,
			CycleCoverEdgesCsvPath:      atspData.cycleCoverEdgesCsvPath,
			AnalysisCsvPath:             atspData.cycleCoverAnalysisCsvPath,
			ThresholdsCsvPath:           atspData.cycleCoverThresholdsCsvPath,
			CycleCoverOverlapMatrixPath: atspData.cycleCoverMsaSupportOverlapCsvPath,
		})
	}

	msaSupportSolutionSummaryPath := filepath.Join(resultsDirectoryName, "msa_support_solution_analysis_summary.csv")
	msaSupportSolutionReportPath := filepath.Join(resultsDirectoryName, "msa_support_solution_analysis_report.md")

	_, err := msaSupportTours.AnalyzeInstances(msaSupportTours.Config{
		Instances:      msaSupportTourConfigs,
		SummaryCsvPath: msaSupportSolutionSummaryPath,
		ReportPath:     msaSupportSolutionReportPath,
		HighThreshold:  0.8,
		Thresholds:     msaSupportTours.DefaultThresholds(),
	})
	if err != nil {
		return err
	}

	cycleCoverSummaryPath := filepath.Join(resultsDirectoryName, "cycle_cover_analysis_summary.csv")
	cycleCoverReportPath := filepath.Join(resultsDirectoryName, "cycle_cover_analysis_report.md")

	cycleCoverAnalyses, err := cycleCover.AnalyzeInstances(cycleCover.Config{
		Instances:      cycleCoverConfigs,
		SummaryCsvPath: cycleCoverSummaryPath,
		ReportPath:     cycleCoverReportPath,
		HighThreshold:  1.0,
		Thresholds:     msaSupportTours.DefaultThresholds(),
	})
	if err != nil {
		return err
	}

	structuralSimilarityReportPath := filepath.Join(finalResultsDirectoryName, "structural_similarity.md")
	if err := saveStructuralSimilarityReport(structuralSimilarityReportPath, cycleCoverAnalyses); err != nil {
		return err
	}

	heuristicOverlapReportPath := filepath.Join(finalResultsDirectoryName, "msa_cycle_cover_overlap.md")
	if err := saveMsaSupportCycleCoverOverlapReport(heuristicOverlapReportPath, cycleCoverAnalyses); err != nil {
		return err
	}

	msaCountScalingReportPath := filepath.Join(finalResultsDirectoryName, "msa_count_scaling.md")
	if err := saveMsaCountScalingReport(msaCountScalingReportPath, atspsData); err != nil {
		return err
	}

	heuristicBoostSummaryPath := filepath.Join(resultsDirectoryName, "heuristic_boosted_edges_summary.csv")
	heuristicBoostRows, err := buildHeuristicBoostSummary(atspsData)
	if err != nil {
		return err
	}
	if err := saveHeuristicBoostSummary(heuristicBoostSummaryPath, heuristicBoostRows); err != nil {
		return err
	}

	finalResultsSummaryPath, finalRows, finalSummarySaved, err := runFinalResultsAnalysis(atspsData, cycleCoverAnalyses, finalResultsDirectoryName)
	if err != nil {
		return err
	}

	finalThreeOptResultsSummaryPath, finalThreeOptRows, finalThreeOptSummarySaved, err := runFinalResultsAnalysis(atspsData, cycleCoverAnalyses, finalThreeOptResultsDirectoryName)
	if err != nil {
		return err
	}

	fmt.Printf("MSA support/solution analysis summary saved to %s\n", msaSupportSolutionSummaryPath)
	fmt.Printf("MSA support/solution analysis report saved to %s\n", msaSupportSolutionReportPath)
	fmt.Printf("Cycle-cover analysis summary saved to %s\n", cycleCoverSummaryPath)
	fmt.Printf("Cycle-cover analysis report saved to %s\n", cycleCoverReportPath)
	fmt.Printf("Structural similarity report saved to %s\n", structuralSimilarityReportPath)
	fmt.Printf("MSA support/cycle-cover overlap report saved to %s\n", heuristicOverlapReportPath)
	fmt.Printf("MSA count scaling report saved to %s\n", msaCountScalingReportPath)
	fmt.Printf("Heuristic boosted-edge summary saved to %s\n", heuristicBoostSummaryPath)
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

func runFinalResultsAnalysis(atspsData []AtspData, cycleCoverAnalyses []cycleCover.InstanceAnalysis, resultsRootPath string) (string, []finalResultsSummaryRow, bool, error) {
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
	if len(cycleCoverAnalyses) != 0 {
		if err := saveStructuralPerformanceLinkReport(filepath.Join(resultsRootPath, "structural_performance_link.md"), finalRows, cycleCoverAnalyses); err != nil {
			return "", nil, false, err
		}
	}
	if err := removeFileIfExists(filepath.Join(resultsRootPath, "summary.csv")); err != nil {
		return "", nil, false, err
	}

	return finalResultsSummaryPath, finalRows, true, nil
}

type heuristicBoostSummaryRow struct {
	instance            string
	heuristic           string
	dimension           int
	boostableEdges      int
	tourEdgeTarget      int
	boostedEdges        int
	boostedEdgeDensity  float64
	boostedToTourTarget float64
}

func buildHeuristicBoostSummary(atspsData []AtspData) ([]heuristicBoostSummaryRow, error) {
	const referenceStrength = 1.0
	heuristicNames := []string{
		heuristicMsaSupport,
		heuristicCycleCover,
	}

	instances := append([]AtspData(nil), atspsData...)
	sort.SliceStable(instances, func(i, j int) bool {
		return instances[i].name < instances[j].name
	})

	rows := make([]heuristicBoostSummaryRow, 0, len(instances)*len(heuristicNames))
	for _, atspData := range instances {
		var cycleCover [][]float64
		var cycleCoverErr error
		cycleCoverReady := false

		for _, heuristic := range heuristicNames {
			heuristicMatrix, err := readMsaSupportMatrixForHeuristic(atspData, heuristic)
			if err != nil {
				return nil, fmt.Errorf("%s/%s: failed to read heuristic matrix: %w", atspData.name, heuristic, err)
			}

			if heuristicUsesCycleCover(heuristic) && !cycleCoverReady {
				cycleCover, _, cycleCoverErr = buildMinimumCycleCoverMatrix(atspData.matrix)
				if cycleCoverErr != nil {
					return nil, fmt.Errorf("%s: failed to build minimum cycle cover: %w", atspData.name, cycleCoverErr)
				}
				cycleCoverReady = true
			}

			modifiers := buildHeuristicModifiers(heuristic, atspData.matrix, heuristicMatrix, cycleCover, referenceStrength)
			rows = append(rows, analyzeHeuristicBoosts(atspData.name, heuristic, modifiers))
		}
	}

	return rows, nil
}

func analyzeHeuristicBoosts(instance, heuristic string, modifiers [][]float64) heuristicBoostSummaryRow {
	dimension := len(modifiers)
	totalDirectedEdges := dimension * (dimension - 1)
	tourEdgeTarget := 0
	if dimension > 1 {
		tourEdgeTarget = dimension - 1
	}
	boostedEdges := 0

	for i := 0; i < dimension; i++ {
		for j := 0; j < len(modifiers[i]); j++ {
			if i == j || modifiers[i][j] <= 1.0 {
				continue
			}

			boostedEdges++
		}
	}

	boostedEdgeDensity := 0.0
	if totalDirectedEdges > 0 {
		boostedEdgeDensity = float64(boostedEdges) / float64(totalDirectedEdges)
	}

	boostedToTourTarget := 0.0
	if tourEdgeTarget > 0 {
		boostedToTourTarget = float64(boostedEdges) / float64(tourEdgeTarget)
	}

	return heuristicBoostSummaryRow{
		instance:            instance,
		heuristic:           heuristic,
		dimension:           dimension,
		boostableEdges:      totalDirectedEdges,
		tourEdgeTarget:      tourEdgeTarget,
		boostedEdges:        boostedEdges,
		boostedEdgeDensity:  boostedEdgeDensity,
		boostedToTourTarget: boostedToTourTarget,
	}
}

func saveHeuristicBoostSummary(path string, rows []heuristicBoostSummaryRow) error {
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		return err
	}

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	header := []string{
		"Instance",
		"Heuristic",
		"Boostable edges",
		"Tour edge target",
		"Boosted edges",
		"Boosted edge density [%]",
		"Boosted edges vs tour target [%]",
	}
	if err := writer.Write(header); err != nil {
		return err
	}

	for _, row := range rows {
		record := []string{
			row.instance,
			row.heuristic,
			strconv.Itoa(row.boostableEdges),
			strconv.Itoa(row.tourEdgeTarget),
			strconv.Itoa(row.boostedEdges),
			fmt.Sprintf("%.2f", 100.0*row.boostedEdgeDensity),
			fmt.Sprintf("%.2f", 100.0*row.boostedToTourTarget),
		}
		if err := writer.Write(record); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func saveExperimentPlots(statistics []ExperimentsDataStatistics, plotTitle, plotPathPrefix string) {
	bestStatistic := statistics[0]

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

		titleSuffix := fmt.Sprintf(" (alpha=%.2f, beta=%.2f, rho=%.2f, heuristicWeight=%.2f)",
			statistic.alpha, statistic.beta, statistic.rho, statistic.heuristicWeight)

		heuristicWeightPlotSuffix := "_heuristicWeight=" + strconv.Itoa(int(100*statistic.heuristicWeight)) + "%"
		plotPath := plotPathPrefix + heuristicWeightPlotSuffix + ".png"

		utilities.SaveLinePlotFromData(lines, plotTitle+titleSuffix, plotPath)
	}
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

func getBestStatisticsFromFiles(resultsFilePaths []string) []ExperimentsDataStatistics {
	topNumber := 3
	bestStatistics := make([]ExperimentsDataStatistics, 0)

	for _, path := range resultsFilePaths {
		statistics, err := readStatistics(path)
		if err != nil {
			fmt.Println(err)
			continue
		}

		statisticsCount := len(statistics)
		if statisticsCount < topNumber {
			topNumber = len(statistics)
		}

		top := statistics[:topNumber]
		bestStatistics = append(bestStatistics, top...)
	}

	return bestStatistics
}

func filterZeroes(data []float64) []float64 {
	var result []float64
	for _, value := range data {
		if value != 0 {
			result = append(result, value)
		}
	}
	return result
}

func flattenMatrix(matrix [][]float64) []float64 {
	var result []float64
	for _, row := range matrix {
		result = append(result, row...)
	}
	return result
}
