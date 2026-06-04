package main

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"atsp_aco_msa/modules/algorithms/heuristics"
	"atsp_aco_msa/modules/algorithms/hungarian"
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	"atsp_aco_msa/modules/analysis/cycleCover"
	"atsp_aco_msa/modules/analysis/msaHeuristicTours"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/parsing"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
	"errors"
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
	alpha, beta, rho, heuristicWeight, msaPatchBias float64
	randomSeed                                      int64
	iterations                                      int
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
	instanceSetTuning     = "tuning"
	instanceSetEvaluation = "evaluation"
	instanceSetAllKnown   = "all-known"
)

const (
	runModeExperiment = "experiment"
	runModeAnalyze    = "analyze"
	runModeAll        = "all"
	runModeFinal      = "final"
	runModeFinal3Opt  = "final+3opt"
)

const (
	analysisScopeAll          = "all"
	analysisScopeGksDeviation = "gks-deviation"
)

const (
	finalNumberOfExperiments               = 50
	finalMsaHeuristicWeight                = 0.9
	finalCycleCoverWeight                  = 0.8
	finalCycleCoverMsaPatchingWeight       = 0.7
	finalCycleCoverMsaPatchingMsaPatchBias = 0.5
	defaultExperimentAlpha                 = 1.0
	defaultExperimentBeta                  = 2.0
	defaultExperimentRho                   = 0.8
	defaultExperimentRunCount              = 30
	defaultBaselineHeuristicWeight         = 0.0
)

const (
	heuristicBaseline              = "baseline"
	heuristicMsaHeuristic          = "msa-heuristic"
	heuristicRandomSparse          = "random-sparse"
	heuristicDistanceRankedSparse  = "distance-ranked-sparse"
	heuristicCycleCover            = "cycle-cover"
	heuristicCycleCoverMsaPatching = "cycle-cover-msa-patching"
)

const finalHeuristicAll = "all"
const finalHeuristicControls = "controls"

var gksDeviationMsaPatchBiases = []float64{0.0, 0.25, 0.5, 0.75, 1.0}
var randomSparseSeeds = []int64{1, 2, 3}

var finalResultsSummaryHeuristics = []string{
	heuristicBaseline,
	heuristicMsaHeuristic,
	heuristicCycleCover,
	heuristicCycleCoverMsaPatching,
}

var msaCountScalingCounts = []int{1, 2, 4, 8, 16, 32, 64, 0}

var tuningInstanceFiles = []string{
	"br17.atsp",
	"ftv33.atsp",
	"p43.atsp",
	"ft53.atsp",
	"ftv64.atsp",
	"crane66_1.atsp",
	"atex5.atsp",
	"ftv90.atsp",
	"crane100_1.atsp",
	"td100_1.atsp",
	"dc112.atsp",
	"ftv120.atsp",
	"dc134.atsp",
	"ftv150.atsp",
	"dc188.atsp",
	"code198.atsp",
}

var evaluationInstanceFiles = []string{
	"atex1.atsp",
	"atex3.atsp",
	"atex4.atsp",
	"ftv35.atsp",
	"ftv38.atsp",
	"ftv44.atsp",
	"ftv47.atsp",
	"ry48p.atsp",
	"ftv55.atsp",
	"crane66_0.atsp",
	"crane66_2.atsp",
	"ft70.atsp",
	"ftv70.atsp",
	"crane100_0.atsp",
	"crane100_2.atsp",
	"ftv100.atsp",
	"ftv110.atsp",
	"dc126.atsp",
	"ftv130.atsp",
	"ftv140.atsp",
	"ftv160.atsp",
	"ftv170.atsp",
	"dc176.atsp",
	"rbg323.atsp",
	"rbg358.atsp",
	"rbg403.atsp",
	"rbg443.atsp",
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

func statisticsCsvHeaderForHeuristic(heuristic string) []string {
	if heuristic == heuristicCycleCoverMsaPatching {
		header := append([]string{}, statisticsCsvHeader[:4]...)
		header = append(header, "MSA patch bias")
		header = append(header, statisticsCsvHeader[4:]...)
		return header
	}

	if heuristic == heuristicRandomSparse {
		header := append([]string{}, statisticsCsvHeader[:4]...)
		header = append(header, "Random seed")
		header = append(header, statisticsCsvHeader[4:]...)
		return header
	}

	return statisticsCsvHeader
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

func saveBestParametersInfo(resultsRootPath, fileName string, bestStatistics []ExperimentsDataStatistics, includeMsaPatchBias bool) {
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
	msaPatchBiasCounts := make(map[float64]int)

	// Count the occurrences
	for params := range uniqueParameters {
		alphaCounts[params.alpha]++
		betaCounts[params.beta]++
		rhoCounts[params.rho]++
		heuristicWeightCounts[params.heuristicWeight]++
		if includeMsaPatchBias {
			msaPatchBiasCounts[params.msaPatchBias]++
		}
	}

	// Markdown content
	markdown := "# Best Parameters Report\n\n"

	// Unique combinations count
	markdown += fmt.Sprintf("Found **%d** best unique parameter combinations.\n\n", len(uniqueParameters))

	// Best parameters
	markdown += "## Best Parameters\n\n"
	if includeMsaPatchBias {
		markdown += "| Alpha | Beta | Rho | Heuristic weight | MSA patch bias | Times used |\n"
		markdown += "|-------|------|-----|------------------|----------------|------------|\n"
	} else {
		markdown += "| Alpha | Beta | Rho | Heuristic weight | Times used |\n"
		markdown += "|-------|------|-----|------------------|------------|\n"
	}

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
		if left.heuristicWeight != right.heuristicWeight {
			return left.heuristicWeight < right.heuristicWeight
		}
		if left.msaPatchBias != right.msaPatchBias {
			return left.msaPatchBias < right.msaPatchBias
		}
		return left.randomSeed < right.randomSeed
	})

	for _, parameters := range sortedParameters {
		timesUsed := uniqueParameters[parameters]
		if includeMsaPatchBias {
			markdown += fmt.Sprintf("| %.2f | %.2f | %.2f | %.2f | %.2f | %d |\n",
				parameters.alpha, parameters.beta, parameters.rho, parameters.heuristicWeight, parameters.msaPatchBias, timesUsed)
		} else {
			markdown += fmt.Sprintf("| %.2f | %.2f | %.2f | %.2f | %d |\n",
				parameters.alpha, parameters.beta, parameters.rho, parameters.heuristicWeight, timesUsed)
		}
	}
	markdown += "\n"

	// Parameter occurrences
	markdown += "## Parameter Values Occurrences\n\n"
	markdown += generateMarkdownCounts("Alpha", alphaCounts)
	markdown += generateMarkdownCounts("Beta", betaCounts)
	markdown += generateMarkdownCounts("Rho", rhoCounts)
	markdown += generateMarkdownCounts("Heuristic weight", heuristicWeightCounts)
	if includeMsaPatchBias {
		markdown += generateMarkdownCounts("MSA patch bias", msaPatchBiasCounts)
	}

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
	if includeMsaPatchBias {
		minMsaPatchBias, maxMsaPatchBias := findMinMax(msaPatchBiasCounts)
		markdown += fmt.Sprintf("- **MSA patch bias**: %.2f - %.2f\n", minMsaPatchBias, maxMsaPatchBias)
	}

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
	includeMsaPatchBias := slices.Equal(header, statisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching))
	includeRandomSeed := slices.Equal(header, statisticsCsvHeaderForHeuristic(heuristicRandomSparse))
	if !slices.Equal(header, statisticsCsvHeader) && !includeMsaPatchBias && !includeRandomSeed {
		return nil, fmt.Errorf("invalid statistics header")
	}
	defaultMsaPatchBias := defaultMsaPatchBiasForStatisticsPath(csvFilePath)

	var statistics []ExperimentsDataStatistics

	// Parse the remaining rows
	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("failed to read records: %w", err)
	}

	for _, record := range records {
		expectedRecordLength := len(statisticsCsvHeader)
		if includeMsaPatchBias {
			expectedRecordLength = len(statisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching))
		} else if includeRandomSeed {
			expectedRecordLength = len(statisticsCsvHeaderForHeuristic(heuristicRandomSparse))
		}
		if len(record) != expectedRecordLength {
			return nil, fmt.Errorf("invalid record length: %d", len(record))
		}

		statistic, err := parseStatisticsRecord(record, includeMsaPatchBias, includeRandomSeed, defaultMsaPatchBias)
		if err != nil {
			return nil, err
		}

		statistics = append(statistics, statistic)
	}

	return statistics, nil
}

func defaultMsaPatchBiasForStatisticsPath(csvFilePath string) float64 {
	if strings.Contains(filepath.Base(csvFilePath), "cycle_cover_msa_patching") {
		return 1.0
	}

	return 0.0
}

func parseStatisticsRecord(record []string, includeMsaPatchBias, includeRandomSeed bool, defaultMsaPatchBias float64) (ExperimentsDataStatistics, error) {
	expectedRecordLength := len(statisticsCsvHeader)
	if includeMsaPatchBias {
		expectedRecordLength = len(statisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching))
	} else if includeRandomSeed {
		expectedRecordLength = len(statisticsCsvHeaderForHeuristic(heuristicRandomSparse))
	}
	if len(record) != expectedRecordLength {
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
	msaPatchBias := defaultMsaPatchBias
	indexOffset := 0
	if includeMsaPatchBias {
		var err error
		msaPatchBias, err = strconv.ParseFloat(record[4], 64)
		if err != nil {
			return ExperimentsDataStatistics{}, fmt.Errorf("invalid MSA patch bias: %w", err)
		}
		indexOffset = 1
	}
	randomSeed := int64(0)
	if includeRandomSeed {
		parsedRandomSeed, err := strconv.ParseInt(record[4+indexOffset], 10, 64)
		if err != nil {
			return ExperimentsDataStatistics{}, fmt.Errorf("invalid random seed: %w", err)
		}
		randomSeed = parsedRandomSeed
		indexOffset++
	}
	iterations, err := strconv.Atoi(record[4+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid iterations: %w", err)
	}
	minBestAtIteration, err := strconv.Atoi(record[5+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min best at iteration: %w", err)
	}
	averageBestAtIteration, err := strconv.ParseFloat(record[6+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg best at iteration: %w", err)
	}
	maxBestAtIteration, err := strconv.Atoi(record[7+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max best at iteration: %w", err)
	}
	minThreeOptImprovementsCount, err := strconv.Atoi(record[8+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min local search improvements: %w", err)
	}
	averageThreeOptImprovementsCount, err := strconv.ParseFloat(record[9+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg local search improvements: %w", err)
	}
	maxThreeOptImprovementsCount, err := strconv.Atoi(record[10+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max local search improvements: %w", err)
	}
	minBestDeviation, err := strconv.ParseFloat(record[11+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min best deviation: %w", err)
	}
	averageBestDeviation, err := strconv.ParseFloat(record[12+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg best deviation: %w", err)
	}
	maxBestDeviation, err := strconv.ParseFloat(record[13+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max best deviation: %w", err)
	}
	successRate, err := strconv.ParseFloat(record[14+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid success rate: %w", err)
	}

	return ExperimentsDataStatistics{
		ExperimentParameters: ExperimentParameters{
			alpha:           alpha,
			beta:            beta,
			rho:             rho,
			heuristicWeight: heuristicWeight,
			msaPatchBias:    msaPatchBias,
			randomSeed:      randomSeed,
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

func saveStatistics(resultCsvPath, heuristic string, statistics []ExperimentsDataStatistics) {
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
	includeMsaPatchBias := heuristic == heuristicCycleCoverMsaPatching
	includeRandomSeed := heuristic == heuristicRandomSparse

	_ = writer.Write(statisticsCsvHeaderForHeuristic(heuristic))

	for _, statistic := range statistics {
		writer.Write(statisticsCsvRecord(statistic, includeMsaPatchBias, includeRandomSeed))
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
		record := append([]string{statistic.heuristic}, statisticsCsvRecord(statistic.statistics, false, false)...)
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
	extendedHeader := append([]string{"Heuristic"}, statisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching)...)
	extendedFormat := slices.Equal(header, extendedHeader)
	if !slices.Equal(header, expectedHeader) && !extendedFormat {
		return nil, fmt.Errorf("invalid heuristic statistics header")
	}

	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	statistics := make([]HeuristicExperimentStatistics, 0, len(records))
	for _, record := range records {
		expectedRecordLength := len(expectedHeader)
		if extendedFormat {
			expectedRecordLength = len(extendedHeader)
		}
		if len(record) != expectedRecordLength {
			return nil, fmt.Errorf("invalid heuristic statistics record length: got %d want %d", len(record), expectedRecordLength)
		}

		parsed, err := parseStatisticsRecord(record[1:], extendedFormat, false, 0.0)
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
		delta := finalAverages[heuristic].averageMinDeviation - finalThreeOptAverages[heuristic].averageMinDeviation
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
		delta := finalThreeOptAverages[heuristic].successRate - finalAverages[heuristic].successRate
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
	patchingPrecision := ratio(structuralTotals.patchingOptimalEdges, structuralTotals.patchingEdges)
	msaRecall := ratio(structuralTotals.msaOptimalEdges, structuralTotals.foundOptimalEdges)
	cycleCoverRecall := ratio(structuralTotals.cycleCoverOptimalEdges, structuralTotals.foundOptimalEdges)
	patchingRecall := ratio(structuralTotals.patchingOptimalEdges, structuralTotals.foundOptimalEdges)

	structuralPrecision := map[string]float64{
		heuristicMsaHeuristic:          msaPrecision,
		heuristicCycleCover:            cycleCoverPrecision,
		heuristicCycleCoverMsaPatching: patchingPrecision,
	}
	structuralRecall := map[string]float64{
		heuristicMsaHeuristic:          msaRecall,
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

		includedHeuristics = append(includedHeuristics, heuristic)
		precisionValues[heuristic] = structuralPrecision[heuristic]
		recallValues[heuristic] = structuralRecall[heuristic]
		deviationValues[heuristic] = metric.averageMinDeviation
		successValues[heuristic] = metric.successRate
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
	builder.WriteString("This table compares the current MSA heuristic, minimum cycle-cover, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.\n\n")
	builder.WriteString("Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.\n\n")
	writeStructuralSimilarityFindings(&builder, totals)
	builder.WriteString("\n")
	writeStructuralSimilarityTable(&builder, rows, totals)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeStructuralSimilarityFindings(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	patchingPrecision := ratio(totals.patchingOptimalEdges, totals.patchingEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)
	patchingRecall := ratio(totals.patchingOptimalEdges, totals.foundOptimalEdges)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **Precision vs found-optimal tours: MSA heuristic %.2f%%, cycle cover %.2f%%, cycle-cover MSA patching %.2f%%.**\n", 100*msaPrecision, 100*cycleCoverPrecision, 100*patchingPrecision)
	fmt.Fprintf(builder, "- **Recall vs found-optimal tours: MSA heuristic %.2f%%, cycle cover %.2f%%, cycle-cover MSA patching %.2f%%.**\n", 100*msaRecall, 100*cycleCoverRecall, 100*patchingRecall)
	fmt.Fprintf(builder, "- **Best-or-tied precision counts: MSA heuristic %d/%d, cycle cover %d/%d, cycle-cover MSA patching %d/%d.**\n",
		totals.msaPrecisionWins,
		totals.instanceCount,
		totals.cycleCoverPrecisionWins,
		totals.instanceCount,
		totals.patchingPrecisionWins,
		totals.instanceCount)
	fmt.Fprintf(builder, "- **Best-or-tied recall counts: MSA heuristic %d/%d, cycle cover %d/%d, cycle-cover MSA patching %d/%d.**\n",
		totals.msaRecallWins,
		totals.instanceCount,
		totals.cycleCoverRecallWins,
		totals.instanceCount,
		totals.patchingRecallWins,
		totals.instanceCount)
}

func writeStructuralSimilarityTable(builder *strings.Builder, rows []cycleCover.InstanceAnalysis, totals structuralSimilaritySummary) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">MSA heuristic</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>\n")
	builder.WriteString("<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>\n")
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
	msaMetrics := metrics.HighMsaHeuristicMetrics
	cycleCoverMetrics := metrics.CycleCoverMetrics
	patchingMetrics := metrics.CycleCoverMsaPatchingMetrics
	precisionHighlights := bestStructuralMetricHighlights(msaMetrics.Precision, cycleCoverMetrics.Precision, patchingMetrics.Precision)
	recallHighlights := bestStructuralMetricHighlights(msaMetrics.Recall, cycleCoverMetrics.Recall, patchingMetrics.Recall)

	writeStructuralSimilarityTableRow(
		builder,
		html.EscapeString(analysis.Instance),
		msaMetrics.Precision,
		msaMetrics.Recall,
		cycleCoverMetrics.Precision,
		cycleCoverMetrics.Recall,
		patchingMetrics.Precision,
		patchingMetrics.Recall,
		precisionHighlights,
		recallHighlights)
}

func writeStructuralSimilarityTotalRow(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	patchingPrecision := ratio(totals.patchingOptimalEdges, totals.patchingEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)
	patchingRecall := ratio(totals.patchingOptimalEdges, totals.foundOptimalEdges)
	precisionHighlights := bestStructuralMetricHighlights(msaPrecision, cycleCoverPrecision, patchingPrecision)
	recallHighlights := bestStructuralMetricHighlights(msaRecall, cycleCoverRecall, patchingRecall)

	writeStructuralSimilarityTableRow(
		builder,
		"<strong>Total</strong>",
		msaPrecision,
		msaRecall,
		cycleCoverPrecision,
		cycleCoverRecall,
		patchingPrecision,
		patchingRecall,
		precisionHighlights,
		recallHighlights)
}

func writeStructuralSimilarityTableRow(builder *strings.Builder, instanceCell string, msaPrecision, msaRecall, cycleCoverPrecision, cycleCoverRecall, patchingPrecision, patchingRecall float64, precisionHighlights, recallHighlights []bool) {
	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		instanceCell,
		finalResultsSummaryMetricCell(100*msaPrecision, precisionHighlights[0]),
		finalResultsSummaryMetricCell(100*msaRecall, recallHighlights[0]),
		finalResultsSummaryMetricCell(100*cycleCoverPrecision, precisionHighlights[1]),
		finalResultsSummaryMetricCell(100*cycleCoverRecall, recallHighlights[1]),
		finalResultsSummaryMetricCell(100*patchingPrecision, precisionHighlights[2]),
		finalResultsSummaryMetricCell(100*patchingRecall, recallHighlights[2]))
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
	patchingEdges           int
	patchingOptimalEdges    int
	msaPrecisionWins        int
	cycleCoverPrecisionWins int
	patchingPrecisionWins   int
	msaRecallWins           int
	cycleCoverRecallWins    int
	patchingRecallWins      int
}

func structuralSimilarityTotals(rows []cycleCover.InstanceAnalysis) structuralSimilaritySummary {
	var totals structuralSimilaritySummary
	totals.instanceCount = len(rows)

	for _, analysis := range rows {
		metrics := analysis.Metrics
		msaMetrics := metrics.HighMsaHeuristicMetrics
		cycleCoverMetrics := metrics.CycleCoverMetrics
		patchingMetrics := metrics.CycleCoverMsaPatchingMetrics
		precisionWins := bestStructuralMetricHighlights(msaMetrics.Precision, cycleCoverMetrics.Precision, patchingMetrics.Precision)
		recallWins := bestStructuralMetricHighlights(msaMetrics.Recall, cycleCoverMetrics.Recall, patchingMetrics.Recall)

		totals.foundOptimalTours += metrics.FoundOptimalTourCount
		totals.foundOptimalEdges += metrics.UniqueFoundOptimalEdgeCount
		totals.tourEdges += analysis.Dimension
		totals.msaEdges += msaMetrics.EdgeCount
		totals.msaOptimalEdges += msaMetrics.OptimalEdgeCount
		totals.cycleCoverEdges += cycleCoverMetrics.EdgeCount
		totals.cycleCoverOptimalEdges += cycleCoverMetrics.OptimalEdgeCount
		totals.patchingEdges += patchingMetrics.EdgeCount
		totals.patchingOptimalEdges += patchingMetrics.OptimalEdgeCount
		if precisionWins[0] {
			totals.msaPrecisionWins++
		}
		if precisionWins[1] {
			totals.cycleCoverPrecisionWins++
		}
		if precisionWins[2] {
			totals.patchingPrecisionWins++
		}
		if recallWins[0] {
			totals.msaRecallWins++
		}
		if recallWins[1] {
			totals.cycleCoverRecallWins++
		}
		if recallWins[2] {
			totals.patchingRecallWins++
		}
	}

	return totals
}

func saveMsaHeuristicCycleCoverOverlapReport(path string, analyses []cycleCover.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows := sortedCycleCoverAnalyses(analyses)
	totals := msaHeuristicCycleCoverOverlapTotals(rows)

	var builder strings.Builder
	builder.WriteString("# MSA Heuristic And Cycle-Cover Overlap\n\n")
	builder.WriteString("This table compares the current MSA heuristic edge set directly with the minimum cycle-cover edge set.\n\n")
	writeMsaHeuristicCycleCoverOverlapFindings(&builder, totals)
	builder.WriteString("\n")
	writeMsaHeuristicCycleCoverOverlapTable(&builder, rows, totals)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

type gksDeviationRow struct {
	instance     string
	msaPatchBias float64
	tourLength   float64
	deviation    float64
}

func saveGksDeviationReport(path string, atspsData []AtspData, msaPatchBiases []float64) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows, err := buildGksDeviationRows(atspsData, msaPatchBiases)
	if err != nil {
		return err
	}

	var builder strings.Builder
	builder.WriteString("# GKS Patching Deviation\n\n")
	builder.WriteString("This table runs the cycle-cover patching heuristic directly, without MMAS. `MSA patch bias = 0.00` is the pure GKS-style cost-based patching variant; positive values bias patch selection toward MSA-supported connector edges.\n\n")
	writeGksDeviationTable(&builder, rows, msaPatchBiases)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildGksDeviationRows(atspsData []AtspData, msaPatchBiases []float64) ([]gksDeviationRow, error) {
	rows := make([]gksDeviationRow, 0, len(atspsData)*len(msaPatchBiases))
	for _, atspData := range atspsData {
		msaHeuristicMatrix, err := msaHeuristic.Read(atspData.msaHeuristicDirectoryPath)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read MSA heuristic: %w", atspData.name, err)
		}

		cycleCoverMatrix, _, err := buildMinimumCycleCoverMatrix(atspData.matrix)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to compute minimum cycle cover: %w", atspData.name, err)
		}

		for _, msaPatchBias := range msaPatchBiases {
			patchingMatrix := heuristics.BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(atspData.matrix, msaHeuristicMatrix, cycleCoverMatrix, msaPatchBias)
			tourLength, err := tourMatrixLength(atspData.matrix, patchingMatrix)
			if err != nil {
				return nil, fmt.Errorf("%s: invalid GKS tour for MSA patch bias %.2f: %w", atspData.name, msaPatchBias, err)
			}

			rows = append(rows, gksDeviationRow{
				instance:     atspData.name,
				msaPatchBias: msaPatchBias,
				tourLength:   tourLength,
				deviation:    100.0 * (tourLength - atspData.knownOptimal) / atspData.knownOptimal,
			})
		}
	}

	sort.Slice(rows, func(i, j int) bool {
		if rows[i].instance != rows[j].instance {
			return rows[i].instance < rows[j].instance
		}
		return rows[i].msaPatchBias < rows[j].msaPatchBias
	})

	return rows, nil
}

func writeGksDeviationTable(builder *strings.Builder, rows []gksDeviationRow, msaPatchBiases []float64) {
	bestBiasesByInstance := bestGksBiasesByInstance(rows)
	averageRows := averageGksDeviationRows(rows, msaPatchBiases)
	averageHighlights := lowestGksDeviationHighlights(averageRows)

	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>MSA patch bias</th><th>Tour length</th><th>Deviation [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		key := gksDeviationRowKey(row.instance, row.msaPatchBias)
		writeGksDeviationTableRow(builder, row.instance, row, bestBiasesByInstance[key])
	}
	for _, row := range averageRows {
		key := gksDeviationRowKey(row.instance, row.msaPatchBias)
		writeGksDeviationTableRow(builder, "Average", row, averageHighlights[key])
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeGksDeviationTableRow(builder *strings.Builder, instance string, row gksDeviationRow, highlightDeviation bool) {
	instanceCell := html.EscapeString(instance)
	if instance == "Average" {
		instanceCell = "<strong>Average</strong>"
	}

	tourLengthCell := fmt.Sprintf("%.2f", row.tourLength)
	if instance == "Average" {
		tourLengthCell = ""
	}

	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		instanceCell,
		row.msaPatchBias,
		tourLengthCell,
		finalResultsSummaryMetricCell(row.deviation, highlightDeviation))
}

func bestGksBiasesByInstance(rows []gksDeviationRow) map[string]bool {
	bestByInstance := make(map[string]float64)
	for _, row := range rows {
		best, ok := bestByInstance[row.instance]
		if !ok || row.deviation < best {
			bestByInstance[row.instance] = row.deviation
		}
	}

	highlights := make(map[string]bool, len(rows))
	for _, row := range rows {
		if math.Abs(row.deviation-bestByInstance[row.instance]) < 1e-9 {
			highlights[gksDeviationRowKey(row.instance, row.msaPatchBias)] = true
		}
	}
	return highlights
}

func averageGksDeviationRows(rows []gksDeviationRow, msaPatchBiases []float64) []gksDeviationRow {
	totals := make(map[float64]float64, len(msaPatchBiases))
	counts := make(map[float64]int, len(msaPatchBiases))
	for _, row := range rows {
		totals[row.msaPatchBias] += row.deviation
		counts[row.msaPatchBias]++
	}

	averageRows := make([]gksDeviationRow, 0, len(msaPatchBiases))
	for _, msaPatchBias := range msaPatchBiases {
		count := counts[msaPatchBias]
		if count == 0 {
			continue
		}

		averageRows = append(averageRows, gksDeviationRow{
			instance:     "Average",
			msaPatchBias: msaPatchBias,
			deviation:    totals[msaPatchBias] / float64(count),
		})
	}
	return averageRows
}

func lowestGksDeviationHighlights(rows []gksDeviationRow) map[string]bool {
	best := math.Inf(1)
	for _, row := range rows {
		if row.deviation < best {
			best = row.deviation
		}
	}

	highlights := make(map[string]bool, len(rows))
	for _, row := range rows {
		if math.Abs(row.deviation-best) < 1e-9 {
			highlights[gksDeviationRowKey(row.instance, row.msaPatchBias)] = true
		}
	}
	return highlights
}

func gksDeviationRowKey(instance string, msaPatchBias float64) string {
	return instance + ":" + fmt.Sprintf("%.6f", msaPatchBias)
}

func tourMatrixLength(distanceMatrix, tourMatrix [][]float64) (float64, error) {
	dimension := len(tourMatrix)
	successors := make([]int, dimension)
	inDegree := make([]int, dimension)
	for vertex := range successors {
		successors[vertex] = -1
	}

	for from, row := range tourMatrix {
		outDegree := 0
		for to, value := range row {
			if value == 0 {
				continue
			}
			if from == to {
				return 0, fmt.Errorf("self-loop at vertex %d", from)
			}
			if from >= len(distanceMatrix) || to >= len(distanceMatrix[from]) {
				return 0, fmt.Errorf("edge %d -> %d is outside distance matrix", from, to)
			}

			outDegree++
			successors[from] = to
			inDegree[to]++
		}
		if outDegree != 1 {
			return 0, fmt.Errorf("vertex %d has out-degree %d", from, outDegree)
		}
	}

	for vertex, degree := range inDegree {
		if degree != 1 {
			return 0, fmt.Errorf("vertex %d has in-degree %d", vertex, degree)
		}
	}

	visited := make([]bool, dimension)
	current := 0
	for step := 0; step < dimension; step++ {
		if current < 0 || current >= dimension {
			return 0, fmt.Errorf("tour left matrix at vertex %d", current)
		}
		if visited[current] {
			return 0, fmt.Errorf("tour returned to vertex %d after %d steps", current, step)
		}
		visited[current] = true
		current = successors[current]
	}
	if current != 0 {
		return 0, fmt.Errorf("tour ended at vertex %d instead of 0", current)
	}

	length := 0.0
	for from, to := range successors {
		length += distanceMatrix[from][to]
	}
	return length, nil
}

func writeMsaHeuristicCycleCoverOverlapFindings(builder *strings.Builder, totals msaHeuristicCycleCoverOverlapSummary) {
	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **%.2f%% of MSA heuristic edges are also cycle-cover edges.**\n", 100*ratio(totals.sharedEdges, totals.msaEdges))
	fmt.Fprintf(builder, "- **%.2f%% of cycle-cover edges are also MSA heuristic edges.**\n", 100*ratio(totals.sharedEdges, totals.cycleCoverEdges))
	fmt.Fprintf(builder,
		"- **Found-optimal edge partition: both %d, only MSA heuristic %d, only cycle cover %d, neither %d.**\n",
		totals.optimalBoth,
		totals.optimalOnlyMsa,
		totals.optimalOnlyCycleCover,
		totals.optimalNeither)
}

func writeMsaHeuristicCycleCoverOverlapTable(builder *strings.Builder, rows []cycleCover.InstanceAnalysis, totals msaHeuristicCycleCoverOverlapSummary) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, analysis := range rows {
		writeMsaHeuristicCycleCoverOverlapRow(builder, analysis)
	}
	writeMsaHeuristicCycleCoverOverlapTotalRow(builder, totals)

	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeMsaHeuristicCycleCoverOverlapRow(builder *strings.Builder, analysis cycleCover.InstanceAnalysis) {
	metrics := analysis.Metrics
	msaEdges := metrics.HighMsaHeuristicMetrics.EdgeCount
	cycleCoverEdges := metrics.CycleCoverMetrics.EdgeCount
	sharedEdges := metrics.CycleCoverHighMsaHeuristicEdges

	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td></tr>\n",
		html.EscapeString(analysis.Instance),
		100*ratio(sharedEdges, msaEdges),
		100*ratio(sharedEdges, cycleCoverEdges),
		metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic,
		metrics.OptimalEdgesInHighMsaHeuristicNotCycleCover,
		metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic)
}

func writeMsaHeuristicCycleCoverOverlapTotalRow(builder *strings.Builder, totals msaHeuristicCycleCoverOverlapSummary) {
	fmt.Fprintf(builder,
		"<tr><td><strong>Total</strong></td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td></tr>\n",
		100*ratio(totals.sharedEdges, totals.msaEdges),
		100*ratio(totals.sharedEdges, totals.cycleCoverEdges),
		totals.optimalBoth,
		totals.optimalOnlyMsa,
		totals.optimalOnlyCycleCover)
}

type msaHeuristicCycleCoverOverlapSummary struct {
	msaEdges              int
	cycleCoverEdges       int
	sharedEdges           int
	optimalBoth           int
	optimalOnlyMsa        int
	optimalOnlyCycleCover int
	optimalNeither        int
}

func msaHeuristicCycleCoverOverlapTotals(rows []cycleCover.InstanceAnalysis) msaHeuristicCycleCoverOverlapSummary {
	var totals msaHeuristicCycleCoverOverlapSummary
	for _, analysis := range rows {
		metrics := analysis.Metrics
		totals.msaEdges += metrics.HighMsaHeuristicMetrics.EdgeCount
		totals.cycleCoverEdges += metrics.CycleCoverMetrics.EdgeCount
		totals.sharedEdges += metrics.CycleCoverHighMsaHeuristicEdges
		totals.optimalBoth += metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic
		totals.optimalOnlyMsa += metrics.OptimalEdgesInHighMsaHeuristicNotCycleCover
		totals.optimalOnlyCycleCover += metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic
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
	builder.WriteString("This table shows how the MSA heuristic signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.\n\n")
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
			boostedEdges := buildPartialMsaHeuristicEdgeSet(msas, selectedRoots)

			row.instanceCount++
			row.boostedEdges += len(boostedEdges)
			row.boostedTargetEdges += maxIntValue(len(msas)-1, 0)

			tours, err := msaHeuristicTours.ReadOptimalTours(atspData.optimalUniqueToursCsvPath)
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
	msas, err := msaHeuristic.ReadMsas(atspData.msaHeuristicDirectoryPath)
	if err == nil && len(msas) == len(atspData.matrix) {
		return msas, nil
	}

	if _, createErr := msaHeuristic.Create(atspData.matrix, atspData.msaHeuristicDirectoryPath); createErr != nil {
		if err != nil {
			return nil, fmt.Errorf("%s: read individual MSAs: %w; create MSA heuristic: %w", atspData.name, err, createErr)
		}
		return nil, fmt.Errorf("%s: create MSA heuristic: %w", atspData.name, createErr)
	}

	msas, err = msaHeuristic.ReadMsas(atspData.msaHeuristicDirectoryPath)
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

func buildPartialMsaHeuristicEdgeSet(msas [][][]float64, selectedRoots []int) map[models.Edge]struct{} {
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
	builder.WriteString("This sanity check compares the final MSA heuristic against deterministic random sparse masks. Each random mask boosts the same number of directed edges as the MSA heuristic and uses the same heuristic weight. The comparison reads MSA from the final results and averages the available final-control random seeds for each instance.\n\n")
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
			missingData = append(missingData, randomSparseControlMissingData{instance: atspData.name, reason: "missing final MSA result"})
			continue
		}

		controlAtspData := withExperimentOutputRoot(atspData, controlResultsRootPath)
		randomStatistics, err := readStatistics(resultFilePathForHeuristic(controlAtspData, heuristicRandomSparse))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, randomSparseControlMissingData{instance: atspData.name, reason: "missing final-control random-sparse result CSV"})
				continue
			}
			return nil, nil, err
		}
		if len(randomStatistics) == 0 {
			missingData = append(missingData, randomSparseControlMissingData{instance: atspData.name, reason: "empty random-sparse result CSV"})
			continue
		}

		randomAverageBestDeviation := 0.0
		randomSuccessRate := 0.0
		randomBestAverageBestDeviation := math.Inf(1)
		randomBestAverageBestDeviationSeed := int64(0)
		for _, statistic := range randomStatistics {
			randomAverageBestDeviation += statistic.averageBestDeviation
			randomSuccessRate += statistic.successRate
			if statistic.averageBestDeviation < randomBestAverageBestDeviation {
				randomBestAverageBestDeviation = statistic.averageBestDeviation
				randomBestAverageBestDeviationSeed = statistic.randomSeed
			}
		}
		randomAverageBestDeviation /= float64(len(randomStatistics))
		randomSuccessRate /= float64(len(randomStatistics))

		rows = append(rows, randomSparseControlRow{
			instance:                           atspData.name,
			msaAverageBestDeviation:            msaMetric.averageMinDeviation,
			randomAverageBestDeviation:         randomAverageBestDeviation,
			randomBestAverageBestDeviation:     randomBestAverageBestDeviation,
			averageBestDeviationDelta:          msaMetric.averageMinDeviation - randomAverageBestDeviation,
			msaSuccessRate:                     msaMetric.successRate,
			randomSuccessRate:                  randomSuccessRate,
			successRateDelta:                   msaMetric.successRate - randomSuccessRate,
			randomSeedCount:                    len(randomStatistics),
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

func readFinalMsaHeuristicControlMetric(atspData AtspData, finalResultsRootPath string) (finalResultsSummaryMetric, bool, error) {
	finalAtspData := withExperimentOutputRoot(atspData, finalResultsRootPath)
	metrics, err := readFinalResultSummaryMetrics(finalAtspData.resultFilePath)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return finalResultsSummaryMetric{}, false, nil
		}
		return finalResultsSummaryMetric{}, false, err
	}

	metric, ok := metrics[heuristicMsaHeuristic]
	return metric, ok, nil
}

func writeRandomSparseControlFindings(builder *strings.Builder, rows []randomSparseControlRow) {
	summary := randomSparseControlSummary(rows)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **MSA had lower average best deviation than the random-sparse mean in %d/%d instances.**\n", summary.msaWins, len(rows))
	fmt.Fprintf(builder, "- **Mean average best deviation: MSA %.2f%%, random sparse %.2f%%, delta %s pp.**\n",
		summary.meanMsaAverageBestDeviation,
		summary.meanRandomAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta))
	fmt.Fprintf(builder, "- **Mean success rate: MSA %.2f%%, random sparse %.2f%%, delta %s pp.**\n",
		summary.meanMsaSuccessRate,
		summary.meanRandomSuccessRate,
		formatSignedFloat(summary.meanSuccessRateDelta))
	fmt.Fprintf(builder, "- MSA also beat the best random seed in %d/%d instances.\n", summary.msaWinsAgainstBestRandomSeed, len(rows))
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
	builder.WriteString("Negative delta means the MSA heuristic had lower average best deviation than the random-sparse mean.\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Random mean avg best dev. [%]</th><th>Best random avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Random success [%]</th><th>Seeds</th></tr>\n")
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
	rows, missingData, err := buildDistanceRankedSparseControlRows(atspsData, finalResultsRootPath, controlResultsRootPath)
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
	builder.WriteString("# Distance-ranked Sparse Control\n\n")
	builder.WriteString("This sanity check compares the final MSA heuristic against a deterministic sparse mask built from the cheapest directed edges. The control boosts the same number of directed edges as the MSA heuristic and uses the same `heuristicWeight=0.90`.\n\n")
	writeDistanceRankedSparseControlFindings(&builder, rows)
	builder.WriteString("\n")
	writeDistanceRankedSparseControlTable(&builder, rows)
	if len(missingData) != 0 {
		writeDistanceRankedSparseControlMissingData(&builder, missingData)
	}

	return true, os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildDistanceRankedSparseControlRows(atspsData []AtspData, finalResultsRootPath, controlResultsRootPath string) ([]distanceRankedSparseControlRow, []distanceRankedSparseControlMissingData, error) {
	rows := make([]distanceRankedSparseControlRow, 0, len(atspsData))
	missingData := make([]distanceRankedSparseControlMissingData, 0)

	for _, atspData := range atspsData {
		msaMetric, ok, err := readFinalMsaHeuristicControlMetric(atspData, finalResultsRootPath)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.name, reason: "missing final MSA result"})
			continue
		}

		controlAtspData := withExperimentOutputRoot(atspData, controlResultsRootPath)
		controlStatistics, err := readStatistics(resultFilePathForHeuristic(controlAtspData, heuristicDistanceRankedSparse))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.name, reason: "missing final-control distance-ranked-sparse result CSV"})
				continue
			}
			return nil, nil, err
		}

		controlStatistic, ok := statisticsForHeuristicWeight(controlStatistics, finalMsaHeuristicWeight)
		if !ok {
			missingData = append(missingData, distanceRankedSparseControlMissingData{instance: atspData.name, reason: "missing distance-ranked-sparse row for heuristic weight 0.90"})
			continue
		}

		rows = append(rows, distanceRankedSparseControlRow{
			instance:                    atspData.name,
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

func writeDistanceRankedSparseControlFindings(builder *strings.Builder, rows []distanceRankedSparseControlRow) {
	summary := distanceRankedSparseControlSummary(rows)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **MSA had lower average best deviation than the distance-ranked sparse control in %d/%d instances.**\n", summary.msaWins, len(rows))
	fmt.Fprintf(builder, "- **Mean average best deviation: MSA %.2f%%, distance-ranked sparse %.2f%%, delta %s pp.**\n",
		summary.meanMsaAverageBestDeviation,
		summary.meanControlAverageBestDeviation,
		formatSignedFloat(summary.meanAverageBestDeviationDelta))
	fmt.Fprintf(builder, "- **Mean success rate: MSA %.2f%%, distance-ranked sparse %.2f%%, delta %s pp.**\n",
		summary.meanMsaSuccessRate,
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

func writeDistanceRankedSparseControlTable(builder *strings.Builder, rows []distanceRankedSparseControlRow) {
	summary := distanceRankedSparseControlSummary(rows)

	builder.WriteString("## Per-instance comparison\n\n")
	builder.WriteString("Negative delta means the MSA heuristic had lower average best deviation than the distance-ranked sparse control.\n\n")
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>MSA avg best dev. [%]</th><th>Distance-ranked avg best dev. [%]</th><th>Delta [pp]</th><th>MSA success [%]</th><th>Distance-ranked success [%]</th></tr>\n")
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

func bestStructuralMetricHighlights(values ...float64) []bool {
	const epsilon = 1e-9
	highlights := make([]bool, len(values))
	best := math.Inf(-1)
	for _, value := range values {
		if value > best {
			best = value
		}
	}
	if best <= 0 {
		return highlights
	}

	for index, value := range values {
		if math.Abs(value-best) < epsilon {
			highlights[index] = true
		}
	}

	return highlights
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
	case heuristicMsaHeuristic:
		return "MSA heuristic"
	case heuristicRandomSparse:
		return "Random sparse"
	case heuristicDistanceRankedSparse:
		return "Distance-ranked sparse"
	case heuristicCycleCover:
		return "Cycle cover"
	case heuristicCycleCoverMsaPatching:
		return "Cycle-cover MSA patching"
	default:
		return heuristic
	}
}

func buildHeuristicModifiers(heuristic string, matrix, msaHeuristic, cycleCover [][]float64, parameters ExperimentParameters) [][]float64 {
	switch heuristic {
	case heuristicBaseline:
		return heuristics.BuildNeutralModifiers(heuristicMatrixDimension(matrix, msaHeuristic, cycleCover))
	case heuristicMsaHeuristic:
		return heuristics.BuildMsaHeuristicModifiers(msaHeuristic, parameters.heuristicWeight)
	case heuristicRandomSparse:
		return heuristics.BuildRandomSparseModifiers(msaHeuristic, parameters.heuristicWeight, parameters.randomSeed)
	case heuristicDistanceRankedSparse:
		return heuristics.BuildDistanceRankedSparseModifiers(matrix, msaHeuristic, parameters.heuristicWeight)
	case heuristicCycleCover:
		return heuristics.BuildCycleCoverModifiers(cycleCover, parameters.heuristicWeight)
	case heuristicCycleCoverMsaPatching:
		return heuristics.BuildCycleCoverMsaPatchingModifiers(matrix, msaHeuristic, cycleCover, parameters.heuristicWeight, parameters.msaPatchBias)
	default:
		return heuristics.BuildNeutralModifiers(heuristicMatrixDimension(matrix, msaHeuristic, cycleCover))
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

func statisticsCsvRecord(statistic ExperimentsDataStatistics, includeMsaPatchBias, includeRandomSeed bool) []string {
	floatFormat := "%.2f"
	record := []string{
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
	if includeMsaPatchBias {
		record = append(record[:4], append([]string{fmt.Sprintf(floatFormat, statistic.msaPatchBias)}, record[4:]...)...)
	} else if includeRandomSeed {
		record = append(record[:4], append([]string{strconv.FormatInt(statistic.randomSeed, 10)}, record[4:]...)...)
	}
	return record
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

func generateParameters(heuristic string) []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0)
	heuristicWeights := []float64{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
	msaPatchBiases := []float64{0.0}
	if heuristic == heuristicCycleCoverMsaPatching {
		msaPatchBiases = []float64{0.0, 0.25, 0.5, 0.75, 1.0}
	}

	for _, alpha := range utilities.GenerateRange(defaultExperimentAlpha, defaultExperimentAlpha, 0.25) {
		for _, beta := range utilities.GenerateRange(defaultExperimentBeta, defaultExperimentBeta, 1.0) {
			for _, rho := range utilities.GenerateRange(defaultExperimentRho, defaultExperimentRho, 0.1) {
				for _, heuristicWeight := range heuristicWeights {
					currentMsaPatchBiases := msaPatchBiases
					if heuristicWeight == 0 {
						currentMsaPatchBiases = []float64{0.0}
					}

					for _, msaPatchBias := range currentMsaPatchBiases {
						parameters = append(parameters,
							ExperimentParameters{
								alpha:           alpha,
								beta:            beta,
								rho:             rho,
								heuristicWeight: heuristicWeight,
								msaPatchBias:    msaPatchBias,
							})
					}
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
			heuristic: heuristicMsaHeuristic,
			parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalMsaHeuristicWeight),
			},
		},
		{
			heuristic: heuristicCycleCover,
			parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalCycleCoverWeight),
			},
		},
		{
			heuristic: heuristicCycleCoverMsaPatching,
			parameters: []ExperimentParameters{
				newPatchingExperimentParameters(finalCycleCoverMsaPatchingWeight, finalCycleCoverMsaPatchingMsaPatchBias),
			},
		},
	}
}

func finalControlExperimentConfigurations() []finalExperimentConfiguration {
	return []finalExperimentConfiguration{
		{
			heuristic:  heuristicRandomSparse,
			parameters: newRandomSparseFinalExperimentParameters(),
		},
		{
			heuristic: heuristicDistanceRankedSparse,
			parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalMsaHeuristicWeight),
			},
		},
	}
}

func selectFinalExperimentConfigurations(finalHeuristic string) ([]finalExperimentConfiguration, error) {
	configurations := finalExperimentConfigurations()
	if finalHeuristic == finalHeuristicAll {
		return configurations, nil
	}
	if finalHeuristic == finalHeuristicControls {
		return finalControlExperimentConfigurations(), nil
	}

	allConfigurations := append(append([]finalExperimentConfiguration{}, configurations...), finalControlExperimentConfigurations()...)
	for _, configuration := range allConfigurations {
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

func newRandomSparseFinalExperimentParameters() []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0, len(randomSparseSeeds))
	for _, randomSeed := range randomSparseSeeds {
		parameter := newDefaultExperimentParameters(finalMsaHeuristicWeight)
		parameter.randomSeed = randomSeed
		parameters = append(parameters, parameter)
	}

	return parameters
}

func newPatchingExperimentParameters(heuristicWeight, msaPatchBias float64) ExperimentParameters {
	parameters := newDefaultExperimentParameters(heuristicWeight)
	parameters.msaPatchBias = msaPatchBias
	return parameters
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

var artifactsDirectoryName = "artifacts"
var resultsDirectoryName = filepath.Join(artifactsDirectoryName, "tuning")
var finalResultsDirectoryName = filepath.Join(artifactsDirectoryName, "final", "no_3opt")
var finalThreeOptResultsDirectoryName = filepath.Join(artifactsDirectoryName, "final", "with_3opt")
var msaHeuristicArtifactsDirectoryName = filepath.Join(artifactsDirectoryName, "msa")
var solutionArtifactsDirectoryName = filepath.Join(artifactsDirectoryName, "solutions")
var resultFileName = "result.csv"

type AtspData struct {
	name         string
	matrix       [][]float64
	knownOptimal float64

	msaHeuristicDirectoryPath string

	msaHeuristicHeatmapPlotPath   string
	msaHeuristicHistogramPlotPath string

	resultFilePath       string
	resultPlotFilePrefix string

	optimalUniqueToursCsvPath               string
	toursHeatmapPlotPath                    string
	toursHistogramPlotPath                  string
	msaHeuristicToursOverlapHeatmapPlotPath string
}

func makeAtspData(name string, matrix [][]float64, knownOptimal float64) AtspData {
	return makeAtspDataInResultsDirectory(name, matrix, knownOptimal, resultsDirectoryName)
}

func makeAtspDataInResultsDirectory(name string, matrix [][]float64, knownOptimal float64, resultsRootPath string) AtspData {
	name = strings.TrimSuffix(name, ".atsp")
	resultsDirectoryPath := filepath.Join(resultsRootPath, name)
	resultsPlotsDirectoryPath := filepath.Join(resultsDirectoryPath, "plots")
	msaHeuristicDirectoryPath := filepath.Join(msaHeuristicArtifactsDirectoryName, name)
	msaHeuristicPlotsDirectoryPath := filepath.Join(msaHeuristicDirectoryPath, "plots")
	solutionsDirectoryPath := filepath.Join(solutionArtifactsDirectoryName, name)
	solutionsPlotsDirectoryPath := filepath.Join(solutionsDirectoryPath, "plots")

	msaHeuristicHeatmapPlotPath := filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_heatmap.png")
	msaHeuristicHistogramPlotPath := filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_histogram.png")

	resultFilePath := filepath.Join(resultsDirectoryPath, resultFileName)
	resultPlotFilePrefix := filepath.Join(resultsPlotsDirectoryPath, "best_result")

	optimalUniqueToursCsvPath := filepath.Join(solutionsDirectoryPath, "solutions.csv")
	toursHeatmapPlotPath := filepath.Join(solutionsPlotsDirectoryPath, "tours_heatmap.png")
	toursHistogramPlotPath := filepath.Join(solutionsPlotsDirectoryPath, "tours_histogram.png")
	msaHeuristicToursOverlapHeatmapPlotPath := filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_tours_overlap_heatmap.png")

	return AtspData{
		name,
		matrix,
		knownOptimal,

		msaHeuristicDirectoryPath,

		msaHeuristicHeatmapPlotPath, msaHeuristicHistogramPlotPath,

		resultFilePath,
		resultPlotFilePrefix,

		optimalUniqueToursCsvPath,
		toursHeatmapPlotPath,
		toursHistogramPlotPath,
		msaHeuristicToursOverlapHeatmapPlotPath,
	}
}

func withExperimentOutputRoot(atspData AtspData, resultsRootPath string) AtspData {
	output := makeAtspDataInResultsDirectory(atspData.name, atspData.matrix, atspData.knownOptimal, resultsRootPath)
	output.msaHeuristicDirectoryPath = atspData.msaHeuristicDirectoryPath
	output.optimalUniqueToursCsvPath = atspData.optimalUniqueToursCsvPath
	return output
}

func selectAtspFiles(atspFilePaths []string, instanceSet string) ([]string, error) {
	switch instanceSet {
	case instanceSetTuning:
		return selectConfiguredAtspFiles(atspFilePaths, tuningInstanceFiles)
	case instanceSetEvaluation:
		return selectConfiguredAtspFiles(atspFilePaths, evaluationInstanceFiles)
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
		return nil, fmt.Errorf("unsupported -instances value %q; use %q, %q, or %q", instanceSet, instanceSetTuning, instanceSetEvaluation, instanceSetAllKnown)
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
		heuristic == heuristicMsaHeuristic ||
		heuristic == heuristicCycleCover ||
		heuristic == heuristicCycleCoverMsaPatching
}

func heuristicUsesCycleCover(heuristic string) bool {
	return heuristic == heuristicCycleCover || heuristic == heuristicCycleCoverMsaPatching
}

func heuristicUsesMsaHeuristic(heuristic string) bool {
	return heuristic == heuristicMsaHeuristic || heuristic == heuristicCycleCoverMsaPatching || heuristicIsSparseControl(heuristic)
}

func heuristicIsRandomSparse(heuristic string) bool {
	return heuristic == heuristicRandomSparse
}

func heuristicIsSparseControl(heuristic string) bool {
	return heuristic == heuristicRandomSparse || heuristic == heuristicDistanceRankedSparse
}

func finalConfigurationsUseMsaHeuristic(configurations []finalExperimentConfiguration) bool {
	for _, configuration := range configurations {
		if heuristicUsesMsaHeuristic(configuration.heuristic) {
			return true
		}
	}

	return false
}

func heuristicFileSuffix(heuristic string) string {
	switch heuristic {
	case heuristicBaseline:
		return "_baseline"
	case heuristicMsaHeuristic:
		return ""
	case heuristicCycleCover:
		return "_cycle_cover"
	case heuristicCycleCoverMsaPatching:
		return "_cycle_cover_msa_patching"
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

func isValidAnalysisScope(scope string) bool {
	return scope == analysisScopeAll || scope == analysisScopeGksDeviation
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

func finalControlsResultsRootPath(finalResultsRootPath string) string {
	return filepath.Join(finalResultsRootPath, "controls")
}

func finalExperimentOutputRootForConfigurations(mode string, configurations []finalExperimentConfiguration) string {
	resultsRootPath := finalExperimentOutputRoot(mode)
	if finalConfigurationsAreSparseControls(configurations) {
		return finalControlsResultsRootPath(resultsRootPath)
	}

	return resultsRootPath
}

func finalConfigurationsAreSparseControls(configurations []finalExperimentConfiguration) bool {
	if len(configurations) == 0 {
		return false
	}

	for _, configuration := range configurations {
		if !heuristicIsSparseControl(configuration.heuristic) {
			return false
		}
	}

	return true
}

func finalExperimentUsesThreeOpt(mode string) bool {
	return mode == runModeFinal3Opt
}

func selectedInstanceSetForMode(mode, requestedInstanceSet string, instanceSetExplicit bool) string {
	if shouldRunFinalExperiments(mode) && !instanceSetExplicit {
		return instanceSetEvaluation
	}

	return requestedInstanceSet
}

func main() {
	instances := flag.String("instances", instanceSetTuning, "ATSP instance set to run: tuning, evaluation, or all-known")
	mode := flag.String("mode", runModeExperiment, "Run mode: experiment, analyze, all, final, or final+3opt")
	analysisScope := flag.String("analysis", analysisScopeAll, "Analysis scope for analyze mode: all or gks-deviation")
	heuristic := flag.String("heuristic", heuristicMsaHeuristic, "ACO heuristic modifier to use in experiment mode: baseline, msa-heuristic, cycle-cover, or cycle-cover-msa-patching")
	finalHeuristic := flag.String("final-heuristic", finalHeuristicAll, "Final-mode heuristic to run: all, controls, baseline, msa-heuristic, random-sparse, distance-ranked-sparse, cycle-cover, or cycle-cover-msa-patching")
	flag.Parse()

	instanceSetExplicit := false
	flag.Visit(func(flag *flag.Flag) {
		if flag.Name == "instances" {
			instanceSetExplicit = true
		}
	})

	selectedHeuristic := *heuristic
	selectedFinalHeuristic := *finalHeuristic

	if !isValidRunMode(*mode) {
		fmt.Printf("Unsupported -mode value %q; use %q, %q, %q, %q, or %q\n", *mode, runModeExperiment, runModeAnalyze, runModeAll, runModeFinal, runModeFinal3Opt)
		return
	}

	if !isValidAnalysisScope(*analysisScope) {
		fmt.Printf("Unsupported -analysis value %q; use %q or %q\n", *analysisScope, analysisScopeAll, analysisScopeGksDeviation)
		return
	}

	if !isValidHeuristic(selectedHeuristic) {
		fmt.Printf("Unsupported -heuristic value %q; use %q, %q, %q, or %q\n", *heuristic, heuristicBaseline, heuristicMsaHeuristic, heuristicCycleCover, heuristicCycleCoverMsaPatching)
		return
	}

	selectedInstances := selectedInstanceSetForMode(*mode, *instances, instanceSetExplicit)
	var finalConfigurations []finalExperimentConfiguration
	if shouldRunFinalExperiments(*mode) {
		configurations, err := selectFinalExperimentConfigurations(selectedFinalHeuristic)
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

		err = runExperimentMode(atspsData, selectedHeuristic)
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

		err = runFinalExperimentMode(atspsData, finalExperimentOutputRootForConfigurations(*mode, finalConfigurations), finalExperimentUsesThreeOpt(*mode), finalConfigurations)
		stopProfiling()
		if err != nil {
			fmt.Println(err)
			return
		}
	}

	if shouldRunAnalysis(*mode) {
		if err := runAnalysisMode(atspsData, *analysisScope); err != nil {
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
	if err := ensureMsaHeuristicArtifacts(atspsData); err != nil {
		return err
	}

	return runExperimentSet(atspsData, resultsDirectoryName, heuristic, generateParameters(heuristic), defaultExperimentRunCount)
}

func runFinalExperimentMode(atspsData []AtspData, resultsRootPath string, useThreeOpt bool, configurations []finalExperimentConfiguration) error {
	needsMsaHeuristic := finalConfigurationsUseMsaHeuristic(configurations)
	if needsMsaHeuristic {
		if err := ensureMsaHeuristicCache(atspsData); err != nil {
			return err
		}
	}

	if err := removeLegacyFinalReports(resultsRootPath); err != nil {
		return err
	}

	for _, atspData := range atspsData {
		finalAtspData := withExperimentOutputRoot(atspData, resultsRootPath)
		if err := runFinalExperimentForInstance(finalAtspData, useThreeOpt, configurations, needsMsaHeuristic); err != nil {
			return err
		}
	}

	return nil
}

func runFinalExperimentForInstance(atspData AtspData, useThreeOpt bool, configurations []finalExperimentConfiguration, needsMsaHeuristic bool) error {
	matrix := atspData.matrix
	knownOptimal := atspData.knownOptimal
	dimension := len(matrix)
	instanceStart := time.Now()
	controlRun := finalConfigurationsAreSparseControls(configurations)
	finalRunName := "final"
	if useThreeOpt {
		finalRunName = "final+3opt"
	}

	if len(configurations) == 0 {
		return fmt.Errorf("no final experiment configurations selected")
	}

	if controlRun {
		if err := removeFileIfExists(atspData.resultFilePath); err != nil {
			return err
		}
	} else {
		if err := removeLegacyFinalResultFiles(atspData); err != nil {
			return err
		}
	}

	var heuristicMatrix [][]float64
	if needsMsaHeuristic {
		var err error
		heuristicMatrix, err = readMsaHeuristicMatrixForHeuristic(atspData, heuristicMsaHeuristic)
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
			heuristicModifiers := buildHeuristicModifiers(config.heuristic, matrix, heuristicMatrix, cycleCover, parameters)
			results := runExperiments(finalNumberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, useThreeOpt)
			data := ExperimentsData{parameters, results}
			experimentData = append(experimentData, data)

			parameterStatistics := calculateStatistics([]ExperimentsData{data})
			if len(parameterStatistics) != 0 {
				statistic := parameterStatistics[0]
				randomSeedLog := ""
				if parameters.randomSeed != 0 {
					randomSeedLog = fmt.Sprintf(" randomSeed=%d", parameters.randomSeed)
				}
				fmt.Printf("\t%s heuristicWeight=%.2f msaPatchBias=%.2f%s iterations=%d runs=%d elapsed=%s min deviation=%.2f avg deviation=%.2f\n",
					config.heuristic,
					parameters.heuristicWeight,
					parameters.msaPatchBias,
					randomSeedLog,
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

		if controlRun {
			saveStatistics(resultFilePathForHeuristic(atspData, config.heuristic), config.heuristic, statistics)
			continue
		}

		finalStatistics = append(finalStatistics, HeuristicExperimentStatistics{
			heuristic:  config.heuristic,
			statistics: statistics[0],
		})

		if err := removeExperimentPlotsForHeuristic(atspData, config.heuristic); err != nil {
			return err
		}
	}

	if !controlRun {
		if err := saveFinalHeuristicStatistics(atspData.resultFilePath, finalStatistics, configurations); err != nil {
			return err
		}
	}

	fmt.Printf("Finished %s %s in %s\n", finalRunName, atspData.name, time.Since(instanceStart).Round(time.Millisecond))
	return nil
}

func removeLegacyFinalReports(resultsRootPath string) error {
	legacyReports := []string{
		filepath.Join(resultsRootPath, bestParametersReportPathForHeuristic(heuristicBaseline)),
		filepath.Join(resultsRootPath, bestParametersReportPathForHeuristic(heuristicMsaHeuristic)),
		filepath.Join(resultsRootPath, bestParametersReportPathForHeuristic(heuristicCycleCover)),
		filepath.Join(resultsRootPath, bestParametersReportPathForHeuristic(heuristicCycleCoverMsaPatching)),
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
		resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching),
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

		heuristicMatrix, err := readMsaHeuristicMatrixForHeuristic(atspData, heuristic)
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
			heuristicModifiers := buildHeuristicModifiers(heuristic, matrix, heuristicMatrix, cycleCover, parameters)
			results := runExperiments(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, false)
			data := ExperimentsData{parameters, results}

			experimentData = append(experimentData, data)

			parameterStatistics := calculateStatistics([]ExperimentsData{data})
			if len(parameterStatistics) != 0 {
				statistic := parameterStatistics[0]
				randomSeedLog := ""
				if parameters.randomSeed != 0 {
					randomSeedLog = fmt.Sprintf(" randomSeed=%d", parameters.randomSeed)
				}
				fmt.Printf("\theuristicWeight=%.2f msaPatchBias=%.2f%s iterations=%d runs=%d elapsed=%s min deviation=%.2f avg deviation=%.2f\n",
					parameters.heuristicWeight,
					parameters.msaPatchBias,
					randomSeedLog,
					parameters.iterations,
					numberOfExperiments,
					time.Since(parameterStart).Round(time.Millisecond),
					statistic.minBestDeviation,
					statistic.averageBestDeviation)
			}
		}

		statistics := calculateStatistics(experimentData)
		if len(statistics) != 0 {
			saveStatistics(resultFilePathForHeuristic(atspData, heuristic), heuristic, statistics)
			if !heuristicIsSparseControl(heuristic) {
				if err := removeExperimentPlotsForHeuristic(atspData, heuristic); err != nil {
					return err
				}
				saveExperimentPlots(statistics, "MMAS deviation per iteration", resultPlotFilePrefixForHeuristic(atspData, heuristic))
			}
		}

		if !heuristicIsSparseControl(heuristic) {
			uniqueOptimalTours, err := msaHeuristicTours.ReadOptimalTours(atspData.optimalUniqueToursCsvPath)
			if err != nil {
				return err
			}

			for _, data := range experimentData {
				for _, result := range data.results {
					if result.deviationPerIteration[result.bestAtIteration] == 0.0 {
						msaHeuristicTours.AddUniqueTour(uniqueOptimalTours, result.bestTour)
					}
				}
			}

			if err := msaHeuristicTours.SaveOptimalToursStatistics(atspData.optimalUniqueToursCsvPath, atspData.msaHeuristicDirectoryPath, uniqueOptimalTours); err != nil {
				return err
			}
		}

		fmt.Printf("Finished %s in %s\n", atspData.name, time.Since(instanceStart).Round(time.Millisecond))
	}

	if !heuristicIsSparseControl(heuristic) {
		bestStatistics := getBestStatisticsFromFiles(resultFilePathsForHeuristic(atspsData, heuristic))
		saveBestParametersInfo(resultsRootPath, bestParametersReportPathForHeuristic(heuristic), bestStatistics, heuristic == heuristicCycleCoverMsaPatching)
	}

	return nil
}

func readMsaHeuristicMatrixForHeuristic(atspData AtspData, heuristic string) ([][]float64, error) {
	return msaHeuristic.Read(atspData.msaHeuristicDirectoryPath)
}

func ensureMsaHeuristicArtifacts(atspsData []AtspData) error {
	for _, atspData := range atspsData {
		name := atspData.name
		matrix := atspData.matrix
		msaHeuristicDirectoryPath := atspData.msaHeuristicDirectoryPath

		msaHeuristicMatrix, err := msaHeuristic.Read(atspData.msaHeuristicDirectoryPath)

		if err != nil {
			start := time.Now()
			msaHeuristicMatrix, err = msaHeuristic.Create(matrix, msaHeuristicDirectoryPath)
			elapsed := time.Since(start)

			fmt.Printf("\tCreating %s took: %d ms\n", msaHeuristicDirectoryPath, elapsed.Milliseconds())

			if err != nil {
				return fmt.Errorf("error saving MSA heuristic: %w", err)
			}
		}

		msaHeuristicHeatmapPlotTitle := name + " MSA heuristic heatmap"

		err = utilities.SaveHeatmapFromMatrix(msaHeuristicMatrix, msaHeuristicHeatmapPlotTitle, atspData.msaHeuristicHeatmapPlotPath)
		if err != nil {
			return err
		}

		dataForHistogram := filterZeroes(flattenMatrix(msaHeuristicMatrix))
		msaHeuristicHistogramPlotTitle := name + " MSA heuristic histogram"

		dimension := len(matrix)
		err = utilities.SaveHistogramFromData(dataForHistogram, dimension-1, msaHeuristicHistogramPlotTitle, atspData.msaHeuristicHistogramPlotPath)
		if err != nil {
			return err
		}
	}

	return nil
}

func ensureMsaHeuristicCache(atspsData []AtspData) error {
	for _, atspData := range atspsData {
		if _, err := msaHeuristic.Read(atspData.msaHeuristicDirectoryPath); err == nil {
			continue
		}

		start := time.Now()
		if _, err := msaHeuristic.Create(atspData.matrix, atspData.msaHeuristicDirectoryPath); err != nil {
			return fmt.Errorf("error saving MSA heuristic: %w", err)
		}

		fmt.Printf("\tCreating %s took: %d ms\n", atspData.msaHeuristicDirectoryPath, time.Since(start).Milliseconds())
	}

	return nil
}

func runAnalysisMode(atspsData []AtspData, analysisScope string) error {
	gksDeviationReportPath := filepath.Join(finalResultsDirectoryName, "gks_deviation.md")
	gksDeviationAtspData, err := loadSelectedAtspData(instanceSetAllKnown)
	if err != nil {
		return err
	}

	if analysisScope == analysisScopeGksDeviation {
		if err := ensureMsaHeuristicCache(gksDeviationAtspData); err != nil {
			return err
		}
		if err := saveGksDeviationReport(gksDeviationReportPath, gksDeviationAtspData, gksDeviationMsaPatchBiases); err != nil {
			return err
		}
		fmt.Printf("GKS deviation report saved to %s using %d all-known instance(s)\n", gksDeviationReportPath, len(gksDeviationAtspData))
		return nil
	}

	if err := ensureMsaHeuristicArtifacts(atspsData); err != nil {
		return err
	}

	msaHeuristicTourConfigs := make([]msaHeuristicTours.InstanceConfig, 0, len(atspsData))
	cycleCoverConfigs := make([]cycleCover.InstanceConfig, 0, len(atspsData))
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

		cycleCoverConfigs = append(cycleCoverConfigs, cycleCover.InstanceConfig{
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

	cycleCoverAnalyses, err := cycleCover.AnalyzeInstances(cycleCover.Config{
		Instances:     cycleCoverConfigs,
		HighThreshold: 1.0,
		MsaPatchBias:  finalCycleCoverMsaPatchingMsaPatchBias,
	})
	if err != nil {
		return err
	}

	structuralSimilarityReportPath := filepath.Join(finalResultsDirectoryName, "structural_similarity.md")
	if err := saveStructuralSimilarityReport(structuralSimilarityReportPath, cycleCoverAnalyses); err != nil {
		return err
	}

	heuristicOverlapReportPath := filepath.Join(finalResultsDirectoryName, "msa_cycle_cover_overlap.md")
	if err := saveMsaHeuristicCycleCoverOverlapReport(heuristicOverlapReportPath, cycleCoverAnalyses); err != nil {
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
	if err := removeFileIfExists(filepath.Join(finalResultsDirectoryName, "random_sparse_control.md")); err != nil {
		return err
	}
	if err := removeFileIfExists(filepath.Join(finalResultsDirectoryName, "distance_ranked_sparse_control.md")); err != nil {
		return err
	}

	if err := ensureMsaHeuristicCache(gksDeviationAtspData); err != nil {
		return err
	}

	if err := saveGksDeviationReport(gksDeviationReportPath, gksDeviationAtspData, gksDeviationMsaPatchBiases); err != nil {
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

	fmt.Printf("Structural similarity report saved to %s\n", structuralSimilarityReportPath)
	fmt.Printf("MSA heuristic/cycle-cover overlap report saved to %s\n", heuristicOverlapReportPath)
	fmt.Printf("MSA count scaling report saved to %s\n", msaCountScalingReportPath)
	if randomSparseControlReportSaved {
		fmt.Printf("Random sparse control report saved to %s\n", randomSparseControlReportPath)
	}
	if distanceRankedSparseControlReportSaved {
		fmt.Printf("Distance-ranked sparse control report saved to %s\n", distanceRankedSparseControlReportPath)
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
