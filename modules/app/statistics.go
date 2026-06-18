package app

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
)

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

	if heuristicUsesRandomSeedColumn(heuristic) {
		header := append([]string{}, statisticsCsvHeader[:4]...)
		header = append(header, "Random seed")
		header = append(header, statisticsCsvHeader[4:]...)
		return header
	}

	return statisticsCsvHeader
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
	includeRandomSeed := heuristicUsesRandomSeedColumn(heuristic)

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
	heuristicWeightIndex := indexOf(header, "Heuristic weight")
	iterationsIndex := indexOf(header, "Iterations")
	averageBestIterationIndex := indexOf(header, "Avg best at iteration")
	averageDeviationIndex := indexOf(header, "Avg best deviation")
	successRateIndex := indexOf(header, "Success rate [%]")
	if heuristicIndex == -1 || heuristicWeightIndex == -1 || iterationsIndex == -1 || averageBestIterationIndex == -1 || averageDeviationIndex == -1 || successRateIndex == -1 {
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

		heuristicWeight, err := strconv.ParseFloat(record[heuristicWeightIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid heuristic weight for %s: %w", record[heuristicIndex], err)
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

		metric := finalResultsSummaryMetric{
			averageMinDeviation:  averageDeviation,
			successRate:          successRate,
			averageBestIteration: averageBestIteration,
			heuristicWeight:      heuristicWeight,
			iterations:           iterations,
		}
		heuristic := record[heuristicIndex]
		if current, ok := metrics[heuristic]; !ok || finalResultsSummaryMetricIsBetter(metric, current) {
			metrics[heuristic] = metric
		}
	}

	return metrics, nil
}

func finalResultsSummaryMetricIsBetter(candidate, current finalResultsSummaryMetric) bool {
	if candidate.averageMinDeviation != current.averageMinDeviation {
		return candidate.averageMinDeviation < current.averageMinDeviation
	}
	if candidate.successRate != current.successRate {
		return candidate.successRate > current.successRate
	}
	if candidate.averageBestIteration != current.averageBestIteration {
		return candidate.averageBestIteration < current.averageBestIteration
	}
	return candidate.heuristicWeight < current.heuristicWeight
}

func indexOf(values []string, value string) int {
	for i, candidate := range values {
		if candidate == value {
			return i
		}
	}

	return -1
}
