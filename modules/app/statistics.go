package app

import (
	"atsp_aco_msa/modules/experiments"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
)

var statisticsCsvHeader = experiments.StatisticsCsvHeader

func statisticsCsvHeaderForHeuristic(heuristic string) []string {
	return experiments.StatisticsCsvHeaderForHeuristic(heuristic)
}

func readStatistics(csvFilePath string) ([]ExperimentsDataStatistics, error) {
	return experiments.ReadStatistics(csvFilePath)
}

func parseStatisticsRecord(record []string, includeMsaPatchBias, includeRandomSeed bool, defaultMsaPatchBias float64) (ExperimentsDataStatistics, error) {
	return experiments.ParseStatisticsRecord(record, includeMsaPatchBias, includeRandomSeed, defaultMsaPatchBias)
}

func saveStatistics(resultCsvPath, heuristic string, statistics []ExperimentsDataStatistics) {
	experiments.SaveStatistics(resultCsvPath, heuristic, statistics)
}

func saveHeuristicStatistics(resultCsvPath string, statistics []HeuristicExperimentStatistics) error {
	return experiments.SaveHeuristicStatistics(resultCsvPath, statistics)
}

func readHeuristicStatistics(resultCsvPath string) ([]HeuristicExperimentStatistics, error) {
	return experiments.ReadHeuristicStatistics(resultCsvPath)
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

func mergeHeuristicStatistics(existingStatistics, newStatistics []HeuristicExperimentStatistics, configurations []finalExperimentConfiguration) []HeuristicExperimentStatistics {
	replaceSet := make(map[string]struct{}, len(configurations))
	for _, configuration := range configurations {
		replaceSet[configuration.Heuristic] = struct{}{}
	}

	knownFinalHeuristics := make(map[string]struct{}, len(finalExperimentConfigurations()))
	for _, configuration := range finalExperimentConfigurations() {
		knownFinalHeuristics[configuration.Heuristic] = struct{}{}
	}

	merged := make([]HeuristicExperimentStatistics, 0, len(existingStatistics)+len(newStatistics))
	for _, statistic := range existingStatistics {
		if _, known := knownFinalHeuristics[statistic.Heuristic]; !known {
			continue
		}
		if _, replace := replaceSet[statistic.Heuristic]; !replace {
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
			AverageMinDeviation:  averageDeviation,
			SuccessRate:          successRate,
			AverageBestIteration: averageBestIteration,
			HeuristicWeight:      heuristicWeight,
			Iterations:           iterations,
		}
		heuristic := record[heuristicIndex]
		if current, ok := metrics[heuristic]; !ok || finalResultsSummaryMetricIsBetter(metric, current) {
			metrics[heuristic] = metric
		}
	}

	return metrics, nil
}

func finalResultsSummaryMetricIsBetter(candidate, current finalResultsSummaryMetric) bool {
	if candidate.AverageMinDeviation != current.AverageMinDeviation {
		return candidate.AverageMinDeviation < current.AverageMinDeviation
	}
	if candidate.SuccessRate != current.SuccessRate {
		return candidate.SuccessRate > current.SuccessRate
	}
	if candidate.AverageBestIteration != current.AverageBestIteration {
		return candidate.AverageBestIteration < current.AverageBestIteration
	}
	return candidate.HeuristicWeight < current.HeuristicWeight
}

func indexOf(values []string, value string) int {
	for i, candidate := range values {
		if candidate == value {
			return i
		}
	}

	return -1
}
