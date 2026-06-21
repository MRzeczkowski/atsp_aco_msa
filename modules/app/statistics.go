package app

import (
	"atsp_aco_msa/modules/experiments"
	"os"
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

func saveEvaluationHeuristicStatistics(resultCsvPath string, statistics []HeuristicExperimentStatistics, configurations []evaluationExperimentConfiguration) error {
	if len(configurations) == len(evaluationExperimentConfigurations()) {
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

func mergeHeuristicStatistics(existingStatistics, newStatistics []HeuristicExperimentStatistics, configurations []evaluationExperimentConfiguration) []HeuristicExperimentStatistics {
	replaceSet := make(map[string]struct{}, len(configurations))
	for _, configuration := range configurations {
		replaceSet[configuration.Heuristic] = struct{}{}
	}

	knownEvaluationHeuristics := make(map[string]struct{}, len(evaluationExperimentConfigurations()))
	for _, configuration := range evaluationExperimentConfigurations() {
		knownEvaluationHeuristics[configuration.Heuristic] = struct{}{}
	}

	merged := make([]HeuristicExperimentStatistics, 0, len(existingStatistics)+len(newStatistics))
	for _, statistic := range existingStatistics {
		if _, known := knownEvaluationHeuristics[statistic.Heuristic]; !known {
			continue
		}
		if _, replace := replaceSet[statistic.Heuristic]; !replace {
			merged = append(merged, statistic)
		}
	}

	merged = append(merged, newStatistics...)
	return merged
}
