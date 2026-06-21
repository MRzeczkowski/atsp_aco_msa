package app

import (
	"atsp_aco_msa/modules/experiments"
	"os"
)

func saveEvaluationHeuristicStatistics(resultCsvPath string, statistics []HeuristicExperimentStatistics, configurations []evaluationExperimentConfiguration) error {
	if len(configurations) == len(evaluationExperimentConfigurations()) {
		return experiments.SaveHeuristicStatistics(resultCsvPath, statistics)
	}

	existingStatistics, err := experiments.ReadHeuristicStatistics(resultCsvPath)
	if err != nil {
		if os.IsNotExist(err) {
			existingStatistics = nil
		} else {
			return err
		}
	}

	mergedStatistics := mergeHeuristicStatistics(existingStatistics, statistics, configurations)
	return experiments.SaveHeuristicStatistics(resultCsvPath, mergedStatistics)
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
