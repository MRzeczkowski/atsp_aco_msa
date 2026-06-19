package app

import "atsp_aco_msa/modules/experiments"

func statisticsCsvRecord(statistic ExperimentsDataStatistics, includeMsaPatchBias, includeRandomSeed bool) []string {
	return experiments.StatisticsCsvRecord(statistic, includeMsaPatchBias, includeRandomSeed)
}

func calculateStatistics(experimentsData []ExperimentsData) []ExperimentsDataStatistics {
	return experiments.CalculateStatistics(experimentsData)
}

func runExperiments(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, heuristicModifiers [][]float64, useThreeOpt bool) []ExperimentResult {
	return experiments.RunExperiments(numberOfRuns, parameters, knownOptimal, matrix, heuristicModifiers, useThreeOpt)
}

func runExperimentsWithRootedHeuristicModifiers(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, heuristicModifiers [][]float64, rootedHeuristicModifiers [][][]float64, useThreeOpt bool) []ExperimentResult {
	return experiments.RunExperimentsWithRootedHeuristicModifiers(numberOfRuns, parameters, knownOptimal, matrix, heuristicModifiers, rootedHeuristicModifiers, useThreeOpt)
}
