package experiments

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"fmt"
	"math"
	"sort"
	"strconv"
)

func StatisticsCsvRecord(statistic ExperimentsDataStatistics, includeMsaPatchBias, includeRandomSeed bool) []string {
	floatFormat := "%.2f"
	record := []string{
		fmt.Sprintf(floatFormat, statistic.Alpha),
		fmt.Sprintf(floatFormat, statistic.Beta),
		fmt.Sprintf(floatFormat, statistic.Rho),
		fmt.Sprintf(floatFormat, statistic.HeuristicWeight),
		strconv.Itoa(statistic.Iterations),
		strconv.Itoa(statistic.MinBestAtIteration),
		fmt.Sprintf(floatFormat, statistic.AverageBestAtIteration),
		strconv.Itoa(statistic.MaxBestAtIteration),
		strconv.Itoa(statistic.MinThreeOptImprovementsCount),
		fmt.Sprintf(floatFormat, statistic.AverageThreeOptImprovementsCount),
		strconv.Itoa(statistic.MaxThreeOptImprovementsCount),
		fmt.Sprintf(floatFormat, statistic.MinBestDeviation),
		fmt.Sprintf(floatFormat, statistic.AverageBestDeviation),
		fmt.Sprintf(floatFormat, statistic.MaxBestDeviation),
		fmt.Sprintf(floatFormat, statistic.SuccessRate),
	}
	if includeMsaPatchBias {
		record = append(record[:4], append([]string{fmt.Sprintf(floatFormat, statistic.MsaPatchBias)}, record[4:]...)...)
	} else if includeRandomSeed {
		record = append(record[:4], append([]string{strconv.FormatInt(statistic.RandomSeed, 10)}, record[4:]...)...)
	}
	return record
}

func CalculateStatistics(experimentsData []ExperimentsData) []ExperimentsDataStatistics {
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
		minDeviationPerIteration := make([]float64, data.Iterations)
		averageDeviationPerIteration := make([]float64, data.Iterations)
		maxDeviationPerIteration := make([]float64, data.Iterations)

		resultsLen := float64(len(data.Results))
		for _, result := range data.Results {
			if result.BestAtIteration < minBestAtIteration {
				minBestAtIteration = result.BestAtIteration
			}

			averageBestAtIteration += float64(result.BestAtIteration)

			if result.BestAtIteration > maxBestAtIteration {
				maxBestAtIteration = result.BestAtIteration
			}

			if result.ThreeOptImprovementsCount < minThreeOptImprovementsCount {
				minThreeOptImprovementsCount = result.ThreeOptImprovementsCount
			}

			averageThreeOptImprovementsCount += float64(result.ThreeOptImprovementsCount)

			if result.ThreeOptImprovementsCount > maxThreeOptImprovementsCount {
				maxThreeOptImprovementsCount = result.ThreeOptImprovementsCount
			}

			bestDeviation := result.DeviationPerIteration[result.BestAtIteration]

			if bestDeviation < minBestDeviation {
				minBestDeviation = bestDeviation
				copy(minDeviationPerIteration, result.DeviationPerIteration)
			}

			averageBestDeviation += bestDeviation
			for i, deviation := range result.DeviationPerIteration {
				averageDeviationPerIteration[i] += deviation / resultsLen
			}

			if bestDeviation > maxBestDeviation {
				maxBestDeviation = bestDeviation
				copy(maxDeviationPerIteration, result.DeviationPerIteration)
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
			MinBestAtIteration:               minBestAtIteration,
			AverageBestAtIteration:           averageBestAtIteration,
			MaxBestAtIteration:               maxBestAtIteration,
			MinThreeOptImprovementsCount:     minThreeOptImprovementsCount,
			AverageThreeOptImprovementsCount: averageThreeOptImprovementsCount,
			MaxThreeOptImprovementsCount:     maxThreeOptImprovementsCount,
			MinBestDeviation:                 minBestDeviation,
			AverageBestDeviation:             averageBestDeviation,
			MaxBestDeviation:                 maxBestDeviation,
			SuccessRate:                      successRate,
			MinDeviationPerIteration:         minDeviationPerIteration,
			AverageDeviationPerIteration:     averageDeviationPerIteration,
			MaxDeviationPerIteration:         maxDeviationPerIteration,
		}
	}

	sort.SliceStable(statistics, func(i, j int) bool {
		if statistics[i].AverageBestDeviation != statistics[j].AverageBestDeviation {
			return statistics[i].AverageBestDeviation < statistics[j].AverageBestDeviation
		}

		if statistics[i].SuccessRate != statistics[j].SuccessRate {
			return statistics[i].SuccessRate > statistics[j].SuccessRate
		}

		return statistics[i].AverageBestAtIteration < statistics[j].AverageBestAtIteration
	})

	return statistics
}

func RunExperiments(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, heuristicModifiers [][]float64, useThreeOpt bool) []ExperimentResult {
	return RunExperimentsWithRootedHeuristicModifiers(numberOfRuns, parameters, knownOptimal, matrix, heuristicModifiers, nil, useThreeOpt)
}

func RunExperimentsWithRootedHeuristicModifiers(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, heuristicModifiers [][]float64, rootedHeuristicModifiers [][][]float64, useThreeOpt bool) []ExperimentResult {
	results := make([]ExperimentResult, numberOfRuns)

	aco := aco.NewACO(
		parameters.Alpha,
		parameters.Beta,
		parameters.Rho,
		parameters.Iterations,
		knownOptimal,
		matrix,
		heuristicModifiers,
		parameters.HeuristicWeight)
	if rootedHeuristicModifiers != nil {
		aco.SetRootedHeuristicModifiers(rootedHeuristicModifiers)
	}
	aco.SetUseThreeOpt(useThreeOpt)

	for i := 0; i < numberOfRuns; i++ {
		aco.Run()

		results[i] = ExperimentResult{
			BestAtIteration:           aco.BestAtIteration,
			ThreeOptImprovementsCount: aco.ThreeOptImprovementsCount,
			BestTour:                  aco.BestTour,
			DeviationPerIteration:     aco.DeviationPerIteration,
		}
	}

	return results
}
