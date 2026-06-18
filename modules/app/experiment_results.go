package app

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"fmt"
	"math"
	"sort"
	"strconv"
)

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
	return runExperimentsWithRootedHeuristicModifiers(numberOfRuns, parameters, knownOptimal, matrix, heuristicModifiers, nil, useThreeOpt)
}

func runExperimentsWithRootedHeuristicModifiers(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, heuristicModifiers [][]float64, rootedHeuristicModifiers [][][]float64, useThreeOpt bool) []ExperimentResult {
	results := make([]ExperimentResult, numberOfRuns)

	aco := aco.NewACO(
		parameters.alpha,
		parameters.beta,
		parameters.rho,
		parameters.iterations,
		knownOptimal,
		matrix,
		heuristicModifiers,
		parameters.heuristicWeight)
	if rootedHeuristicModifiers != nil {
		aco.SetRootedHeuristicModifiers(rootedHeuristicModifiers)
	}
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
