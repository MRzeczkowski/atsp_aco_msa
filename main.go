package main

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"atsp_aco_msa/modules/algorithms/hungarian"
	"atsp_aco_msa/modules/algorithms/msaSupport"
	"atsp_aco_msa/modules/analysis/cycleCover"
	"atsp_aco_msa/modules/analysis/msaSupportTours"
	"atsp_aco_msa/modules/parsing"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
	"flag"
	"fmt"
	"image/color"
	"math"
	"os"
	"path/filepath"
	"runtime/pprof"
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
)

const (
	msaSupportHighSignalThreshold = 1.0
)

const (
	heuristicMsaSupport           = "msa-support"
	heuristicCycleCover           = "cycle-cover"
	heuristicMsaSupportOverlap    = "msa-support-overlap"
	heuristicMsaSupportDifference = "msa-support-difference"
)

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

func saveBestParametersInfo(fileName string, bestStatistics []ExperimentsDataStatistics) {
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
	reportPath := filepath.Join(resultsDirectoryName, fileName)
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

		alpha, _ := strconv.ParseFloat(record[0], 64)
		beta, _ := strconv.ParseFloat(record[1], 64)
		rho, _ := strconv.ParseFloat(record[2], 64)
		heuristicWeight, _ := strconv.ParseFloat(record[3], 64)
		iterations, _ := strconv.Atoi(record[4])
		minBestAtIteration, _ := strconv.Atoi(record[5])
		averageBestAtIteration, _ := strconv.ParseFloat(record[6], 64)
		maxBestAtIteration, _ := strconv.Atoi(record[7])
		minThreeOptImprovementsCount, _ := strconv.Atoi(record[8])
		averageThreeOptImprovementsCount, _ := strconv.ParseFloat(record[9], 64)
		maxThreeOptImprovementsCount, _ := strconv.Atoi(record[10])
		minBestDeviation, _ := strconv.ParseFloat(record[11], 64)
		averageBestDeviation, _ := strconv.ParseFloat(record[12], 64)
		maxBestDeviation, _ := strconv.ParseFloat(record[13], 64)
		successRate, _ := strconv.ParseFloat(record[14], 64)

		statistic := ExperimentsDataStatistics{
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
		}

		statistics = append(statistics, statistic)
	}

	return statistics, nil
}

func saveStatistics(resultCsvPath string, statistics []ExperimentsDataStatistics) {
	file, _ := os.Create(resultCsvPath)
	defer file.Close()

	writer := csv.NewWriter(file)

	_ = writer.Write(statisticsCsvHeader)

	floatFormat := "%.2f"
	for _, statistic := range statistics {

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

		writer.Write(record)
	}

	writer.Flush()
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

func runExperiments(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, heuristicModifiers [][]float64) []ExperimentResult {
	results := make([]ExperimentResult, numberOfRuns)

	aco := aco.NewACO(
		parameters.alpha,
		parameters.beta,
		parameters.rho,
		parameters.iterations,
		knownOptimal,
		matrix,
		heuristicModifiers)

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

func buildHeuristicModifiers(heuristic string, matrix, msaSupport, cycleCover [][]float64, strength float64) [][]float64 {
	switch heuristic {
	case heuristicMsaSupport:
		return buildMsaSupportHeuristicModifiers(msaSupport, strength)
	case heuristicCycleCover:
		return buildCycleCoverHeuristicModifiers(cycleCover, strength)
	case heuristicMsaSupportOverlap:
		return buildMsaSupportCycleCoverMembershipHeuristicModifiers(msaSupport, cycleCover, strength, true)
	case heuristicMsaSupportDifference:
		return buildMsaSupportCycleCoverMembershipHeuristicModifiers(msaSupport, cycleCover, strength, false)
	default:
		return buildNeutralHeuristicModifiers(len(msaSupport))
	}
}

func buildMsaSupportHeuristicModifiers(msaSupport [][]float64, strength float64) [][]float64 {
	dimension := len(msaSupport)
	modifiers := buildNeutralHeuristicModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxMsaSupportSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			msaSupportSignal := msaSupport[i][j] / maxMsaSupportSelections
			if msaSupportSignal >= msaSupportHighSignalThreshold {
				modifiers[i][j] = 1.0 + msaSupportSignal*strength
			}
		}
	}

	return modifiers
}

func buildMsaSupportCycleCoverMembershipHeuristicModifiers(msaSupport, cycleCover [][]float64, strength float64, requireCycleCoverEdge bool) [][]float64 {
	dimension := len(msaSupport)
	modifiers := buildNeutralHeuristicModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxMsaSupportSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			msaSupportSignal := msaSupport[i][j] / maxMsaSupportSelections
			if msaSupportSignal < msaSupportHighSignalThreshold {
				continue
			}

			if matrixContainsEdge(cycleCover, i, j) != requireCycleCoverEdge {
				continue
			}

			modifiers[i][j] = 1.0 + msaSupportSignal*strength
		}
	}

	return modifiers
}

func matrixContainsEdge(matrix [][]float64, from, to int) bool {
	return from >= 0 && from < len(matrix) && to >= 0 && to < len(matrix[from]) && matrix[from][to] != 0
}

func buildCycleCoverHeuristicModifiers(cycleCover [][]float64, strength float64) [][]float64 {
	dimension := len(cycleCover)
	modifiers := buildNeutralHeuristicModifiers(dimension)
	if dimension == 0 || strength == 0 {
		return modifiers
	}

	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i != j && cycleCover[i][j] != 0 {
				modifiers[i][j] = 1.0 + strength
			}
		}
	}

	return modifiers
}

func buildNeutralHeuristicModifiers(dimension int) [][]float64 {
	modifiers := make([][]float64, dimension)
	for i := range modifiers {
		modifiers[i] = make([]float64, dimension)
		for j := range modifiers[i] {
			modifiers[i][j] = 1.0
		}
	}
	return modifiers
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

	for _, alpha := range utilities.GenerateRange(1.0, 1.0, 0.25) {
		for _, beta := range utilities.GenerateRange(2.0, 2.0, 1.0) {
			for _, rho := range utilities.GenerateRange(0.8, 0.8, 0.1) {
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
	name = strings.TrimSuffix(name, ".atsp")
	resultsDirectoryPath := filepath.Join(resultsDirectoryName, name)
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
	return mode == runModeExperiment || mode == runModeAnalyze || mode == runModeAll
}

func isValidHeuristic(heuristic string) bool {
	return heuristic == heuristicMsaSupport ||
		heuristic == heuristicCycleCover ||
		heuristic == heuristicMsaSupportOverlap ||
		heuristic == heuristicMsaSupportDifference
}

func heuristicUsesCycleCover(heuristic string) bool {
	return heuristic == heuristicCycleCover ||
		heuristic == heuristicMsaSupportOverlap ||
		heuristic == heuristicMsaSupportDifference
}

func heuristicFileSuffix(heuristic string) string {
	switch heuristic {
	case heuristicMsaSupport:
		return ""
	case heuristicCycleCover:
		return "_cycle_cover"
	case heuristicMsaSupportOverlap:
		return "_msa_support_overlap"
	case heuristicMsaSupportDifference:
		return "_msa_support_difference"
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
	if heuristic == heuristicMsaSupport {
		paths, _ := filepath.Glob(filepath.Join(resultsDirectoryName, "*", resultFileName))
		return paths
	}

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

func main() {
	instances := flag.String("instances", instanceSetSmoke, "ATSP instance set to run: smoke, tiny, balanced, large, or all-known")
	mode := flag.String("mode", runModeExperiment, "Run mode: experiment, analyze, or all")
	heuristic := flag.String("heuristic", heuristicMsaSupport, "ACO heuristic modifier to use in experiment mode: msa-support, cycle-cover, msa-support-overlap, or msa-support-difference")
	flag.Parse()

	if !isValidRunMode(*mode) {
		fmt.Printf("Unsupported -mode value %q; use %q, %q, or %q\n", *mode, runModeExperiment, runModeAnalyze, runModeAll)
		return
	}

	if !isValidHeuristic(*heuristic) {
		fmt.Printf("Unsupported -heuristic value %q; use %q, %q, %q, or %q\n", *heuristic, heuristicMsaSupport, heuristicCycleCover, heuristicMsaSupportOverlap, heuristicMsaSupportDifference)
		return
	}

	atspsData, err := loadSelectedAtspData(*instances)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Printf("Selected %d ATSP instance(s) with -instances=%s\n", len(atspsData), *instances)

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

	experimentParameters := generateParameters()
	numberOfExperiments := 30
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
			results := runExperiments(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers)
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
	saveBestParametersInfo(bestParametersReportPathForHeuristic(heuristic), bestStatistics)

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

	_, err = cycleCover.AnalyzeInstances(cycleCover.Config{
		Instances:      cycleCoverConfigs,
		SummaryCsvPath: cycleCoverSummaryPath,
		ReportPath:     cycleCoverReportPath,
		HighThreshold:  1.0,
		Thresholds:     msaSupportTours.DefaultThresholds(),
	})
	if err != nil {
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

	fmt.Printf("MSA support/solution analysis summary saved to %s\n", msaSupportSolutionSummaryPath)
	fmt.Printf("MSA support/solution analysis report saved to %s\n", msaSupportSolutionReportPath)
	fmt.Printf("Cycle-cover analysis summary saved to %s\n", cycleCoverSummaryPath)
	fmt.Printf("Cycle-cover analysis report saved to %s\n", cycleCoverReportPath)
	fmt.Printf("Heuristic boosted-edge summary saved to %s\n", heuristicBoostSummaryPath)
	return nil
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
	heuristics := []string{
		heuristicMsaSupport,
		heuristicCycleCover,
		heuristicMsaSupportOverlap,
		heuristicMsaSupportDifference,
	}

	instances := append([]AtspData(nil), atspsData...)
	sort.SliceStable(instances, func(i, j int) bool {
		return instances[i].name < instances[j].name
	})

	rows := make([]heuristicBoostSummaryRow, 0, len(instances)*len(heuristics))
	for _, atspData := range instances {
		var cycleCover [][]float64
		var cycleCoverErr error
		cycleCoverReady := false

		for _, heuristic := range heuristics {
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
