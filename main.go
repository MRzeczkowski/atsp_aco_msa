package main

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"atsp_aco_msa/modules/algorithms/compositeMsa"
	"atsp_aco_msa/modules/algorithms/edmonds"
	"atsp_aco_msa/modules/algorithms/hungarian"
	"atsp_aco_msa/modules/analysis/cmsaTours"
	"atsp_aco_msa/modules/analysis/cycleCover"
	"atsp_aco_msa/modules/models"
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
	alpha, beta, rho, pCmsa float64
	iterations              int
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
	instanceSetBalanced = "balanced"
	instanceSetAllKnown = "all-known"
)

const (
	runModeExperiment = "experiment"
	runModeAnalyze    = "analyze"
	runModeAll        = "all"
)

const (
	cmsaHighSignalThreshold             = 1.0
	cmsaConnectorStrengthFactor         = 0.5
	arborescenceConnectorStrengthFactor = 0.5
)

const (
	heuristicCmsa                    = "cmsa"
	heuristicCmsaa                   = "cmsaa"
	heuristicCmsaAgreement           = "cmsa-agreement"
	heuristicCycleCover              = "cycle-cover"
	heuristicBoth                    = "both"
	heuristicArborescences           = "cycle-cover-arborescences"
	heuristicContractedArborescences = "cycle-cover-contracted-arborescences"
	heuristicSpliceArborescences     = "cycle-cover-splice-arborescences"
	heuristicCmsaOverlap             = "cmsa-overlap"
	heuristicCmsaDifference          = "cmsa-difference"
)

var smokeInstanceFiles = []string{"ftv170.atsp"}

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

var statisticsCsvHeader = []string{
	"Alpha",
	"Beta",
	"Rho",
	"pCmsa",
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
	pCmsaCounts := make(map[float64]int)

	// Count the occurrences
	for params := range uniqueParameters {
		alphaCounts[params.alpha]++
		betaCounts[params.beta]++
		rhoCounts[params.rho]++
		pCmsaCounts[params.pCmsa]++
	}

	// Markdown content
	markdown := "# Best Parameters Report\n\n"

	// Unique combinations count
	markdown += fmt.Sprintf("Found **%d** best unique parameter combinations.\n\n", len(uniqueParameters))

	// Best parameters
	markdown += "## Best Parameters\n\n"
	markdown += "| Alpha | Beta | Rho | pCmsa | Times used |\n"
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
		return left.pCmsa < right.pCmsa
	})

	for _, parameters := range sortedParameters {
		timesUsed := uniqueParameters[parameters]
		markdown += fmt.Sprintf("| %.2f | %.2f | %.2f | %.2f | %d |\n",
			parameters.alpha, parameters.beta, parameters.rho, parameters.pCmsa, timesUsed)
	}
	markdown += "\n"

	// Parameter occurrences
	markdown += "## Parameter Values Occurrences\n\n"
	markdown += generateMarkdownCounts("Alpha", alphaCounts)
	markdown += generateMarkdownCounts("Beta", betaCounts)
	markdown += generateMarkdownCounts("Rho", rhoCounts)
	markdown += generateMarkdownCounts("pCmsa", pCmsaCounts)

	// Parameter ranges
	markdown += "## Parameter Ranges\n\n"
	minAlpha, maxAlpha := findMinMax(alphaCounts)
	minBeta, maxBeta := findMinMax(betaCounts)
	minRho, maxRho := findMinMax(rhoCounts)
	minPCmsa, maxPCmsa := findMinMax(pCmsaCounts)

	markdown += fmt.Sprintf("- **Alpha**: %.2f - %.2f\n", minAlpha, maxAlpha)
	markdown += fmt.Sprintf("- **Beta**: %.2f - %.2f\n", minBeta, maxBeta)
	markdown += fmt.Sprintf("- **Rho**: %.2f - %.2f\n", minRho, maxRho)
	markdown += fmt.Sprintf("- **pCmsa**: %.2f - %.2f\n", minPCmsa, maxPCmsa)

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
		pCmsa, _ := strconv.ParseFloat(record[3], 64)
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
				alpha:      alpha,
				beta:       beta,
				rho:        rho,
				pCmsa:      pCmsa,
				iterations: iterations,
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
			fmt.Sprintf(floatFormat, statistic.pCmsa),
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

func buildHeuristicModifiers(heuristic string, matrix, cmsa, cycleCover [][]float64, strength float64) [][]float64 {
	switch heuristic {
	case heuristicCmsa:
		return buildCmsaHeuristicModifiers(cmsa, strength)
	case heuristicCmsaa:
		return buildCmsaHeuristicModifiers(cmsa, strength)
	case heuristicCmsaAgreement:
		return buildCmsaHeuristicModifiers(cmsa, strength)
	case heuristicCycleCover:
		return buildCycleCoverHeuristicModifiers(cycleCover, strength)
	case heuristicBoth:
		return buildCmsaCycleCoverHeuristicModifiers(cmsa, cycleCover, matrix, strength)
	case heuristicArborescences:
		return buildCycleCoverArborescenceHeuristicModifiers(matrix, cycleCover, strength)
	case heuristicContractedArborescences:
		return buildCycleCoverContractedArborescenceHeuristicModifiers(matrix, cycleCover, strength)
	case heuristicSpliceArborescences:
		return buildCycleCoverSpliceArborescenceHeuristicModifiers(matrix, cycleCover, strength)
	case heuristicCmsaOverlap:
		return buildCmsaCycleCoverMembershipHeuristicModifiers(cmsa, cycleCover, strength, true)
	case heuristicCmsaDifference:
		return buildCmsaCycleCoverMembershipHeuristicModifiers(cmsa, cycleCover, strength, false)
	default:
		return buildNeutralHeuristicModifiers(len(cmsa))
	}
}

func buildCmsaHeuristicModifiers(cmsa [][]float64, strength float64) [][]float64 {
	dimension := len(cmsa)
	modifiers := buildNeutralHeuristicModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxCmsaSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			cmsaSignal := cmsa[i][j] / maxCmsaSelections
			if cmsaSignal >= cmsaHighSignalThreshold {
				modifiers[i][j] = 1.0 + cmsaSignal*strength
			}
		}
	}

	return modifiers
}

func buildCmsaCycleCoverMembershipHeuristicModifiers(cmsa, cycleCover [][]float64, strength float64, requireCycleCoverEdge bool) [][]float64 {
	dimension := len(cmsa)
	modifiers := buildNeutralHeuristicModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxCmsaSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			cmsaSignal := cmsa[i][j] / maxCmsaSelections
			if cmsaSignal < cmsaHighSignalThreshold {
				continue
			}

			if matrixContainsEdge(cycleCover, i, j) != requireCycleCoverEdge {
				continue
			}

			modifiers[i][j] = 1.0 + cmsaSignal*strength
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

func buildCycleCoverArborescenceHeuristicModifiers(matrix, cycleCover [][]float64, strength float64) [][]float64 {
	dimension := len(matrix)
	if dimension <= 1 || strength == 0 {
		return buildNeutralHeuristicModifiers(dimension)
	}

	componentIds := buildCycleCoverComponentIds(cycleCover)
	connectors := selectInterCycleArborescenceConnectors(matrix, componentIds)
	return buildCycleCoverConnectorHeuristicModifiers(dimension, cycleCover, connectors, nil, strength)
}

func buildCycleCoverContractedArborescenceHeuristicModifiers(matrix, cycleCover [][]float64, strength float64) [][]float64 {
	dimension := len(matrix)
	if dimension <= 1 || strength == 0 {
		return buildNeutralHeuristicModifiers(dimension)
	}

	componentIds := buildCycleCoverComponentIds(cycleCover)
	connectors := selectInterCycleContractedArborescenceConnectors(matrix, componentIds)
	return buildCycleCoverConnectorHeuristicModifiers(dimension, cycleCover, connectors, nil, strength)
}

func buildCycleCoverSpliceArborescenceHeuristicModifiers(matrix, cycleCover [][]float64, strength float64) [][]float64 {
	dimension := len(matrix)
	if dimension <= 1 || strength == 0 {
		return buildNeutralHeuristicModifiers(dimension)
	}

	componentIds := buildCycleCoverComponentIds(cycleCover)
	patches := selectInterCycleArborescencePatches(matrix, cycleCover, componentIds)
	return buildCycleCoverPatchHeuristicModifiers(dimension, cycleCover, patches, strength)
}

func buildCycleCoverConnectorHeuristicModifiers(
	dimension int,
	cycleCover [][]float64,
	connectors map[models.Edge]bool,
	splicedCycleCoverEdges map[models.Edge]bool,
	strength float64,
) [][]float64 {
	modifiers := buildNeutralHeuristicModifiers(dimension)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			if matrixContainsEdge(cycleCover, i, j) {
				if !splicedCycleCoverEdges[models.Edge{From: i, To: j}] {
					modifiers[i][j] = 1.0 + strength
				}
				continue
			}

			if connectors[models.Edge{From: i, To: j}] {
				modifiers[i][j] = 1.0 + strength*arborescenceConnectorStrengthFactor
			}
		}
	}

	return modifiers
}

func selectInterCycleContractedArborescenceConnectors(matrix [][]float64, componentIds []int) map[models.Edge]bool {
	connectors := make(map[models.Edge]bool)
	dimension := len(matrix)
	componentCount := cycleCoverComponentCount(componentIds)
	if dimension <= 1 || len(componentIds) != dimension || componentCount <= 1 {
		return connectors
	}

	componentGraph, originalEdgesByComponentEdge := buildCycleCoverComponentGraph(matrix, componentIds)
	componentEdges := selectContractedArborescenceEdges(componentGraph, cycleCoverComponentRoot(componentIds))
	for _, componentEdge := range sortedModelEdges(componentEdges) {
		if originalEdge, ok := originalEdgesByComponentEdge[componentEdge]; ok {
			connectors[originalEdge] = true
		}
	}

	return connectors
}

type cycleCoverPatch struct {
	componentEdge    models.Edge
	insertedForward  models.Edge
	insertedBackward models.Edge
	removedFrom      models.Edge
	removedTo        models.Edge
	delta            float64
	valid            bool
}

func selectInterCycleArborescencePatches(matrix, cycleCover [][]float64, componentIds []int) []cycleCoverPatch {
	dimension := len(matrix)
	componentCount := cycleCoverComponentCount(componentIds)
	if dimension <= 1 || len(cycleCover) != dimension || len(componentIds) != dimension || componentCount <= 1 {
		return nil
	}

	componentEdges := selectInterCycleContractedArborescenceComponentEdges(matrix, componentIds)
	patches := make([]cycleCoverPatch, 0, len(componentEdges))
	for _, componentEdge := range sortedModelEdges(componentEdges) {
		patch, ok := findBestCycleCoverPatch(matrix, cycleCover, componentIds, componentEdge)
		if ok {
			patches = append(patches, patch)
		}
	}

	return patches
}

func selectInterCycleContractedArborescenceComponentEdges(matrix [][]float64, componentIds []int) map[models.Edge]bool {
	componentEdges := make(map[models.Edge]bool)
	dimension := len(matrix)
	componentCount := cycleCoverComponentCount(componentIds)
	if dimension <= 1 || len(componentIds) != dimension || componentCount <= 1 {
		return componentEdges
	}

	componentGraph, _ := buildCycleCoverComponentGraph(matrix, componentIds)
	return selectContractedArborescenceEdges(componentGraph, cycleCoverComponentRoot(componentIds))
}

func selectContractedArborescenceEdges(componentGraph [][]float64, root int) map[models.Edge]bool {
	componentEdges := make(map[models.Edge]bool)
	if len(componentGraph) <= 1 {
		return componentEdges
	}

	vertices, edges, weights := models.ConvertToEdges(componentGraph)
	if root < 0 {
		root = 0
	}

	msa := edmonds.FindMSA(root, vertices, edges, weights)
	msaa := edmonds.FindMSAA(root, vertices, edges, weights)

	for _, edge := range msa {
		componentEdges[edge] = true
	}
	for _, edge := range msaa {
		componentEdges[edge] = true
	}

	return componentEdges
}

func findBestCycleCoverPatch(matrix, cycleCover [][]float64, componentIds []int, componentEdge models.Edge) (cycleCoverPatch, bool) {
	if componentEdge.From == componentEdge.To || len(cycleCover) != len(matrix) || len(componentIds) != len(matrix) {
		return cycleCoverPatch{}, false
	}

	successors := buildCycleCoverSuccessors(cycleCover)
	bestPatch := cycleCoverPatch{}
	for from := 0; from < len(matrix); from++ {
		if from >= len(componentIds) || componentIds[from] != componentEdge.From {
			continue
		}

		fromSuccessor := successors[from]
		if fromSuccessor < 0 {
			continue
		}

		for to := 0; to < len(matrix); to++ {
			if to >= len(componentIds) || componentIds[to] != componentEdge.To {
				continue
			}

			toSuccessor := successors[to]
			if toSuccessor < 0 {
				continue
			}

			insertedForward := models.Edge{From: from, To: toSuccessor}
			insertedBackward := models.Edge{From: to, To: fromSuccessor}
			removedFrom := models.Edge{From: from, To: fromSuccessor}
			removedTo := models.Edge{From: to, To: toSuccessor}
			delta := edgeDistance(matrix, insertedForward.From, insertedForward.To) +
				edgeDistance(matrix, insertedBackward.From, insertedBackward.To) -
				edgeDistance(matrix, removedFrom.From, removedFrom.To) -
				edgeDistance(matrix, removedTo.From, removedTo.To)
			if math.IsInf(delta, 0) || math.IsNaN(delta) {
				continue
			}

			patch := cycleCoverPatch{
				componentEdge:    componentEdge,
				insertedForward:  insertedForward,
				insertedBackward: insertedBackward,
				removedFrom:      removedFrom,
				removedTo:        removedTo,
				delta:            delta,
				valid:            true,
			}
			if isBetterCycleCoverPatch(patch, bestPatch) {
				bestPatch = patch
			}
		}
	}

	return bestPatch, bestPatch.valid
}

func cycleCoverComponentRoot(componentIds []int) int {
	if len(componentIds) == 0 || componentIds[0] < 0 {
		return 0
	}
	return componentIds[0]
}

func buildCycleCoverPatchHeuristicModifiers(dimension int, cycleCover [][]float64, patches []cycleCoverPatch, strength float64) [][]float64 {
	insertedEdges := make(map[models.Edge]bool)
	removedCycleCoverEdges := make(map[models.Edge]bool)
	for _, patch := range patches {
		insertedEdges[patch.insertedForward] = true
		insertedEdges[patch.insertedBackward] = true
		removedCycleCoverEdges[patch.removedFrom] = true
		removedCycleCoverEdges[patch.removedTo] = true
	}

	return buildCycleCoverConnectorHeuristicModifiers(dimension, cycleCover, insertedEdges, removedCycleCoverEdges, strength)
}

func buildCycleCoverSuccessors(cycleCover [][]float64) []int {
	successors := make([]int, len(cycleCover))
	for i := range successors {
		successors[i] = -1
	}

	for from := range cycleCover {
		successors[from] = cycleCoverSuccessor(cycleCover, from)
	}

	return successors
}

func isBetterCycleCoverPatch(candidate, current cycleCoverPatch) bool {
	if !current.valid {
		return true
	}
	if candidate.delta != current.delta {
		return candidate.delta < current.delta
	}
	if candidate.componentEdge != current.componentEdge {
		return edgeLess(candidate.componentEdge, current.componentEdge)
	}
	if candidate.insertedForward != current.insertedForward {
		return edgeLess(candidate.insertedForward, current.insertedForward)
	}
	if candidate.insertedBackward != current.insertedBackward {
		return edgeLess(candidate.insertedBackward, current.insertedBackward)
	}
	if candidate.removedFrom != current.removedFrom {
		return edgeLess(candidate.removedFrom, current.removedFrom)
	}
	return edgeLess(candidate.removedTo, current.removedTo)
}

func sortedModelEdges(edgeSet map[models.Edge]bool) []models.Edge {
	edges := make([]models.Edge, 0, len(edgeSet))
	for edge := range edgeSet {
		edges = append(edges, edge)
	}
	sort.Slice(edges, func(i, j int) bool {
		return edgeLess(edges[i], edges[j])
	})
	return edges
}

func buildCycleCoverComponentGraph(matrix [][]float64, componentIds []int) ([][]float64, map[models.Edge]models.Edge) {
	componentCount := cycleCoverComponentCount(componentIds)
	componentGraph := make([][]float64, componentCount)
	for i := range componentGraph {
		componentGraph[i] = make([]float64, componentCount)
		for j := range componentGraph[i] {
			if i != j {
				componentGraph[i][j] = math.Inf(1)
			}
		}
	}

	originalEdgesByComponentEdge := make(map[models.Edge]models.Edge)
	for from := 0; from < len(matrix); from++ {
		if from >= len(componentIds) || componentIds[from] < 0 {
			continue
		}

		fromComponent := componentIds[from]
		for to := 0; to < len(matrix[from]); to++ {
			if to >= len(componentIds) || componentIds[to] < 0 {
				continue
			}

			toComponent := componentIds[to]
			if fromComponent == toComponent {
				continue
			}

			componentEdge := models.Edge{From: fromComponent, To: toComponent}
			originalEdge := models.Edge{From: from, To: to}
			currentEdge, exists := originalEdgesByComponentEdge[componentEdge]
			currentDistance := componentGraph[fromComponent][toComponent]
			distance := matrix[from][to]
			if !exists || distance < currentDistance || (distance == currentDistance && edgeLess(originalEdge, currentEdge)) {
				componentGraph[fromComponent][toComponent] = distance
				originalEdgesByComponentEdge[componentEdge] = originalEdge
			}
		}
	}

	return componentGraph, originalEdgesByComponentEdge
}

func edgeLess(left, right models.Edge) bool {
	if left.From != right.From {
		return left.From < right.From
	}
	return left.To < right.To
}

func selectInterCycleArborescenceConnectors(matrix [][]float64, componentIds []int) map[models.Edge]bool {
	connectors := make(map[models.Edge]bool)
	dimension := len(matrix)
	if dimension <= 1 || len(componentIds) != dimension || cycleCoverComponentCount(componentIds) <= 1 {
		return connectors
	}

	vertices, edges, weights := models.ConvertToEdges(matrix)
	root := 0
	msa := edmonds.FindMSA(root, vertices, edges, weights)
	msaa := edmonds.FindMSAA(root, vertices, edges, weights)

	for _, edge := range msa {
		if edgeConnectsDifferentComponents(edge, componentIds) {
			connectors[edge] = true
		}
	}
	for _, edge := range msaa {
		if edgeConnectsDifferentComponents(edge, componentIds) {
			connectors[edge] = true
		}
	}

	return connectors
}

func edgeConnectsDifferentComponents(edge models.Edge, componentIds []int) bool {
	if edge.From < 0 || edge.From >= len(componentIds) || edge.To < 0 || edge.To >= len(componentIds) {
		return false
	}

	fromComponent := componentIds[edge.From]
	toComponent := componentIds[edge.To]
	return fromComponent >= 0 && toComponent >= 0 && fromComponent != toComponent
}

func buildCycleCoverComponentIds(cycleCover [][]float64) []int {
	dimension := len(cycleCover)
	componentIds := make([]int, dimension)
	for i := range componentIds {
		componentIds[i] = -1
	}

	componentId := 0
	for start := 0; start < dimension; start++ {
		if componentIds[start] != -1 {
			continue
		}

		current := start
		for current >= 0 && current < dimension && componentIds[current] == -1 {
			componentIds[current] = componentId
			current = cycleCoverSuccessor(cycleCover, current)
		}

		componentId++
	}

	return componentIds
}

func cycleCoverSuccessor(cycleCover [][]float64, vertex int) int {
	for next, value := range cycleCover[vertex] {
		if value != 0 {
			return next
		}
	}

	return -1
}

type cmsaConnectorEdge struct {
	from, to int
}

type cmsaConnectorCandidate struct {
	edge       cmsaConnectorEdge
	cmsaWeight float64
	distance   float64
	valid      bool
}

func selectCmsaCycleCoverConnectors(cmsa, distances [][]float64, cycleCoverComponentIds []int) map[cmsaConnectorEdge]bool {
	connectors := make(map[cmsaConnectorEdge]bool)
	dimension := len(cmsa)
	componentCount := cycleCoverComponentCount(cycleCoverComponentIds)
	if dimension <= 1 || componentCount <= 1 {
		return connectors
	}

	bestOutgoing := make([]cmsaConnectorCandidate, componentCount)
	bestIncoming := make([]cmsaConnectorCandidate, componentCount)
	maxCmsaSelections := float64(dimension - 1)

	for i := 0; i < dimension; i++ {
		fromComponent := cycleCoverComponentIds[i]
		if fromComponent < 0 || fromComponent >= componentCount {
			continue
		}

		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			toComponent := cycleCoverComponentIds[j]
			if toComponent < 0 || toComponent >= componentCount || fromComponent == toComponent {
				continue
			}

			cmsaSignal := cmsa[i][j] / maxCmsaSelections
			if cmsaSignal < cmsaHighSignalThreshold {
				continue
			}

			candidate := cmsaConnectorCandidate{
				edge:       cmsaConnectorEdge{from: i, to: j},
				cmsaWeight: cmsa[i][j],
				distance:   edgeDistance(distances, i, j),
				valid:      true,
			}
			if isBetterCmsaConnectorCandidate(candidate, bestOutgoing[fromComponent]) {
				bestOutgoing[fromComponent] = candidate
			}
			if isBetterCmsaConnectorCandidate(candidate, bestIncoming[toComponent]) {
				bestIncoming[toComponent] = candidate
			}
		}
	}

	for _, candidate := range bestOutgoing {
		if candidate.valid {
			connectors[candidate.edge] = true
		}
	}
	for _, candidate := range bestIncoming {
		if candidate.valid {
			connectors[candidate.edge] = true
		}
	}

	return connectors
}

func cycleCoverComponentCount(componentIds []int) int {
	componentCount := 0
	for _, componentId := range componentIds {
		if componentId >= componentCount {
			componentCount = componentId + 1
		}
	}
	return componentCount
}

func edgeDistance(distances [][]float64, from, to int) float64 {
	if from < len(distances) && to < len(distances[from]) {
		return distances[from][to]
	}
	return math.Inf(1)
}

func isBetterCmsaConnectorCandidate(candidate, current cmsaConnectorCandidate) bool {
	if !current.valid {
		return true
	}
	if candidate.cmsaWeight != current.cmsaWeight {
		return candidate.cmsaWeight > current.cmsaWeight
	}
	if candidate.distance != current.distance {
		return candidate.distance < current.distance
	}
	if candidate.edge.from != current.edge.from {
		return candidate.edge.from < current.edge.from
	}
	return candidate.edge.to < current.edge.to
}

func buildCmsaCycleCoverHeuristicModifiers(cmsa, cycleCover, distances [][]float64, strength float64) [][]float64 {
	dimension := len(cmsa)
	modifiers := buildNeutralHeuristicModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	cycleCoverComponentIds := buildCycleCoverComponentIds(cycleCover)
	cmsaConnectors := selectCmsaCycleCoverConnectors(cmsa, distances, cycleCoverComponentIds)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			cmsaConnectorSignal := 0.0
			if cmsaConnectors[cmsaConnectorEdge{from: i, to: j}] {
				cmsaConnectorSignal = 1.0
			}

			cycleCoverSignal := 0.0
			if cycleCover[i][j] != 0 {
				cycleCoverSignal = 1.0
			}

			combinedSignal := cycleCoverSignal
			if cycleCoverSignal == 0 {
				combinedSignal = cmsaConnectorSignal * cmsaConnectorStrengthFactor
			}
			modifiers[i][j] = 1.0 + combinedSignal*strength
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
				for _, pCmsa := range []float64{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0} {

					parameters = append(parameters,
						ExperimentParameters{
							alpha, beta, rho, pCmsa, 0,
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

	cmsaDirectoryPath,
	cmsaaDirectoryPath,

	cmsaHeatmapPlotPath, cmsaHistogramPlotPath,
	cmsaaHeatmapPlotPath, cmsaaHistogramPlotPath,

	resultFilePath,
	resultPlotFilePrefix,

	optimalUniqueToursCsvPath,
	toursHeatmapPlotPath,
	toursHistogramPlotPath,
	cmsaToursOverlapHeatmapPlotPath,
	cmsaSolutionAnalysisCsvPath,
	cmsaSolutionThresholdsCsvPath,
	cycleCoverEdgesCsvPath,
	cycleCoverAnalysisCsvPath,
	cycleCoverThresholdsCsvPath,
	cycleCoverCmsaOverlapCsvPath string
}

func makeAtspData(name string, matrix [][]float64, knownOptimal float64) AtspData {
	name = strings.TrimSuffix(name, ".atsp")
	resultsDirectoryPath := filepath.Join(resultsDirectoryName, name)
	cmsaDirectoryPath := filepath.Join(resultsDirectoryPath, "cmsa")
	cmsaaDirectoryPath := filepath.Join(resultsDirectoryPath, "cmsaa")
	plotsDirectoryPath := filepath.Join(resultsDirectoryPath, "plots")

	cmsaHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "cmsa_heatmap.png")
	cmsaHistogramPlotPath := filepath.Join(plotsDirectoryPath, "cmsa_histogram.png")
	cmsaaHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "cmsaa_heatmap.png")
	cmsaaHistogramPlotPath := filepath.Join(plotsDirectoryPath, "cmsaa_histogram.png")

	resultFilePath := filepath.Join(resultsDirectoryPath, resultFileName)
	resultPlotFilePrefix := filepath.Join(plotsDirectoryPath, "best_result")

	optimalUniqueToursCsvPath := filepath.Join(resultsDirectoryPath, "solutions.csv")
	toursHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "tours_heatmap.png")
	toursHistogramPlotPath := filepath.Join(plotsDirectoryPath, "tours_histogram.png")
	cmsaToursOverlapHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "cmsa_tours_overlap_heatmap.png")
	cmsaSolutionAnalysisCsvPath := filepath.Join(resultsDirectoryPath, "cmsa_solution_analysis.csv")
	cmsaSolutionThresholdsCsvPath := filepath.Join(resultsDirectoryPath, "cmsa_solution_thresholds.csv")
	cycleCoverEdgesCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_edges.csv")
	cycleCoverAnalysisCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_analysis.csv")
	cycleCoverThresholdsCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_thresholds.csv")
	cycleCoverCmsaOverlapCsvPath := filepath.Join(resultsDirectoryPath, "cycle_cover_cmsa_overlap.csv")

	return AtspData{
		name,
		matrix,
		knownOptimal,

		cmsaDirectoryPath,
		cmsaaDirectoryPath,

		cmsaHeatmapPlotPath, cmsaHistogramPlotPath,
		cmsaaHeatmapPlotPath, cmsaaHistogramPlotPath,

		resultFilePath,
		resultPlotFilePrefix,

		optimalUniqueToursCsvPath,
		toursHeatmapPlotPath,
		toursHistogramPlotPath,
		cmsaToursOverlapHeatmapPlotPath,
		cmsaSolutionAnalysisCsvPath,
		cmsaSolutionThresholdsCsvPath,
		cycleCoverEdgesCsvPath,
		cycleCoverAnalysisCsvPath,
		cycleCoverThresholdsCsvPath,
		cycleCoverCmsaOverlapCsvPath,
	}
}

func selectAtspFiles(atspFilePaths []string, instanceSet string) ([]string, error) {
	switch instanceSet {
	case instanceSetSmoke:
		return selectConfiguredAtspFiles(atspFilePaths, smokeInstanceFiles)
	case instanceSetBalanced:
		return selectConfiguredAtspFiles(atspFilePaths, balancedInstanceFiles)
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
		return nil, fmt.Errorf("unsupported -instances value %q; use %q, %q, or %q", instanceSet, instanceSetSmoke, instanceSetBalanced, instanceSetAllKnown)
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
	return heuristic == heuristicCmsa ||
		heuristic == heuristicCmsaa ||
		heuristic == heuristicCmsaAgreement ||
		heuristic == heuristicCycleCover ||
		heuristic == heuristicBoth ||
		heuristic == heuristicArborescences ||
		heuristic == heuristicContractedArborescences ||
		heuristic == heuristicSpliceArborescences ||
		heuristic == heuristicCmsaOverlap ||
		heuristic == heuristicCmsaDifference
}

func heuristicUsesCycleCover(heuristic string) bool {
	return heuristic == heuristicCycleCover ||
		heuristic == heuristicBoth ||
		heuristic == heuristicArborescences ||
		heuristic == heuristicContractedArborescences ||
		heuristic == heuristicSpliceArborescences ||
		heuristic == heuristicCmsaOverlap ||
		heuristic == heuristicCmsaDifference
}

func heuristicUsesCmsaa(heuristic string) bool {
	return heuristic == heuristicCmsaa ||
		heuristic == heuristicCmsaAgreement
}

func heuristicFileSuffix(heuristic string) string {
	switch heuristic {
	case heuristicCmsa:
		return ""
	case heuristicCmsaa:
		return "_cmsaa"
	case heuristicCmsaAgreement:
		return "_cmsa_agreement"
	case heuristicCycleCover:
		return "_cycle_cover"
	case heuristicBoth:
		return "_both"
	case heuristicArborescences:
		return "_cycle_cover_arborescences"
	case heuristicContractedArborescences:
		return "_cycle_cover_contracted_arborescences"
	case heuristicSpliceArborescences:
		return "_cycle_cover_splice_arborescences"
	case heuristicCmsaOverlap:
		return "_cmsa_overlap"
	case heuristicCmsaDifference:
		return "_cmsa_difference"
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
	if heuristic == heuristicCmsa {
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
	instances := flag.String("instances", instanceSetSmoke, "ATSP instance set to run: smoke, balanced, or all-known")
	mode := flag.String("mode", runModeExperiment, "Run mode: experiment, analyze, or all")
	heuristic := flag.String("heuristic", heuristicCmsa, "ACO heuristic modifier to use in experiment mode: cmsa, cmsaa, cmsa-agreement, cycle-cover, both, cycle-cover-arborescences, cycle-cover-contracted-arborescences, cycle-cover-splice-arborescences, cmsa-overlap, or cmsa-difference")
	flag.Parse()

	if !isValidRunMode(*mode) {
		fmt.Printf("Unsupported -mode value %q; use %q, %q, or %q\n", *mode, runModeExperiment, runModeAnalyze, runModeAll)
		return
	}

	if !isValidHeuristic(*heuristic) {
		fmt.Printf("Unsupported -heuristic value %q; use %q, %q, %q, %q, %q, %q, %q, %q, %q, or %q\n", *heuristic, heuristicCmsa, heuristicCmsaa, heuristicCmsaAgreement, heuristicCycleCover, heuristicBoth, heuristicArborescences, heuristicContractedArborescences, heuristicSpliceArborescences, heuristicCmsaOverlap, heuristicCmsaDifference)
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
	if err := ensureCmsaArtifacts(atspsData); err != nil {
		return err
	}
	if heuristicUsesCmsaa(heuristic) {
		if err := ensureCmsaaArtifacts(atspsData); err != nil {
			return err
		}
	}

	experimentParameters := generateParameters()
	numberOfExperiments := 10
	for _, atspData := range atspsData {
		matrix := atspData.matrix
		knownOptimal := atspData.knownOptimal
		dimension := len(matrix)
		instanceStart := time.Now()

		heuristicMatrix, err := readCompositeMatrixForHeuristic(atspData, heuristic)
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
			heuristicModifiers := buildHeuristicModifiers(heuristic, matrix, heuristicMatrix, cycleCover, parameters.pCmsa)
			results := runExperiments(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers)
			data := ExperimentsData{parameters, results}

			experimentData = append(experimentData, data)

			parameterStatistics := calculateStatistics([]ExperimentsData{data})
			if len(parameterStatistics) != 0 {
				statistic := parameterStatistics[0]
				fmt.Printf("\tpCmsa=%.2f iterations=%d runs=%d elapsed=%s min deviation=%.2f avg deviation=%.2f\n",
					parameters.pCmsa,
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

		uniqueOptimalTours, err := cmsaTours.ReadOptimalTours(atspData.optimalUniqueToursCsvPath)
		if err != nil {
			return err
		}

		for _, data := range experimentData {
			for _, result := range data.results {
				if result.deviationPerIteration[result.bestAtIteration] == 0.0 {
					cmsaTours.AddUniqueTour(uniqueOptimalTours, result.bestTour)
				}
			}
		}

		if err := cmsaTours.SaveOptimalToursStatistics(atspData.optimalUniqueToursCsvPath, atspData.cmsaDirectoryPath, uniqueOptimalTours); err != nil {
			return err
		}

		fmt.Printf("Finished %s in %s\n", atspData.name, time.Since(instanceStart).Round(time.Millisecond))
	}

	bestStatistics := getBestStatisticsFromFiles(resultFilePathsForHeuristic(atspsData, heuristic))
	saveBestParametersInfo(bestParametersReportPathForHeuristic(heuristic), bestStatistics)

	return nil
}

func readCompositeMatrixForHeuristic(atspData AtspData, heuristic string) ([][]float64, error) {
	if heuristic == heuristicCmsaa {
		return compositeMsa.ReadAnti(atspData.cmsaaDirectoryPath)
	}
	if heuristic == heuristicCmsaAgreement {
		cmsa, err := compositeMsa.Read(atspData.cmsaDirectoryPath)
		if err != nil {
			return nil, err
		}

		cmsaa, err := compositeMsa.ReadAnti(atspData.cmsaaDirectoryPath)
		if err != nil {
			return nil, err
		}

		return buildCmsaAgreementMatrix(cmsa, cmsaa)
	}

	return compositeMsa.Read(atspData.cmsaDirectoryPath)
}

func buildCmsaAgreementMatrix(cmsa, cmsaa [][]float64) ([][]float64, error) {
	if len(cmsa) != len(cmsaa) {
		return nil, fmt.Errorf("CMSA and CMSAA dimensions differ: %d != %d", len(cmsa), len(cmsaa))
	}

	dimension := len(cmsa)
	agreement := make([][]float64, dimension)
	for i := 0; i < dimension; i++ {
		if len(cmsa[i]) != dimension {
			return nil, fmt.Errorf("CMSA row %d has length %d, expected %d", i, len(cmsa[i]), dimension)
		}
		if len(cmsaa[i]) != dimension {
			return nil, fmt.Errorf("CMSAA row %d has length %d, expected %d", i, len(cmsaa[i]), dimension)
		}

		agreement[i] = make([]float64, dimension)
		for j := 0; j < dimension; j++ {
			agreement[i][j] = math.Min(cmsa[i][j], cmsaa[i][j])
		}
	}

	return agreement, nil
}

func ensureCmsaArtifacts(atspsData []AtspData) error {
	for _, atspData := range atspsData {
		name := atspData.name
		matrix := atspData.matrix
		cmsaDirectoryPath := atspData.cmsaDirectoryPath

		cmsa, err := compositeMsa.Read(atspData.cmsaDirectoryPath)

		if err != nil {
			start := time.Now()
			cmsa, err = compositeMsa.Create(matrix, cmsaDirectoryPath)
			elapsed := time.Since(start)

			fmt.Printf("\tCreating %s took: %d ms\n", cmsaDirectoryPath, elapsed.Milliseconds())

			if err != nil {
				return fmt.Errorf("error saving CMSA: %w", err)
			}
		}

		cmsaHeatmapPlotTitle := name + " CMSA heatmap"

		err = utilities.SaveHeatmapFromMatrix(cmsa, cmsaHeatmapPlotTitle, atspData.cmsaHeatmapPlotPath)
		if err != nil {
			return err
		}

		dataForHistogram := filterZeroes(flattenMatrix(cmsa))
		cmsaHistogramPlotTitle := name + " CMSA histogram"

		dimension := len(matrix)
		err = utilities.SaveHistogramFromData(dataForHistogram, dimension-1, cmsaHistogramPlotTitle, atspData.cmsaHistogramPlotPath)
		if err != nil {
			return err
		}
	}

	return nil
}

func ensureCmsaaArtifacts(atspsData []AtspData) error {
	for _, atspData := range atspsData {
		name := atspData.name
		matrix := atspData.matrix
		cmsaaDirectoryPath := atspData.cmsaaDirectoryPath

		cmsaa, err := compositeMsa.ReadAnti(cmsaaDirectoryPath)
		if err != nil {
			start := time.Now()
			cmsaa, err = compositeMsa.CreateAnti(matrix, cmsaaDirectoryPath)
			elapsed := time.Since(start)

			fmt.Printf("\tCreating %s took: %d ms\n", cmsaaDirectoryPath, elapsed.Milliseconds())

			if err != nil {
				return fmt.Errorf("error saving CMSAA: %w", err)
			}
		}

		cmsaaHeatmapPlotTitle := name + " CMSAA heatmap"
		if err := utilities.SaveHeatmapFromMatrix(cmsaa, cmsaaHeatmapPlotTitle, atspData.cmsaaHeatmapPlotPath); err != nil {
			return err
		}

		dataForHistogram := filterZeroes(flattenMatrix(cmsaa))
		cmsaaHistogramPlotTitle := name + " CMSAA histogram"
		dimension := len(matrix)
		if err := utilities.SaveHistogramFromData(dataForHistogram, dimension-1, cmsaaHistogramPlotTitle, atspData.cmsaaHistogramPlotPath); err != nil {
			return err
		}
	}

	return nil
}

func runAnalysisMode(atspsData []AtspData) error {
	cmsaTourConfigs := make([]cmsaTours.InstanceConfig, 0, len(atspsData))
	cycleCoverConfigs := make([]cycleCover.InstanceConfig, 0, len(atspsData))
	for _, atspData := range atspsData {
		cmsaTourConfigs = append(cmsaTourConfigs, cmsaTours.InstanceConfig{
			Name:                        atspData.name,
			Dimension:                   len(atspData.matrix),
			CmsaDirectoryPath:           atspData.cmsaDirectoryPath,
			OptimalToursCsvPath:         atspData.optimalUniqueToursCsvPath,
			ToursHeatmapPath:            atspData.toursHeatmapPlotPath,
			ToursHistogramPath:          atspData.toursHistogramPlotPath,
			CmsaToursOverlapHeatmapPath: atspData.cmsaToursOverlapHeatmapPlotPath,
			AnalysisCsvPath:             atspData.cmsaSolutionAnalysisCsvPath,
			ThresholdsCsvPath:           atspData.cmsaSolutionThresholdsCsvPath,
		})

		cycleCoverConfigs = append(cycleCoverConfigs, cycleCover.InstanceConfig{
			Name:                        atspData.name,
			Dimension:                   len(atspData.matrix),
			Matrix:                      atspData.matrix,
			KnownOptimal:                atspData.knownOptimal,
			CmsaDirectoryPath:           atspData.cmsaDirectoryPath,
			OptimalToursCsvPath:         atspData.optimalUniqueToursCsvPath,
			CycleCoverEdgesCsvPath:      atspData.cycleCoverEdgesCsvPath,
			AnalysisCsvPath:             atspData.cycleCoverAnalysisCsvPath,
			ThresholdsCsvPath:           atspData.cycleCoverThresholdsCsvPath,
			CycleCoverOverlapMatrixPath: atspData.cycleCoverCmsaOverlapCsvPath,
		})
	}

	cmsaSolutionSummaryPath := filepath.Join(resultsDirectoryName, "cmsa_solution_analysis_summary.csv")
	cmsaSolutionReportPath := filepath.Join(resultsDirectoryName, "cmsa_solution_analysis_report.md")

	_, err := cmsaTours.AnalyzeInstances(cmsaTours.Config{
		Instances:      cmsaTourConfigs,
		SummaryCsvPath: cmsaSolutionSummaryPath,
		ReportPath:     cmsaSolutionReportPath,
		HighThreshold:  0.8,
		Thresholds:     cmsaTours.DefaultThresholds(),
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
		Thresholds:     cmsaTours.DefaultThresholds(),
	})
	if err != nil {
		return err
	}

	fmt.Printf("CMSA/solution analysis summary saved to %s\n", cmsaSolutionSummaryPath)
	fmt.Printf("CMSA/solution analysis report saved to %s\n", cmsaSolutionReportPath)
	fmt.Printf("Cycle-cover analysis summary saved to %s\n", cycleCoverSummaryPath)
	fmt.Printf("Cycle-cover analysis report saved to %s\n", cycleCoverReportPath)
	return nil
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

		titleSuffix := fmt.Sprintf(" (alpha=%.2f, beta=%.2f, rho=%.2f, pCmsa=%.2f)",
			statistic.alpha, statistic.beta, statistic.rho, statistic.pCmsa)

		pCmsaPlotSuffix := "_pCmsa=" + strconv.Itoa(int(100*statistic.pCmsa)) + "%"
		plotPath := plotPathPrefix + pCmsaPlotSuffix + ".png"

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
