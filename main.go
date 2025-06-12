package main

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"atsp_aco_msa/modules/algorithms/compositeMsa"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/parsing"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"image/color"
	"math"
	"os"
	"path"
	"path/filepath"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"
)

type Edge = models.Edge

type ExperimentsData struct {
	ExperimentParameters
	results []ExperimentResult
}

type ExperimentParameters struct {
	alpha, beta, rho, pCmsa, antsPercentage, localSearchAntsPercentage float64
	antsNumber, localSearchAntsNumber, iterations                      int
}

type ExperimentResult struct {
	bestAtIteration, threeOptImprovementsCount int
	bestTour                                   []int
	computationTime                            int64
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
	averageComputationTime                                                           int64
	minDeviationPerIteration, averageDeviationPerIteration, maxDeviationPerIteration []float64
}

type TourStatistics struct {
	tourId                                                                                       string
	commonalityWithCmsa, minCommonalityWithMsa, averageCommonalityWithMsa, maxCommonalityWithMsa float64
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
	sort.Slice(bestStatistics, func(i, j int) bool {
		return bestStatistics[i].averageBestDeviation < bestStatistics[j].averageBestDeviation
	})

	uniqueParameters := map[ExperimentParameters]bool{}

	for _, statistic := range bestStatistics {
		parameters := statistic.ExperimentParameters
		parameters.antsNumber = 0
		parameters.iterations = 0
		uniqueParameters[parameters] = true
	}

	alphaCounts := make(map[float64]int)
	betaCounts := make(map[float64]int)
	rhoCounts := make(map[float64]int)
	pCmsaCounts := make(map[float64]int)
	antsPercentageCounts := make(map[float64]int)
	localSearchAntsPercentageCounts := make(map[float64]int)

	// Count the occurrences
	for params := range uniqueParameters {
		alphaCounts[params.alpha]++
		betaCounts[params.beta]++
		rhoCounts[params.rho]++
		pCmsaCounts[params.pCmsa]++
		antsPercentageCounts[params.antsPercentage]++
		localSearchAntsPercentageCounts[params.localSearchAntsPercentage]++
	}

	// Markdown content
	markdown := "# Best Parameters Report\n\n"

	// Unique combinations count
	markdown += fmt.Sprintf("Found **%d** best unique parameter combinations.\n\n", len(uniqueParameters))

	// Best parameters
	markdown += "## Best Parameters\n\n"
	markdown += "| Alpha | Beta | Rho | pCmsa | AntsPercentage | LocalSearchAntsPercentage |\n"
	markdown += "|-------|------|-----|-------|----------------|---------------------------|\n"
	for parameters := range uniqueParameters {
		markdown += fmt.Sprintf("| %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |\n",
			parameters.alpha, parameters.beta, parameters.rho,
			parameters.pCmsa, parameters.antsPercentage, parameters.localSearchAntsPercentage)
	}
	markdown += "\n"

	// Parameter occurrences
	markdown += "## Parameter Values Occurrences\n\n"
	markdown += generateMarkdownCounts("Alpha", alphaCounts)
	markdown += generateMarkdownCounts("Beta", betaCounts)
	markdown += generateMarkdownCounts("Rho", rhoCounts)
	markdown += generateMarkdownCounts("pCmsa", pCmsaCounts)
	markdown += generateMarkdownCounts("AntsPercentage", antsPercentageCounts)
	markdown += generateMarkdownCounts("LocalSearchAntsPercentage", localSearchAntsPercentageCounts)

	// Parameter ranges
	markdown += "## Parameter Ranges\n\n"
	minAlpha, maxAlpha := findMinMax(alphaCounts)
	minBeta, maxBeta := findMinMax(betaCounts)
	minRho, maxRho := findMinMax(rhoCounts)
	minPCmsa, maxPCmsa := findMinMax(pCmsaCounts)
	minAntsPercentage, maxAntsPercentage := findMinMax(antsPercentageCounts)
	minLocalSearchAntsPercentage, maxLocalSearchAntsPercentage := findMinMax(localSearchAntsPercentageCounts)

	markdown += fmt.Sprintf("- **Alpha**: %.2f - %.2f\n", minAlpha, maxAlpha)
	markdown += fmt.Sprintf("- **Beta**: %.2f - %.2f\n", minBeta, maxBeta)
	markdown += fmt.Sprintf("- **Rho**: %.2f - %.2f\n", minRho, maxRho)
	markdown += fmt.Sprintf("- **pCmsa**: %.2f - %.2f\n", minPCmsa, maxPCmsa)
	markdown += fmt.Sprintf("- **AntsPercentage**: %.2f - %.2f\n", minAntsPercentage, maxAntsPercentage)
	markdown += fmt.Sprintf("- **LocalSearchAntsPercentage**: %.2f - %.2f\n", minLocalSearchAntsPercentage, maxLocalSearchAntsPercentage)

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

func saveOptimalToursStatistics(optimalUniqueToursCsvPath string, toursStatistics []TourStatistics) {
	header := []string{
		"Tour",
		"Commonality with CMSA",
		"Min commonality with MSA",
		"Avg commonality with MSA",
		"Max commonality with MSA",
	}

	file, _ := os.Create(optimalUniqueToursCsvPath)
	defer file.Close()

	writer := csv.NewWriter(file)

	_ = writer.Write(header)
	floatFormat := "%.2f"

	for _, tourStatistics := range toursStatistics {
		record := []string{
			tourStatistics.tourId,
			fmt.Sprintf(floatFormat, tourStatistics.commonalityWithCmsa),
			fmt.Sprintf(floatFormat, tourStatistics.minCommonalityWithMsa),
			fmt.Sprintf(floatFormat, tourStatistics.averageCommonalityWithMsa),
			fmt.Sprintf(floatFormat, tourStatistics.maxCommonalityWithMsa),
		}

		writer.Write(record)
	}

	writer.Flush()
}

func calculateCommonalityWithMatrix(tourEdges []Edge, matrix [][]float64) float64 {
	commonality := 0.0
	for _, edge := range tourEdges {
		if matrix[edge.From][edge.To] > 0 {
			commonality++
		}
	}
	tourLen := len(tourEdges)

	commonality = 100 * commonality / float64(tourLen)

	return commonality
}

func calculateToursStatistics(cmsaDir string, uniqueOptimalTours map[string][]int) []TourStatistics {

	cmsa, err := compositeMsa.Read(cmsaDir)
	if err != nil {
		fmt.Println(err)
	}

	msas, err := compositeMsa.ReadMsas(cmsaDir)
	if err != nil {
		fmt.Println(err)
	}

	toursStatistics := make([]TourStatistics, len(uniqueOptimalTours))
	i := 0
	for tourId, tour := range uniqueOptimalTours {
		tourEdges := models.ConvertTourToEdges(tour)

		commonalityWithCmsa := calculateCommonalityWithMatrix(tourEdges, cmsa)

		minCommonalityWithMsa := math.MaxFloat64
		averageCommonalityWithMsa := 0.0
		maxCommonalityWithMsa := -math.MaxFloat64

		for _, msa := range msas {
			commonalityWithMsa := calculateCommonalityWithMatrix(tourEdges, msa)

			if commonalityWithMsa < minCommonalityWithMsa {
				minCommonalityWithMsa = commonalityWithMsa
			}

			averageCommonalityWithMsa += commonalityWithMsa

			if commonalityWithMsa > maxCommonalityWithMsa {
				maxCommonalityWithMsa = commonalityWithMsa
			}
		}

		averageCommonalityWithMsa /= float64(len(msas))

		toursStatistics[i] = TourStatistics{tourId, commonalityWithCmsa, minCommonalityWithMsa, averageCommonalityWithMsa, maxCommonalityWithMsa}
		i++
	}

	sort.SliceStable(toursStatistics, func(i, j int) bool {
		if toursStatistics[i].commonalityWithCmsa != toursStatistics[j].commonalityWithCmsa {
			return toursStatistics[i].commonalityWithCmsa > toursStatistics[j].commonalityWithCmsa
		}

		if toursStatistics[i].averageCommonalityWithMsa != toursStatistics[j].averageCommonalityWithMsa {
			return toursStatistics[i].averageCommonalityWithMsa > toursStatistics[j].averageCommonalityWithMsa
		}

		// Tiebreaker so that results don't change unnecessarily.
		return toursStatistics[i].tourId > toursStatistics[j].tourId
	})

	return toursStatistics
}

func normalizeTour(tour []int) []int {
	if len(tour) == 0 {
		return tour
	}

	// Find the smallest element to determine the starting point
	minIndex := 0
	for i, v := range tour {
		if v < tour[minIndex] {
			minIndex = i
		}
	}

	// Rotate the tour to start at the smallest element
	normalized := make([]int, len(tour))
	for i := range tour {
		normalized[i] = tour[(minIndex+i)%len(tour)]
	}

	return normalized
}

func addUniqueTour(savedTours map[string][]int, tour []int) {
	normalizedTour := normalizeTour(tour)
	keyJson, _ := json.Marshal(normalizedTour)
	key := string(keyJson)

	if _, exists := savedTours[key]; exists {
		return
	}

	savedTours[key] = tour
}

func getOptimalTourStatistics(optimalUniqueToursCsvPath string) (map[string][]int, error) {
	result := make(map[string][]int)

	file, err := os.Open(optimalUniqueToursCsvPath)
	if err != nil {
		if os.IsNotExist(err) {
			return result, nil
		}
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	rows, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("failed to read file: %w", err)
	}

	for i, row := range rows {
		if i == 0 {
			// Skip header if present
			continue
		}

		if len(row) < 1 {
			return nil, fmt.Errorf("invalid row format at line %d", i+1)
		}

		var value []int
		if err := json.Unmarshal([]byte(row[0]), &value); err != nil {
			return nil, fmt.Errorf("failed to parse JSON at line %d: %w", i+1, err)
		}

		result[row[0]] = value
	}

	return result, nil
}

func readStatistics(csvFilePath string) ([]ExperimentsDataStatistics, error) {
	file, err := os.Open(csvFilePath)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)

	// Read the header and skip it
	_, err = reader.Read() // Read the first line (header)
	if err != nil {
		return nil, fmt.Errorf("failed to read header: %w", err)
	}

	var statistics []ExperimentsDataStatistics

	// Parse the remaining rows
	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("failed to read records: %w", err)
	}

	for _, record := range records {
		if len(record) < 16 {
			return nil, fmt.Errorf("invalid record length: %d", len(record))
		}

		alpha, _ := strconv.ParseFloat(record[0], 64)
		beta, _ := strconv.ParseFloat(record[1], 64)
		rho, _ := strconv.ParseFloat(record[2], 64)
		pCmsa, _ := strconv.ParseFloat(record[3], 64)
		antsPercentage, _ := strconv.ParseFloat(record[4], 64)
		antsNumber, _ := strconv.Atoi(record[5])
		localSearchAntsPercentage, _ := strconv.ParseFloat(record[6], 64)
		iterations, _ := strconv.Atoi(record[7])
		minBestAtIteration, _ := strconv.Atoi(record[8])
		averageBestAtIteration, _ := strconv.ParseFloat(record[9], 64)
		maxBestAtIteration, _ := strconv.Atoi(record[10])
		minThreeOptImprovementsCount, _ := strconv.Atoi(record[11])
		averageThreeOptImprovementsCount, _ := strconv.ParseFloat(record[12], 64)
		maxThreeOptImprovementsCount, _ := strconv.Atoi(record[13])
		minBestDeviation, _ := strconv.ParseFloat(record[14], 64)
		averageBestDeviation, _ := strconv.ParseFloat(record[15], 64)
		maxBestDeviation, _ := strconv.ParseFloat(record[16], 64)
		successRate, _ := strconv.ParseFloat(record[17], 64)
		averageComputationTime, _ := strconv.ParseInt(record[18], 10, 64)

		statistic := ExperimentsDataStatistics{
			ExperimentParameters: ExperimentParameters{
				alpha:                     alpha,
				beta:                      beta,
				rho:                       rho,
				pCmsa:                     pCmsa,
				antsPercentage:            antsPercentage,
				antsNumber:                antsNumber,
				localSearchAntsPercentage: localSearchAntsPercentage,
				iterations:                iterations,
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
			averageComputationTime:           averageComputationTime,
		}

		statistics = append(statistics, statistic)
	}

	return statistics, nil
}

func saveStatistics(resultCsvPath string, statistics []ExperimentsDataStatistics) {
	header := []string{
		"Alpha",
		"Beta",
		"Rho",
		"pCmsa",
		"Ants fraction",
		"Ants number",
		"Local search ants fraction",
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
		"Avg computation time [ms]"}

	file, _ := os.Create(resultCsvPath)
	defer file.Close()

	writer := csv.NewWriter(file)

	_ = writer.Write(header)

	floatFormat := "%.2f"
	for _, statistic := range statistics {

		record := []string{
			fmt.Sprintf(floatFormat, statistic.alpha),
			fmt.Sprintf(floatFormat, statistic.beta),
			fmt.Sprintf(floatFormat, statistic.rho),
			fmt.Sprintf(floatFormat, statistic.pCmsa),
			fmt.Sprintf(floatFormat, statistic.antsPercentage),
			strconv.Itoa(statistic.antsNumber),
			fmt.Sprintf(floatFormat, statistic.localSearchAntsPercentage),
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
			strconv.FormatInt(statistic.averageComputationTime, 10),
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
		var averageComputationTime int64 = 0
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

			averageComputationTime += result.computationTime
		}

		averageBestAtIteration /= resultsLen
		averageBestDeviation /= resultsLen
		successRate := 100.0 * successCounter / resultsLen
		averageComputationTime /= int64(resultsLen)

		statistics[i] = ExperimentsDataStatistics{
			data.ExperimentParameters,
			minBestAtIteration,
			averageBestAtIteration,
			maxBestAtIteration,
			minThreeOptImprovementsCount,
			averageThreeOptImprovementsCount,
			maxThreeOptImprovementsCount,
			minBestDeviation,
			averageBestDeviation,
			maxBestDeviation,
			successRate,
			averageComputationTime,
			minDeviationPerIteration, averageDeviationPerIteration, maxDeviationPerIteration}
	}

	sort.SliceStable(statistics, func(i, j int) bool {
		if statistics[i].averageBestDeviation != statistics[j].averageBestDeviation {
			return statistics[i].averageBestDeviation < statistics[j].averageBestDeviation
		}

		if statistics[i].successRate != statistics[j].successRate {
			return statistics[i].successRate > statistics[j].successRate
		}

		if statistics[i].localSearchAntsPercentage != statistics[j].localSearchAntsPercentage {
			return statistics[i].localSearchAntsPercentage < statistics[j].localSearchAntsPercentage
		}

		return statistics[i].averageBestAtIteration < statistics[j].averageBestAtIteration
	})

	return statistics
}

func runExperiments(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, cmsa [][]float64) []ExperimentResult {
	results := make([]ExperimentResult, numberOfRuns)

	aco := aco.NewACO(
		parameters.alpha,
		parameters.beta,
		parameters.rho,
		parameters.pCmsa,
		parameters.antsNumber,
		parameters.localSearchAntsNumber,
		parameters.iterations,
		knownOptimal,
		matrix,
		cmsa)

	for i := 0; i < numberOfRuns; i++ {

		start := time.Now()
		aco.Run()
		elapsed := time.Since(start)

		results[i] = ExperimentResult{aco.BestAtIteration, aco.ThreeOptImprovementsCount, aco.BestTour, elapsed.Milliseconds(), aco.DeviationPerIteration}
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

	totalIterations := dimension * iterations

	parameters.antsNumber = int(math.Ceil(float64(dimension) * parameters.antsPercentage))
	parameters.localSearchAntsNumber = int(math.Floor(float64(parameters.antsNumber) * parameters.localSearchAntsPercentage))
	parameters.iterations = totalIterations / parameters.antsNumber
}

func generateParameters() []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0)

	for _, alpha := range utilities.GenerateRange(1.0, 1.0, 0.25) {
		for _, beta := range utilities.GenerateRange(2.0, 2.0, 1.0) {
			for _, rho := range utilities.GenerateRange(0.8, 0.8, 0.1) {
				for _, pCmsa := range utilities.GenerateRange(0.0, 0.0, 0.25) {
					for _, antsPercentage := range utilities.GenerateRange(0.8, 0.8, 0.1) {
						for _, localSearchAntsPercentage := range utilities.GenerateRange(0.0, 0.0, 0.5) {
							parameters = append(parameters,
								ExperimentParameters{
									alpha, beta, rho, pCmsa, antsPercentage, localSearchAntsPercentage, 0, 0, 0,
								})
						}
					}
				}
			}
		}
	}

	return parameters
}

var resultsDirectoryName = "results"
var resultFileName = "result.csv"

type AtspData struct {
	name         string
	matrix       [][]float64
	knownOptimal float64

	cmsaDirectoryPath,

	cmsaHeatmapPlotPath, cmsaHistogramPlotPath,

	resultFilePath,
	resultPlotFilePrefix,

	optimalUniqueToursCsvPath, toursHeatmapPlotPath, toursHistogramPlotPath string
}

func makeAtspData(name string, matrix [][]float64, knownOptimal float64) AtspData {
	name = strings.TrimSuffix(name, ".atsp")
	resultsDirectoryPath := filepath.Join(resultsDirectoryName, name)
	cmsaDirectoryPath := filepath.Join(resultsDirectoryPath, "cmsa")
	plotsDirectoryPath := filepath.Join(resultsDirectoryPath, "plots")

	cmsaHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "cmsa_heatmap.png")
	cmsaHistogramPlotPath := filepath.Join(plotsDirectoryPath, "cmsa_histogram.png")

	resultFilePath := filepath.Join(resultsDirectoryPath, resultFileName)
	resultPlotFilePrefix := filepath.Join(plotsDirectoryPath, "best_result")

	optimalUniqueToursCsvPath := filepath.Join(resultsDirectoryPath, "solutions.csv")
	toursHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "tours_heatmap.png")
	toursHistogramPlotPath := filepath.Join(plotsDirectoryPath, "tours_histogram.png")

	return AtspData{
		name,
		matrix,
		knownOptimal,

		cmsaDirectoryPath,

		cmsaHeatmapPlotPath, cmsaHistogramPlotPath,

		resultFilePath,
		resultPlotFilePrefix,

		optimalUniqueToursCsvPath, toursHeatmapPlotPath, toursHistogramPlotPath,
	}
}

func main() {
	cf, cerr := os.Create("cpu.prof")
	if cerr != nil {
		fmt.Println(cerr)
		return
	}
	pprof.StartCPUProfile(cf)
	defer pprof.StopCPUProfile()

	tsplibDir := "tsplib_files"
	atspFilesPaths, _ := filepath.Glob(filepath.Join(tsplibDir, "*.atsp"))

	atspFilesPaths = utilities.FilterStrings(
		atspFilesPaths,
		func(filePath string) bool {
			var files = []string{
				"atex1.atsp",
				"atex3.atsp",
				"atex4.atsp",
				"atex5.atsp",
				"br17.atsp",
				// "code198.atsp",
				// "crane100_0.atsp",
				// "crane100_1.atsp",
				// "crane100_2.atsp",
				"crane66_0.atsp",
				"crane66_1.atsp",
				"crane66_2.atsp",
				// "dc112.atsp",
				// "dc126.atsp",
				// "dc134.atsp",
				// "dc176.atsp",
				// "dc188.atsp",
				"ft53.atsp",
				"ft70.atsp",
				// "ftv100.atsp",
				// "ftv110.atsp",
				// "ftv120.atsp",
				// "ftv130.atsp",
				// "ftv140.atsp",
				// "ftv150.atsp",
				// "ftv160.atsp",
				// "ftv170.atsp",
				"ftv33.atsp",
				"ftv35.atsp",
				"ftv38.atsp",
				"ftv44.atsp",
				"ftv47.atsp",
				"ftv55.atsp",
				"ftv64.atsp",
				"ftv70.atsp",
				"ftv90.atsp",
				"p43.atsp",
				// "rbg323.atsp",
				// "rbg358.atsp",
				// "rbg403.atsp",
				// "rbg443.atsp",
				"ry48p.atsp",
				// "td100_1.atsp",
			}

			inputFileName := path.Base(filePath)
			for _, file := range files {
				// assumedProblemSize, _ := utilities.ExtractNumber(file)
				if inputFileName == file {
					return true
				}
			}

			return false
		})

	atspsData := make([]AtspData, len(atspFilesPaths))

	for i, atspFilePath := range atspFilesPaths {
		name, matrix, knownOptimal, _ := parsing.ParseTSPLIBFile(atspFilePath)
		atspsData[i] = makeAtspData(name, matrix, knownOptimal)
	}

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
				fmt.Println("\tError saving CMSA: ", err)
				return
			}
		}

		cmsaHeatmapPlotTitle := name + " CMSA heatmap"

		err = utilities.SaveHeatmapFromMatrix(cmsa, cmsaHeatmapPlotTitle, atspData.cmsaHeatmapPlotPath)
		if err != nil {
			fmt.Println(err)
		}

		dataForHistogram := filterZeroes(flattenMatrix(cmsa))
		cmsaHistogramPlotTitle := name + " CMSA histogram"

		dimension := len(matrix)
		err = utilities.SaveHistogramFromData(dataForHistogram, dimension-1, cmsaHistogramPlotTitle, atspData.cmsaHistogramPlotPath)
		if err != nil {
			fmt.Println(err)
		}
	}

	experimentParameters := generateParameters()
	for _, atspData := range atspsData {
		matrix := atspData.matrix
		knownOptimal := atspData.knownOptimal
		dimension := len(matrix)

		cmsa, err := compositeMsa.Read(atspData.cmsaDirectoryPath)
		if err != nil {
			fmt.Println(err)
			return
		}

		numberOfExperiments := 10
		experimentData := make([]ExperimentsData, 0)

		for _, parameters := range experimentParameters {
			setDimensionDependantParameters(dimension, &parameters)
			results := runExperiments(numberOfExperiments, parameters, knownOptimal, matrix, cmsa)
			data := ExperimentsData{parameters, results}

			experimentData = append(experimentData, data)
		}

		statistics := calculateStatistics(experimentData)
		if len(statistics) != 0 {
			saveStatistics(atspData.resultFilePath, statistics)
			saveExperimentPlots(statistics, "MMAS deviation per iteration", atspData.resultPlotFilePrefix)
		}

		uniqueOptimalTours, err := getOptimalTourStatistics(atspData.optimalUniqueToursCsvPath)
		if err != nil {
			fmt.Println(err)
			return
		}

		for _, data := range experimentData {
			for _, result := range data.results {
				if result.deviationPerIteration[result.bestAtIteration] == 0.0 {
					addUniqueTour(uniqueOptimalTours, result.bestTour)
				}
			}
		}

		cmsaDirectoryPath := atspData.cmsaDirectoryPath
		toursStatistics := calculateToursStatistics(cmsaDirectoryPath, uniqueOptimalTours)
		saveOptimalToursStatistics(atspData.optimalUniqueToursCsvPath, toursStatistics)
	}

	for _, atspData := range atspsData {
		name := atspData.name
		dimension := len(atspData.matrix)

		uniqueOptimalTours, err := getOptimalTourStatistics(atspData.optimalUniqueToursCsvPath)

		if err != nil {
			fmt.Println(err)
			return
		}

		toursCount := len(uniqueOptimalTours)
		if toursCount > 0 {
			toursMatrix := make([][]float64, dimension)
			for i := 0; i < dimension; i++ {
				toursMatrix[i] = make([]float64, dimension)
			}

			for _, tour := range uniqueOptimalTours {
				n := len(tour)
				for i := 0; i < n-1; i++ {
					start, end := tour[i], tour[i+1]
					toursMatrix[start][end]++
				}

				last, first := tour[n-1], tour[0]
				toursMatrix[last][first]++
			}

			toursHeatmapPlotTitle := name + " tours heatmap"
			err = utilities.SaveHeatmapFromMatrix(toursMatrix, toursHeatmapPlotTitle, atspData.toursHeatmapPlotPath)
			if err != nil {
				fmt.Println(err)
			}

			if toursCount > 1 {
				dataForHistogram := filterZeroes(flattenMatrix(toursMatrix))
				toursHistogramPlotTitle := name + " tours histogram"

				numberOfBins := countUniqueValues(dataForHistogram)

				err = utilities.SaveHistogramFromData(dataForHistogram, numberOfBins, toursHistogramPlotTitle, atspData.toursHistogramPlotPath)
				if err != nil {
					fmt.Println(err)
				}
			}
		}
	}

	resultsFilePaths, _ := filepath.Glob(filepath.Join(resultsDirectoryName, "*", resultFileName))

	bestStatistics := getBestStatisticsFromFiles(resultsFilePaths)

	saveBestParametersInfo("best_parameters_report.md", bestStatistics)
}

func saveExperimentPlots(statistics []ExperimentsDataStatistics, plotTitle, plotPathPrefix string) {
	bestStatistic := statistics[0]

	for _, statistic := range statistics {
		if statistic.alpha != bestStatistic.alpha ||
			statistic.beta != bestStatistic.beta ||
			statistic.rho != bestStatistic.rho ||
			statistic.antsPercentage != bestStatistic.antsPercentage ||
			statistic.localSearchAntsPercentage != bestStatistic.localSearchAntsPercentage {
			continue
		}

		minDeviationPlotData := utilities.LinePlotData{Name: "min deviation", Color: color.RGBA{G: 255, A: 255}, Values: statistic.minDeviationPerIteration}
		avgDeviationPlotData := utilities.LinePlotData{Name: "avg deviation", Color: color.RGBA{B: 255, A: 255}, Values: statistic.averageDeviationPerIteration}
		maxDeviationPlotData := utilities.LinePlotData{Name: "max deviation", Color: color.RGBA{R: 255, A: 255}, Values: statistic.maxDeviationPerIteration}
		lines := []utilities.LinePlotData{minDeviationPlotData, avgDeviationPlotData, maxDeviationPlotData}

		titleSuffix := fmt.Sprintf(" (alpha=%.2f, beta=%.2f, rho=%.2f, pCmsa=%.2f, antsFraction=%.2f, localSearchAntsFraction=%.2f)",
			statistic.alpha, statistic.beta, statistic.rho, statistic.pCmsa, statistic.antsPercentage, statistic.localSearchAntsPercentage)

		pCmsaPlotSuffix := "_pCmsa=" + strconv.Itoa(int(100*statistic.pCmsa)) + "%"
		plotPath := plotPathPrefix + pCmsaPlotSuffix + ".png"

		utilities.SaveLinePlotFromData(lines, plotTitle+titleSuffix, plotPath)
	}
}

func getBestStatisticsFromFiles(resultsFilePaths []string) []ExperimentsDataStatistics {
	topNumber := 1
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

func countUniqueValues(data []float64) int {
	unique := make(map[float64]struct{})
	for _, value := range data {
		unique[value] = struct{}{}
	}

	return len(unique)
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
