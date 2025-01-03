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
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"time"
)

type Edge = models.Edge

type ExperimentsData struct {
	ExperimentParameters
	results []ExperimentResult
}

type ExperimentParameters struct {
	useLocalSearch                                 bool
	alpha, beta, rho, pBest, pCmsa, antsPercentage float64
	antsNumber, iterations                         int
}

type ExperimentResult struct {
	bestAtIteration       int
	bestTour              []int
	computationTime       int64
	deviationPerIteration []float64
}

type ExperimentsDataStatistics struct {
	ExperimentParameters
	minBestAtIteration, maxBestAtIteration                                                        int
	averageBestAtIteration, minBestDeviation, averageBestDeviation, maxBestDeviation, successRate float64
	averageComputationTime                                                                        int64
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

func saveBestParametersInfo(resultsFolder string, bestStatistics []ExperimentsDataStatistics) {
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
	pBestCounts := make(map[float64]int)
	pCmsaCounts := make(map[float64]int)
	antsPercentageCounts := make(map[float64]int)

	// Count the occurrences
	for params := range uniqueParameters {
		alphaCounts[params.alpha]++
		betaCounts[params.beta]++
		rhoCounts[params.rho]++
		pBestCounts[params.pBest]++
		pCmsaCounts[params.pCmsa]++
		antsPercentageCounts[params.antsPercentage]++
	}

	// Markdown content
	markdown := "# Best Parameters Report\n\n"

	// Unique combinations count
	markdown += fmt.Sprintf("Found **%d** best unique parameter combinations.\n\n", len(uniqueParameters))

	// Best parameters
	markdown += "## Best Parameters\n\n"
	markdown += "| Alpha | Beta | Rho | pBest | pCmsa | AntsPercentage |\n"
	markdown += "|-------|------|-----|-------|-------|----------------|\n"
	for parameters := range uniqueParameters {
		markdown += fmt.Sprintf("| %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |\n",
			parameters.alpha, parameters.beta, parameters.rho,
			parameters.pBest, parameters.pCmsa, parameters.antsPercentage)
	}
	markdown += "\n"

	// Parameter occurrences
	markdown += "## Parameter Values Occurrences\n\n"
	markdown += generateMarkdownCounts("Alpha", alphaCounts)
	markdown += generateMarkdownCounts("Beta", betaCounts)
	markdown += generateMarkdownCounts("Rho", rhoCounts)
	markdown += generateMarkdownCounts("PBest", pBestCounts)
	markdown += generateMarkdownCounts("PCmsa", pCmsaCounts)
	markdown += generateMarkdownCounts("AntsPercentage", antsPercentageCounts)

	// Parameter ranges
	markdown += "## Parameter Ranges\n\n"
	minAlpha, maxAlpha := findMinMax(alphaCounts)
	minBeta, maxBeta := findMinMax(betaCounts)
	minRho, maxRho := findMinMax(rhoCounts)
	minPBest, maxPBest := findMinMax(pBestCounts)
	minPCmsa, maxPCmsa := findMinMax(pCmsaCounts)
	minAntsPercentage, maxAntsPercentage := findMinMax(antsPercentageCounts)

	markdown += fmt.Sprintf("- **Alpha**: %.2f - %.2f\n", minAlpha, maxAlpha)
	markdown += fmt.Sprintf("- **Beta**: %.2f - %.2f\n", minBeta, maxBeta)
	markdown += fmt.Sprintf("- **Rho**: %.2f - %.2f\n", minRho, maxRho)
	markdown += fmt.Sprintf("- **pBest**: %.2f - %.2f\n", minPBest, maxPBest)
	markdown += fmt.Sprintf("- **PCmsa**: %.2f - %.2f\n", minPCmsa, maxPCmsa)
	markdown += fmt.Sprintf("- **AntsPercentage**: %.2f - %.2f\n", minAntsPercentage, maxAntsPercentage)

	// Save to a file
	filename := "best_parameters_report.md"
	reportPath := filepath.Join(resultsFolder, filename)
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
		pBest, _ := strconv.ParseFloat(record[3], 64)
		pCmsa, _ := strconv.ParseFloat(record[4], 64)
		antsPercentage, _ := strconv.ParseFloat(record[5], 64)
		antsNumber, _ := strconv.Atoi(record[6])
		iterations, _ := strconv.Atoi(record[7])
		minBestAtIteration, _ := strconv.Atoi(record[8])
		averageBestAtIteration, _ := strconv.ParseFloat(record[9], 64)
		maxBestAtIteration, _ := strconv.Atoi(record[10])
		minBestDeviation, _ := strconv.ParseFloat(record[11], 64)
		averageBestDeviation, _ := strconv.ParseFloat(record[12], 64)
		maxBestDeviation, _ := strconv.ParseFloat(record[13], 64)
		successRate, _ := strconv.ParseFloat(record[14], 64)
		averageComputationTime, _ := strconv.ParseInt(record[15], 10, 64)

		statistic := ExperimentsDataStatistics{
			ExperimentParameters: ExperimentParameters{
				alpha:          alpha,
				beta:           beta,
				rho:            rho,
				pBest:          pBest,
				pCmsa:          pCmsa,
				antsPercentage: antsPercentage,
				antsNumber:     antsNumber,
				iterations:     iterations,
				useLocalSearch: false, // Assuming default value since it's not in the CSV
			},
			minBestAtIteration:     minBestAtIteration,
			maxBestAtIteration:     maxBestAtIteration,
			averageBestAtIteration: averageBestAtIteration,
			minBestDeviation:       minBestDeviation,
			averageBestDeviation:   averageBestDeviation,
			maxBestDeviation:       maxBestDeviation,
			successRate:            successRate,
			averageComputationTime: averageComputationTime,
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
		"pBest",
		"pCmsa",
		"Ants fraction",
		"Ants number",
		"Iterations",
		"Min best at iteration",
		"Avg best at iteration",
		"Max best at iteration",
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
			fmt.Sprintf(floatFormat, statistic.pBest),
			fmt.Sprintf(floatFormat, statistic.pCmsa),
			fmt.Sprintf(floatFormat, statistic.antsPercentage),
			strconv.Itoa(statistic.antsNumber),
			strconv.Itoa(statistic.iterations),
			strconv.Itoa(statistic.minBestAtIteration),
			fmt.Sprintf(floatFormat, statistic.averageBestAtIteration),
			strconv.Itoa(statistic.maxBestAtIteration),
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
		maxBestAtIteration := -math.MaxInt
		averageBestAtIteration := 0.0
		minBestDeviation := math.MaxFloat64
		averageBestDeviation := 0.0
		maxBestDeviation := -math.MaxFloat64
		successCounter := 0.0
		var averageComputationTime int64 = 0

		for _, result := range data.results {

			if result.bestAtIteration < minBestAtIteration {
				minBestAtIteration = result.bestAtIteration
			}

			if result.bestAtIteration > maxBestAtIteration {
				maxBestAtIteration = result.bestAtIteration
			}

			averageBestAtIteration += float64(result.bestAtIteration)

			bestDeviation := result.deviationPerIteration[result.bestAtIteration]

			if bestDeviation < minBestDeviation {
				minBestDeviation = bestDeviation
			}

			averageBestDeviation += bestDeviation

			if bestDeviation > maxBestDeviation {
				maxBestDeviation = bestDeviation
			}

			if bestDeviation == 0 {
				successCounter++
			}

			averageComputationTime += result.computationTime
		}

		resultsLen := float64(len(data.results))

		averageBestAtIteration /= resultsLen
		averageBestDeviation /= resultsLen
		successRate := 100.0 * successCounter / resultsLen
		averageComputationTime /= int64(resultsLen)

		statistics[i] = ExperimentsDataStatistics{
			data.ExperimentParameters,
			minBestAtIteration,
			maxBestAtIteration,
			averageBestAtIteration,
			minBestDeviation,
			averageBestDeviation,
			maxBestDeviation,
			successRate,
			averageComputationTime}
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

func runExperiments(numberOfRuns int, parameters ExperimentParameters, knownOptimal float64, matrix, cmsa [][]float64) []ExperimentResult {
	results := make([]ExperimentResult, numberOfRuns)

	for i := 0; i < numberOfRuns; i++ {

		aco := aco.NewACO(
			parameters.useLocalSearch,
			parameters.alpha,
			parameters.beta,
			parameters.rho,
			parameters.pBest,
			parameters.pCmsa,
			parameters.antsNumber,
			parameters.iterations,
			knownOptimal,
			matrix,
			cmsa)

		start := time.Now()
		aco.Run()
		elapsed := time.Since(start)

		results[i] = ExperimentResult{aco.BestAtIteration, aco.BestTour, elapsed.Milliseconds(), aco.DeviationPerIteration}
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
		iterations = 10000
	}

	totalIterations := dimension * iterations

	parameters.antsNumber = int(math.Ceil(float64(dimension) * parameters.antsPercentage))
	parameters.iterations = totalIterations / parameters.antsNumber
}

func generateParameters() []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0)

	for _, useLocalSearch := range []bool{true} {
		for _, alpha := range utilities.GenerateRange(1.0, 1.0, 0.25) {
			for _, beta := range utilities.GenerateRange(5.0, 5.0, 1.0) {
				for _, rho := range utilities.GenerateRange(0.8, 0.8, 0.1) {
					for _, pBest := range utilities.GenerateRange(0.05, 0.05, 0.005) {
						for _, pCmsa := range utilities.GenerateRange(0.25, 0.75, 0.25) {
							for _, antsPercentage := range utilities.GenerateRange(0.5, 1.0, 0.1) {
								parameters = append(parameters,
									ExperimentParameters{
										useLocalSearch, alpha, beta, rho, pBest, pCmsa, antsPercentage, 0, 0,
									})
							}
						}
					}
				}
			}
		}
	}

	return parameters
}

var resultsDirectoryName = "results"
var resultFileName = "mmas.csv"
var threeOptResultFileName = "mmas+3opt.csv"

type AtspData struct {
	name         string
	matrix       [][]float64
	knownOptimal float64

	cmsaDirectoryPath,

	cmsaHeatmapPlotPath, cmsaHistogramPlotPath,

	resultFilePath, threeOptResultFilePath,

	optimalUniqueToursCsvPath, toursHeatmapPlotPath, toursHistogramPlotPath string
}

func makeAtspData(name string, matrix [][]float64, knownOptimal float64) AtspData {
	resultsDirectoryPath := filepath.Join(resultsDirectoryName, name)
	cmsaDirectoryPath := filepath.Join(resultsDirectoryPath, "cmsa")
	plotsDirectoryPath := filepath.Join(resultsDirectoryPath, "plots")

	cmsaHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "cmsa_heatmap.png")
	cmsaHistogramPlotPath := filepath.Join(plotsDirectoryPath, "cmsa_histogram.png")

	resultFilePath := filepath.Join(resultsDirectoryPath, resultFileName)
	threeOptResultFilePath := filepath.Join(resultsDirectoryPath, threeOptResultFileName)

	optimalUniqueToursCsvPath := filepath.Join(resultsDirectoryPath, "solutions.csv")
	toursHeatmapPlotPath := filepath.Join(plotsDirectoryPath, "tours_heatmap.png")
	toursHistogramPlotPath := filepath.Join(plotsDirectoryPath, "tours_histogram.png")

	return AtspData{
		name,
		matrix,
		knownOptimal,

		cmsaDirectoryPath,

		cmsaHeatmapPlotPath, cmsaHistogramPlotPath,

		resultFilePath, threeOptResultFilePath,

		optimalUniqueToursCsvPath, toursHeatmapPlotPath, toursHistogramPlotPath,
	}
}

func main() {
	tsplibDir := "tsplib_files"
	atspFilesPaths, _ := filepath.Glob(filepath.Join(tsplibDir, "*.atsp"))

	atspFilesPaths = utilities.FilterStrings(
		atspFilesPaths,
		func(file string) bool {
			var problemSize, _ = utilities.ExtractNumber(file)
			return problemSize != 17 && problemSize == 170
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

		numberOfExperiments := 1
		experimentData := make([]ExperimentsData, 0)
		threeOptExperimentData := make([]ExperimentsData, 0)

		for _, parameters := range experimentParameters {
			setDimensionDependantParameters(dimension, &parameters)
			results := runExperiments(numberOfExperiments, parameters, knownOptimal, matrix, cmsa)
			data := ExperimentsData{parameters, results}

			if parameters.useLocalSearch {
				threeOptExperimentData = append(threeOptExperimentData, data)
			} else {
				experimentData = append(experimentData, data)
			}
		}

		statistics := calculateStatistics(experimentData)
		threeOptStatistics := calculateStatistics(threeOptExperimentData)

		if len(statistics) != 0 {
			saveStatistics(atspData.resultFilePath, statistics)
		}

		if len(threeOptStatistics) != 0 {
			saveStatistics(atspData.threeOptResultFilePath, threeOptStatistics)
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

		for _, data := range threeOptExperimentData {
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

	bestStatistics := make([]ExperimentsDataStatistics, 0)

	topNumber := 1
	resultsFilePaths, _ := filepath.Glob(filepath.Join(resultsDirectoryName, "*", resultFileName))
	threeOptResultsFilePaths, _ := filepath.Glob(filepath.Join(resultsDirectoryName, "*", threeOptResultFileName))

	for _, paths := range [][]string{resultsFilePaths, threeOptResultsFilePaths} {
		for _, path := range paths {
			statistics, err := readStatistics(path)
			if err != nil {
				fmt.Println(err)
				continue
			}

			top := statistics[:topNumber]
			bestStatistics = append(bestStatistics, top...)
		}
	}

	saveBestParametersInfo(resultsDirectoryName, bestStatistics)
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
