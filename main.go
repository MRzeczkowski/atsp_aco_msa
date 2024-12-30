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
	"path"
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

func saveOptimalToursStatistics(optimalUniqueToursCsvPath string, savedTours map[string][]int) {
	header := []string{
		"Tour",
	}

	file, _ := os.Create(optimalUniqueToursCsvPath)
	defer file.Close()

	writer := csv.NewWriter(file)

	_ = writer.Write(header)
	for key := range savedTours {

		record := []string{
			key,
		}

		writer.Write(record)
	}

	writer.Flush()
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
		iterations = 1000
	}

	totalIterations := dimension * iterations

	parameters.antsNumber = int(math.Ceil(float64(dimension) * parameters.antsPercentage))
	parameters.iterations = totalIterations / parameters.antsNumber
}

func generateParameters() []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0)

	for _, useLocalSearch := range []bool{false} {
		for _, alpha := range utilities.GenerateRange(1.0, 1.0, 0.25) {
			for _, beta := range utilities.GenerateRange(5.0, 5.0, 1.0) {
				for _, rho := range utilities.GenerateRange(0.8, 0.8, 0.1) {
					for _, pBest := range utilities.GenerateRange(0.05, 0.05, 0.005) {
						for _, pCmsa := range utilities.GenerateRange(0.0, 1.0, 0.25) {
							for _, antsPercentage := range utilities.GenerateRange(0.1, 1.0, 0.1) {
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

func main() {

	// cf, cerr := os.Create("cpu.prof")
	// if cerr != nil {
	// 	fmt.Println(cerr)
	// 	return
	// }
	// pprof.StartCPUProfile(cf)
	// defer pprof.StopCPUProfile()

	tsplibDir := "tsplib_files"
	atspFilesPaths, _ := filepath.Glob(filepath.Join(tsplibDir, "*.atsp"))

	atspFilesPaths = utilities.FilterStrings(
		atspFilesPaths,
		func(file string) bool {
			var problemSize, _ = utilities.ExtractNumber(file)
			return problemSize >= 33 && problemSize < 40
		})

	resultsFolder := "results"
	numberOfExperiments := 50
	experimentParameters := generateParameters()
	for _, atspFilePath := range atspFilesPaths {

		name, dimension, matrix, knownOptimal, _ := parsing.ParseTSPLIBFile(atspFilePath)
		fmt.Println("Started processing " + name)

		atspResultsDir := filepath.Join(resultsFolder, name)
		cmsaDir := filepath.Join(atspResultsDir, "cmsa")

		cmsa, err := compositeMsa.Read(cmsaDir)

		if err != nil {
			start := time.Now()
			cmsa, err = compositeMsa.Create(matrix, cmsaDir)
			elapsed := time.Since(start)

			fmt.Printf("\tCreating %s took: %d ms\n", cmsaDir, elapsed.Milliseconds())

			if err != nil {
				fmt.Println("\tError saving CMSA: ", err)
				return
			}
		}

		experimentData := make([]ExperimentsData, 0)
		threeOptExperimentData := make([]ExperimentsData, 0)

		start := time.Now()
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
		elapsed := time.Since(start)
		fmt.Printf("\tExperiments took %dms\n", elapsed.Milliseconds())

		start = time.Now()
		statistics := calculateStatistics(experimentData)
		threeOptStatistics := calculateStatistics(threeOptExperimentData)

		if len(statistics) != 0 {
			resultFilePath := filepath.Join(atspResultsDir, "mmas") + ".csv"
			saveStatistics(resultFilePath, statistics)
		}

		if len(threeOptStatistics) != 0 {
			threeOptResultFilePath := filepath.Join(atspResultsDir, "mmas") + "+3opt" + ".csv"
			saveStatistics(threeOptResultFilePath, threeOptStatistics)
		}
		elapsed = time.Since(start)
		fmt.Printf("\tCalculating statistics took %dms\n", elapsed.Milliseconds())

		optimalUniqueToursCsvPath := path.Join(atspResultsDir, "optimal.csv")
		uniqueOptimalTours, err := getOptimalTourStatistics(optimalUniqueToursCsvPath)
		if err != nil {
			fmt.Println(err)
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

		saveOptimalToursStatistics(optimalUniqueToursCsvPath, uniqueOptimalTours)
	}

	// mf, merr := os.Create("mem.prof")
	// if merr != nil {
	// 	fmt.Println(merr)
	// 	return
	// }
	// defer mf.Close()

	// pprof.WriteHeapProfile(mf)
}
