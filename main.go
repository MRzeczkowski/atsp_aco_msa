package main

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"atsp_aco_msa/modules/algorithms/compositeMsa"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/parsing"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
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
	averageBestAtIteration, averageBestDeviation, successRate float64
	averageComputationTime                                    int64
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
		"Avg best at iteration",
		"Avg best deviation",
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
			fmt.Sprintf(floatFormat, statistic.averageBestAtIteration),
			fmt.Sprintf(floatFormat, statistic.averageBestDeviation),
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
		successCounter := 0.0
		averageBestAtIteration := 0.0
		averageBestDeviation := 0.0
		var averageComputationTime int64 = 0

		for _, result := range data.results {
			averageBestAtIteration += float64(result.bestAtIteration)
			bestDeviation := result.deviationPerIteration[result.bestAtIteration]
			averageBestDeviation += bestDeviation
			averageComputationTime += result.computationTime

			if bestDeviation == 0 {
				successCounter++
			}
		}

		resultsLen := float64(len(data.results))

		averageBestAtIteration /= resultsLen
		averageBestDeviation /= resultsLen
		successRate := 100.0 * successCounter / resultsLen
		averageComputationTime /= int64(resultsLen)

		statistics[i] = ExperimentsDataStatistics{
			data.ExperimentParameters, averageBestAtIteration, averageBestDeviation, successRate, averageComputationTime}
	}

	sort.SliceStable(statistics, func(i, j int) bool {
		return statistics[i].averageBestDeviation < statistics[j].averageBestDeviation && statistics[i].successRate < statistics[j].successRate
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
					for _, pBest := range utilities.GenerateRange(0.05, 0.05, 0.01) {
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

	dir := "tsp_files"
	paths, _ := filepath.Glob(filepath.Join(dir, "*.atsp"))

	paths = utilities.FilterStrings(
		paths,
		func(file string) bool {
			var problemSize, _ = utilities.ExtractNumber(file)
			return problemSize < 50
		})

	numberOfExperiments := 50
	experimentParameters := generateParameters()
	for _, path := range paths {

		name, dimension, matrix, knownOptimal, _ := parsing.ParseTSPLIBFile(path)
		fmt.Println("Started processing " + name)

		cmsaCSVPath := filepath.Join("results", name) + "_cmsa.csv"

		cmsa, err := compositeMsa.ReadFromCsv(cmsaCSVPath)

		if err != nil {
			start := time.Now()
			cmsa = compositeMsa.CreateFromData(matrix)
			elapsed := time.Since(start)

			fmt.Printf("Creating %s took: %d ms\n", cmsaCSVPath, elapsed.Milliseconds())

			err := compositeMsa.SaveToCsv(cmsa, cmsaCSVPath)

			if err != nil {
				fmt.Println("Error saving CMSA to file:", cmsaCSVPath, err)
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
			resultFilePath := filepath.Join("results", name) + ".csv"
			saveStatistics(resultFilePath, statistics)
		}

		if len(threeOptStatistics) != 0 {
			threeOptResultFilePath := filepath.Join("results", name) + "+3opt" + ".csv"
			saveStatistics(threeOptResultFilePath, threeOptStatistics)
		}
		elapsed = time.Since(start)
		fmt.Printf("\tCalculating statistics took %dms\n", elapsed.Milliseconds())
	}

	// mf, merr := os.Create("mem.prof")
	// if merr != nil {
	// 	fmt.Println(merr)
	// 	return
	// }
	// defer mf.Close()

	// pprof.WriteHeapProfile(mf)
}
