package main

import (
	"atsp_aco_msa/modules/algorithms/aco"
	"atsp_aco_msa/modules/algorithms/cmsa"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/parsing"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"time"
)

const NumberOfRuns int = 50

type Edge = models.Edge

type ExperimentData struct {
	alpha, beta, alpha1, rho, q0, cmsaP float64
	ExperimentResult
}

type ExperimentResult struct {
	bestAtIteration                                         int
	bestLength, deviation, successRate, commonalityWithCmsa float64
	averageTime                                             int64
}

func (f ExperimentData) ToCSVRow() []string {
	return []string{
		fmt.Sprintf("%.2f", f.alpha),
		fmt.Sprintf("%.2f", f.beta),
		fmt.Sprintf("%.2f", f.alpha1),
		fmt.Sprintf("%.2f", f.rho),
		fmt.Sprintf("%.2f", f.q0),
		fmt.Sprintf("%.2f", f.cmsaP),
		strconv.Itoa(f.ExperimentResult.bestAtIteration),
		fmt.Sprintf("%.2f", f.ExperimentResult.bestLength),
		fmt.Sprintf("%f", f.ExperimentResult.deviation),
		fmt.Sprintf("%.2f", f.ExperimentResult.successRate),
		fmt.Sprintf("%f", f.ExperimentResult.commonalityWithCmsa),
	}
}

var optimalSolutions = map[string]float64{
	"br17":   39,
	"ft53":   6905, // [49,52,50,48,29,28,25,27,26,3,13,11,10,12,14,41,47,42,46,43,45,44,34,32,33,31,30,0,4,2,17,16,15,37,39,38,36,35,40,21,20,24,23,22,19,18,1,8,9,7,6,5,51]
	"ft70":   38673,
	"ftv33":  1286,
	"ftv35":  1473,
	"ftv38":  1530,
	"ftv44":  1613,
	"ftv47":  1776,
	"ftv55":  1608,
	"ftv64":  1839,
	"ftv70":  1950,
	"ftv170": 2755,
	"p43":    5620,
	"rbg323": 1326,
	"rbg358": 1163,
	"rbg403": 2465,
	"rbg443": 2720,
	"ry48p":  14422,
}

func runExperiment(name string, dimension, iterations int, alpha, beta, alpha1, rho, q0, cmsaP float64, matrix, cmsa [][]float64) ExperimentResult {

	var totalBestLength float64
	var totalElapsedTime time.Duration

	bestLength := math.MaxFloat64
	var bestPath []int
	successCounter := 0.0
	bestAtIteration := math.MaxInt

	knownOptimal := optimalSolutions[name]

	ants := dimension

	for i := 0; i < NumberOfRuns; i++ {
		aco := aco.NewACO(alpha, beta, alpha1, rho, q0, cmsaP, ants, iterations, matrix, cmsa)
		start := time.Now()
		aco.Run()
		elapsed := time.Since(start)

		totalBestLength += aco.BestLength
		totalElapsedTime += elapsed

		if aco.BestAtIteration < bestAtIteration && aco.BestLength < bestLength {
			bestAtIteration = aco.BestAtIteration
			bestLength = aco.BestLength
			bestPath = aco.BestPath
		}

		if aco.BestLength == knownOptimal {
			successCounter++
		}
	}

	averageBestLength := totalBestLength / float64(NumberOfRuns)
	averageTime := totalElapsedTime / time.Duration(NumberOfRuns)
	deviation := 100 * (averageBestLength - knownOptimal) / knownOptimal
	successRate := 100 * successCounter / float64(NumberOfRuns)

	bestPathEdges := make([]Edge, len(bestPath))

	for i := 0; i < dimension-1; i++ {
		bestPathEdges[i] = Edge{From: bestPath[i], To: bestPath[i+1]}
	}

	last, first := bestPath[dimension-1], bestPath[0]
	bestPathEdges[last] = Edge{From: bestPath[last], To: bestPath[first]}

	commonalityWithCmsa := 0.0

	for _, edge := range bestPathEdges {
		if cmsa[edge.From][edge.To] > 0 {
			commonalityWithCmsa++
		}
	}

	commonalityWithCmsa = 100 * commonalityWithCmsa / float64(dimension-1)

	return ExperimentResult{bestAtIteration, bestLength, deviation, successRate, commonalityWithCmsa, averageTime.Milliseconds()}
}

func tryFindSolution(path string) {
	name, dimension, matrix, err := parsing.ParseTSPLIBFile(path)

	if err != nil {
		fmt.Println("Error parsing file:", path, err)
		return
	}

	cmsa := cmsa.CreateCMSA(dimension, matrix)

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

	file, err := os.Create(filepath.Join("results", name) + ".csv")
	if err != nil {
		log.Fatalf("Failed to create file: %s", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)

	header := []string{
		"Alpha",
		"Beta",
		"Alpha1",
		"Rho",
		"Q0",
		"CMSA probability",
		"Best at iteration",
		"Best length",
		"Deviation",
		"Success rate",
		"Commonality with CMSA"}

	err = writer.Write(header)
	if err != nil {
		log.Fatalf("Failed to write header: %s", err)
	}

	for _, alpha := range utilities.GenerateRange(1.0, 1.0, 0.25) {
		for _, beta := range utilities.GenerateRange(5.0, 5.0, 1.0) {
			for _, alpha1 := range utilities.GenerateRange(0.9, 0.9, 0.1) {
				for _, rho := range utilities.GenerateRange(0.5, 0.5, 0.25) {
					for _, q0 := range utilities.GenerateRange(0.1, 0.1, 0.25) {
						for _, cmsaP := range utilities.GenerateRange(0.0, 1.0, 0.25) {

							// 1. Analiza grafów
							// 3. Parametry + dopracowanie heurystyki

							result := runExperiment(name, dimension, iterations, alpha, beta, alpha1, q0, rho, cmsaP, matrix, cmsa)

							data := ExperimentData{
								alpha, beta, alpha1, rho, q0, cmsaP, result,
							}

							cswRow := data.ToCSVRow()

							err := writer.Write(cswRow)
							if err != nil {
								log.Fatalf("Failed to write record: %s", err)
							}
						}
					}
				}
			}
		}
	}

	writer.Flush()
	if err := writer.Error(); err != nil {
		log.Fatalf("Error while flushing the data: %s", err)
	}
}

func main() {
	dir := "tsp_files"
	paths, err := filepath.Glob(filepath.Join(dir, "*.atsp"))

	if err != nil {
		fmt.Println("Error finding files:", err)
		return
	}

	if len(paths) == 0 {
		fmt.Println("No files found in the directory.")
		return
	}

	paths = utilities.FilterStrings(
		paths,
		func(file string) bool {
			var problemSize, _ = utilities.ExtractNumber(file)
			return problemSize <= 170
		})

	for _, path := range paths {
		tryFindSolution(path)
	}
}
