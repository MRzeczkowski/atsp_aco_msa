package aco

import (
	"atsp_aco_msa/modules/utilities"
	"math"
	"math/rand"
	"sync"
)

type ACO struct {
	alpha, beta, rho, pDec, pCmsa                       float64
	ants, iterations, currentIteration, BestAtIteration int
	distances, pheromones                               [][]float64
	cmsa                                                [][]float64
	tauMin, tauMax, BestLength                          float64
	BestPath                                            []int
}

func NewACO(alpha, beta, rho, pBest, pCmsa float64, ants, iterations int, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromones := make([][]float64, dimension)

	// Arbitrary high value: `4.3. Pheromone trail initialization`
	for i := range pheromones {
		pheromones[i] = make([]float64, dimension)
		for j := range pheromones[i] {
			pheromones[i][j] = 1.0
		}
	}

	return &ACO{
		alpha:      alpha,
		beta:       beta,
		rho:        rho,
		pDec:       math.Pow(pBest, 1.0/float64(dimension)),
		pCmsa:      pCmsa,
		ants:       ants,
		iterations: iterations,
		distances:  distances,
		cmsa:       cmsa,
		pheromones: pheromones,
		BestLength: math.Inf(1),
	}
}

// Main loop to run MMAS
func (aco *ACO) Run() {
	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {
		paths := make([][]int, aco.ants)
		lengths := make([]float64, aco.ants)

		var wg sync.WaitGroup
		wg.Add(aco.ants)
		for i := 0; i < aco.ants; i++ {
			go func(i int) {
				paths[i], lengths[i] = aco.constructPath(i)
				wg.Done()
			}(i)
		}
		wg.Wait()

		// Find the best path in this iteration
		iterationBestLength := math.Inf(1)
		iterationBestPath := []int{}
		for i := 0; i < aco.ants; i++ {
			path := paths[i]
			length := lengths[i]

			if length < iterationBestLength {
				iterationBestLength = length
				iterationBestPath = append([]int(nil), path...)
			}

			if length < aco.BestLength {
				aco.BestLength = length
				aco.BestPath = append([]int(nil), path...)
				aco.BestAtIteration = aco.currentIteration
			}
		}

		aco.globalPheromoneUpdate(iterationBestPath, iterationBestLength)
		aco.updateLimits()
		aco.clampPheromoneLevels()
	}
}

// Function to construct path for each ant
func (aco *ACO) constructPath(antNumber int) ([]int, float64) {
	dimension := len(aco.distances)
	path := make([]int, dimension)
	visited := make([]bool, dimension)
	current := antNumber % dimension
	path[0] = current
	visited[current] = true

	for i := 1; i < dimension; i++ {
		next := aco.selectNextCity(current, visited)

		path[i] = next
		visited[next] = true

		current = next
	}

	length := aco.pathLength(path)
	return path, length
}

// Function to select the next city for an ant
func (aco *ACO) selectNextCity(current int, visited []bool) int {
	dimension := len(aco.distances)
	probabilities := make([]float64, dimension)
	total := 0.0

	// Apply CMSA logic to bias towards better paths. Other paths will also take part in roulette-wheel selection.
	q := rand.Float64()
	adaptiveCmsaProbability := aco.pCmsa * (1.0 - float64(aco.currentIteration)/float64(aco.iterations))
	if q < adaptiveCmsaProbability {
		for i := 0; i < dimension; i++ {
			if !visited[i] && aco.cmsa[current][i] > 0 {
				probabilities[i] = aco.cmsa[current][i]
				total += probabilities[i]
			}
		}
	}

	for i := 0; i < dimension; i++ {
		if !visited[i] && probabilities[i] == 0.0 {
			pheromone := utilities.FastPow(aco.pheromones[current][i], aco.alpha)

			// Adding 1 to each distance in calculation to avoid division by 0.
			heuristic := 1.0 / (aco.distances[current][i] + 1.0)
			desirability := utilities.FastPow(heuristic, aco.beta)
			probabilities[i] = pheromone * desirability
			total += probabilities[i]
		}
	}

	r := rand.Float64()
	cumulativeProbability := 0.0
	for i := 0; i < dimension; i++ {
		if !visited[i] && probabilities[i] > 0.0 {
			cumulativeProbability += probabilities[i] / total

			if r < cumulativeProbability {
				return i
			}
		}
	}

	return -1 // In case no city is selected. This will never happen.
}

func (aco *ACO) updateLimits() {
	aco.tauMax = 1.0 / ((1 - aco.rho) * aco.BestLength)

	pDec := aco.pDec
	nEffective := float64(len(aco.distances)) / 2.0 // Average possible choices.

	numerator := aco.tauMax * (1.0 - pDec)
	denominator := (nEffective - 1.0) * pDec
	aco.tauMin = numerator / denominator
}

func (aco *ACO) clampPheromoneLevels() {
	for i := range aco.pheromones {
		for j := range aco.pheromones[i] {
			if aco.pheromones[i][j] > aco.tauMax {
				aco.pheromones[i][j] = aco.tauMax
			} else if aco.pheromones[i][j] < aco.tauMin {
				aco.pheromones[i][j] = aco.tauMin
			}
		}
	}
}

// Global pheromones update (best ant)
func (aco *ACO) globalPheromoneUpdate(iterationBestPath []int, iterationBestLength float64) {
	// Evaporate pheromones globally
	for i := range aco.pheromones {
		for j := range aco.pheromones[i] {
			aco.pheromones[i][j] *= (1 - aco.rho)
		}
	}

	pheromoneDeposit := 1.0 / iterationBestLength
	bestPath := iterationBestPath

	// https://sci-hub.se/https://doi.org/10.1016/S0167-739X(00)00043-1
	if aco.currentIteration%10 == 0 {
		pheromoneDeposit = 1.0 / aco.BestLength
		bestPath = aco.BestPath
	}

	// Global update: Only the best path deposits pheromones
	for j := 0; j < len(bestPath)-1; j++ {
		start, end := bestPath[j], bestPath[j+1]
		aco.pheromones[start][end] += aco.rho * pheromoneDeposit
	}

	// Handle the wrap-around from the last to the first node
	last, first := bestPath[len(bestPath)-1], bestPath[0]
	aco.pheromones[last][first] += aco.rho * pheromoneDeposit
}

// Function to calculate the length of a path
func (aco *ACO) pathLength(path []int) float64 {
	sum := 0.0
	p := len(path)

	for i := 0; i < p-1; i++ {
		start, end := path[i], path[i+1]
		sum += aco.distances[start][end]
	}

	if p > 0 {
		last, first := path[p-1], path[0]
		sum += aco.distances[last][first]
	}

	return sum
}
