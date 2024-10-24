package aco

import (
	"atsp_aco_msa/modules/algorithms/threeOpt"
	"atsp_aco_msa/modules/utilities"
	"math"
	"math/rand"
)

type ACO struct {
	useLocalSearch                                      bool
	alpha, beta, rho, pDec, pCmsa                       float64
	ants, iterations, currentIteration, BestAtIteration int
	distances, pheromones                               [][]float64
	cmsa                                                [][]float64
	tauMin, tauMax, BestLength                          float64
	BestTour                                            []int
	reducedThreeOpt                                     *threeOpt.ReducedThreeOpt
	knownOptimal                                        float64
	DeviationPerIteration                               []float64
}

func NewACO(useLocalSearch bool, alpha, beta, rho, pBest, pherCmsa, pCmsa float64, ants, iterations int, knownOptimal float64, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromones := make([][]float64, dimension)

	// Arbitrary high value: `4.3. Pheromone trail initialization`
	for i := range pheromones {
		pheromones[i] = make([]float64, dimension)
		for j := range pheromones[i] {
			pheromones[i][j] = 1.0 + (pherCmsa * cmsa[i][j])
		}
	}

	reducedThreeOpt := threeOpt.NewReducedThreeOpt(distances, 25)

	return &ACO{
		useLocalSearch:        useLocalSearch,
		alpha:                 alpha,
		beta:                  beta,
		rho:                   rho,
		pDec:                  math.Pow(pBest, 1.0/float64(dimension)),
		pCmsa:                 pCmsa,
		ants:                  ants,
		iterations:            iterations,
		distances:             distances,
		cmsa:                  cmsa,
		pheromones:            pheromones,
		BestLength:            math.Inf(1),
		reducedThreeOpt:       reducedThreeOpt,
		knownOptimal:          knownOptimal,
		DeviationPerIteration: make([]float64, iterations),
	}
}

// Main loop to run MMAS
func (aco *ACO) Run() {

	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {
		tours := make([][]int, aco.ants)
		lengths := make([]float64, aco.ants)

		for i := 0; i < aco.ants; i++ {
			tours[i], lengths[i] = aco.constructTour(i)
		}

		iterationBestLength := math.Inf(1)
		iterationBestTour := []int{}
		for i := 0; i < aco.ants; i++ {
			tour := tours[i]
			length := lengths[i]

			if length < iterationBestLength {
				iterationBestLength = length
				iterationBestTour = append([]int(nil), tour...)
			}

			if length < aco.BestLength {
				aco.BestLength = length
				aco.BestTour = append([]int(nil), tour...)
				aco.BestAtIteration = aco.currentIteration
			}
		}

		aco.DeviationPerIteration[aco.currentIteration] = 100 * (iterationBestLength - aco.knownOptimal) / aco.knownOptimal

		aco.globalPheromoneUpdate(iterationBestTour, iterationBestLength)
		aco.updateLimits()
		aco.clampPheromoneLevels()
	}
}

// Function to construct tour for each ant
func (aco *ACO) constructTour(antNumber int) ([]int, float64) {
	dimension := len(aco.distances)
	tour := make([]int, dimension)
	visited := make([]bool, dimension)
	current := antNumber % dimension
	tour[0] = current
	visited[current] = true

	for i := 1; i < dimension; i++ {
		next := aco.selectNextCity(current, visited)

		tour[i] = next
		visited[next] = true

		current = next
	}

	if aco.useLocalSearch {
		aco.reducedThreeOpt.Run(tour)
	}

	length := utilities.TourLength(tour, aco.distances)
	return tour, length
}

// Function to select the next city for an ant
func (aco *ACO) selectNextCity(current int, visited []bool) int {
	dimension := len(aco.distances)
	probabilities := make([]float64, dimension)
	total := 0.0

	// Apply CMSA logic to bias towards better tours. Other tours will also take part in roulette-wheel selection.
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
func (aco *ACO) globalPheromoneUpdate(iterationBestTour []int, iterationBestLength float64) {
	// Evaporate pheromones globally
	for i := range aco.pheromones {
		for j := range aco.pheromones[i] {
			aco.pheromones[i][j] *= (1 - aco.rho)
		}
	}

	var pheromoneDeposit float64
	var bestTour []int

	if aco.currentIteration <= 25 {
		pheromoneDeposit = 1.0 / iterationBestLength
		bestTour = iterationBestTour
	} else {
		var f_gb int

		if aco.currentIteration > 25 && aco.currentIteration <= 75 {
			f_gb = 5
		}

		if aco.currentIteration > 75 && aco.currentIteration <= 125 {
			f_gb = 3
		}

		if aco.currentIteration > 125 && aco.currentIteration <= 250 {
			f_gb = 2
		}

		if aco.currentIteration > 250 {
			f_gb = 1
		}

		if aco.currentIteration%f_gb == 0 {
			pheromoneDeposit = 1.0 / aco.BestLength
			bestTour = aco.BestTour
		} else {
			pheromoneDeposit = 1.0 / iterationBestLength
			bestTour = iterationBestTour
		}
	}

	// Global update: Only the best tour deposits pheromones
	for j := 0; j < len(bestTour)-1; j++ {
		start, end := bestTour[j], bestTour[j+1]
		aco.pheromones[start][end] += aco.rho * pheromoneDeposit
	}

	// Handle the wrap-around from the last to the first node
	last, first := bestTour[len(bestTour)-1], bestTour[0]
	aco.pheromones[last][first] += aco.rho * pheromoneDeposit
}
