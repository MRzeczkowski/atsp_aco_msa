package aco

import (
	"atsp_aco_msa/modules/algorithms/threeOpt"
	"atsp_aco_msa/modules/utilities"
	"math"
	"math/rand"
)

type ACO struct {
	useLocalSearch                                                                 bool
	alpha, beta, rho, pDec, pCmsa                                                  float64
	ants, iterations, currentIteration, BestAtIteration, dimension                 int
	distances, pheromones, desirabilitiesPreCalc, probabilities, cmsaProbabilities [][]float64
	cmsa                                                                           [][]float64
	tauMin, tauMax, BestLength                                                     float64
	BestTour                                                                       []int
	reducedThreeOpt                                                                *threeOpt.ReducedThreeOpt
	knownOptimal                                                                   float64
	DeviationPerIteration                                                          []float64
}

func NewACO(useLocalSearch bool, alpha, beta, rho, pBest, pherCmsa, pCmsa float64, ants, iterations int, knownOptimal float64, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromones := make([][]float64, dimension)
	desirabilitiesPreCalc := make([][]float64, dimension)
	probabilities := make([][]float64, dimension)
	cmsaProbabilities := make([][]float64, dimension)

	// Arbitrary high value: `4.3. Pheromone trail initialization`
	for i := range pheromones {

		pheromones[i] = make([]float64, dimension)
		desirabilitiesPreCalc[i] = make([]float64, dimension)
		probabilities[i] = make([]float64, dimension)
		cmsaProbabilities[i] = make([]float64, dimension)

		for j := range pheromones[i] {
			pheromones[i][j] = 1.0 + (pherCmsa * cmsa[i][j])

			pheromone := math.Pow(pheromones[i][j], alpha)

			// Adding 1 to each distance in calculation to avoid division by 0.
			heuristic := 1.0 / (distances[i][j] + 1.0)
			desirability := math.Pow(heuristic, beta)

			desirabilitiesPreCalc[i][j] = desirability

			probabilities[i][j] = pheromone * desirability
			cmsaProbabilities[i][j] = probabilities[i][j] + cmsa[i][j]
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
		dimension:             dimension,
		distances:             distances,
		cmsa:                  cmsa,
		pheromones:            pheromones,
		desirabilitiesPreCalc: desirabilitiesPreCalc,
		probabilities:         probabilities,
		cmsaProbabilities:     cmsaProbabilities,
		BestLength:            math.Inf(1),
		reducedThreeOpt:       reducedThreeOpt,
		knownOptimal:          knownOptimal,
		DeviationPerIteration: make([]float64, iterations),
	}
}

// Main loop to run MMAS
func (aco *ACO) Run() {

	tours := make([][]int, aco.ants)
	canVisitBits := make([][]float64, aco.ants)
	probabilities := make([][]float64, aco.ants)

	for i := 0; i < aco.ants; i++ {
		tours[i] = make([]int, aco.dimension)
		canVisitBits[i] = make([]float64, aco.dimension)
		probabilities[i] = make([]float64, aco.dimension)
	}

	lengths := make([]float64, aco.ants)

	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {

		for i := 0; i < aco.ants; i++ {
			aco.constructTour(i, tours[i], canVisitBits[i], probabilities[i])
			lengths[i] = utilities.TourLength(tours[i], aco.distances)
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
func (aco *ACO) constructTour(antNumber int, tour []int, canVisitBits []float64, probabilities []float64) {
	current := antNumber % aco.dimension
	tour[0] = current

	for i := 0; i < aco.dimension; i++ {
		canVisitBits[i] = 1.0
	}
	canVisitBits[current] = 0.0

	for i := 1; i < aco.dimension; i++ {
		next := aco.selectNextCity(current, canVisitBits, probabilities)

		tour[i] = next
		canVisitBits[next] = 0.0

		current = next
	}

	if aco.useLocalSearch {
		aco.reducedThreeOpt.Run(tour)
	}
}

// Function to select the next city for an ant
func (aco *ACO) selectNextCity(current int, canVisitBits []float64, probabilities []float64) int {
	total := 0.0

	// Apply CMSA logic to bias towards hopefully better tours. Other tours will also take part in roulette-wheel selection.
	q := rand.Float64()
	adaptiveCmsaProbability := aco.pCmsa * (1.0 - float64(aco.currentIteration)/float64(aco.iterations))

	if q < adaptiveCmsaProbability {
		for i := 0; i < aco.dimension; i++ {
			probabilities[i] = canVisitBits[i] * aco.cmsaProbabilities[current][i]

			total += probabilities[i]
		}
	} else {
		for i := 0; i < aco.dimension; i++ {
			probabilities[i] = canVisitBits[i] * aco.probabilities[current][i]

			total += probabilities[i]
		}
	}

	cumulativeProbability := 0.0
	nextCity := -1
	for i := 0; i < aco.dimension; i++ {
		cumulativeProbability += probabilities[i] / total
		probabilities[i] = 0.0

		if q < cumulativeProbability && nextCity == -1 {
			nextCity = i
		}
	}

	return nextCity
}

func (aco *ACO) updateLimits() {
	aco.tauMax = 1.0 / ((1 - aco.rho) * aco.BestLength)

	pDec := aco.pDec
	nEffective := float64(aco.dimension) / 2.0 // Average possible choices.

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

	evaporationCoefficient := 1 - aco.rho
	// Evaporate pheromones globally
	for i := range aco.pheromones {
		for j := range aco.pheromones[i] {
			aco.pheromones[i][j] *= evaporationCoefficient

			pheromone := math.Pow(aco.pheromones[i][j], aco.alpha)
			aco.probabilities[i][j] = pheromone * aco.desirabilitiesPreCalc[i][j]
			aco.cmsaProbabilities[i][j] = aco.probabilities[i][j] + aco.cmsa[i][j]
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
	for j := 0; j < aco.dimension-1; j++ {
		start, end := bestTour[j], bestTour[j+1]
		aco.pheromones[start][end] += aco.rho * pheromoneDeposit

		pheromone := math.Pow(aco.pheromones[start][end], aco.alpha)
		aco.probabilities[start][end] = pheromone * aco.desirabilitiesPreCalc[start][end]
		aco.cmsaProbabilities[start][end] = aco.probabilities[start][end] + aco.cmsa[start][end]
	}

	// Handle the wrap-around from the last to the first node
	last, first := bestTour[aco.dimension-1], bestTour[0]
	aco.pheromones[last][first] += aco.rho * pheromoneDeposit

	pheromone := math.Pow(aco.pheromones[last][first], aco.alpha)
	aco.probabilities[last][first] = pheromone * aco.desirabilitiesPreCalc[last][first]
	aco.cmsaProbabilities[last][first] = aco.probabilities[last][first] + aco.cmsa[last][first]
}
