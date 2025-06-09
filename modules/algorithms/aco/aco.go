package aco

import (
	"atsp_aco_msa/modules/algorithms/nearestNeighbors"
	"atsp_aco_msa/modules/algorithms/threeOpt"
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"math"

	"pgregory.net/rand"
)

type ACO struct {
	alpha, beta, rho, pDec, pCmsa                                  float64
	ants, localSearchAnts, iterations, currentIteration, dimension int
	distances                                                      [][]float64
	cmsa                                                           [][]float64
	tauMin, tauMax, BestLength                                     float64
	reducedThreeOpt                                                *threeOpt.ReducedThreeOpt
	targetTourLength                                               float64
	neighborsLists                                                 [][]int

	pheromones, desirabilitiesPreCalc, probabilities, cmsaProbabilities [][]float64

	BestAtIteration           int
	BestTour                  []int
	DeviationPerIteration     []float64
	ThreeOptImprovementsCount int
}

func NewACO(alpha, beta, rho, pBest, pCmsa float64, ants, localSearchAnts, iterations int, targetTourLength float64, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromones := make([][]float64, dimension)
	desirabilitiesPreCalc := make([][]float64, dimension)
	probabilities := make([][]float64, dimension)
	cmsaProbabilities := make([][]float64, dimension)

	for i := range pheromones {
		pheromones[i] = make([]float64, dimension)
		desirabilitiesPreCalc[i] = make([]float64, dimension)
		probabilities[i] = make([]float64, dimension)
		cmsaProbabilities[i] = make([]float64, dimension)
	}

	maxLocalSearchNeighborsListSize := 80
	localSearchNeighborsListSize := min(maxLocalSearchNeighborsListSize, dimension)
	localSearchNeighborsLists := nearestNeighbors.BuildNearestNeighborsLists(distances, localSearchNeighborsListSize)

	// We use smaller lists for tour construction than for local search. Just like: https://sci-hub.se/https://doi.org/10.1016/S0167-739X(00)00043-1
	tourConstructionNeighborsListSize := min(maxLocalSearchNeighborsListSize/2, dimension)
	tourConstructionNeighborsLists := make([][]int, dimension)
	for i := 0; i < dimension; i++ {
		tourConstructionNeighborsLists[i] = make([]int, tourConstructionNeighborsListSize)

		for j := 0; j < tourConstructionNeighborsListSize; j++ {
			tourConstructionNeighborsLists[i][j] = localSearchNeighborsLists[i][j]
		}
	}

	reducedThreeOpt := threeOpt.NewReducedThreeOpt(distances, localSearchNeighborsLists)

	return &ACO{
		alpha:            alpha,
		beta:             beta,
		rho:              rho,
		pDec:             math.Pow(pBest, 1.0/float64(dimension)),
		pCmsa:            pCmsa,
		ants:             ants,
		localSearchAnts:  localSearchAnts,
		iterations:       iterations,
		dimension:        dimension,
		distances:        distances,
		cmsa:             cmsa,
		reducedThreeOpt:  reducedThreeOpt,
		targetTourLength: targetTourLength,
		neighborsLists:   tourConstructionNeighborsLists,

		pheromones:            pheromones,
		desirabilitiesPreCalc: desirabilitiesPreCalc,
		probabilities:         probabilities,
		cmsaProbabilities:     cmsaProbabilities,

		BestAtIteration:           math.MaxInt,
		BestLength:                math.MaxFloat64,
		BestTour:                  make([]int, dimension),
		DeviationPerIteration:     make([]float64, iterations),
		ThreeOptImprovementsCount: 0,
	}
}

// Main loop to run MMAS
func (aco *ACO) Run() {
	for i := range aco.pheromones {
		for j := range aco.pheromones[i] {
			// Adding 1 to each distance in calculation to avoid division by 0.
			heuristic := 1.0 / (aco.distances[i][j] + 1.0)
			desirability := math.Pow(heuristic, aco.beta)
			aco.desirabilitiesPreCalc[i][j] = desirability

			aco.setPheromone(i, j, 100.0)
		}
	}

	aco.BestAtIteration = math.MaxInt
	aco.BestLength = math.MaxFloat64
	aco.BestTour = make([]int, aco.dimension)
	aco.DeviationPerIteration = make([]float64, aco.iterations)
	aco.ThreeOptImprovementsCount = 0

	tours := make([][]int, aco.ants)
	canVisitBits := make([][]float64, aco.ants)
	probabilities := make([][]float64, aco.ants)

	for i := 0; i < aco.ants; i++ {
		tours[i] = make([]int, aco.dimension)
		canVisitBits[i] = make([]float64, aco.dimension)
		probabilities[i] = make([]float64, aco.dimension)
	}

	// Restart mechanism state
	noImproveGlobal := 0

	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {

		iterationBestLength := math.MaxFloat64
		iterationBestTour := make([]int, aco.dimension)

		globalImproved := false
		for i := 0; i < aco.ants; i++ {
			aco.constructTour(tours[i], canVisitBits[i], probabilities[i])

			if i < aco.localSearchAnts {
				aco.reducedThreeOpt.Run(tours[i])
			}

			length := utilities.TourLength(tours[i], aco.distances)

			if length < iterationBestLength {
				iterationBestLength = length
				copy(iterationBestTour, tours[i])

				if iterationBestLength < aco.BestLength {
					aco.BestLength = iterationBestLength
					copy(aco.BestTour, iterationBestTour)
					aco.BestAtIteration = aco.currentIteration
					aco.updateLimits()

					globalImproved = true
				}

				currentDeviation := 100 * (iterationBestLength - aco.targetTourLength) / aco.targetTourLength
				aco.DeviationPerIteration[aco.currentIteration] = currentDeviation

				if currentDeviation == 0.0 {
					fmt.Println("DONE!")
					return
				}
			}
		}

		// Update improvement counters
		if !globalImproved {
			noImproveGlobal++
		} else {
			noImproveGlobal = 0
		}

		// Check restart conditions
		if noImproveGlobal >= 50 {
			bf := aco.calculateBranchingFactor(0.1)
			if bf <= 1.1 {
				for i := range aco.pheromones {
					for j := range aco.pheromones[i] {
						aco.setPheromone(i, j, aco.tauMax)
					}
				}
			}
		}

		aco.evaporatePheromones()

		var depositTour []int
		var pheromoneDeposit float64

		var useGlobal bool
		if aco.currentIteration <= 25 {
			useGlobal = false
		} else {
			var phase int

			switch {
			case aco.currentIteration <= 75:
				phase = 5
			case aco.currentIteration <= 125:
				phase = 3
			case aco.currentIteration <= 250:
				phase = 2
			default:
				phase = 1
			}

			useGlobal = aco.currentIteration%phase == 0
		}

		if useGlobal {
			depositTour = aco.BestTour
			pheromoneDeposit = 1.0 / aco.BestLength
		} else {
			depositTour = iterationBestTour
			pheromoneDeposit = 1.0 / iterationBestLength
		}

		aco.depositPheromones(depositTour, pheromoneDeposit)
	}

	aco.ThreeOptImprovementsCount = aco.reducedThreeOpt.Improvements
}

// Function to construct tour for each ant
func (aco *ACO) constructTour(tour []int, canVisitBits []float64, probabilities []float64) {
	current := rand.Intn(aco.dimension)
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
}

// Function to select the next city for an ant
func (aco *ACO) selectNextCity(current int, canVisitBits []float64, probabilities []float64) int {

	q := rand.Float64()
	probabilitiesToUse := aco.probabilities[current]
	// Use CMSA logic to bias towards hopefully better tours. Other tours will also take part in roulette-wheel selection.
	adaptiveCmsaProbability := aco.pCmsa * (1.0 - float64(aco.currentIteration)/float64(aco.iterations))
	if q < adaptiveCmsaProbability {
		probabilitiesToUse = aco.cmsaProbabilities[current]
	}

	total := 0.0
	neighbors := aco.neighborsLists[current]
	for _, i := range neighbors {
		prob := canVisitBits[i] * probabilitiesToUse[i]
		probabilities[i] = prob

		total += prob
	}

	threshold := q * total
	cumulativeProbability := 0.0
	nextCity := -1

	// Check if total is zero (no valid neighbors left), then fallback to all remaining nodes
	if total != 0.0 {
		for _, i := range neighbors {
			if nextCity == -1 && probabilities[i] > 0.0 {
				cumulativeProbability += probabilities[i]
				if threshold < cumulativeProbability {
					nextCity = i
				}
			}

			probabilities[i] = 0.0
		}
	} else {
		for i, probability := range probabilitiesToUse {
			prob := canVisitBits[i] * probability
			probabilities[i] = prob

			total += prob
		}

		threshold = q * total

		for i, probability := range probabilities {
			if nextCity == -1 && probability > 0.0 {
				cumulativeProbability += probability
				if threshold < cumulativeProbability {
					nextCity = i
				}
			}

			probabilities[i] = 0.0
		}
	}

	return nextCity
}

func (aco *ACO) updateLimits() {
	aco.tauMax = 1.0 / ((1 - aco.rho) * aco.BestLength)

	numerator := aco.tauMax * (1.0 - aco.pDec)

	nEffective := float64(aco.dimension) / 2.0 // Average possible choices.
	denominator := (nEffective - 1.0) * aco.pDec

	aco.tauMin = numerator / denominator

	for i := range aco.pheromones {
		for j := range aco.pheromones[i] {
			aco.setPheromone(i, j, aco.clampPheromoneLevel(aco.pheromones[i][j]))
		}
	}
}

func (aco *ACO) setPheromone(i, j int, value float64) {
	aco.pheromones[i][j] = value
	pheromone := math.Pow(aco.pheromones[i][j], aco.alpha)
	aco.probabilities[i][j] = pheromone * aco.desirabilitiesPreCalc[i][j]
	aco.cmsaProbabilities[i][j] = aco.probabilities[i][j] + aco.cmsa[i][j]
}

func (aco *ACO) calculateBranchingFactor(lambda float64) float64 {
	total := 0.0
	for i := 0; i < aco.dimension; i++ {
		// Get all outgoing edges from node i
		edges := aco.pheromones[i]
		tauMin, tauMax := findMinMaxExcludingIndex(edges, i)

		threshold := tauMin + lambda*(tauMax-tauMin)
		count := 0
		for j, tau := range edges {
			if i != j && tau > threshold {
				count++
			}
		}

		total += float64(count)
	}

	return total / float64(aco.dimension)
}

func findMinMaxExcludingIndex(slice []float64, index int) (float64, float64) {
	min := math.MaxFloat64
	max := -math.MaxFloat64

	for i, value := range slice {
		if i == index {
			continue
		}

		if value < min {
			min = value
		}

		if value > max {
			max = value
		}
	}

	return min, max
}

func (aco *ACO) clampPheromoneLevel(pheromone float64) float64 {
	if pheromone > aco.tauMax {
		return aco.tauMax
	} else if pheromone < aco.tauMin {
		return aco.tauMin
	}

	return pheromone
}

func (aco *ACO) evaporatePheromones() {
	// TODO: Evaporate only edges that are in nearest neighbors lists!
	// This will be faster, same as in MMAS paper and will hopefully increase branching factor.
	// There are some issue with it however, since not all edges are updated the branching factor seems to be useless.
	// for i := 0; i < aco.dimension; i++ {
	// 	for _, j := range aco.neighborsLists[i] {

	// Evaporate pheromones globally
	for i := range aco.pheromones {
		for j := range aco.pheromones[i] {
			evaporatedValue := aco.clampPheromoneLevel(aco.pheromones[i][j] * aco.rho)
			aco.setPheromone(i, j, evaporatedValue)
		}
	}
}

func (aco *ACO) depositPheromones(tour []int, pheromoneDeposit float64) {
	for i := 0; i < aco.dimension-1; i++ {
		start, end := tour[i], tour[i+1]

		newValue := aco.clampPheromoneLevel(aco.pheromones[start][end] + pheromoneDeposit)
		aco.setPheromone(start, end, newValue)
	}

	// Handle the wrap-around from the last to the first node
	last, first := tour[aco.dimension-1], tour[0]

	newValue := aco.clampPheromoneLevel(aco.pheromones[last][first] + pheromoneDeposit)
	aco.setPheromone(last, first, newValue)
}
