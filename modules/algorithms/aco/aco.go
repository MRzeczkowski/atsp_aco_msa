package aco

import (
	"atsp_aco_msa/modules/algorithms/nearestNeighbors"
	"atsp_aco_msa/modules/algorithms/threeOpt"
	"atsp_aco_msa/modules/utilities"
	"math"

	"math/rand/v2"
)

type ACO struct {
	alpha, beta, rho, pCmsa                 float64
	iterations, currentIteration, dimension int
	distances                               [][]float64
	hints                                   [][]float64
	tauMin, tauMax, BestLength              float64
	reducedThreeOpt                         *threeOpt.ReducedThreeOpt
	targetTourLength                        float64
	neighborsLists                          [][]int
	rng                                     *rand.Rand

	pheromones, heuristics, desirabilities [][]float64

	BestAtIteration           int
	BestTour                  []int
	DeviationPerIteration     []float64
	ThreeOptImprovementsCount int
}

func NewACO(alpha, beta, rho, pCmsa float64, iterations int, targetTourLength float64, distances, hints [][]float64) *ACO {
	dimension := len(distances)
	pheromones := make([][]float64, dimension)
	heuristics := make([][]float64, dimension)
	desirabilities := make([][]float64, dimension)

	for i := 0; i < dimension; i++ {
		pheromones[i] = make([]float64, dimension)
		heuristics[i] = make([]float64, dimension)
		desirabilities[i] = make([]float64, dimension)
	}

	maxLocalSearchNeighborsListSize := 20
	localSearchNeighborsListSize := min(maxLocalSearchNeighborsListSize, dimension-1)
	localSearchNeighborsLists := nearestNeighbors.BuildNearestNeighborsLists(distances, localSearchNeighborsListSize)

	// We use smaller lists for tour construction than for local search. Just like: https://sci-hub.se/https://doi.org/10.1016/S0167-739X(00)00043-1
	tourConstructionNeighborsListSize := min(maxLocalSearchNeighborsListSize/2, dimension-1)
	tourConstructionNeighborsLists := make([][]int, dimension)
	for i := 0; i < dimension; i++ {
		tourConstructionNeighborsLists[i] = make([]int, tourConstructionNeighborsListSize)

		for j := 0; j < tourConstructionNeighborsListSize; j++ {
			tourConstructionNeighborsLists[i][j] = localSearchNeighborsLists[i][j]
		}
	}

	reducedThreeOpt := threeOpt.NewReducedThreeOpt(distances, localSearchNeighborsLists)

	src := rand.NewPCG(4, 2)
	rng := rand.New(src)

	return &ACO{
		alpha:            alpha,
		beta:             beta,
		rho:              rho,
		pCmsa:            pCmsa,
		iterations:       iterations,
		dimension:        dimension,
		distances:        distances,
		hints:            hints,
		tauMin:           -math.MaxInt,
		tauMax:           math.MaxInt,
		reducedThreeOpt:  reducedThreeOpt,
		targetTourLength: targetTourLength,
		neighborsLists:   tourConstructionNeighborsLists,
		rng:              rng,

		pheromones:     pheromones,
		heuristics:     heuristics,
		desirabilities: desirabilities,

		BestAtIteration:           math.MaxInt,
		BestLength:                math.MaxFloat64,
		BestTour:                  make([]int, dimension),
		DeviationPerIteration:     make([]float64, iterations),
		ThreeOptImprovementsCount: 0,
	}
}

func (aco *ACO) Run() {

	hintsSums := make([]float64, aco.dimension)

	for i := 0; i < aco.dimension; i++ {
		hintsSums[i] = 0

		for j := 0; j < aco.dimension; j++ {
			hintsSums[i] += aco.hints[i][j]
		}
	}

	ants := min(25, aco.dimension)
	initialPheromoneValue := 100000.0 // Arbitrary large value
	for i := 0; i < aco.dimension; i++ {
		for j := 0; j < aco.dimension; j++ {

			// Adding 1 to each distance in calculation to avoid division by 0.
			heuristicBase := 1.0 / (aco.distances[i][j] + 1.0)

			if aco.hints[i][j] != 0 {
				heuristicBase *= 1 + ((aco.hints[i][j] / hintsSums[i]) * aco.pCmsa)
			}

			heuristic := math.Pow(heuristicBase, aco.beta)
			aco.heuristics[i][j] = heuristic

			aco.setPheromone(i, j, initialPheromoneValue)
		}
	}

	aco.BestAtIteration = math.MaxInt
	aco.BestLength = math.MaxFloat64
	aco.BestTour = make([]int, aco.dimension)
	aco.DeviationPerIteration = make([]float64, aco.iterations)
	aco.ThreeOptImprovementsCount = 0

	tours := make([][]int, ants)
	canVisitBits := make([][]float64, ants)
	desirabilities := make([][]float64, ants)

	for i := 0; i < ants; i++ {
		tours[i] = make([]int, aco.dimension)
		canVisitBits[i] = make([]float64, aco.dimension)
		desirabilities[i] = make([]float64, aco.dimension)
	}

	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {

		iterationBestLength := math.MaxFloat64
		iterationBestTour := make([]int, aco.dimension)

		for i := 0; i < ants; i++ {
			aco.constructTour(tours[i], canVisitBits[i], desirabilities[i])

			aco.reducedThreeOpt.Run(tours[i])

			length := utilities.TourLength(tours[i], aco.distances)

			if length < iterationBestLength {
				iterationBestLength = length
				copy(iterationBestTour, tours[i])

				if iterationBestLength < aco.BestLength {
					aco.BestLength = iterationBestLength
					copy(aco.BestTour, iterationBestTour)
					aco.BestAtIteration = aco.currentIteration
				}

				currentDeviation := 100 * (iterationBestLength - aco.targetTourLength) / aco.targetTourLength
				aco.DeviationPerIteration[aco.currentIteration] = currentDeviation

				if currentDeviation == 0.0 {
					return
				}
			}
		}

		aco.updateLimits()
		aco.evaporatePheromones()

		var useGlobal bool

		switch {
		case aco.currentIteration < 25:
			useGlobal = false
		case aco.currentIteration < 75:
			useGlobal = aco.currentIteration%5 == 0
		case aco.currentIteration < 125:
			useGlobal = aco.currentIteration%3 == 0
		case aco.currentIteration < 250:
			useGlobal = aco.currentIteration%2 == 0
		default:
			useGlobal = true
		}

		var depositTour []int
		var pheromoneDeposit float64

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
func (aco *ACO) constructTour(tour []int, canVisitBits []float64, desirabilities []float64) {
	current := aco.rng.IntN(aco.dimension)
	tour[0] = current

	for i := 0; i < aco.dimension; i++ {
		canVisitBits[i] = 1.0
	}
	canVisitBits[current] = 0.0

	for i := 1; i < aco.dimension; i++ {
		next := aco.selectNextCity(current, canVisitBits, desirabilities)

		tour[i] = next
		canVisitBits[next] = 0.0

		current = next
	}
}

// Function to select the next city for an ant
func (aco *ACO) selectNextCity(current int, canVisitBits []float64, desirabilities []float64) int {
	desirabilitiesToUse := aco.desirabilities[current]

	total := 0.0
	neighbors := aco.neighborsLists[current]
	for _, i := range neighbors {
		prob := canVisitBits[i] * desirabilitiesToUse[i]
		desirabilities[i] = prob

		total += prob
	}

	q := aco.rng.Float64()
	threshold := q * total
	cumulativeProbability := 0.0
	nextCity := -1

	// Check if total is zero (no valid neighbors left), then fallback to all remaining nodes
	if total != 0.0 {
		for _, i := range neighbors {
			if nextCity == -1 && desirabilities[i] > 0.0 {
				cumulativeProbability += desirabilities[i]
				if threshold < cumulativeProbability {
					nextCity = i
				}
			}

			desirabilities[i] = 0.0
		}
	} else {
		for i, probability := range desirabilitiesToUse {
			prob := canVisitBits[i] * probability
			desirabilities[i] = prob

			total += prob
		}

		threshold = q * total

		for i, probability := range desirabilities {
			if nextCity == -1 && probability > 0.0 {
				cumulativeProbability += probability
				if threshold < cumulativeProbability {
					nextCity = i
				}
			}

			desirabilities[i] = 0.0
		}
	}

	return nextCity
}

func (aco *ACO) updateLimits() {
	aco.tauMax = 1.0 / ((1.0 - aco.rho) * aco.BestLength)
	aco.tauMin = aco.tauMax / (2.0 * float64(aco.dimension))
}

func (aco *ACO) setPheromone(i, j int, newValue float64) {
	newValue = aco.clampPheromoneLevel(newValue)

	oldValue := aco.pheromones[i][j]
	if oldValue == newValue {
		return
	}

	aco.pheromones[i][j] = newValue
	pheromone := math.Pow(aco.pheromones[i][j], aco.alpha)
	aco.desirabilities[i][j] = pheromone * aco.heuristics[i][j]
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
	for i := 0; i < aco.dimension; i++ {
		for j := 0; j < aco.dimension; j++ {
			evaporatedValue := aco.pheromones[i][j] * aco.rho
			aco.setPheromone(i, j, evaporatedValue)
		}
	}
}

func (aco *ACO) depositPheromones(tour []int, pheromoneDeposit float64) {
	for i := 0; i < aco.dimension-1; i++ {
		start, end := tour[i], tour[i+1]

		newValue := aco.pheromones[start][end] + pheromoneDeposit
		aco.setPheromone(start, end, newValue)
	}

	// Handle the wrap-around from the last to the first node
	last, first := tour[aco.dimension-1], tour[0]

	newValue := aco.pheromones[last][first] + pheromoneDeposit
	aco.setPheromone(last, first, newValue)
}
