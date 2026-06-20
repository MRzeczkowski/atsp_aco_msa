package aco

import (
	"atsp_aco_msa/modules/algorithms/neighbors"
	"atsp_aco_msa/modules/algorithms/threeopt"
	"atsp_aco_msa/modules/utilities"
	"math"

	"math/rand/v2"
)

type ACO struct {
	alpha, beta, rho                        float64
	iterations, currentIteration, dimension int
	distances                               [][]float64
	tauMin, tauMax, BestLength              float64
	reducedThreeOpt                         *threeopt.Reduced
	targetTourLength                        float64
	neighborsLists                          [][]int
	heuristicNeighborsLists                 [][]int
	rootedHeuristicNeighborsLists           [][][]int
	heuristicWeight                         float64
	rng                                     *rand.Rand
	useThreeOpt                             bool

	pheromones, heuristics, desirabilities [][]float64

	BestAtIteration           int
	BestTour                  []int
	DeviationPerIteration     []float64
	ThreeOptImprovementsCount int
}

func NewACO(alpha, beta, rho float64, iterations int, targetTourLength float64, distances, heuristicModifiers [][]float64, heuristicWeight float64) *ACO {
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
	localSearchNeighborsLists := neighbors.BuildLists(distances, localSearchNeighborsListSize)

	// We use smaller lists for tour construction than for local search. Just like: https://sci-hub.se/https://doi.org/10.1016/S0167-739X(00)00043-1
	tourConstructionNeighborsListSize := min(maxLocalSearchNeighborsListSize/2, dimension-1)
	tourConstructionNeighborsLists := neighbors.BuildLists(distances, tourConstructionNeighborsListSize)
	heuristicNeighborsLists := buildConstructionHeuristicNeighborsLists(heuristicModifiers)

	reducedThreeOpt := threeopt.NewReduced(distances, localSearchNeighborsLists)

	src := rand.NewPCG(4, 2)
	rng := rand.New(src)

	return &ACO{
		alpha:                   alpha,
		beta:                    beta,
		rho:                     rho,
		iterations:              iterations,
		dimension:               dimension,
		distances:               distances,
		tauMin:                  -math.MaxInt,
		tauMax:                  math.MaxInt,
		reducedThreeOpt:         reducedThreeOpt,
		targetTourLength:        targetTourLength,
		neighborsLists:          tourConstructionNeighborsLists,
		heuristicNeighborsLists: heuristicNeighborsLists,
		heuristicWeight:         heuristicWeight,
		rng:                     rng,

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

func (aco *ACO) SetUseThreeOpt(enabled bool) {
	aco.useThreeOpt = enabled
}

func (aco *ACO) SetRootedHeuristicModifiers(rootedHeuristicModifiers [][][]float64) {
	aco.rootedHeuristicNeighborsLists = make([][][]int, len(rootedHeuristicModifiers))
	for root := range rootedHeuristicModifiers {
		aco.rootedHeuristicNeighborsLists[root] = buildConstructionHeuristicNeighborsLists(rootedHeuristicModifiers[root])
	}
}

func buildConstructionHeuristicNeighborsLists(heuristicModifiers [][]float64) [][]int {
	neighborsLists := make([][]int, len(heuristicModifiers))

	for i := range heuristicModifiers {
		for j, modifier := range heuristicModifiers[i] {
			if i == j || modifier <= 0.0 {
				continue
			}

			neighborsLists[i] = append(neighborsLists[i], j)
		}
	}

	return neighborsLists
}

func (aco *ACO) Run() {
	ants := min(25, aco.dimension)
	initialPheromoneValue := 100000.0 // Arbitrary large value
	zeroDistanceHeuristicBase := 1000000.0

	for i := 0; i < aco.dimension; i++ {
		for j := 0; j < aco.dimension; j++ {
			if i == j {
				continue
			}

			distance := aco.distances[i][j]
			heuristicBase := zeroDistanceHeuristicBase
			if distance != 0.0 {
				heuristicBase = 1.0 / distance
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
	aco.reducedThreeOpt.Improvements = 0
	defer func() {
		aco.ThreeOptImprovementsCount = aco.reducedThreeOpt.Improvements
	}()

	tours := make([][]int, ants)

	// This is used to implement "branchless programming" i.e. get rid of if statements in tight loops.
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

			if aco.useThreeOpt {
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
				}
			}
		}

		currentDeviation := 100 * (aco.BestLength - aco.targetTourLength) / aco.targetTourLength
		aco.DeviationPerIteration[aco.currentIteration] = currentDeviation

		if currentDeviation == 0.0 {
			return
		}

		aco.updateLimits()
		aco.evaporatePheromones()

		useGlobal := useGlobalBestPheromoneUpdate(aco.currentIteration, aco.iterations)

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
}

func useGlobalBestPheromoneUpdate(currentIteration, iterations int) bool {
	if iterations <= 0 {
		return false
	}

	progress := float64(currentIteration) / float64(iterations)
	switch {
	case progress < 0.05:
		return false
	case progress < 0.30:
		return currentIteration%5 == 0
	case progress < 0.70:
		return currentIteration%3 == 0
	case progress < 0.95:
		return currentIteration%2 == 0
	default:
		return true
	}
}

// Function to construct tour for each ant
func (aco *ACO) constructTour(tour []int, canVisitBits []float64, desirabilities []float64) {
	start := aco.rng.IntN(aco.dimension)
	current := start
	tour[0] = start

	for i := 0; i < aco.dimension; i++ {
		canVisitBits[i] = 1.0
	}
	canVisitBits[current] = 0.0

	for i := 1; i < aco.dimension; i++ {
		next := aco.selectNextCityForStart(start, current, canVisitBits, desirabilities)

		tour[i] = next
		canVisitBits[next] = 0.0

		current = next
	}
}

// Function to select the next city for an ant
func (aco *ACO) selectNextCity(current int, canVisitBits []float64, desirabilities []float64) int {
	return aco.selectNextCityForStart(-1, current, canVisitBits, desirabilities)
}

func (aco *ACO) selectNextCityForStart(start, current int, canVisitBits []float64, desirabilities []float64) int {
	if aco.shouldUseConstructionHeuristic() {
		if nextCity := aco.selectNextCityFromConstructionHeuristic(start, current, canVisitBits, desirabilities); nextCity != -1 {
			return nextCity
		}
	}

	desirabilitiesToUse := aco.desirabilities[current]

	total := 0.0
	neighbors := aco.neighborsLists[current]
	for _, i := range neighbors {
		prob := canVisitBits[i] * desirabilitiesToUse[i]
		desirabilities[i] = prob

		total += prob
	}

	cumulativeProbability := 0.0
	nextCity := -1

	// Check if total is zero (no valid neighbors left), then fallback to all remaining nodes
	if total != 0.0 {
		q := aco.rng.Float64()
		threshold := q * total

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
		bestDesirability := -math.MaxFloat64
		for i, probability := range desirabilitiesToUse {
			desirabilities[i] = 0.0

			if canVisitBits[i] == 0.0 {
				continue
			}

			if probability > bestDesirability {
				bestDesirability = probability
				nextCity = i
			}
		}
	}

	return nextCity
}

func (aco *ACO) shouldUseConstructionHeuristic() bool {
	probability := aco.constructionHeuristicProbability()
	return probability > 0.0 && aco.rng.Float64() < probability
}

func (aco *ACO) constructionHeuristicProbability() float64 {
	if aco.heuristicWeight <= 0.0 || aco.iterations <= 0 {
		return 0.0
	}

	progress := float64(aco.currentIteration) / float64(aco.iterations)
	if progress >= 1.0 {
		return 0.0
	}

	return min(aco.heuristicWeight, 1.0) * (1.0 - progress)
}

func (aco *ACO) selectNextCityFromConstructionHeuristic(start, current int, canVisitBits []float64, desirabilities []float64) int {
	desirabilitiesToUse := aco.desirabilities[current]
	neighbors := aco.constructionHeuristicNeighbors(start, current)

	total := 0.0
	for _, i := range neighbors {
		prob := canVisitBits[i] * desirabilitiesToUse[i]
		desirabilities[i] = prob
		total += prob
	}

	if total == 0.0 {
		for _, i := range neighbors {
			desirabilities[i] = 0.0
		}
		return -1
	}

	threshold := aco.rng.Float64() * total
	cumulativeProbability := 0.0
	nextCity := -1

	for _, i := range neighbors {
		if nextCity == -1 && desirabilities[i] > 0.0 {
			cumulativeProbability += desirabilities[i]
			if threshold < cumulativeProbability {
				nextCity = i
			}
		}

		desirabilities[i] = 0.0
	}

	return nextCity
}

func (aco *ACO) constructionHeuristicNeighbors(start, current int) []int {
	if start >= 0 && start < len(aco.rootedHeuristicNeighborsLists) {
		rootedNeighbors := aco.rootedHeuristicNeighborsLists[start]
		if current >= 0 && current < len(rootedNeighbors) {
			return rootedNeighbors[current]
		}
	}

	if current >= 0 && current < len(aco.heuristicNeighborsLists) {
		return aco.heuristicNeighborsLists[current]
	}

	return nil
}

func (aco *ACO) updateLimits() {
	aco.tauMax = 1.0 / ((1.0 - aco.rho) * aco.BestLength)
	aco.tauMin = aco.tauMax / (2.0 * float64(aco.dimension))
}

func (aco *ACO) setPheromone(i, j int, newValue float64) {
	if i == j {
		return
	}

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
