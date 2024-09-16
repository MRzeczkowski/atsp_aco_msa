package aco

import (
	"atsp_aco_msa/modules/utilities"
	"math"
	"math/rand"
)

type ACO struct {
	alpha, beta, rho, cmsaP                             float64
	ants, iterations, currentIteration, BestAtIteration int
	distances, pheromone                                [][]float64
	cmsa                                                [][]float64
	tauMin, tauMax, BestLength                          float64
	BestPath                                            []int
}

func NewACO(alpha, beta, rho, cmsaP float64, ants, iterations int, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromone := make([][]float64, dimension)

	// Calculate tauMax using nearest neighbor heuristic
	L_nn := nearestNeighborTourLength(distances)
	tauMax := 1.0 / (rho * L_nn)
	// Set tauMin relative to tauMax
	tauMin := tauMax / (2.0 * float64(dimension))

	for i := range pheromone {
		pheromone[i] = make([]float64, dimension)
		for j := range pheromone[i] {
			pheromone[i][j] = tauMax
		}
	}

	return &ACO{
		alpha:      alpha,
		beta:       beta,
		rho:        rho,
		cmsaP:      cmsaP,
		ants:       ants,
		iterations: iterations,
		distances:  distances,
		cmsa:       cmsa,
		pheromone:  pheromone,
		tauMin:     tauMin,
		tauMax:     tauMax,
		BestLength: math.Inf(1),
	}
}

func nearestNeighborTourLength(distances [][]float64) float64 {
	n := len(distances)
	visited := make([]bool, n)
	currentCity := 0
	visited[currentCity] = true
	totalLength := 0.0

	for i := 1; i < n; i++ {
		nextCity := -1
		minDistance := math.Inf(1)

		for j := 0; j < n; j++ {
			if !visited[j] && distances[currentCity][j] < minDistance {
				minDistance = distances[currentCity][j]
				nextCity = j
			}
		}

		totalLength += minDistance
		visited[nextCity] = true
		currentCity = nextCity
	}

	// Return to the starting city
	totalLength += distances[currentCity][0]
	return totalLength
}

// Main loop to run MMAS
func (aco *ACO) Run() {
	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {
		paths := make([][]int, aco.ants)
		lengths := make([]float64, aco.ants)

		for i := 0; i < aco.ants; i++ {
			paths[i], lengths[i] = aco.constructPath(i)
		}

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

		if next == -1 {
			// Select a random unvisited city if none was selected
			unvisited := []int{}

			for j := 0; j < dimension; j++ {
				if !visited[j] {
					unvisited = append(unvisited, j)
				}
			}

			if len(unvisited) == 0 {
				break
			}

			next = unvisited[rand.Intn(len(unvisited))]
		}

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

	// Apply CMSA logic to bias towards better paths, especially in the beginning.
	q := rand.Float64()
	adaptiveCmsaProbability := aco.cmsaP * (1.0 - float64(aco.currentIteration)/float64(aco.iterations))
	if q < adaptiveCmsaProbability {
		for i := 0; i < dimension; i++ {
			if !visited[i] && aco.cmsa[current][i] > 0 {
				probabilities[i] = aco.cmsa[current][i]
				total += probabilities[i]
			}
		}

		// If we use CMSA probabilities, apply roulette-wheel selection.
		if total > 0 {
			r := rand.Float64()
			cumulativeProbability := 0.0
			for i := 0; i < dimension; i++ {
				if !visited[i] && probabilities[i] > 0 {
					cumulativeProbability += probabilities[i] / total
					if r < cumulativeProbability {
						return i
					}
				}
			}
		}
	}

	// Purely probabilistic selection of the next city.
	for i := 0; i < dimension; i++ {
		if !visited[i] {
			pheromone := utilities.FastPow(aco.pheromone[current][i], aco.alpha)
			heuristic := utilities.FastPow(1.0/aco.distances[current][i], aco.beta)
			probabilities[i] = pheromone * heuristic
			total += probabilities[i]
		}
	}

	if total == 0 {
		return -1
	}

	// Roulette-wheel selection.
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

	return -1 // In case no city is selected
}

// Global pheromone update (best ant)
func (aco *ACO) globalPheromoneUpdate(iterationBestPath []int, iterationBestLength float64) {
	// Evaporate pheromones globally
	for i := range aco.pheromone {
		for j := range aco.pheromone[i] {
			aco.pheromone[i][j] *= (1 - aco.rho)
			// Apply pheromone limits
			if aco.pheromone[i][j] < aco.tauMin {
				aco.pheromone[i][j] = aco.tauMin
			}
			if aco.pheromone[i][j] > aco.tauMax {
				aco.pheromone[i][j] = aco.tauMax
			}
		}
	}

	// Decide whether to use best-so-far or iteration-best ant for pheromone update
	// Here, we'll use the best-so-far ant
	pheromoneDeposit := 1.0 / aco.BestLength // Inverse of best tour length
	bestPath := aco.BestPath

	// Optionally, you can switch to iteration-best ant:
	// pheromoneDeposit := 1.0 / iterationBestLength
	// bestPath := iterationBestPath

	// Global update: Only the best path deposits pheromone
	for j := 0; j < len(bestPath)-1; j++ {
		start, end := bestPath[j], bestPath[j+1]
		aco.pheromone[start][end] += aco.rho * pheromoneDeposit

		// Apply pheromone limits
		if aco.pheromone[start][end] > aco.tauMax {
			aco.pheromone[start][end] = aco.tauMax
		}
		if aco.pheromone[start][end] < aco.tauMin {
			aco.pheromone[start][end] = aco.tauMin
		}
	}

	// Handle the wrap-around from the last to the first node
	last, first := bestPath[len(bestPath)-1], bestPath[0]
	aco.pheromone[last][first] += aco.rho * pheromoneDeposit

	// Apply pheromone limits
	if aco.pheromone[last][first] > aco.tauMax {
		aco.pheromone[last][first] = aco.tauMax
	}
	if aco.pheromone[last][first] < aco.tauMin {
		aco.pheromone[last][first] = aco.tauMin
	}
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
