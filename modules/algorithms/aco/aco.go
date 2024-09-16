package aco

import (
	"atsp_aco_msa/modules/utilities"
	"math"
	"math/rand"
)

type ACO struct {
	alpha, beta, alpha1, rho, q0, cmsaP, tau0           float64
	ants, iterations, currentIteration, BestAtIteration int
	distances, pheromone                                [][]float64
	cmsa                                                [][]float64
	BestLength                                          float64
	BestPath                                            []int
}

func NewACO(alpha, beta, alpha1, rho, q0, cmsaP float64, ants, iterations int, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromone := make([][]float64, dimension)

	L_nn := nearestNeighborTourLength(distances)
	tau0 := 1.0 / (float64(dimension) * L_nn)

	for i := range pheromone {
		pheromone[i] = make([]float64, dimension)
		for j := range pheromone[i] {
			pheromone[i][j] = tau0
		}
	}

	return &ACO{
		alpha:      alpha,
		beta:       beta,
		alpha1:     alpha1, // Global pheromone decay parameter
		rho:        rho,    // Local pheromone decay parameter
		q0:         q0,
		cmsaP:      cmsaP,
		tau0:       tau0,
		ants:       ants,
		iterations: iterations,
		distances:  distances,
		cmsa:       cmsa,
		pheromone:  pheromone,
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

// Main loop to run ACO
func (aco *ACO) Run() {
	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {
		paths := make([][]int, aco.ants)
		lengths := make([]float64, aco.ants)

		for i := 0; i < aco.ants; i++ {
			paths[i], lengths[i] = aco.constructPath(i)
		}

		// Find the best path in this iteration
		for i := 0; i < aco.ants; i++ {
			path := paths[i]
			length := lengths[i]

			if length < aco.BestLength {
				aco.BestLength = length
				aco.BestPath = append([]int(nil), path...)
				aco.BestAtIteration = aco.currentIteration
			}
		}

		aco.globalPheromoneUpdate()
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
			break
		}

		path[i] = next
		visited[next] = true

		// Apply local pheromone update after selecting the next city
		aco.localPheromoneUpdate(current, next)

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

	// Exploitation-exploration decision based on q.
	q := rand.Float64()

	// Apply CMSA logic to bias towards better paths, especially in the beginning.
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

	if q <= aco.q0 {
		// Exploitation: Greedy selection of the next city.
		bestCity := -1
		bestValue := -math.MaxFloat64

		for i := 0; i < dimension; i++ {
			if !visited[i] {
				pheromone := utilities.FastPow(aco.pheromone[current][i], aco.alpha)
				heuristic := utilities.FastPow(1.0/aco.distances[current][i], aco.beta)
				value := pheromone * heuristic

				if value > bestValue {
					bestValue = value
					bestCity = i
				}
			}
		}

		return bestCity // Deterministic choice of best city.
	} else {
		// Exploration: Probabilistic selection of the next city.
		for i := 0; i < dimension; i++ {
			if !visited[i] {
				pheromone := utilities.FastPow(aco.pheromone[current][i], aco.alpha)
				heuristic := utilities.FastPow(1.0/aco.distances[current][i], aco.beta)
				probabilities[i] = pheromone * heuristic
				total += probabilities[i]
			}
		}

		// Roulette-wheel selection.
		r := rand.Float64()
		for i, cumulativeProbability := 0, 0.0; i < dimension; i++ {
			if !visited[i] && probabilities[i] > 0.0 {
				probabilities[i] /= total
				cumulativeProbability += probabilities[i]
				if r < cumulativeProbability || math.IsNaN(probabilities[i]) {
					return i
				}
			}
		}
	}

	return -1 // In case no city is selected (shouldn't happen)
}

// Local pheromone update
func (aco *ACO) localPheromoneUpdate(from, to int) {
	aco.pheromone[from][to] = (1-aco.rho)*aco.pheromone[from][to] + aco.rho*aco.tau0
}

// Global pheromone update (best ant)
func (aco *ACO) globalPheromoneUpdate() {
	// Evaporate pheromones globally
	for i := range aco.pheromone {
		for j := range aco.pheromone[i] {
			aco.pheromone[i][j] *= (1 - aco.alpha1)
		}
	}

	// Global update: Only the best path deposits pheromone
	pheromoneDeposit := 1.0 / aco.BestLength // Inverse of best tour length
	for j := 0; j < len(aco.BestPath)-1; j++ {
		start, end := aco.BestPath[j], aco.BestPath[j+1]
		aco.pheromone[start][end] += aco.alpha1 * pheromoneDeposit
	}

	// Handle the wrap-around from the last to the first node
	last, first := aco.BestPath[len(aco.BestPath)-1], aco.BestPath[0]
	aco.pheromone[last][first] += aco.alpha1 * pheromoneDeposit
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
