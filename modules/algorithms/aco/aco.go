package aco

import (
	"atsp_aco_msa/modules/utilities"
	"math"
	"math/rand"
	"sync"
)

type ACO struct {
	alpha, beta, evaporation                            float64
	minPheromone, maxPheromone, exploration, q          float64
	ants, iterations, currentIteration, BestAtIteration int
	distances, pheromone                                [][]float64
	cmsa                                                [][]float64
	BestLength                                          float64
	BestPath                                            []int
}

// NewACO initializes a new ACO instance with initial pheromone levels set to an estimated best value
func NewACO(alpha, beta, evaporation, exploration, q float64, ants, iterations int, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromone := make([][]float64, dimension)
	initialPheromone := 1.0
	for i := range pheromone {
		pheromone[i] = make([]float64, dimension)
		for j := range pheromone[i] {
			pheromone[i][j] = initialPheromone
		}
	}

	return &ACO{
		alpha:        alpha,
		beta:         beta,
		evaporation:  evaporation,
		exploration:  exploration,
		q:            q,
		ants:         ants,
		iterations:   iterations,
		distances:    distances,
		cmsa:         cmsa,
		pheromone:    pheromone,
		BestLength:   math.Inf(1),
		maxPheromone: initialPheromone,
		minPheromone: initialPheromone / (exploration * float64(ants)),
	}
}

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

		aco.updatePheromoneLevels() // Recalculate pheromone limits based on the new best solution
		aco.updatePheromone(paths, lengths)
	}
}

func (aco *ACO) updatePheromoneLevels() {
	aco.maxPheromone = 1.0 / ((1 - aco.evaporation) * aco.BestLength)
	aco.minPheromone = aco.maxPheromone / (aco.exploration * float64(aco.ants))
}

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
		current = next
	}

	length := aco.pathLength(path)
	if length < aco.BestLength {
		aco.BestLength = length
		aco.BestPath = append([]int(nil), path...)
		aco.BestAtIteration = aco.currentIteration
	}

	return path, length
}

func (aco *ACO) selectNextCity(current int, visited []bool) int {
	dimension := len(aco.distances)
	probabilities := make([]float64, dimension)
	total := 0.0

	// This should make ants use better paths in the beginning.
	adaptiveMstProbability := aco.q * (1.0 - float64(aco.currentIteration)/float64(aco.iterations))
	if rand.Float64() < adaptiveMstProbability {
		for i := 0; i < dimension; i++ {
			if !visited[i] && aco.cmsa[current][i] > 0 {
				probabilities[i] = aco.cmsa[current][i]
				total += probabilities[i]
			}
		}
	}

	for i := 0; i < dimension; i++ {
		if !visited[i] && probabilities[i] == 0 {
			pheromone := utilities.FastPow(aco.pheromone[current][i], aco.alpha)
			invDistance := 1.0 / float64(aco.distances[current][i])
			desirability := utilities.FastPow(invDistance, aco.beta)
			probabilities[i] = pheromone * desirability
			total += probabilities[i]
		}
	}

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

	return -1
}

func (aco *ACO) updatePheromone(paths [][]int, lengths []float64) {
	// Find the best path of this iteration
	bestIdx := 0
	for i := 1; i < len(lengths); i++ {
		if lengths[i] < lengths[bestIdx] {
			bestIdx = i
		}
	}

	// Evaporate pheromone first
	for i := range aco.pheromone {
		for j := range aco.pheromone[i] {
			aco.pheromone[i][j] *= (1 - aco.evaporation)
			aco.pheromone[i][j] = math.Max(aco.pheromone[i][j], aco.minPheromone) // Enforce minimum pheromone level
		}
	}

	// Strengthen pheromone trail for the best ant's path
	path := paths[bestIdx]
	delta := 1.0 / lengths[bestIdx]
	for i := 0; i < len(path)-1; i++ {
		start, end := path[i], path[i+1]
		aco.pheromone[start][end] += delta
		aco.pheromone[start][end] = math.Min(aco.pheromone[start][end], aco.maxPheromone) // Enforce maximum pheromone level
	}

	if len(path) > 0 {
		last, first := path[len(path)-1], path[0]
		aco.pheromone[last][first] += delta
		aco.pheromone[last][first] = math.Min(aco.pheromone[last][first], aco.maxPheromone)
	}
}

func (aco *ACO) pathLength(path []int) float64 {
	sum := 0.0
	p := len(path)

	for i := 0; i < p-1; i++ {
		start, end := path[i], path[i+1]
		sum += float64(aco.distances[start][end])
	}

	if p > 0 {
		last, first := path[p-1], path[0]
		sum += float64(aco.distances[last][first])
	}

	return sum
}
