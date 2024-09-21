package aco

import (
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"math"
	"math/rand"
	"slices"
	"sort"
	"strconv"
	"strings"
	"sync"
)

type ACO struct {
	alpha, beta, rho, pDec, pCmsa                       float64
	ants, iterations, currentIteration, BestAtIteration int
	distances, pheromones                               [][]float64
	cmsa                                                [][]float64
	tauMin, tauMax, BestLength                          float64
	BestPath                                            []int
	neighborList                                        [][]int
}

func NewACO(alpha, beta, rho, pBest, pCmsa float64, ants, iterations int, distances, cmsa [][]float64) *ACO {
	dimension := len(distances)
	pheromones := make([][]float64, dimension)

	k := 5
	neighborList := buildNearestNeighborList(distances, k)

	// Arbitrary high value: `4.3. Pheromone trail initialization`
	for i := range pheromones {
		pheromones[i] = make([]float64, dimension)
		for j := range pheromones[i] {
			pheromones[i][j] = 1.0 + (pCmsa * cmsa[i][j])
		}
	}

	return &ACO{
		alpha:        alpha,
		beta:         beta,
		rho:          rho,
		pDec:         math.Pow(pBest, 1.0/float64(dimension)),
		pCmsa:        pCmsa,
		ants:         ants,
		iterations:   iterations,
		distances:    distances,
		cmsa:         cmsa,
		pheromones:   pheromones,
		BestLength:   math.Inf(1),
		neighborList: neighborList,
	}
}

// Main loop to run MMAS
func (aco *ACO) Run() {

	// path := []int{13, 1, 8, 4, 3, 15, 5, 14, 6, 16, 0, 11, 7, 2, 10, 9, 12}
	// path := []int{0, 13, 9, 32, 7, 8, 12, 14, 15, 16, 1, 3, 24, 23, 19, 17, 10, 18, 31, 21, 20, 22, 26, 27, 28, 29, 25, 2, 33, 30, 4, 6, 5, 11}

	// fmt.Println(aco.pathLength(path))
	// aco.threeOpt(path)
	// fmt.Println(aco.pathLength(path))

	for aco.currentIteration = 0; aco.currentIteration < aco.iterations; aco.currentIteration++ {
		paths := make([][]int, aco.ants)
		lengths := make([]float64, aco.ants)

		// for i := 0; i < aco.ants; i++ {
		// 	paths[i], lengths[i] = aco.constructPath(i)
		// }

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

	// Apply Reduced 3-opt to improve the path
	aco.threeOpt(path)
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

// Initialize the don't look bits (one for each city)
func initializeDontLookBits(n int) []bool {
	dontLook := make([]bool, n) // All bits start as false (look at all cities)
	return dontLook
}

// Build the nearest neighbor list for each city
func buildNearestNeighborList(distances [][]float64, k int) [][]int {
	n := len(distances)
	neighborList := make([][]int, n)

	for i := 0; i < n; i++ {
		// Create a list of city indices sorted by distance from city i
		type nodeDist struct {
			id       int
			distance float64
		}

		cityDistances := make([]nodeDist, n)
		for j := 0; j < n; j++ {
			cityDistances[j] = nodeDist{id: j, distance: distances[i][j]}
		}

		// Sort cities by distance from city i
		sort.Slice(cityDistances, func(a, b int) bool {
			return cityDistances[a].distance < cityDistances[b].distance
		})

		// Take only the k nearest neighbors (excluding the city itself)
		neighbors := make([]int, 0, k)
		for j := 1; j <= k && j < n; j++ {
			neighbors = append(neighbors, cityDistances[j].id)
		}

		neighborList[i] = neighbors
	}

	return neighborList
}

// reduced threeOpt optimization that doesn't use segment reversal
func (aco *ACO) threeOpt(path []int) {
	n := len(path)

	improves := true

	// Continue iterating until no improvements are found
	for improves {
		improves = false

	loops:
		for i := 0; i < n-6; i++ {
			for j := i + 3; j < n-4; j++ {
				for k := j + 3; k < n-2; k++ {
					a := path[i]
					b := path[(i+1)%n]
					c := path[j]
					d := path[(j+1)%n]
					e := path[k]
					f := path[(k+1)%n]

					costRemoved := aco.distances[a][b] + aco.distances[c][d] + aco.distances[e][f]

					costAdded := aco.distances[a][d] + aco.distances[e][b] + aco.distances[c][f]

					gain := costAdded - costRemoved

					if gain < 0 {
						beforeMoveLength := aco.pathLength(path)

						var firstSegment []int

						if i < ((k + 1) % n) {
							firstSegment = slices.Concat(path[(k+1)%n:], path[:i+1])
						} else {
							firstSegment = path[(k+1)%n : i+1]
						}

						secondSegment := path[i+1 : j+1]
						thirdSegment := path[j+1 : (k+1)%n]

						newPath := slices.Concat(firstSegment, thirdSegment, secondSegment)

						afterMoveLength := aco.pathLength(newPath)

						for z := 0; z < n; z++ {

							if indexOf(z, newPath) == -1 {
								fmt.Printf("Missing element %d in new path!\n", i)

								fmt.Printf("\tMade move for gain %.2f: i=%d (Node %d), j=%d (Node %d), k=%d (Node %d)\n", gain, i, a, j, c, k, e)

								fmt.Println()

								ab := aco.distances[a][b]
								cd := aco.distances[c][d]
								ef := aco.distances[e][f]

								previousCost := ab + cd + ef

								fmt.Printf("\tPrevious cost = ab + cd + ef : %.0f = %.0f + %.0f + %.0f\n", previousCost, ab, cd, ef)

								ad := aco.distances[a][d]
								eb := aco.distances[e][b]
								cf := aco.distances[c][f]

								newCost := ad + eb + cf

								fmt.Printf("\tNew cost =      ad + eb + cf : %.0f = %.0f + %.0f + %.0f\n", newCost, ad, eb, cf)

								previousWraparoundCost := aco.distances[path[n-1]][path[0]]

								fmt.Println("\tPrevious wraparound cost = ", previousWraparoundCost)

								newWraparoundCost := aco.distances[newPath[n-1]][newPath[0]]

								fmt.Println("\tNew wraparound cost =      ", newWraparoundCost)

								fmt.Println()

								test := make([]string, n)

								for i := 0; i < n; i++ {
									len := len(strconv.Itoa(path[i]))
									test[i] = strings.Repeat(" ", len)
								}

								test[i] = "a" + strings.Repeat(" ", len(strconv.Itoa(path[i]))-1)
								test[(i+1)%n] = "b" + strings.Repeat(" ", len(strconv.Itoa(path[(i+1)%n]))-1)
								test[j] = "c" + strings.Repeat(" ", len(strconv.Itoa(path[j]))-1)
								test[(j+1)%n] = "d" + strings.Repeat(" ", len(strconv.Itoa(path[(j+1)%n]))-1)
								test[k] = "e" + strings.Repeat(" ", len(strconv.Itoa(path[k]))-1)
								test[(k+1)%n] = "f" + strings.Repeat(" ", len(strconv.Itoa(path[(k+1)%n]))-1)
								fmt.Println("\t            ", test)

								fmt.Println("\tBefore Move:", path)
								fmt.Println()

								for i := 0; i < n; i++ {
									len := len(strconv.Itoa(newPath[i]))
									test[i] = strings.Repeat(" ", len)
								}

								aIdx := indexOf(a, newPath)
								bIdx := indexOf(b, newPath)
								cIdx := indexOf(c, newPath)
								dIdx := indexOf(d, newPath)
								eIdx := indexOf(e, newPath)
								fIdx := indexOf(f, newPath)

								if aIdx != -1 {
									test[aIdx] = "a" + strings.Repeat(" ", len(strconv.Itoa(newPath[aIdx]))-1)
								}
								if bIdx != -1 {
									test[bIdx] = "b" + strings.Repeat(" ", len(strconv.Itoa(newPath[bIdx]))-1)
								}
								if cIdx != -1 {
									test[cIdx] = "c" + strings.Repeat(" ", len(strconv.Itoa(newPath[cIdx]))-1)
								}
								if dIdx != -1 {
									test[dIdx] = "d" + strings.Repeat(" ", len(strconv.Itoa(newPath[dIdx]))-1)
								}
								if eIdx != -1 {
									test[eIdx] = "e" + strings.Repeat(" ", len(strconv.Itoa(newPath[eIdx]))-1)
								}
								if fIdx != -1 {
									test[fIdx] = "f" + strings.Repeat(" ", len(strconv.Itoa(newPath[fIdx]))-1)
								}

								fmt.Println("\t            ", test)

								fmt.Println("\tAfter Move: ", newPath)

								fmt.Println()

								fmt.Println("\tPath length before move:", beforeMoveLength)
								fmt.Println("\tPath length after move:", afterMoveLength)

								fmt.Println()
							}
						}

						if afterMoveLength != beforeMoveLength+gain {
							fmt.Println("Error in gain calculation!")

							fmt.Printf("\tMade move for gain %.2f: i=%d (Node %d), j=%d (Node %d), k=%d (Node %d)\n", gain, i, a, j, c, k, e)

							fmt.Println()

							ab := aco.distances[a][b]
							cd := aco.distances[c][d]
							ef := aco.distances[e][f]

							previousCost := ab + cd + ef

							fmt.Printf("\tPrevious cost = ab + cd + ef : %.0f = %.0f + %.0f + %.0f\n", previousCost, ab, cd, ef)

							ad := aco.distances[a][d]
							eb := aco.distances[e][b]
							cf := aco.distances[c][f]

							newCost := ad + eb + cf

							fmt.Printf("\tNew cost =      ad + eb + cf : %.0f = %.0f + %.0f + %.0f\n", newCost, ad, eb, cf)

							previousWraparoundCost := aco.distances[path[n-1]][path[0]]

							fmt.Println("\tPrevious wraparound cost = ", previousWraparoundCost)

							newWraparoundCost := aco.distances[newPath[n-1]][newPath[0]]

							fmt.Println("\tNew wraparound cost =      ", newWraparoundCost)

							fmt.Println()

							test := make([]string, n)

							for i := 0; i < n; i++ {
								len := len(strconv.Itoa(path[i]))
								test[i] = strings.Repeat(" ", len)
							}

							test[i] = "a" + strings.Repeat(" ", len(strconv.Itoa(path[i]))-1)
							test[(i+1)%n] = "b" + strings.Repeat(" ", len(strconv.Itoa(path[(i+1)%n]))-1)
							test[j] = "c" + strings.Repeat(" ", len(strconv.Itoa(path[j]))-1)
							test[(j+1)%n] = "d" + strings.Repeat(" ", len(strconv.Itoa(path[(j+1)%n]))-1)
							test[k] = "e" + strings.Repeat(" ", len(strconv.Itoa(path[k]))-1)
							test[(k+1)%n] = "f" + strings.Repeat(" ", len(strconv.Itoa(path[(k+1)%n]))-1)
							fmt.Println("\t            ", test)

							fmt.Println("\tBefore Move:", path)
							fmt.Println()

							for i := 0; i < n; i++ {
								len := len(strconv.Itoa(newPath[i]))
								test[i] = strings.Repeat(" ", len)
							}

							aIdx := indexOf(a, newPath)
							bIdx := indexOf(b, newPath)
							cIdx := indexOf(c, newPath)
							dIdx := indexOf(d, newPath)
							eIdx := indexOf(e, newPath)
							fIdx := indexOf(f, newPath)

							test[aIdx] = "a" + strings.Repeat(" ", len(strconv.Itoa(newPath[aIdx]))-1)
							test[bIdx] = "b" + strings.Repeat(" ", len(strconv.Itoa(newPath[bIdx]))-1)
							test[cIdx] = "c" + strings.Repeat(" ", len(strconv.Itoa(newPath[cIdx]))-1)
							test[dIdx] = "d" + strings.Repeat(" ", len(strconv.Itoa(newPath[dIdx]))-1)
							test[eIdx] = "e" + strings.Repeat(" ", len(strconv.Itoa(newPath[eIdx]))-1)
							test[fIdx] = "f" + strings.Repeat(" ", len(strconv.Itoa(newPath[fIdx]))-1)
							fmt.Println("\t            ", test)

							fmt.Println("\tAfter Move: ", newPath)

							fmt.Println()

							fmt.Println("\tPath length before move:", beforeMoveLength)
							fmt.Println("\tPath length after move:", afterMoveLength)

							fmt.Println()
						}

						copy(path, newPath)

						improves = true
						break loops // Exit loops after applying a move
					}
				}
			}
		}
	}
}

func indexOf(element int, data []int) int {
	for k, v := range data {
		if element == v {
			return k
		}
	}

	return -1 //not found.
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

	var pheromoneDeposit float64
	var bestPath []int

	if aco.currentIteration <= 25 {
		pheromoneDeposit = 1.0 / iterationBestLength
		bestPath = iterationBestPath
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
			bestPath = aco.BestPath
		} else {
			pheromoneDeposit = 1.0 / iterationBestLength
			bestPath = iterationBestPath
		}
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
