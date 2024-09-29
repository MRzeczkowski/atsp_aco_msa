package threeOpt

import (
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"slices"
	"sort"
	"strconv"
	"strings"
)

var spacing int = 3 // Minimal spacing between indices

type ReducedThreeOpt struct {
	distances      [][]float64
	neighborsLists [][]int
}

func NewReducedThreeOpt(distances [][]float64, k int) *ReducedThreeOpt {

	return &ReducedThreeOpt{
		distances:      distances,
		neighborsLists: buildNearestNeighborsLists(distances, k),
	}
}

// Build the nearest k neighbors list for each city
func buildNearestNeighborsLists(distances [][]float64, k int) [][]int {
	n := len(distances)
	neighborsLists := make([][]int, n)
	k = min(k, n)

	for i := 0; i < n; i++ {

		type nodeDist struct {
			id       int
			distance float64
		}

		cityDistances := make([]nodeDist, n)
		for j := 0; j < n; j++ {
			cityDistances[j] = nodeDist{id: j, distance: distances[i][j]}
		}

		sort.Slice(cityDistances, func(x, y int) bool {
			return cityDistances[x].distance < cityDistances[y].distance
		})

		for j := 0; j < k; j++ {
			neighborsLists[i] = append(neighborsLists[i], cityDistances[j].id)
		}
	}

	return neighborsLists
}

// Reduced 3-opt optimization that doesn't use segment reversal.
// Only one 3-opt move is valid because it can be done without reversal. It's usually referred to as case 7 or a'b'c' and is equivalent to three subsequent 2-opt moves.
// Move is performed by changing segment order from abc to acb.
// Input `tour` is changed in-place.
func (threeOpt *ReducedThreeOpt) Run(tour []int) {

	n := len(tour)

	dontLookBits := make([]bool, n)

	positions := make([]int, n)

	setPositions(positions, tour)

	improves := true

	for improves {
		improves = false

	loops:
		for i := 0; i < n; i++ {
			aIdx := i
			a := tour[aIdx]

			bIdx := (i + 1) % n
			b := tour[bIdx]

			if dontLookBits[a] {
				continue
			}

			for _, d := range threeOpt.neighborsLists[a] {
				dIdx := positions[d]

				j := (dIdx - 1 + n) % n
				c := tour[j]

				distIJ := calculateDistanceInRingBuffer(i, j, n)
				distJI := calculateDistanceInRingBuffer(j, i, n)

				// We don't consider neighbors that are too close.
				if distIJ < spacing || distJI < spacing {
					continue
				}

				for _, f := range threeOpt.neighborsLists[c] {
					fIdx := positions[f]

					k := (fIdx - 1 + n) % n
					e := tour[k]

					distJK := calculateDistanceInRingBuffer(j, k, n)
					distKJ := calculateDistanceInRingBuffer(k, j, n)
					distKI := calculateDistanceInRingBuffer(k, i, n)
					distIK := calculateDistanceInRingBuffer(i, k, n)

					// We don't consider neighbors that are too close.
					if distJK < spacing || distKJ < spacing || distKI < spacing || distIK < spacing {
						continue
					}

					// j must be between i and k
					if !isBetween(i, j, k) {
						continue
					}

					costRemoved := threeOpt.distances[a][b] + threeOpt.distances[c][d] + threeOpt.distances[e][f]

					costAdded := threeOpt.distances[a][d] + threeOpt.distances[e][b] + threeOpt.distances[c][f]

					gain := costAdded - costRemoved

					if gain >= 0 {
						break
					}

					var firstSegment []int
					var secondSegment []int
					var thirdSegment []int

					if aIdx < fIdx || bIdx < fIdx {
						firstSegment = slices.Concat(tour[fIdx:], tour[:bIdx])
					} else {
						firstSegment = tour[fIdx:bIdx]
					}

					if bIdx > dIdx {
						secondSegment = slices.Concat(tour[bIdx:], tour[:dIdx])
					} else {
						secondSegment = tour[bIdx:dIdx]
					}

					if dIdx > fIdx {
						thirdSegment = slices.Concat(tour[dIdx:], tour[:fIdx])
					} else {
						thirdSegment = tour[dIdx:fIdx]
					}

					newTour := slices.Concat(firstSegment, thirdSegment, secondSegment)

					copy(tour, newTour)

					setPositions(positions, tour)

					dontLookBits[a] = false
					dontLookBits[b] = false
					dontLookBits[c] = false
					dontLookBits[d] = false
					dontLookBits[e] = false
					dontLookBits[f] = false

					improves = true
					break loops // Exit after applying a move
				}
			}

			dontLookBits[a] = true
		}
	}
}

func setPositions(positions, tour []int) {
	for i, v := range tour {
		positions[v] = i
	}
}

// Returns true if `x` is between `a` and `b` in cyclic sequence
// with established orientation of nondescending indices
// That is: when one begins a forward traversal of the tour
// at position `a', then position `x` is reached before position `b'.
// Returns true if and only if:
// a < x < b  or  b < a < x  or  x < b < a
func isBetween(a, x, b int) bool {
	if a < b {
		return x > a && x < b
	}

	return x > a || x < b
}

func calculateDistanceInRingBuffer(a, b, n int) int {
	dist := b - a
	if dist < 0 {
		dist += n
	}

	return dist
}

func logReducedThreeOptIssue(tour, newTour []int, distances [][]float64, gain, beforeMoveLength, afterMoveLength float64, i, j, k, a, b, c, d, e, f, n int) {
	fmt.Printf("\tMade move for gain %.2f: i=%d (Node %d), j=%d (Node %d), k=%d (Node %d)\n", gain, i, a, j, c, k, e)

	fmt.Println()

	ab := distances[a][b]
	cd := distances[c][d]
	ef := distances[e][f]

	previousCost := ab + cd + ef

	fmt.Printf("\tPrevious cost = ab + cd + ef : %.0f = %.0f + %.0f + %.0f\n", previousCost, ab, cd, ef)

	ad := distances[a][d]
	eb := distances[e][b]
	cf := distances[c][f]

	newCost := ad + eb + cf

	fmt.Printf("\tNew cost =      ad + eb + cf : %.0f = %.0f + %.0f + %.0f\n", newCost, ad, eb, cf)

	previousWraparoundCost := distances[tour[n-1]][tour[0]]

	fmt.Println("\tPrevious wraparound cost = ", previousWraparoundCost)

	newWraparoundCost := distances[newTour[n-1]][newTour[0]]

	fmt.Println("\tNew wraparound cost =      ", newWraparoundCost)

	fmt.Println()

	test := make([]string, n)

	for i := 0; i < n; i++ {
		len := len(strconv.Itoa(tour[i]))
		test[i] = strings.Repeat(" ", len)
	}

	test[i] = "a" + strings.Repeat(" ", len(strconv.Itoa(tour[i]))-1)
	test[(i+1)%n] = "b" + strings.Repeat(" ", len(strconv.Itoa(tour[(i+1)%n]))-1)
	test[j] = "c" + strings.Repeat(" ", len(strconv.Itoa(tour[j]))-1)
	test[(j+1)%n] = "d" + strings.Repeat(" ", len(strconv.Itoa(tour[(j+1)%n]))-1)
	test[k] = "e" + strings.Repeat(" ", len(strconv.Itoa(tour[k]))-1)
	test[(k+1)%n] = "f" + strings.Repeat(" ", len(strconv.Itoa(tour[(k+1)%n]))-1)
	fmt.Println("\t            ", test)

	fmt.Println("\tBefore Move:", tour)
	fmt.Println()

	for i := 0; i < n; i++ {
		len := len(strconv.Itoa(newTour[i]))
		test[i] = strings.Repeat(" ", len)
	}

	aIdx := utilities.IndexOf(a, newTour)
	bIdx := utilities.IndexOf(b, newTour)
	cIdx := utilities.IndexOf(c, newTour)
	dIdx := utilities.IndexOf(d, newTour)
	eIdx := utilities.IndexOf(e, newTour)
	fIdx := utilities.IndexOf(f, newTour)

	test[aIdx] = "a" + strings.Repeat(" ", len(strconv.Itoa(newTour[aIdx]))-1)
	test[bIdx] = "b" + strings.Repeat(" ", len(strconv.Itoa(newTour[bIdx]))-1)
	test[cIdx] = "c" + strings.Repeat(" ", len(strconv.Itoa(newTour[cIdx]))-1)
	test[dIdx] = "d" + strings.Repeat(" ", len(strconv.Itoa(newTour[dIdx]))-1)
	test[eIdx] = "e" + strings.Repeat(" ", len(strconv.Itoa(newTour[eIdx]))-1)
	test[fIdx] = "f" + strings.Repeat(" ", len(strconv.Itoa(newTour[fIdx]))-1)
	fmt.Println("\t            ", test)

	fmt.Println("\tAfter Move: ", newTour)

	fmt.Println()

	fmt.Println("\tTour length before move:", beforeMoveLength)
	fmt.Println("\tTour length after move:", afterMoveLength)

	fmt.Println()
}
