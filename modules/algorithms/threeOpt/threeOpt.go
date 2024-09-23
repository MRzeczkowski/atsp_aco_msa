package threeOpt

import (
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"slices"
	"strconv"
	"strings"
)

// Reduced 3-opt optimization that doesn't use segment reversal.
// Only one 3-opt move is valid because it can be done without reversal. It's usually referred to as case 7 or a'b'c' and is equivalent to three subsequent 2-opt moves.
// Move is performed by changing segment order from abc to acb.
// Input `tour` is changed in-place.
func ReducedThreeOpt(tour []int, distances [][]float64) {
	n := len(tour)

	dontLookBits := make([]bool, n)

	s := 3 // Minimal spacing between indices

	improves := true

	for improves {
		improves = false

	loops:
		for i := 0; i < n; i++ {

			aIdx := i
			bIdx := (i + 1) % n

			a := tour[aIdx]
			b := tour[bIdx]

			if dontLookBits[a] {
				continue
			}

			for jOffset := s; jOffset < n-s; jOffset++ {
				j := (i + jOffset) % n

				cIdx := j
				dIdx := (j + 1) % n

				c := tour[cIdx]
				d := tour[dIdx]

				for kOffset := jOffset + s; kOffset < n-s+1; kOffset++ {
					k := (i + kOffset) % n

					eIdx := k
					fIdx := (k + 1) % n

					e := tour[eIdx]
					f := tour[fIdx]

					costRemoved := distances[a][b] + distances[c][d] + distances[e][f]

					costAdded := distances[a][d] + distances[e][b] + distances[c][f]

					gain := costAdded - costRemoved

					if gain < 0 {

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
			}

			dontLookBits[a] = true
		}
	}
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
