package threeopt

import "atsp_aco_msa/modules/utilities"

var spacing int = 3 // Minimal spacing between indices

type Reduced struct {
	distances          [][]float64
	neighborsLists     [][]int
	dontLookBits       []bool
	positions, newTour []int
	Improvements       int
}

func NewReduced(distances [][]float64, neighborsLists [][]int) *Reduced {
	n := len(distances)

	return &Reduced{
		distances:      distances,
		neighborsLists: neighborsLists,
		dontLookBits:   make([]bool, n),
		positions:      make([]int, n),
		newTour:        make([]int, n),
		Improvements:   0,
	}
}

// Reduced 3-opt optimization that doesn't use segment reversal.
// Only one 3-opt move is valid because it can be done without reversal. It's usually referred to as case 7 or a'b'c' and is equivalent to three subsequent 2-opt moves.
// Move is performed by changing segment order from abc to acb.
// Input `tour` is changed in-place.
func (opt *Reduced) Run(tour []int) {
	n := len(tour)

	// Initialize total costs and thresholds
	improvementThreshold := 0.001 // 0.1% of initial tour cost
	maxConsecutiveMinorGain := 20 // Allowed minor improvements

	currentLength := utilities.TourLength(tour, opt.distances)
	minorGainThreshold := -currentLength * improvementThreshold
	consecutiveMinorGain := 0

	setPositions(opt.positions, tour)

	for i := 0; i < n; i++ {
		opt.dontLookBits[i] = false
	}

	improves := true

	for improves {
		improves = false

	loops:
		for i := 0; i < n; i++ {
			aIdx := i
			a := tour[aIdx]

			bIdx := (i + 1) % n
			b := tour[bIdx]

			if opt.dontLookBits[a] {
				continue
			}

			distAB := opt.distances[a][b]

			for _, d := range opt.neighborsLists[a] {
				distAD := opt.distances[a][d]

				if distAD >= distAB {
					break
				}

				dIdx := opt.positions[d]

				j := (dIdx - 1 + n) % n
				c := tour[j]

				if !haveCorrectSpacing(i, j, n) {
					continue
				}

				distCD := opt.distances[c][d]

				radius := distAB + distCD - distAD

				for _, f := range opt.neighborsLists[c] {
					distCF := opt.distances[c][f]

					if distCF >= radius {
						break
					}

					fIdx := opt.positions[f]

					k := (fIdx - 1 + n) % n
					e := tour[k]

					if !(isBetween(i, j, k) && haveCorrectSpacing(j, k, n) && haveCorrectSpacing(k, i, n)) {
						continue
					}

					distEF := opt.distances[e][f]

					costRemoved := distAB + distCD + distEF

					distEB := opt.distances[e][b]

					costAdded := distAD + distEB + distCF

					gain := costAdded - costRemoved

					if gain >= 0 {
						continue
					}

					if gain > minorGainThreshold {
						consecutiveMinorGain++
						// Early termination check
						if consecutiveMinorGain > maxConsecutiveMinorGain {
							return // Exit early
						}
						continue
					}
					consecutiveMinorGain = 0 // Reset counter

					pos := 0

					// First Segment
					if aIdx < fIdx || bIdx < fIdx {
						pos += copy(opt.newTour[pos:], tour[fIdx:])
						pos += copy(opt.newTour[pos:], tour[:bIdx])
					} else {
						pos += copy(opt.newTour[pos:], tour[fIdx:bIdx])
					}

					// Third Segment
					if dIdx > fIdx {
						pos += copy(opt.newTour[pos:], tour[dIdx:])
						pos += copy(opt.newTour[pos:], tour[:fIdx])
					} else {
						pos += copy(opt.newTour[pos:], tour[dIdx:fIdx])
					}

					// Second Segment
					if bIdx > dIdx {
						pos += copy(opt.newTour[pos:], tour[bIdx:])
						pos += copy(opt.newTour[pos:], tour[:dIdx])
					} else {
						pos += copy(opt.newTour[pos:], tour[bIdx:dIdx])
					}

					copy(tour, opt.newTour)
					currentLength += gain
					minorGainThreshold = -currentLength * improvementThreshold

					setPositions(opt.positions, tour)

					opt.dontLookBits[a] = false
					opt.dontLookBits[b] = false
					opt.dontLookBits[c] = false
					opt.dontLookBits[d] = false
					opt.dontLookBits[e] = false
					opt.dontLookBits[f] = false

					improves = true
					opt.Improvements++
					break loops // Exit after applying a move
				}
			}

			opt.dontLookBits[a] = true
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

func haveCorrectSpacing(i, j, n int) bool {
	distIJ := calculateDistanceInRingBuffer(i, j, n)
	distJI := calculateDistanceInRingBuffer(j, i, n)

	return distIJ >= spacing && distJI >= spacing
}

func calculateDistanceInRingBuffer(a, b, n int) int {
	dist := b - a
	if dist < 0 {
		dist += n
	}

	return dist
}
