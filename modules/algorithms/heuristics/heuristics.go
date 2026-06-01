package heuristics

import "math"

const msaHeuristicHighSignalThreshold = 1.0

type cyclePatch struct {
	fromA, toB   int
	fromB, toA   int
	costIncrease float64
	msaSupport   float64
	score        float64
	valid        bool
}

func BuildMsaHeuristicModifiers(msaHeuristic [][]float64, strength float64) [][]float64 {
	dimension := len(msaHeuristic)
	modifiers := BuildNeutralModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxMsaHeuristicSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			msaHeuristicSignal := msaHeuristic[i][j] / maxMsaHeuristicSelections
			if msaHeuristicSignal >= msaHeuristicHighSignalThreshold {
				modifiers[i][j] = 1.0 + msaHeuristicSignal*strength
			}
		}
	}

	return modifiers
}

func BuildCycleCoverMsaPatchingModifiers(matrix, msaHeuristic, cycleCover [][]float64, strength float64) [][]float64 {
	patchingMatrix := BuildCycleCoverMsaPatchingMatrix(matrix, msaHeuristic, cycleCover)
	modifiers := BuildNeutralModifiers(len(patchingMatrix))
	if strength == 0 {
		return modifiers
	}

	for i := 0; i < len(patchingMatrix); i++ {
		for j := 0; j < len(patchingMatrix[i]); j++ {
			if i != j && patchingMatrix[i][j] > 0 {
				modifiers[i][j] = 1.0 + patchingMatrix[i][j]*strength
			}
		}
	}

	return modifiers
}

func BuildCycleCoverMsaPatchingMatrix(matrix, msaHeuristic, cycleCover [][]float64) [][]float64 {
	dimension := heuristicDimension(matrix, msaHeuristic, cycleCover)
	if dimension == 0 {
		return newZeroMatrix(dimension)
	}

	successors, ok := cycleCoverSuccessors(cycleCover, dimension)
	if !ok {
		return copyPositiveEdges(cycleCover, dimension)
	}

	usedInPatch := make([]bool, dimension)
	cycles := buildCycles(successors)
	for len(cycles) > 1 {
		// Karp's modified patching process:
		// 1. take a shortest current cycle,
		// 2. patch it to another cycle with a two-edge exchange,
		// 3. never reuse vertices that already served as patch endpoints.
		cycle := shortestCycle(cycles)
		patch := bestCyclePatch(matrix, msaHeuristic, successors, cycle, usedInPatch)
		if !patch.valid {
			return copyPositiveEdges(cycleCover, dimension)
		}

		applyCyclePatch(successors, patch)
		usedInPatch[patch.fromA] = true
		usedInPatch[patch.fromB] = true
		cycles = buildCycles(successors)
	}

	return buildSuccessorMatrix(successors)
}

func BuildCycleCoverModifiers(cycleCover [][]float64, strength float64) [][]float64 {
	dimension := len(cycleCover)
	modifiers := BuildNeutralModifiers(dimension)
	if dimension == 0 || strength == 0 {
		return modifiers
	}

	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i != j && cycleCover[i][j] != 0 {
				modifiers[i][j] = 1.0 + strength
			}
		}
	}

	return modifiers
}

func BuildNeutralModifiers(dimension int) [][]float64 {
	modifiers := make([][]float64, dimension)
	for i := range modifiers {
		modifiers[i] = make([]float64, dimension)
		for j := range modifiers[i] {
			modifiers[i][j] = 1.0
		}
	}
	return modifiers
}

func heuristicDimension(matrices ...[][]float64) int {
	for _, matrix := range matrices {
		if len(matrix) != 0 {
			return len(matrix)
		}
	}

	return 0
}

func newZeroMatrix(dimension int) [][]float64 {
	matrix := make([][]float64, dimension)
	for i := range matrix {
		matrix[i] = make([]float64, dimension)
	}

	return matrix
}

func copyPositiveEdges(matrix [][]float64, dimension int) [][]float64 {
	copyMatrix := newZeroMatrix(dimension)
	for i := 0; i < dimension && i < len(matrix); i++ {
		for j := 0; j < dimension && j < len(matrix[i]); j++ {
			if i != j && matrix[i][j] != 0 {
				copyMatrix[i][j] = 1.0
			}
		}
	}

	return copyMatrix
}

func cycleCoverSuccessors(cycleCover [][]float64, dimension int) ([]int, bool) {
	successors := make([]int, dimension)
	inDegree := make([]int, dimension)
	for vertex := range successors {
		successor := cycleCoverSuccessor(cycleCover, vertex)
		if successor < 0 || successor >= dimension {
			return nil, false
		}
		successors[vertex] = successor
		inDegree[successor]++
	}

	for _, degree := range inDegree {
		if degree != 1 {
			return nil, false
		}
	}

	return successors, true
}

func cycleCoverSuccessor(cycleCover [][]float64, vertex int) int {
	if vertex < 0 || vertex >= len(cycleCover) {
		return -1
	}

	for to, value := range cycleCover[vertex] {
		if vertex != to && value != 0 {
			return to
		}
	}

	return -1
}

func buildCycles(successors []int) [][]int {
	visited := make([]bool, len(successors))
	cycles := make([][]int, 0)

	for start := 0; start < len(successors); start++ {
		if visited[start] {
			continue
		}

		cycle := make([]int, 0)
		current := start
		for !visited[current] {
			visited[current] = true
			cycle = append(cycle, current)
			current = successors[current]
		}
		cycles = append(cycles, cycle)
	}

	return cycles
}

func shortestCycle(cycles [][]int) []int {
	shortest := cycles[0]
	for _, cycle := range cycles[1:] {
		if len(cycle) < len(shortest) {
			shortest = cycle
		}
	}

	return shortest
}

func bestCyclePatch(matrix, msaHeuristic [][]float64, successors []int, cycle []int, usedInPatch []bool) cyclePatch {
	best := cyclePatch{}
	inCycle := make([]bool, len(successors))
	for _, vertex := range cycle {
		inCycle[vertex] = true
	}

	for _, fromA := range cycle {
		if usedInPatch[fromA] {
			continue
		}
		for fromB := range successors {
			if inCycle[fromB] || usedInPatch[fromB] {
				continue
			}

			patch := newCyclePatch(matrix, msaHeuristic, successors, fromA, fromB)
			if betterCyclePatch(patch, best) {
				best = patch
			}
		}
	}

	return best
}

func newCyclePatch(matrix, msaHeuristic [][]float64, successors []int, fromA, fromB int) cyclePatch {
	toA := successors[fromA]
	toB := successors[fromB]
	costIncrease := matrixDistance(matrix, fromA, toB) +
		matrixDistance(matrix, fromB, toA) -
		matrixDistance(matrix, fromA, toA) -
		matrixDistance(matrix, fromB, toB)
	if math.IsInf(costIncrease, 0) || math.IsNaN(costIncrease) {
		return cyclePatch{}
	}

	msaSupport := msaHeuristicSignal(msaHeuristic, fromA, toB) +
		msaHeuristicSignal(msaHeuristic, fromB, toA)
	score := costIncrease / (1.0 + msaSupport)
	if math.IsInf(score, 0) || math.IsNaN(score) {
		return cyclePatch{}
	}

	return cyclePatch{
		fromA:        fromA,
		toB:          toB,
		fromB:        fromB,
		toA:          toA,
		costIncrease: costIncrease,
		msaSupport:   msaSupport,
		score:        score,
		valid:        true,
	}
}

func applyCyclePatch(successors []int, patch cyclePatch) {
	successors[patch.fromA] = patch.toB
	successors[patch.fromB] = patch.toA
}

func msaHeuristicSignal(msaHeuristic [][]float64, from, to int) float64 {
	maxMsaHeuristicSelections := float64(len(msaHeuristic) - 1)
	if maxMsaHeuristicSelections <= 0 ||
		from < 0 || from >= len(msaHeuristic) ||
		to < 0 || to >= len(msaHeuristic[from]) {
		return 0
	}

	return msaHeuristic[from][to] / maxMsaHeuristicSelections
}

func betterCyclePatch(candidate, current cyclePatch) bool {
	if !current.valid {
		return true
	}
	if candidate.score != current.score {
		return candidate.score < current.score
	}
	if candidate.msaSupport != current.msaSupport {
		return candidate.msaSupport > current.msaSupport
	}
	if candidate.costIncrease != current.costIncrease {
		return candidate.costIncrease < current.costIncrease
	}
	if candidate.fromA != current.fromA {
		return candidate.fromA < current.fromA
	}
	if candidate.toB != current.toB {
		return candidate.toB < current.toB
	}
	if candidate.fromB != current.fromB {
		return candidate.fromB < current.fromB
	}
	return candidate.toA < current.toA
}

func matrixDistance(matrix [][]float64, from, to int) float64 {
	if from >= 0 && from < len(matrix) && to >= 0 && to < len(matrix[from]) {
		return matrix[from][to]
	}

	return math.Inf(1)
}

func buildSuccessorMatrix(successors []int) [][]float64 {
	matrix := newZeroMatrix(len(successors))
	for from, to := range successors {
		if from != to && to >= 0 && to < len(successors) {
			matrix[from][to] = 1.0
		}
	}

	return matrix
}
