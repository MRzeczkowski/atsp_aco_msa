package heuristics

import "math"

const msaHeuristicHighSignalThreshold = 1.0

type patchCandidate struct {
	from, to int
	signal   float64
	distance float64
	valid    bool
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
	patchingMatrix := make([][]float64, dimension)
	for i := range patchingMatrix {
		patchingMatrix[i] = make([]float64, dimension)
	}
	if dimension == 0 {
		return patchingMatrix
	}

	for i := 0; i < dimension && i < len(cycleCover); i++ {
		for j := 0; j < dimension && j < len(cycleCover[i]); j++ {
			if i != j && cycleCover[i][j] != 0 {
				patchingMatrix[i][j] = 1.0
			}
		}
	}

	components, componentCount := buildCycleCoverComponents(cycleCover, dimension)
	if componentCount <= 1 || len(msaHeuristic) == 0 {
		return patchingMatrix
	}

	addMsaConnectorHints(patchingMatrix, matrix, msaHeuristic, components, componentCount)

	return patchingMatrix
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

func buildCycleCoverComponents(cycleCover [][]float64, dimension int) ([]int, int) {
	components := make([]int, dimension)
	for i := range components {
		components[i] = -1
	}

	component := 0
	for start := 0; start < dimension; start++ {
		if components[start] != -1 {
			continue
		}

		current := start
		for current >= 0 && current < dimension && components[current] == -1 {
			components[current] = component
			current = cycleCoverSuccessor(cycleCover, current)
		}
		component++
	}

	return components, component
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

func addMsaConnectorHints(patchingMatrix, matrix, msaHeuristic [][]float64, components []int, componentCount int) {
	outgoing := make([]patchCandidate, componentCount)
	incoming := make([]patchCandidate, componentCount)
	dimension := len(components)
	maxMsaHeuristicSelections := float64(dimension - 1)
	if maxMsaHeuristicSelections <= 0 {
		return
	}

	for from := 0; from < dimension && from < len(msaHeuristic); from++ {
		for to := 0; to < dimension && to < len(msaHeuristic[from]); to++ {
			if from == to || components[from] == components[to] {
				continue
			}

			signal := msaHeuristic[from][to] / maxMsaHeuristicSelections
			if signal < msaHeuristicHighSignalThreshold {
				continue
			}

			candidate := patchCandidate{
				from:     from,
				to:       to,
				signal:   signal,
				distance: matrixDistance(matrix, from, to),
				valid:    true,
			}
			if betterPatchCandidate(candidate, outgoing[components[from]]) {
				outgoing[components[from]] = candidate
			}
			if betterPatchCandidate(candidate, incoming[components[to]]) {
				incoming[components[to]] = candidate
			}
		}
	}

	for _, candidate := range outgoing {
		if candidate.valid {
			patchingMatrix[candidate.from][candidate.to] = candidate.signal
		}
	}
	for _, candidate := range incoming {
		if candidate.valid {
			patchingMatrix[candidate.from][candidate.to] = candidate.signal
		}
	}
}

func matrixDistance(matrix [][]float64, from, to int) float64 {
	if from >= 0 && from < len(matrix) && to >= 0 && to < len(matrix[from]) {
		return matrix[from][to]
	}

	return math.Inf(1)
}

func betterPatchCandidate(candidate, current patchCandidate) bool {
	if !current.valid {
		return true
	}
	if candidate.signal != current.signal {
		return candidate.signal > current.signal
	}
	if candidate.distance != current.distance {
		return candidate.distance < current.distance
	}
	if candidate.from != current.from {
		return candidate.from < current.from
	}
	return candidate.to < current.to
}
