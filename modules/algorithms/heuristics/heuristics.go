package heuristics

const msaSupportHighSignalThreshold = 1.0

func BuildMsaSupportModifiers(msaSupport [][]float64, strength float64) [][]float64 {
	dimension := len(msaSupport)
	modifiers := BuildNeutralModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxMsaSupportSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			msaSupportSignal := msaSupport[i][j] / maxMsaSupportSelections
			if msaSupportSignal >= msaSupportHighSignalThreshold {
				modifiers[i][j] = 1.0 + msaSupportSignal*strength
			}
		}
	}

	return modifiers
}

func BuildMsaSupportCycleCoverMembershipModifiers(msaSupport, cycleCover [][]float64, strength float64, requireCycleCoverEdge bool) [][]float64 {
	dimension := len(msaSupport)
	modifiers := BuildNeutralModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxMsaSupportSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			msaSupportSignal := msaSupport[i][j] / maxMsaSupportSelections
			if msaSupportSignal < msaSupportHighSignalThreshold {
				continue
			}

			if matrixContainsEdge(cycleCover, i, j) != requireCycleCoverEdge {
				continue
			}

			modifiers[i][j] = 1.0 + msaSupportSignal*strength
		}
	}

	return modifiers
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

func matrixContainsEdge(matrix [][]float64, from, to int) bool {
	return from >= 0 && from < len(matrix) && to >= 0 && to < len(matrix[from]) && matrix[from][to] != 0
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
