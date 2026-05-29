package heuristics

const msaHeuristicHighSignalThreshold = 1.0

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
