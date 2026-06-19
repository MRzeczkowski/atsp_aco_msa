package app

import "atsp_aco_msa/modules/algorithms/heuristics"

func heuristicDisplayName(heuristic string) string {
	switch heuristic {
	case heuristicBaseline:
		return "Baseline"
	case heuristicStrictMsa:
		return "Strict MSA"
	case heuristicRootedMsa:
		return "Rooted MSA"
	case heuristicRandomSparse:
		return "Random sparse"
	case heuristicDistanceRankedSparse:
		return "Distance-ranked sparse"
	case heuristicShuffledMsa:
		return "Shuffled MSA"
	case heuristicStrictDistanceRanked:
		return "Strict distance-ranked sparse"
	case heuristicStrictShuffledMsa:
		return "Strict shuffled MSA"
	case heuristicRootedDistanceRanked:
		return "Rooted distance-ranked sparse"
	case heuristicRootedShuffledMsa:
		return "Rooted shuffled MSA"
	case heuristicCycleCover:
		return "Cycle cover"
	case heuristicCycleCoverMsaPatching:
		return "Cycle-cover MSA patching"
	default:
		return heuristic
	}
}

func buildHeuristicModifiers(heuristic string, matrix, msaHeuristic, cycleCover [][]float64, parameters ExperimentParameters) [][]float64 {
	switch heuristic {
	case heuristicBaseline:
		return heuristics.BuildNeutralModifiers(heuristicMatrixDimension(matrix, msaHeuristic, cycleCover))
	case heuristicStrictMsa:
		return heuristics.BuildMsaHeuristicModifiers(msaHeuristic)
	case heuristicRandomSparse:
		return heuristics.BuildRandomSparseModifiers(msaHeuristic, parameters.RandomSeed)
	case heuristicDistanceRankedSparse, heuristicStrictDistanceRanked:
		return heuristics.BuildDistanceRankedSparseModifiers(matrix, msaHeuristic)
	case heuristicShuffledMsa, heuristicStrictShuffledMsa:
		return heuristics.BuildShuffledMsaModifiers(msaHeuristic, parameters.RandomSeed)
	case heuristicCycleCover:
		return heuristics.BuildCycleCoverModifiers(cycleCover)
	case heuristicCycleCoverMsaPatching:
		return heuristics.BuildCycleCoverMsaPatchingModifiers(matrix, msaHeuristic, cycleCover, parameters.MsaPatchBias)
	default:
		return heuristics.BuildNeutralModifiers(heuristicMatrixDimension(matrix, msaHeuristic, cycleCover))
	}
}

func buildRootedHeuristicModifiers(heuristic string, matrix [][]float64, rootedMsaHeuristic [][][]float64, parameters ExperimentParameters) [][][]float64 {
	switch heuristic {
	case heuristicRootedMsa:
		return heuristics.BuildRootedMsaHeuristicModifiers(rootedMsaHeuristic)
	case heuristicRandomSparse:
		return heuristics.BuildRootedRandomSparseModifiers(rootedMsaHeuristic, parameters.RandomSeed)
	case heuristicDistanceRankedSparse, heuristicRootedDistanceRanked:
		return heuristics.BuildRootedDistanceRankedSparseModifiers(matrix, rootedMsaHeuristic)
	case heuristicShuffledMsa, heuristicRootedShuffledMsa:
		return heuristics.BuildRootedShuffledMsaModifiers(rootedMsaHeuristic, parameters.RandomSeed)
	default:
		return nil
	}
}

func heuristicMatrixDimension(matrices ...[][]float64) int {
	for _, matrix := range matrices {
		if len(matrix) != 0 {
			return len(matrix)
		}
	}

	return 0
}
