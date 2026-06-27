package app

import (
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"math"
)

func setDimensionDependentParameters(dimension int, parameters *ExperimentParameters) {
	iterations := 100
	if dimension >= 50 {
		iterations = int(math.Round(30.0 * float64(dimension)))
	}
	parameters.Iterations = iterations
}

func generateParameters(heuristic string) []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0)
	heuristicWeights := []float64{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
	if heuristic == heuristicBaseline {
		heuristicWeights = []float64{defaultBaselineHeuristicWeight}
	}
	msaPatchBiases := []float64{0.0}
	if heuristic == heuristicCycleCoverMsaPatching {
		msaPatchBiases = []float64{0.0, 0.25, 0.5, 0.75, 1.0}
	}

	for _, alpha := range utilities.GenerateRange(defaultExperimentAlpha, defaultExperimentAlpha, 0.25) {
		for _, beta := range utilities.GenerateRange(defaultExperimentBeta, defaultExperimentBeta, 1.0) {
			for _, rho := range utilities.GenerateRange(defaultExperimentRho, defaultExperimentRho, 0.1) {
				for _, heuristicWeight := range heuristicWeights {
					for _, msaPatchBias := range msaPatchBiases {
						parameters = append(parameters,
							ExperimentParameters{
								Alpha:           alpha,
								Beta:            beta,
								Rho:             rho,
								HeuristicWeight: heuristicWeight,
								MsaPatchBias:    msaPatchBias,
							})
					}
				}
			}
		}
	}

	return parameters
}

type evaluationExperimentConfiguration struct {
	Heuristic            string
	Parameters           []ExperimentParameters
	SaveAllParameterRows bool
}

func evaluationExperimentConfigurations() []evaluationExperimentConfiguration {
	return []evaluationExperimentConfiguration{
		{
			Heuristic: heuristicBaseline,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(defaultBaselineHeuristicWeight),
			},
		},
		{
			Heuristic: heuristicStrictMsa,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(evaluationStrictMsaHeuristicWeight),
			},
		},
		{
			Heuristic: heuristicRootedMsa,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(evaluationRootedMsaHeuristicWeight),
			},
		},
		{
			Heuristic: heuristicCycleCover,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(evaluationCycleCoverWeight),
			},
		},
		{
			Heuristic: heuristicCycleCoverPatching,
			Parameters: []ExperimentParameters{
				newPatchingExperimentParameters(evaluationCycleCoverPatchingWeight, 0.0),
			},
		},
		{
			Heuristic: heuristicCycleCoverMsaPatching,
			Parameters: []ExperimentParameters{
				newPatchingExperimentParameters(evaluationCycleCoverMsaPatchingWeight, evaluationCycleCoverMsaPatchingMsaPatchBias),
			},
		},
	}
}

func evaluationControlExperimentConfigurations() []evaluationExperimentConfiguration {
	return []evaluationExperimentConfiguration{
		{
			Heuristic:  heuristicRandomSparse,
			Parameters: newRandomSparseEvaluationExperimentParameters(),
		},
		{
			Heuristic: heuristicDistanceRankedSparse,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(evaluationStrictMsaHeuristicWeight),
			},
		},
		{
			Heuristic:  heuristicShuffledMsa,
			Parameters: newShuffledMsaEvaluationExperimentParameters(),
		},
	}
}

func selectEvaluationExperimentConfigurations(evaluationHeuristic string) ([]evaluationExperimentConfiguration, error) {
	configurations := evaluationExperimentConfigurations()
	if evaluationHeuristic == evaluationHeuristicAll {
		return configurations, nil
	}
	if evaluationHeuristic == evaluationHeuristicControls {
		return evaluationControlExperimentConfigurations(), nil
	}

	allConfigurations := append(append([]evaluationExperimentConfiguration{}, configurations...), evaluationControlExperimentConfigurations()...)
	for _, configuration := range allConfigurations {
		if configuration.Heuristic == evaluationHeuristic {
			return []evaluationExperimentConfiguration{configuration}, nil
		}
	}

	return nil, fmt.Errorf("unsupported evaluation heuristic %q", evaluationHeuristic)
}

func newDefaultExperimentParameters(heuristicWeight float64) ExperimentParameters {
	return ExperimentParameters{
		Alpha:           defaultExperimentAlpha,
		Beta:            defaultExperimentBeta,
		Rho:             defaultExperimentRho,
		HeuristicWeight: heuristicWeight,
	}
}

func newRandomSparseEvaluationExperimentParameters() []ExperimentParameters {
	return newSeededControlEvaluationExperimentParameters(randomSparseSeeds)
}

func newShuffledMsaEvaluationExperimentParameters() []ExperimentParameters {
	return newSeededControlEvaluationExperimentParameters(shuffledMsaSeeds)
}

func newSeededControlEvaluationExperimentParameters(seeds []int64) []ExperimentParameters {
	return newSeededControlExperimentParametersForWeights(seeds, []float64{evaluationStrictMsaHeuristicWeight})
}

func newSeededControlExperimentParametersForWeights(seeds []int64, heuristicWeights []float64) []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0, len(seeds)*len(heuristicWeights))
	for _, heuristicWeight := range heuristicWeights {
		for _, randomSeed := range seeds {
			parameter := newDefaultExperimentParameters(heuristicWeight)
			parameter.RandomSeed = randomSeed
			parameters = append(parameters, parameter)
		}
	}

	return parameters
}

func newPatchingExperimentParameters(heuristicWeight, msaPatchBias float64) ExperimentParameters {
	parameters := newDefaultExperimentParameters(heuristicWeight)
	parameters.MsaPatchBias = msaPatchBias
	return parameters
}
