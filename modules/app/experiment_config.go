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

type finalExperimentConfiguration struct {
	Heuristic            string
	Parameters           []ExperimentParameters
	SaveAllParameterRows bool
}

func finalExperimentConfigurations() []finalExperimentConfiguration {
	return []finalExperimentConfiguration{
		{
			Heuristic: heuristicBaseline,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(defaultBaselineHeuristicWeight),
			},
		},
		{
			Heuristic: heuristicStrictMsa,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalStrictMsaHeuristicWeight),
			},
		},
		{
			Heuristic: heuristicRootedMsa,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalRootedMsaHeuristicWeight),
			},
		},
		{
			Heuristic: heuristicCycleCover,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalCycleCoverWeight),
			},
		},
		{
			Heuristic: heuristicCycleCoverMsaPatching,
			Parameters: []ExperimentParameters{
				newPatchingExperimentParameters(finalCycleCoverMsaPatchingWeight, finalCycleCoverMsaPatchingMsaPatchBias),
			},
		},
	}
}

func finalControlExperimentConfigurations() []finalExperimentConfiguration {
	return []finalExperimentConfiguration{
		{
			Heuristic:  heuristicRandomSparse,
			Parameters: newRandomSparseFinalExperimentParameters(),
		},
		{
			Heuristic: heuristicDistanceRankedSparse,
			Parameters: []ExperimentParameters{
				newDefaultExperimentParameters(finalStrictMsaHeuristicWeight),
			},
		},
		{
			Heuristic:  heuristicShuffledMsa,
			Parameters: newShuffledMsaFinalExperimentParameters(),
		},
	}
}

func selectFinalExperimentConfigurations(finalHeuristic string) ([]finalExperimentConfiguration, error) {
	configurations := finalExperimentConfigurations()
	if finalHeuristic == finalHeuristicAll {
		return configurations, nil
	}
	if finalHeuristic == finalHeuristicControls {
		return finalControlExperimentConfigurations(), nil
	}

	allConfigurations := append(append([]finalExperimentConfiguration{}, configurations...), finalControlExperimentConfigurations()...)
	for _, configuration := range allConfigurations {
		if configuration.Heuristic == finalHeuristic {
			return []finalExperimentConfiguration{configuration}, nil
		}
	}

	return nil, fmt.Errorf("unsupported final heuristic %q", finalHeuristic)
}

func newDefaultExperimentParameters(heuristicWeight float64) ExperimentParameters {
	return ExperimentParameters{
		Alpha:           defaultExperimentAlpha,
		Beta:            defaultExperimentBeta,
		Rho:             defaultExperimentRho,
		HeuristicWeight: heuristicWeight,
	}
}

func newDefaultExperimentParametersForWeights(heuristicWeights []float64) []ExperimentParameters {
	parameters := make([]ExperimentParameters, 0, len(heuristicWeights))
	for _, heuristicWeight := range heuristicWeights {
		parameters = append(parameters, newDefaultExperimentParameters(heuristicWeight))
	}
	return parameters
}

func newRandomSparseFinalExperimentParameters() []ExperimentParameters {
	return newSeededControlFinalExperimentParameters(randomSparseSeeds)
}

func newShuffledMsaFinalExperimentParameters() []ExperimentParameters {
	return newSeededControlFinalExperimentParameters(shuffledMsaSeeds)
}

func newSeededControlFinalExperimentParameters(seeds []int64) []ExperimentParameters {
	return newSeededControlExperimentParametersForWeights(seeds, []float64{finalStrictMsaHeuristicWeight})
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
