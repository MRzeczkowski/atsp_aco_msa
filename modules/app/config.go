package app

import (
	"atsp_aco_msa/modules/experiments"
	"atsp_aco_msa/modules/project"
)

type AtspData = project.AtspData
type ExperimentsData = experiments.ExperimentsData
type ExperimentParameters = experiments.ExperimentParameters
type ExperimentsDataStatistics = experiments.ExperimentsDataStatistics
type HeuristicExperimentStatistics = experiments.HeuristicExperimentStatistics

const (
	instanceSetSmoke      = project.InstanceSetSmoke
	instanceSetTuning     = project.InstanceSetTuning
	instanceSetEvaluation = project.InstanceSetEvaluation
	instanceSetAllKnown   = project.InstanceSetAllKnown
)

const (
	runModeExperiment     = "experiment"
	runModeAnalyze        = "analyze"
	runModeAll            = "all"
	runModeEvaluation     = "evaluation"
	runModeEvaluation3Opt = "evaluation+3opt"
	runModeRebuildCache   = "rebuild-cache"
)

const (
	analysisScopeAll          = "all"
	analysisScopeGksDeviation = "gks-deviation"
	analysisScopeTuning       = "tuning"
)

const (
	evaluationNumberOfExperiments               = 50
	evaluationStrictMsaHeuristicWeight          = 0.4
	evaluationRootedMsaHeuristicWeight          = 0.3
	evaluationCycleCoverWeight                  = 0.9
	evaluationCycleCoverMsaPatchingWeight       = 0.9
	evaluationCycleCoverMsaPatchingMsaPatchBias = 0.5
	defaultExperimentAlpha                      = 1.0
	defaultExperimentBeta                       = 2.0
	defaultExperimentRho                        = 0.8
	defaultExperimentRunCount                   = 30
	defaultBaselineHeuristicWeight              = 0.0
)

const (
	heuristicBaseline              = "baseline"
	heuristicStrictMsa             = "strict-msa"
	heuristicRootedMsa             = "rooted-msa"
	heuristicRandomSparse          = "random-sparse"
	heuristicDistanceRankedSparse  = "distance-ranked-sparse"
	heuristicShuffledMsa           = "shuffled-msa"
	heuristicStrictDistanceRanked  = "strict-distance-ranked-sparse"
	heuristicStrictShuffledMsa     = "strict-shuffled-msa"
	heuristicRootedDistanceRanked  = "rooted-distance-ranked-sparse"
	heuristicRootedShuffledMsa     = "rooted-shuffled-msa"
	heuristicCycleCover            = "cycle-cover"
	heuristicCycleCoverMsaPatching = "cycle-cover-msa-patching"
)

const evaluationHeuristicAll = "all"
const evaluationHeuristicControls = "controls"
const experimentHeuristicAll = "all"

var gksDeviationMsaPatchBiases = []float64{0.0, 0.25, 0.5, 0.75, 1.0}
var randomSparseSeeds = []int64{1, 2, 3}
var shuffledMsaSeeds = []int64{101, 102, 103}

var experimentHeuristics = []string{
	heuristicStrictMsa,
	heuristicRootedMsa,
	heuristicCycleCover,
	heuristicCycleCoverMsaPatching,
}

var evaluationResultsSummaryHeuristics = []string{
	heuristicBaseline,
	heuristicStrictMsa,
	heuristicRootedMsa,
	heuristicCycleCover,
	heuristicCycleCoverMsaPatching,
}

var msaCountScalingCounts = []int{1, 2, 4, 8, 16, 32, 64, 0}
