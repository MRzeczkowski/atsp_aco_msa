package app

import (
	"atsp_aco_msa/modules/analysis/reports"
	"atsp_aco_msa/modules/experiments"
	"atsp_aco_msa/modules/project"
)

type AtspData = project.AtspData
type ExperimentsData = experiments.ExperimentsData
type ExperimentParameters = experiments.ExperimentParameters
type ExperimentResult = experiments.ExperimentResult
type ExperimentsDataStatistics = experiments.ExperimentsDataStatistics
type HeuristicExperimentStatistics = experiments.HeuristicExperimentStatistics
type finalResultsSummaryMetric = reports.FinalResultSummaryMetric

type finalResultsSummaryRow struct {
	instance string
	metrics  map[string]finalResultsSummaryMetric
}

const (
	instanceSetSmoke      = project.InstanceSetSmoke
	instanceSetTuning     = project.InstanceSetTuning
	instanceSetEvaluation = project.InstanceSetEvaluation
	instanceSetAllKnown   = project.InstanceSetAllKnown
)

const (
	runModeExperiment   = "experiment"
	runModeAnalyze      = "analyze"
	runModeAll          = "all"
	runModeFinal        = "final"
	runModeFinal3Opt    = "final+3opt"
	runModeRebuildCache = "rebuild-cache"
)

const (
	analysisScopeAll          = "all"
	analysisScopeGksDeviation = "gks-deviation"
	analysisScopeTuning       = "tuning"
)

const (
	finalNumberOfExperiments               = 50
	finalStrictMsaHeuristicWeight          = 0.4
	finalRootedMsaHeuristicWeight          = 0.3
	finalCycleCoverWeight                  = 0.9
	finalCycleCoverMsaPatchingWeight       = 0.9
	finalCycleCoverMsaPatchingMsaPatchBias = 0.5
	defaultExperimentAlpha                 = 1.0
	defaultExperimentBeta                  = 2.0
	defaultExperimentRho                   = 0.8
	defaultExperimentRunCount              = 30
	defaultBaselineHeuristicWeight         = 0.0
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

const finalHeuristicAll = "all"
const finalHeuristicControls = "controls"
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

var finalResultsSummaryHeuristics = []string{
	heuristicBaseline,
	heuristicStrictMsa,
	heuristicRootedMsa,
	heuristicCycleCover,
	heuristicCycleCoverMsaPatching,
}

var msaCountScalingCounts = []int{1, 2, 4, 8, 16, 32, 64, 0}
