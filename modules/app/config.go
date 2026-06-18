package app

import "atsp_aco_msa/modules/project"

type AtspData = project.AtspData

type ExperimentsData struct {
	ExperimentParameters
	results []ExperimentResult
}

type ExperimentParameters struct {
	alpha, beta, rho, heuristicWeight, msaPatchBias float64
	randomSeed                                      int64
	iterations                                      int
}

type ExperimentResult struct {
	bestAtIteration, threeOptImprovementsCount int
	bestTour                                   []int
	deviationPerIteration                      []float64
}

type ExperimentsDataStatistics struct {
	ExperimentParameters
	minBestAtIteration                                                               int
	averageBestAtIteration                                                           float64
	maxBestAtIteration                                                               int
	minThreeOptImprovementsCount                                                     int
	averageThreeOptImprovementsCount                                                 float64
	maxThreeOptImprovementsCount                                                     int
	minBestDeviation, averageBestDeviation, maxBestDeviation, successRate            float64
	minDeviationPerIteration, averageDeviationPerIteration, maxDeviationPerIteration []float64
}

type HeuristicExperimentStatistics struct {
	heuristic  string
	statistics ExperimentsDataStatistics
}

type finalResultsSummaryMetric struct {
	averageMinDeviation  float64
	successRate          float64
	averageBestIteration float64
	heuristicWeight      float64
	iterations           int
}

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
	runModeExperiment = "experiment"
	runModeAnalyze    = "analyze"
	runModeAll        = "all"
	runModeFinal      = "final"
	runModeFinal3Opt  = "final+3opt"
	runModeRebuildMsa = "rebuild-msa"
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
