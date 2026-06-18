package app

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
	instanceSetSmoke      = "smoke"
	instanceSetTuning     = "tuning"
	instanceSetEvaluation = "evaluation"
	// Temporary prototype set for fast MSA impact experiments. Remove with the msa-impact pipeline.
	instanceSetMsaImpact = "msa-impact"
	instanceSetAllKnown  = "all-known"
)

const (
	runModeExperiment = "experiment"
	runModeAnalyze    = "analyze"
	runModeAll        = "all"
	runModeFinal      = "final"
	runModeFinal3Opt  = "final+3opt"
	runModeRebuildMsa = "rebuild-msa"
	// Temporary prototype mode for fast MSA impact experiments. Remove after MSA integration is settled.
	runModeMsaImpact = "msa-impact"
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
	// Temporary prototype run count for msa-impact; final/evaluation runs keep their own budgets.
	msaImpactNumberOfExperiments = 10
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
var msaImpactHeuristicWeights = []float64{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}

var smokeInstanceFiles = []string{
	"br17.atsp",
}

// Temporary prototype instance set for quick MSA impact checks. It is not intended
// as a final benchmark set and should be removed with the msa-impact mode.
var msaImpactInstanceFiles = []string{
	"atex5.atsp",
	"ftv64.atsp",
	"ftv90.atsp",
	"crane100_1.atsp",
	"td100_1.atsp",
	"ry48p.atsp",
	"dc134.atsp",
}

var tuningInstanceFiles = []string{
	"ftv33.atsp",
	"p43.atsp",
	"ft53.atsp",
	"ftv64.atsp",
	"crane66_1.atsp",
	"atex5.atsp",
	"ftv90.atsp",
	"ry48p.atsp",
	"crane100_1.atsp",
	"td100_1.atsp",
	"ftv120.atsp",
	"dc134.atsp",
	"ftv150.atsp",
	"dc188.atsp",
	"rbg323.atsp",
}

var evaluationInstanceFiles = []string{
	"atex1.atsp",
	"atex3.atsp",
	"atex4.atsp",
	"ftv35.atsp",
	"ftv38.atsp",
	"ftv44.atsp",
	"ftv47.atsp",
	"ftv55.atsp",
	"crane66_0.atsp",
	"crane66_2.atsp",
	"ft70.atsp",
	"ftv70.atsp",
	"crane100_0.atsp",
	"crane100_2.atsp",
	"ftv100.atsp",
	"ftv110.atsp",
	"dc112.atsp",
	"dc126.atsp",
	"ftv130.atsp",
	"ftv140.atsp",
	"ftv160.atsp",
	"ftv170.atsp",
	"dc176.atsp",
	"code198.atsp",
	"rbg358.atsp",
	"rbg403.atsp",
	"rbg443.atsp",
}
