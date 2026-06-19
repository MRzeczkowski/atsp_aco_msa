package experiments

type ExperimentsData struct {
	ExperimentParameters
	Results []ExperimentResult
}

type ExperimentParameters struct {
	Alpha, Beta, Rho, HeuristicWeight, MsaPatchBias float64
	RandomSeed                                      int64
	Iterations                                      int
}

type ExperimentResult struct {
	BestAtIteration, ThreeOptImprovementsCount int
	BestTour                                   []int
	DeviationPerIteration                      []float64
}

type ExperimentsDataStatistics struct {
	ExperimentParameters
	MinBestAtIteration                                                               int
	AverageBestAtIteration                                                           float64
	MaxBestAtIteration                                                               int
	MinThreeOptImprovementsCount                                                     int
	AverageThreeOptImprovementsCount                                                 float64
	MaxThreeOptImprovementsCount                                                     int
	MinBestDeviation, AverageBestDeviation, MaxBestDeviation, SuccessRate            float64
	MinDeviationPerIteration, AverageDeviationPerIteration, MaxDeviationPerIteration []float64
}

type HeuristicExperimentStatistics struct {
	Heuristic  string
	Statistics ExperimentsDataStatistics
}
