package app

import (
	"atsp_aco_msa/modules/analysis/tours"
	"atsp_aco_msa/modules/analysis/tuning"
	"atsp_aco_msa/modules/artifacts/cyclecover"
	"atsp_aco_msa/modules/experiments"
	workerpool "atsp_aco_msa/modules/experiments/workers"
	"atsp_aco_msa/modules/project"
	"fmt"
	"sort"
	"sync"
	"time"
)

func runExperimentMode(atspsData []AtspData, heuristics []string, workers int) error {
	if err := ensureMsaHeuristicArtifacts(atspsData, workers, heuristicsUseRootedMsa(heuristics)); err != nil {
		return err
	}
	if heuristicsUseCycleCover(heuristics) {
		if err := ensureCycleCoverCache(atspsData, workers); err != nil {
			return err
		}
	}

	for _, heuristic := range heuristics {
		fmt.Printf("Running tuning heuristic %s\n", heuristic)
		if err := runExperimentSet(atspsData, heuristic, generateParameters(heuristic), defaultExperimentRunCount, workers); err != nil {
			return err
		}
	}

	return saveTuningSummary(project.ResultsDirectoryName, atspsData, heuristics)
}

func saveTuningSummary(resultsRootPath string, atspsData []AtspData, heuristics []string) error {
	heuristicConfigs := make([]tuning.HeuristicConfig, 0, len(heuristics))
	for _, heuristic := range heuristics {
		if heuristicIsSparseControl(heuristic) {
			continue
		}

		resultFiles := make([]tuning.ResultFile, 0, len(atspsData))
		for _, atspData := range atspsData {
			resultFiles = append(resultFiles, tuning.ResultFile{
				Instance: atspData.Name,
				Path:     resultFilePathForHeuristic(atspData, heuristic),
			})
		}

		heuristicConfigs = append(heuristicConfigs, tuning.HeuristicConfig{
			Name:                heuristic,
			DisplayName:         heuristicDisplayName(heuristic),
			IncludeMsaPatchBias: heuristic == heuristicCycleCoverMsaPatching,
			ResultFiles:         resultFiles,
		})
	}

	return tuning.Save(tuning.Config{
		ResultsRootPath: resultsRootPath,
		Heuristics:      heuristicConfigs,
		ReadStatistics:  readTuningSummaryStatistics,
	})
}

func readTuningSummaryStatistics(path string) ([]tuning.Statistic, error) {
	experimentStatistics, err := experiments.ReadStatistics(path)
	if err != nil {
		return nil, err
	}

	statistics := make([]tuning.Statistic, 0, len(experimentStatistics))
	for _, statistic := range experimentStatistics {
		statistics = append(statistics, tuning.Statistic{
			HeuristicWeight:      statistic.HeuristicWeight,
			MsaPatchBias:         statistic.MsaPatchBias,
			AverageBestDeviation: statistic.AverageBestDeviation,
			SuccessRate:          statistic.SuccessRate,
		})
	}
	return statistics, nil
}

func runExperimentSet(atspsData []AtspData, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments, workers int) error {
	workers = max(1, workers)

	instanceRuns := make([]*experimentSetInstanceRun, 0, len(atspsData))
	for _, atspData := range atspsData {
		instanceRun, err := prepareExperimentSetInstance(atspData, heuristic, experimentParameters, numberOfExperiments, workers)
		if err != nil {
			return err
		}
		instanceRuns = append(instanceRuns, instanceRun)
	}

	var logMutex sync.Mutex
	jobs := experimentSetParameterJobs(instanceRuns, experimentParameters, numberOfExperiments, &logMutex)
	if err := workerpool.RunJobs(jobs, workers); err != nil {
		return err
	}

	for _, instanceRun := range instanceRuns {
		if err := finishExperimentSetInstance(instanceRun); err != nil {
			return err
		}
	}

	return nil
}

type experimentSetInstanceRun struct {
	atspData           AtspData
	heuristic          string
	matrix             [][]float64
	knownOptimal       float64
	dimension          int
	instanceStart      time.Time
	heuristicMatrix    [][]float64
	rootedMsaHeuristic [][][]float64
	cycleCover         [][]float64
	experimentData     []ExperimentsData
}

func prepareExperimentSetInstance(atspData AtspData, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments, workers int) (*experimentSetInstanceRun, error) {
	matrix := atspData.Matrix
	knownOptimal := atspData.KnownOptimal
	dimension := len(matrix)
	instanceStart := time.Now()

	instanceRun := &experimentSetInstanceRun{
		atspData:       atspData,
		heuristic:      heuristic,
		matrix:         matrix,
		knownOptimal:   knownOptimal,
		dimension:      dimension,
		instanceStart:  instanceStart,
		experimentData: make([]ExperimentsData, len(experimentParameters)),
	}

	var err error
	if heuristicUsesRootedMsa(heuristic) {
		instanceRun.rootedMsaHeuristic, err = readRootedMsaHeuristics(atspData)
		if err != nil {
			return nil, err
		}
	} else if heuristicUsesMsaHeuristic(heuristic) {
		instanceRun.heuristicMatrix, err = readMsaHeuristicMatrix(atspData)
		if err != nil {
			return nil, err
		}
	}

	if heuristicUsesCycleCover(heuristic) {
		instanceRun.cycleCover, err = cyclecover.Read(atspData.CycleCoverDirectoryPath, dimension)
		if err != nil {
			return nil, err
		}
		cycleCoverCost := cyclecover.Cost(matrix, instanceRun.cycleCover)
		fmt.Printf("\t[%s][%s] Minimum cycle cover cost=%.2f gap=%.2f%%\n",
			atspData.Name,
			heuristic,
			cycleCoverCost,
			100*(knownOptimal-cycleCoverCost)/knownOptimal)
	}

	fmt.Printf("[%s][%s][%s] Starting (dimension=%d, parameters=%d, runs/parameter=%d, workers=%d)\n",
		logTimestamp(instanceStart),
		atspData.Name,
		heuristic,
		dimension,
		len(experimentParameters),
		numberOfExperiments,
		workers)

	return instanceRun, nil
}

func experimentSetParameterJobs(instanceRuns []*experimentSetInstanceRun, experimentParameters []ExperimentParameters, numberOfExperiments int, logMutex *sync.Mutex) []workerpool.Job {
	orderedInstanceRuns := append([]*experimentSetInstanceRun{}, instanceRuns...)
	sort.SliceStable(orderedInstanceRuns, func(i, j int) bool {
		if orderedInstanceRuns[i].dimension != orderedInstanceRuns[j].dimension {
			return orderedInstanceRuns[i].dimension > orderedInstanceRuns[j].dimension
		}

		return orderedInstanceRuns[i].atspData.Name < orderedInstanceRuns[j].atspData.Name
	})

	jobs := make([]workerpool.Job, 0, len(orderedInstanceRuns)*len(experimentParameters))
	for _, instanceRun := range orderedInstanceRuns {
		instanceRun := instanceRun
		for parameterIndex := range experimentParameters {
			parameterIndex := parameterIndex
			jobs = append(jobs, workerpool.Job{
				Label: fmt.Sprintf("%s/%s/parameter-%d", instanceRun.heuristic, instanceRun.atspData.Name, parameterIndex),
				Run: func() error {
					return runExperimentSetParameter(instanceRun, experimentParameters[parameterIndex], parameterIndex, numberOfExperiments, logMutex)
				},
			})
		}
	}

	return jobs
}

func runExperimentSetParameter(instanceRun *experimentSetInstanceRun, parameters ExperimentParameters, parameterIndex, numberOfExperiments int, logMutex *sync.Mutex) error {
	setDimensionDependentParameters(instanceRun.dimension, &parameters)
	parameterStart := time.Now()

	heuristicModifiers := buildHeuristicModifiers(instanceRun.heuristic, instanceRun.matrix, instanceRun.heuristicMatrix, instanceRun.cycleCover, parameters)
	rootedHeuristicModifiers := buildRootedHeuristicModifiers(instanceRun.heuristic, instanceRun.matrix, instanceRun.rootedMsaHeuristic, parameters)
	results := experiments.RunExperimentsWithRootedHeuristicModifiers(numberOfExperiments, parameters, instanceRun.knownOptimal, instanceRun.matrix, heuristicModifiers, rootedHeuristicModifiers, false)
	data := ExperimentsData{ExperimentParameters: parameters, Results: results}
	instanceRun.experimentData[parameterIndex] = data

	parameterStatistics := experiments.CalculateStatistics([]ExperimentsData{data})
	if len(parameterStatistics) == 0 {
		return nil
	}

	statistic := parameterStatistics[0]
	randomSeedLog := ""
	if parameters.RandomSeed != 0 {
		randomSeedLog = fmt.Sprintf(" randomSeed=%d", parameters.RandomSeed)
	}
	logMutex.Lock()
	fmt.Printf("\t[%s][%s] heuristicWeight=%.2f msaPatchBias=%.2f%s iterations=%d runs=%d elapsed=%s min deviation=%.2f avg deviation=%.2f\n",
		instanceRun.atspData.Name,
		instanceRun.heuristic,
		parameters.HeuristicWeight,
		parameters.MsaPatchBias,
		randomSeedLog,
		parameters.Iterations,
		numberOfExperiments,
		time.Since(parameterStart).Round(time.Millisecond),
		statistic.MinBestDeviation,
		statistic.AverageBestDeviation)
	logMutex.Unlock()
	return nil
}

func finishExperimentSetInstance(instanceRun *experimentSetInstanceRun) error {
	atspData := instanceRun.atspData
	heuristic := instanceRun.heuristic
	experimentData := instanceRun.experimentData

	statistics := experiments.CalculateStatistics(experimentData)
	if len(statistics) != 0 {
		experiments.SaveStatistics(resultFilePathForHeuristic(atspData, heuristic), heuristic, statistics)
		if !heuristicIsSparseControl(heuristic) {
			if err := removeExperimentPlotsForHeuristic(atspData, heuristic); err != nil {
				return err
			}
			saveExperimentPlots(statistics, "MMAS deviation per iteration", resultPlotFilePrefixForHeuristic(atspData, heuristic))
		}
	}

	if !heuristicIsSparseControl(heuristic) {
		uniqueOptimalTours, err := tours.ReadOptimal(atspData.OptimalUniqueToursCsvPath)
		if err != nil {
			return err
		}

		for _, data := range experimentData {
			for _, result := range data.Results {
				if result.DeviationPerIteration[result.BestAtIteration] == 0.0 {
					tours.AddUnique(uniqueOptimalTours, result.BestTour)
				}
			}
		}

		if err := tours.SaveStatistics(atspData.OptimalUniqueToursCsvPath, atspData.MsaHeuristicDirectoryPath, uniqueOptimalTours); err != nil {
			return err
		}
	}

	fmt.Printf("[%s][%s][%s] Finished in %s\n", logTimestamp(time.Now()), atspData.Name, heuristic, time.Since(instanceRun.instanceStart).Round(time.Millisecond))
	return nil
}
