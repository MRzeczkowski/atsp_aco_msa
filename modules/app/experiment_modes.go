package app

import (
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	"atsp_aco_msa/modules/analysis/msaHeuristicTours"
	"atsp_aco_msa/modules/analysis/tuningSummary"
	workerpool "atsp_aco_msa/modules/experiments/workers"
	"atsp_aco_msa/modules/project"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"sync"
	"time"
)

func runExperimentMode(atspsData []AtspData, heuristics []string, workers int) error {
	if err := ensureMsaHeuristicArtifacts(atspsData, workers, heuristicsUseRootedMsa(heuristics)); err != nil {
		return err
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
	heuristicConfigs := make([]tuningSummary.HeuristicConfig, 0, len(heuristics))
	for _, heuristic := range heuristics {
		if heuristicIsSparseControl(heuristic) {
			continue
		}

		resultFiles := make([]tuningSummary.ResultFile, 0, len(atspsData))
		for _, atspData := range atspsData {
			resultFiles = append(resultFiles, tuningSummary.ResultFile{
				Instance: atspData.Name,
				Path:     resultFilePathForHeuristic(atspData, heuristic),
			})
		}

		heuristicConfigs = append(heuristicConfigs, tuningSummary.HeuristicConfig{
			Name:                heuristic,
			DisplayName:         heuristicDisplayName(heuristic),
			IncludeMsaPatchBias: heuristic == heuristicCycleCoverMsaPatching,
			ResultFiles:         resultFiles,
		})
	}

	return tuningSummary.Save(tuningSummary.Config{
		ResultsRootPath: resultsRootPath,
		Heuristics:      heuristicConfigs,
		ReadStatistics:  readTuningSummaryStatistics,
	})
}

func readTuningSummaryStatistics(path string) ([]tuningSummary.Statistic, error) {
	experimentStatistics, err := readStatistics(path)
	if err != nil {
		return nil, err
	}

	statistics := make([]tuningSummary.Statistic, 0, len(experimentStatistics))
	for _, statistic := range experimentStatistics {
		statistics = append(statistics, tuningSummary.Statistic{
			HeuristicWeight:      statistic.HeuristicWeight,
			MsaPatchBias:         statistic.MsaPatchBias,
			AverageBestDeviation: statistic.AverageBestDeviation,
			SuccessRate:          statistic.SuccessRate,
		})
	}
	return statistics, nil
}

func runFinalExperimentMode(atspsData []AtspData, resultsRootPath string, useThreeOpt bool, configurations []finalExperimentConfiguration, workers, numberOfExperiments int) error {
	needsMsaHeuristic := finalConfigurationsUseMsaHeuristic(configurations)
	if needsMsaHeuristic {
		if err := ensureMsaHeuristicCache(atspsData, workers, finalConfigurationsUseRootedMsa(configurations)); err != nil {
			return err
		}
	}

	if err := removeLegacyFinalReports(resultsRootPath); err != nil {
		return err
	}

	instanceRuns := make([]*finalExperimentInstanceRun, 0, len(atspsData))
	for _, atspData := range atspsData {
		finalAtspData := project.WithExperimentOutputRoot(atspData, resultsRootPath)
		instanceRun, err := prepareFinalExperimentInstance(finalAtspData, useThreeOpt, configurations, numberOfExperiments)
		if err != nil {
			return err
		}
		instanceRuns = append(instanceRuns, instanceRun)
	}

	if err := workerpool.RunJobs(finalExperimentConfigurationJobs(instanceRuns), workers); err != nil {
		return err
	}

	for _, instanceRun := range instanceRuns {
		if err := finishFinalExperimentInstance(instanceRun); err != nil {
			return err
		}
	}

	return nil
}

type finalExperimentInstanceRun struct {
	atspData                  AtspData
	useThreeOpt               bool
	configurations            []finalExperimentConfiguration
	numberOfExperiments       int
	matrix                    [][]float64
	knownOptimal              float64
	dimension                 int
	instanceStart             time.Time
	controlRun                bool
	finalRunName              string
	statisticsByConfiguration [][]ExperimentsDataStatistics
	cycleCoverOnce            sync.Once
	cycleCover                [][]float64
	cycleCoverErr             error
}

func prepareFinalExperimentInstance(atspData AtspData, useThreeOpt bool, configurations []finalExperimentConfiguration, numberOfExperiments int) (*finalExperimentInstanceRun, error) {
	if len(configurations) == 0 {
		return nil, fmt.Errorf("no final experiment configurations selected")
	}

	controlRun := finalConfigurationsAreSparseControls(configurations)
	if controlRun {
		if err := removeFileIfExists(atspData.ResultFilePath); err != nil {
			return nil, err
		}
	} else {
		if err := removeLegacyFinalResultFiles(atspData); err != nil {
			return nil, err
		}
	}

	finalRunName := "final"
	if useThreeOpt {
		finalRunName = "final+3opt"
	}

	instanceStart := time.Now()
	instanceRun := &finalExperimentInstanceRun{
		atspData:                  atspData,
		useThreeOpt:               useThreeOpt,
		configurations:            configurations,
		numberOfExperiments:       numberOfExperiments,
		matrix:                    atspData.Matrix,
		knownOptimal:              atspData.KnownOptimal,
		dimension:                 len(atspData.Matrix),
		instanceStart:             instanceStart,
		controlRun:                controlRun,
		finalRunName:              finalRunName,
		statisticsByConfiguration: make([][]ExperimentsDataStatistics, len(configurations)),
	}

	fmt.Printf("[%s][%s] Starting %s (dimension=%d, heuristics=%d, parameter sets=%d, runs/parameter=%d)\n",
		logTimestamp(instanceStart),
		atspData.Name,
		finalRunName,
		instanceRun.dimension,
		len(configurations),
		finalExperimentParameterSetCount(configurations),
		numberOfExperiments)

	return instanceRun, nil
}

func finalExperimentConfigurationJobs(instanceRuns []*finalExperimentInstanceRun) []workerpool.Job {
	orderedInstanceRuns := append([]*finalExperimentInstanceRun{}, instanceRuns...)
	sort.SliceStable(orderedInstanceRuns, func(i, j int) bool {
		if orderedInstanceRuns[i].dimension != orderedInstanceRuns[j].dimension {
			return orderedInstanceRuns[i].dimension > orderedInstanceRuns[j].dimension
		}

		return orderedInstanceRuns[i].atspData.Name < orderedInstanceRuns[j].atspData.Name
	})

	jobCount := 0
	for _, instanceRun := range orderedInstanceRuns {
		jobCount += len(instanceRun.configurations)
	}

	jobs := make([]workerpool.Job, 0, jobCount)
	for _, instanceRun := range orderedInstanceRuns {
		instanceRun := instanceRun
		for configurationIndex := range instanceRun.configurations {
			configurationIndex := configurationIndex
			config := instanceRun.configurations[configurationIndex]
			jobs = append(jobs, workerpool.Job{
				Label: fmt.Sprintf("%s/%s/%s", instanceRun.finalRunName, instanceRun.atspData.Name, config.Heuristic),
				Run: func() error {
					return runFinalExperimentConfiguration(instanceRun, configurationIndex)
				},
			})
		}
	}

	return jobs
}

func runFinalExperimentConfiguration(instanceRun *finalExperimentInstanceRun, configurationIndex int) error {
	config := instanceRun.configurations[configurationIndex]
	var heuristicMatrix [][]float64
	var rootedMsaHeuristic [][][]float64
	if heuristicUsesRootedMsa(config.Heuristic) {
		var err error
		rootedMsaHeuristic, err = readRootedMsaHeuristics(instanceRun.atspData)
		if err != nil {
			return err
		}
	} else if heuristicUsesMsaHeuristic(config.Heuristic) {
		var err error
		heuristicMatrix, err = readMsaHeuristicMatrixForHeuristic(instanceRun.atspData, config.Heuristic)
		if err != nil {
			return err
		}
	}

	var cycleCover [][]float64
	if heuristicUsesCycleCover(config.Heuristic) {
		var err error
		cycleCover, err = instanceRun.getCycleCover()
		if err != nil {
			return err
		}
	}

	experimentData, err := runFinalExperimentParameters(
		instanceRun.atspData.Name,
		config.Heuristic,
		config.Parameters,
		instanceRun.numberOfExperiments,
		instanceRun.knownOptimal,
		instanceRun.matrix,
		heuristicMatrix,
		rootedMsaHeuristic,
		cycleCover,
		instanceRun.useThreeOpt,
		instanceRun.dimension,
		1)
	if err != nil {
		return err
	}

	statistics := calculateStatistics(experimentData)
	instanceRun.statisticsByConfiguration[configurationIndex] = statistics
	if len(statistics) != 0 && !instanceRun.controlRun {
		return removeExperimentPlotsForHeuristic(instanceRun.atspData, config.Heuristic)
	}

	return nil
}

func (instanceRun *finalExperimentInstanceRun) getCycleCover() ([][]float64, error) {
	instanceRun.cycleCoverOnce.Do(func() {
		var cycleCoverCost float64
		instanceRun.cycleCover, cycleCoverCost, instanceRun.cycleCoverErr = buildMinimumCycleCoverMatrix(instanceRun.matrix)
		if instanceRun.cycleCoverErr != nil {
			return
		}

		fmt.Printf("\t[%s] Minimum cycle cover cost=%.2f gap=%.2f%%\n",
			instanceRun.atspData.Name,
			cycleCoverCost,
			100*(instanceRun.knownOptimal-cycleCoverCost)/instanceRun.knownOptimal)
	})

	return instanceRun.cycleCover, instanceRun.cycleCoverErr
}

func finishFinalExperimentInstance(instanceRun *finalExperimentInstanceRun) error {
	if instanceRun.controlRun {
		for index, config := range instanceRun.configurations {
			statistics := instanceRun.statisticsByConfiguration[index]
			if len(statistics) != 0 {
				saveStatistics(resultFilePathForHeuristic(instanceRun.atspData, config.Heuristic), config.Heuristic, statistics)
			}
		}
	} else {
		finalStatistics := make([]HeuristicExperimentStatistics, 0, len(instanceRun.configurations))
		for index, config := range instanceRun.configurations {
			statistics := instanceRun.statisticsByConfiguration[index]
			if len(statistics) == 0 {
				continue
			}

			if config.SaveAllParameterRows {
				for _, statistic := range statistics {
					finalStatistics = append(finalStatistics, HeuristicExperimentStatistics{
						Heuristic:  config.Heuristic,
						Statistics: statistic,
					})
				}
			} else {
				finalStatistics = append(finalStatistics, HeuristicExperimentStatistics{
					Heuristic:  config.Heuristic,
					Statistics: statistics[0],
				})
			}
		}

		if err := saveFinalHeuristicStatistics(instanceRun.atspData.ResultFilePath, finalStatistics, instanceRun.configurations); err != nil {
			return err
		}
	}

	fmt.Printf("[%s][%s] Finished %s in %s\n", logTimestamp(time.Now()), instanceRun.atspData.Name, instanceRun.finalRunName, time.Since(instanceRun.instanceStart).Round(time.Millisecond))
	return nil
}

func runFinalExperimentParameters(instanceName, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments int, knownOptimal float64, matrix, heuristicMatrix [][]float64, rootedMsaHeuristic [][][]float64, cycleCover [][]float64, useThreeOpt bool, dimension, workers int) ([]ExperimentsData, error) {
	experimentData := make([]ExperimentsData, len(experimentParameters))
	var logMutex sync.Mutex

	jobs := make([]workerpool.Job, 0, len(experimentParameters))
	for index := range experimentParameters {
		index := index
		jobs = append(jobs, workerpool.Job{
			Label: fmt.Sprintf("%s/%s/parameter-%d", instanceName, heuristic, index),
			Run: func() error {
				parameters := experimentParameters[index]
				setDimensionDependantParameters(dimension, &parameters)
				parameterStart := time.Now()
				heuristicModifiers := buildHeuristicModifiers(heuristic, matrix, heuristicMatrix, cycleCover, parameters)
				rootedHeuristicModifiers := buildRootedHeuristicModifiers(heuristic, matrix, rootedMsaHeuristic, parameters)
				results := runExperimentsWithRootedHeuristicModifiers(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, rootedHeuristicModifiers, useThreeOpt)
				data := ExperimentsData{ExperimentParameters: parameters, Results: results}
				experimentData[index] = data

				parameterStatistics := calculateStatistics([]ExperimentsData{data})
				if len(parameterStatistics) != 0 {
					statistic := parameterStatistics[0]
					randomSeedLog := ""
					if parameters.RandomSeed != 0 {
						randomSeedLog = fmt.Sprintf(" randomSeed=%d", parameters.RandomSeed)
					}
					logMutex.Lock()
					fmt.Printf("\t[%s][%s] heuristicWeight=%.2f msaPatchBias=%.2f%s iterations=%d runs=%d elapsed=%s min deviation=%.2f avg deviation=%.2f\n",
						instanceName,
						heuristic,
						parameters.HeuristicWeight,
						parameters.MsaPatchBias,
						randomSeedLog,
						parameters.Iterations,
						numberOfExperiments,
						time.Since(parameterStart).Round(time.Millisecond),
						statistic.MinBestDeviation,
						statistic.AverageBestDeviation)
					logMutex.Unlock()
				}

				return nil
			},
		})
	}

	err := workerpool.RunJobs(jobs, workers)
	if err != nil {
		return nil, err
	}

	return experimentData, nil
}

func finalExperimentParameterSetCount(configurations []finalExperimentConfiguration) int {
	count := 0
	for _, configuration := range configurations {
		count += len(configuration.Parameters)
	}

	return count
}

func removeLegacyFinalReports(resultsRootPath string) error {
	legacyReports, err := filepath.Glob(filepath.Join(resultsRootPath, "best_parameters_report*.md"))
	if err != nil {
		return err
	}

	for _, path := range legacyReports {
		if err := removeFileIfExists(path); err != nil {
			return err
		}
	}

	return nil
}

func removeLegacyFinalResultFiles(atspData AtspData) error {
	legacyFiles := []string{
		resultFilePathForHeuristic(atspData, heuristicBaseline),
		resultFilePathForHeuristic(atspData, heuristicCycleCover),
		resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching),
		filepath.Join(filepath.Dir(atspData.ResultFilePath), "solutions.csv"),
	}

	for _, path := range legacyFiles {
		if err := removeFileIfExists(path); err != nil {
			return err
		}
	}

	return nil
}

func removeFileIfExists(path string) error {
	if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
		return err
	}

	return nil
}

func runExperimentSet(atspsData []AtspData, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments, workers int) error {
	workers = maxIntValue(1, workers)

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
	cycleCoverOnce     sync.Once
	cycleCoverErr      error
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
	if heuristicUsesRootedMsaForTuning(heuristic) {
		instanceRun.rootedMsaHeuristic, err = readRootedMsaHeuristics(atspData)
		if err != nil {
			return nil, err
		}
	} else if heuristicUsesMsaHeuristic(heuristic) {
		instanceRun.heuristicMatrix, err = readMsaHeuristicMatrixForHeuristic(atspData, heuristic)
		if err != nil {
			return nil, err
		}
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

func (instanceRun *experimentSetInstanceRun) getCycleCover() ([][]float64, error) {
	instanceRun.cycleCoverOnce.Do(func() {
		var cycleCoverCost float64
		instanceRun.cycleCover, cycleCoverCost, instanceRun.cycleCoverErr = buildMinimumCycleCoverMatrix(instanceRun.matrix)
		if instanceRun.cycleCoverErr != nil {
			return
		}

		fmt.Printf("\t[%s][%s] Minimum cycle cover cost=%.2f gap=%.2f%%\n",
			instanceRun.atspData.Name,
			instanceRun.heuristic,
			cycleCoverCost,
			100*(instanceRun.knownOptimal-cycleCoverCost)/instanceRun.knownOptimal)
	})

	return instanceRun.cycleCover, instanceRun.cycleCoverErr
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
	setDimensionDependantParameters(instanceRun.dimension, &parameters)
	parameterStart := time.Now()

	var cycleCover [][]float64
	if heuristicUsesCycleCover(instanceRun.heuristic) {
		var err error
		cycleCover, err = instanceRun.getCycleCover()
		if err != nil {
			return err
		}
	}

	heuristicModifiers := buildHeuristicModifiers(instanceRun.heuristic, instanceRun.matrix, instanceRun.heuristicMatrix, cycleCover, parameters)
	rootedHeuristicModifiers := buildRootedHeuristicModifiers(instanceRun.heuristic, instanceRun.matrix, instanceRun.rootedMsaHeuristic, parameters)
	results := runExperimentsWithRootedHeuristicModifiers(numberOfExperiments, parameters, instanceRun.knownOptimal, instanceRun.matrix, heuristicModifiers, rootedHeuristicModifiers, false)
	data := ExperimentsData{ExperimentParameters: parameters, Results: results}
	instanceRun.experimentData[parameterIndex] = data

	parameterStatistics := calculateStatistics([]ExperimentsData{data})
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

	statistics := calculateStatistics(experimentData)
	if len(statistics) != 0 {
		saveStatistics(resultFilePathForHeuristic(atspData, heuristic), heuristic, statistics)
		if !heuristicIsSparseControl(heuristic) {
			if err := removeExperimentPlotsForHeuristic(atspData, heuristic); err != nil {
				return err
			}
			saveExperimentPlots(statistics, "MMAS deviation per iteration", resultPlotFilePrefixForHeuristic(atspData, heuristic))
		}
	}

	if !heuristicIsSparseControl(heuristic) {
		uniqueOptimalTours, err := msaHeuristicTours.ReadOptimalTours(atspData.OptimalUniqueToursCsvPath)
		if err != nil {
			return err
		}

		for _, data := range experimentData {
			for _, result := range data.Results {
				if result.DeviationPerIteration[result.BestAtIteration] == 0.0 {
					msaHeuristicTours.AddUniqueTour(uniqueOptimalTours, result.BestTour)
				}
			}
		}

		if err := msaHeuristicTours.SaveOptimalToursStatistics(atspData.OptimalUniqueToursCsvPath, atspData.MsaHeuristicDirectoryPath, uniqueOptimalTours); err != nil {
			return err
		}
	}

	fmt.Printf("[%s][%s][%s] Finished in %s\n", logTimestamp(time.Now()), atspData.Name, heuristic, time.Since(instanceRun.instanceStart).Round(time.Millisecond))
	return nil
}

func readMsaHeuristicMatrixForHeuristic(atspData AtspData, heuristic string) ([][]float64, error) {
	return msaHeuristic.Read(atspData.MsaHeuristicDirectoryPath)
}

func readMsaHeuristicMatrixForResultRoot(atspData AtspData, heuristic, resultsRootPath string) ([][]float64, error) {
	return readMsaHeuristicMatrixForHeuristic(atspData, heuristic)
}

func readRootedMsaHeuristics(atspData AtspData) ([][][]float64, error) {
	rootedMsas, err := msaHeuristic.ReadMsas(atspData.MsaHeuristicDirectoryPath)
	if err != nil {
		return nil, err
	}
	if len(rootedMsas) != len(atspData.Matrix) {
		return nil, fmt.Errorf("expected %d rooted MSAs, got %d", len(atspData.Matrix), len(rootedMsas))
	}

	return rootedMsas, nil
}
