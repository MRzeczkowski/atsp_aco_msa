package app

import (
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	"atsp_aco_msa/modules/analysis/msaHeuristicTours"
	"atsp_aco_msa/modules/analysis/tuningSummary"
	"atsp_aco_msa/modules/project"
	"fmt"
	"os"
	"path/filepath"
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

	return runBoundedInstanceJobs(atspsData, workers, func(atspData AtspData) error {
		finalAtspData := project.WithExperimentOutputRoot(atspData, resultsRootPath)
		return runFinalExperimentForInstance(finalAtspData, resultsRootPath, useThreeOpt, configurations, numberOfExperiments)
	})
}

func runFinalExperimentForInstance(atspData AtspData, resultsRootPath string, useThreeOpt bool, configurations []finalExperimentConfiguration, numberOfExperiments int) error {
	return runFinalExperimentForInstanceWithParameterWorkers(atspData, resultsRootPath, useThreeOpt, configurations, numberOfExperiments, 1)
}

func runFinalExperimentForInstanceWithParameterWorkers(atspData AtspData, resultsRootPath string, useThreeOpt bool, configurations []finalExperimentConfiguration, numberOfExperiments, parameterWorkers int) error {
	matrix := atspData.Matrix
	knownOptimal := atspData.KnownOptimal
	dimension := len(matrix)
	instanceStart := time.Now()
	controlRun := finalConfigurationsAreSparseControls(configurations)
	finalRunName := "final"
	if useThreeOpt {
		finalRunName = "final+3opt"
	}

	if len(configurations) == 0 {
		return fmt.Errorf("no final experiment configurations selected")
	}
	if parameterWorkers < 1 {
		return fmt.Errorf("parameter workers must be at least one")
	}

	if controlRun {
		if err := removeFileIfExists(atspData.ResultFilePath); err != nil {
			return err
		}
	} else {
		if err := removeLegacyFinalResultFiles(atspData); err != nil {
			return err
		}
	}

	fmt.Printf("[%s][%s] Starting %s (dimension=%d, heuristics=%d, parameter sets=%d, runs/parameter=%d)\n",
		logTimestamp(instanceStart),
		atspData.Name,
		finalRunName,
		dimension,
		len(configurations),
		finalExperimentParameterSetCount(configurations),
		numberOfExperiments)

	var cycleCover [][]float64
	var cycleCoverErr error
	cycleCoverReady := false
	finalStatistics := make([]HeuristicExperimentStatistics, 0, len(configurations))

	for _, config := range configurations {
		var heuristicMatrix [][]float64
		var rootedMsaHeuristic [][][]float64
		if heuristicUsesRootedMsa(config.Heuristic) {
			var err error
			rootedMsaHeuristic, err = readRootedMsaHeuristics(atspData)
			if err != nil {
				return err
			}
		} else if heuristicUsesMsaHeuristic(config.Heuristic) {
			var err error
			heuristicMatrix, err = readMsaHeuristicMatrixForResultRoot(atspData, config.Heuristic, resultsRootPath)
			if err != nil {
				return err
			}
		}

		if heuristicUsesCycleCover(config.Heuristic) && !cycleCoverReady {
			var cycleCoverCost float64
			cycleCover, cycleCoverCost, cycleCoverErr = buildMinimumCycleCoverMatrix(matrix)
			if cycleCoverErr != nil {
				return cycleCoverErr
			}
			cycleCoverReady = true
			fmt.Printf("\t[%s] Minimum cycle cover cost=%.2f gap=%.2f%%\n",
				atspData.Name,
				cycleCoverCost,
				100*(knownOptimal-cycleCoverCost)/knownOptimal)
		}

		experimentData, err := runFinalExperimentParameters(
			atspData.Name,
			config.Heuristic,
			config.Parameters,
			numberOfExperiments,
			knownOptimal,
			matrix,
			heuristicMatrix,
			rootedMsaHeuristic,
			cycleCover,
			useThreeOpt,
			dimension,
			parameterWorkers)
		if err != nil {
			return err
		}

		statistics := calculateStatistics(experimentData)
		if len(statistics) == 0 {
			continue
		}

		if controlRun {
			saveStatistics(resultFilePathForHeuristic(atspData, config.Heuristic), config.Heuristic, statistics)
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

		if err := removeExperimentPlotsForHeuristic(atspData, config.Heuristic); err != nil {
			return err
		}
	}

	if !controlRun {
		if err := saveFinalHeuristicStatistics(atspData.ResultFilePath, finalStatistics, configurations); err != nil {
			return err
		}
	}

	fmt.Printf("[%s][%s] Finished %s in %s\n", logTimestamp(time.Now()), atspData.Name, finalRunName, time.Since(instanceStart).Round(time.Millisecond))
	return nil
}

func runFinalExperimentParameters(instanceName, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments int, knownOptimal float64, matrix, heuristicMatrix [][]float64, rootedMsaHeuristic [][][]float64, cycleCover [][]float64, useThreeOpt bool, dimension, workers int) ([]ExperimentsData, error) {
	experimentData := make([]ExperimentsData, len(experimentParameters))
	var logMutex sync.Mutex

	err := runBoundedIndexJobs(len(experimentParameters), workers, func(index int) error {
		parameters := experimentParameters[index]
		setDimensionDependantParameters(dimension, &parameters)
		parameterStart := time.Now()
		heuristicModifiers := buildHeuristicModifiers(heuristic, matrix, heuristicMatrix, cycleCover, parameters)
		rootedHeuristicModifiers := buildRootedHeuristicModifiers(heuristic, matrix, rootedMsaHeuristic, parameters)
		results := runExperimentsWithRootedHeuristicModifiers(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, rootedHeuristicModifiers, useThreeOpt)
		data := ExperimentsData{parameters, results}
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
	})
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
	workerGate := make(chan struct{}, workers)

	return runBoundedInstanceJobs(atspsData, min(workers, maxIntValue(1, len(atspsData))), func(atspData AtspData) error {
		return runExperimentSetForInstance(atspData, heuristic, experimentParameters, numberOfExperiments, workerGate)
	})
}

func runExperimentSetForInstance(atspData AtspData, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments int, workerGate chan struct{}) error {
	matrix := atspData.Matrix
	knownOptimal := atspData.KnownOptimal
	dimension := len(matrix)
	instanceStart := time.Now()

	var heuristicMatrix [][]float64
	var rootedMsaHeuristic [][][]float64
	var err error
	if heuristicUsesRootedMsaForTuning(heuristic) {
		rootedMsaHeuristic, err = readRootedMsaHeuristics(atspData)
		if err != nil {
			return err
		}
	} else if heuristicUsesMsaHeuristic(heuristic) {
		heuristicMatrix, err = readMsaHeuristicMatrixForHeuristic(atspData, heuristic)
		if err != nil {
			return err
		}
	}

	fmt.Printf("[%s][%s][%s] Starting (dimension=%d, parameters=%d, runs/parameter=%d, parameter workers=%d)\n",
		logTimestamp(instanceStart),
		atspData.Name,
		heuristic,
		dimension,
		len(experimentParameters),
		numberOfExperiments,
		min(cap(workerGate), len(experimentParameters)))

	var cycleCover [][]float64
	if heuristicUsesCycleCover(heuristic) {
		var cycleCoverCost float64
		cycleCover, cycleCoverCost, err = buildMinimumCycleCoverMatrix(matrix)
		if err != nil {
			return err
		}
		fmt.Printf("\t[%s][%s] Minimum cycle cover cost=%.2f gap=%.2f%%\n",
			atspData.Name,
			heuristic,
			cycleCoverCost,
			100*(knownOptimal-cycleCoverCost)/knownOptimal)
	}

	experimentData := make([]ExperimentsData, len(experimentParameters))
	var logMutex sync.Mutex

	err = runIndexJobsWithSharedWorkers(len(experimentParameters), workerGate, func(index int) error {
		parameters := experimentParameters[index]
		setDimensionDependantParameters(dimension, &parameters)
		parameterStart := time.Now()
		heuristicModifiers := buildHeuristicModifiers(heuristic, matrix, heuristicMatrix, cycleCover, parameters)
		rootedHeuristicModifiers := buildRootedHeuristicModifiers(heuristic, matrix, rootedMsaHeuristic, parameters)
		results := runExperimentsWithRootedHeuristicModifiers(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, rootedHeuristicModifiers, false)
		data := ExperimentsData{parameters, results}
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
				atspData.Name,
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
	})
	if err != nil {
		return err
	}

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

	fmt.Printf("[%s][%s][%s] Finished in %s\n", logTimestamp(time.Now()), atspData.Name, heuristic, time.Since(instanceStart).Round(time.Millisecond))
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
