package app

import (
	"atsp_aco_msa/modules/artifacts/cyclecover"
	"atsp_aco_msa/modules/experiments"
	workerpool "atsp_aco_msa/modules/experiments/workers"
	"atsp_aco_msa/modules/project"
	"fmt"
	"sort"
	"sync"
	"time"
)

func runEvaluationExperimentMode(atspsData []AtspData, resultsRootPath string, useThreeOpt bool, configurations []evaluationExperimentConfiguration, workers, numberOfExperiments int) error {
	needsMsaHeuristic := evaluationConfigurationsUseMsaHeuristic(configurations)
	if needsMsaHeuristic {
		if err := ensureMsaHeuristicCache(atspsData, workers, evaluationConfigurationsUseRootedMsa(configurations)); err != nil {
			return err
		}
	}
	if evaluationConfigurationsUseCycleCover(configurations) {
		if err := ensureCycleCoverCache(atspsData, workers); err != nil {
			return err
		}
	}

	instanceRuns := make([]*evaluationExperimentInstanceRun, 0, len(atspsData))
	for _, atspData := range atspsData {
		evaluationAtspData := project.WithExperimentOutputRoot(atspData, resultsRootPath)
		instanceRun, err := prepareEvaluationExperimentInstance(evaluationAtspData, useThreeOpt, configurations, numberOfExperiments)
		if err != nil {
			return err
		}
		instanceRuns = append(instanceRuns, instanceRun)
	}

	if err := workerpool.RunJobs(evaluationExperimentConfigurationJobs(instanceRuns), workers); err != nil {
		return err
	}

	for _, instanceRun := range instanceRuns {
		if err := finishEvaluationExperimentInstance(instanceRun); err != nil {
			return err
		}
	}

	return nil
}

type evaluationExperimentInstanceRun struct {
	atspData                  AtspData
	useThreeOpt               bool
	configurations            []evaluationExperimentConfiguration
	numberOfExperiments       int
	matrix                    [][]float64
	knownOptimal              float64
	dimension                 int
	instanceStart             time.Time
	controlRun                bool
	evaluationRunName         string
	statisticsByConfiguration [][]ExperimentsDataStatistics
	cycleCover                [][]float64
}

func prepareEvaluationExperimentInstance(atspData AtspData, useThreeOpt bool, configurations []evaluationExperimentConfiguration, numberOfExperiments int) (*evaluationExperimentInstanceRun, error) {
	if len(configurations) == 0 {
		return nil, fmt.Errorf("no evaluation experiment configurations selected")
	}

	controlRun := evaluationConfigurationsAreSparseControls(configurations)

	evaluationRunName := "evaluation"
	if useThreeOpt {
		evaluationRunName = "evaluation+3opt"
	}

	instanceStart := time.Now()
	instanceRun := &evaluationExperimentInstanceRun{
		atspData:                  atspData,
		useThreeOpt:               useThreeOpt,
		configurations:            configurations,
		numberOfExperiments:       numberOfExperiments,
		matrix:                    atspData.Matrix,
		knownOptimal:              atspData.KnownOptimal,
		dimension:                 len(atspData.Matrix),
		instanceStart:             instanceStart,
		controlRun:                controlRun,
		evaluationRunName:         evaluationRunName,
		statisticsByConfiguration: make([][]ExperimentsDataStatistics, len(configurations)),
	}

	if evaluationConfigurationsUseCycleCover(configurations) {
		cycleCoverMatrix, err := cyclecover.Read(atspData.CycleCoverDirectoryPath, instanceRun.dimension)
		if err != nil {
			return nil, err
		}
		instanceRun.cycleCover = cycleCoverMatrix
		cycleCoverCost := cyclecover.Cost(instanceRun.matrix, cycleCoverMatrix)
		fmt.Printf("\t[%s] Minimum cycle cover cost=%.2f gap=%.2f%%\n",
			atspData.Name,
			cycleCoverCost,
			100*(instanceRun.knownOptimal-cycleCoverCost)/instanceRun.knownOptimal)
	}

	fmt.Printf("[%s][%s] Starting %s (dimension=%d, heuristics=%d, parameter sets=%d, runs/parameter=%d)\n",
		logTimestamp(instanceStart),
		atspData.Name,
		evaluationRunName,
		instanceRun.dimension,
		len(configurations),
		evaluationExperimentParameterSetCount(configurations),
		numberOfExperiments)

	return instanceRun, nil
}

func evaluationExperimentConfigurationJobs(instanceRuns []*evaluationExperimentInstanceRun) []workerpool.Job {
	orderedInstanceRuns := append([]*evaluationExperimentInstanceRun{}, instanceRuns...)
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
				Label: fmt.Sprintf("%s/%s/%s", instanceRun.evaluationRunName, instanceRun.atspData.Name, config.Heuristic),
				Run: func() error {
					return runEvaluationExperimentConfiguration(instanceRun, configurationIndex)
				},
			})
		}
	}

	return jobs
}

func runEvaluationExperimentConfiguration(instanceRun *evaluationExperimentInstanceRun, configurationIndex int) error {
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
		heuristicMatrix, err = readMsaHeuristicMatrix(instanceRun.atspData)
		if err != nil {
			return err
		}
	}

	var cycleCover [][]float64
	if heuristicUsesCycleCover(config.Heuristic) {
		cycleCover = instanceRun.cycleCover
	}

	experimentData, err := runEvaluationExperimentParameters(
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

	statistics := experiments.CalculateStatistics(experimentData)
	instanceRun.statisticsByConfiguration[configurationIndex] = statistics
	if len(statistics) != 0 && !instanceRun.controlRun {
		return removeExperimentPlotsForHeuristic(instanceRun.atspData, config.Heuristic)
	}

	return nil
}

func finishEvaluationExperimentInstance(instanceRun *evaluationExperimentInstanceRun) error {
	if instanceRun.controlRun {
		for index, config := range instanceRun.configurations {
			statistics := instanceRun.statisticsByConfiguration[index]
			if len(statistics) != 0 {
				experiments.SaveStatistics(resultFilePathForHeuristic(instanceRun.atspData, config.Heuristic), config.Heuristic, statistics)
			}
		}
	} else {
		evaluationStatistics := make([]HeuristicExperimentStatistics, 0, len(instanceRun.configurations))
		for index, config := range instanceRun.configurations {
			statistics := instanceRun.statisticsByConfiguration[index]
			if len(statistics) == 0 {
				continue
			}

			if config.SaveAllParameterRows {
				for _, statistic := range statistics {
					evaluationStatistics = append(evaluationStatistics, HeuristicExperimentStatistics{
						Heuristic:  config.Heuristic,
						Statistics: statistic,
					})
				}
			} else {
				evaluationStatistics = append(evaluationStatistics, HeuristicExperimentStatistics{
					Heuristic:  config.Heuristic,
					Statistics: statistics[0],
				})
			}
		}

		if err := saveEvaluationHeuristicStatistics(instanceRun.atspData.ResultFilePath, evaluationStatistics, instanceRun.configurations); err != nil {
			return err
		}
	}

	fmt.Printf("[%s][%s] Finished %s in %s\n", logTimestamp(time.Now()), instanceRun.atspData.Name, instanceRun.evaluationRunName, time.Since(instanceRun.instanceStart).Round(time.Millisecond))
	return nil
}

func runEvaluationExperimentParameters(instanceName, heuristic string, experimentParameters []ExperimentParameters, numberOfExperiments int, knownOptimal float64, matrix, heuristicMatrix [][]float64, rootedMsaHeuristic [][][]float64, cycleCover [][]float64, useThreeOpt bool, dimension, workers int) ([]ExperimentsData, error) {
	experimentData := make([]ExperimentsData, len(experimentParameters))
	var logMutex sync.Mutex

	jobs := make([]workerpool.Job, 0, len(experimentParameters))
	for index := range experimentParameters {
		index := index
		jobs = append(jobs, workerpool.Job{
			Label: fmt.Sprintf("%s/%s/parameter-%d", instanceName, heuristic, index),
			Run: func() error {
				parameters := experimentParameters[index]
				setDimensionDependentParameters(dimension, &parameters)
				parameterStart := time.Now()
				heuristicModifiers := buildHeuristicModifiers(heuristic, matrix, heuristicMatrix, cycleCover, parameters)
				rootedHeuristicModifiers := buildRootedHeuristicModifiers(heuristic, matrix, rootedMsaHeuristic, parameters)
				results := experiments.RunExperimentsWithRootedHeuristicModifiers(numberOfExperiments, parameters, knownOptimal, matrix, heuristicModifiers, rootedHeuristicModifiers, useThreeOpt)
				data := ExperimentsData{ExperimentParameters: parameters, Results: results}
				experimentData[index] = data

				parameterStatistics := experiments.CalculateStatistics([]ExperimentsData{data})
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

func evaluationExperimentParameterSetCount(configurations []evaluationExperimentConfiguration) int {
	count := 0
	for _, configuration := range configurations {
		count += len(configuration.Parameters)
	}

	return count
}
