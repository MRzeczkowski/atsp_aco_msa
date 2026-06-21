package app

import (
	"atsp_aco_msa/modules/project"
	"flag"
	"fmt"
	"path/filepath"
	"runtime"
	"strings"
	"time"
)

func isValidRunMode(mode string) bool {
	return mode == runModeExperiment || mode == runModeAnalyze || mode == runModeAll || mode == runModeEvaluation || mode == runModeEvaluation3Opt || mode == runModeRebuildCache
}

func isValidHeuristic(heuristic string) bool {
	return heuristic == heuristicStrictMsa ||
		heuristic == heuristicRootedMsa ||
		heuristic == heuristicCycleCover ||
		heuristic == heuristicCycleCoverMsaPatching
}

func selectExperimentHeuristics(heuristic string, heuristicExplicit bool) ([]string, error) {
	if !heuristicExplicit || heuristic == "" || heuristic == experimentHeuristicAll {
		return append([]string{}, experimentHeuristics...), nil
	}

	if !isValidHeuristic(heuristic) {
		return nil, fmt.Errorf("unsupported -heuristic value %q; omit it or use %q, %q, %q, %q, or %q", heuristic, experimentHeuristicAll, heuristicStrictMsa, heuristicRootedMsa, heuristicCycleCover, heuristicCycleCoverMsaPatching)
	}

	return []string{heuristic}, nil
}

func heuristicUsesCycleCover(heuristic string) bool {
	return heuristic == heuristicCycleCover || heuristic == heuristicCycleCoverMsaPatching
}

func heuristicUsesMsaHeuristic(heuristic string) bool {
	return heuristic == heuristicStrictMsa || heuristic == heuristicRootedMsa || heuristic == heuristicCycleCoverMsaPatching || heuristicIsSparseControl(heuristic)
}

func heuristicUsesRootedMsa(heuristic string) bool {
	return heuristic == heuristicRootedMsa || heuristicIsRootedSparseControl(heuristic)
}

func heuristicUsesRandomSeedColumn(heuristic string) bool {
	return heuristic == heuristicRandomSparse || heuristicIsShuffledMsaControl(heuristic)
}

func heuristicIsSparseControl(heuristic string) bool {
	return heuristic == heuristicRandomSparse || heuristicIsDistanceRankedSparseControl(heuristic) || heuristicIsShuffledMsaControl(heuristic)
}

func heuristicIsDistanceRankedSparseControl(heuristic string) bool {
	return heuristic == heuristicDistanceRankedSparse ||
		heuristic == heuristicStrictDistanceRanked ||
		heuristic == heuristicRootedDistanceRanked
}

func heuristicIsShuffledMsaControl(heuristic string) bool {
	return heuristic == heuristicShuffledMsa ||
		heuristic == heuristicStrictShuffledMsa ||
		heuristic == heuristicRootedShuffledMsa
}

func heuristicIsRootedSparseControl(heuristic string) bool {
	return heuristic == heuristicRootedDistanceRanked || heuristic == heuristicRootedShuffledMsa
}

func evaluationConfigurationsUseMsaHeuristic(configurations []evaluationExperimentConfiguration) bool {
	for _, configuration := range configurations {
		if heuristicUsesMsaHeuristic(configuration.Heuristic) {
			return true
		}
	}

	return false
}

func evaluationConfigurationsUseCycleCover(configurations []evaluationExperimentConfiguration) bool {
	for _, configuration := range configurations {
		if heuristicUsesCycleCover(configuration.Heuristic) {
			return true
		}
	}

	return false
}

func evaluationConfigurationsUseRootedMsa(configurations []evaluationExperimentConfiguration) bool {
	for _, configuration := range configurations {
		if heuristicUsesRootedMsa(configuration.Heuristic) {
			return true
		}
	}

	return false
}

func heuristicsUseCycleCover(heuristics []string) bool {
	for _, heuristic := range heuristics {
		if heuristicUsesCycleCover(heuristic) {
			return true
		}
	}

	return false
}

func heuristicsUseRootedMsa(heuristics []string) bool {
	for _, heuristic := range heuristics {
		if heuristicUsesRootedMsa(heuristic) {
			return true
		}
	}

	return false
}

func heuristicFileSuffix(heuristic string) string {
	switch heuristic {
	case heuristicBaseline:
		return "_baseline"
	case heuristicStrictMsa:
		return ""
	case heuristicCycleCover:
		return "_cycle_cover"
	case heuristicCycleCoverMsaPatching:
		return "_cycle_cover_msa_patching"
	default:
		return "_" + strings.ReplaceAll(heuristic, "-", "_")
	}
}

func resultFilePathForHeuristic(atspData AtspData, heuristic string) string {
	suffix := heuristicFileSuffix(heuristic)
	if suffix == "" {
		return atspData.ResultFilePath
	}

	return strings.TrimSuffix(atspData.ResultFilePath, ".csv") + suffix + ".csv"
}

func resultPlotFilePrefixForHeuristic(atspData AtspData, heuristic string) string {
	return atspData.ResultPlotFilePrefix + heuristicFileSuffix(heuristic)
}

func shouldRunExperiments(mode string) bool {
	return mode == runModeExperiment || mode == runModeAll
}

func shouldRunAnalysis(mode string) bool {
	return mode == runModeAnalyze || mode == runModeAll
}

func shouldRunAnalysisAfterEvaluationExperiments(mode, evaluationHeuristic string) bool {
	return shouldRunEvaluationExperiments(mode) && evaluationHeuristic == evaluationHeuristicAll
}

func isValidAnalysisScope(scope string) bool {
	return scope == analysisScopeAll || scope == analysisScopeGksDeviation || scope == analysisScopeTuning
}

func shouldRunEvaluationExperiments(mode string) bool {
	return mode == runModeEvaluation || mode == runModeEvaluation3Opt
}

func shouldRunCacheRebuild(mode string) bool {
	return mode == runModeRebuildCache
}

func evaluationExperimentOutputRoot(mode string) string {
	if mode == runModeEvaluation3Opt {
		return project.EvaluationThreeOptResultsDirectoryName
	}

	return project.EvaluationResultsDirectoryName
}

func evaluationControlsResultsRootPath(evaluationResultsRootPath string) string {
	return filepath.Join(filepath.Dir(evaluationResultsRootPath), "controls")
}

func evaluationExperimentOutputRootForConfigurations(mode string, configurations []evaluationExperimentConfiguration) string {
	resultsRootPath := evaluationExperimentOutputRoot(mode)
	if evaluationConfigurationsAreSparseControls(configurations) {
		return evaluationControlsResultsRootPath(resultsRootPath)
	}

	return resultsRootPath
}

func evaluationConfigurationsAreSparseControls(configurations []evaluationExperimentConfiguration) bool {
	if len(configurations) == 0 {
		return false
	}

	for _, configuration := range configurations {
		if !heuristicIsSparseControl(configuration.Heuristic) {
			return false
		}
	}

	return true
}

func evaluationExperimentUsesThreeOpt(mode string) bool {
	return mode == runModeEvaluation3Opt
}

func logTimestamp(timestamp time.Time) string {
	return timestamp.Format("2006-01-02 15:04:05")
}

func selectedInstanceSetForMode(mode, requestedInstanceSet string, instanceSetExplicit bool) string {
	if shouldRunEvaluationExperiments(mode) && !instanceSetExplicit {
		return instanceSetEvaluation
	}
	if shouldRunCacheRebuild(mode) && !instanceSetExplicit {
		return instanceSetAllKnown
	}

	return requestedInstanceSet
}

func resolveWorkerCount(requestedWorkers, cpuCount int) (int, error) {
	if requestedWorkers < 0 {
		return 0, fmt.Errorf("workers must be greater than or equal to zero")
	}
	if requestedWorkers > 0 {
		return requestedWorkers, nil
	}
	if cpuCount < 2 {
		return 1, nil
	}

	return max(1, cpuCount/2), nil
}

func configureWorkerCount(requestedWorkers int) (int, error) {
	workers, err := resolveWorkerCount(requestedWorkers, runtime.NumCPU())
	if err != nil {
		return 0, err
	}

	runtime.GOMAXPROCS(workers)
	return workers, nil
}

func Run(args []string) {
	flags := flag.NewFlagSet("atsp_aco_msa", flag.ExitOnError)
	instances := flags.String("instances", instanceSetTuning, "ATSP instance set to run: smoke, tuning, evaluation, or all-known")
	mode := flags.String("mode", runModeExperiment, "Run mode: experiment, analyze, all, evaluation, evaluation+3opt, or rebuild-cache")
	analysisScope := flags.String("analysis", analysisScopeAll, "Analysis scope for analyze mode: all, tuning, or gks-deviation")
	heuristic := flags.String("heuristic", "", "ACO heuristic modifier to use in experiment mode: omit or use all; otherwise strict-msa, rooted-msa, cycle-cover, or cycle-cover-msa-patching")
	evaluationHeuristic := flags.String("evaluation-heuristic", evaluationHeuristicAll, "Evaluation-mode heuristic to run: all, controls, baseline, strict-msa, rooted-msa, random-sparse, distance-ranked-sparse, shuffled-msa, cycle-cover, or cycle-cover-msa-patching")
	workers := flags.Int("workers", 0, "Maximum concurrent instance workers. 0 uses half of available logical CPUs; 1 preserves serial execution")
	flags.Parse(args)

	instanceSetExplicit := false
	heuristicExplicit := false
	flags.Visit(func(flag *flag.Flag) {
		if flag.Name == "instances" {
			instanceSetExplicit = true
		}
		if flag.Name == "heuristic" {
			heuristicExplicit = true
		}
	})

	selectedExperimentHeuristics, err := selectExperimentHeuristics(*heuristic, heuristicExplicit)
	if err != nil {
		fmt.Println(err)
		return
	}
	selectedEvaluationHeuristic := *evaluationHeuristic

	if !isValidRunMode(*mode) {
		fmt.Printf("Unsupported -mode value %q; use %q, %q, %q, %q, %q, or %q\n", *mode, runModeExperiment, runModeAnalyze, runModeAll, runModeEvaluation, runModeEvaluation3Opt, runModeRebuildCache)
		return
	}

	if !isValidAnalysisScope(*analysisScope) {
		fmt.Printf("Unsupported -analysis value %q; use %q, %q, or %q\n", *analysisScope, analysisScopeAll, analysisScopeTuning, analysisScopeGksDeviation)
		return
	}

	effectiveWorkers, err := configureWorkerCount(*workers)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Printf("Using %d worker(s); GOMAXPROCS=%d\n", effectiveWorkers, runtime.GOMAXPROCS(0))

	selectedInstances := selectedInstanceSetForMode(*mode, *instances, instanceSetExplicit)
	var evaluationConfigurations []evaluationExperimentConfiguration
	if shouldRunEvaluationExperiments(*mode) {
		configurations, err := selectEvaluationExperimentConfigurations(selectedEvaluationHeuristic)
		if err != nil {
			fmt.Println(err)
			return
		}
		evaluationConfigurations = configurations
	}

	atspsData, err := project.LoadSelectedAtspData(selectedInstances)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Printf("Selected %d ATSP instance(s) with -instances=%s\n", len(atspsData), selectedInstances)

	if shouldRunCacheRebuild(*mode) {
		err = runRebuildCacheMode(atspsData, effectiveWorkers)
		if err != nil {
			fmt.Println(err)
			return
		}
		return
	}

	if shouldRunExperiments(*mode) {
		err = runExperimentMode(atspsData, selectedExperimentHeuristics, effectiveWorkers)
		if err != nil {
			fmt.Println(err)
			return
		}
	}

	if shouldRunEvaluationExperiments(*mode) {
		err = runEvaluationExperimentMode(atspsData, evaluationExperimentOutputRootForConfigurations(*mode, evaluationConfigurations), evaluationExperimentUsesThreeOpt(*mode), evaluationConfigurations, effectiveWorkers, evaluationNumberOfExperiments)
		if err != nil {
			fmt.Println(err)
			return
		}

		if shouldRunAnalysisAfterEvaluationExperiments(*mode, selectedEvaluationHeuristic) {
			if err := runAnalysisMode(atspsData, *analysisScope, selectedExperimentHeuristics, effectiveWorkers); err != nil {
				fmt.Println(err)
				return
			}
		}
	}

	if shouldRunAnalysis(*mode) {
		if err := runAnalysisMode(atspsData, *analysisScope, selectedExperimentHeuristics, effectiveWorkers); err != nil {
			fmt.Println(err)
			return
		}
	}
}
