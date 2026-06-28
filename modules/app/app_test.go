package app

import (
	"atsp_aco_msa/modules/artifacts/cyclecover"
	"atsp_aco_msa/modules/experiments"
	"atsp_aco_msa/modules/project"
	"atsp_aco_msa/modules/tsplib"
	"os"
	"path/filepath"
	"reflect"
	"strconv"
	"strings"
	"testing"
)

func TestReadStatisticsRequiresFifteenColumns(t *testing.T) {
	dir := t.TempDir()
	validPath := filepath.Join(dir, "valid.csv")
	validRecord := []string{
		"1.00",
		"2.00",
		"0.80",
		"0.50",
		"5000",
		"1",
		"2.00",
		"3",
		"0",
		"0.00",
		"0",
		"1.00",
		"2.00",
		"3.00",
		"40.00",
	}

	validContent := strings.Join(experiments.StatisticsCsvHeader, ",") + "\n" + strings.Join(validRecord, ",") + "\n"
	if err := os.WriteFile(validPath, []byte(validContent), 0644); err != nil {
		t.Fatalf("failed to write valid statistics CSV: %v", err)
	}

	statistics, err := experiments.ReadStatistics(validPath)
	if err != nil {
		t.Fatalf("experiments.ReadStatistics returned unexpected error: %v", err)
	}

	if len(statistics) != 1 {
		t.Fatalf("expected one statistic, got %d", len(statistics))
	}

	if statistics[0].HeuristicWeight != 0.50 || statistics[0].MsaPatchBias != 0.0 || statistics[0].SuccessRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f msaPatchBias=%f successRate=%f", statistics[0].HeuristicWeight, statistics[0].MsaPatchBias, statistics[0].SuccessRate)
	}

	sixteenColumnPath := filepath.Join(dir, "sixteen-column.csv")
	sixteenColumnHeader := append(append([]string{}, experiments.StatisticsCsvHeader...), "Avg computation time [ms]")
	sixteenColumnRecord := append(append([]string{}, validRecord...), "0")
	sixteenColumnContent := strings.Join(sixteenColumnHeader, ",") + "\n" + strings.Join(sixteenColumnRecord, ",") + "\n"
	if err := os.WriteFile(sixteenColumnPath, []byte(sixteenColumnContent), 0644); err != nil {
		t.Fatalf("failed to write 16-column statistics CSV: %v", err)
	}

	if _, err := experiments.ReadStatistics(sixteenColumnPath); err == nil {
		t.Fatal("expected 16-column statistics CSV to be rejected")
	}
}

func TestReadStatisticsAcceptsPatchingMsaPatchBiasColumn(t *testing.T) {
	path := filepath.Join(t.TempDir(), "result_cycle_cover_msa_patching.csv")
	record := []string{
		"1.00",
		"2.00",
		"0.80",
		"0.50",
		"0.25",
		"5000",
		"1",
		"2.00",
		"3",
		"0",
		"0.00",
		"0",
		"1.00",
		"2.00",
		"3.00",
		"40.00",
	}
	content := strings.Join(experiments.StatisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching), ",") + "\n" + strings.Join(record, ",") + "\n"
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write patching statistics CSV: %v", err)
	}

	statistics, err := experiments.ReadStatistics(path)
	if err != nil {
		t.Fatalf("experiments.ReadStatistics returned unexpected error: %v", err)
	}
	if len(statistics) != 1 {
		t.Fatalf("expected one statistic, got %d", len(statistics))
	}
	if statistics[0].HeuristicWeight != 0.50 || statistics[0].MsaPatchBias != 0.25 || statistics[0].SuccessRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f msaPatchBias=%f successRate=%f", statistics[0].HeuristicWeight, statistics[0].MsaPatchBias, statistics[0].SuccessRate)
	}
}

func TestReadStatisticsAcceptsRandomSparseSeedColumn(t *testing.T) {
	path := filepath.Join(t.TempDir(), "result_random_sparse.csv")
	record := []string{
		"1.00",
		"2.00",
		"0.80",
		"0.90",
		"2",
		"5000",
		"1",
		"2.00",
		"3",
		"0",
		"0.00",
		"0",
		"1.00",
		"2.00",
		"3.00",
		"40.00",
	}
	content := strings.Join(experiments.StatisticsCsvHeaderForHeuristic(heuristicRandomSparse), ",") + "\n" + strings.Join(record, ",") + "\n"
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write random-sparse statistics CSV: %v", err)
	}

	statistics, err := experiments.ReadStatistics(path)
	if err != nil {
		t.Fatalf("experiments.ReadStatistics returned unexpected error: %v", err)
	}
	if len(statistics) != 1 {
		t.Fatalf("expected one statistic, got %d", len(statistics))
	}
	if statistics[0].HeuristicWeight != 0.90 || statistics[0].RandomSeed != 2 || statistics[0].SuccessRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f randomSeed=%d successRate=%f", statistics[0].HeuristicWeight, statistics[0].RandomSeed, statistics[0].SuccessRate)
	}
}

func TestReadRootedMsaHeuristicsRequiresOneMsaPerVertex(t *testing.T) {
	root := t.TempDir()
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}, 3, filepath.Join(root, "results"))
	atspData.MsaHeuristicDirectoryPath = filepath.Join(root, "msa", "sample")

	writeTestRootMsaMatrix(t, atspData.MsaHeuristicDirectoryPath, 0, [][]float64{
		{0, 1, 1},
		{0, 0, 0},
		{0, 0, 0},
	})
	writeTestRootMsaMatrix(t, atspData.MsaHeuristicDirectoryPath, 1, [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, 0},
	})

	if _, err := readRootedMsaHeuristics(atspData); err == nil {
		t.Fatal("expected incomplete rooted MSA cache to be rejected")
	}

	writeTestRootMsaMatrix(t, atspData.MsaHeuristicDirectoryPath, 2, [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, 0},
	})

	rootedMsas, err := readRootedMsaHeuristics(atspData)
	if err != nil {
		t.Fatalf("expected complete rooted MSA cache to be accepted: %v", err)
	}
	if len(rootedMsas) != 3 {
		t.Fatalf("expected three rooted MSAs, got %d", len(rootedMsas))
	}
}

func TestSaveHeuristicStatisticsWritesSingleComparisonCsv(t *testing.T) {
	path := filepath.Join(t.TempDir(), "result.csv")
	rows := []HeuristicExperimentStatistics{
		{
			Heuristic: heuristicBaseline,
			Statistics: ExperimentsDataStatistics{
				ExperimentParameters: ExperimentParameters{
					Alpha:           1.0,
					Beta:            2.0,
					Rho:             0.8,
					HeuristicWeight: 0.0,
					Iterations:      500,
				},
				AverageBestDeviation: 3.0,
			},
		},
		{
			Heuristic: heuristicCycleCover,
			Statistics: ExperimentsDataStatistics{
				ExperimentParameters: ExperimentParameters{
					Alpha:           1.0,
					Beta:            2.0,
					Rho:             0.8,
					HeuristicWeight: 0.8,
					Iterations:      500,
				},
				AverageBestDeviation: 1.0,
			},
		},
	}

	if err := experiments.SaveHeuristicStatistics(path, rows); err != nil {
		t.Fatalf("experiments.SaveHeuristicStatistics returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read heuristic statistics CSV: %v", err)
	}

	lines := strings.Split(strings.TrimSpace(string(contentBytes)), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected header plus two rows, got %d lines: %v", len(lines), lines)
	}
	if !strings.HasPrefix(lines[0], "Heuristic,Alpha,Beta,Rho,Heuristic weight") {
		t.Fatalf("unexpected header: %s", lines[0])
	}
	if !strings.HasPrefix(lines[1], heuristicCycleCover+",") {
		t.Fatalf("expected best average-deviation row first, got: %s", lines[1])
	}
	if !strings.HasPrefix(lines[2], heuristicBaseline+",") {
		t.Fatalf("expected baseline row second, got: %s", lines[2])
	}
}

func TestRunEvaluationResultsAnalysisReadsExistingEvaluationResults(t *testing.T) {
	resultsRoot := t.TempDir()
	oldEvaluationResultsDirectoryName := project.EvaluationResultsDirectoryName
	project.EvaluationResultsDirectoryName = filepath.Join(resultsRoot, "evaluation")
	defer func() {
		project.EvaluationResultsDirectoryName = oldEvaluationResultsDirectoryName
	}()

	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	evaluationAtspData := project.WithExperimentOutputRoot(atspData, project.EvaluationResultsDirectoryName)
	rows := []HeuristicExperimentStatistics{
		{
			Heuristic: heuristicBaseline,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 4.25,
				SuccessRate:          10.0,
			},
		},
		{
			Heuristic: heuristicStrictMsa,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 2.50,
				SuccessRate:          20.0,
			},
		},
		{
			Heuristic: heuristicCycleCover,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 1.75,
				SuccessRate:          30.0,
			},
		},
	}
	if err := experiments.SaveHeuristicStatistics(evaluationAtspData.ResultFilePath, rows); err != nil {
		t.Fatalf("experiments.SaveHeuristicStatistics returned unexpected error: %v", err)
	}

	summaryPath, _, saved, err := runEvaluationResultsAnalysis([]AtspData{atspData}, nil, project.EvaluationResultsDirectoryName)
	if err != nil {
		t.Fatalf("runEvaluationResultsAnalysis returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected evaluation results summary to be saved")
	}
	if summaryPath != filepath.Join(project.EvaluationResultsDirectoryName, "summary.md") {
		t.Fatalf("unexpected summary path: %s", summaryPath)
	}
	if _, err := os.Stat(summaryPath); err != nil {
		t.Fatalf("expected evaluation results summary file to exist: %v", err)
	}
}

func TestRunEvaluationResultsAnalysisUsesProvidedResultsRoot(t *testing.T) {
	resultsRoot := t.TempDir()
	evaluationThreeOptRoot := filepath.Join(resultsRoot, "evaluation_3opt")
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	evaluationThreeOptAtspData := project.WithExperimentOutputRoot(atspData, evaluationThreeOptRoot)
	rows := []HeuristicExperimentStatistics{
		{
			Heuristic: heuristicBaseline,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 3.00,
				SuccessRate:          10.0,
			},
		},
		{
			Heuristic: heuristicStrictMsa,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 2.00,
				SuccessRate:          20.0,
			},
		},
		{
			Heuristic: heuristicCycleCover,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 1.00,
				SuccessRate:          30.0,
			},
		},
	}
	if err := experiments.SaveHeuristicStatistics(evaluationThreeOptAtspData.ResultFilePath, rows); err != nil {
		t.Fatalf("experiments.SaveHeuristicStatistics returned unexpected error: %v", err)
	}

	summaryPath, _, saved, err := runEvaluationResultsAnalysis([]AtspData{atspData}, nil, evaluationThreeOptRoot)
	if err != nil {
		t.Fatalf("runEvaluationResultsAnalysis returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected evaluation+3opt results summary to be saved")
	}
	if summaryPath != filepath.Join(evaluationThreeOptRoot, "summary.md") {
		t.Fatalf("unexpected summary path: %s", summaryPath)
	}

	for _, name := range []string{"summary.md", "baseline_comparison.md", "rbg_outlier_summary.md", "convergence_summary.md"} {
		path := filepath.Join(evaluationThreeOptRoot, name)
		if _, err := os.Stat(path); err != nil {
			t.Fatalf("expected %s to exist: %v", path, err)
		}
	}
}

func TestBuildHeuristicModifiersReturnsNeutralMatrixForBaseline(t *testing.T) {
	msaHeuristic := [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}

	modifiers := buildHeuristicModifiers(heuristicBaseline, nil, msaHeuristic, nil, newDefaultExperimentParameters(1.0))
	expected := [][]float64{
		{0, 0, 0},
		{0, 0, 0},
		{0, 0, 0},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected baseline modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildHeuristicModifiersBaselineCanUseDistanceMatrixDimension(t *testing.T) {
	matrix := [][]float64{
		{0, 1},
		{2, 0},
	}

	modifiers := buildHeuristicModifiers(heuristicBaseline, matrix, nil, nil, newDefaultExperimentParameters(1.0))
	expected := [][]float64{
		{0, 0},
		{0, 0},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected baseline modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestCycleCoverBuild(t *testing.T) {
	matrix := [][]float64{
		{0, 5, 1},
		{1, 0, 5},
		{5, 1, 0},
	}

	cycleCoverMatrix, cost, err := cyclecover.Build(matrix)
	if err != nil {
		t.Fatalf("cyclecover.Build returned unexpected error: %v", err)
	}

	expected := [][]float64{
		{0, 0, 1},
		{1, 0, 0},
		{0, 1, 0},
	}
	if !reflect.DeepEqual(cycleCoverMatrix, expected) {
		t.Fatalf("unexpected cycle-cover matrix\nwant: %v\n got: %v", expected, cycleCoverMatrix)
	}
	if cost != 3 {
		t.Fatalf("expected cycle-cover cost 3, got %f", cost)
	}
}

func TestCycleCoverReadOrCreateCreatesCache(t *testing.T) {
	matrix := [][]float64{
		{0, 5, 1},
		{1, 0, 5},
		{5, 1, 0},
	}
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", matrix, 3, t.TempDir())
	atspData.CycleCoverDirectoryPath = filepath.Join(t.TempDir(), "cycle_cover", "sample")

	cycleCoverMatrix, cost, err := cyclecover.ReadOrCreate(atspData.Matrix, atspData.CycleCoverDirectoryPath)
	if err != nil {
		t.Fatalf("cyclecover.ReadOrCreate returned unexpected error: %v", err)
	}

	expected := [][]float64{
		{0, 0, 1},
		{1, 0, 0},
		{0, 1, 0},
	}
	if !reflect.DeepEqual(cycleCoverMatrix, expected) {
		t.Fatalf("unexpected cycle-cover matrix\nwant: %v\n got: %v", expected, cycleCoverMatrix)
	}
	if cost != 3 {
		t.Fatalf("expected cycle-cover cost 3, got %f", cost)
	}
	assertPathExists(t, cyclecover.Path(atspData.CycleCoverDirectoryPath))
}

func TestCycleCoverReadOrCreateReusesValidCache(t *testing.T) {
	matrix := [][]float64{
		{0, 5, 1},
		{1, 0, 5},
		{5, 1, 0},
	}
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", matrix, 3, t.TempDir())
	atspData.CycleCoverDirectoryPath = filepath.Join(t.TempDir(), "cycle_cover", "sample")
	cached := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}
	if err := cyclecover.Save(atspData.CycleCoverDirectoryPath, cached); err != nil {
		t.Fatalf("cyclecover.Save returned unexpected error: %v", err)
	}

	cycleCoverMatrix, cost, err := cyclecover.ReadOrCreate(atspData.Matrix, atspData.CycleCoverDirectoryPath)
	if err != nil {
		t.Fatalf("cyclecover.ReadOrCreate returned unexpected error: %v", err)
	}

	if !reflect.DeepEqual(cycleCoverMatrix, cached) {
		t.Fatalf("expected cached cycle-cover matrix\nwant: %v\n got: %v", cached, cycleCoverMatrix)
	}
	if cost != 15 {
		t.Fatalf("expected cached cycle-cover cost 15, got %f", cost)
	}
}

func TestCycleCoverReadOrCreateRegeneratesInvalidCache(t *testing.T) {
	matrix := [][]float64{
		{0, 5, 1},
		{1, 0, 5},
		{5, 1, 0},
	}
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", matrix, 3, t.TempDir())
	atspData.CycleCoverDirectoryPath = filepath.Join(t.TempDir(), "cycle_cover", "sample")
	if err := cyclecover.Save(atspData.CycleCoverDirectoryPath, [][]float64{
		{0, 0, 0},
		{0, 0, 0},
		{0, 0, 0},
	}); err != nil {
		t.Fatalf("cyclecover.Save returned unexpected error: %v", err)
	}

	cycleCoverMatrix, cost, err := cyclecover.ReadOrCreate(atspData.Matrix, atspData.CycleCoverDirectoryPath)
	if err != nil {
		t.Fatalf("cyclecover.ReadOrCreate returned unexpected error: %v", err)
	}

	expected := [][]float64{
		{0, 0, 1},
		{1, 0, 0},
		{0, 1, 0},
	}
	if !reflect.DeepEqual(cycleCoverMatrix, expected) {
		t.Fatalf("expected regenerated cycle-cover matrix\nwant: %v\n got: %v", expected, cycleCoverMatrix)
	}
	if cost != 3 {
		t.Fatalf("expected regenerated cycle-cover cost 3, got %f", cost)
	}
}

func TestHeuristicSpecificPathsKeepMsaHeuristicBaselinePaths(t *testing.T) {
	atspData := project.MakeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)

	if resultFilePathForHeuristic(atspData, heuristicBaseline) != filepath.Join(project.ResultsDirectoryName, "test", "result_baseline.csv") {
		t.Fatalf("unexpected baseline result path: %s", resultFilePathForHeuristic(atspData, heuristicBaseline))
	}

	if resultFilePathForHeuristic(atspData, heuristicStrictMsa) != filepath.Join(project.ResultsDirectoryName, "test", project.ResultFileName) {
		t.Fatalf("MSA heuristic result path should keep the existing baseline location")
	}

	if resultFilePathForHeuristic(atspData, heuristicRootedMsa) != filepath.Join(project.ResultsDirectoryName, "test", "result_rooted_msa.csv") {
		t.Fatalf("unexpected rooted-MSA result path: %s", resultFilePathForHeuristic(atspData, heuristicRootedMsa))
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCover) != filepath.Join(project.ResultsDirectoryName, "test", "result_cycle_cover.csv") {
		t.Fatalf("unexpected cycle-cover result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCover))
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCoverPatching) != filepath.Join(project.ResultsDirectoryName, "test", "result_cycle_cover_patching.csv") {
		t.Fatalf("unexpected cycle-cover patching result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCoverPatching))
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching) != filepath.Join(project.ResultsDirectoryName, "test", "result_cycle_cover_msa_patching.csv") {
		t.Fatalf("unexpected cycle-cover MSA-patching result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching))
	}
}

func TestWithExperimentOutputRootMovesOutputsButKeepsMsaHeuristicCache(t *testing.T) {
	atspData := project.MakeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)
	output := project.WithExperimentOutputRoot(atspData, project.EvaluationResultsDirectoryName)

	if output.MsaHeuristicDirectoryPath != atspData.MsaHeuristicDirectoryPath {
		t.Fatalf("expected MSA heuristic cache path to stay %s, got %s", atspData.MsaHeuristicDirectoryPath, output.MsaHeuristicDirectoryPath)
	}
	if output.CycleCoverDirectoryPath != atspData.CycleCoverDirectoryPath {
		t.Fatalf("expected cycle-cover cache path to stay %s, got %s", atspData.CycleCoverDirectoryPath, output.CycleCoverDirectoryPath)
	}
	if output.CycleCoverHeatmapPlotPath != atspData.CycleCoverHeatmapPlotPath {
		t.Fatalf("expected cycle-cover heatmap path to stay %s, got %s", atspData.CycleCoverHeatmapPlotPath, output.CycleCoverHeatmapPlotPath)
	}

	expectedResultPath := filepath.Join(project.EvaluationResultsDirectoryName, "test", project.ResultFileName)
	if output.ResultFilePath != expectedResultPath {
		t.Fatalf("unexpected evaluation result path\nwant: %s\n got: %s", expectedResultPath, output.ResultFilePath)
	}

	if output.OptimalUniqueToursCsvPath != atspData.OptimalUniqueToursCsvPath {
		t.Fatalf("expected solutions path to stay %s, got %s", atspData.OptimalUniqueToursCsvPath, output.OptimalUniqueToursCsvPath)
	}
}

func TestEvaluationExperimentConfigurationsUseFixedBalancedComparison(t *testing.T) {
	configs := evaluationExperimentConfigurations()
	expected := []struct {
		heuristic string
		weight    float64
		bias      float64
	}{
		{heuristicBaseline, defaultBaselineHeuristicWeight, 0.0},
		{heuristicStrictMsa, evaluationStrictMsaHeuristicWeight, 0.0},
		{heuristicRootedMsa, evaluationRootedMsaHeuristicWeight, 0.0},
		{heuristicCycleCover, evaluationCycleCoverWeight, 0.0},
		{heuristicCycleCoverPatching, evaluationCycleCoverPatchingWeight, 0.0},
		{heuristicCycleCoverMsaPatching, evaluationCycleCoverMsaPatchingWeight, evaluationCycleCoverMsaPatchingMsaPatchBias},
	}

	if len(configs) != len(expected) {
		t.Fatalf("expected %d evaluation experiment configurations, got %d", len(expected), len(configs))
	}

	for i, config := range configs {
		if config.Heuristic != expected[i].heuristic {
			t.Fatalf("unexpected evaluation heuristic at %d: want %s got %s", i, expected[i].heuristic, config.Heuristic)
		}
		if len(config.Parameters) != 1 {
			t.Fatalf("expected one evaluation parameter set for %s, got %d", config.Heuristic, len(config.Parameters))
		}

		parameters := config.Parameters[0]
		if parameters.Alpha != defaultExperimentAlpha ||
			parameters.Beta != defaultExperimentBeta ||
			parameters.Rho != defaultExperimentRho ||
			parameters.HeuristicWeight != expected[i].weight ||
			parameters.MsaPatchBias != expected[i].bias {
			t.Fatalf("unexpected evaluation parameters for %s: %+v", config.Heuristic, parameters)
		}
	}
}

func TestGenerateParametersUsesMsaPatchBiasOnlyForPatching(t *testing.T) {
	baselineParameters := generateParameters(heuristicBaseline)
	if len(baselineParameters) != 1 {
		t.Fatalf("expected one baseline parameter set, got %d", len(baselineParameters))
	}
	if baselineParameters[0].HeuristicWeight != defaultBaselineHeuristicWeight {
		t.Fatalf("unexpected baseline heuristic weight: %+v", baselineParameters[0])
	}

	msaParameters := generateParameters(heuristicStrictMsa)
	if len(msaParameters) != 10 {
		t.Fatalf("expected 10 MSA heuristic parameter sets, got %d", len(msaParameters))
	}
	for _, parameters := range msaParameters {
		if parameters.HeuristicWeight == 0 {
			t.Fatalf("expected MSA heuristic tuning to skip baseline-equivalent zero weight, got %+v", parameters)
		}
		if parameters.MsaPatchBias != 0 {
			t.Fatalf("expected MSA patch bias to be zero for non-patching heuristic, got %+v", parameters)
		}
	}

	patchingParameters := generateParameters(heuristicCycleCoverMsaPatching)
	if len(patchingParameters) != 50 {
		t.Fatalf("expected 50 patching parameter sets, got %d", len(patchingParameters))
	}

	zeroWeightCount := 0
	zeroBiasCount := 0
	for _, parameters := range patchingParameters {
		if parameters.HeuristicWeight == 0 {
			zeroWeightCount++
		}
		if parameters.MsaPatchBias == 0 {
			zeroBiasCount++
		}
	}
	if zeroWeightCount != 0 {
		t.Fatalf("expected patching tuning to skip baseline-equivalent zero weight, got %d", zeroWeightCount)
	}
	if zeroBiasCount != 10 {
		t.Fatalf("expected one pure patching zero-bias row for each nonzero heuristic weight, got %d", zeroBiasCount)
	}

}

func TestSetDimensionDependentParametersScalesIterationsLinearly(t *testing.T) {
	tests := []struct {
		dimension          int
		expectedIterations int
	}{
		{dimension: 17, expectedIterations: 100},
		{dimension: 49, expectedIterations: 100},
		{dimension: 50, expectedIterations: 1500},
		{dimension: 66, expectedIterations: 1980},
		{dimension: 100, expectedIterations: 3000},
		{dimension: 171, expectedIterations: 5130},
		{dimension: 323, expectedIterations: 9690},
		{dimension: 443, expectedIterations: 13290},
	}

	for _, test := range tests {
		t.Run(strconv.Itoa(test.dimension), func(t *testing.T) {
			parameters := ExperimentParameters{Iterations: 1}
			setDimensionDependentParameters(test.dimension, &parameters)

			if parameters.Iterations != test.expectedIterations {
				t.Fatalf("expected %d iterations, got %d", test.expectedIterations, parameters.Iterations)
			}
		})
	}
}

func TestSelectExperimentHeuristics(t *testing.T) {
	selected, err := selectExperimentHeuristics("", false)
	if err != nil {
		t.Fatalf("selectExperimentHeuristics omitted returned error: %v", err)
	}
	if !reflect.DeepEqual(selected, experimentHeuristics) {
		t.Fatalf("omitted heuristic should select all experiment heuristics\nwant: %v\n got: %v", experimentHeuristics, selected)
	}
	for _, heuristic := range selected {
		if heuristic == heuristicBaseline {
			t.Fatal("default experiment heuristics should not include baseline")
		}
	}

	selected, err = selectExperimentHeuristics(experimentHeuristicAll, true)
	if err != nil {
		t.Fatalf("selectExperimentHeuristics(all) returned error: %v", err)
	}
	if !reflect.DeepEqual(selected, experimentHeuristics) {
		t.Fatalf("all heuristic should select all experiment heuristics\nwant: %v\n got: %v", experimentHeuristics, selected)
	}

	selected, err = selectExperimentHeuristics(heuristicCycleCover, true)
	if err != nil {
		t.Fatalf("selectExperimentHeuristics(cycle-cover) returned error: %v", err)
	}
	if !reflect.DeepEqual(selected, []string{heuristicCycleCover}) {
		t.Fatalf("explicit heuristic should select only that heuristic, got %v", selected)
	}

	selected, err = selectExperimentHeuristics(heuristicRootedMsa, true)
	if err != nil {
		t.Fatalf("selectExperimentHeuristics(rooted-msa) returned error: %v", err)
	}
	if !reflect.DeepEqual(selected, []string{heuristicRootedMsa}) {
		t.Fatalf("explicit rooted MSA heuristic should select only rooted MSA, got %v", selected)
	}

	if _, err := selectExperimentHeuristics(heuristicBaseline, true); err == nil {
		t.Fatal("expected baseline experiment heuristic to be rejected")
	}

	if _, err := selectExperimentHeuristics("unknown", true); err == nil {
		t.Fatal("expected unknown experiment heuristic to be rejected")
	}
}

func TestSelectEvaluationExperimentConfigurations(t *testing.T) {
	allConfigurations, err := selectEvaluationExperimentConfigurations(evaluationHeuristicAll)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(all) returned error: %v", err)
	}
	if len(allConfigurations) != len(evaluationExperimentConfigurations()) {
		t.Fatalf("expected all evaluation configurations, got %d", len(allConfigurations))
	}
	if evaluationConfigurationsAreSparseControls(allConfigurations) {
		t.Fatal("main evaluation all configuration should not be treated as sparse controls")
	}

	controlConfigurations, err := selectEvaluationExperimentConfigurations(evaluationHeuristicControls)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(controls) returned error: %v", err)
	}
	if len(controlConfigurations) != len(evaluationControlExperimentConfigurations()) {
		t.Fatalf("expected all evaluation control configurations, got %d", len(controlConfigurations))
	}
	if !evaluationConfigurationsAreSparseControls(controlConfigurations) {
		t.Fatal("control configurations should be treated as sparse controls")
	}

	cycleCoverConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if len(cycleCoverConfigurations) != 1 || cycleCoverConfigurations[0].Heuristic != heuristicCycleCover {
		t.Fatalf("expected only cycle-cover configuration, got %+v", cycleCoverConfigurations)
	}

	cycleCoverPatchingConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCoverPatching)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover-patching) returned error: %v", err)
	}
	if len(cycleCoverPatchingConfigurations) != 1 || cycleCoverPatchingConfigurations[0].Heuristic != heuristicCycleCoverPatching {
		t.Fatalf("expected only cycle-cover patching configuration, got %+v", cycleCoverPatchingConfigurations)
	}

	patchingConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if len(patchingConfigurations) != 1 || patchingConfigurations[0].Heuristic != heuristicCycleCoverMsaPatching {
		t.Fatalf("expected only cycle-cover MSA-patching configuration, got %+v", patchingConfigurations)
	}

	randomSparseConfigurations, err := selectEvaluationExperimentConfigurations(heuristicRandomSparse)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(random-sparse) returned error: %v", err)
	}
	if len(randomSparseConfigurations) != 1 || randomSparseConfigurations[0].Heuristic != heuristicRandomSparse || len(randomSparseConfigurations[0].Parameters) != len(randomSparseSeeds) {
		t.Fatalf("expected only random-sparse control configuration, got %+v", randomSparseConfigurations)
	}

	shuffledMsaConfigurations, err := selectEvaluationExperimentConfigurations(heuristicShuffledMsa)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(shuffled-msa) returned error: %v", err)
	}
	if len(shuffledMsaConfigurations) != 1 || shuffledMsaConfigurations[0].Heuristic != heuristicShuffledMsa || len(shuffledMsaConfigurations[0].Parameters) != len(shuffledMsaSeeds) {
		t.Fatalf("expected only shuffled-MSA control configuration, got %+v", shuffledMsaConfigurations)
	}

	if _, err := selectEvaluationExperimentConfigurations("unknown"); err == nil {
		t.Fatal("expected invalid evaluation heuristic to be rejected")
	}
}

func TestEvaluationExperimentOutputRootUsesControlsSubdirectoryForSparseControls(t *testing.T) {
	mainConfigurations, err := selectEvaluationExperimentConfigurations(evaluationHeuristicAll)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(all) returned error: %v", err)
	}
	if evaluationExperimentOutputRootForConfigurations(runModeEvaluation, mainConfigurations) != project.EvaluationResultsDirectoryName {
		t.Fatalf("main evaluation run should use %s", project.EvaluationResultsDirectoryName)
	}

	controlConfigurations, err := selectEvaluationExperimentConfigurations(evaluationHeuristicControls)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(controls) returned error: %v", err)
	}
	expectedControlsRoot := filepath.Join(filepath.Dir(project.EvaluationResultsDirectoryName), "controls")
	if evaluationExperimentOutputRootForConfigurations(runModeEvaluation, controlConfigurations) != expectedControlsRoot {
		t.Fatalf("evaluation controls should use %s", expectedControlsRoot)
	}
}

func TestEvaluationConfigurationsUseMsaHeuristicOnlyWhenNeeded(t *testing.T) {
	baselineConfigurations, err := selectEvaluationExperimentConfigurations(heuristicBaseline)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(baseline) returned error: %v", err)
	}
	if evaluationConfigurationsUseMsaHeuristic(baselineConfigurations) {
		t.Fatal("baseline-only evaluation run should not require MSA heuristic")
	}

	cycleCoverConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if evaluationConfigurationsUseMsaHeuristic(cycleCoverConfigurations) {
		t.Fatal("cycle-cover-only evaluation run should not require MSA heuristic")
	}

	cycleCoverPatchingConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCoverPatching)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover-patching) returned error: %v", err)
	}
	if evaluationConfigurationsUseMsaHeuristic(cycleCoverPatchingConfigurations) {
		t.Fatal("cycle-cover patching evaluation run should not require MSA heuristic")
	}

	strictMsaConfigurations, err := selectEvaluationExperimentConfigurations(heuristicStrictMsa)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(strict-msa) returned error: %v", err)
	}
	if !evaluationConfigurationsUseMsaHeuristic(strictMsaConfigurations) {
		t.Fatal("strict MSA evaluation run should require MSA heuristic")
	}

	rootedMsaConfigurations, err := selectEvaluationExperimentConfigurations(heuristicRootedMsa)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(rooted-msa) returned error: %v", err)
	}
	if !evaluationConfigurationsUseMsaHeuristic(rootedMsaConfigurations) {
		t.Fatal("rooted MSA evaluation run should require MSA heuristic")
	}

	patchingConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if !evaluationConfigurationsUseMsaHeuristic(patchingConfigurations) {
		t.Fatal("cycle-cover MSA-patching evaluation run should require MSA heuristic")
	}
}

func TestCycleCoverCacheIsNeededOnlyForCycleCoverHeuristics(t *testing.T) {
	if heuristicsUseCycleCover([]string{heuristicStrictMsa, heuristicRootedMsa}) {
		t.Fatal("MSA-only tuning heuristics should not require cycle-cover cache")
	}
	if !heuristicsUseCycleCover([]string{heuristicStrictMsa, heuristicCycleCover}) {
		t.Fatal("cycle-cover tuning heuristic should require cycle-cover cache")
	}
	if !heuristicsUseCycleCover([]string{heuristicCycleCoverMsaPatching}) {
		t.Fatal("cycle-cover MSA-patching tuning heuristic should require cycle-cover cache")
	}

	baselineConfigurations, err := selectEvaluationExperimentConfigurations(heuristicBaseline)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(baseline) returned error: %v", err)
	}
	if evaluationConfigurationsUseCycleCover(baselineConfigurations) {
		t.Fatal("baseline evaluation run should not require cycle-cover cache")
	}

	strictMsaConfigurations, err := selectEvaluationExperimentConfigurations(heuristicStrictMsa)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(strict-msa) returned error: %v", err)
	}
	if evaluationConfigurationsUseCycleCover(strictMsaConfigurations) {
		t.Fatal("strict MSA evaluation run should not require cycle-cover cache")
	}

	cycleCoverConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if !evaluationConfigurationsUseCycleCover(cycleCoverConfigurations) {
		t.Fatal("cycle-cover evaluation run should require cycle-cover cache")
	}

	cycleCoverPatchingConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCoverPatching)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover-patching) returned error: %v", err)
	}
	if !evaluationConfigurationsUseCycleCover(cycleCoverPatchingConfigurations) {
		t.Fatal("cycle-cover patching evaluation run should require cycle-cover cache")
	}

	patchingConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if !evaluationConfigurationsUseCycleCover(patchingConfigurations) {
		t.Fatal("cycle-cover MSA-patching evaluation run should require cycle-cover cache")
	}
}

func TestSaveEvaluationHeuristicStatisticsMergesSelectedHeuristicIntoExistingResultCsv(t *testing.T) {
	path := filepath.Join(t.TempDir(), "result.csv")
	existing := []HeuristicExperimentStatistics{
		{
			Heuristic:  heuristicBaseline,
			Statistics: makeTestExperimentStatistics(0.0, 5.0, 10.0),
		},
		{
			Heuristic:  heuristicStrictMsa,
			Statistics: makeTestExperimentStatistics(evaluationStrictMsaHeuristicWeight, 3.0, 20.0),
		},
		{
			Heuristic:  heuristicCycleCover,
			Statistics: makeTestExperimentStatistics(evaluationCycleCoverWeight, 7.0, 0.0),
		},
		{
			Heuristic:  "obsolete",
			Statistics: makeTestExperimentStatistics(1.0, 0.5, 100.0),
		},
	}
	if err := experiments.SaveHeuristicStatistics(path, existing); err != nil {
		t.Fatalf("failed to seed evaluation result CSV: %v", err)
	}

	replacement := []HeuristicExperimentStatistics{
		{
			Heuristic:  heuristicCycleCover,
			Statistics: makeTestExperimentStatistics(evaluationCycleCoverWeight, 1.5, 40.0),
		},
	}
	cycleCoverConfigurations, err := selectEvaluationExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectEvaluationExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if err := saveEvaluationHeuristicStatistics(path, replacement, cycleCoverConfigurations); err != nil {
		t.Fatalf("saveEvaluationHeuristicStatistics returned error: %v", err)
	}

	statistics, err := experiments.ReadHeuristicStatistics(path)
	if err != nil {
		t.Fatalf("failed to read merged evaluation result CSV: %v", err)
	}
	if len(statistics) != 3 {
		t.Fatalf("expected three heuristic rows after merge, got %d", len(statistics))
	}

	byHeuristic := make(map[string]HeuristicExperimentStatistics, len(statistics))
	for _, statistic := range statistics {
		byHeuristic[statistic.Heuristic] = statistic
	}

	if byHeuristic[heuristicBaseline].Statistics.AverageBestDeviation != 5.0 {
		t.Fatalf("baseline row was not preserved: %+v", byHeuristic[heuristicBaseline])
	}
	if byHeuristic[heuristicStrictMsa].Statistics.AverageBestDeviation != 3.0 {
		t.Fatalf("strict MSA row was not preserved: %+v", byHeuristic[heuristicStrictMsa])
	}
	if byHeuristic[heuristicCycleCover].Statistics.AverageBestDeviation != 1.5 ||
		byHeuristic[heuristicCycleCover].Statistics.SuccessRate != 40.0 {
		t.Fatalf("cycle-cover row was not replaced: %+v", byHeuristic[heuristicCycleCover])
	}
	if _, ok := byHeuristic["obsolete"]; ok {
		t.Fatal("obsolete heuristic row should not be preserved in evaluation result CSV")
	}
}

func TestEvaluationModesAndExperimentHeuristicsAreValid(t *testing.T) {
	if !isValidRunMode(runModeEvaluation) {
		t.Fatal("evaluation run mode should be valid")
	}
	if !isValidRunMode(runModeEvaluation3Opt) {
		t.Fatal("evaluation+3opt run mode should be valid")
	}
	if !isValidRunMode(runModeRebuildCache) {
		t.Fatal("rebuild-cache run mode should be valid")
	}
	if isValidHeuristic(heuristicBaseline) {
		t.Fatal("baseline heuristic should not be valid in normal experiment mode")
	}
	if isValidHeuristic(heuristicRandomSparse) {
		t.Fatal("random sparse control should not be valid in normal experiment mode")
	}
	if isValidHeuristic(heuristicDistanceRankedSparse) {
		t.Fatal("distance-ranked sparse control should not be valid in normal experiment mode")
	}
	if isValidHeuristic(heuristicShuffledMsa) {
		t.Fatal("shuffled MSA control should not be valid in normal experiment mode")
	}
	if isValidHeuristic(heuristicCycleCoverPatching) {
		t.Fatal("cycle-cover patching should be evaluation-only")
	}
}

func TestEvaluation3OptModeUsesSeparateOutputRootAndThreeOpt(t *testing.T) {
	if evaluationExperimentOutputRoot(runModeEvaluation3Opt) != project.EvaluationThreeOptResultsDirectoryName {
		t.Fatalf("evaluation+3opt should use %s", project.EvaluationThreeOptResultsDirectoryName)
	}
	if evaluationExperimentOutputRoot(runModeEvaluation) != project.EvaluationResultsDirectoryName {
		t.Fatalf("evaluation should use %s", project.EvaluationResultsDirectoryName)
	}
	if !evaluationExperimentUsesThreeOpt(runModeEvaluation3Opt) {
		t.Fatal("evaluation+3opt should enable reduced 3-opt")
	}
	if evaluationExperimentUsesThreeOpt(runModeEvaluation) {
		t.Fatal("evaluation should not enable reduced 3-opt")
	}
}

func TestFullEvaluationRunsTriggerAnalysis(t *testing.T) {
	if !shouldRunAnalysisAfterEvaluationExperiments(runModeEvaluation, evaluationHeuristicAll) {
		t.Fatal("full evaluation run should trigger analysis")
	}
	if !shouldRunAnalysisAfterEvaluationExperiments(runModeEvaluation3Opt, evaluationHeuristicAll) {
		t.Fatal("full evaluation+3opt run should trigger analysis")
	}
	if shouldRunAnalysisAfterEvaluationExperiments(runModeEvaluation, heuristicStrictMsa) {
		t.Fatal("single evaluation heuristic should not trigger analysis")
	}
	if shouldRunAnalysisAfterEvaluationExperiments(runModeEvaluation, evaluationHeuristicControls) {
		t.Fatal("evaluation controls should not trigger analysis")
	}
	if shouldRunAnalysisAfterEvaluationExperiments(runModeExperiment, evaluationHeuristicAll) {
		t.Fatal("experiment mode should not trigger evaluation analysis")
	}
}

func TestAnalysisScopeTuningIsValid(t *testing.T) {
	if !isValidAnalysisScope(analysisScopeTuning) {
		t.Fatal("tuning analysis scope should be valid")
	}
}

func TestRunAnalysisModeTuningRegeneratesTuningSummary(t *testing.T) {
	originalResultsDirectoryName := project.ResultsDirectoryName
	project.ResultsDirectoryName = t.TempDir()
	t.Cleanup(func() {
		project.ResultsDirectoryName = originalResultsDirectoryName
	})

	tuningAtspData := project.MakeAtspDataInResultsDirectory(project.TuningInstanceFiles[0], [][]float64{{0, 1}, {1, 0}}, 2, project.ResultsDirectoryName)
	experiments.SaveStatistics(resultFilePathForHeuristic(tuningAtspData, heuristicCycleCoverMsaPatching), heuristicCycleCoverMsaPatching, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.6, 2.0, 20.0),
	})

	nonTuningAtspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, project.ResultsDirectoryName)
	experiments.SaveStatistics(resultFilePathForHeuristic(nonTuningAtspData, heuristicCycleCoverMsaPatching), heuristicCycleCoverMsaPatching, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.9, 9.0, 90.0),
	})

	if err := runAnalysisMode([]AtspData{nonTuningAtspData}, analysisScopeTuning, []string{heuristicCycleCoverMsaPatching}, 1); err != nil {
		t.Fatalf("runAnalysisMode(tuning) returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(filepath.Join(project.ResultsDirectoryName, "tuning_summary.md"))
	if err != nil {
		t.Fatalf("expected tuning summary to be regenerated: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "# Tuning Summary")
	assertContains(t, content, "Cycle-cover MSA patching")
	assertContains(t, content, "| 0.60 | 0.00 | 2.00 | 2.00 | 1 | 1 | 20.00 |")
	assertDoesNotContain(t, content, "| 0.90 | 0.00 | 9.00 | 9.00 | 1 | 1 | 90.00 |")
}

func TestResolveWorkerCount(t *testing.T) {
	tests := []struct {
		name             string
		requestedWorkers int
		cpuCount         int
		expectedWorkers  int
		expectError      bool
	}{
		{name: "default uses half CPUs", requestedWorkers: 0, cpuCount: 8, expectedWorkers: 4},
		{name: "default keeps at least one worker", requestedWorkers: 0, cpuCount: 1, expectedWorkers: 1},
		{name: "explicit one is serial", requestedWorkers: 1, cpuCount: 8, expectedWorkers: 1},
		{name: "explicit positive value is used", requestedWorkers: 3, cpuCount: 8, expectedWorkers: 3},
		{name: "negative rejected", requestedWorkers: -1, cpuCount: 8, expectError: true},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			workers, err := resolveWorkerCount(test.requestedWorkers, test.cpuCount)
			if test.expectError {
				if err == nil {
					t.Fatal("expected error")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if workers != test.expectedWorkers {
				t.Fatalf("expected %d workers, got %d", test.expectedWorkers, workers)
			}
		})
	}
}

func TestRunEvaluationExperimentParametersRejectsInvalidWorkerCount(t *testing.T) {
	matrix := [][]float64{{0, 1}, {1, 0}}
	_, err := runEvaluationExperimentParameters(
		"sample",
		heuristicBaseline,
		[]ExperimentParameters{newDefaultExperimentParameters(defaultBaselineHeuristicWeight)},
		1,
		2,
		matrix,
		nil,
		nil,
		nil,
		false,
		len(matrix),
		0)
	if err == nil || !strings.Contains(err.Error(), "workers must be at least one") {
		t.Fatalf("expected invalid parameter worker error, got %v", err)
	}
}

func TestRunRebuildCacheModeRebuildsCacheAndPreservesAnalysisPlots(t *testing.T) {
	root := t.TempDir()
	matrix := [][]float64{
		{0, 1, 5, 2},
		{4, 0, 1, 3},
		{2, 4, 0, 1},
		{1, 2, 4, 0},
	}
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", matrix, 0, filepath.Join(root, "results"))
	atspData.MsaHeuristicDirectoryPath = filepath.Join(root, "msa", "sample")
	atspData.MsaHeuristicHeatmapPlotPath = filepath.Join(atspData.MsaHeuristicDirectoryPath, "plots", "msa_heuristic_heatmap.png")
	atspData.MsaHeuristicHistogramPlotPath = filepath.Join(atspData.MsaHeuristicDirectoryPath, "plots", "msa_heuristic_histogram.png")
	atspData.MsaHeuristicToursOverlapHeatmapPlotPath = filepath.Join(root, "solutions", "sample", "plots", "msa_heuristic_tours_overlap_heatmap.png")
	atspData.CycleCoverDirectoryPath = filepath.Join(root, "cycle_cover", "sample")
	atspData.CycleCoverHeatmapPlotPath = filepath.Join(atspData.CycleCoverDirectoryPath, "plots", "cycle_cover_heatmap.png")

	staleRootMsaPath := filepath.Join(atspData.MsaHeuristicDirectoryPath, "msas", "99.csv")
	if err := os.MkdirAll(filepath.Dir(staleRootMsaPath), 0700); err != nil {
		t.Fatalf("failed to create stale root MSA directory: %v", err)
	}
	if err := os.WriteFile(staleRootMsaPath, []byte("0,1\n0,0\n"), 0644); err != nil {
		t.Fatalf("failed to write stale root MSA file: %v", err)
	}
	if err := os.MkdirAll(filepath.Dir(atspData.MsaHeuristicToursOverlapHeatmapPlotPath), 0700); err != nil {
		t.Fatalf("failed to create MSA plots directory: %v", err)
	}
	if err := os.WriteFile(atspData.MsaHeuristicToursOverlapHeatmapPlotPath, []byte("analysis plot"), 0644); err != nil {
		t.Fatalf("failed to write existing analysis plot: %v", err)
	}
	staleCycleCoverPath := filepath.Join(atspData.CycleCoverDirectoryPath, "stale.txt")
	if err := os.MkdirAll(atspData.CycleCoverDirectoryPath, 0700); err != nil {
		t.Fatalf("failed to create stale cycle-cover directory: %v", err)
	}
	if err := os.WriteFile(staleCycleCoverPath, []byte("stale"), 0644); err != nil {
		t.Fatalf("failed to write stale cycle-cover file: %v", err)
	}

	if err := runRebuildCacheMode([]AtspData{atspData}, 1); err != nil {
		t.Fatalf("runRebuildCacheMode returned unexpected error: %v", err)
	}

	if _, err := os.Stat(staleRootMsaPath); !os.IsNotExist(err) {
		t.Fatalf("expected stale root MSA file to be removed, stat error: %v", err)
	}
	if _, err := os.Stat(staleCycleCoverPath); !os.IsNotExist(err) {
		t.Fatalf("expected stale cycle-cover file to be removed, stat error: %v", err)
	}
	assertPathExists(t, filepath.Join(atspData.MsaHeuristicDirectoryPath, "msa_heuristic.csv"))
	assertPathExists(t, atspData.MsaHeuristicHeatmapPlotPath)
	assertPathExists(t, atspData.MsaHeuristicHistogramPlotPath)
	assertPathExists(t, atspData.MsaHeuristicToursOverlapHeatmapPlotPath)
	assertPathExists(t, cyclecover.Path(atspData.CycleCoverDirectoryPath))
	assertPathExists(t, atspData.CycleCoverHeatmapPlotPath)

	rootMsaFiles, err := filepath.Glob(filepath.Join(atspData.MsaHeuristicDirectoryPath, "msas", "*.csv"))
	if err != nil {
		t.Fatalf("failed to glob root MSA files: %v", err)
	}
	if len(rootMsaFiles) != len(matrix) {
		t.Fatalf("expected %d root MSA files, got %d", len(matrix), len(rootMsaFiles))
	}
}

func TestSelectAtspFilesTuning(t *testing.T) {
	paths, err := testTsplibFiles()
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := project.SelectAtspFiles(paths, instanceSetTuning)
	if err != nil {
		t.Fatalf("project.SelectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, project.TuningInstanceFiles) {
		t.Fatalf("tuning selection mismatch\nwant: %v\n got: %v", project.TuningInstanceFiles, selectedFiles)
	}
}

func TestSelectAtspFilesSmoke(t *testing.T) {
	paths, err := testTsplibFiles()
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := project.SelectAtspFiles(paths, instanceSetSmoke)
	if err != nil {
		t.Fatalf("project.SelectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, project.SmokeInstanceFiles) {
		t.Fatalf("smoke selection mismatch\nwant: %v\n got: %v", project.SmokeInstanceFiles, selectedFiles)
	}
}

func TestSelectAtspFilesEvaluation(t *testing.T) {
	paths, err := testTsplibFiles()
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := project.SelectAtspFiles(paths, instanceSetEvaluation)
	if err != nil {
		t.Fatalf("project.SelectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, project.EvaluationInstanceFiles) {
		t.Fatalf("evaluation selection mismatch\nwant: %v\n got: %v", project.EvaluationInstanceFiles, selectedFiles)
	}
}

func TestSelectedAtspFilesHaveKnownOptima(t *testing.T) {
	paths, err := testTsplibFiles()
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	for _, instanceSet := range []string{instanceSetSmoke, instanceSetTuning, instanceSetEvaluation, instanceSetAllKnown} {
		selected, err := project.SelectAtspFiles(paths, instanceSet)
		if err != nil {
			t.Fatalf("project.SelectAtspFiles(%s) returned unexpected error: %v", instanceSet, err)
		}

		for _, selectedPath := range selected {
			name, _, knownOptimal, err := tsplib.ParseFile(selectedPath)
			if err != nil {
				t.Fatalf("failed to parse %s: %v", selectedPath, err)
			}

			if knownOptimal == 0 {
				t.Fatalf("%s selected %s without a known optimum", instanceSet, name)
			}
		}
	}
}

func TestTuningAndEvaluationPartitionKnownInstances(t *testing.T) {
	paths, err := testTsplibFiles()
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	tuning, err := project.SelectAtspFiles(paths, instanceSetTuning)
	if err != nil {
		t.Fatalf("project.SelectAtspFiles(tuning) returned unexpected error: %v", err)
	}
	evaluation, err := project.SelectAtspFiles(paths, instanceSetEvaluation)
	if err != nil {
		t.Fatalf("project.SelectAtspFiles(evaluation) returned unexpected error: %v", err)
	}
	allKnown, err := project.SelectAtspFiles(paths, instanceSetAllKnown)
	if err != nil {
		t.Fatalf("project.SelectAtspFiles(all-known) returned unexpected error: %v", err)
	}
	smoke, err := project.SelectAtspFiles(paths, instanceSetSmoke)
	if err != nil {
		t.Fatalf("project.SelectAtspFiles(smoke) returned unexpected error: %v", err)
	}

	partition := make(map[string]string, len(tuning)+len(evaluation))
	for _, path := range tuning {
		partition[filepath.Base(path)] = instanceSetTuning
	}
	for _, path := range evaluation {
		fileName := filepath.Base(path)
		if previousSet, ok := partition[fileName]; ok {
			t.Fatalf("%s appears in both %s and %s", fileName, previousSet, instanceSetEvaluation)
		}
		partition[fileName] = instanceSetEvaluation
	}

	smokeFiles := make(map[string]struct{}, len(smoke))
	for _, path := range smoke {
		smokeFiles[filepath.Base(path)] = struct{}{}
	}

	expectedKnownCount := len(allKnown) - len(smokeFiles)
	if len(partition) != expectedKnownCount {
		t.Fatalf("tuning/evaluation partition has %d files, all-known excluding smoke has %d", len(partition), expectedKnownCount)
	}
	for _, path := range allKnown {
		fileName := filepath.Base(path)
		if _, ok := smokeFiles[fileName]; ok {
			continue
		}
		if _, ok := partition[fileName]; !ok {
			t.Fatalf("%s is known but missing from tuning/evaluation partition", fileName)
		}
	}
}

func testTsplibFiles() ([]string, error) {
	return filepath.Glob(filepath.Join("..", "..", "tsplib_files", "*.atsp"))
}

func TestSelectedInstanceSetForMode(t *testing.T) {
	if selected := selectedInstanceSetForMode(runModeEvaluation, instanceSetTuning, false); selected != instanceSetEvaluation {
		t.Fatalf("evaluation mode without explicit instances should default to %s, got %s", instanceSetEvaluation, selected)
	}
	if selected := selectedInstanceSetForMode(runModeEvaluation3Opt, instanceSetTuning, false); selected != instanceSetEvaluation {
		t.Fatalf("evaluation+3opt mode without explicit instances should default to %s, got %s", instanceSetEvaluation, selected)
	}
	if selected := selectedInstanceSetForMode(runModeRebuildCache, instanceSetTuning, false); selected != instanceSetAllKnown {
		t.Fatalf("rebuild-cache mode without explicit instances should default to %s, got %s", instanceSetAllKnown, selected)
	}
	if selected := selectedInstanceSetForMode(runModeEvaluation, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("evaluation mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeRebuildCache, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("rebuild-cache mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeExperiment, instanceSetTuning, false); selected != instanceSetTuning {
		t.Fatalf("experiment mode should keep requested instances, got %s", selected)
	}
}

func makeTestExperimentStatistics(heuristicWeight, averageBestDeviation, successRate float64) ExperimentsDataStatistics {
	return ExperimentsDataStatistics{
		ExperimentParameters: ExperimentParameters{
			Alpha:           defaultExperimentAlpha,
			Beta:            defaultExperimentBeta,
			Rho:             defaultExperimentRho,
			HeuristicWeight: heuristicWeight,
			Iterations:      500,
		},
		MinBestAtIteration:               1,
		AverageBestAtIteration:           2.0,
		MaxBestAtIteration:               3,
		MinThreeOptImprovementsCount:     0,
		AverageThreeOptImprovementsCount: 0.0,
		MaxThreeOptImprovementsCount:     0,
		MinBestDeviation:                 averageBestDeviation,
		AverageBestDeviation:             averageBestDeviation,
		MaxBestDeviation:                 averageBestDeviation,
		SuccessRate:                      successRate,
	}
}

func makeTestRandomSparseExperimentStatistics(randomSeed int64, averageBestDeviation, successRate float64) ExperimentsDataStatistics {
	statistics := makeTestExperimentStatistics(evaluationStrictMsaHeuristicWeight, averageBestDeviation, successRate)
	statistics.RandomSeed = randomSeed
	return statistics
}

func writeTestControlStatisticsForWeightSummary(t *testing.T, atspData AtspData, heuristic string) {
	t.Helper()

	firstInstance := strings.Contains(atspData.Name, "sample-a")
	statistics := []ExperimentsDataStatistics{
		makeTestSeededExperimentStatistics(0.4, 1, 4.0, 0.0),
		makeTestSeededExperimentStatistics(0.4, 2, 4.0, 0.0),
		makeTestSeededExperimentStatistics(0.8, 1, 3.0, 0.0),
		makeTestSeededExperimentStatistics(0.8, 2, 3.0, 0.0),
	}
	if !firstInstance {
		statistics = []ExperimentsDataStatistics{
			makeTestSeededExperimentStatistics(0.4, 1, 2.0, 0.0),
			makeTestSeededExperimentStatistics(0.4, 2, 2.0, 0.0),
			makeTestSeededExperimentStatistics(0.8, 1, 4.0, 0.0),
			makeTestSeededExperimentStatistics(0.8, 2, 4.0, 0.0),
		}
	}

	experiments.SaveStatistics(resultFilePathForHeuristic(atspData, heuristic), heuristic, statistics)
}

func makeTestSeededExperimentStatistics(heuristicWeight float64, randomSeed int64, averageBestDeviation, successRate float64) ExperimentsDataStatistics {
	statistics := makeTestExperimentStatistics(heuristicWeight, averageBestDeviation, successRate)
	statistics.RandomSeed = randomSeed
	return statistics
}

func writeTestOptimalToursCsv(t *testing.T, path string, tours []string) {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		t.Fatalf("failed to create optimal tours directory: %v", err)
	}

	lines := []string{"Tour"}
	for _, tour := range tours {
		lines = append(lines, strconv.Quote(tour))
	}
	content := strings.Join(lines, "\n") + "\n"
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write optimal tours CSV: %v", err)
	}
}

func writeTestMsaHeuristicMatrix(t *testing.T, rootPath string, matrix [][]float64) {
	t.Helper()
	if err := os.MkdirAll(rootPath, 0700); err != nil {
		t.Fatalf("failed to create MSA heuristic directory: %v", err)
	}

	rows := make([]string, 0, len(matrix))
	for _, row := range matrix {
		values := make([]string, len(row))
		for i, value := range row {
			values[i] = strconv.FormatFloat(value, 'f', -1, 64)
		}
		rows = append(rows, strings.Join(values, ","))
	}

	path := filepath.Join(rootPath, "msa_heuristic.csv")
	if err := os.WriteFile(path, []byte(strings.Join(rows, "\n")+"\n"), 0644); err != nil {
		t.Fatalf("failed to write test MSA heuristic matrix: %v", err)
	}
}

func writeTestRootMsaMatrix(t *testing.T, rootPath string, root int, matrix [][]float64) {
	t.Helper()
	msasPath := filepath.Join(rootPath, "msas")
	if err := os.MkdirAll(msasPath, 0700); err != nil {
		t.Fatalf("failed to create test root MSA directory: %v", err)
	}

	rows := make([]string, 0, len(matrix))
	for _, row := range matrix {
		values := make([]string, len(row))
		for i, value := range row {
			values[i] = strconv.FormatFloat(value, 'f', -1, 64)
		}
		rows = append(rows, strings.Join(values, ","))
	}

	path := filepath.Join(msasPath, strconv.Itoa(root)+".csv")
	if err := os.WriteFile(path, []byte(strings.Join(rows, "\n")+"\n"), 0644); err != nil {
		t.Fatalf("failed to write test root MSA matrix: %v", err)
	}
}

func assertContains(t *testing.T, content, expected string) {
	t.Helper()
	if !strings.Contains(content, expected) {
		t.Fatalf("expected content to contain:\n%s\n\ncontent:\n%s", expected, content)
	}
}

func assertDoesNotContain(t *testing.T, content, unexpected string) {
	t.Helper()
	if strings.Contains(content, unexpected) {
		t.Fatalf("expected content not to contain:\n%s\n\ncontent:\n%s", unexpected, content)
	}
}

func assertPathExists(t *testing.T, path string) {
	t.Helper()
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("expected %s to exist: %v", path, err)
	}
}
