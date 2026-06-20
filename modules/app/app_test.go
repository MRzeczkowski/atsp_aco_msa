package app

import (
	"atsp_aco_msa/modules/analysis/structure"
	"atsp_aco_msa/modules/artifacts/cyclecover"
	"atsp_aco_msa/modules/models"
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

	validContent := strings.Join(statisticsCsvHeader, ",") + "\n" + strings.Join(validRecord, ",") + "\n"
	if err := os.WriteFile(validPath, []byte(validContent), 0644); err != nil {
		t.Fatalf("failed to write valid statistics CSV: %v", err)
	}

	statistics, err := readStatistics(validPath)
	if err != nil {
		t.Fatalf("readStatistics returned unexpected error: %v", err)
	}

	if len(statistics) != 1 {
		t.Fatalf("expected one statistic, got %d", len(statistics))
	}

	if statistics[0].HeuristicWeight != 0.50 || statistics[0].MsaPatchBias != 0.0 || statistics[0].SuccessRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f msaPatchBias=%f successRate=%f", statistics[0].HeuristicWeight, statistics[0].MsaPatchBias, statistics[0].SuccessRate)
	}

	legacyPath := filepath.Join(dir, "legacy.csv")
	legacyHeader := append(append([]string{}, statisticsCsvHeader...), "Avg computation time [ms]")
	legacyRecord := append(append([]string{}, validRecord...), "0")
	legacyContent := strings.Join(legacyHeader, ",") + "\n" + strings.Join(legacyRecord, ",") + "\n"
	if err := os.WriteFile(legacyPath, []byte(legacyContent), 0644); err != nil {
		t.Fatalf("failed to write legacy statistics CSV: %v", err)
	}

	if _, err := readStatistics(legacyPath); err == nil {
		t.Fatal("expected legacy 16-column statistics CSV to be rejected")
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
	content := strings.Join(statisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching), ",") + "\n" + strings.Join(record, ",") + "\n"
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write patching statistics CSV: %v", err)
	}

	statistics, err := readStatistics(path)
	if err != nil {
		t.Fatalf("readStatistics returned unexpected error: %v", err)
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
	content := strings.Join(statisticsCsvHeaderForHeuristic(heuristicRandomSparse), ",") + "\n" + strings.Join(record, ",") + "\n"
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write random-sparse statistics CSV: %v", err)
	}

	statistics, err := readStatistics(path)
	if err != nil {
		t.Fatalf("readStatistics returned unexpected error: %v", err)
	}
	if len(statistics) != 1 {
		t.Fatalf("expected one statistic, got %d", len(statistics))
	}
	if statistics[0].HeuristicWeight != 0.90 || statistics[0].RandomSeed != 2 || statistics[0].SuccessRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f randomSeed=%d successRate=%f", statistics[0].HeuristicWeight, statistics[0].RandomSeed, statistics[0].SuccessRate)
	}
}

func TestSaveRandomSparseControlReportComparesMsaAgainstRandomSparse(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	finalRoot := filepath.Join(root, "final")
	controlsRoot := finalControlsResultsRootPath(finalRoot)
	first := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	finalFirst := project.WithExperimentOutputRoot(first, finalRoot)
	finalSecond := project.WithExperimentOutputRoot(second, finalRoot)
	controlFirst := project.WithExperimentOutputRoot(first, controlsRoot)
	controlSecond := project.WithExperimentOutputRoot(second, controlsRoot)

	if err := saveHeuristicStatistics(finalFirst.ResultFilePath, []HeuristicExperimentStatistics{
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 2.0, 10.0)},
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(0.8, 0.5, 100.0)},
	}); err != nil {
		t.Fatalf("failed to write first final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlFirst, heuristicRandomSparse), heuristicRandomSparse, []ExperimentsDataStatistics{
		makeTestRandomSparseExperimentStatistics(1, 4.0, 0.0),
		makeTestRandomSparseExperimentStatistics(2, 5.0, 0.0),
		makeTestRandomSparseExperimentStatistics(3, 6.0, 0.0),
	})

	if err := saveHeuristicStatistics(finalSecond.ResultFilePath, []HeuristicExperimentStatistics{
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 4.0, 0.0)},
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(0.8, 0.5, 100.0)},
	}); err != nil {
		t.Fatalf("failed to write second final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlSecond, heuristicRandomSparse), heuristicRandomSparse, []ExperimentsDataStatistics{
		makeTestRandomSparseExperimentStatistics(1, 2.0, 20.0),
		makeTestRandomSparseExperimentStatistics(2, 3.0, 20.0),
		makeTestRandomSparseExperimentStatistics(3, 4.0, 20.0),
	})

	reportPath := filepath.Join(controlsRoot, "random_sparse_control.md")
	saved, err := saveRandomSparseControlReport(reportPath, []AtspData{second, first}, finalRoot, controlsRoot)
	if err != nil {
		t.Fatalf("saveRandomSparseControlReport returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected random sparse control report to be saved")
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read random sparse control report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "Strict MSA had lower average best deviation than the random-sparse mean in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: Strict MSA 3.00%, random sparse 4.00%, delta -1.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func TestSaveDistanceRankedSparseControlReportComparesMsaAgainstDistanceRankedSparse(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	finalRoot := filepath.Join(root, "final")
	controlsRoot := finalControlsResultsRootPath(finalRoot)
	first := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	finalFirst := project.WithExperimentOutputRoot(first, finalRoot)
	finalSecond := project.WithExperimentOutputRoot(second, finalRoot)
	controlFirst := project.WithExperimentOutputRoot(first, controlsRoot)
	controlSecond := project.WithExperimentOutputRoot(second, controlsRoot)

	if err := saveHeuristicStatistics(finalFirst.ResultFilePath, []HeuristicExperimentStatistics{
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 2.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write first final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlFirst, heuristicDistanceRankedSparse), heuristicDistanceRankedSparse, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 4.0, 0.0),
	})

	if err := saveHeuristicStatistics(finalSecond.ResultFilePath, []HeuristicExperimentStatistics{
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 4.0, 0.0)},
	}); err != nil {
		t.Fatalf("failed to write second final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlSecond, heuristicDistanceRankedSparse), heuristicDistanceRankedSparse, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 2.0, 20.0),
	})

	reportPath := filepath.Join(controlsRoot, "distance_ranked_sparse_control.md")
	saved, err := saveDistanceRankedSparseControlReport(reportPath, []AtspData{second, first}, finalRoot, controlsRoot)
	if err != nil {
		t.Fatalf("saveDistanceRankedSparseControlReport returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected distance-ranked sparse control report to be saved")
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read distance-ranked sparse control report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "Strict MSA had lower average best deviation than the distance-ranked sparse control in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: Strict MSA 3.00%, distance-ranked sparse 3.00%, delta +0.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func TestSaveShuffledMsaControlReportComparesMsaAgainstShuffledMsa(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	finalRoot := filepath.Join(root, "final")
	controlsRoot := finalControlsResultsRootPath(finalRoot)
	first := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	finalFirst := project.WithExperimentOutputRoot(first, finalRoot)
	finalSecond := project.WithExperimentOutputRoot(second, finalRoot)
	controlFirst := project.WithExperimentOutputRoot(first, controlsRoot)
	controlSecond := project.WithExperimentOutputRoot(second, controlsRoot)

	if err := saveHeuristicStatistics(finalFirst.ResultFilePath, []HeuristicExperimentStatistics{
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 2.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write first final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlFirst, heuristicShuffledMsa), heuristicShuffledMsa, []ExperimentsDataStatistics{
		makeTestRandomSparseExperimentStatistics(1, 4.0, 0.0),
		makeTestRandomSparseExperimentStatistics(2, 5.0, 0.0),
		makeTestRandomSparseExperimentStatistics(3, 6.0, 0.0),
	})

	if err := saveHeuristicStatistics(finalSecond.ResultFilePath, []HeuristicExperimentStatistics{
		{Heuristic: heuristicStrictMsa, Statistics: makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 4.0, 0.0)},
	}); err != nil {
		t.Fatalf("failed to write second final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlSecond, heuristicShuffledMsa), heuristicShuffledMsa, []ExperimentsDataStatistics{
		makeTestRandomSparseExperimentStatistics(1, 2.0, 20.0),
		makeTestRandomSparseExperimentStatistics(2, 3.0, 20.0),
		makeTestRandomSparseExperimentStatistics(3, 4.0, 20.0),
	})

	reportPath := filepath.Join(controlsRoot, "shuffled_msa_control.md")
	saved, err := saveShuffledMsaControlReport(reportPath, []AtspData{second, first}, finalRoot, controlsRoot)
	if err != nil {
		t.Fatalf("saveShuffledMsaControlReport returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected shuffled MSA control report to be saved")
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read shuffled MSA control report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "Strict MSA had lower average best deviation than the shuffled MSA mean in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: Strict MSA 3.00%, shuffled MSA 4.00%, delta -1.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func TestReadMsaHeuristicMatrixForResultRootUsesCompositeMsaForMatrixFallback(t *testing.T) {
	root := t.TempDir()
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}, 3, filepath.Join(root, "results"))
	atspData.MsaHeuristicDirectoryPath = filepath.Join(root, "msa", "sample")

	compositeMsa := [][]float64{
		{0, 2, 0},
		{0, 0, 2},
		{2, 0, 0},
	}
	writeTestMsaHeuristicMatrix(t, atspData.MsaHeuristicDirectoryPath, compositeMsa)
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
	writeTestRootMsaMatrix(t, atspData.MsaHeuristicDirectoryPath, 2, [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, 0},
	})

	normalMatrix, err := readMsaHeuristicMatrixForResultRoot(atspData, heuristicStrictMsa, project.FinalResultsDirectoryName)
	if err != nil {
		t.Fatalf("read normal MSA matrix: %v", err)
	}
	if !reflect.DeepEqual(normalMatrix, compositeMsa) {
		t.Fatalf("expected normal run to use composite MSA, got %v", normalMatrix)
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

	if err := saveHeuristicStatistics(path, rows); err != nil {
		t.Fatalf("saveHeuristicStatistics returned unexpected error: %v", err)
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

func TestSaveFinalResultsSummaryWritesMarkdownTableWithHighlightedFindings(t *testing.T) {
	resultsRoot := t.TempDir()
	firstAtspData := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	secondAtspData := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
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
			Heuristic: heuristicRootedMsa,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 3.00,
				SuccessRate:          15.0,
			},
		},
		{
			Heuristic: heuristicCycleCover,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 1.75,
				SuccessRate:          30.0,
			},
		},
		{
			Heuristic: heuristicCycleCoverMsaPatching,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 2.00,
				SuccessRate:          25.0,
			},
		},
	}
	secondRows := []HeuristicExperimentStatistics{
		{
			Heuristic: heuristicBaseline,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 6.75,
				SuccessRate:          20.0,
			},
		},
		{
			Heuristic: heuristicStrictMsa,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 3.50,
				SuccessRate:          40.0,
			},
		},
		{
			Heuristic: heuristicRootedMsa,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 4.00,
				SuccessRate:          35.0,
			},
		},
		{
			Heuristic: heuristicCycleCover,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 2.25,
				SuccessRate:          60.0,
			},
		},
		{
			Heuristic: heuristicCycleCoverMsaPatching,
			Statistics: ExperimentsDataStatistics{
				AverageBestDeviation: 2.50,
				SuccessRate:          55.0,
			},
		},
	}

	if err := saveHeuristicStatistics(firstAtspData.ResultFilePath, rows); err != nil {
		t.Fatalf("saveHeuristicStatistics returned unexpected error: %v", err)
	}
	if err := saveHeuristicStatistics(secondAtspData.ResultFilePath, secondRows); err != nil {
		t.Fatalf("saveHeuristicStatistics returned unexpected error: %v", err)
	}

	summaryPath := filepath.Join(resultsRoot, "summary.md")
	if err := saveFinalResultsSummary([]AtspData{firstAtspData, secondAtspData}, summaryPath); err != nil {
		t.Fatalf("saveFinalResultsSummary returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(summaryPath)
	if err != nil {
		t.Fatalf("failed to read final summary Markdown: %v", err)
	}

	lines := strings.Split(strings.TrimSpace(string(contentBytes)), "\n")
	expected := []string{
		"# Final Results Summary",
		"",
		"## Findings",
		"",
		"- **Cycle cover has the lowest average best deviation overall: 2.00%.**",
		"- **Cycle cover has the highest average success rate overall: 45.00%.**",
		"- **Best-or-tied average best deviation counts: Baseline 0/2, Strict MSA 0/2, Rooted MSA 0/2, Cycle cover 2/2, Cycle-cover MSA patching 0/2.**",
		"",
		"<table>",
		"<thead>",
		"<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">Baseline</th><th colspan=\"2\">Strict MSA</th><th colspan=\"2\">Rooted MSA</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>",
		"<tr><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th></tr>",
		"</thead>",
		"<tbody>",
		"<tr><td>sample-a</td><td align=\"right\">4.25</td><td align=\"right\">10.00</td><td align=\"right\">2.50</td><td align=\"right\">20.00</td><td align=\"right\">3.00</td><td align=\"right\">15.00</td><td align=\"right\"><strong>1.75</strong></td><td align=\"right\"><strong>30.00</strong></td><td align=\"right\">2.00</td><td align=\"right\">25.00</td></tr>",
		"<tr><td>sample-b</td><td align=\"right\">6.75</td><td align=\"right\">20.00</td><td align=\"right\">3.50</td><td align=\"right\">40.00</td><td align=\"right\">4.00</td><td align=\"right\">35.00</td><td align=\"right\"><strong>2.25</strong></td><td align=\"right\"><strong>60.00</strong></td><td align=\"right\">2.50</td><td align=\"right\">55.00</td></tr>",
		"<tr><td><strong>Average</strong></td><td align=\"right\">5.50</td><td align=\"right\">15.00</td><td align=\"right\">3.00</td><td align=\"right\">30.00</td><td align=\"right\">3.50</td><td align=\"right\">25.00</td><td align=\"right\"><strong>2.00</strong></td><td align=\"right\"><strong>45.00</strong></td><td align=\"right\">2.25</td><td align=\"right\">40.00</td></tr>",
		"</tbody>",
		"</table>",
	}
	if !reflect.DeepEqual(lines, expected) {
		t.Fatalf("unexpected final summary Markdown\nwant: %v\n got: %v", expected, lines)
	}
}

func TestRunFinalResultsAnalysisReadsExistingFinalResults(t *testing.T) {
	resultsRoot := t.TempDir()
	oldFinalResultsDirectoryName := project.FinalResultsDirectoryName
	project.FinalResultsDirectoryName = filepath.Join(resultsRoot, "final")
	defer func() {
		project.FinalResultsDirectoryName = oldFinalResultsDirectoryName
	}()

	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	finalAtspData := project.WithExperimentOutputRoot(atspData, project.FinalResultsDirectoryName)
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
	if err := saveHeuristicStatistics(finalAtspData.ResultFilePath, rows); err != nil {
		t.Fatalf("saveHeuristicStatistics returned unexpected error: %v", err)
	}

	summaryPath, _, saved, err := runFinalResultsAnalysis([]AtspData{atspData}, nil, project.FinalResultsDirectoryName)
	if err != nil {
		t.Fatalf("runFinalResultsAnalysis returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected final results summary to be saved")
	}
	if summaryPath != filepath.Join(project.FinalResultsDirectoryName, "summary.md") {
		t.Fatalf("unexpected summary path: %s", summaryPath)
	}
	if _, err := os.Stat(summaryPath); err != nil {
		t.Fatalf("expected final results summary file to exist: %v", err)
	}
}

func TestRunFinalResultsAnalysisUsesProvidedResultsRoot(t *testing.T) {
	resultsRoot := t.TempDir()
	finalThreeOptRoot := filepath.Join(resultsRoot, "final_3opt")
	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	finalThreeOptAtspData := project.WithExperimentOutputRoot(atspData, finalThreeOptRoot)
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
	if err := saveHeuristicStatistics(finalThreeOptAtspData.ResultFilePath, rows); err != nil {
		t.Fatalf("saveHeuristicStatistics returned unexpected error: %v", err)
	}

	summaryPath, _, saved, err := runFinalResultsAnalysis([]AtspData{atspData}, nil, finalThreeOptRoot)
	if err != nil {
		t.Fatalf("runFinalResultsAnalysis returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected final+3opt results summary to be saved")
	}
	if summaryPath != filepath.Join(finalThreeOptRoot, "summary.md") {
		t.Fatalf("unexpected summary path: %s", summaryPath)
	}

	for _, name := range []string{"summary.md", "pairwise_performance.md", "convergence_summary.md"} {
		path := filepath.Join(finalThreeOptRoot, name)
		if _, err := os.Stat(path); err != nil {
			t.Fatalf("expected %s to exist: %v", path, err)
		}
	}
}

func TestSaveFinalThreeOptComparisonReportShowsHiddenHeuristicEffect(t *testing.T) {
	path := filepath.Join(t.TempDir(), "comparison_to_final.md")
	finalRows := []finalResultsSummaryRow{
		{
			instance: "sample",
			metrics: map[string]finalResultsSummaryMetric{
				heuristicBaseline: {
					AverageMinDeviation:  10.0,
					SuccessRate:          10.0,
					AverageBestIteration: 50.0,
					Iterations:           100,
				},
				heuristicStrictMsa: {
					AverageMinDeviation:  8.0,
					SuccessRate:          12.0,
					AverageBestIteration: 60.0,
					Iterations:           100,
				},
				heuristicCycleCover: {
					AverageMinDeviation:  7.0,
					SuccessRate:          11.0,
					AverageBestIteration: 30.0,
					Iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  7.5,
					SuccessRate:          13.0,
					AverageBestIteration: 35.0,
					Iterations:           100,
				},
			},
		},
	}
	finalThreeOptRows := []finalResultsSummaryRow{
		{
			instance: "sample",
			metrics: map[string]finalResultsSummaryMetric{
				heuristicBaseline: {
					AverageMinDeviation:  1.0,
					SuccessRate:          50.0,
					AverageBestIteration: 20.0,
					Iterations:           100,
				},
				heuristicStrictMsa: {
					AverageMinDeviation:  0.9,
					SuccessRate:          51.0,
					AverageBestIteration: 19.0,
					Iterations:           100,
				},
				heuristicCycleCover: {
					AverageMinDeviation:  0.8,
					SuccessRate:          52.0,
					AverageBestIteration: 19.0,
					Iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  0.85,
					SuccessRate:          53.0,
					AverageBestIteration: 18.0,
					Iterations:           100,
				},
			},
		},
	}

	if err := saveFinalThreeOptComparisonReport(path, finalRows, finalThreeOptRows); err != nil {
		t.Fatalf("saveFinalThreeOptComparisonReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read comparison report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "# Reduced 3-Opt Impact")
	assertContains(t, content, "Cycle-cover MSA patching +2.50 -> +0.15 pp")
	assertContains(t, content, "<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>")
	assertContains(t, content, "<tr><td>Strict MSA</td><td align=\"right\">+2.00</td><td align=\"right\">+0.10</td><td align=\"right\">5.00</td></tr>")
	assertContains(t, content, "<tr><td>Cycle cover</td><td align=\"right\">+3.00</td><td align=\"right\">+0.20</td><td align=\"right\">6.67</td></tr>")
	assertContains(t, content, "<tr><td>Cycle-cover MSA patching</td><td align=\"right\">+2.50</td><td align=\"right\">+0.15</td><td align=\"right\">6.00</td></tr>")
}

func TestSaveStructuralSimilarityReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "structural_similarity.md")
	if err := saveStructuralSimilarityReport(path, sampleStructuralAnalyses()); err != nil {
		t.Fatalf("saveStructuralSimilarityReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read structural similarity report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "- **Precision vs found-optimal tours: MSA heuristic 71.43%, cycle cover 66.67%, cycle-cover MSA patching 63.64%.**")
	assertContains(t, content, "- **Recall vs found-optimal tours: MSA heuristic 35.71%, cycle cover 42.86%, cycle-cover MSA patching 50.00%.**")
	assertContains(t, content, "<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">MSA heuristic</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>")
	assertContains(t, content, "<tr><td>a</td><td align=\"right\">50.00</td><td align=\"right\">25.00</td><td align=\"right\"><strong>75.00</strong></td><td align=\"right\"><strong>75.00</strong></td><td align=\"right\">60.00</td><td align=\"right\"><strong>75.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Total</strong></td><td align=\"right\"><strong>71.43</strong></td><td align=\"right\">35.71</td><td align=\"right\">66.67</td><td align=\"right\">42.86</td><td align=\"right\">63.64</td><td align=\"right\"><strong>50.00</strong></td></tr>")
}

func TestSaveMsaHeuristicCycleCoverOverlapReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "msa_cycle_cover_overlap.md")
	if err := saveMsaHeuristicCycleCoverOverlapReport(path, sampleStructuralAnalyses()); err != nil {
		t.Fatalf("saveMsaHeuristicCycleCoverOverlapReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read MSA heuristic/cycle-cover overlap report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "- **42.86% of MSA heuristic edges are also cycle-cover edges.**")
	assertContains(t, content, "- **33.33% of cycle-cover edges are also MSA heuristic edges.**")
	assertContains(t, content, "- **Found-optimal edge partition: both 3, only MSA heuristic 2, only cycle cover 3, neither 6.**")
	assertContains(t, content, "<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>")
	assertContains(t, content, "<tr><td>a</td><td align=\"right\">50.00</td><td align=\"right\">25.00</td><td align=\"right\">1</td><td align=\"right\">0</td><td align=\"right\">2</td></tr>")
	assertContains(t, content, "<tr><td><strong>Total</strong></td><td align=\"right\">42.86</td><td align=\"right\">33.33</td><td align=\"right\">3</td><td align=\"right\">2</td><td align=\"right\">3</td></tr>")
}

func TestTourMatrixLengthValidatesSingleTourAndCalculatesLength(t *testing.T) {
	distanceMatrix := [][]float64{
		{0, 2, 8},
		{4, 0, 3},
		{5, 6, 0},
	}
	tourMatrix := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}

	length, err := tourMatrixLength(distanceMatrix, tourMatrix)
	if err != nil {
		t.Fatalf("tourMatrixLength returned unexpected error: %v", err)
	}
	if length != 10 {
		t.Fatalf("expected tour length 10, got %.2f", length)
	}

	disconnected := [][]float64{
		{0, 1, 0},
		{1, 0, 0},
		{0, 0, 0},
	}
	if _, err := tourMatrixLength(distanceMatrix, disconnected); err == nil {
		t.Fatal("expected invalid tour matrix to be rejected")
	}
}

func TestSaveGksDeviationReportShowsMsaPatchBiasDeviation(t *testing.T) {
	path := filepath.Join(t.TempDir(), "gks_deviation.md")
	rows := []gksDeviationRow{
		{instance: "sample-a", msaPatchBias: 0.0, tourLength: 110, deviation: 10},
		{instance: "sample-a", msaPatchBias: 0.5, tourLength: 105, deviation: 5},
		{instance: "sample-b", msaPatchBias: 0.0, tourLength: 240, deviation: 20},
		{instance: "sample-b", msaPatchBias: 0.5, tourLength: 260, deviation: 30},
	}

	var builder strings.Builder
	builder.WriteString("# GKS Patching Deviation\n\n")
	writeGksDeviationTable(&builder, rows, []float64{0.0, 0.5})
	if err := os.WriteFile(path, []byte(builder.String()), 0644); err != nil {
		t.Fatalf("failed to write test GKS report: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read GKS report: %v", err)
	}
	content := string(contentBytes)

	assertContains(t, content, "<tr><th>Instance</th><th>MSA patch bias</th><th>Tour length</th><th>Deviation [%]</th></tr>")
	assertContains(t, content, "<tr><td>sample-a</td><td align=\"right\">0.50</td><td align=\"right\">105.00</td><td align=\"right\"><strong>5.00</strong></td></tr>")
	assertContains(t, content, "<tr><td>sample-b</td><td align=\"right\">0.00</td><td align=\"right\">240.00</td><td align=\"right\"><strong>20.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Average</strong></td><td align=\"right\">0.00</td><td align=\"right\"></td><td align=\"right\"><strong>15.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Average</strong></td><td align=\"right\">0.50</td><td align=\"right\"></td><td align=\"right\">17.50</td></tr>")
}

func TestSaveFinalPairwisePerformanceReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "pairwise_performance.md")
	if err := saveFinalPairwisePerformanceReport(path, sampleFinalSummaryRows()); err != nil {
		t.Fatalf("saveFinalPairwisePerformanceReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read pairwise performance report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "<tr><td>Strict MSA vs Baseline</td><td align=\"right\">-1.50</td><td align=\"right\">2</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">+15.00</td></tr>")
	assertContains(t, content, "<tr><td>Cycle-cover MSA patching vs Baseline</td><td align=\"right\">-2.10</td><td align=\"right\">2</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">+20.00</td></tr>")
	assertContains(t, content, "<tr><td>Strict MSA vs Cycle cover</td><td align=\"right\">+0.25</td><td align=\"right\">0</td><td align=\"right\">1</td><td align=\"right\">1</td><td align=\"right\">-5.00</td></tr>")
}

func TestSaveFinalConvergenceSummaryReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "convergence_summary.md")
	if err := saveFinalConvergenceSummaryReport(path, sampleFinalSummaryRows()); err != nil {
		t.Fatalf("saveFinalConvergenceSummaryReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read convergence summary report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "- **Average best-iteration position: Baseline 75.00%, Strict MSA 50.00%, Cycle cover 35.00%, Cycle-cover MSA patching 27.50%.**")
	assertContains(t, content, "<tr><td>a</td><td align=\"right\">80.00</td><td align=\"right\">60.00</td><td></td><td align=\"right\">40.00</td><td align=\"right\"><strong>35.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Average</strong></td><td align=\"right\">75.00</td><td align=\"right\">50.00</td><td></td><td align=\"right\">35.00</td><td align=\"right\"><strong>27.50</strong></td></tr>")
}

func TestSaveStructuralPerformanceLinkReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "structural_performance_link.md")
	if err := saveStructuralPerformanceLinkReport(path, sampleFinalSummaryRows(), sampleStructuralAnalyses()); err != nil {
		t.Fatalf("saveStructuralPerformanceLinkReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read structural performance link report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "<tr><td>Strict MSA</td><td align=\"right\"><strong>71.43</strong></td><td align=\"right\">35.71</td><td align=\"right\">3.00</td><td align=\"right\">20.00</td></tr>")
	assertContains(t, content, "<tr><td>Cycle cover</td><td align=\"right\">66.67</td><td align=\"right\">42.86</td><td align=\"right\">2.75</td><td align=\"right\"><strong>25.00</strong></td></tr>")
	assertContains(t, content, "<tr><td>Cycle-cover MSA patching</td><td align=\"right\">63.64</td><td align=\"right\"><strong>50.00</strong></td><td align=\"right\"><strong>2.40</strong></td><td align=\"right\"><strong>25.00</strong></td></tr>")
}

func TestSelectEvenlySpacedRootIndexes(t *testing.T) {
	actual := selectEvenlySpacedRootIndexes(10, 4)
	expected := []int{0, 3, 6, 9}
	if !reflect.DeepEqual(actual, expected) {
		t.Fatalf("unexpected root selection\nwant: %v\n got: %v", expected, actual)
	}

	actual = selectEvenlySpacedRootIndexes(3, 64)
	expected = []int{0, 1, 2}
	if !reflect.DeepEqual(actual, expected) {
		t.Fatalf("unexpected capped root selection\nwant: %v\n got: %v", expected, actual)
	}
}

func TestBuildPartialMsaHeuristicEdgeSetUsesEligibilityByRoot(t *testing.T) {
	msas := [][][]float64{
		{
			{0, 1, 0, 0},
			{0, 0, 1, 0},
			{0, 0, 0, 1},
			{0, 0, 0, 0},
		},
		emptyMatrix(4),
		{
			{0, 1, 0, 0},
			{0, 0, 0, 1},
			{0, 0, 0, 0},
			{1, 0, 0, 0},
		},
		emptyMatrix(4),
	}

	actual := buildPartialMsaHeuristicEdgeSet(msas, []int{0, 2})
	expected := map[models.Edge]struct{}{
		models.Edge{From: 0, To: 1}: {},
		models.Edge{From: 1, To: 2}: {},
		models.Edge{From: 3, To: 0}: {},
	}
	if !reflect.DeepEqual(actual, expected) {
		t.Fatalf("unexpected partial MSA heuristic edge set\nwant: %v\n got: %v", expected, actual)
	}
}

func emptyMatrix(size int) [][]float64 {
	matrix := make([][]float64, size)
	for i := range matrix {
		matrix[i] = make([]float64, size)
	}

	return matrix
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

	if resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching) != filepath.Join(project.ResultsDirectoryName, "test", "result_cycle_cover_msa_patching.csv") {
		t.Fatalf("unexpected cycle-cover MSA-patching result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching))
	}
}

func TestWithExperimentOutputRootMovesOutputsButKeepsMsaHeuristicCache(t *testing.T) {
	atspData := project.MakeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)
	output := project.WithExperimentOutputRoot(atspData, project.FinalResultsDirectoryName)

	if output.MsaHeuristicDirectoryPath != atspData.MsaHeuristicDirectoryPath {
		t.Fatalf("expected MSA heuristic cache path to stay %s, got %s", atspData.MsaHeuristicDirectoryPath, output.MsaHeuristicDirectoryPath)
	}
	if output.CycleCoverDirectoryPath != atspData.CycleCoverDirectoryPath {
		t.Fatalf("expected cycle-cover cache path to stay %s, got %s", atspData.CycleCoverDirectoryPath, output.CycleCoverDirectoryPath)
	}
	if output.CycleCoverHeatmapPlotPath != atspData.CycleCoverHeatmapPlotPath {
		t.Fatalf("expected cycle-cover heatmap path to stay %s, got %s", atspData.CycleCoverHeatmapPlotPath, output.CycleCoverHeatmapPlotPath)
	}

	expectedResultPath := filepath.Join(project.FinalResultsDirectoryName, "test", project.ResultFileName)
	if output.ResultFilePath != expectedResultPath {
		t.Fatalf("unexpected final result path\nwant: %s\n got: %s", expectedResultPath, output.ResultFilePath)
	}

	if output.OptimalUniqueToursCsvPath != atspData.OptimalUniqueToursCsvPath {
		t.Fatalf("expected solutions path to stay %s, got %s", atspData.OptimalUniqueToursCsvPath, output.OptimalUniqueToursCsvPath)
	}
}

func TestFinalExperimentConfigurationsUseFixedBalancedComparison(t *testing.T) {
	configs := finalExperimentConfigurations()
	expected := []struct {
		heuristic string
		weight    float64
		bias      float64
	}{
		{heuristicBaseline, defaultBaselineHeuristicWeight, 0.0},
		{heuristicStrictMsa, finalStrictMsaHeuristicWeight, 0.0},
		{heuristicRootedMsa, finalRootedMsaHeuristicWeight, 0.0},
		{heuristicCycleCover, finalCycleCoverWeight, 0.0},
		{heuristicCycleCoverMsaPatching, finalCycleCoverMsaPatchingWeight, finalCycleCoverMsaPatchingMsaPatchBias},
	}

	if len(configs) != len(expected) {
		t.Fatalf("expected %d final experiment configurations, got %d", len(expected), len(configs))
	}

	for i, config := range configs {
		if config.Heuristic != expected[i].heuristic {
			t.Fatalf("unexpected final heuristic at %d: want %s got %s", i, expected[i].heuristic, config.Heuristic)
		}
		if len(config.Parameters) != 1 {
			t.Fatalf("expected one final parameter set for %s, got %d", config.Heuristic, len(config.Parameters))
		}

		parameters := config.Parameters[0]
		if parameters.Alpha != defaultExperimentAlpha ||
			parameters.Beta != defaultExperimentBeta ||
			parameters.Rho != defaultExperimentRho ||
			parameters.HeuristicWeight != expected[i].weight ||
			parameters.MsaPatchBias != expected[i].bias {
			t.Fatalf("unexpected final parameters for %s: %+v", config.Heuristic, parameters)
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

func TestSelectFinalExperimentConfigurations(t *testing.T) {
	allConfigurations, err := selectFinalExperimentConfigurations(finalHeuristicAll)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(all) returned error: %v", err)
	}
	if len(allConfigurations) != len(finalExperimentConfigurations()) {
		t.Fatalf("expected all final configurations, got %d", len(allConfigurations))
	}
	if finalConfigurationsAreSparseControls(allConfigurations) {
		t.Fatal("main final all configuration should not be treated as sparse controls")
	}

	controlConfigurations, err := selectFinalExperimentConfigurations(finalHeuristicControls)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(controls) returned error: %v", err)
	}
	if len(controlConfigurations) != len(finalControlExperimentConfigurations()) {
		t.Fatalf("expected all final control configurations, got %d", len(controlConfigurations))
	}
	if !finalConfigurationsAreSparseControls(controlConfigurations) {
		t.Fatal("control configurations should be treated as sparse controls")
	}

	cycleCoverConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if len(cycleCoverConfigurations) != 1 || cycleCoverConfigurations[0].Heuristic != heuristicCycleCover {
		t.Fatalf("expected only cycle-cover configuration, got %+v", cycleCoverConfigurations)
	}

	patchingConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if len(patchingConfigurations) != 1 || patchingConfigurations[0].Heuristic != heuristicCycleCoverMsaPatching {
		t.Fatalf("expected only cycle-cover MSA-patching configuration, got %+v", patchingConfigurations)
	}

	randomSparseConfigurations, err := selectFinalExperimentConfigurations(heuristicRandomSparse)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(random-sparse) returned error: %v", err)
	}
	if len(randomSparseConfigurations) != 1 || randomSparseConfigurations[0].Heuristic != heuristicRandomSparse || len(randomSparseConfigurations[0].Parameters) != len(randomSparseSeeds) {
		t.Fatalf("expected only random-sparse control configuration, got %+v", randomSparseConfigurations)
	}

	shuffledMsaConfigurations, err := selectFinalExperimentConfigurations(heuristicShuffledMsa)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(shuffled-msa) returned error: %v", err)
	}
	if len(shuffledMsaConfigurations) != 1 || shuffledMsaConfigurations[0].Heuristic != heuristicShuffledMsa || len(shuffledMsaConfigurations[0].Parameters) != len(shuffledMsaSeeds) {
		t.Fatalf("expected only shuffled-MSA control configuration, got %+v", shuffledMsaConfigurations)
	}

	if _, err := selectFinalExperimentConfigurations("unknown"); err == nil {
		t.Fatal("expected invalid final heuristic to be rejected")
	}
}

func TestFinalExperimentOutputRootUsesControlsSubdirectoryForSparseControls(t *testing.T) {
	mainConfigurations, err := selectFinalExperimentConfigurations(finalHeuristicAll)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(all) returned error: %v", err)
	}
	if finalExperimentOutputRootForConfigurations(runModeFinal, mainConfigurations) != project.FinalResultsDirectoryName {
		t.Fatalf("main final run should use %s", project.FinalResultsDirectoryName)
	}

	controlConfigurations, err := selectFinalExperimentConfigurations(finalHeuristicControls)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(controls) returned error: %v", err)
	}
	expectedControlsRoot := filepath.Join(filepath.Dir(project.FinalResultsDirectoryName), "controls")
	if finalExperimentOutputRootForConfigurations(runModeFinal, controlConfigurations) != expectedControlsRoot {
		t.Fatalf("final controls should use %s", expectedControlsRoot)
	}
}

func TestFinalConfigurationsUseMsaHeuristicOnlyWhenNeeded(t *testing.T) {
	baselineConfigurations, err := selectFinalExperimentConfigurations(heuristicBaseline)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(baseline) returned error: %v", err)
	}
	if finalConfigurationsUseMsaHeuristic(baselineConfigurations) {
		t.Fatal("baseline-only final run should not require MSA heuristic")
	}

	cycleCoverConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if finalConfigurationsUseMsaHeuristic(cycleCoverConfigurations) {
		t.Fatal("cycle-cover-only final run should not require MSA heuristic")
	}

	strictMsaConfigurations, err := selectFinalExperimentConfigurations(heuristicStrictMsa)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(strict-msa) returned error: %v", err)
	}
	if !finalConfigurationsUseMsaHeuristic(strictMsaConfigurations) {
		t.Fatal("strict MSA final run should require MSA heuristic")
	}

	rootedMsaConfigurations, err := selectFinalExperimentConfigurations(heuristicRootedMsa)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(rooted-msa) returned error: %v", err)
	}
	if !finalConfigurationsUseMsaHeuristic(rootedMsaConfigurations) {
		t.Fatal("rooted MSA final run should require MSA heuristic")
	}

	patchingConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if !finalConfigurationsUseMsaHeuristic(patchingConfigurations) {
		t.Fatal("cycle-cover MSA-patching final run should require MSA heuristic")
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

	baselineConfigurations, err := selectFinalExperimentConfigurations(heuristicBaseline)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(baseline) returned error: %v", err)
	}
	if finalConfigurationsUseCycleCover(baselineConfigurations) {
		t.Fatal("baseline final run should not require cycle-cover cache")
	}

	strictMsaConfigurations, err := selectFinalExperimentConfigurations(heuristicStrictMsa)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(strict-msa) returned error: %v", err)
	}
	if finalConfigurationsUseCycleCover(strictMsaConfigurations) {
		t.Fatal("strict MSA final run should not require cycle-cover cache")
	}

	cycleCoverConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if !finalConfigurationsUseCycleCover(cycleCoverConfigurations) {
		t.Fatal("cycle-cover final run should require cycle-cover cache")
	}

	patchingConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if !finalConfigurationsUseCycleCover(patchingConfigurations) {
		t.Fatal("cycle-cover MSA-patching final run should require cycle-cover cache")
	}
}

func TestSaveFinalHeuristicStatisticsMergesSelectedHeuristicIntoExistingResultCsv(t *testing.T) {
	path := filepath.Join(t.TempDir(), "result.csv")
	existing := []HeuristicExperimentStatistics{
		{
			Heuristic:  heuristicBaseline,
			Statistics: makeTestExperimentStatistics(0.0, 5.0, 10.0),
		},
		{
			Heuristic:  heuristicStrictMsa,
			Statistics: makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, 3.0, 20.0),
		},
		{
			Heuristic:  heuristicCycleCover,
			Statistics: makeTestExperimentStatistics(finalCycleCoverWeight, 7.0, 0.0),
		},
		{
			Heuristic:  "obsolete",
			Statistics: makeTestExperimentStatistics(1.0, 0.5, 100.0),
		},
	}
	if err := saveHeuristicStatistics(path, existing); err != nil {
		t.Fatalf("failed to seed final result CSV: %v", err)
	}

	replacement := []HeuristicExperimentStatistics{
		{
			Heuristic:  heuristicCycleCover,
			Statistics: makeTestExperimentStatistics(finalCycleCoverWeight, 1.5, 40.0),
		},
	}
	cycleCoverConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCover)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover) returned error: %v", err)
	}
	if err := saveFinalHeuristicStatistics(path, replacement, cycleCoverConfigurations); err != nil {
		t.Fatalf("saveFinalHeuristicStatistics returned error: %v", err)
	}

	statistics, err := readHeuristicStatistics(path)
	if err != nil {
		t.Fatalf("failed to read merged final result CSV: %v", err)
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
		t.Fatal("obsolete heuristic row should not be preserved in final result CSV")
	}
}

func TestFinalModesAndExperimentHeuristicsAreValid(t *testing.T) {
	if !isValidRunMode(runModeFinal) {
		t.Fatal("final run mode should be valid")
	}
	if !isValidRunMode(runModeFinal3Opt) {
		t.Fatal("final+3opt run mode should be valid")
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
}

func TestFinal3OptModeUsesSeparateOutputRootAndThreeOpt(t *testing.T) {
	if finalExperimentOutputRoot(runModeFinal3Opt) != project.FinalThreeOptResultsDirectoryName {
		t.Fatalf("final+3opt should use %s", project.FinalThreeOptResultsDirectoryName)
	}
	if finalExperimentOutputRoot(runModeFinal) != project.FinalResultsDirectoryName {
		t.Fatalf("final should use %s", project.FinalResultsDirectoryName)
	}
	if !finalExperimentUsesThreeOpt(runModeFinal3Opt) {
		t.Fatal("final+3opt should enable reduced 3-opt")
	}
	if finalExperimentUsesThreeOpt(runModeFinal) {
		t.Fatal("final should not enable reduced 3-opt")
	}
}

func TestFullFinalRunsTriggerAnalysis(t *testing.T) {
	if !shouldRunAnalysisAfterFinalExperiments(runModeFinal, finalHeuristicAll) {
		t.Fatal("full final run should trigger analysis")
	}
	if !shouldRunAnalysisAfterFinalExperiments(runModeFinal3Opt, finalHeuristicAll) {
		t.Fatal("full final+3opt run should trigger analysis")
	}
	if shouldRunAnalysisAfterFinalExperiments(runModeFinal, heuristicStrictMsa) {
		t.Fatal("single final heuristic should not trigger analysis")
	}
	if shouldRunAnalysisAfterFinalExperiments(runModeFinal, finalHeuristicControls) {
		t.Fatal("final controls should not trigger analysis")
	}
	if shouldRunAnalysisAfterFinalExperiments(runModeExperiment, finalHeuristicAll) {
		t.Fatal("experiment mode should not trigger final analysis")
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

	atspData := project.MakeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, project.ResultsDirectoryName)
	saveStatistics(resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching), heuristicCycleCoverMsaPatching, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.6, 2.0, 20.0),
	})

	if err := runAnalysisMode([]AtspData{atspData}, analysisScopeTuning, []string{heuristicCycleCoverMsaPatching}, 1); err != nil {
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

func TestRunFinalExperimentParametersRejectsInvalidWorkerCount(t *testing.T) {
	matrix := [][]float64{{0, 1}, {1, 0}}
	_, err := runFinalExperimentParameters(
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
	if selected := selectedInstanceSetForMode(runModeFinal, instanceSetTuning, false); selected != instanceSetEvaluation {
		t.Fatalf("final mode without explicit instances should default to %s, got %s", instanceSetEvaluation, selected)
	}
	if selected := selectedInstanceSetForMode(runModeFinal3Opt, instanceSetTuning, false); selected != instanceSetEvaluation {
		t.Fatalf("final+3opt mode without explicit instances should default to %s, got %s", instanceSetEvaluation, selected)
	}
	if selected := selectedInstanceSetForMode(runModeRebuildCache, instanceSetTuning, false); selected != instanceSetAllKnown {
		t.Fatalf("rebuild-cache mode without explicit instances should default to %s, got %s", instanceSetAllKnown, selected)
	}
	if selected := selectedInstanceSetForMode(runModeFinal, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("final mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeRebuildCache, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("rebuild-cache mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeExperiment, instanceSetTuning, false); selected != instanceSetTuning {
		t.Fatalf("experiment mode should keep requested instances, got %s", selected)
	}
}

func sampleStructuralAnalyses() []structure.InstanceAnalysis {
	return []structure.InstanceAnalysis{
		{
			Instance:  "b",
			Dimension: 5,
			Metrics: structure.InstanceMetrics{
				FoundOptimalTourCount:       2,
				UniqueFoundOptimalEdgeCount: 10,
				HighMsaHeuristicMetrics: structure.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 4,
					Precision:        0.8,
					Recall:           0.4,
				},
				CycleCoverMetrics: structure.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 3,
					Precision:        0.6,
					Recall:           0.3,
				},
				CycleCoverMsaPatchingMetrics: structure.EdgeSetMetrics{
					EdgeCount:        6,
					OptimalEdgeCount: 4,
					Precision:        4.0 / 6.0,
					Recall:           0.4,
				},
				CycleCoverHighMsaHeuristicEdges:             2,
				OptimalEdgesInCycleCoverAndHighMsaHeuristic: 2,
				OptimalEdgesInHighMsaHeuristicNotCycleCover: 2,
				OptimalEdgesInCycleCoverNotHighMsaHeuristic: 1,
				OptimalEdgesInNeitherCycleCoverNorHigh:      5,
			},
		},
		{
			Instance:  "a",
			Dimension: 4,
			Metrics: structure.InstanceMetrics{
				FoundOptimalTourCount:       1,
				UniqueFoundOptimalEdgeCount: 4,
				HighMsaHeuristicMetrics: structure.EdgeSetMetrics{
					EdgeCount:        2,
					OptimalEdgeCount: 1,
					Precision:        0.5,
					Recall:           0.25,
				},
				CycleCoverMetrics: structure.EdgeSetMetrics{
					EdgeCount:        4,
					OptimalEdgeCount: 3,
					Precision:        0.75,
					Recall:           0.75,
				},
				CycleCoverMsaPatchingMetrics: structure.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 3,
					Precision:        0.6,
					Recall:           0.75,
				},
				CycleCoverHighMsaHeuristicEdges:             1,
				OptimalEdgesInCycleCoverAndHighMsaHeuristic: 1,
				OptimalEdgesInHighMsaHeuristicNotCycleCover: 0,
				OptimalEdgesInCycleCoverNotHighMsaHeuristic: 2,
				OptimalEdgesInNeitherCycleCoverNorHigh:      1,
			},
		},
	}
}

func sampleFinalSummaryRows() []finalResultsSummaryRow {
	return []finalResultsSummaryRow{
		{
			instance: "a",
			metrics: map[string]finalResultsSummaryMetric{
				heuristicBaseline: {
					AverageMinDeviation:  4.0,
					SuccessRate:          10.0,
					AverageBestIteration: 80.0,
					Iterations:           100,
				},
				heuristicStrictMsa: {
					AverageMinDeviation:  2.5,
					SuccessRate:          20.0,
					AverageBestIteration: 60.0,
					Iterations:           100,
				},
				heuristicCycleCover: {
					AverageMinDeviation:  2.0,
					SuccessRate:          30.0,
					AverageBestIteration: 40.0,
					Iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  1.8,
					SuccessRate:          25.0,
					AverageBestIteration: 35.0,
					Iterations:           100,
				},
			},
		},
		{
			instance: "b",
			metrics: map[string]finalResultsSummaryMetric{
				heuristicBaseline: {
					AverageMinDeviation:  5.0,
					SuccessRate:          0.0,
					AverageBestIteration: 70.0,
					Iterations:           100,
				},
				heuristicStrictMsa: {
					AverageMinDeviation:  3.5,
					SuccessRate:          20.0,
					AverageBestIteration: 40.0,
					Iterations:           100,
				},
				heuristicCycleCover: {
					AverageMinDeviation:  3.5,
					SuccessRate:          20.0,
					AverageBestIteration: 30.0,
					Iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  3.0,
					SuccessRate:          25.0,
					AverageBestIteration: 20.0,
					Iterations:           100,
				},
			},
		},
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
	statistics := makeTestExperimentStatistics(finalStrictMsaHeuristicWeight, averageBestDeviation, successRate)
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

	saveStatistics(resultFilePathForHeuristic(atspData, heuristic), heuristic, statistics)
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

func assertPathExists(t *testing.T, path string) {
	t.Helper()
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("expected %s to exist: %v", path, err)
	}
}
