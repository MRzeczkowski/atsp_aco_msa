package main

import (
	"atsp_aco_msa/modules/analysis/structuralComparison"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/parsing"
	"os"
	"path/filepath"
	"reflect"
	"strconv"
	"strings"
	"sync/atomic"
	"testing"
	"time"
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

	if statistics[0].heuristicWeight != 0.50 || statistics[0].msaPatchBias != 0.0 || statistics[0].successRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f msaPatchBias=%f successRate=%f", statistics[0].heuristicWeight, statistics[0].msaPatchBias, statistics[0].successRate)
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
	if statistics[0].heuristicWeight != 0.50 || statistics[0].msaPatchBias != 0.25 || statistics[0].successRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f msaPatchBias=%f successRate=%f", statistics[0].heuristicWeight, statistics[0].msaPatchBias, statistics[0].successRate)
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
	if statistics[0].heuristicWeight != 0.90 || statistics[0].randomSeed != 2 || statistics[0].successRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f randomSeed=%d successRate=%f", statistics[0].heuristicWeight, statistics[0].randomSeed, statistics[0].successRate)
	}
}

func TestSaveRandomSparseControlReportComparesMsaAgainstRandomSparse(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	finalRoot := filepath.Join(root, "final")
	controlsRoot := finalControlsResultsRootPath(finalRoot)
	first := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	finalFirst := withExperimentOutputRoot(first, finalRoot)
	finalSecond := withExperimentOutputRoot(second, finalRoot)
	controlFirst := withExperimentOutputRoot(first, controlsRoot)
	controlSecond := withExperimentOutputRoot(second, controlsRoot)

	if err := saveHeuristicStatistics(finalFirst.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 2.0, 10.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.8, 0.5, 100.0)},
	}); err != nil {
		t.Fatalf("failed to write first final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlFirst, heuristicRandomSparse), heuristicRandomSparse, []ExperimentsDataStatistics{
		makeTestRandomSparseExperimentStatistics(1, 4.0, 0.0),
		makeTestRandomSparseExperimentStatistics(2, 5.0, 0.0),
		makeTestRandomSparseExperimentStatistics(3, 6.0, 0.0),
	})

	if err := saveHeuristicStatistics(finalSecond.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 4.0, 0.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.8, 0.5, 100.0)},
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
	assertContains(t, content, "MSA had lower average best deviation than the random-sparse mean in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: MSA 3.00%, random sparse 4.00%, delta -1.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func TestSaveDistanceRankedSparseControlReportComparesMsaAgainstDistanceRankedSparse(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	finalRoot := filepath.Join(root, "final")
	controlsRoot := finalControlsResultsRootPath(finalRoot)
	first := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	finalFirst := withExperimentOutputRoot(first, finalRoot)
	finalSecond := withExperimentOutputRoot(second, finalRoot)
	controlFirst := withExperimentOutputRoot(first, controlsRoot)
	controlSecond := withExperimentOutputRoot(second, controlsRoot)

	if err := saveHeuristicStatistics(finalFirst.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 2.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write first final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlFirst, heuristicDistanceRankedSparse), heuristicDistanceRankedSparse, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(finalMsaHeuristicWeight, 4.0, 0.0),
	})

	if err := saveHeuristicStatistics(finalSecond.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 4.0, 0.0)},
	}); err != nil {
		t.Fatalf("failed to write second final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlSecond, heuristicDistanceRankedSparse), heuristicDistanceRankedSparse, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(finalMsaHeuristicWeight, 2.0, 20.0),
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
	assertContains(t, content, "MSA had lower average best deviation than the distance-ranked sparse control in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: MSA 3.00%, distance-ranked sparse 3.00%, delta +0.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func TestSaveShuffledMsaControlReportComparesMsaAgainstShuffledMsa(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	finalRoot := filepath.Join(root, "final")
	controlsRoot := finalControlsResultsRootPath(finalRoot)
	first := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	finalFirst := withExperimentOutputRoot(first, finalRoot)
	finalSecond := withExperimentOutputRoot(second, finalRoot)
	controlFirst := withExperimentOutputRoot(first, controlsRoot)
	controlSecond := withExperimentOutputRoot(second, controlsRoot)

	if err := saveHeuristicStatistics(finalFirst.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 2.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write first final MSA result: %v", err)
	}
	saveStatistics(resultFilePathForHeuristic(controlFirst, heuristicShuffledMsa), heuristicShuffledMsa, []ExperimentsDataStatistics{
		makeTestRandomSparseExperimentStatistics(1, 4.0, 0.0),
		makeTestRandomSparseExperimentStatistics(2, 5.0, 0.0),
		makeTestRandomSparseExperimentStatistics(3, 6.0, 0.0),
	})

	if err := saveHeuristicStatistics(finalSecond.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 4.0, 0.0)},
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
	assertContains(t, content, "MSA had lower average best deviation than the shuffled MSA mean in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: MSA 3.00%, shuffled MSA 4.00%, delta -1.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

// Temporary msa-impact prototype tests. Remove with the msa-impact pipeline.
func TestSaveMsaImpactSummaryComparesMsaAgainstBaseline(t *testing.T) {
	originalMsaImpactResultsDirectoryName := msaImpactResultsDirectoryName
	msaImpactResultsDirectoryName = t.TempDir()
	t.Cleanup(func() {
		msaImpactResultsDirectoryName = originalMsaImpactResultsDirectoryName
	})

	sourceRoot := filepath.Join(t.TempDir(), "results")
	first := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	impactFirst := withExperimentOutputRoot(first, msaImpactResultsDirectoryName)
	impactSecond := withExperimentOutputRoot(second, msaImpactResultsDirectoryName)

	if err := saveHeuristicStatistics(impactFirst.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicBaseline, statistics: makeTestExperimentStatistics(defaultBaselineHeuristicWeight, 5.0, 10.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 3.0, 20.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.8, 4.0, 0.0)},
	}); err != nil {
		t.Fatalf("failed to write first impact result: %v", err)
	}
	if err := saveHeuristicStatistics(impactSecond.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicBaseline, statistics: makeTestExperimentStatistics(defaultBaselineHeuristicWeight, 2.0, 30.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 4.0, 10.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.2, 1.0, 50.0)},
	}); err != nil {
		t.Fatalf("failed to write second impact result: %v", err)
	}

	reportPath := filepath.Join(msaImpactResultsDirectoryName, "msa_impact_summary.md")
	if err := saveMsaImpactSummary(reportPath, []AtspData{second, first}); err != nil {
		t.Fatalf("saveMsaImpactSummary returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read MSA impact summary: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "# MSA Impact Summary")
	assertContains(t, content, "MSA heuristic improved average best deviation in 2/2 instances, tied in 0, and was worse in 0.")
	assertContains(t, content, "Mean average best deviation: baseline 3.50%, MSA heuristic 2.00%, delta -1.50 pp.")
	assertContains(t, content, "Best MSA heuristic weights: 0.20 (1), 0.40 (1).")
	assertContains(t, content, "<tr><td>sample-a</td><td align=\"right\">5.00</td><td align=\"right\">3.00</td><td align=\"right\">0.40</td><td align=\"right\">-2.00</td><td align=\"right\">10.00</td><td align=\"right\">20.00</td></tr>")
	assertContains(t, content, "<tr><td>sample-b</td><td align=\"right\">2.00</td><td align=\"right\">1.00</td><td align=\"right\">0.20</td><td align=\"right\">-1.00</td><td align=\"right\">30.00</td><td align=\"right\">50.00</td></tr>")
}

func TestSaveMsaImpactStructureReportCombinesStructureWithPerformance(t *testing.T) {
	originalMsaImpactResultsDirectoryName := msaImpactResultsDirectoryName
	msaImpactResultsDirectoryName = t.TempDir()
	t.Cleanup(func() {
		msaImpactResultsDirectoryName = originalMsaImpactResultsDirectoryName
	})

	sourceRoot := filepath.Join(t.TempDir(), "results")
	msaRoot := filepath.Join(t.TempDir(), "msa")
	first := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{
		{0, 1, 1, 1},
		{1, 0, 1, 1},
		{1, 1, 0, 1},
		{1, 1, 1, 0},
	}, 4, sourceRoot)
	second := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}, 3, sourceRoot)
	first.msaHeuristicDirectoryPath = filepath.Join(msaRoot, "sample-a")
	second.msaHeuristicDirectoryPath = filepath.Join(msaRoot, "sample-b")

	writeTestMsaHeuristicMatrix(t, first.msaHeuristicDirectoryPath, [][]float64{
		{0, 3, 0, 0},
		{0, 0, 3, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
	})
	writeTestRootMsaMatrix(t, first.msaHeuristicDirectoryPath, 0, [][]float64{
		{0, 1, 1, 1},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
	})
	writeTestRootMsaMatrix(t, first.msaHeuristicDirectoryPath, 1, [][]float64{
		{0, 0, 1, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 0, 0},
	})
	writeTestMsaHeuristicMatrix(t, second.msaHeuristicDirectoryPath, [][]float64{
		{0, 2, 0},
		{0, 0, 2},
		{2, 0, 0},
	})
	writeTestRootMsaMatrix(t, second.msaHeuristicDirectoryPath, 0, [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, 0},
	})

	impactFirst := withExperimentOutputRoot(first, msaImpactResultsDirectoryName)
	impactSecond := withExperimentOutputRoot(second, msaImpactResultsDirectoryName)
	if err := saveHeuristicStatistics(impactFirst.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicBaseline, statistics: makeTestExperimentStatistics(defaultBaselineHeuristicWeight, 5.0, 10.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 3.0, 20.0)},
	}); err != nil {
		t.Fatalf("failed to write first impact result: %v", err)
	}
	if err := saveHeuristicStatistics(impactSecond.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicBaseline, statistics: makeTestExperimentStatistics(defaultBaselineHeuristicWeight, 2.0, 30.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 4.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write second impact result: %v", err)
	}

	reportPath := filepath.Join(msaImpactResultsDirectoryName, "msa_impact_structure.md")
	if err := saveMsaImpactStructureReport(reportPath, []AtspData{first, second}); err != nil {
		t.Fatalf("saveMsaImpactStructureReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read MSA impact structure report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "# MSA Impact Structure")
	assertContains(t, content, "Average boosted-edge ratio against n-1: 1.08.")
	assertContains(t, content, "Average missing outgoing vertices: 1.00; average missing incoming vertices: 1.00.")
	assertContains(t, content, "Average missing endpoints: improved cases 4.00, tied/worse cases 0.00.")
	assertContains(t, content, "<tr><td>sample-a</td><td align=\"right\">4</td><td align=\"right\">2</td><td align=\"right\">0.67</td><td align=\"right\">2</td><td align=\"right\">2</td><td align=\"right\">-2.00</td><td align=\"right\">+10.00</td></tr>")
	assertContains(t, content, "<tr><td>sample-b</td><td align=\"right\">3</td><td align=\"right\">3</td><td align=\"right\">1.50</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">+2.00</td><td align=\"right\">-20.00</td></tr>")
}

func TestReadMsaHeuristicMatrixForResultRootUsesDefaultMsaForMsaImpact(t *testing.T) {
	root := t.TempDir()
	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}, 3, filepath.Join(root, "results"))
	atspData.msaHeuristicDirectoryPath = filepath.Join(root, "msa", "sample")

	compositeMsa := [][]float64{
		{0, 2, 0},
		{0, 0, 2},
		{2, 0, 0},
	}
	writeTestMsaHeuristicMatrix(t, atspData.msaHeuristicDirectoryPath, compositeMsa)
	writeTestRootMsaMatrix(t, atspData.msaHeuristicDirectoryPath, 0, [][]float64{
		{0, 1, 1},
		{0, 0, 0},
		{0, 0, 0},
	})
	writeTestRootMsaMatrix(t, atspData.msaHeuristicDirectoryPath, 1, [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, 0},
	})

	normalMatrix, err := readMsaHeuristicMatrixForResultRoot(atspData, heuristicMsaHeuristic, finalResultsDirectoryName)
	if err != nil {
		t.Fatalf("read normal MSA matrix: %v", err)
	}
	if !reflect.DeepEqual(normalMatrix, compositeMsa) {
		t.Fatalf("expected normal run to use composite MSA, got %v", normalMatrix)
	}

	impactMatrix, err := readMsaHeuristicMatrixForResultRoot(atspData, heuristicMsaHeuristic, msaImpactResultsDirectoryName)
	if err != nil {
		t.Fatalf("read MSA impact matrix: %v", err)
	}
	if !reflect.DeepEqual(impactMatrix, compositeMsa) {
		t.Fatalf("expected MSA impact to use composite MSA, got %v", impactMatrix)
	}

	controlMatrix, err := readMsaHeuristicMatrixForResultRoot(atspData, heuristicRandomSparse, msaImpactControlsDirectoryName)
	if err != nil {
		t.Fatalf("read MSA impact control matrix: %v", err)
	}
	if !reflect.DeepEqual(controlMatrix, compositeMsa) {
		t.Fatalf("expected MSA impact controls to use composite MSA, got %v", controlMatrix)
	}
}

func TestReadFinalMsaHeuristicControlMetricUsesBestWeightForMsaImpact(t *testing.T) {
	originalMsaImpactResultsDirectoryName := msaImpactResultsDirectoryName
	root := t.TempDir()
	msaImpactResultsDirectoryName = filepath.Join(root, "msa-impact")
	t.Cleanup(func() {
		msaImpactResultsDirectoryName = originalMsaImpactResultsDirectoryName
	})

	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, filepath.Join(root, "source"))
	impactAtspData := withExperimentOutputRoot(atspData, msaImpactResultsDirectoryName)
	if err := saveHeuristicStatistics(impactAtspData.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.4, 5.0, 0.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.8, 2.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write MSA impact result: %v", err)
	}

	impactMetric, ok, err := readFinalMsaHeuristicControlMetric(atspData, msaImpactResultsDirectoryName)
	if err != nil {
		t.Fatalf("read MSA impact control metric: %v", err)
	}
	if !ok || impactMetric.heuristicWeight != 0.8 || impactMetric.averageMinDeviation != 2.0 {
		t.Fatalf("expected best MSA impact metric, got metric=%+v ok=%t", impactMetric, ok)
	}

	finalRoot := filepath.Join(root, "final")
	finalAtspData := withExperimentOutputRoot(atspData, finalRoot)
	if err := saveHeuristicStatistics(finalAtspData.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 6.0, 0.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.8, 1.0, 100.0)},
	}); err != nil {
		t.Fatalf("failed to write final result: %v", err)
	}

	finalMetric, ok, err := readFinalMsaHeuristicControlMetric(atspData, finalRoot)
	if err != nil {
		t.Fatalf("read final control metric: %v", err)
	}
	if !ok || finalMetric.heuristicWeight != finalMsaHeuristicWeight || finalMetric.averageMinDeviation != 6.0 {
		t.Fatalf("expected fixed-weight final metric, got metric=%+v ok=%t", finalMetric, ok)
	}
}

func TestSaveMsaImpactControlWeightSummaryReportsPValuePerWeight(t *testing.T) {
	originalMsaImpactResultsDirectoryName := msaImpactResultsDirectoryName
	originalMsaImpactControlsDirectoryName := msaImpactControlsDirectoryName
	root := t.TempDir()
	msaImpactResultsDirectoryName = filepath.Join(root, "results")
	msaImpactControlsDirectoryName = filepath.Join(root, "controls")
	t.Cleanup(func() {
		msaImpactResultsDirectoryName = originalMsaImpactResultsDirectoryName
		msaImpactControlsDirectoryName = originalMsaImpactControlsDirectoryName
	})

	sourceRoot := filepath.Join(t.TempDir(), "source")
	first := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	impactFirst := withExperimentOutputRoot(first, msaImpactResultsDirectoryName)
	impactSecond := withExperimentOutputRoot(second, msaImpactResultsDirectoryName)
	controlFirst := withExperimentOutputRoot(first, msaImpactControlsDirectoryName)
	controlSecond := withExperimentOutputRoot(second, msaImpactControlsDirectoryName)

	if err := saveHeuristicStatistics(impactFirst.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.4, 2.0, 10.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.8, 5.0, 0.0)},
	}); err != nil {
		t.Fatalf("failed to write first MSA impact result: %v", err)
	}
	if err := saveHeuristicStatistics(impactSecond.resultFilePath, []HeuristicExperimentStatistics{
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.4, 1.0, 20.0)},
		{heuristic: heuristicMsaHeuristic, statistics: makeTestExperimentStatistics(0.8, 1.0, 30.0)},
	}); err != nil {
		t.Fatalf("failed to write second MSA impact result: %v", err)
	}

	writeTestControlStatisticsForWeightSummary(t, controlFirst, heuristicShuffledMsa)
	writeTestControlStatisticsForWeightSummary(t, controlSecond, heuristicShuffledMsa)
	saveStatistics(resultFilePathForHeuristic(controlFirst, heuristicDistanceRankedSparse), heuristicDistanceRankedSparse, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.4, 4.0, 0.0),
		makeTestExperimentStatistics(0.8, 3.0, 0.0),
	})
	saveStatistics(resultFilePathForHeuristic(controlSecond, heuristicDistanceRankedSparse), heuristicDistanceRankedSparse, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.4, 2.0, 0.0),
		makeTestExperimentStatistics(0.8, 4.0, 0.0),
	})

	reportPath := filepath.Join(msaImpactControlsDirectoryName, "msa_weight_control_summary.md")
	saved, err := saveMsaImpactControlWeightSummary(reportPath, []AtspData{second, first})
	if err != nil {
		t.Fatalf("saveMsaImpactControlWeightSummary returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected MSA impact control weight summary to be saved")
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read MSA impact control weight summary: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "# MSA Weight Control Summary")
	if strings.Contains(content, "## Random Sparse") {
		t.Fatal("MSA impact control summary should not include the removed random-sparse control")
	}
	assertContains(t, content, "## Distance-ranked Sparse")
	assertContains(t, content, "## Shuffled MSA")
	assertContains(t, content, "<tr><td align=\"right\">0.40</td><td align=\"right\">1</td><td align=\"right\">1</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">1.000000")
	assertContains(t, content, "<tr><td align=\"right\">0.80</td><td align=\"right\">1</td><td align=\"right\">1</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">1.000000")
}

func TestSaveMsaDistanceRankedCategoryReportComparesEdgeCategories(t *testing.T) {
	root := t.TempDir()
	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{
		{0, 1, 2, 100},
		{100, 0, 50, 100},
		{100, 100, 0, 60},
		{3, 100, 100, 0},
	}, 10, filepath.Join(root, "results"))
	atspData.msaHeuristicDirectoryPath = filepath.Join(root, "msa", "sample")
	atspData.optimalUniqueToursCsvPath = filepath.Join(root, "solutions", "sample", "solutions.csv")

	writeTestMsaHeuristicMatrix(t, atspData.msaHeuristicDirectoryPath, [][]float64{
		{0, 3, 0, 0},
		{0, 0, 3, 0},
		{0, 0, 0, 3},
		{0, 0, 0, 0},
	})
	writeTestRootMsaMatrix(t, atspData.msaHeuristicDirectoryPath, 0, [][]float64{
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
		{0, 0, 0, 0},
	})
	writeTestRootMsaMatrix(t, atspData.msaHeuristicDirectoryPath, 1, [][]float64{
		{0, 1, 1, 1},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
	})
	writeTestOptimalToursCsv(t, atspData.optimalUniqueToursCsvPath, []string{"[0,1,2,3]"})

	reportPath := filepath.Join(root, "msa_distance_ranked_edge_categories.md")
	if err := saveMsaDistanceRankedCategoryReport(reportPath, []AtspData{atspData}); err != nil {
		t.Fatalf("saveMsaDistanceRankedCategoryReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read MSA distance-ranked category report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "# MSA vs Distance-ranked Edge Categories")
	assertContains(t, content, "MSA-only precision 100.00% vs distance-only precision 50.00%.")
	assertContains(t, content, "MSA-only recall 50.00% vs distance-only recall 25.00%.")
	assertContains(t, content, "<tr><td>Both</td><td align=\"right\">1</td><td align=\"right\">1</td><td align=\"right\">100.00</td><td align=\"right\">25.00</td></tr>")
	assertContains(t, content, "<tr><td>MSA only</td><td align=\"right\">2</td><td align=\"right\">2</td><td align=\"right\">100.00</td><td align=\"right\">50.00</td></tr>")
	assertContains(t, content, "<tr><td>Distance-ranked only</td><td align=\"right\">2</td><td align=\"right\">1</td><td align=\"right\">50.00</td><td align=\"right\">25.00</td></tr>")
}

func TestSaveHeuristicStatisticsWritesSingleComparisonCsv(t *testing.T) {
	path := filepath.Join(t.TempDir(), "result.csv")
	rows := []HeuristicExperimentStatistics{
		{
			heuristic: heuristicBaseline,
			statistics: ExperimentsDataStatistics{
				ExperimentParameters: ExperimentParameters{
					alpha:           1.0,
					beta:            2.0,
					rho:             0.8,
					heuristicWeight: 0.0,
					iterations:      500,
				},
				averageBestDeviation: 3.0,
			},
		},
		{
			heuristic: heuristicCycleCover,
			statistics: ExperimentsDataStatistics{
				ExperimentParameters: ExperimentParameters{
					alpha:           1.0,
					beta:            2.0,
					rho:             0.8,
					heuristicWeight: 0.8,
					iterations:      500,
				},
				averageBestDeviation: 1.0,
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
	firstAtspData := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	secondAtspData := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	rows := []HeuristicExperimentStatistics{
		{
			heuristic: heuristicBaseline,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 4.25,
				successRate:          10.0,
			},
		},
		{
			heuristic: heuristicMsaHeuristic,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 2.50,
				successRate:          20.0,
			},
		},
		{
			heuristic: heuristicCycleCover,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 1.75,
				successRate:          30.0,
			},
		},
		{
			heuristic: heuristicCycleCoverMsaPatching,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 2.00,
				successRate:          25.0,
			},
		},
	}
	secondRows := []HeuristicExperimentStatistics{
		{
			heuristic: heuristicBaseline,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 6.75,
				successRate:          20.0,
			},
		},
		{
			heuristic: heuristicMsaHeuristic,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 3.50,
				successRate:          40.0,
			},
		},
		{
			heuristic: heuristicCycleCover,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 2.25,
				successRate:          60.0,
			},
		},
		{
			heuristic: heuristicCycleCoverMsaPatching,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 2.50,
				successRate:          55.0,
			},
		},
	}

	if err := saveHeuristicStatistics(firstAtspData.resultFilePath, rows); err != nil {
		t.Fatalf("saveHeuristicStatistics returned unexpected error: %v", err)
	}
	if err := saveHeuristicStatistics(secondAtspData.resultFilePath, secondRows); err != nil {
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
		"- **Best-or-tied average best deviation counts: Baseline 0/2, MSA heuristic 0/2, Cycle cover 2/2, Cycle-cover MSA patching 0/2.**",
		"",
		"<table>",
		"<thead>",
		"<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">Baseline</th><th colspan=\"2\">MSA heuristic</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>",
		"<tr><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th></tr>",
		"</thead>",
		"<tbody>",
		"<tr><td>sample-a</td><td align=\"right\">4.25</td><td align=\"right\">10.00</td><td align=\"right\">2.50</td><td align=\"right\">20.00</td><td align=\"right\"><strong>1.75</strong></td><td align=\"right\"><strong>30.00</strong></td><td align=\"right\">2.00</td><td align=\"right\">25.00</td></tr>",
		"<tr><td>sample-b</td><td align=\"right\">6.75</td><td align=\"right\">20.00</td><td align=\"right\">3.50</td><td align=\"right\">40.00</td><td align=\"right\"><strong>2.25</strong></td><td align=\"right\"><strong>60.00</strong></td><td align=\"right\">2.50</td><td align=\"right\">55.00</td></tr>",
		"<tr><td><strong>Average</strong></td><td align=\"right\">5.50</td><td align=\"right\">15.00</td><td align=\"right\">3.00</td><td align=\"right\">30.00</td><td align=\"right\"><strong>2.00</strong></td><td align=\"right\"><strong>45.00</strong></td><td align=\"right\">2.25</td><td align=\"right\">40.00</td></tr>",
		"</tbody>",
		"</table>",
	}
	if !reflect.DeepEqual(lines, expected) {
		t.Fatalf("unexpected final summary Markdown\nwant: %v\n got: %v", expected, lines)
	}
}

func TestRunFinalResultsAnalysisReadsExistingFinalResults(t *testing.T) {
	resultsRoot := t.TempDir()
	oldFinalResultsDirectoryName := finalResultsDirectoryName
	finalResultsDirectoryName = filepath.Join(resultsRoot, "final")
	defer func() {
		finalResultsDirectoryName = oldFinalResultsDirectoryName
	}()

	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	finalAtspData := withExperimentOutputRoot(atspData, finalResultsDirectoryName)
	rows := []HeuristicExperimentStatistics{
		{
			heuristic: heuristicBaseline,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 4.25,
				successRate:          10.0,
			},
		},
		{
			heuristic: heuristicMsaHeuristic,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 2.50,
				successRate:          20.0,
			},
		},
		{
			heuristic: heuristicCycleCover,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 1.75,
				successRate:          30.0,
			},
		},
	}
	if err := saveHeuristicStatistics(finalAtspData.resultFilePath, rows); err != nil {
		t.Fatalf("saveHeuristicStatistics returned unexpected error: %v", err)
	}

	summaryPath, _, saved, err := runFinalResultsAnalysis([]AtspData{atspData}, nil, finalResultsDirectoryName)
	if err != nil {
		t.Fatalf("runFinalResultsAnalysis returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected final results summary to be saved")
	}
	if summaryPath != filepath.Join(finalResultsDirectoryName, "summary.md") {
		t.Fatalf("unexpected summary path: %s", summaryPath)
	}
	if _, err := os.Stat(summaryPath); err != nil {
		t.Fatalf("expected final results summary file to exist: %v", err)
	}
}

func TestRunFinalResultsAnalysisUsesProvidedResultsRoot(t *testing.T) {
	resultsRoot := t.TempDir()
	finalThreeOptRoot := filepath.Join(resultsRoot, "final_3opt")
	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	finalThreeOptAtspData := withExperimentOutputRoot(atspData, finalThreeOptRoot)
	rows := []HeuristicExperimentStatistics{
		{
			heuristic: heuristicBaseline,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 3.00,
				successRate:          10.0,
			},
		},
		{
			heuristic: heuristicMsaHeuristic,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 2.00,
				successRate:          20.0,
			},
		},
		{
			heuristic: heuristicCycleCover,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 1.00,
				successRate:          30.0,
			},
		},
	}
	if err := saveHeuristicStatistics(finalThreeOptAtspData.resultFilePath, rows); err != nil {
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
					averageMinDeviation:  10.0,
					successRate:          10.0,
					averageBestIteration: 50.0,
					iterations:           100,
				},
				heuristicMsaHeuristic: {
					averageMinDeviation:  8.0,
					successRate:          12.0,
					averageBestIteration: 60.0,
					iterations:           100,
				},
				heuristicCycleCover: {
					averageMinDeviation:  7.0,
					successRate:          11.0,
					averageBestIteration: 30.0,
					iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					averageMinDeviation:  7.5,
					successRate:          13.0,
					averageBestIteration: 35.0,
					iterations:           100,
				},
			},
		},
	}
	finalThreeOptRows := []finalResultsSummaryRow{
		{
			instance: "sample",
			metrics: map[string]finalResultsSummaryMetric{
				heuristicBaseline: {
					averageMinDeviation:  1.0,
					successRate:          50.0,
					averageBestIteration: 20.0,
					iterations:           100,
				},
				heuristicMsaHeuristic: {
					averageMinDeviation:  0.9,
					successRate:          51.0,
					averageBestIteration: 19.0,
					iterations:           100,
				},
				heuristicCycleCover: {
					averageMinDeviation:  0.8,
					successRate:          52.0,
					averageBestIteration: 19.0,
					iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					averageMinDeviation:  0.85,
					successRate:          53.0,
					averageBestIteration: 18.0,
					iterations:           100,
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
	assertContains(t, content, "<tr><td>MSA heuristic</td><td align=\"right\">+2.00</td><td align=\"right\">+0.10</td><td align=\"right\">5.00</td></tr>")
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
	assertContains(t, content, "<tr><td>MSA heuristic vs Baseline</td><td align=\"right\">-1.50</td><td align=\"right\">2</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">+15.00</td></tr>")
	assertContains(t, content, "<tr><td>Cycle-cover MSA patching vs Baseline</td><td align=\"right\">-2.10</td><td align=\"right\">2</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">+20.00</td></tr>")
	assertContains(t, content, "<tr><td>MSA heuristic vs Cycle cover</td><td align=\"right\">+0.25</td><td align=\"right\">0</td><td align=\"right\">1</td><td align=\"right\">1</td><td align=\"right\">-5.00</td></tr>")
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
	assertContains(t, content, "- **Average best-iteration position: Baseline 75.00%, MSA heuristic 50.00%, Cycle cover 35.00%, Cycle-cover MSA patching 27.50%.**")
	assertContains(t, content, "<tr><td>a</td><td align=\"right\">80.00</td><td align=\"right\">60.00</td><td align=\"right\">40.00</td><td align=\"right\"><strong>35.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Average</strong></td><td align=\"right\">75.00</td><td align=\"right\">50.00</td><td align=\"right\">35.00</td><td align=\"right\"><strong>27.50</strong></td></tr>")
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
	assertContains(t, content, "<tr><td>MSA heuristic</td><td align=\"right\"><strong>71.43</strong></td><td align=\"right\">35.71</td><td align=\"right\">3.00</td><td align=\"right\">20.00</td></tr>")
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
		{1, 1, 1},
		{1, 1, 1},
		{1, 1, 1},
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
		{1, 1},
		{1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected baseline modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildMinimumCycleCoverMatrix(t *testing.T) {
	matrix := [][]float64{
		{0, 5, 1},
		{1, 0, 5},
		{5, 1, 0},
	}

	cycleCover, cost, err := buildMinimumCycleCoverMatrix(matrix)
	if err != nil {
		t.Fatalf("buildMinimumCycleCoverMatrix returned unexpected error: %v", err)
	}

	expected := [][]float64{
		{0, 0, 1},
		{1, 0, 0},
		{0, 1, 0},
	}
	if !reflect.DeepEqual(cycleCover, expected) {
		t.Fatalf("unexpected cycle-cover matrix\nwant: %v\n got: %v", expected, cycleCover)
	}
	if cost != 3 {
		t.Fatalf("expected cycle-cover cost 3, got %f", cost)
	}
}

func TestHeuristicSpecificPathsKeepMsaHeuristicBaselinePaths(t *testing.T) {
	atspData := makeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)

	if resultFilePathForHeuristic(atspData, heuristicBaseline) != filepath.Join(resultsDirectoryName, "test", "result_baseline.csv") {
		t.Fatalf("unexpected baseline result path: %s", resultFilePathForHeuristic(atspData, heuristicBaseline))
	}

	if resultFilePathForHeuristic(atspData, heuristicMsaHeuristic) != filepath.Join(resultsDirectoryName, "test", resultFileName) {
		t.Fatalf("MSA heuristic result path should keep the existing baseline location")
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCover) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover.csv") {
		t.Fatalf("unexpected cycle-cover result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCover))
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover_msa_patching.csv") {
		t.Fatalf("unexpected cycle-cover MSA-patching result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching))
	}
}

func TestWithExperimentOutputRootMovesOutputsButKeepsMsaHeuristicCache(t *testing.T) {
	atspData := makeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)
	output := withExperimentOutputRoot(atspData, finalResultsDirectoryName)

	if output.msaHeuristicDirectoryPath != atspData.msaHeuristicDirectoryPath {
		t.Fatalf("expected MSA heuristic cache path to stay %s, got %s", atspData.msaHeuristicDirectoryPath, output.msaHeuristicDirectoryPath)
	}

	expectedResultPath := filepath.Join(finalResultsDirectoryName, "test", resultFileName)
	if output.resultFilePath != expectedResultPath {
		t.Fatalf("unexpected final result path\nwant: %s\n got: %s", expectedResultPath, output.resultFilePath)
	}

	if output.optimalUniqueToursCsvPath != atspData.optimalUniqueToursCsvPath {
		t.Fatalf("expected solutions path to stay %s, got %s", atspData.optimalUniqueToursCsvPath, output.optimalUniqueToursCsvPath)
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
		{heuristicMsaHeuristic, finalMsaHeuristicWeight, 0.0},
		{heuristicCycleCover, finalCycleCoverWeight, 0.0},
		{heuristicCycleCoverMsaPatching, finalCycleCoverMsaPatchingWeight, finalCycleCoverMsaPatchingMsaPatchBias},
	}

	if len(configs) != len(expected) {
		t.Fatalf("expected %d final experiment configurations, got %d", len(expected), len(configs))
	}

	for i, config := range configs {
		if config.heuristic != expected[i].heuristic {
			t.Fatalf("unexpected final heuristic at %d: want %s got %s", i, expected[i].heuristic, config.heuristic)
		}
		if len(config.parameters) != 1 {
			t.Fatalf("expected one final parameter set for %s, got %d", config.heuristic, len(config.parameters))
		}

		parameters := config.parameters[0]
		if parameters.alpha != defaultExperimentAlpha ||
			parameters.beta != defaultExperimentBeta ||
			parameters.rho != defaultExperimentRho ||
			parameters.heuristicWeight != expected[i].weight ||
			parameters.msaPatchBias != expected[i].bias {
			t.Fatalf("unexpected final parameters for %s: %+v", config.heuristic, parameters)
		}
	}
}

func TestGenerateParametersUsesMsaPatchBiasOnlyForPatching(t *testing.T) {
	baselineParameters := generateParameters(heuristicBaseline)
	if len(baselineParameters) != 1 {
		t.Fatalf("expected one baseline parameter set, got %d", len(baselineParameters))
	}
	if baselineParameters[0].heuristicWeight != defaultBaselineHeuristicWeight {
		t.Fatalf("unexpected baseline heuristic weight: %+v", baselineParameters[0])
	}

	msaParameters := generateParameters(heuristicMsaHeuristic)
	if len(msaParameters) != 10 {
		t.Fatalf("expected 10 MSA heuristic parameter sets, got %d", len(msaParameters))
	}
	for _, parameters := range msaParameters {
		if parameters.heuristicWeight == 0 {
			t.Fatalf("expected MSA heuristic tuning to skip baseline-equivalent zero weight, got %+v", parameters)
		}
		if parameters.msaPatchBias != 0 {
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
		if parameters.heuristicWeight == 0 {
			zeroWeightCount++
		}
		if parameters.msaPatchBias == 0 {
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

func TestSetDimensionDependantParametersScalesIterationsLinearly(t *testing.T) {
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
			parameters := ExperimentParameters{iterations: 1}
			setDimensionDependantParameters(test.dimension, &parameters)

			if parameters.iterations != test.expectedIterations {
				t.Fatalf("expected %d iterations, got %d", test.expectedIterations, parameters.iterations)
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
	if len(cycleCoverConfigurations) != 1 || cycleCoverConfigurations[0].heuristic != heuristicCycleCover {
		t.Fatalf("expected only cycle-cover configuration, got %+v", cycleCoverConfigurations)
	}

	patchingConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if len(patchingConfigurations) != 1 || patchingConfigurations[0].heuristic != heuristicCycleCoverMsaPatching {
		t.Fatalf("expected only cycle-cover MSA-patching configuration, got %+v", patchingConfigurations)
	}

	randomSparseConfigurations, err := selectFinalExperimentConfigurations(heuristicRandomSparse)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(random-sparse) returned error: %v", err)
	}
	if len(randomSparseConfigurations) != 1 || randomSparseConfigurations[0].heuristic != heuristicRandomSparse || len(randomSparseConfigurations[0].parameters) != len(randomSparseSeeds) {
		t.Fatalf("expected only random-sparse control configuration, got %+v", randomSparseConfigurations)
	}

	shuffledMsaConfigurations, err := selectFinalExperimentConfigurations(heuristicShuffledMsa)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(shuffled-msa) returned error: %v", err)
	}
	if len(shuffledMsaConfigurations) != 1 || shuffledMsaConfigurations[0].heuristic != heuristicShuffledMsa || len(shuffledMsaConfigurations[0].parameters) != len(shuffledMsaSeeds) {
		t.Fatalf("expected only shuffled-MSA control configuration, got %+v", shuffledMsaConfigurations)
	}

	if _, err := selectFinalExperimentConfigurations("unknown"); err == nil {
		t.Fatal("expected invalid final heuristic to be rejected")
	}
}

func TestMsaImpactExperimentConfigurationsSweepMsaHeuristicWeights(t *testing.T) {
	configurations := msaImpactExperimentConfigurations()
	if len(configurations) != 2 {
		t.Fatalf("expected baseline and MSA impact configurations, got %+v", configurations)
	}

	msaConfiguration := configurations[1]
	if msaConfiguration.heuristic != heuristicMsaHeuristic {
		t.Fatalf("expected second configuration to be MSA heuristic, got %+v", msaConfiguration)
	}
	if !msaConfiguration.saveAllParameterRows {
		t.Fatal("MSA impact should save every swept parameter row")
	}
	if len(msaConfiguration.parameters) != len(msaImpactHeuristicWeights) {
		t.Fatalf("expected %d MSA impact weights, got %d", len(msaImpactHeuristicWeights), len(msaConfiguration.parameters))
	}
	for i, expectedWeight := range msaImpactHeuristicWeights {
		if msaConfiguration.parameters[i].heuristicWeight != expectedWeight {
			t.Fatalf("unexpected MSA impact weight at %d: want %.2f got %.2f", i, expectedWeight, msaConfiguration.parameters[i].heuristicWeight)
		}
	}
}

func TestMsaImpactControlExperimentConfigurationsUseSelectedControlWeight(t *testing.T) {
	const selectedWeight = 0.3

	configurations := msaImpactControlExperimentConfigurations(selectedWeight)
	if len(configurations) != 2 {
		t.Fatalf("expected two MSA impact control configurations, got %+v", configurations)
	}

	if configurations[0].heuristic != heuristicDistanceRankedSparse || len(configurations[0].parameters) != 1 {
		t.Fatalf("unexpected distance-ranked MSA impact control configuration: %+v", configurations[0])
	}
	if configurations[0].parameters[0].heuristicWeight != selectedWeight {
		t.Fatalf("distance-ranked control should use selected weight %.2f, got %+v", selectedWeight, configurations[0].parameters[0])
	}

	if configurations[1].heuristic != heuristicShuffledMsa || len(configurations[1].parameters) != len(shuffledMsaSeeds) {
		t.Fatalf("unexpected shuffled-MSA impact control configuration: %+v", configurations[1])
	}
	for i, parameter := range configurations[1].parameters {
		if parameter.heuristicWeight != selectedWeight || parameter.randomSeed != shuffledMsaSeeds[i] {
			t.Fatalf("unexpected shuffled-MSA parameter at %d: %+v", i, parameter)
		}
	}
}

func TestFinalExperimentOutputRootUsesControlsSubdirectoryForSparseControls(t *testing.T) {
	mainConfigurations, err := selectFinalExperimentConfigurations(finalHeuristicAll)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(all) returned error: %v", err)
	}
	if finalExperimentOutputRootForConfigurations(runModeFinal, mainConfigurations) != finalResultsDirectoryName {
		t.Fatalf("main final run should use %s", finalResultsDirectoryName)
	}

	controlConfigurations, err := selectFinalExperimentConfigurations(finalHeuristicControls)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(controls) returned error: %v", err)
	}
	expectedControlsRoot := filepath.Join(filepath.Dir(finalResultsDirectoryName), "controls")
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

	msaHeuristicConfigurations, err := selectFinalExperimentConfigurations(heuristicMsaHeuristic)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(msa-heuristic) returned error: %v", err)
	}
	if !finalConfigurationsUseMsaHeuristic(msaHeuristicConfigurations) {
		t.Fatal("MSA heuristic final run should require MSA heuristic")
	}

	patchingConfigurations, err := selectFinalExperimentConfigurations(heuristicCycleCoverMsaPatching)
	if err != nil {
		t.Fatalf("selectFinalExperimentConfigurations(cycle-cover-msa-patching) returned error: %v", err)
	}
	if !finalConfigurationsUseMsaHeuristic(patchingConfigurations) {
		t.Fatal("cycle-cover MSA-patching final run should require MSA heuristic")
	}
}

func TestSaveFinalHeuristicStatisticsMergesSelectedHeuristicIntoExistingResultCsv(t *testing.T) {
	path := filepath.Join(t.TempDir(), "result.csv")
	existing := []HeuristicExperimentStatistics{
		{
			heuristic:  heuristicBaseline,
			statistics: makeTestExperimentStatistics(0.0, 5.0, 10.0),
		},
		{
			heuristic:  heuristicMsaHeuristic,
			statistics: makeTestExperimentStatistics(finalMsaHeuristicWeight, 3.0, 20.0),
		},
		{
			heuristic:  heuristicCycleCover,
			statistics: makeTestExperimentStatistics(finalCycleCoverWeight, 7.0, 0.0),
		},
		{
			heuristic:  "obsolete",
			statistics: makeTestExperimentStatistics(1.0, 0.5, 100.0),
		},
	}
	if err := saveHeuristicStatistics(path, existing); err != nil {
		t.Fatalf("failed to seed final result CSV: %v", err)
	}

	replacement := []HeuristicExperimentStatistics{
		{
			heuristic:  heuristicCycleCover,
			statistics: makeTestExperimentStatistics(finalCycleCoverWeight, 1.5, 40.0),
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
		byHeuristic[statistic.heuristic] = statistic
	}

	if byHeuristic[heuristicBaseline].statistics.averageBestDeviation != 5.0 {
		t.Fatalf("baseline row was not preserved: %+v", byHeuristic[heuristicBaseline])
	}
	if byHeuristic[heuristicMsaHeuristic].statistics.averageBestDeviation != 3.0 {
		t.Fatalf("MSA heuristic row was not preserved: %+v", byHeuristic[heuristicMsaHeuristic])
	}
	if byHeuristic[heuristicCycleCover].statistics.averageBestDeviation != 1.5 ||
		byHeuristic[heuristicCycleCover].statistics.successRate != 40.0 {
		t.Fatalf("cycle-cover row was not replaced: %+v", byHeuristic[heuristicCycleCover])
	}
	if _, ok := byHeuristic["obsolete"]; ok {
		t.Fatal("obsolete heuristic row should not be preserved in final result CSV")
	}
}

func TestFinalModeAndBaselineHeuristicAreValid(t *testing.T) {
	if !isValidRunMode(runModeFinal) {
		t.Fatal("final run mode should be valid")
	}
	if !isValidRunMode(runModeFinal3Opt) {
		t.Fatal("final+3opt run mode should be valid")
	}
	if !isValidRunMode(runModeRebuildMsa) {
		t.Fatal("rebuild-msa run mode should be valid")
	}
	if !isValidRunMode(runModeMsaImpact) {
		t.Fatal("MSA impact run mode should be valid")
	}
	if !isValidHeuristic(heuristicBaseline) {
		t.Fatal("baseline heuristic should be valid")
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
	if finalExperimentOutputRoot(runModeFinal3Opt) != finalThreeOptResultsDirectoryName {
		t.Fatalf("final+3opt should use %s", finalThreeOptResultsDirectoryName)
	}
	if finalExperimentOutputRoot(runModeFinal) != finalResultsDirectoryName {
		t.Fatalf("final should use %s", finalResultsDirectoryName)
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
	if shouldRunAnalysisAfterFinalExperiments(runModeFinal, heuristicMsaHeuristic) {
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
	originalResultsDirectoryName := resultsDirectoryName
	resultsDirectoryName = t.TempDir()
	t.Cleanup(func() {
		resultsDirectoryName = originalResultsDirectoryName
	})

	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsDirectoryName)
	saveStatistics(resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching), heuristicCycleCoverMsaPatching, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.6, 2.0, 20.0),
	})

	if err := runAnalysisMode([]AtspData{atspData}, analysisScopeTuning, []string{heuristicCycleCoverMsaPatching}, 1); err != nil {
		t.Fatalf("runAnalysisMode(tuning) returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(filepath.Join(resultsDirectoryName, "tuning_summary.md"))
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

func TestRunBoundedInstanceJobsRespectsWorkerLimit(t *testing.T) {
	atspsData := []AtspData{
		makeAtspDataInResultsDirectory("a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		makeAtspDataInResultsDirectory("b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		makeAtspDataInResultsDirectory("c.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		makeAtspDataInResultsDirectory("d.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		makeAtspDataInResultsDirectory("e.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
	}

	var activeWorkers int32
	var maxActiveWorkers int32
	var processedJobs int32

	err := runBoundedInstanceJobs(atspsData, 2, func(atspData AtspData) error {
		currentActiveWorkers := atomic.AddInt32(&activeWorkers, 1)
		for {
			currentMax := atomic.LoadInt32(&maxActiveWorkers)
			if currentActiveWorkers <= currentMax || atomic.CompareAndSwapInt32(&maxActiveWorkers, currentMax, currentActiveWorkers) {
				break
			}
		}

		time.Sleep(10 * time.Millisecond)
		atomic.AddInt32(&processedJobs, 1)
		atomic.AddInt32(&activeWorkers, -1)
		return nil
	})
	if err != nil {
		t.Fatalf("runBoundedInstanceJobs returned unexpected error: %v", err)
	}

	if processedJobs != int32(len(atspsData)) {
		t.Fatalf("expected %d processed jobs, got %d", len(atspsData), processedJobs)
	}
	if maxActiveWorkers > 2 {
		t.Fatalf("expected at most two active workers, got %d", maxActiveWorkers)
	}
}

func TestRunBoundedIndexJobsRespectsWorkerLimitAndKeepsIndexedResults(t *testing.T) {
	result := make([]int, 6)
	var activeWorkers int32
	var maxActiveWorkers int32

	err := runBoundedIndexJobs(len(result), 3, func(index int) error {
		currentActiveWorkers := atomic.AddInt32(&activeWorkers, 1)
		for {
			currentMax := atomic.LoadInt32(&maxActiveWorkers)
			if currentActiveWorkers <= currentMax || atomic.CompareAndSwapInt32(&maxActiveWorkers, currentMax, currentActiveWorkers) {
				break
			}
		}

		time.Sleep(10 * time.Millisecond)
		result[index] = index + 100
		atomic.AddInt32(&activeWorkers, -1)
		return nil
	})
	if err != nil {
		t.Fatalf("runBoundedIndexJobs returned unexpected error: %v", err)
	}

	if maxActiveWorkers > 3 {
		t.Fatalf("expected at most three active workers, got %d", maxActiveWorkers)
	}
	for index, value := range result {
		if value != index+100 {
			t.Fatalf("expected indexed result %d to be %d, got %d", index, index+100, value)
		}
	}
}

func TestRunBoundedIndexJobsReturnsFirstError(t *testing.T) {
	err := runBoundedIndexJobs(3, 1, func(index int) error {
		if index == 1 {
			return os.ErrInvalid
		}
		return nil
	})
	if err != os.ErrInvalid {
		t.Fatalf("expected %v, got %v", os.ErrInvalid, err)
	}
}

func TestRunFinalExperimentForInstanceWithParameterWorkersRejectsInvalidWorkerCount(t *testing.T) {
	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir())
	configurations := []finalExperimentConfiguration{
		{
			heuristic: heuristicBaseline,
			parameters: []ExperimentParameters{
				newDefaultExperimentParameters(defaultBaselineHeuristicWeight),
			},
		},
	}

	err := runFinalExperimentForInstanceWithParameterWorkers(atspData, t.TempDir(), false, configurations, 1, 0)
	if err == nil || !strings.Contains(err.Error(), "parameter workers must be at least one") {
		t.Fatalf("expected invalid parameter worker error, got %v", err)
	}
}

func TestRunBoundedInstanceJobsReturnsFirstError(t *testing.T) {
	atspsData := []AtspData{
		makeAtspDataInResultsDirectory("ok.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
		makeAtspDataInResultsDirectory("bad.atsp", [][]float64{{0, 1}, {1, 0}}, 2, t.TempDir()),
	}

	err := runBoundedInstanceJobs(atspsData, 1, func(atspData AtspData) error {
		if atspData.name == "bad" {
			return os.ErrInvalid
		}
		return nil
	})
	if err != os.ErrInvalid {
		t.Fatalf("expected %v, got %v", os.ErrInvalid, err)
	}
}

func TestRunBoundedInstanceJobsHandlesEmptyInput(t *testing.T) {
	called := false
	err := runBoundedInstanceJobs(nil, 2, func(atspData AtspData) error {
		called = true
		return nil
	})
	if err != nil {
		t.Fatalf("runBoundedInstanceJobs returned unexpected error: %v", err)
	}
	if called {
		t.Fatal("job should not be called for empty input")
	}
}

func TestRunRebuildMsaModeRemovesStaleArtifactsAndRecreatesFiles(t *testing.T) {
	root := t.TempDir()
	matrix := [][]float64{
		{0, 1, 5, 2},
		{4, 0, 1, 3},
		{2, 4, 0, 1},
		{1, 2, 4, 0},
	}
	atspData := makeAtspDataInResultsDirectory("sample.atsp", matrix, 0, filepath.Join(root, "results"))
	atspData.msaHeuristicDirectoryPath = filepath.Join(root, "msa", "sample")
	atspData.msaHeuristicHeatmapPlotPath = filepath.Join(atspData.msaHeuristicDirectoryPath, "plots", "msa_heuristic_heatmap.png")
	atspData.msaHeuristicHistogramPlotPath = filepath.Join(atspData.msaHeuristicDirectoryPath, "plots", "msa_heuristic_histogram.png")
	atspData.msaThinnessScoresCsvPath = filepath.Join(atspData.msaHeuristicDirectoryPath, "msa_thinness_scores.csv")
	atspData.msaThinnessHistogramPlotPath = filepath.Join(atspData.msaHeuristicDirectoryPath, "plots", "msa_thinness_histogram.png")

	stalePath := filepath.Join(atspData.msaHeuristicDirectoryPath, "stale.txt")
	if err := os.MkdirAll(atspData.msaHeuristicDirectoryPath, 0700); err != nil {
		t.Fatalf("failed to create stale MSA directory: %v", err)
	}
	if err := os.WriteFile(stalePath, []byte("stale"), 0644); err != nil {
		t.Fatalf("failed to write stale file: %v", err)
	}

	if err := runRebuildMsaMode([]AtspData{atspData}, 1); err != nil {
		t.Fatalf("runRebuildMsaMode returned unexpected error: %v", err)
	}

	if _, err := os.Stat(stalePath); !os.IsNotExist(err) {
		t.Fatalf("expected stale file to be removed, stat error: %v", err)
	}
	assertPathExists(t, filepath.Join(atspData.msaHeuristicDirectoryPath, "msa_heuristic.csv"))
	assertPathExists(t, atspData.msaHeuristicHeatmapPlotPath)
	assertPathExists(t, atspData.msaHeuristicHistogramPlotPath)
	assertPathExists(t, atspData.msaThinnessScoresCsvPath)
	assertPathExists(t, atspData.msaThinnessHistogramPlotPath)

	thinnessCsv, err := os.ReadFile(atspData.msaThinnessScoresCsvPath)
	if err != nil {
		t.Fatalf("failed to read MSA thinness scores CSV: %v", err)
	}
	assertContains(t, string(thinnessCsv), "Root,Thinness score (branch surplus),Max outgoing degree,Branching vertices,Total cost")

	rootMsaFiles, err := filepath.Glob(filepath.Join(atspData.msaHeuristicDirectoryPath, "msas", "*.csv"))
	if err != nil {
		t.Fatalf("failed to glob root MSA files: %v", err)
	}
	if len(rootMsaFiles) != len(matrix) {
		t.Fatalf("expected %d root MSA files, got %d", len(matrix), len(rootMsaFiles))
	}
}

func TestStartCPUProfileIsOptional(t *testing.T) {
	stopProfiling, err := startCPUProfile("")
	if err != nil {
		t.Fatalf("startCPUProfile returned unexpected error: %v", err)
	}
	stopProfiling()
}

func TestStartCPUProfileWritesRequestedPath(t *testing.T) {
	profilePath := filepath.Join(t.TempDir(), "profiles", "cpu.prof")
	stopProfiling, err := startCPUProfile(profilePath)
	if err != nil {
		t.Fatalf("startCPUProfile returned unexpected error: %v", err)
	}
	for i := 0; i < 10000; i++ {
		_ = i * i
	}
	stopProfiling()

	if _, err := os.Stat(profilePath); err != nil {
		t.Fatalf("expected CPU profile to be written: %v", err)
	}
}

func TestSelectAtspFilesTuning(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := selectAtspFiles(paths, instanceSetTuning)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, tuningInstanceFiles) {
		t.Fatalf("tuning selection mismatch\nwant: %v\n got: %v", tuningInstanceFiles, selectedFiles)
	}
}

func TestSelectAtspFilesSmoke(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := selectAtspFiles(paths, instanceSetSmoke)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, smokeInstanceFiles) {
		t.Fatalf("smoke selection mismatch\nwant: %v\n got: %v", smokeInstanceFiles, selectedFiles)
	}
}

func TestSelectAtspFilesMsaImpact(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := selectAtspFiles(paths, instanceSetMsaImpact)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, msaImpactInstanceFiles) {
		t.Fatalf("MSA impact selection mismatch\nwant: %v\n got: %v", msaImpactInstanceFiles, selectedFiles)
	}
}

func TestSelectAtspFilesEvaluation(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := selectAtspFiles(paths, instanceSetEvaluation)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, evaluationInstanceFiles) {
		t.Fatalf("evaluation selection mismatch\nwant: %v\n got: %v", evaluationInstanceFiles, selectedFiles)
	}
}

func TestSelectedAtspFilesHaveKnownOptima(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	for _, instanceSet := range []string{instanceSetSmoke, instanceSetMsaImpact, instanceSetTuning, instanceSetEvaluation, instanceSetAllKnown} {
		selected, err := selectAtspFiles(paths, instanceSet)
		if err != nil {
			t.Fatalf("selectAtspFiles(%s) returned unexpected error: %v", instanceSet, err)
		}

		for _, selectedPath := range selected {
			name, _, knownOptimal, err := parsing.ParseTSPLIBFile(selectedPath)
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
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	tuning, err := selectAtspFiles(paths, instanceSetTuning)
	if err != nil {
		t.Fatalf("selectAtspFiles(tuning) returned unexpected error: %v", err)
	}
	evaluation, err := selectAtspFiles(paths, instanceSetEvaluation)
	if err != nil {
		t.Fatalf("selectAtspFiles(evaluation) returned unexpected error: %v", err)
	}
	allKnown, err := selectAtspFiles(paths, instanceSetAllKnown)
	if err != nil {
		t.Fatalf("selectAtspFiles(all-known) returned unexpected error: %v", err)
	}
	smoke, err := selectAtspFiles(paths, instanceSetSmoke)
	if err != nil {
		t.Fatalf("selectAtspFiles(smoke) returned unexpected error: %v", err)
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

func TestSelectedInstanceSetForMode(t *testing.T) {
	if selected := selectedInstanceSetForMode(runModeFinal, instanceSetTuning, false); selected != instanceSetEvaluation {
		t.Fatalf("final mode without explicit instances should default to %s, got %s", instanceSetEvaluation, selected)
	}
	if selected := selectedInstanceSetForMode(runModeFinal3Opt, instanceSetTuning, false); selected != instanceSetEvaluation {
		t.Fatalf("final+3opt mode without explicit instances should default to %s, got %s", instanceSetEvaluation, selected)
	}
	if selected := selectedInstanceSetForMode(runModeMsaImpact, instanceSetTuning, false); selected != instanceSetMsaImpact {
		t.Fatalf("MSA impact mode without explicit instances should default to %s, got %s", instanceSetMsaImpact, selected)
	}
	if selected := selectedInstanceSetForMode(runModeRebuildMsa, instanceSetTuning, false); selected != instanceSetAllKnown {
		t.Fatalf("rebuild-msa mode without explicit instances should default to %s, got %s", instanceSetAllKnown, selected)
	}
	if selected := selectedInstanceSetForMode(runModeFinal, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("final mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeMsaImpact, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("MSA impact mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeRebuildMsa, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("rebuild-msa mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeExperiment, instanceSetTuning, false); selected != instanceSetTuning {
		t.Fatalf("experiment mode should keep requested instances, got %s", selected)
	}
}

func sampleStructuralAnalyses() []structuralComparison.InstanceAnalysis {
	return []structuralComparison.InstanceAnalysis{
		{
			Instance:  "b",
			Dimension: 5,
			Metrics: structuralComparison.InstanceMetrics{
				FoundOptimalTourCount:       2,
				UniqueFoundOptimalEdgeCount: 10,
				HighMsaHeuristicMetrics: structuralComparison.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 4,
					Precision:        0.8,
					Recall:           0.4,
				},
				CycleCoverMetrics: structuralComparison.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 3,
					Precision:        0.6,
					Recall:           0.3,
				},
				CycleCoverMsaPatchingMetrics: structuralComparison.EdgeSetMetrics{
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
			Metrics: structuralComparison.InstanceMetrics{
				FoundOptimalTourCount:       1,
				UniqueFoundOptimalEdgeCount: 4,
				HighMsaHeuristicMetrics: structuralComparison.EdgeSetMetrics{
					EdgeCount:        2,
					OptimalEdgeCount: 1,
					Precision:        0.5,
					Recall:           0.25,
				},
				CycleCoverMetrics: structuralComparison.EdgeSetMetrics{
					EdgeCount:        4,
					OptimalEdgeCount: 3,
					Precision:        0.75,
					Recall:           0.75,
				},
				CycleCoverMsaPatchingMetrics: structuralComparison.EdgeSetMetrics{
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
					averageMinDeviation:  4.0,
					successRate:          10.0,
					averageBestIteration: 80.0,
					iterations:           100,
				},
				heuristicMsaHeuristic: {
					averageMinDeviation:  2.5,
					successRate:          20.0,
					averageBestIteration: 60.0,
					iterations:           100,
				},
				heuristicCycleCover: {
					averageMinDeviation:  2.0,
					successRate:          30.0,
					averageBestIteration: 40.0,
					iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					averageMinDeviation:  1.8,
					successRate:          25.0,
					averageBestIteration: 35.0,
					iterations:           100,
				},
			},
		},
		{
			instance: "b",
			metrics: map[string]finalResultsSummaryMetric{
				heuristicBaseline: {
					averageMinDeviation:  5.0,
					successRate:          0.0,
					averageBestIteration: 70.0,
					iterations:           100,
				},
				heuristicMsaHeuristic: {
					averageMinDeviation:  3.5,
					successRate:          20.0,
					averageBestIteration: 40.0,
					iterations:           100,
				},
				heuristicCycleCover: {
					averageMinDeviation:  3.5,
					successRate:          20.0,
					averageBestIteration: 30.0,
					iterations:           100,
				},
				heuristicCycleCoverMsaPatching: {
					averageMinDeviation:  3.0,
					successRate:          25.0,
					averageBestIteration: 20.0,
					iterations:           100,
				},
			},
		},
	}
}

func makeTestExperimentStatistics(heuristicWeight, averageBestDeviation, successRate float64) ExperimentsDataStatistics {
	return ExperimentsDataStatistics{
		ExperimentParameters: ExperimentParameters{
			alpha:           defaultExperimentAlpha,
			beta:            defaultExperimentBeta,
			rho:             defaultExperimentRho,
			heuristicWeight: heuristicWeight,
			iterations:      500,
		},
		minBestAtIteration:               1,
		averageBestAtIteration:           2.0,
		maxBestAtIteration:               3,
		minThreeOptImprovementsCount:     0,
		averageThreeOptImprovementsCount: 0.0,
		maxThreeOptImprovementsCount:     0,
		minBestDeviation:                 averageBestDeviation,
		averageBestDeviation:             averageBestDeviation,
		maxBestDeviation:                 averageBestDeviation,
		successRate:                      successRate,
	}
}

func makeTestRandomSparseExperimentStatistics(randomSeed int64, averageBestDeviation, successRate float64) ExperimentsDataStatistics {
	statistics := makeTestExperimentStatistics(finalMsaHeuristicWeight, averageBestDeviation, successRate)
	statistics.randomSeed = randomSeed
	return statistics
}

func writeTestControlStatisticsForWeightSummary(t *testing.T, atspData AtspData, heuristic string) {
	t.Helper()

	firstInstance := strings.Contains(atspData.name, "sample-a")
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
	statistics.randomSeed = randomSeed
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
