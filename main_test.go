package main

import (
	"atsp_aco_msa/modules/analysis/cycleCover"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/parsing"
	"os"
	"path/filepath"
	"reflect"
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
	if err := saveStructuralSimilarityReport(path, sampleCycleCoverAnalyses()); err != nil {
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
	if err := saveMsaHeuristicCycleCoverOverlapReport(path, sampleCycleCoverAnalyses()); err != nil {
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
	if err := saveStructuralPerformanceLinkReport(path, sampleFinalSummaryRows(), sampleCycleCoverAnalyses()); err != nil {
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
	msaParameters := generateParameters(heuristicMsaHeuristic)
	if len(msaParameters) != 11 {
		t.Fatalf("expected 11 MSA heuristic parameter sets, got %d", len(msaParameters))
	}
	for _, parameters := range msaParameters {
		if parameters.msaPatchBias != 0 {
			t.Fatalf("expected MSA patch bias to be zero for non-patching heuristic, got %+v", parameters)
		}
	}

	patchingParameters := generateParameters(heuristicCycleCoverMsaPatching)
	if len(patchingParameters) != 51 {
		t.Fatalf("expected 51 patching parameter sets, got %d", len(patchingParameters))
	}

	zeroWeightCount := 0
	for _, parameters := range patchingParameters {
		if parameters.heuristicWeight == 0 {
			zeroWeightCount++
			if parameters.msaPatchBias != 0 {
				t.Fatalf("expected zero-weight patching parameter to use zero MSA patch bias, got %+v", parameters)
			}
		}
	}
	if zeroWeightCount != 1 {
		t.Fatalf("expected one zero-weight patching parameter set, got %d", zeroWeightCount)
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

	if _, err := selectFinalExperimentConfigurations("unknown"); err == nil {
		t.Fatal("expected invalid final heuristic to be rejected")
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
	expectedControlsRoot := filepath.Join(finalResultsDirectoryName, "controls")
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
	if !isValidHeuristic(heuristicBaseline) {
		t.Fatal("baseline heuristic should be valid")
	}
	if isValidHeuristic(heuristicRandomSparse) {
		t.Fatal("random sparse control should not be valid in normal experiment mode")
	}
	if isValidHeuristic(heuristicDistanceRankedSparse) {
		t.Fatal("distance-ranked sparse control should not be valid in normal experiment mode")
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

	for _, instanceSet := range []string{instanceSetTuning, instanceSetEvaluation, instanceSetAllKnown} {
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

	if len(partition) != len(allKnown) {
		t.Fatalf("tuning/evaluation partition has %d files, all-known has %d", len(partition), len(allKnown))
	}
	for _, path := range allKnown {
		fileName := filepath.Base(path)
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
	if selected := selectedInstanceSetForMode(runModeFinal, instanceSetTuning, true); selected != instanceSetTuning {
		t.Fatalf("final mode should respect explicit instances, got %s", selected)
	}
	if selected := selectedInstanceSetForMode(runModeExperiment, instanceSetTuning, false); selected != instanceSetTuning {
		t.Fatalf("experiment mode should keep requested instances, got %s", selected)
	}
}

func sampleCycleCoverAnalyses() []cycleCover.InstanceAnalysis {
	return []cycleCover.InstanceAnalysis{
		{
			Instance:  "b",
			Dimension: 5,
			Metrics: cycleCover.InstanceMetrics{
				FoundOptimalTourCount:       2,
				UniqueFoundOptimalEdgeCount: 10,
				HighMsaHeuristicMetrics: cycleCover.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 4,
					Precision:        0.8,
					Recall:           0.4,
				},
				CycleCoverMetrics: cycleCover.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 3,
					Precision:        0.6,
					Recall:           0.3,
				},
				CycleCoverMsaPatchingMetrics: cycleCover.EdgeSetMetrics{
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
			Metrics: cycleCover.InstanceMetrics{
				FoundOptimalTourCount:       1,
				UniqueFoundOptimalEdgeCount: 4,
				HighMsaHeuristicMetrics: cycleCover.EdgeSetMetrics{
					EdgeCount:        2,
					OptimalEdgeCount: 1,
					Precision:        0.5,
					Recall:           0.25,
				},
				CycleCoverMetrics: cycleCover.EdgeSetMetrics{
					EdgeCount:        4,
					OptimalEdgeCount: 3,
					Precision:        0.75,
					Recall:           0.75,
				},
				CycleCoverMsaPatchingMetrics: cycleCover.EdgeSetMetrics{
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

func assertContains(t *testing.T, content, expected string) {
	t.Helper()
	if !strings.Contains(content, expected) {
		t.Fatalf("expected content to contain:\n%s\n\ncontent:\n%s", expected, content)
	}
}
