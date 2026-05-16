package main

import (
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

	if statistics[0].heuristicWeight != 0.50 || statistics[0].successRate != 40.00 {
		t.Fatalf("unexpected statistic values: heuristicWeight=%f successRate=%f", statistics[0].heuristicWeight, statistics[0].successRate)
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
			heuristic: heuristicMsaSupport,
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
	secondRows := []HeuristicExperimentStatistics{
		{
			heuristic: heuristicBaseline,
			statistics: ExperimentsDataStatistics{
				averageBestDeviation: 6.75,
				successRate:          20.0,
			},
		},
		{
			heuristic: heuristicMsaSupport,
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
		"- **Best-or-tied average best deviation counts: Baseline 0/2, MSA support 0/2, Cycle cover 2/2.**",
		"",
		"<table>",
		"<thead>",
		"<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">Baseline</th><th colspan=\"2\">MSA support</th><th colspan=\"2\">Cycle cover</th></tr>",
		"<tr><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th></tr>",
		"</thead>",
		"<tbody>",
		"<tr><td>sample-a</td><td align=\"right\">4.25</td><td align=\"right\">10.00</td><td align=\"right\">2.50</td><td align=\"right\">20.00</td><td align=\"right\"><strong>1.75</strong></td><td align=\"right\"><strong>30.00</strong></td></tr>",
		"<tr><td>sample-b</td><td align=\"right\">6.75</td><td align=\"right\">20.00</td><td align=\"right\">3.50</td><td align=\"right\">40.00</td><td align=\"right\"><strong>2.25</strong></td><td align=\"right\"><strong>60.00</strong></td></tr>",
		"<tr><td><strong>Average</strong></td><td align=\"right\">5.50</td><td align=\"right\">15.00</td><td align=\"right\">3.00</td><td align=\"right\">30.00</td><td align=\"right\"><strong>2.00</strong></td><td align=\"right\"><strong>45.00</strong></td></tr>",
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
			heuristic: heuristicMsaSupport,
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

	summaryPath, saved, err := runFinalResultsAnalysis([]AtspData{atspData})
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

func TestBuildMsaSupportHeuristicModifiersBoostsOnlyEdgesUsedByEveryMsa(t *testing.T) {
	msaSupport := [][]float64{
		{0, 3, 2, 0},
		{0, 0, 3, 1},
		{1, 0, 0, 2},
		{0, 0, 0, 0},
	}

	modifiers := buildMsaSupportHeuristicModifiers(msaSupport, 0.5)
	expected := [][]float64{
		{1, 1.5, 1, 1},
		{1, 1, 1.5, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected MSA support modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildMsaSupportHeuristicModifiersReturnsNeutralMatrixWhenStrengthIsZero(t *testing.T) {
	msaSupport := [][]float64{
		{0, 3},
		{3, 0},
	}

	modifiers := buildMsaSupportHeuristicModifiers(msaSupport, 0)
	expected := [][]float64{
		{1, 1},
		{1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected neutral modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildMsaSupportCycleCoverMembershipHeuristicModifiersSplitsOverlapAndDifference(t *testing.T) {
	msaSupport := [][]float64{
		{0, 3, 3, 0},
		{2, 0, 0, 0},
		{0, 0, 0, 3},
		{3, 0, 2, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 0},
	}

	overlapModifiers := buildMsaSupportCycleCoverMembershipHeuristicModifiers(msaSupport, cycleCover, 0.5, true)
	expectedOverlap := [][]float64{
		{1, 1.5, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1.5},
		{1, 1, 1, 1},
	}

	if !reflect.DeepEqual(overlapModifiers, expectedOverlap) {
		t.Fatalf("unexpected MSA support-overlap modifiers\nwant: %v\n got: %v", expectedOverlap, overlapModifiers)
	}

	differenceModifiers := buildMsaSupportCycleCoverMembershipHeuristicModifiers(msaSupport, cycleCover, 0.5, false)
	expectedDifference := [][]float64{
		{1, 1, 1.5, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1.5, 1, 1, 1},
	}

	if !reflect.DeepEqual(differenceModifiers, expectedDifference) {
		t.Fatalf("unexpected MSA support-difference modifiers\nwant: %v\n got: %v", expectedDifference, differenceModifiers)
	}
}

func TestBuildCycleCoverHeuristicModifiersBoostsOnlyCycleCoverEdges(t *testing.T) {
	cycleCover := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}

	modifiers := buildCycleCoverHeuristicModifiers(cycleCover, 0.4)
	expected := [][]float64{
		{1, 1.4, 1},
		{1, 1, 1.4},
		{1.4, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected cycle-cover modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildHeuristicModifiersReturnsNeutralMatrixForBaseline(t *testing.T) {
	msaSupport := [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}

	modifiers := buildHeuristicModifiers(heuristicBaseline, nil, msaSupport, nil, 1.0)
	expected := [][]float64{
		{1, 1, 1},
		{1, 1, 1},
		{1, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected baseline modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestAnalyzeHeuristicBoostsSummarizesModifierMatrix(t *testing.T) {
	modifiers := [][]float64{
		{1, 2, 1.5},
		{1, 1, 1},
		{1.25, 1, 1},
	}

	row := analyzeHeuristicBoosts("test", "example", modifiers)

	if row.dimension != 3 {
		t.Fatalf("expected dimension 3, got %d", row.dimension)
	}
	if row.boostableEdges != 6 {
		t.Fatalf("expected 6 boostable edges, got %d", row.boostableEdges)
	}
	if row.tourEdgeTarget != 2 {
		t.Fatalf("expected tour edge target 2, got %d", row.tourEdgeTarget)
	}
	if row.boostedEdges != 3 {
		t.Fatalf("expected 3 boosted edges, got %d", row.boostedEdges)
	}
	if row.boostedEdgeDensity != 0.5 {
		t.Fatalf("expected density 0.5, got %f", row.boostedEdgeDensity)
	}
	if row.boostedToTourTarget != 1.5 {
		t.Fatalf("expected boosted-to-tour-target ratio 1.5, got %f", row.boostedToTourTarget)
	}
}

func TestAnalyzeHeuristicBoostsHandlesNeutralMatrix(t *testing.T) {
	modifiers := [][]float64{
		{1, 1},
		{1, 1},
	}

	row := analyzeHeuristicBoosts("test", "neutral", modifiers)

	if row.boostedEdges != 0 {
		t.Fatalf("expected no boosted edges, got %d", row.boostedEdges)
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

func TestHeuristicSpecificPathsKeepMsaSupportBaselinePaths(t *testing.T) {
	atspData := makeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)

	if resultFilePathForHeuristic(atspData, heuristicBaseline) != filepath.Join(resultsDirectoryName, "test", "result_baseline.csv") {
		t.Fatalf("unexpected baseline result path: %s", resultFilePathForHeuristic(atspData, heuristicBaseline))
	}

	if resultFilePathForHeuristic(atspData, heuristicMsaSupport) != filepath.Join(resultsDirectoryName, "test", resultFileName) {
		t.Fatalf("MSA support result path should keep the existing baseline location")
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCover) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover.csv") {
		t.Fatalf("unexpected cycle-cover result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCover))
	}

	if resultFilePathForHeuristic(atspData, heuristicMsaSupportOverlap) != filepath.Join(resultsDirectoryName, "test", "result_msa_support_overlap.csv") {
		t.Fatalf("unexpected MSA support-overlap result path: %s", resultFilePathForHeuristic(atspData, heuristicMsaSupportOverlap))
	}

	if resultFilePathForHeuristic(atspData, heuristicMsaSupportDifference) != filepath.Join(resultsDirectoryName, "test", "result_msa_support_difference.csv") {
		t.Fatalf("unexpected MSA support-difference result path: %s", resultFilePathForHeuristic(atspData, heuristicMsaSupportDifference))
	}
}

func TestWithExperimentOutputRootMovesOutputsButKeepsMsaSupportCache(t *testing.T) {
	atspData := makeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)
	output := withExperimentOutputRoot(atspData, finalResultsDirectoryName)

	if output.msaSupportDirectoryPath != atspData.msaSupportDirectoryPath {
		t.Fatalf("expected MSA support cache path to stay %s, got %s", atspData.msaSupportDirectoryPath, output.msaSupportDirectoryPath)
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
	}{
		{heuristicBaseline, defaultBaselineHeuristicWeight},
		{heuristicMsaSupport, finalMsaSupportWeight},
		{heuristicCycleCover, finalCycleCoverWeight},
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
			parameters.heuristicWeight != expected[i].weight {
			t.Fatalf("unexpected final parameters for %s: %+v", config.heuristic, parameters)
		}
	}
}

func TestFinalModeAndBaselineHeuristicAreValid(t *testing.T) {
	if !isValidRunMode(runModeFinal) {
		t.Fatal("final run mode should be valid")
	}
	if !isValidHeuristic(heuristicBaseline) {
		t.Fatal("baseline heuristic should be valid")
	}
}

func TestMsaSupportOverlapAndDifferenceHeuristicsUseCycleCover(t *testing.T) {
	if !isValidHeuristic(heuristicMsaSupportOverlap) || !isValidHeuristic(heuristicMsaSupportDifference) {
		t.Fatal("MSA support-overlap and MSA support-difference should be valid heuristics")
	}

	if !heuristicUsesCycleCover(heuristicMsaSupportOverlap) || !heuristicUsesCycleCover(heuristicMsaSupportDifference) {
		t.Fatal("MSA support-overlap and MSA support-difference should require cycle cover")
	}
}

func TestSelectAtspFilesSmoke(t *testing.T) {
	paths := []string{
		filepath.Join("tsplib_files", "ft53.atsp"),
	}
	for _, smokeInstanceFile := range smokeInstanceFiles {
		paths = append(paths, filepath.Join("tsplib_files", smokeInstanceFile))
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

func TestSelectAtspFilesBalanced(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := selectAtspFiles(paths, instanceSetBalanced)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, balancedInstanceFiles) {
		t.Fatalf("balanced selection mismatch\nwant: %v\n got: %v", balancedInstanceFiles, selectedFiles)
	}
}

func TestSelectAtspFilesTiny(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := selectAtspFiles(paths, instanceSetTiny)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, tinyInstanceFiles) {
		t.Fatalf("tiny selection mismatch\nwant: %v\n got: %v", tinyInstanceFiles, selectedFiles)
	}
}

func TestSelectAtspFilesLarge(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	selected, err := selectAtspFiles(paths, instanceSetLarge)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	selectedFiles := make([]string, len(selected))
	for i, selectedPath := range selected {
		selectedFiles[i] = filepath.Base(selectedPath)
	}

	if !reflect.DeepEqual(selectedFiles, largeInstanceFiles) {
		t.Fatalf("large selection mismatch\nwant: %v\n got: %v", largeInstanceFiles, selectedFiles)
	}
}

func TestSelectedAtspFilesHaveKnownOptima(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	for _, instanceSet := range []string{instanceSetSmoke, instanceSetTiny, instanceSetBalanced, instanceSetLarge, instanceSetAllKnown} {
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
