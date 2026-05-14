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

	if statistics[0].pCmsa != 0.50 || statistics[0].successRate != 40.00 {
		t.Fatalf("unexpected statistic values: pCmsa=%f successRate=%f", statistics[0].pCmsa, statistics[0].successRate)
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

func TestBuildCmsaHeuristicModifiersBoostsOnlyEdgesUsedByEveryMsa(t *testing.T) {
	cmsa := [][]float64{
		{0, 3, 2, 0},
		{0, 0, 3, 1},
		{1, 0, 0, 2},
		{0, 0, 0, 0},
	}

	modifiers := buildCmsaHeuristicModifiers(cmsa, 0.5)
	expected := [][]float64{
		{1, 1.5, 1, 1},
		{1, 1, 1.5, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected CMSA modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCmsaHeuristicModifiersReturnsNeutralMatrixWhenStrengthIsZero(t *testing.T) {
	cmsa := [][]float64{
		{0, 3},
		{3, 0},
	}

	modifiers := buildCmsaHeuristicModifiers(cmsa, 0)
	expected := [][]float64{
		{1, 1},
		{1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected neutral modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCmsaCycleCoverMembershipHeuristicModifiersSplitsOverlapAndDifference(t *testing.T) {
	cmsa := [][]float64{
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

	overlapModifiers := buildCmsaCycleCoverMembershipHeuristicModifiers(cmsa, cycleCover, 0.5, true)
	expectedOverlap := [][]float64{
		{1, 1.5, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1.5},
		{1, 1, 1, 1},
	}

	if !reflect.DeepEqual(overlapModifiers, expectedOverlap) {
		t.Fatalf("unexpected CMSA-overlap modifiers\nwant: %v\n got: %v", expectedOverlap, overlapModifiers)
	}

	differenceModifiers := buildCmsaCycleCoverMembershipHeuristicModifiers(cmsa, cycleCover, 0.5, false)
	expectedDifference := [][]float64{
		{1, 1, 1.5, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1.5, 1, 1, 1},
	}

	if !reflect.DeepEqual(differenceModifiers, expectedDifference) {
		t.Fatalf("unexpected CMSA-difference modifiers\nwant: %v\n got: %v", expectedDifference, differenceModifiers)
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

func TestHeuristicSpecificPathsKeepCmsaBaselinePaths(t *testing.T) {
	atspData := makeAtspData("test.atsp", [][]float64{{0, 1}, {1, 0}}, 2)

	if resultFilePathForHeuristic(atspData, heuristicCmsa) != filepath.Join(resultsDirectoryName, "test", resultFileName) {
		t.Fatalf("CMSA result path should keep the existing baseline location")
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCover) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover.csv") {
		t.Fatalf("unexpected cycle-cover result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCover))
	}

	if resultFilePathForHeuristic(atspData, heuristicCmsaOverlap) != filepath.Join(resultsDirectoryName, "test", "result_cmsa_overlap.csv") {
		t.Fatalf("unexpected CMSA-overlap result path: %s", resultFilePathForHeuristic(atspData, heuristicCmsaOverlap))
	}

	if resultFilePathForHeuristic(atspData, heuristicCmsaDifference) != filepath.Join(resultsDirectoryName, "test", "result_cmsa_difference.csv") {
		t.Fatalf("unexpected CMSA-difference result path: %s", resultFilePathForHeuristic(atspData, heuristicCmsaDifference))
	}
}

func TestCmsaOverlapAndDifferenceHeuristicsUseCycleCover(t *testing.T) {
	if !isValidHeuristic(heuristicCmsaOverlap) || !isValidHeuristic(heuristicCmsaDifference) {
		t.Fatal("CMSA-overlap and CMSA-difference should be valid heuristics")
	}

	if !heuristicUsesCycleCover(heuristicCmsaOverlap) || !heuristicUsesCycleCover(heuristicCmsaDifference) {
		t.Fatal("CMSA-overlap and CMSA-difference should require cycle cover")
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

func TestSelectedAtspFilesHaveKnownOptima(t *testing.T) {
	paths, err := filepath.Glob(filepath.Join("tsplib_files", "*.atsp"))
	if err != nil {
		t.Fatalf("failed to glob ATSP files: %v", err)
	}

	for _, instanceSet := range []string{instanceSetSmoke, instanceSetBalanced, instanceSetAllKnown} {
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
