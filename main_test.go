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

func TestBuildCmsaHeuristicModifiersPreservesCurrentCmsaRule(t *testing.T) {
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

func TestBuildCycleCoverComponentIds(t *testing.T) {
	cycleCover := [][]float64{
		{0, 1, 0, 0, 0},
		{0, 0, 1, 0, 0},
		{1, 0, 0, 0, 0},
		{0, 0, 0, 0, 1},
		{0, 0, 0, 1, 0},
	}

	componentIds := buildCycleCoverComponentIds(cycleCover)
	expected := []int{0, 0, 0, 1, 1}

	if !reflect.DeepEqual(componentIds, expected) {
		t.Fatalf("unexpected cycle-cover component ids\nwant: %v\n got: %v", expected, componentIds)
	}
}

func TestSelectCmsaCycleCoverConnectorsKeepsOneIncomingAndOutgoingPerComponent(t *testing.T) {
	cmsa := [][]float64{
		{0, 0, 5, 0, 4, 0},
		{0, 0, 0, 5, 0, 0},
		{0, 0, 0, 0, 5, 0},
		{0, 0, 0, 0, 0, 5},
		{5, 0, 0, 0, 0, 0},
		{0, 5, 0, 0, 0, 0},
	}
	distances := [][]float64{
		{0, 1, 9, 0, 1, 0},
		{1, 0, 0, 1, 0, 0},
		{0, 0, 0, 1, 1, 0},
		{0, 0, 1, 0, 0, 9},
		{9, 0, 0, 0, 0, 1},
		{0, 1, 0, 0, 1, 0},
	}
	componentIds := []int{0, 0, 1, 1, 2, 2}

	connectors := selectCmsaCycleCoverConnectors(cmsa, distances, componentIds)
	expected := map[cmsaConnectorEdge]bool{
		{from: 1, to: 3}: true,
		{from: 2, to: 4}: true,
		{from: 5, to: 1}: true,
	}

	if !reflect.DeepEqual(connectors, expected) {
		t.Fatalf("unexpected CMSA connectors\nwant: %v\n got: %v", expected, connectors)
	}
}

func TestSelectCmsaCycleCoverConnectorsReturnsNoneForSingleCycle(t *testing.T) {
	cmsa := [][]float64{
		{0, 3, 3, 3},
		{3, 0, 3, 3},
		{3, 3, 0, 3},
		{3, 3, 3, 0},
	}
	componentIds := []int{0, 0, 0, 0}

	connectors := selectCmsaCycleCoverConnectors(cmsa, cmsa, componentIds)
	if len(connectors) != 0 {
		t.Fatalf("expected no connectors for a single cycle-cover component, got %v", connectors)
	}
}

func TestBuildCmsaCycleCoverHeuristicModifiersBoostsCycleCoverEdgesMoreThanHighCmsaConnectors(t *testing.T) {
	cmsa := [][]float64{
		{0, 0, 4, 4, 0},
		{0, 0, 4, 0, 0},
		{0, 0, 0, 2, 0},
		{0, 0, 0, 0, 0},
		{4, 0, 0, 0, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0, 0},
		{0, 0, 1, 0, 0},
		{1, 0, 0, 0, 0},
		{0, 0, 0, 0, 1},
		{0, 0, 0, 1, 0},
	}

	modifiers := buildCmsaCycleCoverHeuristicModifiers(cmsa, cycleCover, cmsa, 0.5)
	expected := [][]float64{
		{1, 1.5, 1, 1.25, 1},
		{1, 1, 1.5, 1, 1},
		{1.5, 1, 1, 1, 1},
		{1, 1, 1, 1, 1.5},
		{1.25, 1, 1, 1.5, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected combined modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCmsaCycleCoverHeuristicModifiersBoostsOnlySelectedConnectors(t *testing.T) {
	cmsa := [][]float64{
		{0, 0, 5, 0, 4, 0},
		{0, 0, 0, 5, 0, 0},
		{0, 0, 0, 0, 5, 0},
		{0, 0, 0, 0, 0, 5},
		{5, 0, 0, 0, 0, 0},
		{0, 5, 0, 0, 0, 0},
	}
	distances := [][]float64{
		{0, 1, 9, 0, 1, 0},
		{1, 0, 0, 1, 0, 0},
		{0, 0, 0, 1, 1, 0},
		{0, 0, 1, 0, 0, 9},
		{9, 0, 0, 0, 0, 1},
		{0, 1, 0, 0, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 1, 0},
	}

	modifiers := buildCmsaCycleCoverHeuristicModifiers(cmsa, cycleCover, distances, 0.5)

	if modifiers[0][1] != 1.5 {
		t.Fatalf("expected cycle-cover edge to receive full boost, got %f", modifiers[0][1])
	}
	if modifiers[1][3] != 1.25 || modifiers[2][4] != 1.25 || modifiers[5][1] != 1.25 {
		t.Fatalf("expected selected CMSA connectors to receive weak boost")
	}
	if modifiers[0][2] != 1.0 || modifiers[3][5] != 1.0 || modifiers[4][0] != 1.0 {
		t.Fatalf("expected unselected high-CMSA crossing edges to stay neutral")
	}
	if modifiers[0][4] != 1.0 {
		t.Fatalf("expected below-threshold CMSA crossing edge to stay neutral")
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

	if resultFilePathForHeuristic(atspData, heuristicBoth) != filepath.Join(resultsDirectoryName, "test", "result_both.csv") {
		t.Fatalf("unexpected combined result path: %s", resultFilePathForHeuristic(atspData, heuristicBoth))
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
		filepath.Join("tsplib_files", "ftv170.atsp"),
	}

	selected, err := selectAtspFiles(paths, instanceSetSmoke)
	if err != nil {
		t.Fatalf("selectAtspFiles returned unexpected error: %v", err)
	}

	if len(selected) != 1 || filepath.Base(selected[0]) != "ftv170.atsp" {
		t.Fatalf("expected smoke to select only ftv170.atsp, got %v", selected)
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
