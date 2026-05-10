package main

import (
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

func TestBuildCmsaAgreementMatrixUsesMinimumSupport(t *testing.T) {
	cmsa := [][]float64{
		{0, 3, 1},
		{2, 0, 4},
		{5, 1, 0},
	}
	cmsaa := [][]float64{
		{0, 1, 2},
		{3, 0, 1},
		{4, 2, 0},
	}

	agreement, err := buildCmsaAgreementMatrix(cmsa, cmsaa)
	if err != nil {
		t.Fatalf("buildCmsaAgreementMatrix returned unexpected error: %v", err)
	}

	expected := [][]float64{
		{0, 1, 1},
		{2, 0, 1},
		{4, 1, 0},
	}
	if !reflect.DeepEqual(agreement, expected) {
		t.Fatalf("unexpected CMSA agreement matrix\nwant: %v\n got: %v", expected, agreement)
	}
}

func TestBuildCmsaAgreementMatrixRejectsMismatchedDimensions(t *testing.T) {
	cmsa := [][]float64{
		{0, 1},
		{1, 0},
	}
	cmsaa := [][]float64{
		{0, 1},
	}

	if _, err := buildCmsaAgreementMatrix(cmsa, cmsaa); err == nil {
		t.Fatal("expected mismatched CMSA/CMSAA dimensions to be rejected")
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

func TestBuildCycleCoverArborescenceHeuristicModifiersBoostsCycleCoverAndInterCycleConnectors(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 50, 50},
		{1, 0, 2, 50},
		{50, 2, 0, 1},
		{50, 50, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 0},
	}

	modifiers := buildCycleCoverArborescenceHeuristicModifiers(matrix, cycleCover, 0.5)
	expected := [][]float64{
		{1, 1.5, 1, 1},
		{1.5, 1, 1.25, 1},
		{1, 1.25, 1, 1.5},
		{1, 1, 1.5, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected cycle-cover/arborescence modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCycleCoverComponentGraphKeepsCheapestOriginalInterComponentEdges(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 2, 4},
		{1, 0, 2, 3},
		{5, 6, 0, 1},
		{7, 8, 1, 0},
	}
	componentIds := []int{0, 0, 1, 1}

	componentGraph, originalEdges := buildCycleCoverComponentGraph(matrix, componentIds)
	expectedGraph := [][]float64{
		{0, 2},
		{5, 0},
	}
	expectedEdges := map[models.Edge]models.Edge{
		{From: 0, To: 1}: {From: 0, To: 2},
		{From: 1, To: 0}: {From: 2, To: 0},
	}

	if !reflect.DeepEqual(componentGraph, expectedGraph) {
		t.Fatalf("unexpected component graph\nwant: %v\n got: %v", expectedGraph, componentGraph)
	}
	if !reflect.DeepEqual(originalEdges, expectedEdges) {
		t.Fatalf("unexpected component edge mapping\nwant: %v\n got: %v", expectedEdges, originalEdges)
	}
}

func TestBuildCycleCoverContractedArborescenceHeuristicModifiersUsesComponentGraphConnectors(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 8, 6, 9, 10},
		{1, 0, 2, 7, 11, 12},
		{6, 5, 0, 1, 10, 8},
		{9, 6, 1, 0, 3, 9},
		{7, 8, 5, 6, 0, 1},
		{4, 9, 7, 5, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 1, 0},
	}

	modifiers := buildCycleCoverContractedArborescenceHeuristicModifiers(matrix, cycleCover, 0.5)
	expected := [][]float64{
		{1, 1.5, 1, 1, 1, 1},
		{1.5, 1, 1.25, 1, 1, 1},
		{1, 1, 1, 1.5, 1, 1},
		{1, 1, 1.5, 1, 1.25, 1},
		{1, 1, 1, 1, 1, 1.5},
		{1.25, 1, 1, 1, 1.5, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected contracted cycle-cover/arborescence modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCycleCoverSpliceArborescenceHeuristicModifiersUsesArborescenceGuidedPatches(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 8, 6, 9, 10},
		{1, 0, 2, 7, 11, 12},
		{6, 5, 0, 1, 10, 8},
		{9, 6, 1, 0, 3, 9},
		{7, 8, 5, 6, 0, 1},
		{4, 9, 7, 5, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 1, 0},
	}

	modifiers := buildCycleCoverSpliceArborescenceHeuristicModifiers(matrix, cycleCover, 0.5)
	expected := [][]float64{
		{1, 1, 1, 1.25, 1, 1},
		{1, 1, 1, 1, 1.25, 1},
		{1, 1.25, 1, 1, 1, 1},
		{1, 1, 1, 1, 1.25, 1},
		{1, 1, 1, 1, 1, 1.5},
		{1.25, 1, 1.25, 1, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected patch-aware cycle-cover/arborescence modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestFindBestCycleCoverPatchUsesDirectedTwoEdgeExchange(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 20, 50},
		{1, 0, 2, 40},
		{20, 50, 0, 1},
		{2, 40, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 0},
	}
	componentIds := []int{0, 0, 1, 1}

	patch, ok := findBestCycleCoverPatch(matrix, cycleCover, componentIds, models.Edge{From: 0, To: 1})
	if !ok {
		t.Fatal("expected a patch between two cycle-cover components")
	}

	expected := cycleCoverPatch{
		componentEdge:    models.Edge{From: 0, To: 1},
		insertedForward:  models.Edge{From: 1, To: 2},
		insertedBackward: models.Edge{From: 3, To: 0},
		removedFrom:      models.Edge{From: 1, To: 0},
		removedTo:        models.Edge{From: 3, To: 2},
		delta:            2,
		valid:            true,
	}

	if !reflect.DeepEqual(patch, expected) {
		t.Fatalf("unexpected best patch\nwant: %v\n got: %v", expected, patch)
	}
}

func TestSelectInterCycleArborescencePatchesUsesMsaAndMsaaComponentPairs(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 8, 6, 9, 10},
		{1, 0, 2, 7, 11, 12},
		{6, 5, 0, 1, 10, 8},
		{9, 6, 1, 0, 3, 9},
		{7, 8, 5, 6, 0, 1},
		{4, 9, 7, 5, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 1, 0},
	}
	componentIds := []int{0, 0, 1, 1, 2, 2}

	patches := selectInterCycleArborescencePatches(matrix, cycleCover, componentIds)
	expected := []cycleCoverPatch{
		{
			componentEdge:    models.Edge{From: 0, To: 1},
			insertedForward:  models.Edge{From: 0, To: 3},
			insertedBackward: models.Edge{From: 2, To: 1},
			removedFrom:      models.Edge{From: 0, To: 1},
			removedTo:        models.Edge{From: 2, To: 3},
			delta:            9,
			valid:            true,
		},
		{
			componentEdge:    models.Edge{From: 1, To: 2},
			insertedForward:  models.Edge{From: 3, To: 4},
			insertedBackward: models.Edge{From: 5, To: 2},
			removedFrom:      models.Edge{From: 3, To: 2},
			removedTo:        models.Edge{From: 5, To: 4},
			delta:            8,
			valid:            true,
		},
		{
			componentEdge:    models.Edge{From: 2, To: 0},
			insertedForward:  models.Edge{From: 5, To: 0},
			insertedBackward: models.Edge{From: 1, To: 4},
			removedFrom:      models.Edge{From: 5, To: 4},
			removedTo:        models.Edge{From: 1, To: 0},
			delta:            13,
			valid:            true,
		},
	}

	if !reflect.DeepEqual(patches, expected) {
		t.Fatalf("unexpected arborescence-guided patches\nwant: %v\n got: %v", expected, patches)
	}
}

func TestEdgeConnectsDifferentComponentsRejectsInternalEdges(t *testing.T) {
	componentIds := []int{0, 0, 0, 1, 1}

	if edgeConnectsDifferentComponents(models.Edge{From: 0, To: 2}, componentIds) {
		t.Fatal("edge inside one cycle-cover component should not be used as connector")
	}

	if !edgeConnectsDifferentComponents(models.Edge{From: 2, To: 3}, componentIds) {
		t.Fatal("edge between cycle-cover components should be used as connector")
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

	if resultFilePathForHeuristic(atspData, heuristicCmsaa) != filepath.Join(resultsDirectoryName, "test", "result_cmsaa.csv") {
		t.Fatalf("unexpected CMSAA result path: %s", resultFilePathForHeuristic(atspData, heuristicCmsaa))
	}

	if resultFilePathForHeuristic(atspData, heuristicCmsaAgreement) != filepath.Join(resultsDirectoryName, "test", "result_cmsa_agreement.csv") {
		t.Fatalf("unexpected CMSA agreement result path: %s", resultFilePathForHeuristic(atspData, heuristicCmsaAgreement))
	}

	if resultFilePathForHeuristic(atspData, heuristicCycleCover) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover.csv") {
		t.Fatalf("unexpected cycle-cover result path: %s", resultFilePathForHeuristic(atspData, heuristicCycleCover))
	}

	if resultFilePathForHeuristic(atspData, heuristicBoth) != filepath.Join(resultsDirectoryName, "test", "result_both.csv") {
		t.Fatalf("unexpected combined result path: %s", resultFilePathForHeuristic(atspData, heuristicBoth))
	}

	if resultFilePathForHeuristic(atspData, heuristicArborescences) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover_arborescences.csv") {
		t.Fatalf("unexpected arborescence result path: %s", resultFilePathForHeuristic(atspData, heuristicArborescences))
	}

	if resultFilePathForHeuristic(atspData, heuristicContractedArborescences) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover_contracted_arborescences.csv") {
		t.Fatalf("unexpected contracted arborescence result path: %s", resultFilePathForHeuristic(atspData, heuristicContractedArborescences))
	}

	if resultFilePathForHeuristic(atspData, heuristicSpliceArborescences) != filepath.Join(resultsDirectoryName, "test", "result_cycle_cover_splice_arborescences.csv") {
		t.Fatalf("unexpected splice arborescence result path: %s", resultFilePathForHeuristic(atspData, heuristicSpliceArborescences))
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

func TestCmsaaHeuristicIsValidAndDoesNotUseCycleCover(t *testing.T) {
	if !isValidHeuristic(heuristicCmsaa) {
		t.Fatal("CMSAA heuristic should be valid")
	}

	if heuristicUsesCycleCover(heuristicCmsaa) {
		t.Fatal("CMSAA heuristic should not require cycle cover")
	}
}

func TestCmsaAgreementHeuristicIsValidAndUsesCmsaa(t *testing.T) {
	if !isValidHeuristic(heuristicCmsaAgreement) {
		t.Fatal("CMSA agreement heuristic should be valid")
	}

	if heuristicUsesCycleCover(heuristicCmsaAgreement) {
		t.Fatal("CMSA agreement heuristic should not require cycle cover")
	}

	if !heuristicUsesCmsaa(heuristicCmsaAgreement) {
		t.Fatal("CMSA agreement heuristic should require CMSAA artifacts")
	}
}

func TestArborescenceHeuristicUsesCycleCover(t *testing.T) {
	if !isValidHeuristic(heuristicArborescences) {
		t.Fatal("arborescence heuristic should be valid")
	}
	if !isValidHeuristic(heuristicContractedArborescences) {
		t.Fatal("contracted arborescence heuristic should be valid")
	}
	if !isValidHeuristic(heuristicSpliceArborescences) {
		t.Fatal("splice arborescence heuristic should be valid")
	}

	if !heuristicUsesCycleCover(heuristicArborescences) {
		t.Fatal("arborescence heuristic should require cycle cover")
	}
	if !heuristicUsesCycleCover(heuristicContractedArborescences) {
		t.Fatal("contracted arborescence heuristic should require cycle cover")
	}
	if !heuristicUsesCycleCover(heuristicSpliceArborescences) {
		t.Fatal("splice arborescence heuristic should require cycle cover")
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
