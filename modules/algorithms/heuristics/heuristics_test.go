package heuristics

import (
	"reflect"
	"testing"
)

func TestBuildMsaHeuristicModifiersBoostsOnlyHighSupportEdges(t *testing.T) {
	msaHeuristic := [][]float64{
		{0, 3, 2, 0},
		{0, 0, 3, 1},
		{1, 0, 0, 2},
		{0, 0, 0, 0},
	}

	modifiers := BuildMsaHeuristicModifiers(msaHeuristic, 0.5)
	expected := [][]float64{
		{1, 1.5, 1, 1},
		{1, 1, 1.5, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected MSA heuristic modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildMsaHeuristicModifiersReturnsNeutralMatrixWhenStrengthIsZero(t *testing.T) {
	msaHeuristic := [][]float64{
		{0, 3},
		{3, 0},
	}

	modifiers := BuildMsaHeuristicModifiers(msaHeuristic, 0)
	expected := [][]float64{
		{1, 1},
		{1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected neutral modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCycleCoverModifiersBoostsOnlyCycleCoverEdges(t *testing.T) {
	cycleCover := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}

	modifiers := BuildCycleCoverModifiers(cycleCover, 0.4)
	expected := [][]float64{
		{1, 1.4, 1},
		{1, 1, 1.4},
		{1.4, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected cycle-cover modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCycleCoverMsaPatchingModifiersUsesCycleCoverAndSparseMsaConnectors(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 9, 6},
		{1, 0, 4, 9},
		{9, 9, 0, 1},
		{2, 9, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 0},
	}
	msaHeuristic := [][]float64{
		{0, 3, 3, 3},
		{0, 0, 3, 0},
		{3, 0, 0, 3},
		{3, 3, 0, 0},
	}

	modifiers := BuildCycleCoverMsaPatchingModifiers(matrix, msaHeuristic, cycleCover, 0.5)
	expected := [][]float64{
		{1, 1.5, 1, 1},
		{1.5, 1, 1.5, 1},
		{1, 1, 1, 1.5},
		{1.5, 1, 1.5, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected cycle-cover MSA-patching modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCycleCoverMsaPatchingMatrixIgnoresIntraCycleMsaEdges(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 2},
		{1, 0, 3},
		{2, 3, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}
	msaHeuristic := [][]float64{
		{0, 2, 2},
		{2, 0, 2},
		{2, 2, 0},
	}

	patchingMatrix := BuildCycleCoverMsaPatchingMatrix(matrix, msaHeuristic, cycleCover)
	expected := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}

	if !reflect.DeepEqual(patchingMatrix, expected) {
		t.Fatalf("unexpected cycle-cover MSA-patching matrix\nwant: %v\n got: %v", expected, patchingMatrix)
	}
}

func TestBuildNeutralModifiers(t *testing.T) {
	modifiers := BuildNeutralModifiers(3)
	expected := [][]float64{
		{1, 1, 1},
		{1, 1, 1},
		{1, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected baseline modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}
