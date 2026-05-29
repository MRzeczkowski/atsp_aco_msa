package heuristics

import (
	"reflect"
	"testing"
)

func TestBuildMsaHeuristicModifiersBoostsOnlyHighSupportEdges(t *testing.T) {
	msaSupport := [][]float64{
		{0, 3, 2, 0},
		{0, 0, 3, 1},
		{1, 0, 0, 2},
		{0, 0, 0, 0},
	}

	modifiers := BuildMsaHeuristicModifiers(msaSupport, 0.5)
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
	msaSupport := [][]float64{
		{0, 3},
		{3, 0},
	}

	modifiers := BuildMsaHeuristicModifiers(msaSupport, 0)
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
