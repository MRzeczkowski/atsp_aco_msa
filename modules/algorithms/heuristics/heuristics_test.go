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

func TestBuildCycleCoverMsaPatchingModifiersUsesMsaSupportInPatchScore(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 20, 3},
		{1, 0, 5, 20},
		{20, 3, 0, 1},
		{5, 20, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 0},
	}
	msaHeuristic := [][]float64{
		{0, 0, 0, 0},
		{0, 0, 2, 0},
		{0, 0, 0, 0},
		{2, 0, 0, 0},
	}

	modifiers := BuildCycleCoverMsaPatchingModifiers(matrix, msaHeuristic, cycleCover, 0.5, 1.0)
	expected := [][]float64{
		{1, 1.5, 1, 1},
		{1, 1, 1.5, 1},
		{1, 1, 1, 1.5},
		{1.5, 1, 1, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected cycle-cover MSA-patching modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCycleCoverMsaPatchingModifiersCanDisableMsaPatchBias(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 20, 3},
		{1, 0, 5, 20},
		{20, 3, 0, 1},
		{5, 20, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 0},
	}
	msaHeuristic := [][]float64{
		{0, 0, 0, 0},
		{0, 0, 2, 0},
		{0, 0, 0, 0},
		{2, 0, 0, 0},
	}

	modifiers := BuildCycleCoverMsaPatchingModifiers(matrix, msaHeuristic, cycleCover, 0.5, 0.0)
	expected := [][]float64{
		{1, 1, 1, 1.5},
		{1.5, 1, 1, 1},
		{1, 1.5, 1, 1},
		{1, 1, 1.5, 1},
	}

	if !reflect.DeepEqual(modifiers, expected) {
		t.Fatalf("unexpected cost-only patching modifiers\nwant: %v\n got: %v", expected, modifiers)
	}
}

func TestBuildCycleCoverMsaPatchingMatrixAppliesTwoEdgePatch(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 20, 1},
		{1, 0, 20, 20},
		{20, 1, 0, 1},
		{20, 20, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 0},
	}
	msaHeuristic := [][]float64{
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
	}

	patchingMatrix := BuildCycleCoverMsaPatchingMatrix(matrix, msaHeuristic, cycleCover)
	expected := [][]float64{
		{0, 0, 0, 1},
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
	}

	if !reflect.DeepEqual(patchingMatrix, expected) {
		t.Fatalf("unexpected patching matrix\nwant: %v\n got: %v", expected, patchingMatrix)
	}
}

func TestLargestCycleSelectsLargestCycleDeterministically(t *testing.T) {
	cycles := [][]int{
		{0, 1},
		{2, 3, 4},
		{5, 6, 7},
	}

	selected := largestCycle(cycles)

	if !reflect.DeepEqual(selected, []int{2, 3, 4}) {
		t.Fatalf("expected first largest cycle to be selected, got %v", selected)
	}
}

func TestBuildCycleCoverMsaPatchingMatrixReturnsSingleTourWithValidDegrees(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 2, 3, 4, 5},
		{1, 0, 2, 3, 4, 5},
		{2, 3, 0, 1, 4, 5},
		{2, 3, 1, 0, 4, 5},
		{2, 3, 4, 5, 0, 1},
		{2, 3, 4, 5, 1, 0},
	}
	cycleCover := [][]float64{
		{0, 1, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 1},
		{0, 0, 0, 0, 1, 0},
	}
	msaHeuristic := [][]float64{
		{0, 0, 5, 0, 0, 0},
		{0, 0, 0, 0, 5, 0},
		{0, 0, 0, 0, 5, 0},
		{5, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 5, 0, 0, 0, 0},
	}

	patchingMatrix := BuildCycleCoverMsaPatchingMatrix(matrix, msaHeuristic, cycleCover)

	assertSingleTourMatrix(t, patchingMatrix)
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

func TestBestCyclePatchUsesOnlyVerticesFromSelectedCycle(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 6, 6, 6, 6, 6, 6},
		{1, 0, 20, 3, 20, 20, 20, 20},
		{6, 6, 0, 1, 6, 0, 6, 6},
		{6, 6, 6, 0, 1, 6, 0, 6},
		{6, 6, 1, 6, 0, 6, 6, 0},
		{6, 6, 0, 6, 6, 0, 1, 6},
		{6, 6, 6, 0, 6, 6, 0, 1},
		{6, 6, 6, 6, 0, 1, 6, 0},
	}
	successors := []int{1, 0, 3, 4, 2, 6, 7, 5}
	selectedCycle := []int{0, 1}
	usedInPatch := make([]bool, len(successors))
	msaHeuristic := newZeroMatrix(len(successors))

	patch := bestCyclePatch(matrix, msaHeuristic, successors, selectedCycle, usedInPatch, 1.0)

	if !patch.valid {
		t.Fatal("expected a valid patch")
	}
	if patch.fromA != 1 || patch.fromB != 2 {
		t.Fatalf("expected best patch from selected cycle to be 1/2, got %d/%d", patch.fromA, patch.fromB)
	}
}

func TestBestCyclePatchDoesNotReusePatchVertices(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 1, 1},
		{1, 0, 20, 2},
		{1, 20, 0, 1},
		{1, 2, 1, 0},
	}
	successors := []int{1, 0, 3, 2}
	selectedCycle := []int{0, 1}
	usedInPatch := []bool{true, false, true, false}
	msaHeuristic := newZeroMatrix(len(successors))

	patch := bestCyclePatch(matrix, msaHeuristic, successors, selectedCycle, usedInPatch, 1.0)

	if !patch.valid {
		t.Fatal("expected a valid patch")
	}
	if patch.fromA != 1 || patch.fromB != 3 {
		t.Fatalf("expected only unused vertices 1/3 to be patched, got %d/%d", patch.fromA, patch.fromB)
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

func assertSingleTourMatrix(t *testing.T, matrix [][]float64) {
	t.Helper()

	dimension := len(matrix)
	successors := make([]int, dimension)
	inDegree := make([]int, dimension)
	for i := range successors {
		successors[i] = -1
	}

	for from, row := range matrix {
		outDegree := 0
		for to, value := range row {
			if from == to || value == 0 {
				continue
			}

			outDegree++
			successors[from] = to
			inDegree[to]++
		}
		if outDegree != 1 {
			t.Fatalf("vertex %d has out-degree %d, want 1", from, outDegree)
		}
	}

	for vertex, degree := range inDegree {
		if degree != 1 {
			t.Fatalf("vertex %d has in-degree %d, want 1", vertex, degree)
		}
	}

	visited := make([]bool, dimension)
	current := 0
	for step := 0; step < dimension; step++ {
		if current < 0 || current >= dimension {
			t.Fatalf("tour left matrix at vertex %d", current)
		}
		if visited[current] {
			t.Fatalf("tour returned to vertex %d after only %d steps", current, step)
		}
		visited[current] = true
		current = successors[current]
	}

	if current != 0 {
		t.Fatalf("tour ended at %d, want 0", current)
	}
}
