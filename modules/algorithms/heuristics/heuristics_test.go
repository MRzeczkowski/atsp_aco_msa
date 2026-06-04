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

func TestBuildRandomSparseModifiersBoostsSameNumberOfEdgesAsMsaHeuristic(t *testing.T) {
	msaHeuristic := [][]float64{
		{0, 4, 3, 0, 0},
		{0, 0, 4, 1, 0},
		{1, 0, 0, 4, 0},
		{0, 0, 0, 0, 4},
		{4, 0, 0, 0, 0},
	}

	msaModifiers := BuildMsaHeuristicModifiers(msaHeuristic, 0.9)
	randomModifiers := BuildRandomSparseModifiers(msaHeuristic, 0.9, 17)

	expectedBoostedEdges := boostedModifierEdgeCount(msaModifiers)
	actualBoostedEdges := boostedModifierEdgeCount(randomModifiers)
	if actualBoostedEdges != expectedBoostedEdges {
		t.Fatalf("unexpected boosted-edge count: want %d, got %d", expectedBoostedEdges, actualBoostedEdges)
	}

	for i := range randomModifiers {
		if randomModifiers[i][i] != 1.0 {
			t.Fatalf("random sparse modifier boosted self-loop %d/%d", i, i)
		}
		for j := range randomModifiers[i] {
			if randomModifiers[i][j] != 1.0 && randomModifiers[i][j] != 1.9 {
				t.Fatalf("unexpected modifier at %d/%d: %f", i, j, randomModifiers[i][j])
			}
		}
	}
}

func TestBuildRandomSparseModifiersIsDeterministicForSeed(t *testing.T) {
	msaHeuristic := [][]float64{
		{0, 4, 3, 0, 0},
		{0, 0, 4, 1, 0},
		{1, 0, 0, 4, 0},
		{0, 0, 0, 0, 4},
		{4, 0, 0, 0, 0},
	}

	first := BuildRandomSparseModifiers(msaHeuristic, 0.5, 42)
	second := BuildRandomSparseModifiers(msaHeuristic, 0.5, 42)

	if !reflect.DeepEqual(first, second) {
		t.Fatalf("expected deterministic random sparse modifiers\nfirst:  %v\nsecond: %v", first, second)
	}
}

func TestBuildRandomSparseModifiersReturnsNeutralMatrixWhenStrengthIsZero(t *testing.T) {
	msaHeuristic := [][]float64{
		{0, 2, 0},
		{0, 0, 2},
		{2, 0, 0},
	}

	modifiers := BuildRandomSparseModifiers(msaHeuristic, 0, 1)
	expected := BuildNeutralModifiers(3)

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

func TestBestGreedyCyclePatchCanPatchAnyPairOfCycles(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 10, 10, 10, 10},
		{1, 0, 10, 10, 10, 10},
		{10, 10, 0, 1, 10, 10},
		{10, 10, 1, 0, 10, 1},
		{10, 10, 1, 10, 0, 1},
		{10, 10, 10, 10, 1, 0},
	}
	successors := []int{1, 0, 3, 2, 5, 4}
	cycles := buildCycles(successors)
	msaHeuristic := newZeroMatrix(len(successors))

	patch := bestGreedyCyclePatch(matrix, msaHeuristic, successors, cycles, 1.0)

	if !patch.valid {
		t.Fatal("expected a valid patch")
	}
	if patch.fromA != 3 || patch.fromB != 4 {
		t.Fatalf("expected greedy patch from cycles 2/3 and 4/5 to be 3/4, got %d/%d", patch.fromA, patch.fromB)
	}
}

func TestPatchCostCacheMatchesFullScanAfterPatchUpdate(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 10, 10, 10, 10},
		{1, 0, 3, 10, 4, 10},
		{10, 4, 0, 1, 10, 10},
		{3, 10, 1, 0, 10, 5},
		{10, 10, 3, 10, 0, 1},
		{10, 4, 10, 10, 1, 0},
	}
	successors := []int{1, 0, 3, 2, 5, 4}
	cycles := buildCycles(successors)
	cycleIndex := vertexCycleIndexes(cycles, len(successors))
	msaHeuristic := newZeroMatrix(len(successors))
	cache := newPatchCostCache(matrix, msaHeuristic, successors, cycleIndex, 1.0)

	patch := cache.bestPatch()
	expected := bestGreedyCyclePatch(matrix, msaHeuristic, successors, cycles, 1.0)

	if !samePatchEndpoints(patch, expected.fromA, expected.fromB) {
		t.Fatalf("cached patch does not match full scan before update: cache=%+v scan=%+v", patch, expected)
	}

	fromCycle := cycleIndex[patch.fromA]
	toCycle := cycleIndex[patch.fromB]
	fromCycleVertices := cycles[fromCycle]
	toCycleVertices := cycles[toCycle]

	applyCyclePatch(successors, patch)
	for _, vertex := range toCycleVertices {
		cycleIndex[vertex] = fromCycle
	}
	cycles[fromCycle] = append(cycles[fromCycle], toCycleVertices...)
	cycles[toCycle] = nil
	cache.updateAfterPatch(patch, fromCycleVertices, toCycleVertices, fromCycle)

	currentCycles := buildCycles(successors)
	patch = cache.bestPatch()
	expected = bestGreedyCyclePatch(matrix, msaHeuristic, successors, currentCycles, 1.0)

	if !samePatchEndpoints(patch, expected.fromA, expected.fromB) {
		t.Fatalf("cached patch does not match full scan after update: cache=%+v scan=%+v", patch, expected)
	}
}

func TestCycleCoverMsaPatchingMatrixCacheMatchesFullScanReference(t *testing.T) {
	for dimension := 4; dimension <= 10; dimension++ {
		matrix := deterministicAsymmetricMatrix(dimension)
		msaHeuristic := deterministicMsaHeuristicMatrix(dimension)
		cycleCover := deterministicCycleCoverMatrix(dimension)

		cached := BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(matrix, msaHeuristic, cycleCover, 0.75)
		fullScan := buildCycleCoverMsaPatchingMatrixByFullScan(matrix, msaHeuristic, cycleCover, 0.75)

		if !reflect.DeepEqual(cached, fullScan) {
			t.Fatalf("cached patching differs from full scan for dimension %d\ncached: %v\nfull:   %v", dimension, cached, fullScan)
		}
	}
}

func TestInvalidCyclePatchNeverBeatsValidPatch(t *testing.T) {
	valid := cyclePatch{fromA: 1, fromB: 2, score: 10, valid: true}
	invalid := cyclePatch{}

	if betterCyclePatch(invalid, valid) {
		t.Fatal("invalid patch should not beat a valid patch")
	}
	if !betterCyclePatch(valid, invalid) {
		t.Fatal("valid patch should beat an invalid patch")
	}
}

func buildCycleCoverMsaPatchingMatrixByFullScan(matrix, msaHeuristic, cycleCover [][]float64, msaPatchBias float64) [][]float64 {
	dimension := heuristicDimension(matrix, msaHeuristic, cycleCover)
	successors, ok := cycleCoverSuccessors(cycleCover, dimension)
	if !ok {
		return copyPositiveEdges(cycleCover, dimension)
	}

	cycles := buildCycles(successors)
	for len(cycles) > 1 {
		patch := bestGreedyCyclePatch(matrix, msaHeuristic, successors, cycles, msaPatchBias)
		if !patch.valid {
			return copyPositiveEdges(cycleCover, dimension)
		}

		applyCyclePatch(successors, patch)
		cycles = buildCycles(successors)
	}

	return buildSuccessorMatrix(successors)
}

func deterministicAsymmetricMatrix(dimension int) [][]float64 {
	matrix := newZeroMatrix(dimension)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			matrix[i][j] = float64(((i+3)*(j+5)+7*i+j*j)%23 + 1)
		}
	}
	return matrix
}

func deterministicMsaHeuristicMatrix(dimension int) [][]float64 {
	matrix := newZeroMatrix(dimension)
	maxSupport := dimension - 1
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			matrix[i][j] = float64(((i+1)*(j+2) + i + 2*j) % maxSupport)
		}
	}
	return matrix
}

func deterministicCycleCoverMatrix(dimension int) [][]float64 {
	matrix := newZeroMatrix(dimension)
	for start := 0; start < dimension; {
		cycleLength := 2
		if dimension-start == 3 || dimension-start >= 5 {
			cycleLength = 3
		}

		for offset := 0; offset < cycleLength; offset++ {
			from := start + offset
			to := start + ((offset + 1) % cycleLength)
			matrix[from][to] = 1
		}
		start += cycleLength
	}
	return matrix
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
