package hungarian

import (
	"atsp_aco_msa/modules/models"
	"math"
	"reflect"
	"testing"
)

const floatTolerance = 1e-9

func TestSolveReferenceExamples(t *testing.T) {
	tests := []struct {
		name           string
		costs          [][]float64
		wantAssignment []int
		wantCost       float64
	}{
		{
			name: "oddg three by three",
			costs: [][]float64{
				{11, 6, 12},
				{12, 4, 6},
				{8, 12, 11},
			},
			wantAssignment: []int{1, 2, 0},
			wantCost:       20,
		},
		{
			name: "oddg six by six",
			costs: [][]float64{
				{13, 13, 19, 50, 33, 38},
				{73, 33, 71, 77, 97, 95},
				{20, 8, 56, 55, 64, 35},
				{26, 25, 72, 32, 55, 77},
				{83, 40, 69, 3, 53, 49},
				{67, 20, 44, 29, 86, 61},
			},
			wantAssignment: []int{4, 1, 5, 0, 3, 2},
			wantCost:       174,
		},
		{
			name: "tdedecko cost example",
			costs: [][]float64{
				{4, 2, 8},
				{4, 3, 7},
				{3, 1, 6},
			},
			wantAssignment: []int{1, 0, 2},
			wantCost:       12,
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			assignment, cost, err := Solve(test.costs)
			if err != nil {
				t.Fatalf("Solve returned error: %v", err)
			}

			assertAssignment(t, assignment, test.wantAssignment)
			assertFloat(t, cost, test.wantCost)
			assertValidAssignment(t, assignment)
		})
	}
}

func TestSolveMatchesBruteForceOnSmallMatrices(t *testing.T) {
	for n := 1; n <= 7; n++ {
		for matrixIndex := 0; matrixIndex < 15; matrixIndex++ {
			costs := deterministicCostMatrix(n, matrixIndex)

			assignment, cost, err := Solve(costs)
			if err != nil {
				t.Fatalf("Solve returned error for n=%d matrix=%d: %v", n, matrixIndex, err)
			}

			wantCost := bruteForceAssignmentCost(costs)
			assertValidAssignment(t, assignment)
			assertFloat(t, cost, wantCost)
		}
	}
}

func TestSolveKeepsValidZeroCostEdges(t *testing.T) {
	costs := [][]float64{
		{9, 0, 7},
		{5, 9, 0},
		{0, 6, 9},
	}

	assignment, cost, err := Solve(costs)
	if err != nil {
		t.Fatalf("Solve returned error: %v", err)
	}

	assertAssignment(t, assignment, []int{1, 2, 0})
	assertFloat(t, cost, 0)
}

func TestSolveDeterministicWhenSeveralOptimaExist(t *testing.T) {
	costs := [][]float64{
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
		{1, 1, 1, 1},
	}

	firstAssignment, firstCost, err := Solve(costs)
	if err != nil {
		t.Fatalf("Solve returned error: %v", err)
	}

	for i := 0; i < 20; i++ {
		assignment, cost, err := Solve(costs)
		if err != nil {
			t.Fatalf("Solve returned error on repeat %d: %v", i, err)
		}
		assertAssignment(t, assignment, firstAssignment)
		assertFloat(t, cost, firstCost)
	}
}

func TestMinimumCycleCoverForbidsSelfLoops(t *testing.T) {
	distances := [][]float64{
		{0, 5, 1},
		{1, 0, 5},
		{5, 1, 0},
	}

	edges, cost, err := MinimumCycleCover(distances)
	if err != nil {
		t.Fatalf("MinimumCycleCover returned error: %v", err)
	}

	assertEdges(t, edges, []models.Edge{
		{From: 0, To: 2},
		{From: 1, To: 0},
		{From: 2, To: 1},
	})
	assertFloat(t, cost, 3)
	assertCycleCover(t, edges, len(distances))
}

func TestMinimumCycleCoverKeepsZeroDistanceEdges(t *testing.T) {
	distances := [][]float64{
		{0, 0, 5},
		{2, 0, 0},
		{0, 3, 0},
	}

	edges, cost, err := MinimumCycleCover(distances)
	if err != nil {
		t.Fatalf("MinimumCycleCover returned error: %v", err)
	}

	assertEdges(t, edges, []models.Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
		{From: 2, To: 0},
	})
	assertFloat(t, cost, 0)
	assertCycleCover(t, edges, len(distances))
}

func TestMinimumCycleCoverForAllZeroDistances(t *testing.T) {
	distances := [][]float64{
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
	}

	edges, cost, err := MinimumCycleCover(distances)
	if err != nil {
		t.Fatalf("MinimumCycleCover returned error: %v", err)
	}

	assertFloat(t, cost, 0)
	assertCycleCover(t, edges, len(distances))
}

func TestMinimumCycleCoverMatchesBruteForceDerangement(t *testing.T) {
	for n := 2; n <= 7; n++ {
		for matrixIndex := 0; matrixIndex < 15; matrixIndex++ {
			distances := deterministicCostMatrix(n, matrixIndex+100)

			edges, cost, err := MinimumCycleCover(distances)
			if err != nil {
				t.Fatalf("MinimumCycleCover returned error for n=%d matrix=%d: %v", n, matrixIndex, err)
			}

			assertCycleCover(t, edges, n)
			assertFloat(t, cost, bruteForceDerangementCost(distances))
		}
	}
}

func TestAssignmentCycles(t *testing.T) {
	cycles, err := AssignmentCycles([]int{1, 2, 0, 4, 3})
	if err != nil {
		t.Fatalf("AssignmentCycles returned error: %v", err)
	}

	want := [][]int{
		{0, 1, 2},
		{3, 4},
	}
	if !reflect.DeepEqual(cycles, want) {
		t.Fatalf("expected cycles %v, got %v", want, cycles)
	}
}

func TestValidationRejectsInvalidInputs(t *testing.T) {
	tests := []struct {
		name  string
		costs [][]float64
	}{
		{name: "empty"},
		{name: "not square", costs: [][]float64{{1, 2}, {3}}},
		{name: "negative", costs: [][]float64{{0, -1}, {1, 0}}},
		{name: "nan", costs: [][]float64{{0, math.NaN()}, {1, 0}}},
		{name: "infinity", costs: [][]float64{{0, math.Inf(1)}, {1, 0}}},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			if _, _, err := Solve(test.costs); err == nil {
				t.Fatal("expected Solve to reject invalid input")
			}
		})
	}
}

func TestAssignmentCyclesRejectsInvalidAssignment(t *testing.T) {
	if _, err := AssignmentCycles([]int{1, 1, 0}); err == nil {
		t.Fatal("expected duplicate column assignment to be rejected")
	}

	if _, err := AssignmentCycles([]int{1, 3, 0}); err == nil {
		t.Fatal("expected out-of-range assignment to be rejected")
	}
}

func deterministicCostMatrix(n, offset int) [][]float64 {
	matrix := make([][]float64, n)
	seed := uint64(1469598103934665603 + offset*1099511628211 + n)

	for i := range matrix {
		matrix[i] = make([]float64, n)
		for j := range matrix[i] {
			seed = seed*2862933555777941757 + 3037000493
			matrix[i][j] = float64(seed % 23)
		}
	}

	return matrix
}

func bruteForceAssignmentCost(costs [][]float64) float64 {
	n := len(costs)
	used := make([]bool, n)
	best := math.Inf(1)

	var search func(row int, cost float64)
	search = func(row int, cost float64) {
		if cost >= best {
			return
		}
		if row == n {
			best = cost
			return
		}

		for column := 0; column < n; column++ {
			if used[column] {
				continue
			}

			used[column] = true
			search(row+1, cost+costs[row][column])
			used[column] = false
		}
	}

	search(0, 0)
	return best
}

func bruteForceDerangementCost(costs [][]float64) float64 {
	n := len(costs)
	used := make([]bool, n)
	best := math.Inf(1)

	var search func(row int, cost float64)
	search = func(row int, cost float64) {
		if cost >= best {
			return
		}
		if row == n {
			best = cost
			return
		}

		for column := 0; column < n; column++ {
			if used[column] || row == column {
				continue
			}

			used[column] = true
			search(row+1, cost+costs[row][column])
			used[column] = false
		}
	}

	search(0, 0)
	return best
}

func assertAssignment(t *testing.T, got, want []int) {
	t.Helper()
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("expected assignment %v, got %v", want, got)
	}
}

func assertValidAssignment(t *testing.T, assignment []int) {
	t.Helper()
	if err := validateAssignment(assignment); err != nil {
		t.Fatalf("invalid assignment %v: %v", assignment, err)
	}
}

func assertCycleCover(t *testing.T, edges []models.Edge, dimension int) {
	t.Helper()

	if len(edges) != dimension {
		t.Fatalf("expected %d edges, got %d", dimension, len(edges))
	}

	inDegree := make([]int, dimension)
	outDegree := make([]int, dimension)
	for _, edge := range edges {
		if edge.From == edge.To {
			t.Fatalf("cycle cover contains self-loop %v", edge)
		}
		outDegree[edge.From]++
		inDegree[edge.To]++
	}

	for i := 0; i < dimension; i++ {
		if inDegree[i] != 1 || outDegree[i] != 1 {
			t.Fatalf("vertex %d has in-degree %d and out-degree %d", i, inDegree[i], outDegree[i])
		}
	}
}

func assertEdges(t *testing.T, got, want []models.Edge) {
	t.Helper()
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("expected edges %v, got %v", want, got)
	}
}

func assertFloat(t *testing.T, got, want float64) {
	t.Helper()
	if math.Abs(got-want) > floatTolerance {
		t.Fatalf("expected %v, got %v", want, got)
	}
}
