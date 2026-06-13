package msaHeuristic

import (
	"os"
	"path/filepath"
	"testing"
)

func TestCreateKeepsOriginalWeightsWhenOnlySomeMsasCached(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 1},
		{1, 0, 100},
		{100, 100, 0},
	}

	cleanDir := filepath.Join(t.TempDir(), "clean")
	expected, err := Create(matrix, cleanDir)
	if err != nil {
		t.Fatalf("create clean MSA heuristic: %v", err)
	}

	cachedRootDir := filepath.Join(t.TempDir(), "cached")
	cachedMsaDir := filepath.Join(cachedRootDir, "msas")
	rootZeroMsa, err := readFromCsv(filepath.Join(cleanDir, "msas", "0.csv"))
	if err != nil {
		t.Fatalf("read seeded MSA: %v", err)
	}
	if err := os.MkdirAll(cachedMsaDir, os.ModePerm); err != nil {
		t.Fatalf("create partial MSA cache directory: %v", err)
	}
	if err := saveToCsv(rootZeroMsa, filepath.Join(cachedMsaDir, "0.csv")); err != nil {
		t.Fatalf("seed partial MSA cache: %v", err)
	}

	actual, err := Create(matrix, cachedRootDir)
	if err != nil {
		t.Fatalf("create partially cached MSA heuristic: %v", err)
	}

	if !compareMatrices(actual, expected) {
		t.Fatalf("expected partially cached MSA heuristic %v to equal clean MSA heuristic %v", actual, expected)
	}
}

func TestReadMsasSortsCacheByNumericRoot(t *testing.T) {
	rootDir := t.TempDir()
	msaDir := filepath.Join(rootDir, "msas")
	if err := os.MkdirAll(msaDir, os.ModePerm); err != nil {
		t.Fatalf("create MSA cache directory: %v", err)
	}

	if err := saveToCsv([][]float64{{10}}, filepath.Join(msaDir, "10.csv")); err != nil {
		t.Fatalf("write root 10 MSA: %v", err)
	}
	if err := saveToCsv([][]float64{{2}}, filepath.Join(msaDir, "2.csv")); err != nil {
		t.Fatalf("write root 2 MSA: %v", err)
	}
	if err := saveToCsv([][]float64{{0}}, filepath.Join(msaDir, "0.csv")); err != nil {
		t.Fatalf("write root 0 MSA: %v", err)
	}

	msas, err := ReadMsas(rootDir)
	if err != nil {
		t.Fatalf("read MSAs: %v", err)
	}

	expected := []float64{0, 2, 10}
	if len(msas) != len(expected) {
		t.Fatalf("expected %d MSAs, got %d", len(expected), len(msas))
	}
	for i, expectedValue := range expected {
		if len(msas[i]) != 1 || len(msas[i][0]) != 1 || msas[i][0][0] != expectedValue {
			t.Fatalf("unexpected MSA order at index %d: %v", i, msas)
		}
	}
}

func TestReadMsaThinnessScoresSortsByRootAndScoresEachMsa(t *testing.T) {
	rootDir := t.TempDir()
	msaDir := filepath.Join(rootDir, "msas")
	if err := os.MkdirAll(msaDir, os.ModePerm); err != nil {
		t.Fatalf("create MSA cache directory: %v", err)
	}

	if err := saveToCsv([][]float64{
		{0, 1, 1},
		{0, 0, 0},
		{0, 0, 0},
	}, filepath.Join(msaDir, "10.csv")); err != nil {
		t.Fatalf("write root 10 MSA: %v", err)
	}
	if err := saveToCsv([][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, 0},
	}, filepath.Join(msaDir, "2.csv")); err != nil {
		t.Fatalf("write root 2 MSA: %v", err)
	}

	matrix := [][]float64{
		{0, 5, 7},
		{1, 0, 11},
		{13, 17, 0},
	}
	scores, err := ReadMsaThinnessScores(rootDir, matrix)
	if err != nil {
		t.Fatalf("read MSA thinness scores: %v", err)
	}

	expected := []MsaThinnessScore{
		{Root: 2, BranchSurplus: 0, MaxOutgoingDegree: 1, BranchingVertices: 0, TotalCost: 16},
		{Root: 10, BranchSurplus: 1, MaxOutgoingDegree: 2, BranchingVertices: 1, TotalCost: 12},
	}
	if len(scores) != len(expected) {
		t.Fatalf("expected %d scores, got %d", len(expected), len(scores))
	}
	for i := range expected {
		if scores[i] != expected[i] {
			t.Fatalf("unexpected score at %d: want %+v got %+v", i, expected[i], scores[i])
		}
	}
}

func TestSelectThinnestMsaPrefersFewerBranchingVertices(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 1, 1},
		{1, 0, 1, 1},
		{1, 1, 0, 1},
		{1, 1, 1, 0},
	}
	branching := [][]float64{
		{0, 1, 1, 1},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
	}
	chain := [][]float64{
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
		{0, 0, 0, 0},
	}

	selected := SelectThinnestMsa([][][]float64{branching, chain}, matrix)
	expected := [][]float64{
		{0, 3, 0, 0},
		{0, 0, 3, 0},
		{0, 0, 0, 3},
		{0, 0, 0, 0},
	}

	if !compareMatrices(selected, expected) {
		t.Fatalf("expected chain MSA scaled to full support, got %v", selected)
	}
}

func TestSelectThinnestMsaDoesNotPreferHighDegreeHub(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 1, 1, 1},
		{1, 0, 1, 1, 1},
		{1, 1, 0, 1, 1},
		{1, 1, 1, 0, 1},
		{1, 1, 1, 1, 0},
	}
	hub := [][]float64{
		{0, 1, 1, 1, 1},
		{0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0},
	}
	twoSmallBranches := [][]float64{
		{0, 1, 1, 0, 0},
		{0, 0, 0, 1, 0},
		{0, 0, 0, 0, 1},
		{0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0},
	}

	selected := SelectThinnestMsa([][][]float64{hub, twoSmallBranches}, matrix)
	expected := [][]float64{
		{0, 4, 4, 0, 0},
		{0, 0, 0, 4, 0},
		{0, 0, 0, 0, 4},
		{0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0},
	}

	if !compareMatrices(selected, expected) {
		t.Fatalf("expected lower-surplus MSA scaled to full support, got %v", selected)
	}
}

func TestSelectThinnestMsaTieBreaksByCost(t *testing.T) {
	matrix := [][]float64{
		{0, 10, 1, 10},
		{10, 0, 10, 1},
		{10, 1, 0, 10},
		{10, 10, 10, 0},
	}
	expensive := [][]float64{
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
		{0, 0, 0, 0},
	}
	cheap := [][]float64{
		{0, 0, 1, 0},
		{0, 0, 0, 1},
		{0, 1, 0, 0},
		{0, 0, 0, 0},
	}

	selected := SelectThinnestMsa([][][]float64{expensive, cheap}, matrix)
	expected := [][]float64{
		{0, 0, 3, 0},
		{0, 0, 0, 3},
		{0, 3, 0, 0},
		{0, 0, 0, 0},
	}

	if !compareMatrices(selected, expected) {
		t.Fatalf("expected cheaper MSA scaled to full support, got %v", selected)
	}
}

func TestSelectThinnestMsaTieBreaksByRootOrder(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}
	first := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 0, 0},
	}
	second := [][]float64{
		{0, 0, 1},
		{1, 0, 0},
		{0, 0, 0},
	}

	selected := SelectThinnestMsa([][][]float64{first, second}, matrix)
	expected := [][]float64{
		{0, 2, 0},
		{0, 0, 2},
		{0, 0, 0},
	}

	if !compareMatrices(selected, expected) {
		t.Fatalf("expected first equal-scoring MSA scaled to full support, got %v", selected)
	}
}

func compareMatrices(a, b [][]float64) bool {
	if len(a) != len(b) {
		return false
	}

	for i := range a {
		if len(a[i]) != len(b[i]) {
			return false
		}

		for j := range a[i] {
			if a[i][j] != b[i][j] {
				return false
			}
		}
	}

	return true
}
