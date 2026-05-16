package msaSupport

import (
	"atsp_aco_msa/modules/algorithms/edmonds"
	"atsp_aco_msa/modules/models"
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
		t.Fatalf("create clean MSA support: %v", err)
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
		t.Fatalf("create partially cached MSA support: %v", err)
	}

	if !compareMatrices(actual, expected) {
		t.Fatalf("expected partially cached MSA support %v to equal clean MSA support %v", actual, expected)
	}
}

func TestCreateAntiBuildsCompositeFromAllRootAntiArborescences(t *testing.T) {
	matrix := [][]float64{
		{0, 1, 8, 6},
		{4, 0, 2, 7},
		{6, 5, 0, 1},
		{9, 6, 3, 0},
	}
	rootDir := t.TempDir()

	actual, err := CreateAnti(matrix, rootDir)
	if err != nil {
		t.Fatalf("create anti-arborescence support: %v", err)
	}

	expected := expectedComposite(matrix, edmonds.FindMSAA)
	if !compareMatrices(actual, expected) {
		t.Fatalf("unexpected anti-arborescence support\nwant: %v\n got: %v", expected, actual)
	}

	readBack, err := ReadAnti(rootDir)
	if err != nil {
		t.Fatalf("read anti-arborescence support: %v", err)
	}
	if !compareMatrices(readBack, expected) {
		t.Fatalf("unexpected read anti-arborescence support\nwant: %v\n got: %v", expected, readBack)
	}

	if _, err := os.Stat(filepath.Join(rootDir, "msaas", "0.csv")); err != nil {
		t.Fatalf("expected antiarborescence cache file: %v", err)
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

func expectedComposite(matrix [][]float64, find finder) [][]float64 {
	vertices, edges, weights := models.ConvertToEdges(matrix)
	dimension := len(matrix)
	composite := make([][]float64, dimension)
	for i := range composite {
		composite[i] = make([]float64, dimension)
	}

	for root := 0; root < dimension; root++ {
		for _, edge := range find(root, vertices, edges, weights) {
			composite[edge.From][edge.To]++
		}
	}

	return composite
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
