package compositeMsa

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
		t.Fatalf("create clean CMSA: %v", err)
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
		t.Fatalf("create partially cached CMSA: %v", err)
	}

	if !compareMatrices(actual, expected) {
		t.Fatalf("expected partially cached CMSA %v to equal clean CMSA %v", actual, expected)
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
