package tours

import (
	"atsp_aco_msa/modules/artifacts/cyclecover"
	"os"
	"path/filepath"
	"testing"
)

func TestAddUniqueNormalizesRotation(t *testing.T) {
	tours := make(map[string][]int)

	AddUnique(tours, []int{2, 3, 0, 1})
	AddUnique(tours, []int{0, 1, 2, 3})
	AddUnique(tours, []int{0, 2, 3, 1})

	if len(tours) != 2 {
		t.Fatalf("expected 2 unique tours, got %d", len(tours))
	}
	if _, ok := tours["[0,1,2,3]"]; !ok {
		t.Fatalf("expected normalized tour key [0,1,2,3], got %v", tours)
	}
	if _, ok := tours["[0,2,3,1]"]; !ok {
		t.Fatalf("expected normalized tour key [0,2,3,1], got %v", tours)
	}
}

func TestReadOptimalMissingFile(t *testing.T) {
	tours, err := ReadOptimal(filepath.Join(t.TempDir(), "missing.csv"))
	if err != nil {
		t.Fatalf("ReadOptimal returned unexpected error: %v", err)
	}
	if len(tours) != 0 {
		t.Fatalf("expected no tours, got %d", len(tours))
	}
}

func TestReadOptimal(t *testing.T) {
	path := filepath.Join(t.TempDir(), "solutions.csv")
	data := "Tour,Commonality with MSA heuristic,Min commonality with MSA,Avg commonality with MSA,Max commonality with MSA\n" +
		"\"[0,1,2]\",100.00,100.00,100.00,100.00\n"
	if err := os.WriteFile(path, []byte(data), 0644); err != nil {
		t.Fatalf("failed to write test CSV: %v", err)
	}

	tours, err := ReadOptimal(path)
	if err != nil {
		t.Fatalf("ReadOptimal returned unexpected error: %v", err)
	}

	tour := tours["[0,1,2]"]
	if len(tour) != 3 || tour[0] != 0 || tour[1] != 1 || tour[2] != 2 {
		t.Fatalf("unexpected tour data: %v", tours)
	}
}

func TestAnalyzeInstancesWritesTourPlots(t *testing.T) {
	dir := t.TempDir()
	msaHeuristicDir := filepath.Join(dir, "msa_heuristic")
	if err := os.MkdirAll(msaHeuristicDir, 0700); err != nil {
		t.Fatalf("failed to create MSA heuristic dir: %v", err)
	}

	if err := os.WriteFile(filepath.Join(msaHeuristicDir, "msa_heuristic.csv"), []byte("0,2,0\n0,0,2\n2,0,0\n"), 0644); err != nil {
		t.Fatalf("failed to write MSA heuristic: %v", err)
	}

	cycleCoverDir := filepath.Join(dir, "cycle_cover")
	if err := cyclecover.Save(cycleCoverDir, [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}); err != nil {
		t.Fatalf("failed to write cycle cover: %v", err)
	}

	solutionsPath := filepath.Join(dir, "solutions.csv")
	solutions := "Tour,Commonality with MSA heuristic,Min commonality with MSA,Avg commonality with MSA,Max commonality with MSA\n" +
		"\"[0,1,2]\",100.00,100.00,100.00,100.00\n"
	if err := os.WriteFile(solutionsPath, []byte(solutions), 0644); err != nil {
		t.Fatalf("failed to write solutions CSV: %v", err)
	}

	toursHeatmapPath := filepath.Join(dir, "plots", "tours_heatmap.png")
	toursHistogramPath := filepath.Join(dir, "plots", "tours_histogram.png")
	overlapHeatmapPath := filepath.Join(dir, "plots", "msa_heuristic_tours_overlap_heatmap.png")
	cycleCoverOverlapHeatmapPath := filepath.Join(dir, "plots", "cycle_cover_tours_overlap_heatmap.png")

	if err := AnalyzeInstances(Config{
		Instances: []InstanceConfig{
			{
				Name:                                "synthetic",
				Dimension:                           3,
				MsaHeuristicDirectoryPath:           msaHeuristicDir,
				CycleCoverDirectoryPath:             cycleCoverDir,
				OptimalToursCsvPath:                 solutionsPath,
				ToursHeatmapPath:                    toursHeatmapPath,
				ToursHistogramPath:                  toursHistogramPath,
				MsaHeuristicToursOverlapHeatmapPath: overlapHeatmapPath,
				CycleCoverToursOverlapHeatmapPath:   cycleCoverOverlapHeatmapPath,
			},
		},
	}); err != nil {
		t.Fatalf("AnalyzeInstances returned unexpected error: %v", err)
	}

	assertFileExists(t, toursHeatmapPath)
	assertFileExists(t, overlapHeatmapPath)
	assertFileExists(t, cycleCoverOverlapHeatmapPath)
	if _, err := os.Stat(toursHistogramPath); !os.IsNotExist(err) {
		t.Fatalf("expected no histogram for a single tour, got stat error %v", err)
	}
}

func TestBuildCycleCoverOverlap(t *testing.T) {
	cycleCoverMatrix := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}
	toursMatrix := [][]float64{
		{0, 2, 1},
		{1, 0, 1},
		{2, 0, 0},
	}

	overlapMatrix := BuildCycleCoverOverlap(cycleCoverMatrix, toursMatrix, 2)
	expected := [][]float64{
		{0, 1, 0},
		{0, 0, 0.5},
		{1, 0, 0},
	}

	for i := range expected {
		for j := range expected[i] {
			if overlapMatrix[i][j] != expected[i][j] {
				t.Fatalf("unexpected overlap at %d,%d: want %.2f, got %.2f", i, j, expected[i][j], overlapMatrix[i][j])
			}
		}
	}
}

func assertFileExists(t *testing.T, path string) {
	t.Helper()
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected %s to exist: %v", path, err)
	}
	if info.Size() == 0 {
		t.Fatalf("expected %s to be non-empty", path)
	}
}
