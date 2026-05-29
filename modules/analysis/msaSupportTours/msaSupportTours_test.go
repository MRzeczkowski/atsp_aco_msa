package msaSupportTours

import (
	"os"
	"path/filepath"
	"testing"
)

func TestAddUniqueTourNormalizesRotation(t *testing.T) {
	tours := make(map[string][]int)

	AddUniqueTour(tours, []int{2, 3, 0, 1})
	AddUniqueTour(tours, []int{0, 1, 2, 3})
	AddUniqueTour(tours, []int{0, 2, 3, 1})

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

func TestReadOptimalToursMissingFile(t *testing.T) {
	tours, err := ReadOptimalTours(filepath.Join(t.TempDir(), "missing.csv"))
	if err != nil {
		t.Fatalf("ReadOptimalTours returned unexpected error: %v", err)
	}
	if len(tours) != 0 {
		t.Fatalf("expected no tours, got %d", len(tours))
	}
}

func TestReadOptimalTours(t *testing.T) {
	path := filepath.Join(t.TempDir(), "solutions.csv")
	data := "Tour,Commonality with MSA support,Min commonality with MSA,Avg commonality with MSA,Max commonality with MSA\n" +
		"\"[0,1,2]\",100.00,100.00,100.00,100.00\n"
	if err := os.WriteFile(path, []byte(data), 0644); err != nil {
		t.Fatalf("failed to write test CSV: %v", err)
	}

	tours, err := ReadOptimalTours(path)
	if err != nil {
		t.Fatalf("ReadOptimalTours returned unexpected error: %v", err)
	}

	tour := tours["[0,1,2]"]
	if len(tour) != 3 || tour[0] != 0 || tour[1] != 1 || tour[2] != 2 {
		t.Fatalf("unexpected tour data: %v", tours)
	}
}

func TestAnalyzeInstancesWritesTourPlots(t *testing.T) {
	dir := t.TempDir()
	msaSupportDir := filepath.Join(dir, "msa_support")
	if err := os.MkdirAll(msaSupportDir, 0700); err != nil {
		t.Fatalf("failed to create MSA support dir: %v", err)
	}

	if err := os.WriteFile(filepath.Join(msaSupportDir, "msa_support.csv"), []byte("0,2,0\n0,0,2\n2,0,0\n"), 0644); err != nil {
		t.Fatalf("failed to write MSA support: %v", err)
	}

	solutionsPath := filepath.Join(dir, "solutions.csv")
	solutions := "Tour,Commonality with MSA support,Min commonality with MSA,Avg commonality with MSA,Max commonality with MSA\n" +
		"\"[0,1,2]\",100.00,100.00,100.00,100.00\n"
	if err := os.WriteFile(solutionsPath, []byte(solutions), 0644); err != nil {
		t.Fatalf("failed to write solutions CSV: %v", err)
	}

	toursHeatmapPath := filepath.Join(dir, "plots", "tours_heatmap.png")
	toursHistogramPath := filepath.Join(dir, "plots", "tours_histogram.png")
	overlapHeatmapPath := filepath.Join(dir, "plots", "msa_support_tours_overlap_heatmap.png")

	if err := AnalyzeInstances(Config{
		Instances: []InstanceConfig{
			{
				Name:                              "synthetic",
				Dimension:                         3,
				MsaSupportDirectoryPath:           msaSupportDir,
				OptimalToursCsvPath:               solutionsPath,
				ToursHeatmapPath:                  toursHeatmapPath,
				ToursHistogramPath:                toursHistogramPath,
				MsaSupportToursOverlapHeatmapPath: overlapHeatmapPath,
			},
		},
	}); err != nil {
		t.Fatalf("AnalyzeInstances returned unexpected error: %v", err)
	}

	assertFileExists(t, toursHeatmapPath)
	assertFileExists(t, overlapHeatmapPath)
	if _, err := os.Stat(toursHistogramPath); !os.IsNotExist(err) {
		t.Fatalf("expected no histogram for a single tour, got stat error %v", err)
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
