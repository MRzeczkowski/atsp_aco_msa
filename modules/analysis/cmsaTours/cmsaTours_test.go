package cmsaTours

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestCalculateMetricsPerfectOverlap(t *testing.T) {
	cmsa := newTestMatrix(4)
	for _, edge := range []Edge{{0, 1}, {1, 2}, {2, 3}, {3, 0}} {
		cmsa[edge.From][edge.To] = 3
	}

	metrics := CalculateMetrics(cmsa, map[string][]int{"tour": {0, 1, 2, 3}}, 0.8)

	assertInt(t, "found optimal tours", metrics.FoundOptimalTourCount, 1)
	assertInt(t, "unique found-optimal edges", metrics.UniqueFoundOptimalEdgeCount, 4)
	assertFloat(t, "found-optimal edge density", metrics.FoundOptimalEdgeDensity, 4.0/12.0)
	assertInt(t, "CMSA positive edges", metrics.CmsaPositiveEdgeCount, 4)
	assertFloat(t, "CMSA density", metrics.CmsaPositiveDensity, 4.0/12.0)
	assertFloat(t, "coverage", metrics.OptimalEdgeCmsaCoverage, 1.0)
	assertFloat(t, "precision", metrics.HighCmsaPrecision, 1.0)
	assertFloat(t, "recall", metrics.HighCmsaRecall, 1.0)
	assertFloat(t, "lift", metrics.HighCmsaLift, 3.0)
	assertInt(t, "top CMSA edges", metrics.TopCmsaEdgeCount, 3)
	assertFloat(t, "top CMSA precision", metrics.TopCmsaPrecision, 1.0)
	assertFloat(t, "top CMSA recall", metrics.TopCmsaRecall, 3.0/4.0)
	assertFloat(t, "top CMSA lift", metrics.TopCmsaLift, 3.0)
	assertFloat(t, "top CMSA average tour coverage", metrics.TopCmsaAverageTourCoverage, 3.0/4.0)
	assertFloat(t, "top CMSA max tour coverage", metrics.TopCmsaMaxTourCoverage, 3.0/4.0)
	assertFloat(t, "avg optimal CMSA", metrics.AverageNormalizedCmsaOnOptimalEdges, 1.0)
	assertFloat(t, "median optimal CMSA", metrics.MedianNormalizedCmsaOnOptimalEdges, 1.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 0)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighCmsaEdges, 0)
	assertInt(t, "false-positive top edges", metrics.TopCmsaFalsePositiveEdges, 0)
}

func TestCalculateMetricsNoOverlap(t *testing.T) {
	cmsa := newTestMatrix(4)
	for _, edge := range []Edge{{0, 2}, {1, 3}, {2, 0}, {3, 1}} {
		cmsa[edge.From][edge.To] = 3
	}

	metrics := CalculateMetrics(cmsa, map[string][]int{"tour": {0, 1, 2, 3}}, 0.8)

	assertFloat(t, "coverage", metrics.OptimalEdgeCmsaCoverage, 0.0)
	assertFloat(t, "precision", metrics.HighCmsaPrecision, 0.0)
	assertFloat(t, "recall", metrics.HighCmsaRecall, 0.0)
	assertFloat(t, "lift", metrics.HighCmsaLift, 0.0)
	assertFloat(t, "top CMSA precision", metrics.TopCmsaPrecision, 0.0)
	assertFloat(t, "top CMSA recall", metrics.TopCmsaRecall, 0.0)
	assertFloat(t, "top CMSA average tour coverage", metrics.TopCmsaAverageTourCoverage, 0.0)
	assertFloat(t, "avg optimal CMSA", metrics.AverageNormalizedCmsaOnOptimalEdges, 0.0)
	assertFloat(t, "median optimal CMSA", metrics.MedianNormalizedCmsaOnOptimalEdges, 0.0)
	assertFloat(t, "avg not-found positive CMSA", metrics.AverageNormalizedCmsaOnNotFoundPositive, 1.0)
	assertFloat(t, "median not-found positive CMSA", metrics.MedianNormalizedCmsaOnNotFoundPositive, 1.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 4)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighCmsaEdges, 4)
}

func TestCalculateMetricsPartialOverlap(t *testing.T) {
	cmsa := newTestMatrix(4)
	cmsa[0][1] = 3
	cmsa[1][2] = 1
	cmsa[0][2] = 3
	cmsa[1][3] = 3
	cmsa[2][0] = 1

	metrics := CalculateMetrics(cmsa, map[string][]int{"tour": {0, 1, 2, 3}}, 0.8)

	assertInt(t, "CMSA positive edges", metrics.CmsaPositiveEdgeCount, 5)
	assertInt(t, "high CMSA edges", metrics.HighCmsaEdgeCount, 3)
	assertFloat(t, "coverage", metrics.OptimalEdgeCmsaCoverage, 0.5)
	assertFloat(t, "precision", metrics.HighCmsaPrecision, 1.0/3.0)
	assertFloat(t, "recall", metrics.HighCmsaRecall, 0.25)
	assertFloat(t, "lift", metrics.HighCmsaLift, 1.0)
	assertInt(t, "top CMSA edges", metrics.TopCmsaEdgeCount, 3)
	assertFloat(t, "top CMSA precision", metrics.TopCmsaPrecision, 1.0/3.0)
	assertFloat(t, "top CMSA recall", metrics.TopCmsaRecall, 0.25)
	assertFloat(t, "top CMSA lift", metrics.TopCmsaLift, 1.0)
	assertFloat(t, "top CMSA average tour coverage", metrics.TopCmsaAverageTourCoverage, 0.25)
	assertFloat(t, "avg optimal CMSA", metrics.AverageNormalizedCmsaOnOptimalEdges, 1.0/3.0)
	assertFloat(t, "median optimal CMSA", metrics.MedianNormalizedCmsaOnOptimalEdges, 1.0/6.0)
	assertFloat(t, "avg not-found positive CMSA", metrics.AverageNormalizedCmsaOnNotFoundPositive, 7.0/9.0)
	assertFloat(t, "median not-found positive CMSA", metrics.MedianNormalizedCmsaOnNotFoundPositive, 1.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 2)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighCmsaEdges, 2)
}

func TestCalculateMetricsNoFoundOptimalTours(t *testing.T) {
	cmsa := newTestMatrix(4)
	cmsa[0][1] = 3
	cmsa[1][2] = 1

	metrics := CalculateMetrics(cmsa, nil, 0.8)

	assertInt(t, "found optimal tours", metrics.FoundOptimalTourCount, 0)
	assertInt(t, "unique found-optimal edges", metrics.UniqueFoundOptimalEdgeCount, 0)
	assertInt(t, "CMSA positive edges", metrics.CmsaPositiveEdgeCount, 2)
	assertInt(t, "high CMSA edges", metrics.HighCmsaEdgeCount, 1)
	assertFloat(t, "coverage", metrics.OptimalEdgeCmsaCoverage, 0.0)
	assertFloat(t, "precision", metrics.HighCmsaPrecision, 0.0)
	assertFloat(t, "recall", metrics.HighCmsaRecall, 0.0)
	assertFloat(t, "lift", metrics.HighCmsaLift, 0.0)
	assertInt(t, "top CMSA edges", metrics.TopCmsaEdgeCount, 2)
	assertFloat(t, "top CMSA precision", metrics.TopCmsaPrecision, 0.0)
	assertFloat(t, "top CMSA recall", metrics.TopCmsaRecall, 0.0)
	assertFloat(t, "avg not-found positive CMSA", metrics.AverageNormalizedCmsaOnNotFoundPositive, 0.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 0)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighCmsaEdges, 0)
}

func TestCalculateThresholdMetricsOrdering(t *testing.T) {
	cmsa := newTestMatrix(4)
	cmsa[0][1] = 3
	cmsa[1][2] = 1

	metrics := CalculateThresholdMetrics(cmsa, map[string][]int{"tour": {0, 1, 2, 3}}, []float64{0.8, 0.1, 1.0})

	if len(metrics) != 3 {
		t.Fatalf("expected 3 threshold rows, got %d", len(metrics))
	}

	assertFloat(t, "first threshold", metrics[0].Threshold, 0.1)
	assertFloat(t, "second threshold", metrics[1].Threshold, 0.8)
	assertFloat(t, "third threshold", metrics[2].Threshold, 1.0)
	assertInt(t, "0.1 high edges", metrics[0].HighCmsaEdgeCount, 2)
	assertInt(t, "0.8 high edges", metrics[1].HighCmsaEdgeCount, 1)
}

func TestDefaultThresholds(t *testing.T) {
	thresholds := DefaultThresholds()
	if len(thresholds) != 10 {
		t.Fatalf("expected 10 thresholds, got %d", len(thresholds))
	}

	for i, threshold := range thresholds {
		assertFloat(t, "threshold", threshold, float64(i+1)/10.0)
	}
}

func TestAnalyzeInstancesNoFoundToursSmoke(t *testing.T) {
	dir := t.TempDir()
	cmsaDir := filepath.Join(dir, "cmsa")
	if err := os.MkdirAll(cmsaDir, 0700); err != nil {
		t.Fatalf("failed to create CMSA dir: %v", err)
	}

	if err := os.WriteFile(filepath.Join(cmsaDir, "cmsa.csv"), []byte("0,1,0\n0,0,1\n1,0,0\n"), 0644); err != nil {
		t.Fatalf("failed to write CMSA: %v", err)
	}

	solutionsPath := filepath.Join(dir, "solutions.csv")
	if err := os.WriteFile(solutionsPath, []byte("Tour,Commonality with CMSA,Min commonality with MSA,Avg commonality with MSA,Max commonality with MSA\n"), 0644); err != nil {
		t.Fatalf("failed to write solutions CSV: %v", err)
	}

	summaryPath := filepath.Join(dir, "summary.csv")
	reportPath := filepath.Join(dir, "report.md")
	_, err := AnalyzeInstances(Config{
		Instances: []InstanceConfig{
			{
				Name:                "synthetic",
				Dimension:           3,
				CmsaDirectoryPath:   cmsaDir,
				OptimalToursCsvPath: solutionsPath,
				AnalysisCsvPath:     filepath.Join(dir, "cmsa_solution_analysis.csv"),
				ThresholdsCsvPath:   filepath.Join(dir, "cmsa_solution_thresholds.csv"),
			},
		},
		SummaryCsvPath: summaryPath,
		ReportPath:     reportPath,
		HighThreshold:  0.8,
		Thresholds:     []float64{0.8, 0.1},
	})
	if err != nil {
		t.Fatalf("AnalyzeInstances returned unexpected error: %v", err)
	}

	summary, err := os.ReadFile(summaryPath)
	if err != nil {
		t.Fatalf("failed to read summary: %v", err)
	}
	if !strings.Contains(string(summary), "synthetic,0,0,0.00,0.800000,0") {
		t.Fatalf("summary does not contain expected no-tour row:\n%s", string(summary))
	}

	thresholds, err := os.ReadFile(filepath.Join(dir, "cmsa_solution_thresholds.csv"))
	if err != nil {
		t.Fatalf("failed to read thresholds: %v", err)
	}
	if !strings.Contains(string(thresholds), "0.100000,3,50.00") {
		t.Fatalf("thresholds do not contain expected sorted 0.1 row:\n%s", string(thresholds))
	}

	if _, err := os.Stat(reportPath); err != nil {
		t.Fatalf("expected report to be written: %v", err)
	}
}

func newTestMatrix(dimension int) [][]float64 {
	matrix := make([][]float64, dimension)
	for i := range matrix {
		matrix[i] = make([]float64, dimension)
	}
	return matrix
}

func assertFloat(t *testing.T, name string, got, want float64) {
	t.Helper()
	if math.Abs(got-want) > 1e-9 {
		t.Fatalf("%s mismatch: got %.12f want %.12f", name, got, want)
	}
}

func assertInt(t *testing.T, name string, got, want int) {
	t.Helper()
	if got != want {
		t.Fatalf("%s mismatch: got %d want %d", name, got, want)
	}
}
