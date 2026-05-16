package msaSupportTours

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestCalculateMetricsPerfectOverlap(t *testing.T) {
	msaSupport := newTestMatrix(4)
	for _, edge := range []Edge{{0, 1}, {1, 2}, {2, 3}, {3, 0}} {
		msaSupport[edge.From][edge.To] = 3
	}

	metrics := CalculateMetrics(msaSupport, map[string][]int{"tour": {0, 1, 2, 3}}, 0.8)

	assertInt(t, "found optimal tours", metrics.FoundOptimalTourCount, 1)
	assertInt(t, "unique found-optimal edges", metrics.UniqueFoundOptimalEdgeCount, 4)
	assertFloat(t, "found-optimal edge density", metrics.FoundOptimalEdgeDensity, 4.0/12.0)
	assertInt(t, "MSA support positive edges", metrics.MsaSupportPositiveEdgeCount, 4)
	assertFloat(t, "MSA support density", metrics.MsaSupportPositiveDensity, 4.0/12.0)
	assertFloat(t, "coverage", metrics.OptimalEdgeMsaSupportCoverage, 1.0)
	assertFloat(t, "precision", metrics.HighMsaSupportPrecision, 1.0)
	assertFloat(t, "recall", metrics.HighMsaSupportRecall, 1.0)
	assertFloat(t, "lift", metrics.HighMsaSupportLift, 3.0)
	assertInt(t, "top MSA support edges", metrics.TopMsaSupportEdgeCount, 3)
	assertFloat(t, "top MSA support precision", metrics.TopMsaSupportPrecision, 1.0)
	assertFloat(t, "top MSA support recall", metrics.TopMsaSupportRecall, 3.0/4.0)
	assertFloat(t, "top MSA support lift", metrics.TopMsaSupportLift, 3.0)
	assertFloat(t, "top MSA support average tour coverage", metrics.TopMsaSupportAverageTourCoverage, 3.0/4.0)
	assertFloat(t, "top MSA support max tour coverage", metrics.TopMsaSupportMaxTourCoverage, 3.0/4.0)
	assertFloat(t, "avg optimal MSA support", metrics.AverageNormalizedMsaSupportOnOptimalEdges, 1.0)
	assertFloat(t, "median optimal MSA support", metrics.MedianNormalizedMsaSupportOnOptimalEdges, 1.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 0)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighMsaSupportEdges, 0)
	assertInt(t, "false-positive top edges", metrics.TopMsaSupportFalsePositiveEdges, 0)
}

func TestCalculateMetricsNoOverlap(t *testing.T) {
	msaSupport := newTestMatrix(4)
	for _, edge := range []Edge{{0, 2}, {1, 3}, {2, 0}, {3, 1}} {
		msaSupport[edge.From][edge.To] = 3
	}

	metrics := CalculateMetrics(msaSupport, map[string][]int{"tour": {0, 1, 2, 3}}, 0.8)

	assertFloat(t, "coverage", metrics.OptimalEdgeMsaSupportCoverage, 0.0)
	assertFloat(t, "precision", metrics.HighMsaSupportPrecision, 0.0)
	assertFloat(t, "recall", metrics.HighMsaSupportRecall, 0.0)
	assertFloat(t, "lift", metrics.HighMsaSupportLift, 0.0)
	assertFloat(t, "top MSA support precision", metrics.TopMsaSupportPrecision, 0.0)
	assertFloat(t, "top MSA support recall", metrics.TopMsaSupportRecall, 0.0)
	assertFloat(t, "top MSA support average tour coverage", metrics.TopMsaSupportAverageTourCoverage, 0.0)
	assertFloat(t, "avg optimal MSA support", metrics.AverageNormalizedMsaSupportOnOptimalEdges, 0.0)
	assertFloat(t, "median optimal MSA support", metrics.MedianNormalizedMsaSupportOnOptimalEdges, 0.0)
	assertFloat(t, "avg not-found positive MSA support", metrics.AverageNormalizedMsaSupportOnNotFoundPositive, 1.0)
	assertFloat(t, "median not-found positive MSA support", metrics.MedianNormalizedMsaSupportOnNotFoundPositive, 1.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 4)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighMsaSupportEdges, 4)
}

func TestCalculateMetricsPartialOverlap(t *testing.T) {
	msaSupport := newTestMatrix(4)
	msaSupport[0][1] = 3
	msaSupport[1][2] = 1
	msaSupport[0][2] = 3
	msaSupport[1][3] = 3
	msaSupport[2][0] = 1

	metrics := CalculateMetrics(msaSupport, map[string][]int{"tour": {0, 1, 2, 3}}, 0.8)

	assertInt(t, "MSA support positive edges", metrics.MsaSupportPositiveEdgeCount, 5)
	assertInt(t, "high MSA support edges", metrics.HighMsaSupportEdgeCount, 3)
	assertFloat(t, "coverage", metrics.OptimalEdgeMsaSupportCoverage, 0.5)
	assertFloat(t, "precision", metrics.HighMsaSupportPrecision, 1.0/3.0)
	assertFloat(t, "recall", metrics.HighMsaSupportRecall, 0.25)
	assertFloat(t, "lift", metrics.HighMsaSupportLift, 1.0)
	assertInt(t, "top MSA support edges", metrics.TopMsaSupportEdgeCount, 3)
	assertFloat(t, "top MSA support precision", metrics.TopMsaSupportPrecision, 1.0/3.0)
	assertFloat(t, "top MSA support recall", metrics.TopMsaSupportRecall, 0.25)
	assertFloat(t, "top MSA support lift", metrics.TopMsaSupportLift, 1.0)
	assertFloat(t, "top MSA support average tour coverage", metrics.TopMsaSupportAverageTourCoverage, 0.25)
	assertFloat(t, "avg optimal MSA support", metrics.AverageNormalizedMsaSupportOnOptimalEdges, 1.0/3.0)
	assertFloat(t, "median optimal MSA support", metrics.MedianNormalizedMsaSupportOnOptimalEdges, 1.0/6.0)
	assertFloat(t, "avg not-found positive MSA support", metrics.AverageNormalizedMsaSupportOnNotFoundPositive, 7.0/9.0)
	assertFloat(t, "median not-found positive MSA support", metrics.MedianNormalizedMsaSupportOnNotFoundPositive, 1.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 2)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighMsaSupportEdges, 2)
}

func TestCalculateMetricsNoFoundOptimalTours(t *testing.T) {
	msaSupport := newTestMatrix(4)
	msaSupport[0][1] = 3
	msaSupport[1][2] = 1

	metrics := CalculateMetrics(msaSupport, nil, 0.8)

	assertInt(t, "found optimal tours", metrics.FoundOptimalTourCount, 0)
	assertInt(t, "unique found-optimal edges", metrics.UniqueFoundOptimalEdgeCount, 0)
	assertInt(t, "MSA support positive edges", metrics.MsaSupportPositiveEdgeCount, 2)
	assertInt(t, "high MSA support edges", metrics.HighMsaSupportEdgeCount, 1)
	assertFloat(t, "coverage", metrics.OptimalEdgeMsaSupportCoverage, 0.0)
	assertFloat(t, "precision", metrics.HighMsaSupportPrecision, 0.0)
	assertFloat(t, "recall", metrics.HighMsaSupportRecall, 0.0)
	assertFloat(t, "lift", metrics.HighMsaSupportLift, 0.0)
	assertInt(t, "top MSA support edges", metrics.TopMsaSupportEdgeCount, 2)
	assertFloat(t, "top MSA support precision", metrics.TopMsaSupportPrecision, 0.0)
	assertFloat(t, "top MSA support recall", metrics.TopMsaSupportRecall, 0.0)
	assertFloat(t, "avg not-found positive MSA support", metrics.AverageNormalizedMsaSupportOnNotFoundPositive, 0.0)
	assertInt(t, "missing optimal edges", metrics.MissingFoundOptimalEdges, 0)
	assertInt(t, "false-positive high edges", metrics.FalsePositiveHighMsaSupportEdges, 0)
}

func TestCalculateThresholdMetricsOrdering(t *testing.T) {
	msaSupport := newTestMatrix(4)
	msaSupport[0][1] = 3
	msaSupport[1][2] = 1

	metrics := CalculateThresholdMetrics(msaSupport, map[string][]int{"tour": {0, 1, 2, 3}}, []float64{0.8, 0.1, 1.0})

	if len(metrics) != 3 {
		t.Fatalf("expected 3 threshold rows, got %d", len(metrics))
	}

	assertFloat(t, "first threshold", metrics[0].Threshold, 0.1)
	assertFloat(t, "second threshold", metrics[1].Threshold, 0.8)
	assertFloat(t, "third threshold", metrics[2].Threshold, 1.0)
	assertInt(t, "0.1 high edges", metrics[0].HighMsaSupportEdgeCount, 2)
	assertInt(t, "0.8 high edges", metrics[1].HighMsaSupportEdgeCount, 1)
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
	msaSupportDir := filepath.Join(dir, "msa_support")
	if err := os.MkdirAll(msaSupportDir, 0700); err != nil {
		t.Fatalf("failed to create MSA support dir: %v", err)
	}

	if err := os.WriteFile(filepath.Join(msaSupportDir, "msa_support.csv"), []byte("0,1,0\n0,0,1\n1,0,0\n"), 0644); err != nil {
		t.Fatalf("failed to write MSA support: %v", err)
	}

	solutionsPath := filepath.Join(dir, "solutions.csv")
	if err := os.WriteFile(solutionsPath, []byte("Tour,Commonality with MSA support,Min commonality with MSA,Avg commonality with MSA,Max commonality with MSA\n"), 0644); err != nil {
		t.Fatalf("failed to write solutions CSV: %v", err)
	}

	summaryPath := filepath.Join(dir, "summary.csv")
	reportPath := filepath.Join(dir, "report.md")
	_, err := AnalyzeInstances(Config{
		Instances: []InstanceConfig{
			{
				Name:                    "synthetic",
				Dimension:               3,
				MsaSupportDirectoryPath: msaSupportDir,
				OptimalToursCsvPath:     solutionsPath,
				AnalysisCsvPath:         filepath.Join(dir, "msa_support_solution_analysis.csv"),
				ThresholdsCsvPath:       filepath.Join(dir, "msa_support_solution_thresholds.csv"),
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

	thresholds, err := os.ReadFile(filepath.Join(dir, "msa_support_solution_thresholds.csv"))
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
