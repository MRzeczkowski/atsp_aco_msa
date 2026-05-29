package cycleCover

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestCalculateAnalysisStructuralMetrics(t *testing.T) {
	msaSupport := newTestMatrix(4)
	msaSupport[0][1] = 3
	msaSupport[0][2] = 3

	tours := map[string][]int{"tour": {0, 1, 2, 3}}
	cycleCoverEdges := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
		{From: 2, To: 3},
		{From: 3, To: 0},
	}

	analysis := calculateAnalysis("test", 4, msaSupport, tours, cycleCoverEdges, 1.0)
	metrics := analysis.Metrics

	assertFloat(t, "cycle-cover precision", metrics.CycleCoverMetrics.Precision, 1)
	assertFloat(t, "cycle-cover recall", metrics.CycleCoverMetrics.Recall, 1)
	assertFloat(t, "high-MSA support precision", metrics.HighMsaSupportMetrics.Precision, 0.5)
	assertFloat(t, "high-MSA support recall", metrics.HighMsaSupportMetrics.Recall, 0.25)

	if metrics.CycleCoverHighMsaEdges != 1 {
		t.Fatalf("expected one shared cycle-cover/MSA edge, got %d", metrics.CycleCoverHighMsaEdges)
	}
	if metrics.OptimalEdgesInCycleCoverAndHighMsaSupport != 1 {
		t.Fatalf("expected one optimal edge in both sets, got %d", metrics.OptimalEdgesInCycleCoverAndHighMsaSupport)
	}
	if metrics.OptimalEdgesInCycleCoverNotHighMsaSupport != 3 {
		t.Fatalf("expected three optimal edges only in cycle cover, got %d", metrics.OptimalEdgesInCycleCoverNotHighMsaSupport)
	}
}

func TestCalculateAnalysisOptimalEdgePartition(t *testing.T) {
	msaSupport := newTestMatrix(4)
	msaSupport[0][1] = 3
	msaSupport[1][2] = 3

	tours := map[string][]int{"tour": {0, 1, 2, 3}}
	cycleCoverEdges := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 3},
		{From: 2, To: 0},
		{From: 3, To: 2},
	}

	analysis := calculateAnalysis("test", 4, msaSupport, tours, cycleCoverEdges, 1.0)
	metrics := analysis.Metrics

	if metrics.OptimalEdgesInCycleCoverAndHighMsaSupport != 1 {
		t.Fatalf("expected one optimal edge in both sets, got %d", metrics.OptimalEdgesInCycleCoverAndHighMsaSupport)
	}
	if metrics.OptimalEdgesInHighMsaSupportNotCycleCover != 1 {
		t.Fatalf("expected one optimal edge only in MSA support, got %d", metrics.OptimalEdgesInHighMsaSupportNotCycleCover)
	}
	if metrics.OptimalEdgesInCycleCoverNotHighMsaSupport != 0 {
		t.Fatalf("expected no optimal edge only in cycle cover, got %d", metrics.OptimalEdgesInCycleCoverNotHighMsaSupport)
	}
	if metrics.OptimalEdgesInNeitherCycleCoverNorHigh != 2 {
		t.Fatalf("expected two optimal edges in neither set, got %d", metrics.OptimalEdgesInNeitherCycleCoverNorHigh)
	}
}

func TestEdgeSetMetricsWithoutFoundOptimalTours(t *testing.T) {
	edgeSet := map[Edge]struct{}{
		{From: 0, To: 1}: {},
		{From: 1, To: 0}: {},
	}

	metrics := calculateEdgeSetMetrics(edgeSet, nil)

	if metrics.EdgeCount != 2 {
		t.Fatalf("expected two edges, got %d", metrics.EdgeCount)
	}
	assertFloat(t, "precision", metrics.Precision, 0)
	assertFloat(t, "recall", metrics.Recall, 0)
}

func TestAnalyzeInstanceCalculatesExpectedMetrics(t *testing.T) {
	dir := t.TempDir()
	msaSupportDir := filepath.Join(dir, "msa_support")
	if err := os.MkdirAll(msaSupportDir, 0700); err != nil {
		t.Fatal(err)
	}

	writeFile(t, filepath.Join(msaSupportDir, "msa_support.csv"), strings.Join([]string{
		"0,2,0",
		"0,0,2",
		"2,0,0",
	}, "\n"))
	writeFile(t, filepath.Join(dir, "solutions.csv"), strings.Join([]string{
		"Tour,Commonality with MSA support",
		`"[0,1,2]",100`,
	}, "\n"))

	config := InstanceConfig{
		Name:                    "test",
		Dimension:               3,
		Matrix:                  [][]float64{{0, 1, 5}, {5, 0, 1}, {1, 5, 0}},
		MsaSupportDirectoryPath: msaSupportDir,
		OptimalToursCsvPath:     filepath.Join(dir, "solutions.csv"),
	}

	analysis, err := AnalyzeInstance(config, 1.0)
	if err != nil {
		t.Fatalf("AnalyzeInstance returned error: %v", err)
	}

	if analysis.Metrics.CycleCoverMetrics.EdgeCount != 3 {
		t.Fatalf("expected three cycle-cover edges, got %d", analysis.Metrics.CycleCoverMetrics.EdgeCount)
	}
	assertFloat(t, "cycle-cover precision", analysis.Metrics.CycleCoverMetrics.Precision, 1)
	assertFloat(t, "cycle-cover recall", analysis.Metrics.CycleCoverMetrics.Recall, 1)
}

func newTestMatrix(dimension int) [][]float64 {
	matrix := make([][]float64, dimension)
	for i := range matrix {
		matrix[i] = make([]float64, dimension)
	}
	return matrix
}

func writeFile(t *testing.T, path, content string) {
	t.Helper()
	if err := os.WriteFile(path, []byte(content), 0600); err != nil {
		t.Fatal(err)
	}
}

func assertFloat(t *testing.T, name string, got, want float64) {
	t.Helper()
	const tolerance = 1e-9
	if got < want-tolerance || got > want+tolerance {
		t.Fatalf("%s: expected %v, got %v", name, want, got)
	}
}
