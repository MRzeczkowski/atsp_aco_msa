package cycleCover

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestCalculateAnalysisStructuralMetrics(t *testing.T) {
	msaHeuristic := newTestMatrix(4)
	msaHeuristic[0][1] = 3
	msaHeuristic[0][2] = 3

	tours := map[string][]int{"tour": {0, 1, 2, 3}}
	cycleCoverEdges := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
		{From: 2, To: 3},
		{From: 3, To: 0},
	}

	matrix := [][]float64{
		{0, 1, 2, 3},
		{1, 0, 1, 2},
		{2, 1, 0, 1},
		{1, 2, 1, 0},
	}

	analysis := calculateAnalysis("test", 4, matrix, msaHeuristic, tours, cycleCoverEdges, 1.0)
	metrics := analysis.Metrics

	assertFloat(t, "cycle-cover precision", metrics.CycleCoverMetrics.Precision, 1)
	assertFloat(t, "cycle-cover recall", metrics.CycleCoverMetrics.Recall, 1)
	assertFloat(t, "high-MSA heuristic precision", metrics.HighMsaHeuristicMetrics.Precision, 0.5)
	assertFloat(t, "high-MSA heuristic recall", metrics.HighMsaHeuristicMetrics.Recall, 0.25)
	assertFloat(t, "cycle-cover MSA-patching precision", metrics.CycleCoverMsaPatchingMetrics.Precision, 1)
	assertFloat(t, "cycle-cover MSA-patching recall", metrics.CycleCoverMsaPatchingMetrics.Recall, 1)

	if metrics.CycleCoverHighMsaHeuristicEdges != 1 {
		t.Fatalf("expected one shared cycle-cover/MSA edge, got %d", metrics.CycleCoverHighMsaHeuristicEdges)
	}
	if metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic != 1 {
		t.Fatalf("expected one optimal edge in both sets, got %d", metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic)
	}
	if metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic != 3 {
		t.Fatalf("expected three optimal edges only in cycle cover, got %d", metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic)
	}
}

func TestCalculateAnalysisOptimalEdgePartition(t *testing.T) {
	msaHeuristic := newTestMatrix(4)
	msaHeuristic[0][1] = 3
	msaHeuristic[1][2] = 3

	tours := map[string][]int{"tour": {0, 1, 2, 3}}
	cycleCoverEdges := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 3},
		{From: 2, To: 0},
		{From: 3, To: 2},
	}

	matrix := [][]float64{
		{0, 1, 2, 3},
		{1, 0, 1, 2},
		{1, 2, 0, 1},
		{2, 3, 1, 0},
	}

	analysis := calculateAnalysis("test", 4, matrix, msaHeuristic, tours, cycleCoverEdges, 1.0)
	metrics := analysis.Metrics

	if metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic != 1 {
		t.Fatalf("expected one optimal edge in both sets, got %d", metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic)
	}
	if metrics.OptimalEdgesInHighMsaHeuristicNotCycleCover != 1 {
		t.Fatalf("expected one optimal edge only in MSA heuristic, got %d", metrics.OptimalEdgesInHighMsaHeuristicNotCycleCover)
	}
	if metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic != 0 {
		t.Fatalf("expected no optimal edge only in cycle cover, got %d", metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic)
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
	msaHeuristicDir := filepath.Join(dir, "msa_heuristic")
	if err := os.MkdirAll(msaHeuristicDir, 0700); err != nil {
		t.Fatal(err)
	}

	writeFile(t, filepath.Join(msaHeuristicDir, "msa_heuristic.csv"), strings.Join([]string{
		"0,2,0",
		"0,0,2",
		"2,0,0",
	}, "\n"))
	writeFile(t, filepath.Join(dir, "solutions.csv"), strings.Join([]string{
		"Tour,Commonality with MSA heuristic",
		`"[0,1,2]",100`,
	}, "\n"))

	config := InstanceConfig{
		Name:                      "test",
		Dimension:                 3,
		Matrix:                    [][]float64{{0, 1, 5}, {5, 0, 1}, {1, 5, 0}},
		MsaHeuristicDirectoryPath: msaHeuristicDir,
		OptimalToursCsvPath:       filepath.Join(dir, "solutions.csv"),
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
