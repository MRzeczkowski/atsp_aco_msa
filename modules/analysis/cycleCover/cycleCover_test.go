package cycleCover

import (
	"encoding/csv"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

func TestCalculateAnalysisCycleCoverGateFiltersFalsePositive(t *testing.T) {
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

	analysis := calculateAnalysis("test", 4, 4, msaSupport, tours, cycleCoverEdges, 4, []float64{0.8}, 1.0)
	metrics := analysis.Metrics

	assertFloat(t, "cycle-cover precision", metrics.CycleCoverMetrics.Precision, 1)
	assertFloat(t, "cycle-cover recall", metrics.CycleCoverMetrics.Recall, 1)
	assertFloat(t, "high-MSA support precision", metrics.HighMsaSupportMetrics.Precision, 0.5)
	assertFloat(t, "high-MSA support recall", metrics.HighMsaSupportMetrics.Recall, 0.25)
	assertFloat(t, "gated high-MSA support precision", metrics.CycleCoverHighMsaSupportMetrics.Precision, 1)
	assertFloat(t, "gated high-MSA support recall", metrics.CycleCoverHighMsaSupportMetrics.Recall, 0.25)
	assertFloat(t, "cycle-cover positive-MSA support share", metrics.CycleCoverPositiveMsaSupportShare, 0.25)
	assertFloat(t, "cycle-cover high-MSA support share", metrics.CycleCoverHighMsaSupportShare, 0.25)
	assertFloat(t, "precision gain", metrics.HighMsaSupportPrecisionGainFromCycleCover, 0.5)
	assertFloat(t, "recall loss", metrics.HighMsaSupportRecallLossFromCycleCover, 0)

	if metrics.CycleCoverEdgesWithPositiveMsaSupport != 1 {
		t.Fatalf("expected one cycle-cover edge with positive MSA support, got %d", metrics.CycleCoverEdgesWithPositiveMsaSupport)
	}
	if metrics.HighMsaSupportEdgesRemovedByCycleCover != 1 {
		t.Fatalf("expected one high-MSA support edge to be removed, got %d", metrics.HighMsaSupportEdgesRemovedByCycleCover)
	}
	if metrics.HighMsaSupportOptimalEdgesRemovedByCycleCover != 0 {
		t.Fatalf("expected no optimal high-MSA support edge to be removed, got %d", metrics.HighMsaSupportOptimalEdgesRemovedByCycleCover)
	}
}

func TestCalculateAnalysisCycleCoverGateCanLoseRecall(t *testing.T) {
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

	analysis := calculateAnalysis("test", 4, 4, msaSupport, tours, cycleCoverEdges, 4, []float64{1.0}, 1.0)
	metrics := analysis.Metrics

	assertFloat(t, "high-MSA support recall", metrics.HighMsaSupportMetrics.Recall, 0.5)
	assertFloat(t, "gated high-MSA support recall", metrics.CycleCoverHighMsaSupportMetrics.Recall, 0.25)
	assertFloat(t, "recall loss", metrics.HighMsaSupportRecallLossFromCycleCover, 0.25)

	if metrics.HighMsaSupportOptimalEdgesRemovedByCycleCover != 1 {
		t.Fatalf("expected one optimal high-MSA support edge to be removed, got %d", metrics.HighMsaSupportOptimalEdgesRemovedByCycleCover)
	}
}

func TestThresholdMetricsAreSorted(t *testing.T) {
	msaSupport := newTestMatrix(4)
	msaSupport[0][1] = 3
	msaSupport[1][2] = 2
	msaSupport[2][3] = 1

	optimalEdges := buildTourEdgeSet(map[string][]int{"tour": {0, 1, 2, 3}})
	cycleCoverSet := buildEdgeSet([]Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
		{From: 2, To: 3},
		{From: 3, To: 0},
	})

	metrics := calculateThresholdMetrics(msaSupport, optimalEdges, cycleCoverSet, []float64{0.8, 0.1, 1.0}, 12)

	got := []float64{metrics[0].Threshold, metrics[1].Threshold, metrics[2].Threshold}
	want := []float64{0.1, 0.8, 1.0}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("expected sorted thresholds %v, got %v", want, got)
	}
}

func TestMetricsWithoutFoundOptimalToursDoNotReportFalsePositives(t *testing.T) {
	edgeSet := map[Edge]struct{}{
		{From: 0, To: 1}: {},
		{From: 1, To: 0}: {},
	}

	metrics := calculateEdgeSetMetrics(edgeSet, nil, 2)

	if metrics.FalsePositiveEdges != 0 {
		t.Fatalf("expected no false-positive count without found optimal tours, got %d", metrics.FalsePositiveEdges)
	}
	assertFloat(t, "precision", metrics.Precision, 0)
	assertFloat(t, "recall", metrics.Recall, 0)
	assertFloat(t, "lift", metrics.Lift, 0)
}

func TestAnalyzeInstanceWritesExpectedFiles(t *testing.T) {
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
		Name:                        "test",
		Dimension:                   3,
		Matrix:                      [][]float64{{0, 1, 5}, {5, 0, 1}, {1, 5, 0}},
		KnownOptimal:                3,
		MsaSupportDirectoryPath:     msaSupportDir,
		OptimalToursCsvPath:         filepath.Join(dir, "solutions.csv"),
		CycleCoverEdgesCsvPath:      filepath.Join(dir, "cycle_cover_edges.csv"),
		AnalysisCsvPath:             filepath.Join(dir, "cycle_cover_analysis.csv"),
		ThresholdsCsvPath:           filepath.Join(dir, "cycle_cover_thresholds.csv"),
		CycleCoverOverlapMatrixPath: filepath.Join(dir, "cycle_cover_msa_support_overlap.csv"),
	}

	analysis, err := AnalyzeInstance(config, []float64{1.0, 0.5}, 1.0)
	if err != nil {
		t.Fatalf("AnalyzeInstance returned error: %v", err)
	}

	if analysis.Metrics.CycleCoverCost != 3 {
		t.Fatalf("expected cycle-cover cost 3, got %v", analysis.Metrics.CycleCoverCost)
	}

	assertCsvRowCount(t, config.CycleCoverEdgesCsvPath, 4)
	assertCsvRowCount(t, config.AnalysisCsvPath, 53)
	assertCsvRowCount(t, config.ThresholdsCsvPath, 3)
	assertCsvRowCount(t, config.CycleCoverOverlapMatrixPath, 3)
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

func assertCsvRowCount(t *testing.T, path string, want int) {
	t.Helper()

	file, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer file.Close()

	rows, err := csv.NewReader(file).ReadAll()
	if err != nil {
		t.Fatal(err)
	}
	if len(rows) != want {
		t.Fatalf("%s: expected %d rows, got %d", path, want, len(rows))
	}
}

func assertFloat(t *testing.T, name string, got, want float64) {
	t.Helper()
	const tolerance = 1e-9
	if got < want-tolerance || got > want+tolerance {
		t.Fatalf("%s: expected %v, got %v", name, want, got)
	}
}
