package reports

import (
	"atsp_aco_msa/modules/analysis/structure"
	"atsp_aco_msa/modules/models"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

func TestSaveStructuralSimilarityReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "structural_similarity.md")
	if err := SaveStructuralSimilarity(path, sampleStructuralAnalyses()); err != nil {
		t.Fatalf("SaveStructuralSimilarity returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read structural similarity report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "- **Precision vs found-optimal tours: MSA heuristic 71.43%, cycle cover 66.67%, cycle-cover MSA patching 63.64%.**")
	assertContains(t, content, "- **Recall vs found-optimal tours: MSA heuristic 35.71%, cycle cover 42.86%, cycle-cover MSA patching 50.00%.**")
	assertContains(t, content, "<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">MSA heuristic</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>")
	assertContains(t, content, "<tr><td>a</td><td align=\"right\">50.00</td><td align=\"right\">25.00</td><td align=\"right\"><strong>75.00</strong></td><td align=\"right\"><strong>75.00</strong></td><td align=\"right\">60.00</td><td align=\"right\"><strong>75.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Total</strong></td><td align=\"right\"><strong>71.43</strong></td><td align=\"right\">35.71</td><td align=\"right\">66.67</td><td align=\"right\">42.86</td><td align=\"right\">63.64</td><td align=\"right\"><strong>50.00</strong></td></tr>")
}

func TestSaveMsaHeuristicCycleCoverOverlapReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "msa_cycle_cover_overlap.md")
	if err := SaveMsaHeuristicCycleCoverOverlap(path, sampleStructuralAnalyses()); err != nil {
		t.Fatalf("SaveMsaHeuristicCycleCoverOverlap returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read MSA heuristic/cycle-cover overlap report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "- **42.86% of MSA heuristic edges are also cycle-cover edges.**")
	assertContains(t, content, "- **33.33% of cycle-cover edges are also MSA heuristic edges.**")
	assertContains(t, content, "- **Found-optimal edge partition: both 3, only MSA heuristic 2, only cycle cover 3, neither 6.**")
	assertContains(t, content, "<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>")
	assertContains(t, content, "<tr><td>a</td><td align=\"right\">50.00</td><td align=\"right\">25.00</td><td align=\"right\">1</td><td align=\"right\">0</td><td align=\"right\">2</td></tr>")
	assertContains(t, content, "<tr><td><strong>Total</strong></td><td align=\"right\">42.86</td><td align=\"right\">33.33</td><td align=\"right\">3</td><td align=\"right\">2</td><td align=\"right\">3</td></tr>")
}

func TestTourMatrixLengthValidatesSingleTourAndCalculatesLength(t *testing.T) {
	distanceMatrix := [][]float64{
		{0, 2, 8},
		{4, 0, 3},
		{5, 6, 0},
	}
	tourMatrix := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 0},
	}

	length, err := tourMatrixLength(distanceMatrix, tourMatrix)
	if err != nil {
		t.Fatalf("tourMatrixLength returned unexpected error: %v", err)
	}
	if length != 10 {
		t.Fatalf("expected tour length 10, got %.2f", length)
	}

	disconnected := [][]float64{
		{0, 1, 0},
		{1, 0, 0},
		{0, 0, 0},
	}
	if _, err := tourMatrixLength(distanceMatrix, disconnected); err == nil {
		t.Fatal("expected invalid tour matrix to be rejected")
	}
}

func TestSaveGksDeviationReportShowsMsaPatchBiasDeviation(t *testing.T) {
	path := filepath.Join(t.TempDir(), "gks_deviation.md")
	rows := []gksDeviationRow{
		{instance: "sample-a", msaPatchBias: 0.0, tourLength: 110, deviation: 10},
		{instance: "sample-a", msaPatchBias: 0.5, tourLength: 105, deviation: 5},
		{instance: "sample-b", msaPatchBias: 0.0, tourLength: 240, deviation: 20},
		{instance: "sample-b", msaPatchBias: 0.5, tourLength: 260, deviation: 30},
	}

	var builder strings.Builder
	builder.WriteString("# GKS Patching Deviation\n\n")
	writeGksDeviationTable(&builder, rows, []float64{0.0, 0.5})
	if err := os.WriteFile(path, []byte(builder.String()), 0644); err != nil {
		t.Fatalf("failed to write test GKS report: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read GKS report: %v", err)
	}
	content := string(contentBytes)

	assertContains(t, content, "<tr><th>Instance</th><th>MSA patch bias</th><th>Tour length</th><th>Deviation [%]</th></tr>")
	assertContains(t, content, "<tr><td>sample-a</td><td align=\"right\">0.50</td><td align=\"right\">105.00</td><td align=\"right\"><strong>5.00</strong></td></tr>")
	assertContains(t, content, "<tr><td>sample-b</td><td align=\"right\">0.00</td><td align=\"right\">240.00</td><td align=\"right\"><strong>20.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Average</strong></td><td align=\"right\">0.00</td><td align=\"right\"></td><td align=\"right\"><strong>15.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Average</strong></td><td align=\"right\">0.50</td><td align=\"right\"></td><td align=\"right\">17.50</td></tr>")
}

func TestSelectEvenlySpacedRootIndexes(t *testing.T) {
	actual := selectEvenlySpacedRootIndexes(10, 4)
	expected := []int{0, 3, 6, 9}
	if !reflect.DeepEqual(actual, expected) {
		t.Fatalf("unexpected root selection\nwant: %v\n got: %v", expected, actual)
	}

	actual = selectEvenlySpacedRootIndexes(3, 64)
	expected = []int{0, 1, 2}
	if !reflect.DeepEqual(actual, expected) {
		t.Fatalf("unexpected capped root selection\nwant: %v\n got: %v", expected, actual)
	}
}

func TestBuildPartialMsaHeuristicEdgeSetUsesEligibilityByRoot(t *testing.T) {
	msas := [][][]float64{
		{
			{0, 1, 0, 0},
			{0, 0, 1, 0},
			{0, 0, 0, 1},
			{0, 0, 0, 0},
		},
		emptyMatrix(4),
		{
			{0, 1, 0, 0},
			{0, 0, 0, 1},
			{0, 0, 0, 0},
			{1, 0, 0, 0},
		},
		emptyMatrix(4),
	}

	actual := buildPartialMsaHeuristicEdgeSet(msas, []int{0, 2})
	expected := map[models.Edge]struct{}{
		models.Edge{From: 0, To: 1}: {},
		models.Edge{From: 1, To: 2}: {},
		models.Edge{From: 3, To: 0}: {},
	}
	if !reflect.DeepEqual(actual, expected) {
		t.Fatalf("unexpected partial MSA heuristic edge set\nwant: %v\n got: %v", expected, actual)
	}
}

func sampleStructuralAnalyses() []structure.InstanceAnalysis {
	return []structure.InstanceAnalysis{
		{
			Instance:  "b",
			Dimension: 5,
			Metrics: structure.InstanceMetrics{
				FoundOptimalTourCount:       2,
				UniqueFoundOptimalEdgeCount: 10,
				HighMsaHeuristicMetrics: structure.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 4,
					Precision:        0.8,
					Recall:           0.4,
				},
				CycleCoverMetrics: structure.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 3,
					Precision:        0.6,
					Recall:           0.3,
				},
				CycleCoverMsaPatchingMetrics: structure.EdgeSetMetrics{
					EdgeCount:        6,
					OptimalEdgeCount: 4,
					Precision:        4.0 / 6.0,
					Recall:           0.4,
				},
				CycleCoverHighMsaHeuristicEdges:             2,
				OptimalEdgesInCycleCoverAndHighMsaHeuristic: 2,
				OptimalEdgesInHighMsaHeuristicNotCycleCover: 2,
				OptimalEdgesInCycleCoverNotHighMsaHeuristic: 1,
				OptimalEdgesInNeitherCycleCoverNorHigh:      5,
			},
		},
		{
			Instance:  "a",
			Dimension: 4,
			Metrics: structure.InstanceMetrics{
				FoundOptimalTourCount:       1,
				UniqueFoundOptimalEdgeCount: 4,
				HighMsaHeuristicMetrics: structure.EdgeSetMetrics{
					EdgeCount:        2,
					OptimalEdgeCount: 1,
					Precision:        0.5,
					Recall:           0.25,
				},
				CycleCoverMetrics: structure.EdgeSetMetrics{
					EdgeCount:        4,
					OptimalEdgeCount: 3,
					Precision:        0.75,
					Recall:           0.75,
				},
				CycleCoverMsaPatchingMetrics: structure.EdgeSetMetrics{
					EdgeCount:        5,
					OptimalEdgeCount: 3,
					Precision:        0.6,
					Recall:           0.75,
				},
				CycleCoverHighMsaHeuristicEdges:             1,
				OptimalEdgesInCycleCoverAndHighMsaHeuristic: 1,
				OptimalEdgesInHighMsaHeuristicNotCycleCover: 0,
				OptimalEdgesInCycleCoverNotHighMsaHeuristic: 2,
				OptimalEdgesInNeitherCycleCoverNorHigh:      1,
			},
		},
	}
}

func emptyMatrix(size int) [][]float64 {
	matrix := make([][]float64, size)
	for i := range matrix {
		matrix[i] = make([]float64, size)
	}

	return matrix
}

func assertContains(t *testing.T, content, expected string) {
	t.Helper()
	if !strings.Contains(content, expected) {
		t.Fatalf("expected content to contain:\n%s\n\nactual content:\n%s", expected, content)
	}
}
