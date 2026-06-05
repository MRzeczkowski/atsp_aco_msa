package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestSaveTuningSummaryRanksParametersAcrossInstances(t *testing.T) {
	root := t.TempDir()
	first := makeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, root)
	second := makeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, root)

	saveStatistics(resultFilePathForHeuristic(first, heuristicMsaHeuristic), heuristicMsaHeuristic, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.50, 2.0, 10.0),
		makeTestExperimentStatistics(0.90, 1.0, 80.0),
	})
	saveStatistics(resultFilePathForHeuristic(second, heuristicMsaHeuristic), heuristicMsaHeuristic, []ExperimentsDataStatistics{
		makeTestExperimentStatistics(0.50, 1.0, 20.0),
		makeTestExperimentStatistics(0.90, 1.2, 20.0),
	})

	if err := saveTuningSummary(root, []AtspData{second, first}, []string{heuristicMsaHeuristic}); err != nil {
		t.Fatalf("saveTuningSummary returned unexpected error: %v", err)
	}

	content := readTuningSummary(t, root)
	assertContains(t, content, "# Tuning Summary")
	assertContains(t, content, "## MSA heuristic")
	assertContains(t, content, "Instances used: 2/2")
	assertContains(t, content, "| 0.90 | 1.10 | 1.10 | 1 | 2 | 50.00 |")
	assertContains(t, content, "| 0.50 | 1.50 | 1.50 | 1 | 1 | 15.00 |")

	bestRow := strings.Index(content, "| 0.90 | 1.10 | 1.10 | 1 | 2 | 50.00 |")
	secondRow := strings.Index(content, "| 0.50 | 1.50 | 1.50 | 1 | 1 | 15.00 |")
	if bestRow < 0 || secondRow < 0 || bestRow > secondRow {
		t.Fatalf("expected rows to be sorted by mean result, got:\n%s", content)
	}
}

func TestSaveTuningSummaryIncludesPatchingMsaPatchBias(t *testing.T) {
	root := t.TempDir()
	atspData := makeAtspDataInResultsDirectory("sample.atsp", [][]float64{{0, 1}, {1, 0}}, 2, root)

	lowerBias := makeTestExperimentStatistics(0.70, 1.0, 60.0)
	lowerBias.msaPatchBias = 0.25
	higherBias := makeTestExperimentStatistics(0.70, 2.0, 20.0)
	higherBias.msaPatchBias = 0.75
	saveStatistics(resultFilePathForHeuristic(atspData, heuristicCycleCoverMsaPatching), heuristicCycleCoverMsaPatching, []ExperimentsDataStatistics{
		lowerBias,
		higherBias,
	})

	if err := saveTuningSummary(root, []AtspData{atspData}, []string{heuristicCycleCoverMsaPatching}); err != nil {
		t.Fatalf("saveTuningSummary returned unexpected error: %v", err)
	}

	content := readTuningSummary(t, root)
	assertContains(t, content, "## Cycle-cover MSA patching")
	assertContains(t, content, "| Heuristic weight | MSA patch bias | Mean result [%] | Median result [%] | Exact best | Near best | Mean success [%] |")
	assertContains(t, content, "| 0.70 | 0.25 | 1.00 | 1.00 | 1 | 1 | 60.00 |")
	assertContains(t, content, "| 0.70 | 0.75 | 2.00 | 2.00 | 0 | 0 | 20.00 |")
}

func TestSaveTuningSummaryRemovesLegacyBestParametersReports(t *testing.T) {
	root := t.TempDir()
	atspData := makeAtspDataInResultsDirectory("missing.atsp", [][]float64{{0, 1}, {1, 0}}, 2, root)
	legacyPath := filepath.Join(root, "best_parameters_report.md")
	if err := os.WriteFile(legacyPath, []byte("legacy"), 0644); err != nil {
		t.Fatalf("failed to write legacy report: %v", err)
	}

	if err := saveTuningSummary(root, []AtspData{atspData}, []string{heuristicMsaHeuristic}); err != nil {
		t.Fatalf("saveTuningSummary returned unexpected error: %v", err)
	}

	if _, err := os.Stat(legacyPath); !os.IsNotExist(err) {
		t.Fatalf("expected legacy best-parameters report to be removed, stat err=%v", err)
	}
	content := readTuningSummary(t, root)
	assertContains(t, content, "Instances used: 0/1")
	assertContains(t, content, "Missing result files: missing")
}

func readTuningSummary(t *testing.T, root string) string {
	t.Helper()
	content, err := os.ReadFile(filepath.Join(root, tuningSummaryReportFileName))
	if err != nil {
		t.Fatalf("failed to read tuning summary: %v", err)
	}
	return string(content)
}
