package tuning

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestSaveRanksParametersAcrossInstances(t *testing.T) {
	root := t.TempDir()
	statisticsByPath := map[string][]Statistic{
		"sample-a.csv": {
			{HeuristicWeight: 0.50, AverageBestDeviation: 2.0, SuccessRate: 10.0},
			{HeuristicWeight: 0.90, AverageBestDeviation: 1.0, SuccessRate: 80.0},
		},
		"sample-b.csv": {
			{HeuristicWeight: 0.50, AverageBestDeviation: 1.0, SuccessRate: 20.0},
			{HeuristicWeight: 0.90, AverageBestDeviation: 1.2, SuccessRate: 20.0},
		},
	}

	err := Save(Config{
		ResultsRootPath: root,
		Heuristics: []HeuristicConfig{
			{
				Name:        "msa-heuristic",
				DisplayName: "MSA heuristic",
				ResultFiles: []ResultFile{
					{Instance: "sample-b", Path: "sample-b.csv"},
					{Instance: "sample-a", Path: "sample-a.csv"},
				},
			},
		},
		ReadStatistics: fakeReader(statisticsByPath),
	})
	if err != nil {
		t.Fatalf("Save returned unexpected error: %v", err)
	}

	content := readSummary(t, root)
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

func TestSaveIncludesPatchingMsaPatchBias(t *testing.T) {
	root := t.TempDir()
	err := Save(Config{
		ResultsRootPath: root,
		Heuristics: []HeuristicConfig{
			{
				Name:                "cycle-cover-msa-patching",
				DisplayName:         "Cycle-cover MSA patching",
				IncludeMsaPatchBias: true,
				ResultFiles: []ResultFile{
					{Instance: "sample", Path: "sample.csv"},
				},
			},
		},
		ReadStatistics: fakeReader(map[string][]Statistic{
			"sample.csv": {
				{HeuristicWeight: 0.70, MsaPatchBias: 0.25, AverageBestDeviation: 1.0, SuccessRate: 60.0},
				{HeuristicWeight: 0.70, MsaPatchBias: 0.75, AverageBestDeviation: 2.0, SuccessRate: 20.0},
			},
		}),
	})
	if err != nil {
		t.Fatalf("Save returned unexpected error: %v", err)
	}

	content := readSummary(t, root)
	assertContains(t, content, "## Cycle-cover MSA patching")
	assertContains(t, content, "| Heuristic weight | MSA patch bias | Mean result [%] | Median result [%] | Exact best | Near best | Mean success [%] |")
	assertContains(t, content, "| 0.70 | 0.25 | 1.00 | 1.00 | 1 | 1 | 60.00 |")
	assertContains(t, content, "| 0.70 | 0.75 | 2.00 | 2.00 | 0 | 0 | 20.00 |")
}

func TestSaveRemovesLegacyBestParametersReports(t *testing.T) {
	root := t.TempDir()
	legacyPath := filepath.Join(root, "best_parameters_report.md")
	if err := os.WriteFile(legacyPath, []byte("legacy"), 0644); err != nil {
		t.Fatalf("failed to write legacy report: %v", err)
	}

	err := Save(Config{
		ResultsRootPath: root,
		Heuristics: []HeuristicConfig{
			{
				Name:        "msa-heuristic",
				DisplayName: "MSA heuristic",
				ResultFiles: []ResultFile{
					{Instance: "missing", Path: "missing.csv"},
				},
			},
		},
		ReadStatistics: fakeReader(nil),
	})
	if err != nil {
		t.Fatalf("Save returned unexpected error: %v", err)
	}

	if _, err := os.Stat(legacyPath); !os.IsNotExist(err) {
		t.Fatalf("expected legacy best-parameters report to be removed, stat err=%v", err)
	}
	content := readSummary(t, root)
	assertContains(t, content, "Instances used: 0/1")
	assertContains(t, content, "Missing result files: missing")
}

func fakeReader(statisticsByPath map[string][]Statistic) func(path string) ([]Statistic, error) {
	return func(path string) ([]Statistic, error) {
		statistics, ok := statisticsByPath[path]
		if !ok {
			return nil, os.ErrNotExist
		}

		return statistics, nil
	}
}

func readSummary(t *testing.T, root string) string {
	t.Helper()
	content, err := os.ReadFile(filepath.Join(root, ReportFileName))
	if err != nil {
		t.Fatalf("failed to read tuning summary: %v", err)
	}
	return string(content)
}

func assertContains(t *testing.T, content, expected string) {
	t.Helper()
	if !strings.Contains(content, expected) {
		t.Fatalf("expected content to contain:\n%s\n\ncontent:\n%s", expected, content)
	}
}
