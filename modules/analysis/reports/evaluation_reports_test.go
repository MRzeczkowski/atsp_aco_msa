package reports

import (
	"atsp_aco_msa/modules/experiments"
	"atsp_aco_msa/modules/project"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

const (
	testHeuristicBaseline              = "baseline"
	testHeuristicStrictMsa             = "strict-msa"
	testHeuristicRootedMsa             = "rooted-msa"
	testHeuristicCycleCover            = "cycle-cover"
	testHeuristicCycleCoverMsaPatching = "cycle-cover-msa-patching"
)

func TestSaveEvaluationResultsSummaryWritesMarkdownTableWithHighlightedFindings(t *testing.T) {
	resultsRoot := t.TempDir()
	firstAtspData := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	secondAtspData := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, resultsRoot)
	rows := []experiments.HeuristicExperimentStatistics{
		{
			Heuristic: testHeuristicBaseline,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 4.25,
				SuccessRate:          10.0,
			},
		},
		{
			Heuristic: testHeuristicStrictMsa,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 2.50,
				SuccessRate:          20.0,
			},
		},
		{
			Heuristic: testHeuristicRootedMsa,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 3.00,
				SuccessRate:          15.0,
			},
		},
		{
			Heuristic: testHeuristicCycleCover,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 1.75,
				SuccessRate:          30.0,
			},
		},
		{
			Heuristic: testHeuristicCycleCoverMsaPatching,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 2.00,
				SuccessRate:          25.0,
			},
		},
	}
	secondRows := []experiments.HeuristicExperimentStatistics{
		{
			Heuristic: testHeuristicBaseline,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 6.75,
				SuccessRate:          20.0,
			},
		},
		{
			Heuristic: testHeuristicStrictMsa,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 3.50,
				SuccessRate:          40.0,
			},
		},
		{
			Heuristic: testHeuristicRootedMsa,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 4.00,
				SuccessRate:          35.0,
			},
		},
		{
			Heuristic: testHeuristicCycleCover,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 2.25,
				SuccessRate:          60.0,
			},
		},
		{
			Heuristic: testHeuristicCycleCoverMsaPatching,
			Statistics: experiments.ExperimentsDataStatistics{
				AverageBestDeviation: 2.50,
				SuccessRate:          55.0,
			},
		},
	}

	if err := experiments.SaveHeuristicStatistics(firstAtspData.ResultFilePath, rows); err != nil {
		t.Fatalf("SaveHeuristicStatistics returned unexpected error: %v", err)
	}
	if err := experiments.SaveHeuristicStatistics(secondAtspData.ResultFilePath, secondRows); err != nil {
		t.Fatalf("SaveHeuristicStatistics returned unexpected error: %v", err)
	}

	summaryPath := filepath.Join(resultsRoot, "summary.md")
	summaryRows, err := ReadEvaluationResultsSummaryRows([]project.AtspData{firstAtspData, secondAtspData})
	if err != nil {
		t.Fatalf("ReadEvaluationResultsSummaryRows returned unexpected error: %v", err)
	}
	if err := SaveEvaluationResultsSummaryRows(summaryRows, summaryPath, evaluationReportsTestConfig()); err != nil {
		t.Fatalf("SaveEvaluationResultsSummaryRows returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(summaryPath)
	if err != nil {
		t.Fatalf("failed to read evaluation summary Markdown: %v", err)
	}

	lines := strings.Split(strings.TrimSpace(string(contentBytes)), "\n")
	expected := []string{
		"# Evaluation Results Summary",
		"",
		"## Findings",
		"",
		"- **Cycle cover has the lowest average best deviation overall: 2.00%.**",
		"- **Cycle cover has the highest average success rate overall: 45.00%.**",
		"- **Best-or-tied average best deviation counts: Baseline 0/2, Strict MSA 0/2, Rooted MSA 0/2, Cycle cover 2/2, Cycle-cover MSA patching 0/2.**",
		"",
		"<table>",
		"<thead>",
		"<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">Baseline</th><th colspan=\"2\">Strict MSA</th><th colspan=\"2\">Rooted MSA</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>",
		"<tr><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th><th>Avg best dev. [%]</th><th>Success [%]</th></tr>",
		"</thead>",
		"<tbody>",
		"<tr><td>sample-a</td><td align=\"right\">4.25</td><td align=\"right\">10.00</td><td align=\"right\">2.50</td><td align=\"right\">20.00</td><td align=\"right\">3.00</td><td align=\"right\">15.00</td><td align=\"right\"><strong>1.75</strong></td><td align=\"right\"><strong>30.00</strong></td><td align=\"right\">2.00</td><td align=\"right\">25.00</td></tr>",
		"<tr><td>sample-b</td><td align=\"right\">6.75</td><td align=\"right\">20.00</td><td align=\"right\">3.50</td><td align=\"right\">40.00</td><td align=\"right\">4.00</td><td align=\"right\">35.00</td><td align=\"right\"><strong>2.25</strong></td><td align=\"right\"><strong>60.00</strong></td><td align=\"right\">2.50</td><td align=\"right\">55.00</td></tr>",
		"<tr><td><strong>Average</strong></td><td align=\"right\">5.50</td><td align=\"right\">15.00</td><td align=\"right\">3.00</td><td align=\"right\">30.00</td><td align=\"right\">3.50</td><td align=\"right\">25.00</td><td align=\"right\"><strong>2.00</strong></td><td align=\"right\"><strong>45.00</strong></td><td align=\"right\">2.25</td><td align=\"right\">40.00</td></tr>",
		"</tbody>",
		"</table>",
	}
	if !reflect.DeepEqual(lines, expected) {
		t.Fatalf("unexpected evaluation summary Markdown\nwant: %v\n got: %v", expected, lines)
	}
}

func TestSaveEvaluationThreeOptComparisonReportShowsHiddenHeuristicEffect(t *testing.T) {
	path := filepath.Join(t.TempDir(), "comparison_to_evaluation.md")
	evaluationRows := []EvaluationResultsSummaryRow{
		{
			Instance: "sample",
			Metrics: map[string]EvaluationResultSummaryMetric{
				testHeuristicBaseline: {
					AverageMinDeviation:  10.0,
					SuccessRate:          10.0,
					AverageBestIteration: 50.0,
					Iterations:           100,
				},
				testHeuristicStrictMsa: {
					AverageMinDeviation:  8.0,
					SuccessRate:          12.0,
					AverageBestIteration: 60.0,
					Iterations:           100,
				},
				testHeuristicCycleCover: {
					AverageMinDeviation:  7.0,
					SuccessRate:          11.0,
					AverageBestIteration: 30.0,
					Iterations:           100,
				},
				testHeuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  7.5,
					SuccessRate:          13.0,
					AverageBestIteration: 35.0,
					Iterations:           100,
				},
			},
		},
	}
	evaluationThreeOptRows := []EvaluationResultsSummaryRow{
		{
			Instance: "sample",
			Metrics: map[string]EvaluationResultSummaryMetric{
				testHeuristicBaseline: {
					AverageMinDeviation:  1.0,
					SuccessRate:          50.0,
					AverageBestIteration: 20.0,
					Iterations:           100,
				},
				testHeuristicStrictMsa: {
					AverageMinDeviation:  0.9,
					SuccessRate:          51.0,
					AverageBestIteration: 19.0,
					Iterations:           100,
				},
				testHeuristicCycleCover: {
					AverageMinDeviation:  0.8,
					SuccessRate:          52.0,
					AverageBestIteration: 19.0,
					Iterations:           100,
				},
				testHeuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  0.85,
					SuccessRate:          53.0,
					AverageBestIteration: 18.0,
					Iterations:           100,
				},
			},
		},
	}

	if err := SaveEvaluationThreeOptComparisonReport(path, evaluationRows, evaluationThreeOptRows, evaluationReportsTestConfig()); err != nil {
		t.Fatalf("SaveEvaluationThreeOptComparisonReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read comparison report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "# Reduced 3-Opt Impact")
	assertContains(t, content, "Cycle-cover MSA patching +2.50 -> +0.15 pp")
	assertContains(t, content, "<tr><th>Heuristic</th><th>Avg best dev. without 3-opt [%]</th><th>Avg best dev. with 3-opt [%]</th><th>Success without 3-opt [%]</th><th>Success with 3-opt [%]</th></tr>")
	assertContains(t, content, "<tr><td>Strict MSA</td><td align=\"right\">+2.00</td><td align=\"right\">+0.10</td><td align=\"right\">5.00</td></tr>")
	assertContains(t, content, "<tr><td>Cycle cover</td><td align=\"right\">+3.00</td><td align=\"right\">+0.20</td><td align=\"right\">6.67</td></tr>")
	assertContains(t, content, "<tr><td>Cycle-cover MSA patching</td><td align=\"right\">+2.50</td><td align=\"right\">+0.15</td><td align=\"right\">6.00</td></tr>")
}

func TestSaveEvaluationPairwisePerformanceReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "pairwise_performance.md")
	if err := SaveEvaluationPairwisePerformanceReport(path, sampleEvaluationSummaryRows(), evaluationReportsTestConfig()); err != nil {
		t.Fatalf("SaveEvaluationPairwisePerformanceReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read pairwise performance report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "<tr><td>Strict MSA vs Baseline</td><td align=\"right\">-1.50</td><td align=\"right\">2</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">+15.00</td></tr>")
	assertContains(t, content, "<tr><td>Cycle-cover MSA patching vs Baseline</td><td align=\"right\">-2.10</td><td align=\"right\">2</td><td align=\"right\">0</td><td align=\"right\">0</td><td align=\"right\">+20.00</td></tr>")
	assertContains(t, content, "<tr><td>Strict MSA vs Cycle cover</td><td align=\"right\">+0.25</td><td align=\"right\">0</td><td align=\"right\">1</td><td align=\"right\">1</td><td align=\"right\">-5.00</td></tr>")
}

func TestSaveEvaluationConvergenceSummaryReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "convergence_summary.md")
	if err := SaveEvaluationConvergenceSummaryReport(path, sampleEvaluationSummaryRows(), evaluationReportsTestConfig()); err != nil {
		t.Fatalf("SaveEvaluationConvergenceSummaryReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read convergence summary report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "- **Average best-iteration position: Baseline 75.00%, Strict MSA 50.00%, Cycle cover 35.00%, Cycle-cover MSA patching 27.50%.**")
	assertContains(t, content, "<tr><td>a</td><td align=\"right\">80.00</td><td align=\"right\">60.00</td><td></td><td align=\"right\">40.00</td><td align=\"right\"><strong>35.00</strong></td></tr>")
	assertContains(t, content, "<tr><td><strong>Average</strong></td><td align=\"right\">75.00</td><td align=\"right\">50.00</td><td></td><td align=\"right\">35.00</td><td align=\"right\"><strong>27.50</strong></td></tr>")
}

func TestSaveStructuralPerformanceLinkReport(t *testing.T) {
	path := filepath.Join(t.TempDir(), "structural_performance_link.md")
	if err := SaveStructuralPerformanceLinkReport(path, sampleEvaluationSummaryRows(), sampleStructuralAnalyses(), evaluationReportsTestConfig()); err != nil {
		t.Fatalf("SaveStructuralPerformanceLinkReport returned unexpected error: %v", err)
	}

	contentBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("failed to read structural performance link report: %v", err)
	}

	content := string(contentBytes)
	assertContains(t, content, "<tr><td>Strict MSA</td><td align=\"right\"><strong>71.43</strong></td><td align=\"right\">35.71</td><td align=\"right\">3.00</td><td align=\"right\">20.00</td></tr>")
	assertContains(t, content, "<tr><td>Cycle cover</td><td align=\"right\">66.67</td><td align=\"right\">42.86</td><td align=\"right\">2.75</td><td align=\"right\"><strong>25.00</strong></td></tr>")
	assertContains(t, content, "<tr><td>Cycle-cover MSA patching</td><td align=\"right\">63.64</td><td align=\"right\"><strong>50.00</strong></td><td align=\"right\"><strong>2.40</strong></td><td align=\"right\"><strong>25.00</strong></td></tr>")
}

func evaluationReportsTestConfig() EvaluationReportsConfig {
	return EvaluationReportsConfig{
		Heuristics: []string{
			testHeuristicBaseline,
			testHeuristicStrictMsa,
			testHeuristicRootedMsa,
			testHeuristicCycleCover,
			testHeuristicCycleCoverMsaPatching,
		},
		BaselineHeuristic:              testHeuristicBaseline,
		StrictMsaHeuristic:             testHeuristicStrictMsa,
		CycleCoverHeuristic:            testHeuristicCycleCover,
		CycleCoverMsaPatchingHeuristic: testHeuristicCycleCoverMsaPatching,
		DisplayName:                    evaluationReportTestDisplayName,
	}
}

func evaluationReportTestDisplayName(heuristic string) string {
	switch heuristic {
	case testHeuristicBaseline:
		return "Baseline"
	case testHeuristicStrictMsa:
		return "Strict MSA"
	case testHeuristicRootedMsa:
		return "Rooted MSA"
	case testHeuristicCycleCover:
		return "Cycle cover"
	case testHeuristicCycleCoverMsaPatching:
		return "Cycle-cover MSA patching"
	default:
		return heuristic
	}
}

func sampleEvaluationSummaryRows() []EvaluationResultsSummaryRow {
	return []EvaluationResultsSummaryRow{
		{
			Instance: "a",
			Metrics: map[string]EvaluationResultSummaryMetric{
				testHeuristicBaseline: {
					AverageMinDeviation:  4.0,
					SuccessRate:          10.0,
					AverageBestIteration: 80.0,
					Iterations:           100,
				},
				testHeuristicStrictMsa: {
					AverageMinDeviation:  2.5,
					SuccessRate:          20.0,
					AverageBestIteration: 60.0,
					Iterations:           100,
				},
				testHeuristicCycleCover: {
					AverageMinDeviation:  2.0,
					SuccessRate:          30.0,
					AverageBestIteration: 40.0,
					Iterations:           100,
				},
				testHeuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  1.8,
					SuccessRate:          25.0,
					AverageBestIteration: 35.0,
					Iterations:           100,
				},
			},
		},
		{
			Instance: "b",
			Metrics: map[string]EvaluationResultSummaryMetric{
				testHeuristicBaseline: {
					AverageMinDeviation:  5.0,
					SuccessRate:          0.0,
					AverageBestIteration: 70.0,
					Iterations:           100,
				},
				testHeuristicStrictMsa: {
					AverageMinDeviation:  3.5,
					SuccessRate:          20.0,
					AverageBestIteration: 40.0,
					Iterations:           100,
				},
				testHeuristicCycleCover: {
					AverageMinDeviation:  3.5,
					SuccessRate:          20.0,
					AverageBestIteration: 30.0,
					Iterations:           100,
				},
				testHeuristicCycleCoverMsaPatching: {
					AverageMinDeviation:  3.0,
					SuccessRate:          25.0,
					AverageBestIteration: 20.0,
					Iterations:           100,
				},
			},
		},
	}
}
