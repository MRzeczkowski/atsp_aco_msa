package reports

import (
	"atsp_aco_msa/modules/experiments"
	"atsp_aco_msa/modules/project"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const (
	testHeuristicRandomSparse               = "random-sparse"
	testHeuristicDistanceRankedSparse       = "distance-ranked-sparse"
	testHeuristicShuffledMsa                = "shuffled-msa"
	testEvaluationStrictMsaHeuristicWeight  = 0.4
	testEvaluationControlsResultsFolderName = "controls"
)

func TestSaveRandomSparseControlReportComparesMsaAgainstRandomSparse(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	evaluationRoot := filepath.Join(root, "evaluation")
	controlsRoot := filepath.Join(root, testEvaluationControlsResultsFolderName)
	first := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	evaluationFirst := project.WithExperimentOutputRoot(first, evaluationRoot)
	evaluationSecond := project.WithExperimentOutputRoot(second, evaluationRoot)
	controlFirst := project.WithExperimentOutputRoot(first, controlsRoot)
	controlSecond := project.WithExperimentOutputRoot(second, controlsRoot)

	if err := experiments.SaveHeuristicStatistics(evaluationFirst.ResultFilePath, []experiments.HeuristicExperimentStatistics{
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 2.0, 10.0)},
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(0.8, 0.5, 100.0)},
	}); err != nil {
		t.Fatalf("failed to write first evaluation MSA result: %v", err)
	}
	experiments.SaveStatistics(controlReportTestResultFilePathForHeuristic(controlFirst, testHeuristicRandomSparse), testHeuristicRandomSparse, []experiments.ExperimentsDataStatistics{
		makeControlReportTestRandomSparseStatistics(1, 4.0, 0.0),
		makeControlReportTestRandomSparseStatistics(2, 5.0, 0.0),
		makeControlReportTestRandomSparseStatistics(3, 6.0, 0.0),
	})

	if err := experiments.SaveHeuristicStatistics(evaluationSecond.ResultFilePath, []experiments.HeuristicExperimentStatistics{
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 4.0, 0.0)},
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(0.8, 0.5, 100.0)},
	}); err != nil {
		t.Fatalf("failed to write second evaluation MSA result: %v", err)
	}
	experiments.SaveStatistics(controlReportTestResultFilePathForHeuristic(controlSecond, testHeuristicRandomSparse), testHeuristicRandomSparse, []experiments.ExperimentsDataStatistics{
		makeControlReportTestRandomSparseStatistics(1, 2.0, 20.0),
		makeControlReportTestRandomSparseStatistics(2, 3.0, 20.0),
		makeControlReportTestRandomSparseStatistics(3, 4.0, 20.0),
	})

	reportPath := filepath.Join(controlsRoot, "random_sparse_control.md")
	saved, err := SaveRandomSparseControlReport(reportPath, []project.AtspData{second, first}, evaluationRoot, controlsRoot, controlReportsTestConfig())
	if err != nil {
		t.Fatalf("SaveRandomSparseControlReport returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected random sparse control report to be saved")
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read random sparse control report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "Strict MSA had lower average best deviation than the random-sparse mean in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: Strict MSA 3.00%, random sparse 4.00%, delta -1.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func TestSaveDistanceRankedSparseControlReportComparesMsaAgainstDistanceRankedSparse(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	evaluationRoot := filepath.Join(root, "evaluation")
	controlsRoot := filepath.Join(root, testEvaluationControlsResultsFolderName)
	first := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	evaluationFirst := project.WithExperimentOutputRoot(first, evaluationRoot)
	evaluationSecond := project.WithExperimentOutputRoot(second, evaluationRoot)
	controlFirst := project.WithExperimentOutputRoot(first, controlsRoot)
	controlSecond := project.WithExperimentOutputRoot(second, controlsRoot)

	if err := experiments.SaveHeuristicStatistics(evaluationFirst.ResultFilePath, []experiments.HeuristicExperimentStatistics{
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 2.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write first evaluation MSA result: %v", err)
	}
	experiments.SaveStatistics(controlReportTestResultFilePathForHeuristic(controlFirst, testHeuristicDistanceRankedSparse), testHeuristicDistanceRankedSparse, []experiments.ExperimentsDataStatistics{
		makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 4.0, 0.0),
	})

	if err := experiments.SaveHeuristicStatistics(evaluationSecond.ResultFilePath, []experiments.HeuristicExperimentStatistics{
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 4.0, 0.0)},
	}); err != nil {
		t.Fatalf("failed to write second evaluation MSA result: %v", err)
	}
	experiments.SaveStatistics(controlReportTestResultFilePathForHeuristic(controlSecond, testHeuristicDistanceRankedSparse), testHeuristicDistanceRankedSparse, []experiments.ExperimentsDataStatistics{
		makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 2.0, 20.0),
	})

	reportPath := filepath.Join(controlsRoot, "distance_ranked_sparse_control.md")
	saved, err := SaveDistanceRankedSparseControlReport(reportPath, []project.AtspData{second, first}, evaluationRoot, controlsRoot, controlReportsTestConfig())
	if err != nil {
		t.Fatalf("SaveDistanceRankedSparseControlReport returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected distance-ranked sparse control report to be saved")
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read distance-ranked sparse control report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "Strict MSA had lower average best deviation than the distance-ranked sparse control in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: Strict MSA 3.00%, distance-ranked sparse 3.00%, delta +0.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func TestSaveShuffledMsaControlReportComparesMsaAgainstShuffledMsa(t *testing.T) {
	root := t.TempDir()
	sourceRoot := filepath.Join(root, "results")
	evaluationRoot := filepath.Join(root, "evaluation")
	controlsRoot := filepath.Join(root, testEvaluationControlsResultsFolderName)
	first := project.MakeAtspDataInResultsDirectory("sample-a.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	second := project.MakeAtspDataInResultsDirectory("sample-b.atsp", [][]float64{{0, 1}, {1, 0}}, 2, sourceRoot)
	evaluationFirst := project.WithExperimentOutputRoot(first, evaluationRoot)
	evaluationSecond := project.WithExperimentOutputRoot(second, evaluationRoot)
	controlFirst := project.WithExperimentOutputRoot(first, controlsRoot)
	controlSecond := project.WithExperimentOutputRoot(second, controlsRoot)

	if err := experiments.SaveHeuristicStatistics(evaluationFirst.ResultFilePath, []experiments.HeuristicExperimentStatistics{
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 2.0, 10.0)},
	}); err != nil {
		t.Fatalf("failed to write first evaluation MSA result: %v", err)
	}
	experiments.SaveStatistics(controlReportTestResultFilePathForHeuristic(controlFirst, testHeuristicShuffledMsa), testHeuristicShuffledMsa, []experiments.ExperimentsDataStatistics{
		makeControlReportTestRandomSparseStatistics(1, 4.0, 0.0),
		makeControlReportTestRandomSparseStatistics(2, 5.0, 0.0),
		makeControlReportTestRandomSparseStatistics(3, 6.0, 0.0),
	})

	if err := experiments.SaveHeuristicStatistics(evaluationSecond.ResultFilePath, []experiments.HeuristicExperimentStatistics{
		{Heuristic: testHeuristicStrictMsa, Statistics: makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, 4.0, 0.0)},
	}); err != nil {
		t.Fatalf("failed to write second evaluation MSA result: %v", err)
	}
	experiments.SaveStatistics(controlReportTestResultFilePathForHeuristic(controlSecond, testHeuristicShuffledMsa), testHeuristicShuffledMsa, []experiments.ExperimentsDataStatistics{
		makeControlReportTestRandomSparseStatistics(1, 2.0, 20.0),
		makeControlReportTestRandomSparseStatistics(2, 3.0, 20.0),
		makeControlReportTestRandomSparseStatistics(3, 4.0, 20.0),
	})

	reportPath := filepath.Join(controlsRoot, "shuffled_msa_control.md")
	saved, err := SaveShuffledMsaControlReport(reportPath, []project.AtspData{second, first}, evaluationRoot, controlsRoot, controlReportsTestConfig())
	if err != nil {
		t.Fatalf("SaveShuffledMsaControlReport returned unexpected error: %v", err)
	}
	if !saved {
		t.Fatal("expected shuffled MSA control report to be saved")
	}

	contentBytes, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("failed to read shuffled MSA control report: %v", err)
	}
	content := string(contentBytes)
	assertContains(t, content, "Strict MSA had lower average best deviation than the shuffled MSA mean in 1/2 instances.")
	assertContains(t, content, "Mean average best deviation: Strict MSA 3.00%, shuffled MSA 4.00%, delta -1.00 pp.")
	assertContains(t, content, "<td>sample-a</td>")
	assertContains(t, content, "<td>sample-b</td>")
}

func controlReportsTestConfig() ControlReportsConfig {
	return ControlReportsConfig{
		StrictMsaHeuristic:                 testHeuristicStrictMsa,
		RandomSparseHeuristic:              testHeuristicRandomSparse,
		DistanceRankedSparseHeuristic:      testHeuristicDistanceRankedSparse,
		ShuffledMsaHeuristic:               testHeuristicShuffledMsa,
		EvaluationStrictMsaHeuristicWeight: testEvaluationStrictMsaHeuristicWeight,
		ReadStatistics:                     experiments.ReadStatistics,
		ReadHeuristicStatistics:            experiments.ReadHeuristicStatistics,
		ResultFilePathForHeuristic:         controlReportTestResultFilePathForHeuristic,
	}
}

func controlReportTestResultFilePathForHeuristic(atspData project.AtspData, heuristic string) string {
	if heuristic == testHeuristicStrictMsa {
		return atspData.ResultFilePath
	}

	return strings.TrimSuffix(atspData.ResultFilePath, ".csv") + "_" + strings.ReplaceAll(heuristic, "-", "_") + ".csv"
}

func makeControlReportTestExperimentStatistics(heuristicWeight, averageBestDeviation, successRate float64) experiments.ExperimentsDataStatistics {
	return experiments.ExperimentsDataStatistics{
		ExperimentParameters: experiments.ExperimentParameters{
			Alpha:           1.0,
			Beta:            2.0,
			Rho:             0.8,
			HeuristicWeight: heuristicWeight,
			Iterations:      500,
		},
		MinBestAtIteration:               1,
		AverageBestAtIteration:           2.0,
		MaxBestAtIteration:               3,
		MinThreeOptImprovementsCount:     0,
		AverageThreeOptImprovementsCount: 0.0,
		MaxThreeOptImprovementsCount:     0,
		MinBestDeviation:                 averageBestDeviation,
		AverageBestDeviation:             averageBestDeviation,
		MaxBestDeviation:                 averageBestDeviation,
		SuccessRate:                      successRate,
	}
}

func makeControlReportTestRandomSparseStatistics(randomSeed int64, averageBestDeviation, successRate float64) experiments.ExperimentsDataStatistics {
	statistics := makeControlReportTestExperimentStatistics(testEvaluationStrictMsaHeuristicWeight, averageBestDeviation, successRate)
	statistics.RandomSeed = randomSeed
	return statistics
}
