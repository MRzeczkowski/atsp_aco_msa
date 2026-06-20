package app

import (
	"atsp_aco_msa/modules/analysis/structure"
	"sort"
)

func maxIntValue(left, right int) int {
	if left > right {
		return left
	}

	return right
}

func ratio(numerator, denominator int) float64 {
	if denominator == 0 {
		return 0
	}

	return float64(numerator) / float64(denominator)
}

type structuralSimilaritySummary struct {
	instanceCount           int
	foundOptimalTours       int
	foundOptimalEdges       int
	tourEdges               int
	msaEdges                int
	msaOptimalEdges         int
	cycleCoverEdges         int
	cycleCoverOptimalEdges  int
	patchingEdges           int
	patchingOptimalEdges    int
	msaPrecisionWins        int
	cycleCoverPrecisionWins int
	patchingPrecisionWins   int
	msaRecallWins           int
	cycleCoverRecallWins    int
	patchingRecallWins      int
}

func structuralSimilarityTotals(rows []structure.InstanceAnalysis) structuralSimilaritySummary {
	var totals structuralSimilaritySummary
	totals.instanceCount = len(rows)

	for _, analysis := range rows {
		metrics := analysis.Metrics
		msaMetrics := metrics.HighMsaHeuristicMetrics
		cycleCoverMetrics := metrics.CycleCoverMetrics
		patchingMetrics := metrics.CycleCoverMsaPatchingMetrics

		totals.foundOptimalTours += metrics.FoundOptimalTourCount
		totals.foundOptimalEdges += metrics.UniqueFoundOptimalEdgeCount
		totals.tourEdges += analysis.Dimension
		totals.msaEdges += msaMetrics.EdgeCount
		totals.msaOptimalEdges += msaMetrics.OptimalEdgeCount
		totals.cycleCoverEdges += cycleCoverMetrics.EdgeCount
		totals.cycleCoverOptimalEdges += cycleCoverMetrics.OptimalEdgeCount
		totals.patchingEdges += patchingMetrics.EdgeCount
		totals.patchingOptimalEdges += patchingMetrics.OptimalEdgeCount
	}

	return totals
}

func sortedStructuralAnalyses(analyses []structure.InstanceAnalysis) []structure.InstanceAnalysis {
	rows := append([]structure.InstanceAnalysis(nil), analyses...)
	sort.SliceStable(rows, func(i, j int) bool {
		return rows[i].Instance < rows[j].Instance
	})

	return rows
}

func filterAnalysesWithFoundOptimalEdges(analyses []structure.InstanceAnalysis) []structure.InstanceAnalysis {
	rows := make([]structure.InstanceAnalysis, 0, len(analyses))
	for _, analysis := range analyses {
		if analysis.Metrics.UniqueFoundOptimalEdgeCount == 0 {
			continue
		}

		rows = append(rows, analysis)
	}

	return rows
}
