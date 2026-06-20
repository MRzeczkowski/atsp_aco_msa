package reports

import (
	"atsp_aco_msa/modules/analysis/structure"
	"fmt"
	"math"
	"sort"
)

func maxIntValue(left, right int) int {
	if left > right {
		return left
	}

	return right
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

func bestStructuralMetricHighlights(values ...float64) []bool {
	const epsilon = 1e-9
	highlights := make([]bool, len(values))
	best := math.Inf(-1)
	for _, value := range values {
		if value > best {
			best = value
		}
	}
	if best <= 0 {
		return highlights
	}

	for index, value := range values {
		if math.Abs(value-best) < epsilon {
			highlights[index] = true
		}
	}

	return highlights
}

func ratio(numerator, denominator int) float64 {
	if denominator == 0 {
		return 0
	}

	return float64(numerator) / float64(denominator)
}

func metricCell(value float64, bold bool) string {
	valueText := fmt.Sprintf("%.2f", value)
	if bold {
		return "<strong>" + valueText + "</strong>"
	}

	return valueText
}
