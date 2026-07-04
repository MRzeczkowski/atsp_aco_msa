package reports

import (
	"atsp_aco_msa/modules/analysis/structure"
	"fmt"
	"html"
	"os"
	"path/filepath"
	"strings"
)

func SaveStructuralSimilarity(path string, analyses []structure.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows := filterAnalysesWithFoundOptimalEdges(sortedStructuralAnalyses(analyses))
	totals := structuralSimilarityTotals(rows)

	var builder strings.Builder
	builder.WriteString("# Structural Similarity To Found Optimal Tours\n\n")
	builder.WriteString("This table compares the current MSA heuristic, minimum cycle-cover, GKS patching, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.\n\n")
	builder.WriteString("Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.\n\n")
	writeStructuralSimilarityFindings(&builder, totals)
	builder.WriteString("\n")
	writeStructuralSimilarityTable(&builder, rows, totals)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeStructuralSimilarityFindings(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	cycleCoverPatchingPrecision := ratio(totals.cycleCoverPatchingOptimalEdges, totals.cycleCoverPatchingEdges)
	cycleCoverMsaPatchingPrecision := ratio(totals.cycleCoverMsaPatchingOptimalEdges, totals.cycleCoverMsaPatchingEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)
	cycleCoverPatchingRecall := ratio(totals.cycleCoverPatchingOptimalEdges, totals.foundOptimalEdges)
	cycleCoverMsaPatchingRecall := ratio(totals.cycleCoverMsaPatchingOptimalEdges, totals.foundOptimalEdges)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **Precision vs found-optimal tours: MSA heuristic %.2f%%, cycle cover %.2f%%, GKS patching %.2f%%, cycle-cover MSA patching %.2f%%.**\n", 100*msaPrecision, 100*cycleCoverPrecision, 100*cycleCoverPatchingPrecision, 100*cycleCoverMsaPatchingPrecision)
	fmt.Fprintf(builder, "- **Recall vs found-optimal tours: MSA heuristic %.2f%%, cycle cover %.2f%%, GKS patching %.2f%%, cycle-cover MSA patching %.2f%%.**\n", 100*msaRecall, 100*cycleCoverRecall, 100*cycleCoverPatchingRecall, 100*cycleCoverMsaPatchingRecall)
	fmt.Fprintf(builder, "- **Best-or-tied precision counts: MSA heuristic %d/%d, cycle cover %d/%d, GKS patching %d/%d, cycle-cover MSA patching %d/%d.**\n",
		totals.msaPrecisionWins,
		totals.instanceCount,
		totals.cycleCoverPrecisionWins,
		totals.instanceCount,
		totals.cycleCoverPatchingPrecisionWins,
		totals.instanceCount,
		totals.cycleCoverMsaPatchingPrecisionWins,
		totals.instanceCount)
	fmt.Fprintf(builder, "- **Best-or-tied recall counts: MSA heuristic %d/%d, cycle cover %d/%d, GKS patching %d/%d, cycle-cover MSA patching %d/%d.**\n",
		totals.msaRecallWins,
		totals.instanceCount,
		totals.cycleCoverRecallWins,
		totals.instanceCount,
		totals.cycleCoverPatchingRecallWins,
		totals.instanceCount,
		totals.cycleCoverMsaPatchingRecallWins,
		totals.instanceCount)
}

func writeStructuralSimilarityTable(builder *strings.Builder, rows []structure.InstanceAnalysis, totals structuralSimilaritySummary) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">MSA heuristic</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">GKS patching</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>\n")
	builder.WriteString("<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, analysis := range rows {
		writeStructuralSimilarityRow(builder, analysis)
	}
	writeStructuralSimilarityTotalRow(builder, totals)

	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeStructuralSimilarityRow(builder *strings.Builder, analysis structure.InstanceAnalysis) {
	metrics := analysis.Metrics
	msaMetrics := metrics.HighMsaHeuristicMetrics
	cycleCoverMetrics := metrics.CycleCoverMetrics
	cycleCoverPatchingMetrics := metrics.CycleCoverPatchingMetrics
	cycleCoverMsaPatchingMetrics := metrics.CycleCoverMsaPatchingMetrics
	precisionHighlights := bestStructuralMetricHighlights(msaMetrics.Precision, cycleCoverMetrics.Precision, cycleCoverPatchingMetrics.Precision, cycleCoverMsaPatchingMetrics.Precision)
	recallHighlights := bestStructuralMetricHighlights(msaMetrics.Recall, cycleCoverMetrics.Recall, cycleCoverPatchingMetrics.Recall, cycleCoverMsaPatchingMetrics.Recall)

	writeStructuralSimilarityTableRow(
		builder,
		html.EscapeString(analysis.Instance),
		msaMetrics.Precision,
		msaMetrics.Recall,
		cycleCoverMetrics.Precision,
		cycleCoverMetrics.Recall,
		cycleCoverPatchingMetrics.Precision,
		cycleCoverPatchingMetrics.Recall,
		cycleCoverMsaPatchingMetrics.Precision,
		cycleCoverMsaPatchingMetrics.Recall,
		precisionHighlights,
		recallHighlights)
}

func writeStructuralSimilarityTotalRow(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	cycleCoverPatchingPrecision := ratio(totals.cycleCoverPatchingOptimalEdges, totals.cycleCoverPatchingEdges)
	cycleCoverMsaPatchingPrecision := ratio(totals.cycleCoverMsaPatchingOptimalEdges, totals.cycleCoverMsaPatchingEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)
	cycleCoverPatchingRecall := ratio(totals.cycleCoverPatchingOptimalEdges, totals.foundOptimalEdges)
	cycleCoverMsaPatchingRecall := ratio(totals.cycleCoverMsaPatchingOptimalEdges, totals.foundOptimalEdges)
	precisionHighlights := bestStructuralMetricHighlights(msaPrecision, cycleCoverPrecision, cycleCoverPatchingPrecision, cycleCoverMsaPatchingPrecision)
	recallHighlights := bestStructuralMetricHighlights(msaRecall, cycleCoverRecall, cycleCoverPatchingRecall, cycleCoverMsaPatchingRecall)

	writeStructuralSimilarityTableRow(
		builder,
		"<strong>Total</strong>",
		msaPrecision,
		msaRecall,
		cycleCoverPrecision,
		cycleCoverRecall,
		cycleCoverPatchingPrecision,
		cycleCoverPatchingRecall,
		cycleCoverMsaPatchingPrecision,
		cycleCoverMsaPatchingRecall,
		precisionHighlights,
		recallHighlights)
}

func writeStructuralSimilarityTableRow(builder *strings.Builder, instanceCell string, msaPrecision, msaRecall, cycleCoverPrecision, cycleCoverRecall, cycleCoverPatchingPrecision, cycleCoverPatchingRecall, cycleCoverMsaPatchingPrecision, cycleCoverMsaPatchingRecall float64, precisionHighlights, recallHighlights []bool) {
	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		instanceCell,
		metricCell(100*msaPrecision, precisionHighlights[0]),
		metricCell(100*msaRecall, recallHighlights[0]),
		metricCell(100*cycleCoverPrecision, precisionHighlights[1]),
		metricCell(100*cycleCoverRecall, recallHighlights[1]),
		metricCell(100*cycleCoverPatchingPrecision, precisionHighlights[2]),
		metricCell(100*cycleCoverPatchingRecall, recallHighlights[2]),
		metricCell(100*cycleCoverMsaPatchingPrecision, precisionHighlights[3]),
		metricCell(100*cycleCoverMsaPatchingRecall, recallHighlights[3]))
}

type structuralSimilaritySummary struct {
	instanceCount                      int
	foundOptimalTours                  int
	foundOptimalEdges                  int
	tourEdges                          int
	msaEdges                           int
	msaOptimalEdges                    int
	cycleCoverEdges                    int
	cycleCoverOptimalEdges             int
	cycleCoverPatchingEdges            int
	cycleCoverPatchingOptimalEdges     int
	cycleCoverMsaPatchingEdges         int
	cycleCoverMsaPatchingOptimalEdges  int
	msaPrecisionWins                   int
	cycleCoverPrecisionWins            int
	cycleCoverPatchingPrecisionWins    int
	cycleCoverMsaPatchingPrecisionWins int
	msaRecallWins                      int
	cycleCoverRecallWins               int
	cycleCoverPatchingRecallWins       int
	cycleCoverMsaPatchingRecallWins    int
}

func structuralSimilarityTotals(rows []structure.InstanceAnalysis) structuralSimilaritySummary {
	var totals structuralSimilaritySummary
	totals.instanceCount = len(rows)

	for _, analysis := range rows {
		metrics := analysis.Metrics
		msaMetrics := metrics.HighMsaHeuristicMetrics
		cycleCoverMetrics := metrics.CycleCoverMetrics
		cycleCoverPatchingMetrics := metrics.CycleCoverPatchingMetrics
		cycleCoverMsaPatchingMetrics := metrics.CycleCoverMsaPatchingMetrics
		precisionWins := bestStructuralMetricHighlights(msaMetrics.Precision, cycleCoverMetrics.Precision, cycleCoverPatchingMetrics.Precision, cycleCoverMsaPatchingMetrics.Precision)
		recallWins := bestStructuralMetricHighlights(msaMetrics.Recall, cycleCoverMetrics.Recall, cycleCoverPatchingMetrics.Recall, cycleCoverMsaPatchingMetrics.Recall)

		totals.foundOptimalTours += metrics.FoundOptimalTourCount
		totals.foundOptimalEdges += metrics.UniqueFoundOptimalEdgeCount
		totals.tourEdges += analysis.Dimension
		totals.msaEdges += msaMetrics.EdgeCount
		totals.msaOptimalEdges += msaMetrics.OptimalEdgeCount
		totals.cycleCoverEdges += cycleCoverMetrics.EdgeCount
		totals.cycleCoverOptimalEdges += cycleCoverMetrics.OptimalEdgeCount
		totals.cycleCoverPatchingEdges += cycleCoverPatchingMetrics.EdgeCount
		totals.cycleCoverPatchingOptimalEdges += cycleCoverPatchingMetrics.OptimalEdgeCount
		totals.cycleCoverMsaPatchingEdges += cycleCoverMsaPatchingMetrics.EdgeCount
		totals.cycleCoverMsaPatchingOptimalEdges += cycleCoverMsaPatchingMetrics.OptimalEdgeCount
		if precisionWins[0] {
			totals.msaPrecisionWins++
		}
		if precisionWins[1] {
			totals.cycleCoverPrecisionWins++
		}
		if precisionWins[2] {
			totals.cycleCoverPatchingPrecisionWins++
		}
		if precisionWins[3] {
			totals.cycleCoverMsaPatchingPrecisionWins++
		}
		if recallWins[0] {
			totals.msaRecallWins++
		}
		if recallWins[1] {
			totals.cycleCoverRecallWins++
		}
		if recallWins[2] {
			totals.cycleCoverPatchingRecallWins++
		}
		if recallWins[3] {
			totals.cycleCoverMsaPatchingRecallWins++
		}
	}

	return totals
}
