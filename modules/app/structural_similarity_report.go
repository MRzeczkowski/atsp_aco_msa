package app

import (
	"atsp_aco_msa/modules/analysis/structure"
	"fmt"
	"html"
	"os"
	"path/filepath"
	"strings"
)

func saveStructuralSimilarityReport(path string, analyses []structure.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows := filterAnalysesWithFoundOptimalEdges(sortedStructuralAnalyses(analyses))
	totals := structuralSimilarityTotals(rows)

	var builder strings.Builder
	builder.WriteString("# Structural Similarity To Found Optimal Tours\n\n")
	builder.WriteString("This table compares the current MSA heuristic, minimum cycle-cover, and cycle-cover MSA-patching edge sets against the found optimal tours saved in `solutions.csv`.\n\n")
	builder.WriteString("Instances without found optimal tours are omitted because precision and recall cannot be interpreted without a reference edge set.\n\n")
	writeStructuralSimilarityFindings(&builder, totals)
	builder.WriteString("\n")
	writeStructuralSimilarityTable(&builder, rows, totals)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeStructuralSimilarityFindings(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	patchingPrecision := ratio(totals.patchingOptimalEdges, totals.patchingEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)
	patchingRecall := ratio(totals.patchingOptimalEdges, totals.foundOptimalEdges)

	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **Precision vs found-optimal tours: MSA heuristic %.2f%%, cycle cover %.2f%%, cycle-cover MSA patching %.2f%%.**\n", 100*msaPrecision, 100*cycleCoverPrecision, 100*patchingPrecision)
	fmt.Fprintf(builder, "- **Recall vs found-optimal tours: MSA heuristic %.2f%%, cycle cover %.2f%%, cycle-cover MSA patching %.2f%%.**\n", 100*msaRecall, 100*cycleCoverRecall, 100*patchingRecall)
	fmt.Fprintf(builder, "- **Best-or-tied precision counts: MSA heuristic %d/%d, cycle cover %d/%d, cycle-cover MSA patching %d/%d.**\n",
		totals.msaPrecisionWins,
		totals.instanceCount,
		totals.cycleCoverPrecisionWins,
		totals.instanceCount,
		totals.patchingPrecisionWins,
		totals.instanceCount)
	fmt.Fprintf(builder, "- **Best-or-tied recall counts: MSA heuristic %d/%d, cycle cover %d/%d, cycle-cover MSA patching %d/%d.**\n",
		totals.msaRecallWins,
		totals.instanceCount,
		totals.cycleCoverRecallWins,
		totals.instanceCount,
		totals.patchingRecallWins,
		totals.instanceCount)
}

func writeStructuralSimilarityTable(builder *strings.Builder, rows []structure.InstanceAnalysis, totals structuralSimilaritySummary) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th rowspan=\"2\">Instance</th><th colspan=\"2\">MSA heuristic</th><th colspan=\"2\">Cycle cover</th><th colspan=\"2\">Cycle-cover MSA patching</th></tr>\n")
	builder.WriteString("<tr><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>\n")
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
	patchingMetrics := metrics.CycleCoverMsaPatchingMetrics
	precisionHighlights := bestStructuralMetricHighlights(msaMetrics.Precision, cycleCoverMetrics.Precision, patchingMetrics.Precision)
	recallHighlights := bestStructuralMetricHighlights(msaMetrics.Recall, cycleCoverMetrics.Recall, patchingMetrics.Recall)

	writeStructuralSimilarityTableRow(
		builder,
		html.EscapeString(analysis.Instance),
		msaMetrics.Precision,
		msaMetrics.Recall,
		cycleCoverMetrics.Precision,
		cycleCoverMetrics.Recall,
		patchingMetrics.Precision,
		patchingMetrics.Recall,
		precisionHighlights,
		recallHighlights)
}

func writeStructuralSimilarityTotalRow(builder *strings.Builder, totals structuralSimilaritySummary) {
	msaPrecision := ratio(totals.msaOptimalEdges, totals.msaEdges)
	cycleCoverPrecision := ratio(totals.cycleCoverOptimalEdges, totals.cycleCoverEdges)
	patchingPrecision := ratio(totals.patchingOptimalEdges, totals.patchingEdges)
	msaRecall := ratio(totals.msaOptimalEdges, totals.foundOptimalEdges)
	cycleCoverRecall := ratio(totals.cycleCoverOptimalEdges, totals.foundOptimalEdges)
	patchingRecall := ratio(totals.patchingOptimalEdges, totals.foundOptimalEdges)
	precisionHighlights := bestStructuralMetricHighlights(msaPrecision, cycleCoverPrecision, patchingPrecision)
	recallHighlights := bestStructuralMetricHighlights(msaRecall, cycleCoverRecall, patchingRecall)

	writeStructuralSimilarityTableRow(
		builder,
		"<strong>Total</strong>",
		msaPrecision,
		msaRecall,
		cycleCoverPrecision,
		cycleCoverRecall,
		patchingPrecision,
		patchingRecall,
		precisionHighlights,
		recallHighlights)
}

func writeStructuralSimilarityTableRow(builder *strings.Builder, instanceCell string, msaPrecision, msaRecall, cycleCoverPrecision, cycleCoverRecall, patchingPrecision, patchingRecall float64, precisionHighlights, recallHighlights []bool) {
	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		instanceCell,
		finalResultsSummaryMetricCell(100*msaPrecision, precisionHighlights[0]),
		finalResultsSummaryMetricCell(100*msaRecall, recallHighlights[0]),
		finalResultsSummaryMetricCell(100*cycleCoverPrecision, precisionHighlights[1]),
		finalResultsSummaryMetricCell(100*cycleCoverRecall, recallHighlights[1]),
		finalResultsSummaryMetricCell(100*patchingPrecision, precisionHighlights[2]),
		finalResultsSummaryMetricCell(100*patchingRecall, recallHighlights[2]))
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
		precisionWins := bestStructuralMetricHighlights(msaMetrics.Precision, cycleCoverMetrics.Precision, patchingMetrics.Precision)
		recallWins := bestStructuralMetricHighlights(msaMetrics.Recall, cycleCoverMetrics.Recall, patchingMetrics.Recall)

		totals.foundOptimalTours += metrics.FoundOptimalTourCount
		totals.foundOptimalEdges += metrics.UniqueFoundOptimalEdgeCount
		totals.tourEdges += analysis.Dimension
		totals.msaEdges += msaMetrics.EdgeCount
		totals.msaOptimalEdges += msaMetrics.OptimalEdgeCount
		totals.cycleCoverEdges += cycleCoverMetrics.EdgeCount
		totals.cycleCoverOptimalEdges += cycleCoverMetrics.OptimalEdgeCount
		totals.patchingEdges += patchingMetrics.EdgeCount
		totals.patchingOptimalEdges += patchingMetrics.OptimalEdgeCount
		if precisionWins[0] {
			totals.msaPrecisionWins++
		}
		if precisionWins[1] {
			totals.cycleCoverPrecisionWins++
		}
		if precisionWins[2] {
			totals.patchingPrecisionWins++
		}
		if recallWins[0] {
			totals.msaRecallWins++
		}
		if recallWins[1] {
			totals.cycleCoverRecallWins++
		}
		if recallWins[2] {
			totals.patchingRecallWins++
		}
	}

	return totals
}
