package reports

import (
	"atsp_aco_msa/modules/analysis/structure"
	"fmt"
	"html"
	"os"
	"path/filepath"
	"strings"
)

func SaveMsaHeuristicCycleCoverOverlap(path string, analyses []structure.InstanceAnalysis) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows := sortedStructuralAnalyses(analyses)
	totals := msaHeuristicCycleCoverOverlapTotals(rows)

	var builder strings.Builder
	builder.WriteString("# MSA Heuristic And Cycle-Cover Overlap\n\n")
	builder.WriteString("This table compares the current MSA heuristic edge set directly with the minimum cycle-cover edge set.\n\n")
	writeMsaHeuristicCycleCoverOverlapFindings(&builder, totals)
	builder.WriteString("\n")
	writeMsaHeuristicCycleCoverOverlapTable(&builder, rows, totals)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func writeMsaHeuristicCycleCoverOverlapFindings(builder *strings.Builder, totals msaHeuristicCycleCoverOverlapSummary) {
	builder.WriteString("## Findings\n\n")
	fmt.Fprintf(builder, "- **%.2f%% of MSA heuristic edges are also cycle-cover edges.**\n", 100*ratio(totals.sharedEdges, totals.msaEdges))
	fmt.Fprintf(builder, "- **%.2f%% of cycle-cover edges are also MSA heuristic edges.**\n", 100*ratio(totals.sharedEdges, totals.cycleCoverEdges))
	fmt.Fprintf(builder,
		"- **Found-optimal edge partition: both %d, only MSA heuristic %d, only cycle cover %d, neither %d.**\n",
		totals.optimalBoth,
		totals.optimalOnlyMsa,
		totals.optimalOnlyCycleCover,
		totals.optimalNeither)
}

func writeMsaHeuristicCycleCoverOverlapTable(builder *strings.Builder, rows []structure.InstanceAnalysis, totals msaHeuristicCycleCoverOverlapSummary) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>MSA in CC [%]</th><th>CC in MSA [%]</th><th>Optimal both</th><th>Optimal only MSA</th><th>Optimal only CC</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, analysis := range rows {
		writeMsaHeuristicCycleCoverOverlapRow(builder, analysis)
	}
	writeMsaHeuristicCycleCoverOverlapTotalRow(builder, totals)

	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeMsaHeuristicCycleCoverOverlapRow(builder *strings.Builder, analysis structure.InstanceAnalysis) {
	metrics := analysis.Metrics
	msaEdges := metrics.HighMsaHeuristicMetrics.EdgeCount
	cycleCoverEdges := metrics.CycleCoverMetrics.EdgeCount
	sharedEdges := metrics.CycleCoverHighMsaHeuristicEdges

	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td></tr>\n",
		html.EscapeString(analysis.Instance),
		100*ratio(sharedEdges, msaEdges),
		100*ratio(sharedEdges, cycleCoverEdges),
		metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic,
		metrics.OptimalEdgesInHighMsaHeuristicNotCycleCover,
		metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic)
}

func writeMsaHeuristicCycleCoverOverlapTotalRow(builder *strings.Builder, totals msaHeuristicCycleCoverOverlapSummary) {
	fmt.Fprintf(builder,
		"<tr><td><strong>Total</strong></td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%d</td><td align=\"right\">%d</td><td align=\"right\">%d</td></tr>\n",
		100*ratio(totals.sharedEdges, totals.msaEdges),
		100*ratio(totals.sharedEdges, totals.cycleCoverEdges),
		totals.optimalBoth,
		totals.optimalOnlyMsa,
		totals.optimalOnlyCycleCover)
}

type msaHeuristicCycleCoverOverlapSummary struct {
	msaEdges              int
	cycleCoverEdges       int
	sharedEdges           int
	optimalBoth           int
	optimalOnlyMsa        int
	optimalOnlyCycleCover int
	optimalNeither        int
}

func msaHeuristicCycleCoverOverlapTotals(rows []structure.InstanceAnalysis) msaHeuristicCycleCoverOverlapSummary {
	var totals msaHeuristicCycleCoverOverlapSummary
	for _, analysis := range rows {
		metrics := analysis.Metrics
		totals.msaEdges += metrics.HighMsaHeuristicMetrics.EdgeCount
		totals.cycleCoverEdges += metrics.CycleCoverMetrics.EdgeCount
		totals.sharedEdges += metrics.CycleCoverHighMsaHeuristicEdges
		totals.optimalBoth += metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic
		totals.optimalOnlyMsa += metrics.OptimalEdgesInHighMsaHeuristicNotCycleCover
		totals.optimalOnlyCycleCover += metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic
		totals.optimalNeither += metrics.OptimalEdgesInNeitherCycleCoverNorHigh
	}

	return totals
}
