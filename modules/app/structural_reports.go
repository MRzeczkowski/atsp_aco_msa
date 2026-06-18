package app

import (
	"atsp_aco_msa/modules/algorithms/heuristics"
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	"atsp_aco_msa/modules/analysis/msaHeuristicTours"
	"atsp_aco_msa/modules/analysis/structuralComparison"
	"atsp_aco_msa/modules/models"
	"fmt"
	"html"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

func saveStructuralSimilarityReport(path string, analyses []structuralComparison.InstanceAnalysis) error {
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

func writeStructuralSimilarityTable(builder *strings.Builder, rows []structuralComparison.InstanceAnalysis, totals structuralSimilaritySummary) {
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

func writeStructuralSimilarityRow(builder *strings.Builder, analysis structuralComparison.InstanceAnalysis) {
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

func structuralSimilarityTotals(rows []structuralComparison.InstanceAnalysis) structuralSimilaritySummary {
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

func saveMsaHeuristicCycleCoverOverlapReport(path string, analyses []structuralComparison.InstanceAnalysis) error {
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

type gksDeviationRow struct {
	instance     string
	msaPatchBias float64
	tourLength   float64
	deviation    float64
}

func saveGksDeviationReport(path string, atspsData []AtspData, msaPatchBiases []float64) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows, err := buildGksDeviationRows(atspsData, msaPatchBiases)
	if err != nil {
		return err
	}

	var builder strings.Builder
	builder.WriteString("# GKS Patching Deviation\n\n")
	builder.WriteString("This table runs the cycle-cover patching heuristic directly, without MMAS. `MSA patch bias = 0.00` is the pure GKS-style cost-based patching variant; positive values bias patch selection toward MSA-supported connector edges.\n\n")
	writeGksDeviationTable(&builder, rows, msaPatchBiases)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildGksDeviationRows(atspsData []AtspData, msaPatchBiases []float64) ([]gksDeviationRow, error) {
	rows := make([]gksDeviationRow, 0, len(atspsData)*len(msaPatchBiases))
	for _, atspData := range atspsData {
		msaHeuristicMatrix, err := msaHeuristic.Read(atspData.msaHeuristicDirectoryPath)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read MSA heuristic: %w", atspData.name, err)
		}

		cycleCoverMatrix, _, err := buildMinimumCycleCoverMatrix(atspData.matrix)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to compute minimum cycle cover: %w", atspData.name, err)
		}

		for _, msaPatchBias := range msaPatchBiases {
			patchingMatrix := heuristics.BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(atspData.matrix, msaHeuristicMatrix, cycleCoverMatrix, msaPatchBias)
			tourLength, err := tourMatrixLength(atspData.matrix, patchingMatrix)
			if err != nil {
				return nil, fmt.Errorf("%s: invalid GKS tour for MSA patch bias %.2f: %w", atspData.name, msaPatchBias, err)
			}

			rows = append(rows, gksDeviationRow{
				instance:     atspData.name,
				msaPatchBias: msaPatchBias,
				tourLength:   tourLength,
				deviation:    100.0 * (tourLength - atspData.knownOptimal) / atspData.knownOptimal,
			})
		}
	}

	sort.Slice(rows, func(i, j int) bool {
		if rows[i].instance != rows[j].instance {
			return rows[i].instance < rows[j].instance
		}
		return rows[i].msaPatchBias < rows[j].msaPatchBias
	})

	return rows, nil
}

func writeGksDeviationTable(builder *strings.Builder, rows []gksDeviationRow, msaPatchBiases []float64) {
	bestBiasesByInstance := bestGksBiasesByInstance(rows)
	averageRows := averageGksDeviationRows(rows, msaPatchBiases)
	averageHighlights := lowestGksDeviationHighlights(averageRows)

	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>Instance</th><th>MSA patch bias</th><th>Tour length</th><th>Deviation [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")
	for _, row := range rows {
		key := gksDeviationRowKey(row.instance, row.msaPatchBias)
		writeGksDeviationTableRow(builder, row.instance, row, bestBiasesByInstance[key])
	}
	for _, row := range averageRows {
		key := gksDeviationRowKey(row.instance, row.msaPatchBias)
		writeGksDeviationTableRow(builder, "Average", row, averageHighlights[key])
	}
	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

func writeGksDeviationTableRow(builder *strings.Builder, instance string, row gksDeviationRow, highlightDeviation bool) {
	instanceCell := html.EscapeString(instance)
	if instance == "Average" {
		instanceCell = "<strong>Average</strong>"
	}

	tourLengthCell := fmt.Sprintf("%.2f", row.tourLength)
	if instance == "Average" {
		tourLengthCell = ""
	}

	fmt.Fprintf(builder,
		"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%s</td><td align=\"right\">%s</td></tr>\n",
		instanceCell,
		row.msaPatchBias,
		tourLengthCell,
		finalResultsSummaryMetricCell(row.deviation, highlightDeviation))
}

func bestGksBiasesByInstance(rows []gksDeviationRow) map[string]bool {
	bestByInstance := make(map[string]float64)
	for _, row := range rows {
		best, ok := bestByInstance[row.instance]
		if !ok || row.deviation < best {
			bestByInstance[row.instance] = row.deviation
		}
	}

	highlights := make(map[string]bool, len(rows))
	for _, row := range rows {
		if math.Abs(row.deviation-bestByInstance[row.instance]) < 1e-9 {
			highlights[gksDeviationRowKey(row.instance, row.msaPatchBias)] = true
		}
	}
	return highlights
}

func averageGksDeviationRows(rows []gksDeviationRow, msaPatchBiases []float64) []gksDeviationRow {
	totals := make(map[float64]float64, len(msaPatchBiases))
	counts := make(map[float64]int, len(msaPatchBiases))
	for _, row := range rows {
		totals[row.msaPatchBias] += row.deviation
		counts[row.msaPatchBias]++
	}

	averageRows := make([]gksDeviationRow, 0, len(msaPatchBiases))
	for _, msaPatchBias := range msaPatchBiases {
		count := counts[msaPatchBias]
		if count == 0 {
			continue
		}

		averageRows = append(averageRows, gksDeviationRow{
			instance:     "Average",
			msaPatchBias: msaPatchBias,
			deviation:    totals[msaPatchBias] / float64(count),
		})
	}
	return averageRows
}

func lowestGksDeviationHighlights(rows []gksDeviationRow) map[string]bool {
	best := math.Inf(1)
	for _, row := range rows {
		if row.deviation < best {
			best = row.deviation
		}
	}

	highlights := make(map[string]bool, len(rows))
	for _, row := range rows {
		if math.Abs(row.deviation-best) < 1e-9 {
			highlights[gksDeviationRowKey(row.instance, row.msaPatchBias)] = true
		}
	}
	return highlights
}

func gksDeviationRowKey(instance string, msaPatchBias float64) string {
	return instance + ":" + fmt.Sprintf("%.6f", msaPatchBias)
}

func tourMatrixLength(distanceMatrix, tourMatrix [][]float64) (float64, error) {
	dimension := len(tourMatrix)
	successors := make([]int, dimension)
	inDegree := make([]int, dimension)
	for vertex := range successors {
		successors[vertex] = -1
	}

	for from, row := range tourMatrix {
		outDegree := 0
		for to, value := range row {
			if value == 0 {
				continue
			}
			if from == to {
				return 0, fmt.Errorf("self-loop at vertex %d", from)
			}
			if from >= len(distanceMatrix) || to >= len(distanceMatrix[from]) {
				return 0, fmt.Errorf("edge %d -> %d is outside distance matrix", from, to)
			}

			outDegree++
			successors[from] = to
			inDegree[to]++
		}
		if outDegree != 1 {
			return 0, fmt.Errorf("vertex %d has out-degree %d", from, outDegree)
		}
	}

	for vertex, degree := range inDegree {
		if degree != 1 {
			return 0, fmt.Errorf("vertex %d has in-degree %d", vertex, degree)
		}
	}

	visited := make([]bool, dimension)
	current := 0
	for step := 0; step < dimension; step++ {
		if current < 0 || current >= dimension {
			return 0, fmt.Errorf("tour left matrix at vertex %d", current)
		}
		if visited[current] {
			return 0, fmt.Errorf("tour returned to vertex %d after %d steps", current, step)
		}
		visited[current] = true
		current = successors[current]
	}
	if current != 0 {
		return 0, fmt.Errorf("tour ended at vertex %d instead of 0", current)
	}

	length := 0.0
	for from, to := range successors {
		length += distanceMatrix[from][to]
	}
	return length, nil
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

func writeMsaHeuristicCycleCoverOverlapTable(builder *strings.Builder, rows []structuralComparison.InstanceAnalysis, totals msaHeuristicCycleCoverOverlapSummary) {
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

func writeMsaHeuristicCycleCoverOverlapRow(builder *strings.Builder, analysis structuralComparison.InstanceAnalysis) {
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

func msaHeuristicCycleCoverOverlapTotals(rows []structuralComparison.InstanceAnalysis) msaHeuristicCycleCoverOverlapSummary {
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

func saveMsaCountScalingReport(path string, atspsData []AtspData) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}

	rows, err := buildMsaCountScalingRows(atspsData, msaCountScalingCounts)
	if err != nil {
		return err
	}

	var builder strings.Builder
	builder.WriteString("# MSA Count Scaling\n\n")
	builder.WriteString("This table shows how the MSA heuristic signal changes when it is built from an increasing number of individual root MSAs. For each requested count, roots are selected deterministically and spread across the available root indices. If the requested count is larger than an instance dimension, all roots are used for that instance.\n\n")
	if len(rows) > 0 {
		fmt.Fprintf(&builder, "The table uses pooled totals over %d selected instances. Precision and recall use the %d instances with at least one found optimal tour in `solutions.csv`; boosted-edge density uses every selected instance. The boosted-edge target is `n - 1`, matching the edge count of a single arborescence.\n\n", rows[0].instanceCount, rows[0].referenceInstanceCount)
	}
	writeMsaCountScalingFindings(&builder, rows)
	builder.WriteString("\n")
	writeMsaCountScalingTable(&builder, rows)

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildMsaCountScalingRows(atspsData []AtspData, requestedCounts []int) ([]msaCountScalingRow, error) {
	rows := make([]msaCountScalingRow, 0, len(requestedCounts))
	for _, requestedCount := range requestedCounts {
		row := msaCountScalingRow{requestedCount: requestedCount}

		for _, atspData := range atspsData {
			msas, err := readOrCreateIndividualMsas(atspData)
			if err != nil {
				return nil, err
			}
			if len(msas) == 0 {
				continue
			}

			selectedRootCount := requestedCount
			if selectedRootCount == 0 || selectedRootCount > len(msas) {
				selectedRootCount = len(msas)
			}
			selectedRoots := selectEvenlySpacedRootIndexes(len(msas), selectedRootCount)
			boostedEdges := buildPartialMsaHeuristicEdgeSet(msas, selectedRoots)

			row.instanceCount++
			row.boostedEdges += len(boostedEdges)
			row.boostedTargetEdges += maxIntValue(len(msas)-1, 0)

			tours, err := msaHeuristicTours.ReadOptimalTours(atspData.optimalUniqueToursCsvPath)
			if err != nil {
				return nil, fmt.Errorf("%s: read found optimal tours: %w", atspData.name, err)
			}
			optimalEdges := buildAnalysisTourEdgeSet(tours)
			if len(optimalEdges) == 0 {
				continue
			}

			overlap := countEdgeSetIntersection(boostedEdges, optimalEdges)
			row.referenceInstanceCount++
			row.referenceBoostedEdges += len(boostedEdges)
			row.optimalEdges += len(optimalEdges)
			row.overlapEdges += overlap
		}

		rows = append(rows, row)
	}

	return rows, nil
}

func readOrCreateIndividualMsas(atspData AtspData) ([][][]float64, error) {
	msas, err := msaHeuristic.ReadMsas(atspData.msaHeuristicDirectoryPath)
	if err == nil && len(msas) == len(atspData.matrix) {
		return msas, nil
	}

	if _, createErr := msaHeuristic.Create(atspData.matrix, atspData.msaHeuristicDirectoryPath); createErr != nil {
		if err != nil {
			return nil, fmt.Errorf("%s: read individual MSAs: %w; create MSA heuristic: %w", atspData.name, err, createErr)
		}
		return nil, fmt.Errorf("%s: create MSA heuristic: %w", atspData.name, createErr)
	}

	msas, err = msaHeuristic.ReadMsas(atspData.msaHeuristicDirectoryPath)
	if err != nil {
		return nil, fmt.Errorf("%s: read individual MSAs: %w", atspData.name, err)
	}

	return msas, nil
}

func writeMsaCountScalingFindings(builder *strings.Builder, rows []msaCountScalingRow) {
	if len(rows) == 0 {
		return
	}

	first := rows[0]
	last := rows[len(rows)-1]
	fmt.Fprintf(builder, "## Findings\n\n")
	fmt.Fprintf(builder, "- **Precision changes from %.2f%% with %s MSA to %.2f%% with %s MSAs.**\n",
		100*first.precision(),
		first.label(),
		100*last.precision(),
		last.label())
	fmt.Fprintf(builder, "- **Recall changes from %.2f%% with %s MSA to %.2f%% with %s MSAs.**\n",
		100*first.recall(),
		first.label(),
		100*last.recall(),
		last.label())
}

func writeMsaCountScalingTable(builder *strings.Builder, rows []msaCountScalingRow) {
	builder.WriteString("<table>\n")
	builder.WriteString("<thead>\n")
	builder.WriteString("<tr><th>MSAs used</th><th>Boosted / n-1 [%]</th><th>Precision [%]</th><th>Recall [%]</th></tr>\n")
	builder.WriteString("</thead>\n")
	builder.WriteString("<tbody>\n")

	for _, row := range rows {
		fmt.Fprintf(builder,
			"<tr><td>%s</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td><td align=\"right\">%.2f</td></tr>\n",
			html.EscapeString(row.label()),
			100*row.boostedToTarget(),
			100*row.precision(),
			100*row.recall())
	}

	builder.WriteString("</tbody>\n")
	builder.WriteString("</table>\n")
}

type msaCountScalingRow struct {
	requestedCount         int
	instanceCount          int
	referenceInstanceCount int
	boostedEdges           int
	boostedTargetEdges     int
	referenceBoostedEdges  int
	optimalEdges           int
	overlapEdges           int
}

func (row msaCountScalingRow) label() string {
	if row.requestedCount == 0 {
		return "all"
	}

	return strconv.Itoa(row.requestedCount)
}

func (row msaCountScalingRow) boostedToTarget() float64 {
	return ratio(row.boostedEdges, row.boostedTargetEdges)
}

func (row msaCountScalingRow) precision() float64 {
	return ratio(row.overlapEdges, row.referenceBoostedEdges)
}

func (row msaCountScalingRow) recall() float64 {
	return ratio(row.overlapEdges, row.optimalEdges)
}

func buildPartialMsaHeuristicEdgeSet(msas [][][]float64, selectedRoots []int) map[models.Edge]struct{} {
	selectedRootSet := make(map[int]struct{}, len(selectedRoots))
	support := make(map[models.Edge]int)
	dimension := 0
	if len(msas) != 0 {
		dimension = len(msas[0])
	}

	for _, root := range selectedRoots {
		selectedRootSet[root] = struct{}{}
		for from := 0; from < len(msas[root]); from++ {
			for to, value := range msas[root][from] {
				if from == to || value <= 0 {
					continue
				}

				support[models.Edge{From: from, To: to}]++
			}
		}
	}

	edges := make(map[models.Edge]struct{})
	for edge, count := range support {
		eligibleRoots := len(selectedRoots)
		if _, ok := selectedRootSet[edge.To]; ok {
			eligibleRoots--
		}
		if eligibleRoots > 0 && count == eligibleRoots && edge.From >= 0 && edge.From < dimension && edge.To >= 0 && edge.To < dimension {
			edges[edge] = struct{}{}
		}
	}

	return edges
}

func selectEvenlySpacedRootIndexes(total, count int) []int {
	if count <= 0 || count >= total {
		roots := make([]int, total)
		for i := range roots {
			roots[i] = i
		}
		return roots
	}
	if count == 1 {
		return []int{0}
	}

	roots := make([]int, count)
	for i := 0; i < count; i++ {
		roots[i] = i * (total - 1) / (count - 1)
	}

	return roots
}

func buildAnalysisTourEdgeSet(tours map[string][]int) map[models.Edge]struct{} {
	edges := make(map[models.Edge]struct{})
	tourIds := make([]string, 0, len(tours))
	for tourId := range tours {
		tourIds = append(tourIds, tourId)
	}
	sort.Strings(tourIds)

	for _, tourId := range tourIds {
		for _, edge := range models.ConvertTourToEdges(tours[tourId]) {
			edges[edge] = struct{}{}
		}
	}

	return edges
}

func countEdgeSetIntersection(left, right map[models.Edge]struct{}) int {
	count := 0
	for edge := range left {
		if _, ok := right[edge]; ok {
			count++
		}
	}

	return count
}

func maxIntValue(left, right int) int {
	if left > right {
		return left
	}

	return right
}

func sortedStructuralAnalyses(analyses []structuralComparison.InstanceAnalysis) []structuralComparison.InstanceAnalysis {
	rows := append([]structuralComparison.InstanceAnalysis(nil), analyses...)
	sort.SliceStable(rows, func(i, j int) bool {
		return rows[i].Instance < rows[j].Instance
	})

	return rows
}

func filterAnalysesWithFoundOptimalEdges(analyses []structuralComparison.InstanceAnalysis) []structuralComparison.InstanceAnalysis {
	rows := make([]structuralComparison.InstanceAnalysis, 0, len(analyses))
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
