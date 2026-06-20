package app

import (
	"atsp_aco_msa/modules/analysis/tours"
	"atsp_aco_msa/modules/models"
	"fmt"
	"html"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

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

			optimalTours, err := tours.ReadOptimal(atspData.OptimalUniqueToursCsvPath)
			if err != nil {
				return nil, fmt.Errorf("%s: read found optimal tours: %w", atspData.Name, err)
			}
			optimalEdges := buildAnalysisTourEdgeSet(optimalTours)
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
