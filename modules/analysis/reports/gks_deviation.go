package reports

import (
	"atsp_aco_msa/modules/algorithms/heuristics"
	"atsp_aco_msa/modules/artifacts/cyclecover"
	"atsp_aco_msa/modules/artifacts/msaheuristic"
	"atsp_aco_msa/modules/project"
	"fmt"
	"html"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

type gksDeviationRow struct {
	instance     string
	msaPatchBias float64
	tourLength   float64
	deviation    float64
}

func SaveGksDeviation(path string, atspsData []project.AtspData, msaPatchBiases []float64) error {
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

func buildGksDeviationRows(atspsData []project.AtspData, msaPatchBiases []float64) ([]gksDeviationRow, error) {
	rows := make([]gksDeviationRow, 0, len(atspsData)*len(msaPatchBiases))
	for _, atspData := range atspsData {
		msaHeuristicMatrix, err := msaheuristic.Read(atspData.MsaHeuristicDirectoryPath)
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read MSA Heuristic: %w", atspData.Name, err)
		}

		cycleCoverMatrix, err := cyclecover.Read(atspData.CycleCoverDirectoryPath, len(atspData.Matrix))
		if err != nil {
			return nil, fmt.Errorf("%s: failed to read cycle-cover cache: %w", atspData.Name, err)
		}

		for _, msaPatchBias := range msaPatchBiases {
			patchingMatrix := heuristics.BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(atspData.Matrix, msaHeuristicMatrix, cycleCoverMatrix, msaPatchBias)
			tourLength, err := tourMatrixLength(atspData.Matrix, patchingMatrix)
			if err != nil {
				return nil, fmt.Errorf("%s: invalid GKS tour for MSA patch bias %.2f: %w", atspData.Name, msaPatchBias, err)
			}

			rows = append(rows, gksDeviationRow{
				instance:     atspData.Name,
				msaPatchBias: msaPatchBias,
				tourLength:   tourLength,
				deviation:    100.0 * (tourLength - atspData.KnownOptimal) / atspData.KnownOptimal,
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
		metricCell(row.deviation, highlightDeviation))
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
