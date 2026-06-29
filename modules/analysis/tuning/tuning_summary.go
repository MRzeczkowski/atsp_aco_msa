package tuning

import (
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

const (
	ReportFileName = "tuning_summary.md"
)

type Statistic struct {
	HeuristicWeight      float64
	MsaPatchBias         float64
	AverageBestDeviation float64
	SuccessRate          float64
}

type ResultFile struct {
	Instance string
	Path     string
}

type HeuristicConfig struct {
	Name                string
	DisplayName         string
	IncludeMsaPatchBias bool
	ResultFiles         []ResultFile
}

type Config struct {
	ResultsRootPath string
	Heuristics      []HeuristicConfig
	ReadStatistics  func(path string) ([]Statistic, error)
}

type summaryKey struct {
	heuristicWeight float64
	msaPatchBias    float64
}

type summaryAggregate struct {
	key          summaryKey
	results      []float64
	successRates []float64
}

type summaryRow struct {
	key          summaryKey
	meanResult   float64
	medianResult float64
	meanSuccess  float64
}

func Save(config Config) error {
	if config.ReadStatistics == nil {
		return fmt.Errorf("statistics reader is required")
	}

	reportPath := filepath.Join(config.ResultsRootPath, ReportFileName)
	if err := os.MkdirAll(filepath.Dir(reportPath), 0700); err != nil {
		return err
	}

	var builder strings.Builder
	builder.WriteString("# Tuning Summary\n\n")
	builder.WriteString("Result means average best deviation from `result.csv`; lower is better.\n\n")

	for _, heuristic := range config.Heuristics {
		if err := appendHeuristicSummary(&builder, heuristic, config.ReadStatistics); err != nil {
			return err
		}
	}

	return os.WriteFile(reportPath, []byte(builder.String()), 0644)
}

func appendHeuristicSummary(builder *strings.Builder, heuristic HeuristicConfig, readStatistics func(path string) ([]Statistic, error)) error {
	rows, instanceCount, missingInstances, err := buildSummaryRows(heuristic, readStatistics)
	if err != nil {
		return err
	}

	displayName := heuristic.DisplayName
	if displayName == "" {
		displayName = heuristic.Name
	}

	fmt.Fprintf(builder, "## %s\n\n", displayName)
	fmt.Fprintf(builder, "Instances used: %d/%d\n\n", instanceCount, len(heuristic.ResultFiles))
	if len(missingInstances) != 0 {
		fmt.Fprintf(builder, "Missing result files: %s\n\n", strings.Join(missingInstances, ", "))
	}
	if len(rows) == 0 {
		builder.WriteString("No tuning results found.\n\n")
		return nil
	}

	if heuristic.IncludeMsaPatchBias {
		builder.WriteString("| Heuristic weight | MSA patch bias | Mean result [%] | Median result [%] | Mean success [%] |\n")
		builder.WriteString("|------------------|----------------|-----------------|-------------------|------------------|\n")
	} else {
		builder.WriteString("| Heuristic weight | Mean result [%] | Median result [%] | Mean success [%] |\n")
		builder.WriteString("|------------------|-----------------|-------------------|------------------|\n")
	}

	for _, row := range rows {
		if heuristic.IncludeMsaPatchBias {
			fmt.Fprintf(builder, "| %.2f | %.2f | %.2f | %.2f | %.2f |\n",
				row.key.heuristicWeight,
				row.key.msaPatchBias,
				row.meanResult,
				row.medianResult,
				row.meanSuccess)
		} else {
			fmt.Fprintf(builder, "| %.2f | %.2f | %.2f | %.2f |\n",
				row.key.heuristicWeight,
				row.meanResult,
				row.medianResult,
				row.meanSuccess)
		}
	}

	builder.WriteString("\n")
	return nil
}

func buildSummaryRows(heuristic HeuristicConfig, readStatistics func(path string) ([]Statistic, error)) ([]summaryRow, int, []string, error) {
	aggregates := map[summaryKey]*summaryAggregate{}
	missingInstances := make([]string, 0)
	instanceCount := 0

	for _, resultFile := range heuristic.ResultFiles {
		statistics, err := readStatistics(resultFile.Path)
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingInstances = append(missingInstances, resultFile.Instance)
				continue
			}
			return nil, 0, nil, err
		}
		if len(statistics) == 0 {
			missingInstances = append(missingInstances, resultFile.Instance)
			continue
		}

		instanceCount++
		for _, statistic := range statistics {
			key := summaryKey{
				heuristicWeight: statistic.HeuristicWeight,
				msaPatchBias:    statistic.MsaPatchBias,
			}
			if !heuristic.IncludeMsaPatchBias {
				key.msaPatchBias = 0
			}

			aggregate, ok := aggregates[key]
			if !ok {
				aggregate = &summaryAggregate{key: key}
				aggregates[key] = aggregate
			}

			aggregate.results = append(aggregate.results, statistic.AverageBestDeviation)
			aggregate.successRates = append(aggregate.successRates, statistic.SuccessRate)
		}
	}

	sort.Strings(missingInstances)
	rows := make([]summaryRow, 0, len(aggregates))
	for _, aggregate := range aggregates {
		rows = append(rows, summaryRow{
			key:          aggregate.key,
			meanResult:   average(aggregate.results),
			medianResult: median(aggregate.results),
			meanSuccess:  average(aggregate.successRates),
		})
	}

	sortSummaryRows(rows)
	return rows, instanceCount, missingInstances, nil
}

func sortSummaryRows(rows []summaryRow) {
	sort.SliceStable(rows, func(i, j int) bool {
		left, right := rows[i], rows[j]
		if !floatEqual(left.meanResult, right.meanResult) {
			return left.meanResult < right.meanResult
		}
		if !floatEqual(left.medianResult, right.medianResult) {
			return left.medianResult < right.medianResult
		}
		if !floatEqual(left.meanSuccess, right.meanSuccess) {
			return left.meanSuccess > right.meanSuccess
		}
		if !floatEqual(left.key.heuristicWeight, right.key.heuristicWeight) {
			return left.key.heuristicWeight < right.key.heuristicWeight
		}
		return left.key.msaPatchBias < right.key.msaPatchBias
	})
}

func average(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}

	sum := 0.0
	for _, value := range values {
		sum += value
	}
	return sum / float64(len(values))
}

func median(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}

	sortedValues := append([]float64(nil), values...)
	sort.Float64s(sortedValues)
	middle := len(sortedValues) / 2
	if len(sortedValues)%2 == 1 {
		return sortedValues[middle]
	}

	return (sortedValues[middle-1] + sortedValues[middle]) / 2
}

func floatEqual(left, right float64) bool {
	return math.Abs(left-right) < 1e-9
}
