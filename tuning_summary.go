package main

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
	tuningSummaryReportFileName = "tuning_summary.md"
	tuningNearBestTolerance     = 0.25
)

type tuningSummaryKey struct {
	heuristicWeight float64
	msaPatchBias    float64
}

type tuningSummaryAggregate struct {
	key            tuningSummaryKey
	results        []float64
	successRates   []float64
	exactBestCount int
	nearBestCount  int
}

type tuningSummaryRow struct {
	key            tuningSummaryKey
	meanResult     float64
	medianResult   float64
	meanSuccess    float64
	exactBestCount int
	nearBestCount  int
}

func saveTuningSummary(resultsRootPath string, atspsData []AtspData, heuristics []string) error {
	if err := removeLegacyBestParametersReports(resultsRootPath); err != nil {
		return err
	}

	reportPath := filepath.Join(resultsRootPath, tuningSummaryReportFileName)
	if err := os.MkdirAll(filepath.Dir(reportPath), 0700); err != nil {
		return err
	}

	var builder strings.Builder
	builder.WriteString("# Tuning Summary\n\n")
	builder.WriteString("Result means average best deviation from `result.csv`; lower is better. ")
	builder.WriteString(fmt.Sprintf("Near-best means within %.2f percentage points of the best result for the same instance.\n\n", tuningNearBestTolerance))

	for _, heuristic := range heuristics {
		if heuristicIsSparseControl(heuristic) {
			continue
		}

		if err := appendTuningHeuristicSummary(&builder, atspsData, heuristic); err != nil {
			return err
		}
	}

	return os.WriteFile(reportPath, []byte(builder.String()), 0644)
}

func appendTuningHeuristicSummary(builder *strings.Builder, atspsData []AtspData, heuristic string) error {
	rows, instanceCount, missingInstances, err := buildTuningSummaryRows(atspsData, heuristic)
	if err != nil {
		return err
	}

	fmt.Fprintf(builder, "## %s\n\n", heuristicDisplayName(heuristic))
	fmt.Fprintf(builder, "Instances used: %d/%d\n\n", instanceCount, len(atspsData))
	if len(missingInstances) != 0 {
		fmt.Fprintf(builder, "Missing result files: %s\n\n", strings.Join(missingInstances, ", "))
	}
	if len(rows) == 0 {
		builder.WriteString("No tuning results found.\n\n")
		return nil
	}

	includeMsaPatchBias := heuristic == heuristicCycleCoverMsaPatching
	if includeMsaPatchBias {
		builder.WriteString("| Heuristic weight | MSA patch bias | Mean result [%] | Median result [%] | Exact best | Near best | Mean success [%] |\n")
		builder.WriteString("|------------------|----------------|-----------------|-------------------|------------|-----------|------------------|\n")
	} else {
		builder.WriteString("| Heuristic weight | Mean result [%] | Median result [%] | Exact best | Near best | Mean success [%] |\n")
		builder.WriteString("|------------------|-----------------|-------------------|------------|-----------|------------------|\n")
	}

	for _, row := range rows {
		if includeMsaPatchBias {
			fmt.Fprintf(builder, "| %.2f | %.2f | %.2f | %.2f | %d | %d | %.2f |\n",
				row.key.heuristicWeight,
				row.key.msaPatchBias,
				row.meanResult,
				row.medianResult,
				row.exactBestCount,
				row.nearBestCount,
				row.meanSuccess)
		} else {
			fmt.Fprintf(builder, "| %.2f | %.2f | %.2f | %d | %d | %.2f |\n",
				row.key.heuristicWeight,
				row.meanResult,
				row.medianResult,
				row.exactBestCount,
				row.nearBestCount,
				row.meanSuccess)
		}
	}

	builder.WriteString("\n")
	return nil
}

func buildTuningSummaryRows(atspsData []AtspData, heuristic string) ([]tuningSummaryRow, int, []string, error) {
	aggregates := map[tuningSummaryKey]*tuningSummaryAggregate{}
	missingInstances := make([]string, 0)
	instanceCount := 0

	for _, atspData := range atspsData {
		statistics, err := readStatistics(resultFilePathForHeuristic(atspData, heuristic))
		if err != nil {
			if errors.Is(err, os.ErrNotExist) {
				missingInstances = append(missingInstances, atspData.name)
				continue
			}
			return nil, 0, nil, err
		}
		if len(statistics) == 0 {
			missingInstances = append(missingInstances, atspData.name)
			continue
		}

		instanceCount++
		instanceBestResult := bestTuningResult(statistics)
		for _, statistic := range statistics {
			key := tuningSummaryKey{
				heuristicWeight: statistic.heuristicWeight,
				msaPatchBias:    statistic.msaPatchBias,
			}
			if heuristic != heuristicCycleCoverMsaPatching {
				key.msaPatchBias = 0
			}

			aggregate, ok := aggregates[key]
			if !ok {
				aggregate = &tuningSummaryAggregate{key: key}
				aggregates[key] = aggregate
			}

			aggregate.results = append(aggregate.results, statistic.averageBestDeviation)
			aggregate.successRates = append(aggregate.successRates, statistic.successRate)
			if tuningFloatEqual(statistic.averageBestDeviation, instanceBestResult) {
				aggregate.exactBestCount++
			}
			if statistic.averageBestDeviation <= instanceBestResult+tuningNearBestTolerance {
				aggregate.nearBestCount++
			}
		}
	}

	sort.Strings(missingInstances)
	rows := make([]tuningSummaryRow, 0, len(aggregates))
	for _, aggregate := range aggregates {
		rows = append(rows, tuningSummaryRow{
			key:            aggregate.key,
			meanResult:     averageFloat64(aggregate.results),
			medianResult:   medianFloat64(aggregate.results),
			meanSuccess:    averageFloat64(aggregate.successRates),
			exactBestCount: aggregate.exactBestCount,
			nearBestCount:  aggregate.nearBestCount,
		})
	}

	sortTuningSummaryRows(rows)
	return rows, instanceCount, missingInstances, nil
}

func bestTuningResult(statistics []ExperimentsDataStatistics) float64 {
	best := math.Inf(1)
	for _, statistic := range statistics {
		if statistic.averageBestDeviation < best {
			best = statistic.averageBestDeviation
		}
	}
	return best
}

func sortTuningSummaryRows(rows []tuningSummaryRow) {
	sort.SliceStable(rows, func(i, j int) bool {
		left, right := rows[i], rows[j]
		if !tuningFloatEqual(left.meanResult, right.meanResult) {
			return left.meanResult < right.meanResult
		}
		if !tuningFloatEqual(left.medianResult, right.medianResult) {
			return left.medianResult < right.medianResult
		}
		if left.exactBestCount != right.exactBestCount {
			return left.exactBestCount > right.exactBestCount
		}
		if left.nearBestCount != right.nearBestCount {
			return left.nearBestCount > right.nearBestCount
		}
		if !tuningFloatEqual(left.meanSuccess, right.meanSuccess) {
			return left.meanSuccess > right.meanSuccess
		}
		if !tuningFloatEqual(left.key.heuristicWeight, right.key.heuristicWeight) {
			return left.key.heuristicWeight < right.key.heuristicWeight
		}
		return left.key.msaPatchBias < right.key.msaPatchBias
	})
}

func averageFloat64(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}

	sum := 0.0
	for _, value := range values {
		sum += value
	}
	return sum / float64(len(values))
}

func medianFloat64(values []float64) float64 {
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

func tuningFloatEqual(left, right float64) bool {
	return math.Abs(left-right) < 1e-9
}

func removeLegacyBestParametersReports(resultsRootPath string) error {
	matches, err := filepath.Glob(filepath.Join(resultsRootPath, "best_parameters_report*.md"))
	if err != nil {
		return err
	}

	for _, match := range matches {
		if err := removeFileIfExists(match); err != nil {
			return err
		}
	}

	return nil
}
