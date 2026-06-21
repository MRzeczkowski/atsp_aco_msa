package reports

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
)

type EvaluationResultSummaryMetric struct {
	AverageMinDeviation  float64
	SuccessRate          float64
	AverageBestIteration float64
	HeuristicWeight      float64
	Iterations           int
}

func ReadEvaluationResultSummaryMetrics(resultCsvPath string) (map[string]EvaluationResultSummaryMetric, error) {
	file, err := os.Open(resultCsvPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		return nil, err
	}

	heuristicIndex := indexOf(header, "Heuristic")
	heuristicWeightIndex := indexOf(header, "Heuristic weight")
	iterationsIndex := indexOf(header, "Iterations")
	averageBestIterationIndex := indexOf(header, "Avg best at iteration")
	averageDeviationIndex := indexOf(header, "Avg best deviation")
	successRateIndex := indexOf(header, "Success rate [%]")
	if heuristicIndex == -1 || heuristicWeightIndex == -1 || iterationsIndex == -1 || averageBestIterationIndex == -1 || averageDeviationIndex == -1 || successRateIndex == -1 {
		return nil, fmt.Errorf("missing required summary columns")
	}

	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	metrics := make(map[string]EvaluationResultSummaryMetric, len(records))
	for _, record := range records {
		if len(record) != len(header) {
			return nil, fmt.Errorf("invalid record length: got %d want %d", len(record), len(header))
		}

		heuristicWeight, err := strconv.ParseFloat(record[heuristicWeightIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid heuristic weight for %s: %w", record[heuristicIndex], err)
		}
		iterations, err := strconv.Atoi(record[iterationsIndex])
		if err != nil {
			return nil, fmt.Errorf("invalid iterations for %s: %w", record[heuristicIndex], err)
		}
		averageBestIteration, err := strconv.ParseFloat(record[averageBestIterationIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid average best iteration for %s: %w", record[heuristicIndex], err)
		}
		averageDeviation, err := strconv.ParseFloat(record[averageDeviationIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid average deviation for %s: %w", record[heuristicIndex], err)
		}
		successRate, err := strconv.ParseFloat(record[successRateIndex], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid success rate for %s: %w", record[heuristicIndex], err)
		}

		metric := EvaluationResultSummaryMetric{
			AverageMinDeviation:  averageDeviation,
			SuccessRate:          successRate,
			AverageBestIteration: averageBestIteration,
			HeuristicWeight:      heuristicWeight,
			Iterations:           iterations,
		}
		heuristic := record[heuristicIndex]
		if current, ok := metrics[heuristic]; !ok || evaluationResultSummaryMetricIsBetter(metric, current) {
			metrics[heuristic] = metric
		}
	}

	return metrics, nil
}

func evaluationResultSummaryMetricIsBetter(candidate, current EvaluationResultSummaryMetric) bool {
	if candidate.AverageMinDeviation != current.AverageMinDeviation {
		return candidate.AverageMinDeviation < current.AverageMinDeviation
	}
	if candidate.SuccessRate != current.SuccessRate {
		return candidate.SuccessRate > current.SuccessRate
	}
	if candidate.AverageBestIteration != current.AverageBestIteration {
		return candidate.AverageBestIteration < current.AverageBestIteration
	}
	return candidate.HeuristicWeight < current.HeuristicWeight
}

func indexOf(values []string, value string) int {
	for i, candidate := range values {
		if candidate == value {
			return i
		}
	}

	return -1
}
