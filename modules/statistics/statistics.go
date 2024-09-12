package statistics

import (
	"atsp_aco_msa/modules/models"
	"math"
)

type Edge = models.Edge

type GraphStats struct {
	minWeight, maxWeight, avgWeight, stdDevWeight, skewness, kurtosis float64
}

func CalculateEdgesStats(edges []Edge, weights map[Edge]float64) (stats GraphStats) {
	var totalWeight float64
	var count int
	minWeight := math.MaxFloat64
	maxWeight := -math.MaxFloat64

	for _, edge := range edges {
		weight := weights[edge]
		totalWeight += weight

		if weight < minWeight {
			minWeight = weight
		}

		if weight > maxWeight {
			maxWeight = weight
		}

		count++
	}

	avgWeight := totalWeight / float64(count)

	var sumOfSquares, sumOfCubes, sumOfFourthPowers float64
	for _, edge := range edges {
		diff := weights[edge] - avgWeight
		square := diff * diff
		cube := square * diff
		fourthPower := cube * diff

		sumOfSquares += square
		sumOfCubes += cube
		sumOfFourthPowers += fourthPower
	}

	n := float64(count)
	stdDevWeight := math.Sqrt(sumOfSquares / n)
	skewness := (sumOfCubes / n) / math.Pow(stdDevWeight, 3)
	kurtosis := (sumOfFourthPowers/n)/math.Pow(stdDevWeight, 4) - 3

	return GraphStats{
		minWeight:    minWeight,
		maxWeight:    maxWeight,
		avgWeight:    avgWeight,
		stdDevWeight: stdDevWeight,
		skewness:     skewness,
		kurtosis:     kurtosis,
	}
}

func CalculateMatrixStats(matrix [][]float64) (stats GraphStats) {
	var totalWeight float64
	var count int
	minWeight := math.MaxFloat64
	maxWeight := -math.MaxFloat64

	for i := range matrix {
		for j := range matrix[i] {
			if i != j {
				weight := matrix[i][j]
				totalWeight += weight

				if weight < minWeight {
					minWeight = weight
				}

				if weight > maxWeight {
					maxWeight = weight
				}

				count++
			}
		}
	}

	avgWeight := totalWeight / float64(count)

	var sumOfSquares, sumOfCubes, sumOfFourthPowers float64
	for i := range matrix {
		for j := range matrix[i] {
			if i != j {
				diff := matrix[i][j] - avgWeight
				square := diff * diff
				cube := square * diff
				fourthPower := cube * diff

				sumOfSquares += square
				sumOfCubes += cube
				sumOfFourthPowers += fourthPower
			}
		}
	}

	n := float64(count)
	stdDevWeight := math.Sqrt(sumOfSquares / n)
	skewness := (sumOfCubes / n) / math.Pow(stdDevWeight, 3)
	kurtosis := (sumOfFourthPowers/n)/math.Pow(stdDevWeight, 4) - 3

	return GraphStats{
		minWeight:    minWeight,
		maxWeight:    maxWeight,
		avgWeight:    avgWeight,
		stdDevWeight: stdDevWeight,
		skewness:     skewness,
		kurtosis:     kurtosis,
	}
}

// Function to count leaves in a tree represented by edges
func CountLeaves(edges []Edge) float64 {
	if len(edges) == 0 {
		return 0
	}

	// Create a map to hold child nodes for each node
	children := make(map[int][]int)
	allNodes := make(map[int]bool)

	// Populate the children map
	for _, edge := range edges {
		children[edge.From] = append(children[edge.From], edge.To)
		allNodes[edge.From] = true
		allNodes[edge.To] = true
	}

	leafCount := 0.0
	for node := range allNodes {
		if len(children[node]) == 0 {
			leafCount++
		}
	}

	return leafCount
}
