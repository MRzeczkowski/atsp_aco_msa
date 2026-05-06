package hungarian

import (
	"atsp_aco_msa/modules/models"
	"errors"
	"fmt"
	"math"
)

// MinimumCycleCover returns a minimum-cost directed cycle cover for an ATSP
// distance matrix. Self-loops are forbidden even when the diagonal is cheap.
func MinimumCycleCover(distances [][]float64) ([]models.Edge, float64, error) {
	if err := validateCostMatrix(distances); err != nil {
		return nil, 0, err
	}

	n := len(distances)
	if n < 2 {
		return nil, 0, errors.New("minimum cycle cover requires at least two vertices")
	}

	costs := make([][]float64, n)
	for i := range costs {
		costs[i] = make([]float64, n)
		copy(costs[i], distances[i])
	}

	forbiddenCost, err := forbiddenSelfLoopCost(distances)
	if err != nil {
		return nil, 0, err
	}

	for i := 0; i < n; i++ {
		costs[i][i] = forbiddenCost
	}

	assignment, _, err := Solve(costs)
	if err != nil {
		return nil, 0, err
	}

	edges := make([]models.Edge, n)
	totalCost := 0.0
	for from, to := range assignment {
		if from == to {
			return nil, 0, errors.New("minimum cycle cover selected a forbidden self-loop")
		}

		edges[from] = models.Edge{From: from, To: to}
		totalCost += distances[from][to]
	}

	return edges, totalCost, nil
}

// AssignmentCycles decomposes an assignment into directed cycles. Cycles are
// returned in deterministic order by their smallest unvisited start vertex.
func AssignmentCycles(assignment []int) ([][]int, error) {
	if err := validateAssignment(assignment); err != nil {
		return nil, err
	}

	visited := make([]bool, len(assignment))
	cycles := make([][]int, 0)

	for start := range assignment {
		if visited[start] {
			continue
		}

		cycle := make([]int, 0)
		current := start
		for !visited[current] {
			visited[current] = true
			cycle = append(cycle, current)
			current = assignment[current]
		}

		cycles = append(cycles, cycle)
	}

	return cycles, nil
}

func forbiddenSelfLoopCost(distances [][]float64) (float64, error) {
	n := len(distances)
	maxOffDiagonal := 0.0

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			if distances[i][j] > maxOffDiagonal {
				maxOffDiagonal = distances[i][j]
			}
		}
	}

	upperBound := float64(n) * maxOffDiagonal
	if math.IsInf(upperBound, 0) || math.IsNaN(upperBound) {
		return 0, errors.New("could not compute finite forbidden self-loop cost")
	}
	if upperBound == 0 {
		return 1, nil
	}

	forbiddenCost := math.Nextafter(upperBound, math.Inf(1))
	if math.IsInf(forbiddenCost, 0) || math.IsNaN(forbiddenCost) {
		return 0, errors.New("could not compute finite forbidden self-loop cost")
	}

	return forbiddenCost, nil
}

func validateAssignment(assignment []int) error {
	n := len(assignment)
	if n == 0 {
		return errors.New("assignment is empty")
	}

	seen := make([]bool, n)
	for row, column := range assignment {
		if column < 0 || column >= n {
			return fmt.Errorf("assignment[%d] is out of range: %d", row, column)
		}
		if seen[column] {
			return fmt.Errorf("assignment is not one-to-one: column %d appears more than once", column)
		}
		seen[column] = true
	}

	return nil
}
