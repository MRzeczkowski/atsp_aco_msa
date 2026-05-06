package hungarian

import (
	"errors"
	"fmt"
	"math"
)

// Solve returns a minimum-cost assignment for a square non-negative cost matrix.
// Costs must be finite. The returned assignment maps each row i to column
// assignment[i].
func Solve(costs [][]float64) ([]int, float64, error) {
	if err := validateCostMatrix(costs); err != nil {
		return nil, 0, err
	}

	n := len(costs)

	// This is the O(n^3) primal-dual Hungarian algorithm. It scans rows and
	// columns in index order, so ties are resolved deterministically.
	rowPotential := make([]float64, n+1)
	columnPotential := make([]float64, n+1)
	matchedRowByColumn := make([]int, n+1)
	previousColumn := make([]int, n+1)

	for row := 1; row <= n; row++ {
		matchedRowByColumn[0] = row
		currentColumn := 0

		minReducedCost := make([]float64, n+1)
		usedColumn := make([]bool, n+1)
		for i := range minReducedCost {
			minReducedCost[i] = math.Inf(1)
		}

		for {
			usedColumn[currentColumn] = true
			currentRow := matchedRowByColumn[currentColumn]
			delta := math.Inf(1)
			nextColumn := 0

			for column := 1; column <= n; column++ {
				if usedColumn[column] {
					continue
				}

				reducedCost := costs[currentRow-1][column-1] - rowPotential[currentRow] - columnPotential[column]
				if reducedCost < minReducedCost[column] {
					minReducedCost[column] = reducedCost
					previousColumn[column] = currentColumn
				}
				if minReducedCost[column] < delta {
					delta = minReducedCost[column]
					nextColumn = column
				}
			}

			if nextColumn == 0 || math.IsInf(delta, 1) {
				return nil, 0, errors.New("could not find an augmenting path")
			}

			for column := 0; column <= n; column++ {
				if usedColumn[column] {
					rowPotential[matchedRowByColumn[column]] += delta
					columnPotential[column] -= delta
				} else {
					minReducedCost[column] -= delta
				}
			}

			currentColumn = nextColumn
			if matchedRowByColumn[currentColumn] == 0 {
				break
			}
		}

		for {
			nextColumn := previousColumn[currentColumn]
			matchedRowByColumn[currentColumn] = matchedRowByColumn[nextColumn]
			currentColumn = nextColumn
			if currentColumn == 0 {
				break
			}
		}
	}

	assignment := make([]int, n)
	for column := 1; column <= n; column++ {
		row := matchedRowByColumn[column]
		if row == 0 {
			continue
		}
		assignment[row-1] = column - 1
	}
	if err := validateAssignment(assignment); err != nil {
		return nil, 0, fmt.Errorf("constructed invalid assignment: %w", err)
	}

	totalCost := AssignmentCost(costs, assignment)
	if math.IsInf(totalCost, 0) || math.IsNaN(totalCost) {
		return nil, 0, errors.New("assignment has invalid cost")
	}

	return assignment, totalCost, nil
}

// AssignmentCost returns the total cost of an assignment.
func AssignmentCost(costs [][]float64, assignment []int) float64 {
	total := 0.0
	for row, column := range assignment {
		total += costs[row][column]
	}
	return total
}

func validateCostMatrix(costs [][]float64) error {
	n := len(costs)
	if n == 0 {
		return errors.New("cost matrix is empty")
	}

	for row := 0; row < n; row++ {
		if len(costs[row]) != n {
			return fmt.Errorf("cost matrix must be square: row %d has length %d, expected %d", row, len(costs[row]), n)
		}
		for column := 0; column < n; column++ {
			cost := costs[row][column]
			if math.IsNaN(cost) {
				return fmt.Errorf("cost matrix contains NaN at (%d,%d)", row, column)
			}
			if math.IsInf(cost, 0) {
				return fmt.Errorf("cost matrix contains infinity at (%d,%d)", row, column)
			}
			if cost < 0 {
				return fmt.Errorf("cost matrix contains negative value at (%d,%d): %v", row, column, cost)
			}
		}
	}

	return nil
}
