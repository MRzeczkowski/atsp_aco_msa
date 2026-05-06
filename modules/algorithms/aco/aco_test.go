package aco

import (
	"math"
	"testing"
)

func TestRunKeepsSelfLoopDesirabilityZero(t *testing.T) {
	distances := [][]float64{
		{0, 1, 2},
		{3, 0, 4},
		{5, 6, 0},
	}
	hints := [][]float64{
		{0, 0, 0},
		{0, 0, 0},
		{0, 0, 0},
	}

	aco := NewACO(1.0, 2.0, 0.8, 0.0, 1, 100.0, distances, hints)
	aco.Run()

	for i := range distances {
		if aco.pheromones[i][i] != 0 {
			t.Fatalf("expected self-loop pheromone at %d,%d to remain zero, got %f", i, i, aco.pheromones[i][i])
		}
		if aco.heuristics[i][i] != 0 {
			t.Fatalf("expected self-loop heuristic at %d,%d to remain zero, got %f", i, i, aco.heuristics[i][i])
		}
		if aco.desirabilities[i][i] != 0 {
			t.Fatalf("expected self-loop desirability at %d,%d to remain zero, got %f", i, i, aco.desirabilities[i][i])
		}
	}
}

func TestSelectNextCityUsesDeterministicFallback(t *testing.T) {
	aco := &ACO{
		neighborsLists: [][]int{
			{1},
		},
		desirabilities: [][]float64{
			{0, 0, 2, 5},
		},
	}
	canVisitBits := []float64{0, 0, 1, 1}
	desirabilities := make([]float64, 4)

	nextCity := aco.selectNextCity(0, canVisitBits, desirabilities)

	if nextCity != 3 {
		t.Fatalf("expected deterministic fallback to choose city 3, got %d", nextCity)
	}
}

func TestNewACOKeepsTourConstructionNeighborsDistanceOnly(t *testing.T) {
	distances := make([][]float64, 12)
	hints := make([][]float64, 12)
	for i := range distances {
		distances[i] = make([]float64, 12)
		hints[i] = make([]float64, 12)

		for j := range distances[i] {
			distances[i][j] = 1.0
		}
	}

	for j := 1; j < 12; j++ {
		distances[0][j] = float64(j)
	}
	hints[0][11] = 11.0

	withoutCmsa := NewACO(1.0, 1.0, 0.8, 0.0, 1, 100.0, distances, hints)
	withCmsa := NewACO(1.0, 1.0, 0.8, 1.0, 1, 100.0, distances, hints)

	if !equalIntSlices(withoutCmsa.neighborsLists[0], withCmsa.neighborsLists[0]) {
		t.Fatalf("expected pCmsa not to change candidate list, got without=%v with=%v", withoutCmsa.neighborsLists[0], withCmsa.neighborsLists[0])
	}

	for _, neighbor := range withCmsa.neighborsLists[0] {
		if neighbor == 11 {
			t.Fatalf("expected distance-only candidate list to exclude far CMSA edge, got %v", withCmsa.neighborsLists[0])
		}
	}
}

func equalIntSlices(left, right []int) bool {
	if len(left) != len(right) {
		return false
	}

	for i := range left {
		if left[i] != right[i] {
			return false
		}
	}

	return true
}

func TestRunHandlesZeroDistanceHeuristic(t *testing.T) {
	distances := [][]float64{
		{0, 0, 2},
		{1, 0, 3},
		{4, 5, 0},
	}
	hints := [][]float64{
		{0, 0, 0},
		{0, 0, 0},
		{0, 0, 0},
	}

	aco := NewACO(1.0, 1.0, 0.8, 0.0, 1, 100.0, distances, hints)
	aco.Run()

	if aco.heuristics[1][0] != 1.0 {
		t.Fatalf("expected positive edge heuristic to use 1/distance, got %f", aco.heuristics[1][0])
	}
	if aco.heuristics[0][2] != 0.5 {
		t.Fatalf("expected positive edge heuristic to use 1/distance, got %f", aco.heuristics[0][2])
	}
	if math.IsInf(aco.heuristics[0][1], 0) || math.IsNaN(aco.heuristics[0][1]) {
		t.Fatalf("expected zero edge heuristic to stay finite, got %f", aco.heuristics[0][1])
	}
	if aco.heuristics[0][1] <= aco.heuristics[1][0] {
		t.Fatalf("expected zero edge heuristic to be preferred over the shortest positive edge, got zero=%f positive=%f", aco.heuristics[0][1], aco.heuristics[1][0])
	}
}
