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
	heuristicModifiers := newNeutralHeuristicModifiers(len(distances))

	aco := NewACO(1.0, 2.0, 0.8, 1, 100.0, distances, heuristicModifiers)
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
	for i := range distances {
		distances[i] = make([]float64, 12)

		for j := range distances[i] {
			distances[i][j] = 1.0
		}
	}

	for j := 1; j < 12; j++ {
		distances[0][j] = float64(j)
	}
	neutralModifiers := newNeutralHeuristicModifiers(len(distances))
	strongFarEdgeModifiers := newNeutralHeuristicModifiers(len(distances))
	strongFarEdgeModifiers[0][11] = 2.0

	withoutModifier := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, neutralModifiers)
	withModifier := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, strongFarEdgeModifiers)

	if !equalIntSlices(withoutModifier.neighborsLists[0], withModifier.neighborsLists[0]) {
		t.Fatalf("expected heuristic modifiers not to change candidate list, got without=%v with=%v", withoutModifier.neighborsLists[0], withModifier.neighborsLists[0])
	}

	for _, neighbor := range withModifier.neighborsLists[0] {
		if neighbor == 11 {
			t.Fatalf("expected distance-only candidate list to exclude far modified edge, got %v", withModifier.neighborsLists[0])
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
	heuristicModifiers := newNeutralHeuristicModifiers(len(distances))

	aco := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, heuristicModifiers)
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

func TestRunAppliesHeuristicModifierMatrix(t *testing.T) {
	distances := [][]float64{
		{0, 2, 4},
		{3, 0, 5},
		{6, 7, 0},
	}
	heuristicModifiers := newNeutralHeuristicModifiers(len(distances))
	heuristicModifiers[0][1] = 3.0

	aco := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, heuristicModifiers)
	aco.Run()

	expected := 1.5
	if aco.heuristics[0][1] != expected {
		t.Fatalf("expected modifier to be applied before beta, got %f", aco.heuristics[0][1])
	}
}

func newNeutralHeuristicModifiers(dimension int) [][]float64 {
	modifiers := make([][]float64, dimension)
	for i := range modifiers {
		modifiers[i] = make([]float64, dimension)
		for j := range modifiers[i] {
			modifiers[i][j] = 1.0
		}
	}
	return modifiers
}
