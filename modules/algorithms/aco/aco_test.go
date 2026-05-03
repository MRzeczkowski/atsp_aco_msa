package aco

import "testing"

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
