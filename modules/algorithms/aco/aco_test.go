package aco

import (
	"math"
	"math/rand/v2"
	"testing"
)

func TestRunKeepsSelfLoopDesirabilityZero(t *testing.T) {
	distances := [][]float64{
		{0, 1, 2},
		{3, 0, 4},
		{5, 6, 0},
	}
	heuristicModifiers := newNeutralHeuristicModifiers(len(distances))

	aco := NewACO(1.0, 2.0, 0.8, 1, 100.0, distances, heuristicModifiers, 0.0)
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
	strongFarEdgeModifiers[0][11] = 1.0

	withoutModifier := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, neutralModifiers, 0.0)
	withModifier := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, strongFarEdgeModifiers, 0.0)

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

	aco := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, heuristicModifiers, 0.0)
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

func TestRunDoesNotBakeHeuristicModifierIntoVisibility(t *testing.T) {
	distances := [][]float64{
		{0, 2, 4},
		{3, 0, 5},
		{6, 7, 0},
	}
	heuristicModifiers := newNeutralHeuristicModifiers(len(distances))
	heuristicModifiers[0][1] = 1.0

	aco := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, heuristicModifiers, 1.0)
	aco.Run()

	expected := 0.5
	if aco.heuristics[0][1] != expected {
		t.Fatalf("expected static visibility to use only 1/distance, got %f", aco.heuristics[0][1])
	}
}

func TestBuildConstructionHeuristicNeighborsListsUsesBoostedEdges(t *testing.T) {
	modifiers := [][]float64{
		{0.0, 1.0, 0.0},
		{1.0, 0.0, 1.0},
		{0.0, 0.0, 0.0},
	}

	neighborsLists := buildConstructionHeuristicNeighborsLists(modifiers)

	if !equalIntSlices(neighborsLists[0], []int{1}) {
		t.Fatalf("unexpected boosted neighbors for row 0: %v", neighborsLists[0])
	}
	if !equalIntSlices(neighborsLists[1], []int{0, 2}) {
		t.Fatalf("unexpected boosted neighbors for row 1: %v", neighborsLists[1])
	}
	if len(neighborsLists[2]) != 0 {
		t.Fatalf("expected no boosted neighbors for row 2, got %v", neighborsLists[2])
	}
}

func TestConstructionHeuristicProbabilityDecays(t *testing.T) {
	aco := &ACO{
		iterations:      100,
		heuristicWeight: 0.8,
	}

	aco.currentIteration = 0
	if probability := aco.constructionHeuristicProbability(); probability != 0.8 {
		t.Fatalf("expected initial probability 0.8, got %f", probability)
	}

	aco.currentIteration = 50
	if probability := aco.constructionHeuristicProbability(); probability != 0.4 {
		t.Fatalf("expected halfway probability 0.4, got %f", probability)
	}

	aco.currentIteration = 100
	if probability := aco.constructionHeuristicProbability(); probability != 0.0 {
		t.Fatalf("expected final probability 0, got %f", probability)
	}
}

func TestSelectNextCityCanUseConstructionHeuristicOutsideDistanceNeighbors(t *testing.T) {
	aco := &ACO{
		iterations:      10,
		heuristicWeight: 1.0,
		neighborsLists: [][]int{
			{1},
		},
		heuristicNeighborsLists: [][]int{
			{2},
		},
		desirabilities: [][]float64{
			{0.0, 100.0, 1.0, 5.0},
		},
		rng: rand.New(rand.NewPCG(1, 2)),
	}
	canVisitBits := []float64{0.0, 1.0, 1.0, 1.0}
	desirabilities := make([]float64, 4)

	nextCity := aco.selectNextCity(0, canVisitBits, desirabilities)

	if nextCity != 2 {
		t.Fatalf("expected construction heuristic to choose boosted city 2, got %d", nextCity)
	}
}

func TestSelectNextCityFallsBackWhenConstructionHeuristicHasNoFeasibleEdge(t *testing.T) {
	aco := &ACO{
		iterations:      10,
		heuristicWeight: 1.0,
		neighborsLists: [][]int{
			{1},
		},
		heuristicNeighborsLists: [][]int{
			{2},
		},
		desirabilities: [][]float64{
			{0.0, 100.0, 1.0, 5.0},
		},
		rng: rand.New(rand.NewPCG(1, 2)),
	}
	canVisitBits := []float64{0.0, 1.0, 0.0, 1.0}
	desirabilities := make([]float64, 4)

	nextCity := aco.selectNextCity(0, canVisitBits, desirabilities)

	if nextCity != 1 {
		t.Fatalf("expected normal distance-neighbor selection to choose city 1, got %d", nextCity)
	}
}

func TestNewACODisablesThreeOptByDefault(t *testing.T) {
	distances := [][]float64{
		{0, 1, 2},
		{3, 0, 4},
		{5, 6, 0},
	}
	heuristicModifiers := newNeutralHeuristicModifiers(len(distances))

	aco := NewACO(1.0, 1.0, 0.8, 1, 100.0, distances, heuristicModifiers, 0.0)
	if aco.useThreeOpt {
		t.Fatal("expected reduced 3-opt to be disabled by default")
	}

	aco.SetUseThreeOpt(true)
	if !aco.useThreeOpt {
		t.Fatal("expected SetUseThreeOpt(true) to enable reduced 3-opt")
	}
}

func TestUseGlobalBestPheromoneUpdateUsesIterationBudgetPercentages(t *testing.T) {
	tests := []struct {
		name             string
		currentIteration int
		iterations       int
		expected         bool
	}{
		{name: "before first boundary uses iteration best", currentIteration: 249, iterations: 5000, expected: false},
		{name: "first boundary enters every fifth phase", currentIteration: 250, iterations: 5000, expected: true},
		{name: "every fifth phase skips non-multiple", currentIteration: 251, iterations: 5000, expected: false},
		{name: "second boundary enters every third phase", currentIteration: 1500, iterations: 5000, expected: true},
		{name: "every third phase skips non-multiple", currentIteration: 1501, iterations: 5000, expected: false},
		{name: "third boundary enters every second phase", currentIteration: 3500, iterations: 5000, expected: true},
		{name: "every second phase skips odd iteration", currentIteration: 3501, iterations: 5000, expected: false},
		{name: "final boundary uses global best only", currentIteration: 4750, iterations: 5000, expected: true},
		{name: "final phase ignores modulo", currentIteration: 4751, iterations: 5000, expected: true},
		{name: "non-positive budget uses iteration best", currentIteration: 0, iterations: 0, expected: false},
		{name: "same schedule scales to short runs", currentIteration: 95, iterations: 100, expected: true},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			actual := useGlobalBestPheromoneUpdate(test.currentIteration, test.iterations)
			if actual != test.expected {
				t.Fatalf("expected %t, got %t", test.expected, actual)
			}
		})
	}
}

func newNeutralHeuristicModifiers(dimension int) [][]float64 {
	modifiers := make([][]float64, dimension)
	for i := range modifiers {
		modifiers[i] = make([]float64, dimension)
	}
	return modifiers
}
