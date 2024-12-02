package edmonds

import (
	"testing"
)

func TestFindMSA(t *testing.T) {
	// Example graph
	V := []int{0, 1, 2, 3}
	E := []Edge{
		{From: 0, To: 1},
		{From: 0, To: 2},
		{From: 1, To: 2},
		{From: 2, To: 3},
		{From: 1, To: 3},
	}
	w := map[Edge]float64{
		{From: 0, To: 1}: 1.0,
		{From: 0, To: 2}: 2.0,
		{From: 1, To: 2}: 1.5,
		{From: 2, To: 3}: 1.0,
		{From: 1, To: 3}: 3.0,
	}

	// Expected result
	expected := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
		{From: 2, To: 3},
	}

	// Call the function
	result := FindMSA(V, E, 0, w)

	// Compare results
	if !compareEdges(result, expected) {
		t.Errorf("Expected %v, got %v", expected, result)
	}
}

func TestFindMSAWithCycle(t *testing.T) {
	V := []int{0, 1, 2, 3}
	E := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
		{From: 2, To: 1}, // Cycle edge
		{From: 2, To: 3},
	}
	w := map[Edge]float64{
		{From: 0, To: 1}: 1.0,
		{From: 1, To: 2}: 2.0,
		{From: 2, To: 1}: 0.5, // Cycle weight
		{From: 2, To: 3}: 1.0,
	}

	expected := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
		{From: 2, To: 3},
	}

	result := FindMSA(V, E, 0, w)
	if !compareEdges(result, expected) {
		t.Errorf("Expected %v, got %v", expected, result)
	}
}

func compareEdges(a, b []Edge) bool {
	if len(a) != len(b) {
		return false
	}

	m := make(map[Edge]int)
	for _, edge := range a {
		m[edge]++
	}
	for _, edge := range b {
		m[edge]--
		if m[edge] < 0 {
			return false
		}
	}
	return true
}
