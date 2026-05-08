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
	result := FindMSA(0, V, E, w)

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

	result := FindMSA(0, V, E, w)
	if !compareEdges(result, expected) {
		t.Errorf("Expected %v, got %v", expected, result)
	}
}

func TestFindMSAIgnoresSelfLoops(t *testing.T) {
	V := []int{0, 1, 2}
	E := []Edge{
		{From: 0, To: 1},
		{From: 0, To: 2},
		{From: 1, To: 2},
		{From: 1, To: 1},
		{From: 2, To: 2},
	}
	w := map[Edge]float64{
		{From: 0, To: 1}: 1.0,
		{From: 0, To: 2}: 5.0,
		{From: 1, To: 2}: 1.0,
		{From: 1, To: 1}: 0.0,
		{From: 2, To: 2}: 0.0,
	}

	expected := []Edge{
		{From: 0, To: 1},
		{From: 1, To: 2},
	}

	result := FindMSA(0, V, E, w)
	if !compareEdges(result, expected) {
		t.Errorf("Expected %v, got %v", expected, result)
	}
}

func TestFindMSAA(t *testing.T) {
	V := []int{0, 1, 2, 3}
	E := []Edge{
		{From: 1, To: 0},
		{From: 2, To: 0},
		{From: 2, To: 1},
		{From: 3, To: 1},
		{From: 3, To: 2},
	}
	w := map[Edge]float64{
		{From: 1, To: 0}: 1.0,
		{From: 2, To: 0}: 5.0,
		{From: 2, To: 1}: 1.0,
		{From: 3, To: 1}: 3.0,
		{From: 3, To: 2}: 1.0,
	}

	expected := []Edge{
		{From: 1, To: 0},
		{From: 2, To: 1},
		{From: 3, To: 2},
	}

	result := FindMSAA(0, V, E, w)
	if !compareEdges(result, expected) {
		t.Errorf("Expected %v, got %v", expected, result)
	}
}

func TestFindMSAAWithCycle(t *testing.T) {
	V := []int{0, 1, 2, 3}
	E := []Edge{
		{From: 1, To: 0},
		{From: 2, To: 1},
		{From: 1, To: 2},
		{From: 3, To: 2},
	}
	w := map[Edge]float64{
		{From: 1, To: 0}: 1.0,
		{From: 2, To: 1}: 2.0,
		{From: 1, To: 2}: 0.5,
		{From: 3, To: 2}: 1.0,
	}

	expected := []Edge{
		{From: 1, To: 0},
		{From: 2, To: 1},
		{From: 3, To: 2},
	}

	result := FindMSAA(0, V, E, w)
	if !compareEdges(result, expected) {
		t.Errorf("Expected %v, got %v", expected, result)
	}
}

func TestFindMSAARejectsEdgesLeavingRoot(t *testing.T) {
	V := []int{0, 1, 2}
	E := []Edge{
		{From: 0, To: 1},
		{From: 0, To: 2},
		{From: 1, To: 0},
		{From: 1, To: 2},
		{From: 2, To: 0},
		{From: 2, To: 1},
	}
	w := map[Edge]float64{
		{From: 0, To: 1}: 0.0,
		{From: 0, To: 2}: 0.0,
		{From: 1, To: 0}: 3.0,
		{From: 1, To: 2}: 1.0,
		{From: 2, To: 0}: 1.0,
		{From: 2, To: 1}: 3.0,
	}

	expected := []Edge{
		{From: 1, To: 2},
		{From: 2, To: 0},
	}

	result := FindMSAA(0, V, E, w)
	if !compareEdges(result, expected) {
		t.Errorf("Expected %v, got %v", expected, result)
	}

	for _, edge := range result {
		if edge.From == 0 {
			t.Fatalf("antiarborescence rooted at 0 must not contain an edge leaving the root: %v", result)
		}
	}
}

func TestMakeSuperNodeUsesUnusedNegativeVertex(t *testing.T) {
	vertices := []int{0, 1, 2, -1, -2}

	superNode := makeSuperNode(vertices)

	if superNode != -3 {
		t.Fatalf("expected -3, got %d", superNode)
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
