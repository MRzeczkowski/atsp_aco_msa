package models

import "testing"

func TestConvertToEdgesExcludesSelfLoops(t *testing.T) {
	matrix := [][]float64{
		{0, 2, 3},
		{4, 0, 5},
		{6, 7, 0},
	}

	vertices, edges, weights := ConvertToEdges(matrix)

	if len(vertices) != len(matrix) {
		t.Fatalf("expected %d vertices, got %d", len(matrix), len(vertices))
	}

	expectedEdgesCount := len(matrix) * (len(matrix) - 1)
	if len(edges) != expectedEdgesCount {
		t.Fatalf("expected %d edges, got %d", expectedEdgesCount, len(edges))
	}

	for i := range matrix {
		edge := Edge{From: i, To: i}
		if _, exists := weights[edge]; exists {
			t.Fatalf("expected self-loop %v to be excluded", edge)
		}
	}

	for _, edge := range edges {
		if edge.From == edge.To {
			t.Fatalf("expected no self-loops, got %v", edge)
		}
		if weights[edge] != matrix[edge.From][edge.To] {
			t.Fatalf("expected weight %.2f for %v, got %.2f", matrix[edge.From][edge.To], edge, weights[edge])
		}
	}
}
