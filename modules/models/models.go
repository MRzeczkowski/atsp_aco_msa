package models

type Edge struct {
	From, To int
}

func ConvertToEdges(matrix [][]float64) ([]int, []Edge, map[Edge]float64) {
	var vertices []int
	edges := make([]Edge, 0)
	weights := make(map[Edge]float64)

	for i := range matrix {
		vertices = append(vertices, i)

		for j := range matrix[i] {
			edge := Edge{From: i, To: j}
			edges = append(edges, edge)
			weights[edge] = matrix[i][j]
		}
	}

	return vertices, edges, weights
}

func ConvertToMatrix(edges []Edge, size int) [][]float64 {
	matrix := make([][]float64, size)
	for i := range matrix {
		matrix[i] = make([]float64, size)
	}

	for _, edge := range edges {
		matrix[edge.From][edge.To] = 1 // Use 1 to indicate an edge.
	}

	return matrix
}
