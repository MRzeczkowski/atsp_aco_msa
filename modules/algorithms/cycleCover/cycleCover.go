package cycleCover

import (
	"atsp_aco_msa/modules/algorithms/hungarian"
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
)

const FileName = "cycle_cover.csv"

func CsvPath(rootPath string) string {
	return filepath.Join(rootPath, FileName)
}

func Read(rootPath string, dimension int) ([][]float64, error) {
	file, err := os.Open(CsvPath(rootPath))
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	matrix := make([][]float64, len(records))
	for i, record := range records {
		matrix[i] = make([]float64, len(record))
		for j, value := range record {
			matrix[i][j], err = strconv.ParseFloat(value, 64)
			if err != nil {
				return nil, err
			}
		}
	}

	if err := Validate(matrix, dimension); err != nil {
		return nil, err
	}

	return matrix, nil
}

func ReadOrCreate(matrix [][]float64, rootPath string) ([][]float64, float64, error) {
	cycleCover, err := Read(rootPath, len(matrix))
	if err == nil {
		return cycleCover, MatrixCost(matrix, cycleCover), nil
	}

	return Create(matrix, rootPath)
}

func Rebuild(matrix [][]float64, rootPath string) ([][]float64, float64, error) {
	if err := os.RemoveAll(rootPath); err != nil {
		return nil, 0, err
	}

	return Create(matrix, rootPath)
}

func Create(matrix [][]float64, rootPath string) ([][]float64, float64, error) {
	cycleCover, cost, err := BuildMatrix(matrix)
	if err != nil {
		return nil, 0, err
	}
	if err := Save(rootPath, cycleCover); err != nil {
		return nil, 0, err
	}

	return cycleCover, cost, nil
}

func BuildMatrix(matrix [][]float64) ([][]float64, float64, error) {
	edges, cost, err := hungarian.MinimumCycleCover(matrix)
	if err != nil {
		return nil, 0, err
	}

	cycleCover := make([][]float64, len(matrix))
	for i := range cycleCover {
		cycleCover[i] = make([]float64, len(matrix))
	}

	for _, edge := range edges {
		cycleCover[edge.From][edge.To] = 1.0
	}

	return cycleCover, cost, nil
}

func Save(rootPath string, matrix [][]float64) error {
	if err := os.MkdirAll(rootPath, os.ModePerm); err != nil {
		return err
	}

	file, err := os.Create(CsvPath(rootPath))
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	for _, row := range matrix {
		record := make([]string, len(row))
		for i, value := range row {
			record[i] = strconv.FormatFloat(value, 'f', -1, 64)
		}
		if err := writer.Write(record); err != nil {
			return err
		}
	}

	return nil
}

func Validate(matrix [][]float64, dimension int) error {
	if len(matrix) != dimension {
		return fmt.Errorf("cycle-cover dimension %d does not match matrix dimension %d", len(matrix), dimension)
	}

	inDegree := make([]int, dimension)
	outDegree := make([]int, dimension)
	for from, row := range matrix {
		if len(row) != dimension {
			return fmt.Errorf("cycle-cover row %d has length %d, want %d", from, len(row), dimension)
		}

		for to, value := range row {
			if value == 0 {
				continue
			}
			if value < 0 {
				return fmt.Errorf("cycle-cover edge %d -> %d has negative value %.2f", from, to, value)
			}
			if from == to {
				return fmt.Errorf("cycle-cover contains self-loop at vertex %d", from)
			}

			outDegree[from]++
			inDegree[to]++
		}
	}

	for vertex := 0; vertex < dimension; vertex++ {
		if inDegree[vertex] != 1 || outDegree[vertex] != 1 {
			return fmt.Errorf("cycle-cover vertex %d has in-degree %d and out-degree %d", vertex, inDegree[vertex], outDegree[vertex])
		}
	}

	return nil
}

func MatrixCost(matrix, cycleCover [][]float64) float64 {
	cost := 0.0
	for i := 0; i < len(cycleCover); i++ {
		for j, value := range cycleCover[i] {
			if value != 0 {
				cost += matrix[i][j]
			}
		}
	}

	return cost
}
