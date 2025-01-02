package parsing

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

var optimalSolutions = map[string]float64{
	"br17":   39,
	"ft53":   6905,
	"ft70":   38673,
	"ftv33":  1286,
	"ftv35":  1473,
	"ftv38":  1530,
	"ftv44":  1613,
	"ftv47":  1776,
	"ftv55":  1608,
	"ftv64":  1839,
	"ftv70":  1950,
	"ftv170": 2755,
	"p43":    5620,
	"rbg323": 1326,
	"rbg358": 1163,
	"rbg403": 2465,
	"rbg443": 2720,
	"ry48p":  14422,
}

func ParseTSPLIBFile(path string) (name string, matrix [][]float64, knownOptimal float64, err error) {
	file, err := os.Open(path)
	if err != nil {
		return "", nil, 0, err
	}

	defer file.Close()

	scanner := bufio.NewScanner(file)
	readMatrix := false

	// Initialize variables to keep track of matrix data as it's read
	var dimension int
	var valuesInMatrix []float64

	for scanner.Scan() {
		line := scanner.Text()
		if strings.Contains(line, "EOF") {
			break
		}

		if strings.HasPrefix(line, "NAME") {
			parts := strings.Split(line, ":")
			name = strings.TrimSpace(parts[1])
		}

		if strings.HasPrefix(line, "DIMENSION") {
			parts := strings.Split(line, ":")
			dimension, err = strconv.Atoi(strings.TrimSpace(parts[1]))
			if err != nil {
				return "", nil, 0, err
			}

			matrix = make([][]float64, dimension)
			for i := range matrix {
				matrix[i] = make([]float64, dimension)
			}
		}

		if readMatrix {
			// Add to a continuous list of values
			rowValues := strings.Fields(line)
			for _, val := range rowValues {
				num, err := strconv.ParseFloat(val, 64)
				if err != nil {
					return "", nil, 0, err
				}

				valuesInMatrix = append(valuesInMatrix, num)
			}
		}

		if strings.HasPrefix(line, "EDGE_WEIGHT_SECTION") {
			readMatrix = true
		}
	}

	if err := scanner.Err(); err != nil {
		return "", nil, 0, err
	}

	// Populate the matrix from the list of values
	if len(valuesInMatrix) != dimension*dimension {
		return "", nil, 0, fmt.Errorf("the total numbers in matrix (%d) does not match expected dimension squared (%d)", len(valuesInMatrix), dimension*dimension)
	}

	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			matrix[i][j] = valuesInMatrix[i*dimension+j]
		}
	}

	return name, matrix, optimalSolutions[name], nil
}
