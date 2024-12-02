package compositeMsa

import (
	"atsp_aco_msa/modules/algorithms/edmonds"
	"atsp_aco_msa/modules/models"
	"encoding/csv"
	"os"
	"strconv"
)

type Edge = models.Edge

func ReadFromCsv(path string) ([][]float64, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)

	rows, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	data := make([][]float64, len(rows))
	for i, row := range rows {
		data[i] = make([]float64, len(row))
		for j, value := range row {
			data[i][j], err = strconv.ParseFloat(value, 64)
			if err != nil {
				return nil, err
			}
		}
	}

	return data, nil
}

func CreateFromData(matrix [][]float64) [][]float64 {
	vertices, edges, weights := models.ConvertToEdges(matrix)

	dimension := len(matrix)

	cmsa := make([][]float64, dimension)
	for i := range dimension {
		cmsa[i] = make([]float64, dimension)
	}

	msas := make([][]Edge, dimension)
	occurrences := make(map[Edge]float64)

	for i := 0; i < dimension; i++ {

		msa := edmonds.FindMSA(vertices, edges, i, weights)

		msas[i] = msa

		for _, edge := range msa {
			occurrences[edge]++
		}
	}

	for _, msa := range msas {
		for _, edge := range msa {
			cmsa[edge.From][edge.To] = occurrences[edge]
		}
	}

	return cmsa
}

func SaveToCsv(cmsa [][]float64, path string) error {

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	for _, row := range cmsa {
		strRow := make([]string, len(row))
		for i, value := range row {
			strRow[i] = strconv.FormatFloat(value, 'f', -1, 64)
		}
		err := writer.Write(strRow)
		if err != nil {
			return err
		}
	}

	return nil
}
