package compositeMsa

import (
	"atsp_aco_msa/modules/algorithms/edmonds"
	"atsp_aco_msa/modules/models"
	"encoding/csv"
	"fmt"
	"os"
	"path"
	"strconv"
)

type Edge = models.Edge

func Read(rootCmsaPath string) ([][]float64, error) {
	cmsaPath := path.Join(rootCmsaPath, "cmsa.csv")
	return readFromCsv(cmsaPath)
}

func Create(matrix [][]float64, rootCmsaPath string) ([][]float64, error) {
	vertices, edges, weights := models.ConvertToEdges(matrix)

	dimension := len(matrix)

	cmsa := make([][]float64, dimension)
	for i := range dimension {
		cmsa[i] = make([]float64, dimension)
	}

	msaRootPath := path.Join(rootCmsaPath, "msas")
	err := os.MkdirAll(msaRootPath, os.ModePerm)
	if err != nil {
		return nil, err
	}

	msas := make([][]Edge, dimension)
	occurrences := make(map[Edge]float64)

	for i := 0; i < dimension; i++ {

		msaCsvFileName := fmt.Sprintf("%d.csv", i)
		msaPath := path.Join(msaRootPath, msaCsvFileName)

		var msa []Edge
		msaMatrix, _ := readFromCsv(msaPath)

		if msaMatrix != nil {
			_, msa, weights = models.ConvertToEdges(msaMatrix)

			for _, edge := range msa {
				if weights[edge] == 0 {
					continue
				}

				occurrences[edge]++
			}
		} else {
			msa = edmonds.FindMSA(i, vertices, edges, weights)
			msaMatrix = models.ConvertToMatrix(msa, dimension)

			err := saveToCsv(msaMatrix, msaPath)
			if err != nil {
				return nil, err
			}

			for _, edge := range msa {
				occurrences[edge]++
			}
		}

		msas[i] = msa
	}

	for _, msa := range msas {
		for _, edge := range msa {
			cmsa[edge.From][edge.To] = occurrences[edge]
		}
	}

	cmsaPath := path.Join(rootCmsaPath, "cmsa.csv")
	err = saveToCsv(cmsa, cmsaPath)
	if err != nil {
		return nil, err
	}

	return cmsa, nil
}

func readFromCsv(path string) ([][]float64, error) {
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

func saveToCsv(matrix [][]float64, path string) error {

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	for _, row := range matrix {
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
