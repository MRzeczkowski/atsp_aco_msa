package msaSupport

import (
	"atsp_aco_msa/modules/algorithms/edmonds"
	"atsp_aco_msa/modules/models"
	"encoding/csv"
	"fmt"
	"os"
	"path"
	"path/filepath"
	"slices"
	"strconv"
)

type Edge = models.Edge
type finder func(root int, vertices []int, edges []Edge, weights map[Edge]float64) []Edge

func Read(rootMsaSupportPath string) ([][]float64, error) {
	msaSupportPath := path.Join(rootMsaSupportPath, "msa_support.csv")
	return readFromCsv(msaSupportPath)
}

func ReadAnti(rootAntiSupportPath string) ([][]float64, error) {
	antiSupportPath := path.Join(rootAntiSupportPath, "msaa_support.csv")
	return readFromCsv(antiSupportPath)
}

func ReadMsas(rootMsaSupportPath string) ([][][]float64, error) {

	msaRootPath := path.Join(rootMsaSupportPath, "msas")
	msasPaths, err := filepath.Glob(filepath.Join(msaRootPath, "*.csv"))
	if err != nil {
		return nil, err
	}

	slices.Sort(msasPaths)

	msas := make([][][]float64, len(msasPaths))

	for i, msaPath := range msasPaths {
		msa, err := readFromCsv(msaPath)
		if err != nil {
			return nil, err
		}

		msas[i] = msa
	}

	return msas, nil
}

func Create(matrix [][]float64, rootMsaSupportPath string) ([][]float64, error) {
	return createComposite(matrix, rootMsaSupportPath, "msa_support.csv", "msas", edmonds.FindMSA)
}

func CreateAnti(matrix [][]float64, rootAntiSupportPath string) ([][]float64, error) {
	return createComposite(matrix, rootAntiSupportPath, "msaa_support.csv", "msaas", edmonds.FindMSAA)
}

func createComposite(matrix [][]float64, rootPath, compositeFileName, arborescencesDirectoryName string, find finder) ([][]float64, error) {
	vertices, edges, weights := models.ConvertToEdges(matrix)

	dimension := len(matrix)

	composite := make([][]float64, dimension)
	for i := range dimension {
		composite[i] = make([]float64, dimension)
	}

	arborescencesRootPath := path.Join(rootPath, arborescencesDirectoryName)
	err := os.MkdirAll(arborescencesRootPath, os.ModePerm)
	if err != nil {
		return nil, err
	}

	arborescences := make([][]Edge, dimension)
	occurrences := make(map[Edge]float64)

	for i := 0; i < dimension; i++ {

		arborescenceCsvFileName := fmt.Sprintf("%d.csv", i)
		arborescencePath := path.Join(arborescencesRootPath, arborescenceCsvFileName)

		var arborescence []Edge
		arborescenceMatrix, _ := readFromCsv(arborescencePath)

		if arborescenceMatrix != nil {
			_, arborescence, arborescenceMatrixWeights := models.ConvertToEdges(arborescenceMatrix)

			for _, edge := range arborescence {
				if arborescenceMatrixWeights[edge] == 0 {
					continue
				}

				occurrences[edge]++
			}
		} else {
			arborescence = find(i, vertices, edges, weights)
			arborescenceMatrix = models.ConvertToMatrix(arborescence, dimension)

			err := saveToCsv(arborescenceMatrix, arborescencePath)
			if err != nil {
				return nil, err
			}

			for _, edge := range arborescence {
				occurrences[edge]++
			}
		}

		arborescences[i] = arborescence
	}

	for _, arborescence := range arborescences {
		for _, edge := range arborescence {
			composite[edge.From][edge.To] = occurrences[edge]
		}
	}

	compositePath := path.Join(rootPath, compositeFileName)
	err = saveToCsv(composite, compositePath)
	if err != nil {
		return nil, err
	}

	return composite, nil
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
