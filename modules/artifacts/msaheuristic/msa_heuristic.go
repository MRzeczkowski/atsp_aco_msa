package msaheuristic

import (
	"atsp_aco_msa/modules/algorithms/edmonds"
	"atsp_aco_msa/modules/models"
	"encoding/csv"
	"fmt"
	"os"
	"path"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

type Edge = models.Edge
type finder func(root int, vertices []int, edges []Edge, weights map[Edge]float64) []Edge

func Read(rootMsaHeuristicPath string) ([][]float64, error) {
	msaHeuristicPath := path.Join(rootMsaHeuristicPath, "msa_heuristic.csv")
	return readFromCsv(msaHeuristicPath)
}

func ReadMsas(rootMsaHeuristicPath string) ([][][]float64, error) {
	rootedMsas, err := readRootedMsas(rootMsaHeuristicPath)
	if err != nil {
		return nil, err
	}

	msas := make([][][]float64, len(rootedMsas))
	for i, rootedMsa := range rootedMsas {
		msas[i] = rootedMsa.matrix
	}

	return msas, nil
}

type rootedMsa struct {
	root   int
	matrix [][]float64
}

type rootedMsaPath struct {
	root int
	path string
}

func readRootedMsas(rootMsaHeuristicPath string) ([]rootedMsa, error) {
	rootedPaths, err := readRootedMsaPaths(rootMsaHeuristicPath)
	if err != nil {
		return nil, err
	}

	rootedMsas := make([]rootedMsa, len(rootedPaths))
	for i, rootedPath := range rootedPaths {
		msa, err := readFromCsv(rootedPath.path)
		if err != nil {
			return nil, err
		}

		rootedMsas[i] = rootedMsa{root: rootedPath.root, matrix: msa}
	}

	return rootedMsas, nil
}

func readRootedMsaPaths(rootMsaHeuristicPath string) ([]rootedMsaPath, error) {
	msaRootPath := path.Join(rootMsaHeuristicPath, "msas")
	msasPaths, err := filepath.Glob(filepath.Join(msaRootPath, "*.csv"))
	if err != nil {
		return nil, err
	}

	rootedPaths := make([]rootedMsaPath, len(msasPaths))
	for i, msaPath := range msasPaths {
		root, err := strconv.Atoi(strings.TrimSuffix(filepath.Base(msaPath), filepath.Ext(msaPath)))
		if err != nil {
			return nil, fmt.Errorf("parse MSA root from %s: %w", msaPath, err)
		}
		rootedPaths[i] = rootedMsaPath{root: root, path: msaPath}
	}

	sort.SliceStable(rootedPaths, func(i, j int) bool {
		return rootedPaths[i].root < rootedPaths[j].root
	})

	return rootedPaths, nil
}

func Create(matrix [][]float64, rootMsaHeuristicPath string) ([][]float64, error) {
	return createComposite(matrix, rootMsaHeuristicPath, "msa_heuristic.csv", "msas", edmonds.FindMSA)
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
