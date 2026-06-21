package app

import (
	"atsp_aco_msa/modules/artifacts/msaheuristic"
	"fmt"
	"os"
)

func removeFileIfExists(path string) error {
	if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
		return err
	}

	return nil
}

func readMsaHeuristicMatrixForHeuristic(atspData AtspData, heuristic string) ([][]float64, error) {
	return msaheuristic.Read(atspData.MsaHeuristicDirectoryPath)
}

func readMsaHeuristicMatrixForResultRoot(atspData AtspData, heuristic, resultsRootPath string) ([][]float64, error) {
	return readMsaHeuristicMatrixForHeuristic(atspData, heuristic)
}

func readRootedMsaHeuristics(atspData AtspData) ([][][]float64, error) {
	rootedMsas, err := msaheuristic.ReadMsas(atspData.MsaHeuristicDirectoryPath)
	if err != nil {
		return nil, err
	}
	if len(rootedMsas) != len(atspData.Matrix) {
		return nil, fmt.Errorf("expected %d rooted MSAs, got %d", len(atspData.Matrix), len(rootedMsas))
	}

	return rootedMsas, nil
}
