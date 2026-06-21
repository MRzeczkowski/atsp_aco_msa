package app

import (
	"atsp_aco_msa/modules/artifacts/msaheuristic"
	"fmt"
)

func readMsaHeuristicMatrix(atspData AtspData) ([][]float64, error) {
	return msaheuristic.Read(atspData.MsaHeuristicDirectoryPath)
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
