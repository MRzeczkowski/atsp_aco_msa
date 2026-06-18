package app

import (
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"os"
	"time"
)

func runRebuildMsaMode(atspsData []AtspData, workers int) error {
	return runBoundedInstanceJobs(atspsData, workers, func(atspData AtspData) error {
		start := time.Now()
		fmt.Printf("[%s][%s] Rebuilding MSA artifacts in %s\n", logTimestamp(start), atspData.name, atspData.msaHeuristicDirectoryPath)

		if err := os.RemoveAll(atspData.msaHeuristicDirectoryPath); err != nil {
			return fmt.Errorf("%s: remove MSA artifacts: %w", atspData.name, err)
		}

		msaHeuristicMatrix, err := msaHeuristic.Create(atspData.matrix, atspData.msaHeuristicDirectoryPath)
		if err != nil {
			return fmt.Errorf("%s: create MSA artifacts: %w", atspData.name, err)
		}

		if err := saveMsaHeuristicPlots(atspData, msaHeuristicMatrix); err != nil {
			return fmt.Errorf("%s: save MSA plots: %w", atspData.name, err)
		}

		fmt.Printf("[%s][%s] Rebuilt MSA artifacts in %s\n", logTimestamp(time.Now()), atspData.name, time.Since(start).Round(time.Millisecond))
		return nil
	})
}

func ensureMsaHeuristicArtifacts(atspsData []AtspData, workers int, requireRooted bool) error {
	return runBoundedInstanceJobs(atspsData, workers, func(atspData AtspData) error {
		name := atspData.name
		matrix := atspData.matrix
		msaHeuristicDirectoryPath := atspData.msaHeuristicDirectoryPath

		msaHeuristicMatrix, err := msaHeuristic.Read(atspData.msaHeuristicDirectoryPath)

		if err != nil {
			start := time.Now()
			msaHeuristicMatrix, err = msaHeuristic.Create(matrix, msaHeuristicDirectoryPath)
			elapsed := time.Since(start)

			fmt.Printf("\t[%s] Creating %s took: %d ms\n", name, msaHeuristicDirectoryPath, elapsed.Milliseconds())

			if err != nil {
				return fmt.Errorf("error saving MSA heuristic: %w", err)
			}
		}

		if requireRooted {
			if _, err := readOrCreateIndividualMsas(atspData); err != nil {
				return err
			}
			msaHeuristicMatrix, err = msaHeuristic.Read(atspData.msaHeuristicDirectoryPath)
			if err != nil {
				return err
			}
		}

		return saveMsaHeuristicPlots(atspData, msaHeuristicMatrix)
	})
}

func saveMsaHeuristicPlots(atspData AtspData, msaHeuristicMatrix [][]float64) error {
	msaHeuristicHeatmapPlotTitle := atspData.name + " MSA heuristic heatmap"
	if err := utilities.SaveHeatmapFromMatrix(msaHeuristicMatrix, msaHeuristicHeatmapPlotTitle, atspData.msaHeuristicHeatmapPlotPath); err != nil {
		return err
	}

	dataForHistogram := filterZeroes(flattenMatrix(msaHeuristicMatrix))
	msaHeuristicHistogramPlotTitle := atspData.name + " MSA heuristic histogram"
	dimension := len(atspData.matrix)
	if err := utilities.SaveHistogramFromData(dataForHistogram, dimension-1, msaHeuristicHistogramPlotTitle, atspData.msaHeuristicHistogramPlotPath); err != nil {
		return err
	}

	return nil
}

func ensureMsaHeuristicCache(atspsData []AtspData, workers int, requireRooted bool) error {
	return runBoundedInstanceJobs(atspsData, workers, func(atspData AtspData) error {
		if _, err := msaHeuristic.Read(atspData.msaHeuristicDirectoryPath); err == nil {
			if requireRooted {
				_, err = readOrCreateIndividualMsas(atspData)
			}
			return err
		}

		start := time.Now()
		if _, err := msaHeuristic.Create(atspData.matrix, atspData.msaHeuristicDirectoryPath); err != nil {
			return fmt.Errorf("error saving MSA heuristic: %w", err)
		}

		if requireRooted {
			if _, err := readOrCreateIndividualMsas(atspData); err != nil {
				return err
			}
		}

		fmt.Printf("\t[%s] Creating %s took: %d ms\n", atspData.name, atspData.msaHeuristicDirectoryPath, time.Since(start).Milliseconds())
		return nil
	})
}

func filterZeroes(data []float64) []float64 {
	var result []float64
	for _, value := range data {
		if value != 0 {
			result = append(result, value)
		}
	}
	return result
}

func countUniqueFloatValues(data []float64) int {
	uniqueValues := make(map[float64]struct{}, len(data))
	for _, value := range data {
		uniqueValues[value] = struct{}{}
	}
	return len(uniqueValues)
}

func flattenMatrix(matrix [][]float64) []float64 {
	var result []float64
	for _, row := range matrix {
		result = append(result, row...)
	}
	return result
}
