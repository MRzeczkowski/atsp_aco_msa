package app

import (
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	workerpool "atsp_aco_msa/modules/experiments/workers"
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"os"
	"time"
)

func runRebuildCacheMode(atspsData []AtspData, workers int) error {
	return workerpool.RunJobs(instanceJobsByDescendingDimension(atspsData, "rebuild-cache", func(atspData AtspData) error {
		start := time.Now()
		fmt.Printf("[%s][%s] Rebuilding cached MSA artifacts in %s\n", logTimestamp(start), atspData.Name, atspData.MsaHeuristicDirectoryPath)

		if err := os.RemoveAll(atspData.MsaHeuristicDirectoryPath); err != nil {
			return fmt.Errorf("%s: remove MSA artifacts: %w", atspData.Name, err)
		}

		msaHeuristicMatrix, err := msaHeuristic.Create(atspData.Matrix, atspData.MsaHeuristicDirectoryPath)
		if err != nil {
			return fmt.Errorf("%s: create MSA artifacts: %w", atspData.Name, err)
		}

		if err := saveMsaHeuristicPlots(atspData, msaHeuristicMatrix); err != nil {
			return fmt.Errorf("%s: save MSA plots: %w", atspData.Name, err)
		}

		fmt.Printf("[%s][%s] Rebuilt cached MSA artifacts in %s\n", logTimestamp(time.Now()), atspData.Name, time.Since(start).Round(time.Millisecond))
		return nil
	}), workers)
}

func ensureMsaHeuristicArtifacts(atspsData []AtspData, workers int, requireRooted bool) error {
	return workerpool.RunJobs(instanceJobsByDescendingDimension(atspsData, "ensure-msa-artifacts", func(atspData AtspData) error {
		name := atspData.Name
		matrix := atspData.Matrix
		msaHeuristicDirectoryPath := atspData.MsaHeuristicDirectoryPath

		msaHeuristicMatrix, err := msaHeuristic.Read(atspData.MsaHeuristicDirectoryPath)

		if err != nil {
			start := time.Now()
			msaHeuristicMatrix, err = msaHeuristic.Create(matrix, msaHeuristicDirectoryPath)
			elapsed := time.Since(start)

			fmt.Printf("\t[%s] Creating %s took: %d ms\n", name, msaHeuristicDirectoryPath, elapsed.Milliseconds())

			if err != nil {
				return fmt.Errorf("error saving MSA Heuristic: %w", err)
			}
		}

		if requireRooted {
			if _, err := readOrCreateIndividualMsas(atspData); err != nil {
				return err
			}
			msaHeuristicMatrix, err = msaHeuristic.Read(atspData.MsaHeuristicDirectoryPath)
			if err != nil {
				return err
			}
		}

		return saveMsaHeuristicPlots(atspData, msaHeuristicMatrix)
	}), workers)
}

func saveMsaHeuristicPlots(atspData AtspData, msaHeuristicMatrix [][]float64) error {
	msaHeuristicHeatmapPlotTitle := atspData.Name + " MSA heuristic heatmap"
	if err := utilities.SaveHeatmapFromMatrix(msaHeuristicMatrix, msaHeuristicHeatmapPlotTitle, atspData.MsaHeuristicHeatmapPlotPath); err != nil {
		return err
	}

	dataForHistogram := filterZeroes(flattenMatrix(msaHeuristicMatrix))
	msaHeuristicHistogramPlotTitle := atspData.Name + " MSA heuristic histogram"
	dimension := len(atspData.Matrix)
	if err := utilities.SaveHistogramFromData(dataForHistogram, dimension-1, msaHeuristicHistogramPlotTitle, atspData.MsaHeuristicHistogramPlotPath); err != nil {
		return err
	}

	return nil
}

func ensureMsaHeuristicCache(atspsData []AtspData, workers int, requireRooted bool) error {
	return workerpool.RunJobs(instanceJobsByDescendingDimension(atspsData, "ensure-msa-cache", func(atspData AtspData) error {
		if _, err := msaHeuristic.Read(atspData.MsaHeuristicDirectoryPath); err == nil {
			if requireRooted {
				_, err = readOrCreateIndividualMsas(atspData)
			}
			return err
		}

		start := time.Now()
		if _, err := msaHeuristic.Create(atspData.Matrix, atspData.MsaHeuristicDirectoryPath); err != nil {
			return fmt.Errorf("error saving MSA Heuristic: %w", err)
		}

		if requireRooted {
			if _, err := readOrCreateIndividualMsas(atspData); err != nil {
				return err
			}
		}

		fmt.Printf("\t[%s] Creating %s took: %d ms\n", atspData.Name, atspData.MsaHeuristicDirectoryPath, time.Since(start).Milliseconds())
		return nil
	}), workers)
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
