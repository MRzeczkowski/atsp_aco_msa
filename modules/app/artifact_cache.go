package app

import (
	"atsp_aco_msa/modules/artifacts/cyclecover"
	"atsp_aco_msa/modules/artifacts/msaheuristic"
	workerpool "atsp_aco_msa/modules/experiments/workers"
	"atsp_aco_msa/modules/utilities"
	"fmt"
	"os"
	"path/filepath"
	"time"
)

func runRebuildCacheMode(atspsData []AtspData, workers int) error {
	return workerpool.RunJobs(instanceJobsByDescendingDimension(atspsData, "rebuild-cache", func(atspData AtspData) error {
		start := time.Now()
		fmt.Printf("[%s][%s] Rebuilding cached MSA artifacts in %s\n", logTimestamp(start), atspData.Name, atspData.MsaHeuristicDirectoryPath)

		if err := removeMsaHeuristicCache(atspData); err != nil {
			return fmt.Errorf("%s: remove MSA cache: %w", atspData.Name, err)
		}

		msaHeuristicMatrix, err := msaheuristic.Create(atspData.Matrix, atspData.MsaHeuristicDirectoryPath)
		if err != nil {
			return fmt.Errorf("%s: create MSA artifacts: %w", atspData.Name, err)
		}

		if err := saveMsaHeuristicPlots(atspData, msaHeuristicMatrix); err != nil {
			return fmt.Errorf("%s: save MSA plots: %w", atspData.Name, err)
		}

		cycleCoverMatrix, cycleCoverCost, err := cyclecover.Rebuild(atspData.Matrix, atspData.CycleCoverDirectoryPath)
		if err != nil {
			return fmt.Errorf("%s: rebuild cycle-cover cache: %w", atspData.Name, err)
		}
		if err := saveCycleCoverPlots(atspData, cycleCoverMatrix); err != nil {
			return fmt.Errorf("%s: save cycle-cover plots: %w", atspData.Name, err)
		}
		fmt.Printf("[%s][%s] Rebuilt cached cycle cover in %s (cost=%.2f)\n", logTimestamp(time.Now()), atspData.Name, cyclecover.Path(atspData.CycleCoverDirectoryPath), cycleCoverCost)

		fmt.Printf("[%s][%s] Rebuilt cached artifacts in %s\n", logTimestamp(time.Now()), atspData.Name, time.Since(start).Round(time.Millisecond))
		return nil
	}), workers)
}

func removeMsaHeuristicCache(atspData AtspData) error {
	if err := removeFileIfExists(filepath.Join(atspData.MsaHeuristicDirectoryPath, "msa_heuristic.csv")); err != nil {
		return err
	}

	return os.RemoveAll(filepath.Join(atspData.MsaHeuristicDirectoryPath, "msas"))
}

func ensureCycleCoverCache(atspsData []AtspData, workers int) error {
	return workerpool.RunJobs(instanceJobsByDescendingDimension(atspsData, "ensure-cycle-cover-cache", func(atspData AtspData) error {
		cycleCoverMatrix, _, err := cyclecover.ReadOrCreate(atspData.Matrix, atspData.CycleCoverDirectoryPath)
		if err != nil {
			return err
		}

		return saveCycleCoverPlots(atspData, cycleCoverMatrix)
	}), workers)
}

func ensureMsaHeuristicArtifacts(atspsData []AtspData, workers int, requireRooted bool) error {
	return workerpool.RunJobs(instanceJobsByDescendingDimension(atspsData, "ensure-msa-artifacts", func(atspData AtspData) error {
		name := atspData.Name
		matrix := atspData.Matrix
		msaHeuristicDirectoryPath := atspData.MsaHeuristicDirectoryPath

		msaHeuristicMatrix, err := msaheuristic.Read(atspData.MsaHeuristicDirectoryPath)

		if err != nil {
			start := time.Now()
			msaHeuristicMatrix, err = msaheuristic.Create(matrix, msaHeuristicDirectoryPath)
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
			msaHeuristicMatrix, err = msaheuristic.Read(atspData.MsaHeuristicDirectoryPath)
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

func saveCycleCoverPlots(atspData AtspData, cycleCoverMatrix [][]float64) error {
	cycleCoverHeatmapPlotTitle := atspData.Name + " cycle cover heatmap"
	return utilities.SaveHeatmapFromMatrix(cycleCoverMatrix, cycleCoverHeatmapPlotTitle, atspData.CycleCoverHeatmapPlotPath)
}

func ensureMsaHeuristicCache(atspsData []AtspData, workers int, requireRooted bool) error {
	return workerpool.RunJobs(instanceJobsByDescendingDimension(atspsData, "ensure-msa-cache", func(atspData AtspData) error {
		if _, err := msaheuristic.Read(atspData.MsaHeuristicDirectoryPath); err == nil {
			if requireRooted {
				_, err = readOrCreateIndividualMsas(atspData)
			}
			return err
		}

		start := time.Now()
		if _, err := msaheuristic.Create(atspData.Matrix, atspData.MsaHeuristicDirectoryPath); err != nil {
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
