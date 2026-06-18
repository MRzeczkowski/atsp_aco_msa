package app

import (
	"atsp_aco_msa/modules/parsing"
	"fmt"
	"path/filepath"
	"sort"
	"strings"
)

var artifactsDirectoryName = "artifacts"
var resultsDirectoryName = filepath.Join(artifactsDirectoryName, "tuning")
var finalResultsDirectoryName = filepath.Join(artifactsDirectoryName, "final", "no_3opt")
var finalThreeOptResultsDirectoryName = filepath.Join(artifactsDirectoryName, "final", "with_3opt")
var msaHeuristicArtifactsDirectoryName = filepath.Join(artifactsDirectoryName, "msa")
var solutionArtifactsDirectoryName = filepath.Join(artifactsDirectoryName, "solutions")
var resultFileName = "result.csv"

type AtspData struct {
	name         string
	matrix       [][]float64
	knownOptimal float64

	msaHeuristicDirectoryPath string

	msaHeuristicHeatmapPlotPath   string
	msaHeuristicHistogramPlotPath string

	resultFilePath       string
	resultPlotFilePrefix string

	optimalUniqueToursCsvPath               string
	toursHeatmapPlotPath                    string
	toursHistogramPlotPath                  string
	msaHeuristicToursOverlapHeatmapPlotPath string
}

func makeAtspData(name string, matrix [][]float64, knownOptimal float64) AtspData {
	return makeAtspDataInResultsDirectory(name, matrix, knownOptimal, resultsDirectoryName)
}

func makeAtspDataInResultsDirectory(name string, matrix [][]float64, knownOptimal float64, resultsRootPath string) AtspData {
	name = strings.TrimSuffix(name, ".atsp")
	resultsDirectoryPath := filepath.Join(resultsRootPath, name)
	resultsPlotsDirectoryPath := filepath.Join(resultsDirectoryPath, "plots")
	msaHeuristicDirectoryPath := filepath.Join(msaHeuristicArtifactsDirectoryName, name)
	msaHeuristicPlotsDirectoryPath := filepath.Join(msaHeuristicDirectoryPath, "plots")
	solutionsDirectoryPath := filepath.Join(solutionArtifactsDirectoryName, name)
	solutionsPlotsDirectoryPath := filepath.Join(solutionsDirectoryPath, "plots")

	msaHeuristicHeatmapPlotPath := filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_heatmap.png")
	msaHeuristicHistogramPlotPath := filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_histogram.png")

	resultFilePath := filepath.Join(resultsDirectoryPath, resultFileName)
	resultPlotFilePrefix := filepath.Join(resultsPlotsDirectoryPath, "best_result")

	optimalUniqueToursCsvPath := filepath.Join(solutionsDirectoryPath, "solutions.csv")
	toursHeatmapPlotPath := filepath.Join(solutionsPlotsDirectoryPath, "tours_heatmap.png")
	toursHistogramPlotPath := filepath.Join(solutionsPlotsDirectoryPath, "tours_histogram.png")
	msaHeuristicToursOverlapHeatmapPlotPath := filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_tours_overlap_heatmap.png")

	return AtspData{
		name,
		matrix,
		knownOptimal,

		msaHeuristicDirectoryPath,

		msaHeuristicHeatmapPlotPath, msaHeuristicHistogramPlotPath,

		resultFilePath,
		resultPlotFilePrefix,

		optimalUniqueToursCsvPath,
		toursHeatmapPlotPath,
		toursHistogramPlotPath,
		msaHeuristicToursOverlapHeatmapPlotPath,
	}
}

func withExperimentOutputRoot(atspData AtspData, resultsRootPath string) AtspData {
	output := makeAtspDataInResultsDirectory(atspData.name, atspData.matrix, atspData.knownOptimal, resultsRootPath)
	output.msaHeuristicDirectoryPath = atspData.msaHeuristicDirectoryPath
	output.optimalUniqueToursCsvPath = atspData.optimalUniqueToursCsvPath
	return output
}

func selectAtspFiles(atspFilePaths []string, instanceSet string) ([]string, error) {
	switch instanceSet {
	case instanceSetSmoke:
		return selectConfiguredAtspFiles(atspFilePaths, smokeInstanceFiles)
	case instanceSetTuning:
		return selectConfiguredAtspFiles(atspFilePaths, tuningInstanceFiles)
	case instanceSetEvaluation:
		return selectConfiguredAtspFiles(atspFilePaths, evaluationInstanceFiles)
	case instanceSetAllKnown:
		selected := make([]string, 0, len(atspFilePaths))
		for _, atspFilePath := range atspFilePaths {
			if parsing.HasKnownOptimalSolution(filepath.Base(atspFilePath)) {
				selected = append(selected, atspFilePath)
			}
		}

		sort.Strings(selected)
		if len(selected) == 0 {
			return nil, fmt.Errorf("no ATSP files with known optima were found")
		}

		return selected, nil
	default:
		return nil, fmt.Errorf("unsupported -instances value %q; use %q, %q, %q, or %q", instanceSet, instanceSetSmoke, instanceSetTuning, instanceSetEvaluation, instanceSetAllKnown)
	}
}

func selectConfiguredAtspFiles(atspFilePaths, configuredFiles []string) ([]string, error) {
	pathByFileName := make(map[string]string, len(atspFilePaths))
	for _, atspFilePath := range atspFilePaths {
		pathByFileName[filepath.Base(atspFilePath)] = atspFilePath
	}

	selected := make([]string, 0, len(configuredFiles))
	for _, fileName := range configuredFiles {
		atspFilePath, ok := pathByFileName[fileName]
		if !ok {
			return nil, fmt.Errorf("configured ATSP instance %q was not found", fileName)
		}

		selected = append(selected, atspFilePath)
	}

	return selected, nil
}

func loadSelectedAtspData(instances string) ([]AtspData, error) {
	tsplibDir := "tsplib_files"
	atspFilesPaths, err := filepath.Glob(filepath.Join(tsplibDir, "*.atsp"))
	if err != nil {
		return nil, err
	}

	atspFilesPaths, err = selectAtspFiles(atspFilesPaths, instances)
	if err != nil {
		return nil, err
	}

	atspsData := make([]AtspData, len(atspFilesPaths))
	for i, atspFilePath := range atspFilesPaths {
		name, matrix, knownOptimal, err := parsing.ParseTSPLIBFile(atspFilePath)
		if err != nil {
			return nil, err
		}
		if knownOptimal == 0 {
			return nil, fmt.Errorf("missing known optimal solution for %s", atspFilePath)
		}

		atspsData[i] = makeAtspData(name, matrix, knownOptimal)
	}

	return atspsData, nil
}
