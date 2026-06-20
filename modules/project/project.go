package project

import (
	"atsp_aco_msa/modules/parsing"
	"fmt"
	"path/filepath"
	"sort"
	"strings"
)

var ArtifactsDirectoryName = "artifacts"
var ResultsDirectoryName = filepath.Join(ArtifactsDirectoryName, "tuning")
var FinalResultsDirectoryName = filepath.Join(ArtifactsDirectoryName, "final", "no_3opt")
var FinalThreeOptResultsDirectoryName = filepath.Join(ArtifactsDirectoryName, "final", "with_3opt")
var MsaHeuristicArtifactsDirectoryName = filepath.Join(ArtifactsDirectoryName, "msa")
var CycleCoverArtifactsDirectoryName = filepath.Join(ArtifactsDirectoryName, "cycle_cover")
var SolutionArtifactsDirectoryName = filepath.Join(ArtifactsDirectoryName, "solutions")
var ResultFileName = "result.csv"

type AtspData struct {
	Name         string
	Matrix       [][]float64
	KnownOptimal float64

	MsaHeuristicDirectoryPath string

	CycleCoverDirectoryPath string

	MsaHeuristicHeatmapPlotPath   string
	MsaHeuristicHistogramPlotPath string

	ResultFilePath       string
	ResultPlotFilePrefix string

	OptimalUniqueToursCsvPath               string
	ToursHeatmapPlotPath                    string
	ToursHistogramPlotPath                  string
	MsaHeuristicToursOverlapHeatmapPlotPath string
}

func MakeAtspData(name string, matrix [][]float64, knownOptimal float64) AtspData {
	return MakeAtspDataInResultsDirectory(name, matrix, knownOptimal, ResultsDirectoryName)
}

func MakeAtspDataInResultsDirectory(name string, matrix [][]float64, knownOptimal float64, resultsRootPath string) AtspData {
	name = strings.TrimSuffix(name, ".atsp")
	resultsDirectoryPath := filepath.Join(resultsRootPath, name)
	resultsPlotsDirectoryPath := filepath.Join(resultsDirectoryPath, "plots")
	msaHeuristicDirectoryPath := filepath.Join(MsaHeuristicArtifactsDirectoryName, name)
	msaHeuristicPlotsDirectoryPath := filepath.Join(msaHeuristicDirectoryPath, "plots")
	cycleCoverDirectoryPath := filepath.Join(CycleCoverArtifactsDirectoryName, name)
	solutionsDirectoryPath := filepath.Join(SolutionArtifactsDirectoryName, name)
	solutionsPlotsDirectoryPath := filepath.Join(solutionsDirectoryPath, "plots")

	return AtspData{
		Name:                                    name,
		Matrix:                                  matrix,
		KnownOptimal:                            knownOptimal,
		MsaHeuristicDirectoryPath:               msaHeuristicDirectoryPath,
		CycleCoverDirectoryPath:                 cycleCoverDirectoryPath,
		MsaHeuristicHeatmapPlotPath:             filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_heatmap.png"),
		MsaHeuristicHistogramPlotPath:           filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_histogram.png"),
		ResultFilePath:                          filepath.Join(resultsDirectoryPath, ResultFileName),
		ResultPlotFilePrefix:                    filepath.Join(resultsPlotsDirectoryPath, "best_result"),
		OptimalUniqueToursCsvPath:               filepath.Join(solutionsDirectoryPath, "solutions.csv"),
		ToursHeatmapPlotPath:                    filepath.Join(solutionsPlotsDirectoryPath, "tours_heatmap.png"),
		ToursHistogramPlotPath:                  filepath.Join(solutionsPlotsDirectoryPath, "tours_histogram.png"),
		MsaHeuristicToursOverlapHeatmapPlotPath: filepath.Join(msaHeuristicPlotsDirectoryPath, "msa_heuristic_tours_overlap_heatmap.png"),
	}
}

func WithExperimentOutputRoot(atspData AtspData, resultsRootPath string) AtspData {
	output := MakeAtspDataInResultsDirectory(atspData.Name, atspData.Matrix, atspData.KnownOptimal, resultsRootPath)
	output.MsaHeuristicDirectoryPath = atspData.MsaHeuristicDirectoryPath
	output.CycleCoverDirectoryPath = atspData.CycleCoverDirectoryPath
	output.OptimalUniqueToursCsvPath = atspData.OptimalUniqueToursCsvPath
	return output
}

const (
	InstanceSetSmoke      = "smoke"
	InstanceSetTuning     = "tuning"
	InstanceSetEvaluation = "evaluation"
	InstanceSetAllKnown   = "all-known"
)

var SmokeInstanceFiles = []string{
	"br17.atsp",
}

var TuningInstanceFiles = []string{
	"ftv33.atsp",
	"p43.atsp",
	"ft53.atsp",
	"ftv64.atsp",
	"crane66_1.atsp",
	"atex5.atsp",
	"ftv90.atsp",
	"ry48p.atsp",
	"crane100_1.atsp",
	"td100_1.atsp",
	"ftv120.atsp",
	"dc134.atsp",
	"ftv150.atsp",
	"dc188.atsp",
	"rbg323.atsp",
}

var EvaluationInstanceFiles = []string{
	"atex1.atsp",
	"atex3.atsp",
	"atex4.atsp",
	"ftv35.atsp",
	"ftv38.atsp",
	"ftv44.atsp",
	"ftv47.atsp",
	"ftv55.atsp",
	"crane66_0.atsp",
	"crane66_2.atsp",
	"ft70.atsp",
	"ftv70.atsp",
	"crane100_0.atsp",
	"crane100_2.atsp",
	"ftv100.atsp",
	"ftv110.atsp",
	"dc112.atsp",
	"dc126.atsp",
	"ftv130.atsp",
	"ftv140.atsp",
	"ftv160.atsp",
	"ftv170.atsp",
	"dc176.atsp",
	"code198.atsp",
	"rbg358.atsp",
	"rbg403.atsp",
	"rbg443.atsp",
}

func SelectAtspFiles(atspFilePaths []string, instanceSet string) ([]string, error) {
	switch instanceSet {
	case InstanceSetSmoke:
		return selectConfiguredAtspFiles(atspFilePaths, SmokeInstanceFiles)
	case InstanceSetTuning:
		return selectConfiguredAtspFiles(atspFilePaths, TuningInstanceFiles)
	case InstanceSetEvaluation:
		return selectConfiguredAtspFiles(atspFilePaths, EvaluationInstanceFiles)
	case InstanceSetAllKnown:
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
		return nil, fmt.Errorf("unsupported -instances value %q; use %q, %q, %q, or %q", instanceSet, InstanceSetSmoke, InstanceSetTuning, InstanceSetEvaluation, InstanceSetAllKnown)
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

func LoadSelectedAtspData(instances string) ([]AtspData, error) {
	tsplibDir := "tsplib_files"
	atspFilesPaths, err := filepath.Glob(filepath.Join(tsplibDir, "*.atsp"))
	if err != nil {
		return nil, err
	}

	atspFilesPaths, err = SelectAtspFiles(atspFilesPaths, instances)
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

		atspsData[i] = MakeAtspData(name, matrix, knownOptimal)
	}

	return atspsData, nil
}
