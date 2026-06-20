package structuralComparison

import (
	"atsp_aco_msa/modules/algorithms/cycleCover"
	"atsp_aco_msa/modules/algorithms/heuristics"
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	"atsp_aco_msa/modules/analysis/solutionTours"
	"atsp_aco_msa/modules/models"
	"fmt"
	"sort"
)

type Edge = models.Edge

type InstanceConfig struct {
	Name                      string
	Dimension                 int
	Matrix                    [][]float64
	CycleCoverDirectoryPath   string
	MsaHeuristicDirectoryPath string
	OptimalToursCsvPath       string
}

type Config struct {
	Instances     []InstanceConfig
	HighThreshold float64
	MsaPatchBias  float64
}

type InstanceAnalysis struct {
	Instance  string
	Dimension int
	Metrics   InstanceMetrics
}

type InstanceMetrics struct {
	FoundOptimalTourCount       int
	UniqueFoundOptimalEdgeCount int

	CycleCoverMetrics               EdgeSetMetrics
	HighMsaHeuristicMetrics         EdgeSetMetrics
	CycleCoverMsaPatchingMetrics    EdgeSetMetrics
	CycleCoverHighMsaHeuristicEdges int

	OptimalEdgesInCycleCoverAndHighMsaHeuristic int
	OptimalEdgesInCycleCoverNotHighMsaHeuristic int
	OptimalEdgesInHighMsaHeuristicNotCycleCover int
	OptimalEdgesInNeitherCycleCoverNorHigh      int
}

type EdgeSetMetrics struct {
	EdgeCount        int
	OptimalEdgeCount int
	Precision        float64
	Recall           float64
}

func AnalyzeInstances(config Config) ([]InstanceAnalysis, error) {
	highThreshold := config.HighThreshold
	if highThreshold == 0 {
		highThreshold = 1.0
	}
	msaPatchBias := config.MsaPatchBias

	instances := append([]InstanceConfig(nil), config.Instances...)
	sort.SliceStable(instances, func(i, j int) bool {
		return instances[i].Name < instances[j].Name
	})

	analyses := make([]InstanceAnalysis, 0, len(instances))
	for _, instance := range instances {
		analysis, err := analyzeInstance(instance, highThreshold, msaPatchBias)
		if err != nil {
			return nil, err
		}
		analyses = append(analyses, analysis)
	}

	return analyses, nil
}

func AnalyzeInstance(config InstanceConfig, highThreshold float64) (InstanceAnalysis, error) {
	return analyzeInstance(config, highThreshold, 1.0)
}

func analyzeInstance(config InstanceConfig, highThreshold, msaPatchBias float64) (InstanceAnalysis, error) {
	if highThreshold == 0 {
		highThreshold = 1.0
	}

	if err := validateSquareMatrix(config.Matrix); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid distance matrix: %w", config.Name, err)
	}
	if config.Dimension != 0 && config.Dimension != len(config.Matrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: configured dimension %d does not match matrix dimension %d", config.Name, config.Dimension, len(config.Matrix))
	}

	msaHeuristicMatrix, err := msaHeuristic.Read(config.MsaHeuristicDirectoryPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read MSA heuristic: %w", config.Name, err)
	}
	if err := validateSquareMatrix(msaHeuristicMatrix); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid MSA heuristic: %w", config.Name, err)
	}
	if len(msaHeuristicMatrix) != len(config.Matrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: MSA heuristic dimension %d does not match matrix dimension %d", config.Name, len(msaHeuristicMatrix), len(config.Matrix))
	}

	uniqueOptimalTours, err := solutionTours.ReadOptimalTours(config.OptimalToursCsvPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read found optimal tours: %w", config.Name, err)
	}
	if err := validateTours(uniqueOptimalTours, len(config.Matrix)); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid found optimal tours: %w", config.Name, err)
	}

	cycleCoverMatrix, err := cycleCover.Read(config.CycleCoverDirectoryPath, len(config.Matrix))
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read cycle cover: %w", config.Name, err)
	}
	if err := validateSquareMatrix(cycleCoverMatrix); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid cycle cover: %w", config.Name, err)
	}
	if len(cycleCoverMatrix) != len(config.Matrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: cycle-cover dimension %d does not match matrix dimension %d", config.Name, len(cycleCoverMatrix), len(config.Matrix))
	}
	cycleCoverEdges := buildMatrixEdges(cycleCoverMatrix)
	if err := validateCycleCover(cycleCoverEdges, len(config.Matrix)); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid cycle cover: %w", config.Name, err)
	}

	analysis := calculateAnalysis(config.Name, len(config.Matrix), config.Matrix, msaHeuristicMatrix, uniqueOptimalTours, cycleCoverMatrix, cycleCoverEdges, highThreshold, msaPatchBias)
	return analysis, nil
}

func calculateAnalysis(instance string, dimension int, matrix, msaHeuristic [][]float64, uniqueOptimalTours map[string][]int, cycleCoverMatrix [][]float64, cycleCoverEdges []Edge, highThreshold, msaPatchBias float64) InstanceAnalysis {
	optimalEdges := buildTourEdgeSet(uniqueOptimalTours)
	cycleCoverSet := buildEdgeSet(cycleCoverEdges)
	highMsaHeuristicSet := buildMsaHeuristicThresholdSet(msaHeuristic, highThreshold, true)
	cycleCoverHighMsaHeuristicSet := intersectEdgeSets(cycleCoverSet, highMsaHeuristicSet)
	cycleCoverMsaPatchingSet := buildMatrixEdgeSet(heuristics.BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(matrix, msaHeuristic, cycleCoverMatrix, msaPatchBias))

	metrics := InstanceMetrics{
		FoundOptimalTourCount:           len(uniqueOptimalTours),
		UniqueFoundOptimalEdgeCount:     len(optimalEdges),
		CycleCoverMetrics:               calculateEdgeSetMetrics(cycleCoverSet, optimalEdges),
		HighMsaHeuristicMetrics:         calculateEdgeSetMetrics(highMsaHeuristicSet, optimalEdges),
		CycleCoverMsaPatchingMetrics:    calculateEdgeSetMetrics(cycleCoverMsaPatchingSet, optimalEdges),
		CycleCoverHighMsaHeuristicEdges: len(cycleCoverHighMsaHeuristicSet),
	}

	metrics.OptimalEdgesInCycleCoverAndHighMsaHeuristic = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaHeuristicSet, true, true)
	metrics.OptimalEdgesInCycleCoverNotHighMsaHeuristic = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaHeuristicSet, true, false)
	metrics.OptimalEdgesInHighMsaHeuristicNotCycleCover = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaHeuristicSet, false, true)
	metrics.OptimalEdgesInNeitherCycleCoverNorHigh = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaHeuristicSet, false, false)

	return InstanceAnalysis{
		Instance:  instance,
		Dimension: dimension,
		Metrics:   metrics,
	}
}

func calculateEdgeSetMetrics(edgeSet, optimalEdges map[Edge]struct{}) EdgeSetMetrics {
	optimalEdgeCount := countIntersection(edgeSet, optimalEdges)
	metrics := EdgeSetMetrics{
		EdgeCount:        len(edgeSet),
		OptimalEdgeCount: optimalEdgeCount,
	}

	if len(edgeSet) > 0 && len(optimalEdges) > 0 {
		metrics.Precision = float64(optimalEdgeCount) / float64(len(edgeSet))
	}
	if len(optimalEdges) > 0 {
		metrics.Recall = float64(optimalEdgeCount) / float64(len(optimalEdges))
	}

	return metrics
}

func buildTourEdgeSet(uniqueOptimalTours map[string][]int) map[Edge]struct{} {
	edges := make(map[Edge]struct{})
	tourIds := sortedTourIds(uniqueOptimalTours)
	for _, tourId := range tourIds {
		for _, edge := range models.ConvertTourToEdges(uniqueOptimalTours[tourId]) {
			edges[edge] = struct{}{}
		}
	}
	return edges
}

func buildEdgeSet(edges []Edge) map[Edge]struct{} {
	set := make(map[Edge]struct{}, len(edges))
	for _, edge := range edges {
		set[edge] = struct{}{}
	}
	return set
}

func buildMatrixEdgeSet(matrix [][]float64) map[Edge]struct{} {
	set := make(map[Edge]struct{})
	for from := 0; from < len(matrix); from++ {
		for to, value := range matrix[from] {
			if from != to && value != 0 {
				set[Edge{From: from, To: to}] = struct{}{}
			}
		}
	}
	return set
}

func buildMatrixEdges(matrix [][]float64) []Edge {
	edges := make([]Edge, 0)
	for from := 0; from < len(matrix); from++ {
		for to, value := range matrix[from] {
			if from != to && value != 0 {
				edges = append(edges, Edge{From: from, To: to})
			}
		}
	}
	return edges
}

func buildEdgeMatrix(dimension int, edges []Edge) [][]float64 {
	matrix := make([][]float64, dimension)
	for i := range matrix {
		matrix[i] = make([]float64, dimension)
	}
	for _, edge := range edges {
		if edge.From >= 0 && edge.From < dimension && edge.To >= 0 && edge.To < dimension {
			matrix[edge.From][edge.To] = 1.0
		}
	}
	return matrix
}

func buildMsaHeuristicThresholdSet(msaHeuristic [][]float64, threshold float64, includeEqual bool) map[Edge]struct{} {
	dimension := len(msaHeuristic)
	set := make(map[Edge]struct{})
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			value := normalizedMsaHeuristicValue(msaHeuristic[i][j], dimension)
			if includeEqual {
				if value >= threshold {
					set[Edge{From: i, To: j}] = struct{}{}
				}
			} else if value > threshold {
				set[Edge{From: i, To: j}] = struct{}{}
			}
		}
	}
	return set
}

func intersectEdgeSets(left, right map[Edge]struct{}) map[Edge]struct{} {
	if len(left) > len(right) {
		left, right = right, left
	}

	intersection := make(map[Edge]struct{})
	for edge := range left {
		if _, ok := right[edge]; ok {
			intersection[edge] = struct{}{}
		}
	}
	return intersection
}

func countIntersection(left, right map[Edge]struct{}) int {
	count := 0
	for edge := range left {
		if _, ok := right[edge]; ok {
			count++
		}
	}
	return count
}

func countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaHeuristicSet map[Edge]struct{}, inCycleCover, inHighMsaHeuristic bool) int {
	count := 0
	for edge := range optimalEdges {
		_, cycleCoverContains := cycleCoverSet[edge]
		_, highMsaHeuristicContains := highMsaHeuristicSet[edge]
		if cycleCoverContains == inCycleCover && highMsaHeuristicContains == inHighMsaHeuristic {
			count++
		}
	}
	return count
}

func validateSquareMatrix(matrix [][]float64) error {
	dimension := len(matrix)
	if dimension == 0 {
		return fmt.Errorf("matrix is empty")
	}

	for i, row := range matrix {
		if len(row) != dimension {
			return fmt.Errorf("row %d has length %d, expected %d", i, len(row), dimension)
		}
	}

	return nil
}

func validateTours(tours map[string][]int, dimension int) error {
	for _, tourId := range sortedTourIds(tours) {
		tour := tours[tourId]
		if len(tour) != dimension {
			return fmt.Errorf("tour %s has length %d, expected %d", tourId, len(tour), dimension)
		}

		seen := make([]bool, dimension)
		for _, city := range tour {
			if city < 0 || city >= dimension {
				return fmt.Errorf("tour %s contains city %d outside [0, %d)", tourId, city, dimension)
			}
			if seen[city] {
				return fmt.Errorf("tour %s contains city %d more than once", tourId, city)
			}
			seen[city] = true
		}
	}

	return nil
}

func validateCycleCover(edges []Edge, dimension int) error {
	if len(edges) != dimension {
		return fmt.Errorf("has %d edges, expected %d", len(edges), dimension)
	}

	inDegree := make([]int, dimension)
	outDegree := make([]int, dimension)
	for _, edge := range edges {
		if edge.From < 0 || edge.From >= dimension || edge.To < 0 || edge.To >= dimension {
			return fmt.Errorf("edge %v is outside [0, %d)", edge, dimension)
		}
		if edge.From == edge.To {
			return fmt.Errorf("contains self-loop %v", edge)
		}
		outDegree[edge.From]++
		inDegree[edge.To]++
	}

	for vertex := 0; vertex < dimension; vertex++ {
		if inDegree[vertex] != 1 || outDegree[vertex] != 1 {
			return fmt.Errorf("vertex %d has in-degree %d and out-degree %d", vertex, inDegree[vertex], outDegree[vertex])
		}
	}

	return nil
}

func sortedTourIds(tours map[string][]int) []string {
	tourIds := make([]string, 0, len(tours))
	for tourId := range tours {
		tourIds = append(tourIds, tourId)
	}
	sort.Strings(tourIds)
	return tourIds
}

func normalizedMsaHeuristicValue(value float64, dimension int) float64 {
	if dimension <= 1 {
		return 0
	}
	return value / float64(dimension-1)
}
