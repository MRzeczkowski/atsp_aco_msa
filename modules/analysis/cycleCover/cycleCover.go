package cycleCover

import (
	"atsp_aco_msa/modules/algorithms/hungarian"
	"atsp_aco_msa/modules/algorithms/msaSupport"
	"atsp_aco_msa/modules/analysis/msaSupportTours"
	"atsp_aco_msa/modules/models"
	"fmt"
	"sort"
)

type Edge = models.Edge

type InstanceConfig struct {
	Name                    string
	Dimension               int
	Matrix                  [][]float64
	MsaSupportDirectoryPath string
	OptimalToursCsvPath     string
}

type Config struct {
	Instances     []InstanceConfig
	HighThreshold float64
}

type InstanceAnalysis struct {
	Instance  string
	Dimension int
	Metrics   InstanceMetrics
}

type InstanceMetrics struct {
	FoundOptimalTourCount       int
	UniqueFoundOptimalEdgeCount int

	CycleCoverMetrics      EdgeSetMetrics
	HighMsaSupportMetrics  EdgeSetMetrics
	CycleCoverHighMsaEdges int

	OptimalEdgesInCycleCoverAndHighMsaSupport int
	OptimalEdgesInCycleCoverNotHighMsaSupport int
	OptimalEdgesInHighMsaSupportNotCycleCover int
	OptimalEdgesInNeitherCycleCoverNorHigh    int
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

	instances := append([]InstanceConfig(nil), config.Instances...)
	sort.SliceStable(instances, func(i, j int) bool {
		return instances[i].Name < instances[j].Name
	})

	analyses := make([]InstanceAnalysis, 0, len(instances))
	for _, instance := range instances {
		analysis, err := AnalyzeInstance(instance, highThreshold)
		if err != nil {
			return nil, err
		}
		analyses = append(analyses, analysis)
	}

	return analyses, nil
}

func AnalyzeInstance(config InstanceConfig, highThreshold float64) (InstanceAnalysis, error) {
	if highThreshold == 0 {
		highThreshold = 1.0
	}

	if err := validateSquareMatrix(config.Matrix); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid distance matrix: %w", config.Name, err)
	}
	if config.Dimension != 0 && config.Dimension != len(config.Matrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: configured dimension %d does not match matrix dimension %d", config.Name, config.Dimension, len(config.Matrix))
	}

	msaSupportMatrix, err := msaSupport.Read(config.MsaSupportDirectoryPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read MSA support: %w", config.Name, err)
	}
	if err := validateSquareMatrix(msaSupportMatrix); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid MSA support: %w", config.Name, err)
	}
	if len(msaSupportMatrix) != len(config.Matrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: MSA support dimension %d does not match matrix dimension %d", config.Name, len(msaSupportMatrix), len(config.Matrix))
	}

	uniqueOptimalTours, err := msaSupportTours.ReadOptimalTours(config.OptimalToursCsvPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read found optimal tours: %w", config.Name, err)
	}
	if err := validateTours(uniqueOptimalTours, len(config.Matrix)); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid found optimal tours: %w", config.Name, err)
	}

	cycleCoverEdges, _, err := hungarian.MinimumCycleCover(config.Matrix)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to compute minimum cycle cover: %w", config.Name, err)
	}
	if err := validateCycleCover(cycleCoverEdges, len(config.Matrix)); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid minimum cycle cover: %w", config.Name, err)
	}

	analysis := calculateAnalysis(config.Name, len(config.Matrix), msaSupportMatrix, uniqueOptimalTours, cycleCoverEdges, highThreshold)
	return analysis, nil
}

func calculateAnalysis(instance string, dimension int, msaSupport [][]float64, uniqueOptimalTours map[string][]int, cycleCoverEdges []Edge, highThreshold float64) InstanceAnalysis {
	optimalEdges := buildTourEdgeSet(uniqueOptimalTours)
	cycleCoverSet := buildEdgeSet(cycleCoverEdges)
	highMsaSupportSet := buildMsaSupportThresholdSet(msaSupport, highThreshold, true)
	cycleCoverHighMsaSupportSet := intersectEdgeSets(cycleCoverSet, highMsaSupportSet)

	metrics := InstanceMetrics{
		FoundOptimalTourCount:       len(uniqueOptimalTours),
		UniqueFoundOptimalEdgeCount: len(optimalEdges),
		CycleCoverMetrics:           calculateEdgeSetMetrics(cycleCoverSet, optimalEdges),
		HighMsaSupportMetrics:       calculateEdgeSetMetrics(highMsaSupportSet, optimalEdges),
		CycleCoverHighMsaEdges:      len(cycleCoverHighMsaSupportSet),
	}

	metrics.OptimalEdgesInCycleCoverAndHighMsaSupport = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaSupportSet, true, true)
	metrics.OptimalEdgesInCycleCoverNotHighMsaSupport = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaSupportSet, true, false)
	metrics.OptimalEdgesInHighMsaSupportNotCycleCover = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaSupportSet, false, true)
	metrics.OptimalEdgesInNeitherCycleCoverNorHigh = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaSupportSet, false, false)

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

func buildMsaSupportThresholdSet(msaSupport [][]float64, threshold float64, includeEqual bool) map[Edge]struct{} {
	dimension := len(msaSupport)
	set := make(map[Edge]struct{})
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			value := normalizedMsaSupportValue(msaSupport[i][j], dimension)
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

func countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highMsaSupportSet map[Edge]struct{}, inCycleCover, inHighMsaSupport bool) int {
	count := 0
	for edge := range optimalEdges {
		_, cycleCoverContains := cycleCoverSet[edge]
		_, highMsaSupportContains := highMsaSupportSet[edge]
		if cycleCoverContains == inCycleCover && highMsaSupportContains == inHighMsaSupport {
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

func normalizedMsaSupportValue(value float64, dimension int) float64 {
	if dimension <= 1 {
		return 0
	}
	return value / float64(dimension-1)
}
