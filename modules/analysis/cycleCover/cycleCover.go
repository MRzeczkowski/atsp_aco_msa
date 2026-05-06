package cycleCover

import (
	"atsp_aco_msa/modules/algorithms/compositeMsa"
	"atsp_aco_msa/modules/algorithms/hungarian"
	"atsp_aco_msa/modules/analysis/cmsaTours"
	"atsp_aco_msa/modules/models"
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

type Edge = models.Edge

type InstanceConfig struct {
	Name                        string
	Dimension                   int
	Matrix                      [][]float64
	KnownOptimal                float64
	CmsaDirectoryPath           string
	OptimalToursCsvPath         string
	CycleCoverEdgesCsvPath      string
	AnalysisCsvPath             string
	ThresholdsCsvPath           string
	CycleCoverOverlapMatrixPath string
}

type Config struct {
	Instances      []InstanceConfig
	SummaryCsvPath string
	ReportPath     string
	HighThreshold  float64
	Thresholds     []float64
}

type InstanceAnalysis struct {
	Instance          string
	Dimension         int
	Metrics           InstanceMetrics
	ThresholdsMetrics []ThresholdMetrics
}

type InstanceMetrics struct {
	FoundOptimalTourCount       int
	UniqueFoundOptimalEdgeCount int
	FoundOptimalEdgeDensity     float64

	CycleCoverCost               float64
	CycleCoverGapToKnownOptimal  float64
	CycleCoverCycleCount         int
	CycleCoverMinCycleLength     int
	CycleCoverAverageCycleLength float64
	CycleCoverMaxCycleLength     int

	CycleCoverMetrics             EdgeSetMetrics
	HighCmsaThreshold             float64
	HighCmsaMetrics               EdgeSetMetrics
	CycleCoverPositiveCmsaMetrics EdgeSetMetrics
	CycleCoverHighCmsaMetrics     EdgeSetMetrics

	CycleCoverEdgesWithPositiveCmsa int
	CycleCoverEdgesWithHighCmsa     int

	HighCmsaEdgesRemovedByCycleCover        int
	HighCmsaOptimalEdgesRemovedByCycleCover int
	HighCmsaPrecisionGainFromCycleCover     float64
	HighCmsaRecallLossFromCycleCover        float64

	OptimalEdgesInCycleCoverAndHighCmsa    int
	OptimalEdgesInCycleCoverNotHighCmsa    int
	OptimalEdgesInHighCmsaNotCycleCover    int
	OptimalEdgesInNeitherCycleCoverNorHigh int
}

type EdgeSetMetrics struct {
	EdgeCount          int
	EdgeDensity        float64
	OptimalEdgeCount   int
	FalsePositiveEdges int
	Precision          float64
	Recall             float64
	Lift               float64
}

type ThresholdMetrics struct {
	Threshold                               float64
	HighCmsa                                EdgeSetMetrics
	CycleCoverHighCmsa                      EdgeSetMetrics
	HighCmsaEdgesRemovedByCycleCover        int
	HighCmsaOptimalEdgesRemovedByCycleCover int
	HighCmsaPrecisionGainFromCycleCover     float64
	HighCmsaRecallLossFromCycleCover        float64
	CycleCoverShareOfHighCmsaEdges          float64
	CycleCoverShareOfHighCmsaOptimalEdges   float64
}

func AnalyzeInstances(config Config) ([]InstanceAnalysis, error) {
	thresholds := append([]float64(nil), config.Thresholds...)
	if len(thresholds) == 0 {
		thresholds = cmsaTours.DefaultThresholds()
	}
	sort.Float64s(thresholds)

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
		analysis, err := AnalyzeInstance(instance, thresholds, highThreshold)
		if err != nil {
			return nil, err
		}
		analyses = append(analyses, analysis)
	}

	if config.SummaryCsvPath != "" {
		if err := saveSummary(config.SummaryCsvPath, analyses); err != nil {
			return nil, err
		}
	}

	if config.ReportPath != "" {
		if err := saveReport(config.ReportPath, analyses, highThreshold); err != nil {
			return nil, err
		}
	}

	return analyses, nil
}

func AnalyzeInstance(config InstanceConfig, thresholds []float64, highThreshold float64) (InstanceAnalysis, error) {
	if len(thresholds) == 0 {
		thresholds = cmsaTours.DefaultThresholds()
	}
	thresholds = append([]float64(nil), thresholds...)
	sort.Float64s(thresholds)

	if highThreshold == 0 {
		highThreshold = 1.0
	}

	if err := validateSquareMatrix(config.Matrix); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid distance matrix: %w", config.Name, err)
	}
	if config.Dimension != 0 && config.Dimension != len(config.Matrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: configured dimension %d does not match matrix dimension %d", config.Name, config.Dimension, len(config.Matrix))
	}

	cmsa, err := compositeMsa.Read(config.CmsaDirectoryPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read CMSA: %w", config.Name, err)
	}
	if err := validateSquareMatrix(cmsa); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid CMSA: %w", config.Name, err)
	}
	if len(cmsa) != len(config.Matrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: CMSA dimension %d does not match matrix dimension %d", config.Name, len(cmsa), len(config.Matrix))
	}

	uniqueOptimalTours, err := cmsaTours.ReadOptimalTours(config.OptimalToursCsvPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read found optimal tours: %w", config.Name, err)
	}
	if err := validateTours(uniqueOptimalTours, len(config.Matrix)); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid found optimal tours: %w", config.Name, err)
	}

	cycleCoverEdges, cycleCoverCost, err := hungarian.MinimumCycleCover(config.Matrix)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to compute minimum cycle cover: %w", config.Name, err)
	}
	if err := validateCycleCover(cycleCoverEdges, len(config.Matrix)); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid minimum cycle cover: %w", config.Name, err)
	}

	analysis := calculateAnalysis(config.Name, len(config.Matrix), config.KnownOptimal, cmsa, uniqueOptimalTours, cycleCoverEdges, cycleCoverCost, thresholds, highThreshold)

	if err := saveCycleCoverEdges(config.CycleCoverEdgesCsvPath, config.Matrix, cmsa, buildTourEdgeSet(uniqueOptimalTours), buildEdgeSet(cycleCoverEdges), highThreshold); err != nil {
		return InstanceAnalysis{}, err
	}
	if err := saveInstanceAnalysis(config.AnalysisCsvPath, analysis); err != nil {
		return InstanceAnalysis{}, err
	}
	if err := saveThresholdAnalysis(config.ThresholdsCsvPath, analysis.ThresholdsMetrics); err != nil {
		return InstanceAnalysis{}, err
	}
	if config.CycleCoverOverlapMatrixPath != "" {
		matrix := buildCycleCoverCmsaOverlapMatrix(cmsa, buildEdgeSet(cycleCoverEdges))
		if err := saveMatrixCsv(config.CycleCoverOverlapMatrixPath, matrix); err != nil {
			return InstanceAnalysis{}, err
		}
	}

	return analysis, nil
}

func calculateAnalysis(instance string, dimension int, knownOptimal float64, cmsa [][]float64, uniqueOptimalTours map[string][]int, cycleCoverEdges []Edge, cycleCoverCost float64, thresholds []float64, highThreshold float64) InstanceAnalysis {
	optimalEdges := buildTourEdgeSet(uniqueOptimalTours)
	cycleCoverSet := buildEdgeSet(cycleCoverEdges)
	positiveCmsaSet := buildCmsaThresholdSet(cmsa, 0, false)
	highCmsaSet := buildCmsaThresholdSet(cmsa, highThreshold, true)
	cycleCoverPositiveCmsaSet := intersectEdgeSets(cycleCoverSet, positiveCmsaSet)
	cycleCoverHighCmsaSet := intersectEdgeSets(cycleCoverSet, highCmsaSet)

	totalDirectedEdges := dimension * (dimension - 1)
	cycleLengths := cycleLengths(cycleCoverEdges, dimension)

	metrics := InstanceMetrics{
		FoundOptimalTourCount:           len(uniqueOptimalTours),
		UniqueFoundOptimalEdgeCount:     len(optimalEdges),
		CycleCoverCost:                  cycleCoverCost,
		CycleCoverCycleCount:            len(cycleLengths),
		CycleCoverMinCycleLength:        minInt(cycleLengths),
		CycleCoverAverageCycleLength:    averageInt(cycleLengths),
		CycleCoverMaxCycleLength:        maxInt(cycleLengths),
		CycleCoverMetrics:               calculateEdgeSetMetrics(cycleCoverSet, optimalEdges, totalDirectedEdges),
		HighCmsaThreshold:               highThreshold,
		HighCmsaMetrics:                 calculateEdgeSetMetrics(highCmsaSet, optimalEdges, totalDirectedEdges),
		CycleCoverPositiveCmsaMetrics:   calculateEdgeSetMetrics(cycleCoverPositiveCmsaSet, optimalEdges, totalDirectedEdges),
		CycleCoverHighCmsaMetrics:       calculateEdgeSetMetrics(cycleCoverHighCmsaSet, optimalEdges, totalDirectedEdges),
		CycleCoverEdgesWithPositiveCmsa: len(cycleCoverPositiveCmsaSet),
		CycleCoverEdgesWithHighCmsa:     len(cycleCoverHighCmsaSet),
	}

	if totalDirectedEdges > 0 {
		metrics.FoundOptimalEdgeDensity = float64(len(optimalEdges)) / float64(totalDirectedEdges)
	}
	if knownOptimal > 0 {
		metrics.CycleCoverGapToKnownOptimal = (knownOptimal - cycleCoverCost) / knownOptimal
	}

	metrics.HighCmsaEdgesRemovedByCycleCover = metrics.HighCmsaMetrics.EdgeCount - metrics.CycleCoverHighCmsaMetrics.EdgeCount
	metrics.HighCmsaOptimalEdgesRemovedByCycleCover = metrics.HighCmsaMetrics.OptimalEdgeCount - metrics.CycleCoverHighCmsaMetrics.OptimalEdgeCount
	metrics.HighCmsaPrecisionGainFromCycleCover = metrics.CycleCoverHighCmsaMetrics.Precision - metrics.HighCmsaMetrics.Precision
	metrics.HighCmsaRecallLossFromCycleCover = metrics.HighCmsaMetrics.Recall - metrics.CycleCoverHighCmsaMetrics.Recall

	metrics.OptimalEdgesInCycleCoverAndHighCmsa = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highCmsaSet, true, true)
	metrics.OptimalEdgesInCycleCoverNotHighCmsa = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highCmsaSet, true, false)
	metrics.OptimalEdgesInHighCmsaNotCycleCover = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highCmsaSet, false, true)
	metrics.OptimalEdgesInNeitherCycleCoverNorHigh = countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highCmsaSet, false, false)

	thresholdMetrics := calculateThresholdMetrics(cmsa, optimalEdges, cycleCoverSet, thresholds, totalDirectedEdges)

	return InstanceAnalysis{
		Instance:          instance,
		Dimension:         dimension,
		Metrics:           metrics,
		ThresholdsMetrics: thresholdMetrics,
	}
}

func calculateThresholdMetrics(cmsa [][]float64, optimalEdges, cycleCoverSet map[Edge]struct{}, thresholds []float64, totalDirectedEdges int) []ThresholdMetrics {
	metrics := make([]ThresholdMetrics, 0, len(thresholds))
	sortedThresholds := append([]float64(nil), thresholds...)
	sort.Float64s(sortedThresholds)

	for _, threshold := range sortedThresholds {
		highCmsaSet := buildCmsaThresholdSet(cmsa, threshold, true)
		cycleCoverHighCmsaSet := intersectEdgeSets(cycleCoverSet, highCmsaSet)

		highCmsaMetrics := calculateEdgeSetMetrics(highCmsaSet, optimalEdges, totalDirectedEdges)
		cycleCoverHighCmsaMetrics := calculateEdgeSetMetrics(cycleCoverHighCmsaSet, optimalEdges, totalDirectedEdges)

		metric := ThresholdMetrics{
			Threshold:                               threshold,
			HighCmsa:                                highCmsaMetrics,
			CycleCoverHighCmsa:                      cycleCoverHighCmsaMetrics,
			HighCmsaEdgesRemovedByCycleCover:        highCmsaMetrics.EdgeCount - cycleCoverHighCmsaMetrics.EdgeCount,
			HighCmsaOptimalEdgesRemovedByCycleCover: highCmsaMetrics.OptimalEdgeCount - cycleCoverHighCmsaMetrics.OptimalEdgeCount,
			HighCmsaPrecisionGainFromCycleCover:     cycleCoverHighCmsaMetrics.Precision - highCmsaMetrics.Precision,
			HighCmsaRecallLossFromCycleCover:        highCmsaMetrics.Recall - cycleCoverHighCmsaMetrics.Recall,
		}

		if highCmsaMetrics.EdgeCount > 0 {
			metric.CycleCoverShareOfHighCmsaEdges = float64(cycleCoverHighCmsaMetrics.EdgeCount) / float64(highCmsaMetrics.EdgeCount)
		}
		if highCmsaMetrics.OptimalEdgeCount > 0 {
			metric.CycleCoverShareOfHighCmsaOptimalEdges = float64(cycleCoverHighCmsaMetrics.OptimalEdgeCount) / float64(highCmsaMetrics.OptimalEdgeCount)
		}

		metrics = append(metrics, metric)
	}

	return metrics
}

func calculateEdgeSetMetrics(edgeSet, optimalEdges map[Edge]struct{}, totalDirectedEdges int) EdgeSetMetrics {
	optimalEdgeCount := countIntersection(edgeSet, optimalEdges)
	metrics := EdgeSetMetrics{
		EdgeCount:        len(edgeSet),
		OptimalEdgeCount: optimalEdgeCount,
	}
	if len(optimalEdges) > 0 {
		metrics.FalsePositiveEdges = len(edgeSet) - optimalEdgeCount
	}

	if totalDirectedEdges > 0 {
		metrics.EdgeDensity = float64(len(edgeSet)) / float64(totalDirectedEdges)
	}
	if len(edgeSet) > 0 && len(optimalEdges) > 0 {
		metrics.Precision = float64(optimalEdgeCount) / float64(len(edgeSet))
	}
	if len(optimalEdges) > 0 {
		metrics.Recall = float64(optimalEdgeCount) / float64(len(optimalEdges))
	}
	if totalDirectedEdges > 0 && len(optimalEdges) > 0 && metrics.Precision > 0 {
		randomPrecision := float64(len(optimalEdges)) / float64(totalDirectedEdges)
		metrics.Lift = metrics.Precision / randomPrecision
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

func buildCmsaThresholdSet(cmsa [][]float64, threshold float64, includeEqual bool) map[Edge]struct{} {
	dimension := len(cmsa)
	set := make(map[Edge]struct{})
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			value := normalizedCmsaValue(cmsa[i][j], dimension)
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

func countOptimalEdgesByMembership(optimalEdges, cycleCoverSet, highCmsaSet map[Edge]struct{}, inCycleCover, inHighCmsa bool) int {
	count := 0
	for edge := range optimalEdges {
		_, cycleCoverContains := cycleCoverSet[edge]
		_, highCmsaContains := highCmsaSet[edge]
		if cycleCoverContains == inCycleCover && highCmsaContains == inHighCmsa {
			count++
		}
	}
	return count
}

func cycleLengths(edges []Edge, dimension int) []int {
	assignment := make([]int, dimension)
	for _, edge := range edges {
		assignment[edge.From] = edge.To
	}

	cycles, err := hungarian.AssignmentCycles(assignment)
	if err != nil {
		return nil
	}

	lengths := make([]int, len(cycles))
	for i, cycle := range cycles {
		lengths[i] = len(cycle)
	}
	sort.Ints(lengths)
	return lengths
}

func buildCycleCoverCmsaOverlapMatrix(cmsa [][]float64, cycleCoverSet map[Edge]struct{}) [][]float64 {
	dimension := len(cmsa)
	matrix := make([][]float64, dimension)
	for i := range matrix {
		matrix[i] = make([]float64, dimension)
	}

	for edge := range cycleCoverSet {
		matrix[edge.From][edge.To] = normalizedCmsaValue(cmsa[edge.From][edge.To], dimension)
	}

	return matrix
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

func sortedEdges(edgeSet map[Edge]struct{}) []Edge {
	edges := make([]Edge, 0, len(edgeSet))
	for edge := range edgeSet {
		edges = append(edges, edge)
	}
	sort.SliceStable(edges, func(i, j int) bool {
		if edges[i].From != edges[j].From {
			return edges[i].From < edges[j].From
		}
		return edges[i].To < edges[j].To
	})
	return edges
}

func normalizedCmsaValue(value float64, dimension int) float64 {
	if dimension <= 1 {
		return 0
	}
	return value / float64(dimension-1)
}

func minInt(values []int) int {
	if len(values) == 0 {
		return 0
	}
	return values[0]
}

func maxInt(values []int) int {
	if len(values) == 0 {
		return 0
	}
	return values[len(values)-1]
}

func averageInt(values []int) float64 {
	if len(values) == 0 {
		return 0
	}

	sum := 0
	for _, value := range values {
		sum += value
	}
	return float64(sum) / float64(len(values))
}

func ensureParentDirectory(path string) error {
	if path == "" {
		return nil
	}

	directory := filepath.Dir(path)
	if directory == "." {
		return nil
	}

	return os.MkdirAll(directory, 0700)
}

func formatFloat(value float64) string {
	return strconv.FormatFloat(value, 'f', 6, 64)
}

func formatPercent(value float64) string {
	return strconv.FormatFloat(100*value, 'f', 2, 64)
}

func saveCycleCoverEdges(path string, distances, cmsa [][]float64, optimalEdges, cycleCoverSet map[Edge]struct{}, highThreshold float64) error {
	if path == "" {
		return nil
	}
	if err := ensureParentDirectory(path); err != nil {
		return err
	}

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	header := []string{
		"From",
		"To",
		"Distance",
		"CMSA",
		"Normalized CMSA",
		"High CMSA",
		"Found-optimal edge",
	}
	if err := writer.Write(header); err != nil {
		return err
	}

	dimension := len(cmsa)
	for _, edge := range sortedEdges(cycleCoverSet) {
		normalizedCmsa := normalizedCmsaValue(cmsa[edge.From][edge.To], dimension)
		_, isOptimal := optimalEdges[edge]
		row := []string{
			strconv.Itoa(edge.From),
			strconv.Itoa(edge.To),
			formatFloat(distances[edge.From][edge.To]),
			formatFloat(cmsa[edge.From][edge.To]),
			formatFloat(normalizedCmsa),
			strconv.FormatBool(normalizedCmsa >= highThreshold),
			strconv.FormatBool(isOptimal),
		}
		if err := writer.Write(row); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func saveInstanceAnalysis(path string, analysis InstanceAnalysis) error {
	if path == "" {
		return nil
	}
	if err := ensureParentDirectory(path); err != nil {
		return err
	}

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	if err := writer.Write([]string{"Metric", "Value"}); err != nil {
		return err
	}

	metrics := analysis.Metrics
	rows := [][]string{
		{"Instance", analysis.Instance},
		{"Dimension", strconv.Itoa(analysis.Dimension)},
		{"Unique found-optimal tours", strconv.Itoa(metrics.FoundOptimalTourCount)},
		{"Unique found-optimal edges", strconv.Itoa(metrics.UniqueFoundOptimalEdgeCount)},
		{"Found-optimal edge density [%]", formatPercent(metrics.FoundOptimalEdgeDensity)},
		{"Cycle cover cost", formatFloat(metrics.CycleCoverCost)},
		{"Cycle cover gap to known optimum [%]", formatPercent(metrics.CycleCoverGapToKnownOptimal)},
		{"Cycle cover cycles", strconv.Itoa(metrics.CycleCoverCycleCount)},
		{"Cycle cover min cycle length", strconv.Itoa(metrics.CycleCoverMinCycleLength)},
		{"Cycle cover avg cycle length", formatFloat(metrics.CycleCoverAverageCycleLength)},
		{"Cycle cover max cycle length", strconv.Itoa(metrics.CycleCoverMaxCycleLength)},
	}
	rows = append(rows, edgeSetRows("Cycle cover", metrics.CycleCoverMetrics)...)
	rows = append(rows, []string{"High CMSA threshold", formatFloat(metrics.HighCmsaThreshold)})
	rows = append(rows, edgeSetRows("High CMSA", metrics.HighCmsaMetrics)...)
	rows = append(rows, edgeSetRows("Cycle cover positive-CMSA", metrics.CycleCoverPositiveCmsaMetrics)...)
	rows = append(rows, edgeSetRows("Cycle cover high-CMSA", metrics.CycleCoverHighCmsaMetrics)...)
	rows = append(rows,
		[]string{"Cycle cover edges with positive CMSA", strconv.Itoa(metrics.CycleCoverEdgesWithPositiveCmsa)},
		[]string{"Cycle cover edges with high CMSA", strconv.Itoa(metrics.CycleCoverEdgesWithHighCmsa)},
		[]string{"High-CMSA edges removed by cycle-cover gate", strconv.Itoa(metrics.HighCmsaEdgesRemovedByCycleCover)},
		[]string{"High-CMSA found-optimal edges removed by cycle-cover gate", strconv.Itoa(metrics.HighCmsaOptimalEdgesRemovedByCycleCover)},
		[]string{"High-CMSA precision gain from cycle-cover gate", formatFloat(metrics.HighCmsaPrecisionGainFromCycleCover)},
		[]string{"High-CMSA recall loss from cycle-cover gate", formatFloat(metrics.HighCmsaRecallLossFromCycleCover)},
		[]string{"Found-optimal edges in both cycle cover and high CMSA", strconv.Itoa(metrics.OptimalEdgesInCycleCoverAndHighCmsa)},
		[]string{"Found-optimal edges in cycle cover but not high CMSA", strconv.Itoa(metrics.OptimalEdgesInCycleCoverNotHighCmsa)},
		[]string{"Found-optimal edges in high CMSA but not cycle cover", strconv.Itoa(metrics.OptimalEdgesInHighCmsaNotCycleCover)},
		[]string{"Found-optimal edges in neither cycle cover nor high CMSA", strconv.Itoa(metrics.OptimalEdgesInNeitherCycleCoverNorHigh)},
	)

	for _, row := range rows {
		if err := writer.Write(row); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func edgeSetRows(prefix string, metrics EdgeSetMetrics) [][]string {
	return [][]string{
		{prefix + " edges", strconv.Itoa(metrics.EdgeCount)},
		{prefix + " edge density [%]", formatPercent(metrics.EdgeDensity)},
		{prefix + " found-optimal edges", strconv.Itoa(metrics.OptimalEdgeCount)},
		{prefix + " edges not seen in found-optimal tours", strconv.Itoa(metrics.FalsePositiveEdges)},
		{prefix + " precision vs found-optimal edges [%]", formatPercent(metrics.Precision)},
		{prefix + " recall vs found-optimal edges [%]", formatPercent(metrics.Recall)},
		{prefix + " lift vs found-optimal edge density", formatFloat(metrics.Lift)},
	}
}

func saveThresholdAnalysis(path string, metrics []ThresholdMetrics) error {
	if path == "" {
		return nil
	}
	if err := ensureParentDirectory(path); err != nil {
		return err
	}

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	header := []string{
		"Threshold",
		"High CMSA edges",
		"High CMSA precision [%]",
		"High CMSA recall [%]",
		"High CMSA lift",
		"Cycle-cover high-CMSA edges",
		"Cycle-cover high-CMSA precision [%]",
		"Cycle-cover high-CMSA recall [%]",
		"Cycle-cover high-CMSA lift",
		"High-CMSA edges removed by cycle cover",
		"High-CMSA found-optimal edges removed by cycle cover",
		"Precision gain from cycle cover",
		"Recall loss from cycle cover",
		"Cycle-cover share of high-CMSA edges [%]",
		"Cycle-cover share of high-CMSA found-optimal edges [%]",
	}
	if err := writer.Write(header); err != nil {
		return err
	}

	for _, metric := range metrics {
		row := []string{
			formatFloat(metric.Threshold),
			strconv.Itoa(metric.HighCmsa.EdgeCount),
			formatPercent(metric.HighCmsa.Precision),
			formatPercent(metric.HighCmsa.Recall),
			formatFloat(metric.HighCmsa.Lift),
			strconv.Itoa(metric.CycleCoverHighCmsa.EdgeCount),
			formatPercent(metric.CycleCoverHighCmsa.Precision),
			formatPercent(metric.CycleCoverHighCmsa.Recall),
			formatFloat(metric.CycleCoverHighCmsa.Lift),
			strconv.Itoa(metric.HighCmsaEdgesRemovedByCycleCover),
			strconv.Itoa(metric.HighCmsaOptimalEdgesRemovedByCycleCover),
			formatFloat(metric.HighCmsaPrecisionGainFromCycleCover),
			formatFloat(metric.HighCmsaRecallLossFromCycleCover),
			formatPercent(metric.CycleCoverShareOfHighCmsaEdges),
			formatPercent(metric.CycleCoverShareOfHighCmsaOptimalEdges),
		}
		if err := writer.Write(row); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func saveSummary(path string, analyses []InstanceAnalysis) error {
	if err := ensureParentDirectory(path); err != nil {
		return err
	}

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	header := []string{
		"Instance",
		"Dimension",
		"Unique found-optimal tours",
		"Unique found-optimal edges",
		"Found-optimal edge density [%]",
		"Cycle cover cost",
		"Cycle cover gap to known optimum [%]",
		"Cycle cover cycles",
		"Cycle cover precision [%]",
		"Cycle cover recall [%]",
		"Cycle cover lift",
		"High CMSA threshold",
		"High CMSA edges",
		"High CMSA precision [%]",
		"High CMSA recall [%]",
		"High CMSA lift",
		"Cycle-cover high-CMSA edges",
		"Cycle-cover high-CMSA precision [%]",
		"Cycle-cover high-CMSA recall [%]",
		"Cycle-cover high-CMSA lift",
		"High-CMSA edges removed by cycle-cover gate",
		"High-CMSA found-optimal edges removed by cycle-cover gate",
		"High-CMSA precision gain from cycle-cover gate",
		"High-CMSA recall loss from cycle-cover gate",
		"Found-optimal edges in both cycle cover and high CMSA",
		"Found-optimal edges in cycle cover but not high CMSA",
		"Found-optimal edges in high CMSA but not cycle cover",
		"Found-optimal edges in neither cycle cover nor high CMSA",
	}
	if err := writer.Write(header); err != nil {
		return err
	}

	for _, analysis := range analyses {
		metrics := analysis.Metrics
		row := []string{
			analysis.Instance,
			strconv.Itoa(analysis.Dimension),
			strconv.Itoa(metrics.FoundOptimalTourCount),
			strconv.Itoa(metrics.UniqueFoundOptimalEdgeCount),
			formatPercent(metrics.FoundOptimalEdgeDensity),
			formatFloat(metrics.CycleCoverCost),
			formatPercent(metrics.CycleCoverGapToKnownOptimal),
			strconv.Itoa(metrics.CycleCoverCycleCount),
			formatPercent(metrics.CycleCoverMetrics.Precision),
			formatPercent(metrics.CycleCoverMetrics.Recall),
			formatFloat(metrics.CycleCoverMetrics.Lift),
			formatFloat(metrics.HighCmsaThreshold),
			strconv.Itoa(metrics.HighCmsaMetrics.EdgeCount),
			formatPercent(metrics.HighCmsaMetrics.Precision),
			formatPercent(metrics.HighCmsaMetrics.Recall),
			formatFloat(metrics.HighCmsaMetrics.Lift),
			strconv.Itoa(metrics.CycleCoverHighCmsaMetrics.EdgeCount),
			formatPercent(metrics.CycleCoverHighCmsaMetrics.Precision),
			formatPercent(metrics.CycleCoverHighCmsaMetrics.Recall),
			formatFloat(metrics.CycleCoverHighCmsaMetrics.Lift),
			strconv.Itoa(metrics.HighCmsaEdgesRemovedByCycleCover),
			strconv.Itoa(metrics.HighCmsaOptimalEdgesRemovedByCycleCover),
			formatFloat(metrics.HighCmsaPrecisionGainFromCycleCover),
			formatFloat(metrics.HighCmsaRecallLossFromCycleCover),
			strconv.Itoa(metrics.OptimalEdgesInCycleCoverAndHighCmsa),
			strconv.Itoa(metrics.OptimalEdgesInCycleCoverNotHighCmsa),
			strconv.Itoa(metrics.OptimalEdgesInHighCmsaNotCycleCover),
			strconv.Itoa(metrics.OptimalEdgesInNeitherCycleCoverNorHigh),
		}
		if err := writer.Write(row); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func saveReport(path string, analyses []InstanceAnalysis, highThreshold float64) error {
	if err := ensureParentDirectory(path); err != nil {
		return err
	}

	var instancesWithTours int
	var cycleCoverPrecisionSum, cycleCoverRecallSum, highCmsaPrecisionSum, highCmsaRecallSum float64
	var gatedPrecisionSum, gatedRecallSum, precisionGainSum, recallLossSum float64

	for _, analysis := range analyses {
		if analysis.Metrics.FoundOptimalTourCount == 0 {
			continue
		}
		instancesWithTours++
		metrics := analysis.Metrics
		cycleCoverPrecisionSum += metrics.CycleCoverMetrics.Precision
		cycleCoverRecallSum += metrics.CycleCoverMetrics.Recall
		highCmsaPrecisionSum += metrics.HighCmsaMetrics.Precision
		highCmsaRecallSum += metrics.HighCmsaMetrics.Recall
		gatedPrecisionSum += metrics.CycleCoverHighCmsaMetrics.Precision
		gatedRecallSum += metrics.CycleCoverHighCmsaMetrics.Recall
		precisionGainSum += metrics.HighCmsaPrecisionGainFromCycleCover
		recallLossSum += metrics.HighCmsaRecallLossFromCycleCover
	}

	var builder strings.Builder
	builder.WriteString("# Cycle-Cover/CMSA Analysis\n\n")
	builder.WriteString(fmt.Sprintf("Main high-CMSA threshold: %.1f\n\n", highThreshold))
	builder.WriteString("## Summary\n\n")
	builder.WriteString(fmt.Sprintf("- Instances analyzed: %d\n", len(analyses)))
	builder.WriteString(fmt.Sprintf("- Instances with found optimal tours: %d\n", instancesWithTours))
	if instancesWithTours > 0 {
		count := float64(instancesWithTours)
		builder.WriteString(fmt.Sprintf("- Average cycle-cover precision: %s%%\n", formatPercent(cycleCoverPrecisionSum/count)))
		builder.WriteString(fmt.Sprintf("- Average cycle-cover recall: %s%%\n", formatPercent(cycleCoverRecallSum/count)))
		builder.WriteString(fmt.Sprintf("- Average high-CMSA precision: %s%%\n", formatPercent(highCmsaPrecisionSum/count)))
		builder.WriteString(fmt.Sprintf("- Average high-CMSA recall: %s%%\n", formatPercent(highCmsaRecallSum/count)))
		builder.WriteString(fmt.Sprintf("- Average cycle-cover-gated high-CMSA precision: %s%%\n", formatPercent(gatedPrecisionSum/count)))
		builder.WriteString(fmt.Sprintf("- Average cycle-cover-gated high-CMSA recall: %s%%\n", formatPercent(gatedRecallSum/count)))
		builder.WriteString(fmt.Sprintf("- Average precision gain from cycle-cover gate: %s\n", formatFloat(precisionGainSum/count)))
		builder.WriteString(fmt.Sprintf("- Average recall loss from cycle-cover gate: %s\n", formatFloat(recallLossSum/count)))
	}
	builder.WriteString("\n")

	builder.WriteString("## Notes\n\n")
	builder.WriteString("- `Cycle cover` is the minimum assignment/cycle-cover solution computed from the original ATSP matrix with self-loops forbidden.\n")
	builder.WriteString("- `High CMSA` uses normalized `CMSA / (dimension - 1)` and the threshold shown above.\n")
	builder.WriteString("- `Cycle-cover high-CMSA` is the strict intersection of high-CMSA edges and minimum-cycle-cover edges.\n")
	builder.WriteString("- Precision/recall use the union of tours recorded in `solutions.csv`; absent edges are not proven non-optimal.\n")
	builder.WriteString("- Precision gain and recall loss compare `cycle-cover high-CMSA` against high-CMSA alone.\n\n")

	builder.WriteString("## Instances\n\n")
	builder.WriteString("| Instance | Tours | Opt edges | CC gap % | CC cycles | CC precision % | CC recall % | High precision % | High recall % | Gated precision % | Gated recall % | Precision gain | Recall loss |\n")
	builder.WriteString("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
	for _, analysis := range analyses {
		metrics := analysis.Metrics
		builder.WriteString(fmt.Sprintf("| %s | %d | %d | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |\n",
			analysis.Instance,
			metrics.FoundOptimalTourCount,
			metrics.UniqueFoundOptimalEdgeCount,
			formatPercent(metrics.CycleCoverGapToKnownOptimal),
			metrics.CycleCoverCycleCount,
			formatPercent(metrics.CycleCoverMetrics.Precision),
			formatPercent(metrics.CycleCoverMetrics.Recall),
			formatPercent(metrics.HighCmsaMetrics.Precision),
			formatPercent(metrics.HighCmsaMetrics.Recall),
			formatPercent(metrics.CycleCoverHighCmsaMetrics.Precision),
			formatPercent(metrics.CycleCoverHighCmsaMetrics.Recall),
			formatFloat(metrics.HighCmsaPrecisionGainFromCycleCover),
			formatFloat(metrics.HighCmsaRecallLossFromCycleCover),
		))
	}

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func saveMatrixCsv(path string, matrix [][]float64) error {
	if path == "" {
		return nil
	}
	if err := ensureParentDirectory(path); err != nil {
		return err
	}

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	for _, row := range matrix {
		record := make([]string, len(row))
		for i, value := range row {
			record[i] = formatFloat(value)
		}
		if err := writer.Write(record); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}
