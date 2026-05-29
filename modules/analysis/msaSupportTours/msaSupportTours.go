package msaSupportTours

import (
	"atsp_aco_msa/modules/algorithms/msaSupport"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

type Edge = models.Edge

type InstanceConfig struct {
	Name                              string
	Dimension                         int
	MsaSupportDirectoryPath           string
	OptimalToursCsvPath               string
	ToursHeatmapPath                  string
	ToursHistogramPath                string
	MsaSupportToursOverlapHeatmapPath string
	AnalysisCsvPath                   string
	ThresholdsCsvPath                 string
}

type Config struct {
	Instances      []InstanceConfig
	SummaryCsvPath string
	ReportPath     string
	HighThreshold  float64
	Thresholds     []float64
}

type TourStatistics struct {
	tourId                                                                                             string
	commonalityWithMsaSupport, minCommonalityWithMsa, averageCommonalityWithMsa, maxCommonalityWithMsa float64
}

type InstanceAnalysis struct {
	Instance          string
	Dimension         int
	Metrics           InstanceMetrics
	ThresholdsMetrics []ThresholdMetrics
}

type InstanceMetrics struct {
	FoundOptimalTourCount                         int
	UniqueFoundOptimalEdgeCount                   int
	FoundOptimalEdgeDensity                       float64
	MsaSupportPositiveEdgeCount                   int
	MsaSupportPositiveDensity                     float64
	HighMsaSupportThreshold                       float64
	HighMsaSupportEdgeCount                       int
	HighMsaSupportDensity                         float64
	OptimalEdgeMsaSupportCoverage                 float64
	HighMsaSupportPrecision                       float64
	HighMsaSupportRecall                          float64
	HighMsaSupportLift                            float64
	TopMsaSupportEdgeCount                        int
	TopMsaSupportDensity                          float64
	TopMsaSupportPrecision                        float64
	TopMsaSupportRecall                           float64
	TopMsaSupportLift                             float64
	TopMsaSupportFalsePositiveEdges               int
	TopMsaSupportMinTourCoverage                  float64
	TopMsaSupportAverageTourCoverage              float64
	TopMsaSupportMedianTourCoverage               float64
	TopMsaSupportMaxTourCoverage                  float64
	AverageNormalizedMsaSupportOnOptimalEdges     float64
	MedianNormalizedMsaSupportOnOptimalEdges      float64
	AverageNormalizedMsaSupportOnNotFoundPositive float64
	MedianNormalizedMsaSupportOnNotFoundPositive  float64
	MissingFoundOptimalEdges                      int
	FalsePositiveHighMsaSupportEdges              int
}

type ThresholdMetrics struct {
	Threshold               float64
	HighMsaSupportEdgeCount int
	HighMsaSupportDensity   float64
	Precision               float64
	Recall                  float64
	Lift                    float64
	FalsePositiveHighEdges  int
}

func DefaultThresholds() []float64 {
	thresholds := make([]float64, 0, 10)
	for i := 1; i <= 10; i++ {
		thresholds = append(thresholds, float64(i)/10.0)
	}
	return thresholds
}

func AnalyzeInstances(config Config) ([]InstanceAnalysis, error) {
	thresholds := append([]float64(nil), config.Thresholds...)
	if len(thresholds) == 0 {
		thresholds = DefaultThresholds()
	}
	sort.Float64s(thresholds)

	highThreshold := config.HighThreshold
	if highThreshold == 0 {
		highThreshold = 0.8
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
		thresholds = DefaultThresholds()
	}
	thresholds = append([]float64(nil), thresholds...)
	sort.Float64s(thresholds)

	if highThreshold == 0 {
		highThreshold = 0.8
	}

	uniqueOptimalTours, err := ReadOptimalTours(config.OptimalToursCsvPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read optimal tours: %w", config.Name, err)
	}

	msaSupportMatrix, err := msaSupport.Read(config.MsaSupportDirectoryPath)
	if err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: failed to read MSA support: %w", config.Name, err)
	}
	if err := validateSquareMatrix(msaSupportMatrix); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid MSA support: %w", config.Name, err)
	}
	if config.Dimension != 0 && config.Dimension != len(msaSupportMatrix) {
		return InstanceAnalysis{}, fmt.Errorf("%s: configured dimension %d does not match MSA support dimension %d", config.Name, config.Dimension, len(msaSupportMatrix))
	}
	if err := validateTours(uniqueOptimalTours, len(msaSupportMatrix)); err != nil {
		return InstanceAnalysis{}, fmt.Errorf("%s: invalid found optimal tours: %w", config.Name, err)
	}

	metrics := CalculateMetrics(msaSupportMatrix, uniqueOptimalTours, highThreshold)
	thresholdMetrics := CalculateThresholdMetrics(msaSupportMatrix, uniqueOptimalTours, thresholds)

	if err := saveInstanceAnalysis(config.AnalysisCsvPath, config.Name, config.Dimension, metrics); err != nil {
		return InstanceAnalysis{}, err
	}

	if err := saveThresholdAnalysis(config.ThresholdsCsvPath, thresholdMetrics); err != nil {
		return InstanceAnalysis{}, err
	}

	toursCount := len(uniqueOptimalTours)
	if toursCount > 0 {
		toursMatrix := BuildToursMatrix(uniqueOptimalTours, config.Dimension)
		toursHeatmapPlotTitle := config.Name + " tours heatmap"
		if err := utilities.SaveHeatmapFromMatrix(toursMatrix, toursHeatmapPlotTitle, config.ToursHeatmapPath); err != nil {
			return InstanceAnalysis{}, err
		}

		overlapMatrix := BuildMsaSupportToursOverlapMatrix(msaSupportMatrix, toursMatrix, toursCount)
		overlapHeatmapPlotTitle := config.Name + " MSA support/tours overlap heatmap"
		if err := utilities.SaveHeatmapFromMatrix(overlapMatrix, overlapHeatmapPlotTitle, config.MsaSupportToursOverlapHeatmapPath); err != nil {
			return InstanceAnalysis{}, err
		}

		if toursCount > 1 {
			dataForHistogram := filterZeroes(flattenMatrix(toursMatrix))
			if len(dataForHistogram) > 0 {
				toursHistogramPlotTitle := config.Name + " tours histogram"
				numberOfBins := countUniqueValues(dataForHistogram)

				if err := utilities.SaveHistogramFromData(dataForHistogram, numberOfBins, toursHistogramPlotTitle, config.ToursHistogramPath); err != nil {
					return InstanceAnalysis{}, err
				}
			}
		} else if err := removeIfExists(config.ToursHistogramPath); err != nil {
			return InstanceAnalysis{}, err
		}
	} else if err := removeStaleTourPlots(config); err != nil {
		return InstanceAnalysis{}, err
	}

	return InstanceAnalysis{
		Instance:          config.Name,
		Dimension:         config.Dimension,
		Metrics:           metrics,
		ThresholdsMetrics: thresholdMetrics,
	}, nil
}

func CalculateMetrics(msaSupport [][]float64, uniqueOptimalTours map[string][]int, highThreshold float64) InstanceMetrics {
	thresholdMetrics := calculateMetrics(msaSupport, uniqueOptimalTours, highThreshold)
	return thresholdMetrics
}

func CalculateThresholdMetrics(msaSupport [][]float64, uniqueOptimalTours map[string][]int, thresholds []float64) []ThresholdMetrics {
	metrics := make([]ThresholdMetrics, 0, len(thresholds))
	sortedThresholds := append([]float64(nil), thresholds...)
	sort.Float64s(sortedThresholds)

	for _, threshold := range sortedThresholds {
		instanceMetrics := calculateMetrics(msaSupport, uniqueOptimalTours, threshold)
		metrics = append(metrics, ThresholdMetrics{
			Threshold:               threshold,
			HighMsaSupportEdgeCount: instanceMetrics.HighMsaSupportEdgeCount,
			HighMsaSupportDensity:   instanceMetrics.HighMsaSupportDensity,
			Precision:               instanceMetrics.HighMsaSupportPrecision,
			Recall:                  instanceMetrics.HighMsaSupportRecall,
			Lift:                    instanceMetrics.HighMsaSupportLift,
			FalsePositiveHighEdges:  instanceMetrics.FalsePositiveHighMsaSupportEdges,
		})
	}

	return metrics
}

func calculateMetrics(msaSupport [][]float64, uniqueOptimalTours map[string][]int, highThreshold float64) InstanceMetrics {
	dimension := len(msaSupport)
	totalDirectedEdges := dimension * (dimension - 1)
	tourEdges := buildTourEdgeSet(uniqueOptimalTours)
	foundOptimalEdgeCount := len(tourEdges)
	hasFoundOptimalEdges := foundOptimalEdgeCount > 0

	positiveMsaSupportEdges := 0
	highMsaSupportEdges := 0
	highMsaSupportOptimalEdges := 0
	positiveOptimalEdges := 0
	falsePositiveHighEdges := 0
	optimalValues := make([]float64, 0, foundOptimalEdgeCount)
	notFoundPositiveValues := make([]float64, 0)
	topMsaSupportEdges := buildTopMsaSupportEdgeSet(msaSupport, dimension-1)
	topMsaSupportOptimalEdges := countIntersection(topMsaSupportEdges, tourEdges)
	topMsaSupportFalsePositiveEdges := 0
	if hasFoundOptimalEdges {
		topMsaSupportFalsePositiveEdges = len(topMsaSupportEdges) - topMsaSupportOptimalEdges
	}
	topMsaSupportTourCoverages := calculateTourCoverages(uniqueOptimalTours, topMsaSupportEdges)

	for i := 0; i < dimension; i++ {
		for j := 0; j < len(msaSupport[i]); j++ {
			if i == j {
				continue
			}

			value := normalizedMsaSupportValue(msaSupport[i][j], dimension)
			edge := Edge{From: i, To: j}
			_, isOptimalEdge := tourEdges[edge]

			if isOptimalEdge {
				optimalValues = append(optimalValues, value)
				if msaSupport[i][j] > 0 {
					positiveOptimalEdges++
				}
			}

			if msaSupport[i][j] <= 0 {
				continue
			}

			positiveMsaSupportEdges++
			if hasFoundOptimalEdges && !isOptimalEdge {
				notFoundPositiveValues = append(notFoundPositiveValues, value)
			}

			if value >= highThreshold {
				highMsaSupportEdges++
				if isOptimalEdge {
					highMsaSupportOptimalEdges++
				} else if hasFoundOptimalEdges {
					falsePositiveHighEdges++
				}
			}
		}
	}

	metrics := InstanceMetrics{
		FoundOptimalTourCount:                         len(uniqueOptimalTours),
		UniqueFoundOptimalEdgeCount:                   foundOptimalEdgeCount,
		MsaSupportPositiveEdgeCount:                   positiveMsaSupportEdges,
		HighMsaSupportThreshold:                       highThreshold,
		HighMsaSupportEdgeCount:                       highMsaSupportEdges,
		TopMsaSupportEdgeCount:                        len(topMsaSupportEdges),
		TopMsaSupportFalsePositiveEdges:               topMsaSupportFalsePositiveEdges,
		TopMsaSupportMinTourCoverage:                  minFloat(topMsaSupportTourCoverages),
		TopMsaSupportAverageTourCoverage:              average(topMsaSupportTourCoverages),
		TopMsaSupportMedianTourCoverage:               median(topMsaSupportTourCoverages),
		TopMsaSupportMaxTourCoverage:                  maxFloat(topMsaSupportTourCoverages),
		AverageNormalizedMsaSupportOnOptimalEdges:     average(optimalValues),
		MedianNormalizedMsaSupportOnOptimalEdges:      median(optimalValues),
		AverageNormalizedMsaSupportOnNotFoundPositive: average(notFoundPositiveValues),
		MedianNormalizedMsaSupportOnNotFoundPositive:  median(notFoundPositiveValues),
		MissingFoundOptimalEdges:                      foundOptimalEdgeCount - positiveOptimalEdges,
		FalsePositiveHighMsaSupportEdges:              falsePositiveHighEdges,
	}

	if totalDirectedEdges > 0 {
		metrics.FoundOptimalEdgeDensity = float64(foundOptimalEdgeCount) / float64(totalDirectedEdges)
		metrics.MsaSupportPositiveDensity = float64(positiveMsaSupportEdges) / float64(totalDirectedEdges)
		metrics.HighMsaSupportDensity = float64(highMsaSupportEdges) / float64(totalDirectedEdges)
		metrics.TopMsaSupportDensity = float64(len(topMsaSupportEdges)) / float64(totalDirectedEdges)
	}

	if foundOptimalEdgeCount > 0 {
		metrics.OptimalEdgeMsaSupportCoverage = float64(positiveOptimalEdges) / float64(foundOptimalEdgeCount)
		metrics.HighMsaSupportRecall = float64(highMsaSupportOptimalEdges) / float64(foundOptimalEdgeCount)
		metrics.TopMsaSupportRecall = float64(topMsaSupportOptimalEdges) / float64(foundOptimalEdgeCount)
	}

	if highMsaSupportEdges > 0 && foundOptimalEdgeCount > 0 {
		metrics.HighMsaSupportPrecision = float64(highMsaSupportOptimalEdges) / float64(highMsaSupportEdges)
	}
	if len(topMsaSupportEdges) > 0 && foundOptimalEdgeCount > 0 {
		metrics.TopMsaSupportPrecision = float64(topMsaSupportOptimalEdges) / float64(len(topMsaSupportEdges))
	}

	if totalDirectedEdges > 0 && foundOptimalEdgeCount > 0 && metrics.HighMsaSupportPrecision > 0 {
		randomPrecision := float64(foundOptimalEdgeCount) / float64(totalDirectedEdges)
		metrics.HighMsaSupportLift = metrics.HighMsaSupportPrecision / randomPrecision
	}
	if totalDirectedEdges > 0 && foundOptimalEdgeCount > 0 && metrics.TopMsaSupportPrecision > 0 {
		randomPrecision := float64(foundOptimalEdgeCount) / float64(totalDirectedEdges)
		metrics.TopMsaSupportLift = metrics.TopMsaSupportPrecision / randomPrecision
	}

	return metrics
}

func ReadOptimalTours(optimalUniqueToursCsvPath string) (map[string][]int, error) {
	result := make(map[string][]int)

	file, err := os.Open(optimalUniqueToursCsvPath)
	if err != nil {
		if os.IsNotExist(err) {
			return result, nil
		}
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	rows, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("failed to read file: %w", err)
	}

	for i, row := range rows {
		if i == 0 {
			continue
		}

		if len(row) < 1 {
			return nil, fmt.Errorf("invalid row format at line %d", i+1)
		}

		var value []int
		if err := json.Unmarshal([]byte(row[0]), &value); err != nil {
			return nil, fmt.Errorf("failed to parse JSON at line %d: %w", i+1, err)
		}

		result[row[0]] = value
	}

	return result, nil
}

func AddUniqueTour(savedTours map[string][]int, tour []int) {
	normalizedTour := normalizeTour(tour)
	keyJson, _ := json.Marshal(normalizedTour)
	key := string(keyJson)

	if _, exists := savedTours[key]; exists {
		return
	}

	savedTours[key] = tour
}

func SaveOptimalToursStatistics(optimalUniqueToursCsvPath, msaSupportDir string, uniqueOptimalTours map[string][]int) error {
	toursStatistics, err := calculateToursStatistics(msaSupportDir, uniqueOptimalTours)
	if err != nil {
		return err
	}

	return saveOptimalToursStatistics(optimalUniqueToursCsvPath, toursStatistics)
}

func BuildToursMatrix(uniqueOptimalTours map[string][]int, dimension int) [][]float64 {
	toursMatrix := make([][]float64, dimension)
	for i := 0; i < dimension; i++ {
		toursMatrix[i] = make([]float64, dimension)
	}

	tourIds := sortedTourIds(uniqueOptimalTours)
	for _, tourId := range tourIds {
		tour := uniqueOptimalTours[tourId]
		n := len(tour)
		for i := 0; i < n-1; i++ {
			start, end := tour[i], tour[i+1]
			toursMatrix[start][end]++
		}

		last, first := tour[n-1], tour[0]
		toursMatrix[last][first]++
	}

	return toursMatrix
}

func BuildMsaSupportToursOverlapMatrix(msaSupport, toursMatrix [][]float64, toursCount int) [][]float64 {
	dimension := len(msaSupport)
	overlapMatrix := make([][]float64, dimension)

	maxMsaSupportSelections := float64(dimension - 1)
	maxTourSelections := float64(toursCount)

	for i := 0; i < dimension; i++ {
		overlapMatrix[i] = make([]float64, dimension)
		for j := 0; j < dimension; j++ {
			msaSupportFrequency := 0.0
			if maxMsaSupportSelections > 0 {
				msaSupportFrequency = msaSupport[i][j] / maxMsaSupportSelections
			}

			toursFrequency := 0.0
			if maxTourSelections > 0 {
				toursFrequency = toursMatrix[i][j] / maxTourSelections
			}

			overlapMatrix[i][j] = msaSupportFrequency * toursFrequency
		}
	}

	return overlapMatrix
}

func saveOptimalToursStatistics(optimalUniqueToursCsvPath string, toursStatistics []TourStatistics) error {
	header := []string{
		"Tour",
		"Commonality with MSA support",
		"Min commonality with MSA",
		"Avg commonality with MSA",
		"Max commonality with MSA",
	}

	file, err := os.Create(optimalUniqueToursCsvPath)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)

	if err := writer.Write(header); err != nil {
		return err
	}

	floatFormat := "%.2f"
	for _, tourStatistics := range toursStatistics {
		record := []string{
			tourStatistics.tourId,
			fmt.Sprintf(floatFormat, tourStatistics.commonalityWithMsaSupport),
			fmt.Sprintf(floatFormat, tourStatistics.minCommonalityWithMsa),
			fmt.Sprintf(floatFormat, tourStatistics.averageCommonalityWithMsa),
			fmt.Sprintf(floatFormat, tourStatistics.maxCommonalityWithMsa),
		}

		if err := writer.Write(record); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func calculateToursStatistics(msaSupportDir string, uniqueOptimalTours map[string][]int) ([]TourStatistics, error) {
	if len(uniqueOptimalTours) == 0 {
		return nil, nil
	}

	msaSupportMatrix, err := msaSupport.Read(msaSupportDir)
	if err != nil {
		return nil, err
	}

	msas, err := msaSupport.ReadMsas(msaSupportDir)
	if err != nil {
		return nil, err
	}

	toursStatistics := make([]TourStatistics, 0, len(uniqueOptimalTours))
	for _, tourId := range sortedTourIds(uniqueOptimalTours) {
		tour := uniqueOptimalTours[tourId]
		tourEdges := models.ConvertTourToEdges(tour)

		commonalityWithMsaSupport := calculateCommonalityWithMatrix(tourEdges, msaSupportMatrix)

		minCommonalityWithMsa := math.MaxFloat64
		averageCommonalityWithMsa := 0.0
		maxCommonalityWithMsa := -math.MaxFloat64

		for _, msa := range msas {
			commonalityWithMsa := calculateCommonalityWithMatrix(tourEdges, msa)

			if commonalityWithMsa < minCommonalityWithMsa {
				minCommonalityWithMsa = commonalityWithMsa
			}

			averageCommonalityWithMsa += commonalityWithMsa

			if commonalityWithMsa > maxCommonalityWithMsa {
				maxCommonalityWithMsa = commonalityWithMsa
			}
		}

		if len(msas) > 0 {
			averageCommonalityWithMsa /= float64(len(msas))
		} else {
			minCommonalityWithMsa = 0
			maxCommonalityWithMsa = 0
		}

		toursStatistics = append(toursStatistics, TourStatistics{
			tourId,
			commonalityWithMsaSupport,
			minCommonalityWithMsa,
			averageCommonalityWithMsa,
			maxCommonalityWithMsa,
		})
	}

	sort.SliceStable(toursStatistics, func(i, j int) bool {
		if toursStatistics[i].commonalityWithMsaSupport != toursStatistics[j].commonalityWithMsaSupport {
			return toursStatistics[i].commonalityWithMsaSupport > toursStatistics[j].commonalityWithMsaSupport
		}

		if toursStatistics[i].averageCommonalityWithMsa != toursStatistics[j].averageCommonalityWithMsa {
			return toursStatistics[i].averageCommonalityWithMsa > toursStatistics[j].averageCommonalityWithMsa
		}

		return toursStatistics[i].tourId > toursStatistics[j].tourId
	})

	return toursStatistics, nil
}

func calculateCommonalityWithMatrix(tourEdges []Edge, matrix [][]float64) float64 {
	commonality := 0.0
	for _, edge := range tourEdges {
		if matrix[edge.From][edge.To] > 0 {
			commonality++
		}
	}

	if len(tourEdges) == 0 {
		return 0
	}

	return 100 * commonality / float64(len(tourEdges))
}

func saveInstanceAnalysis(path, instance string, dimension int, metrics InstanceMetrics) error {
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

	rows := [][]string{
		{"Instance", instance},
		{"Dimension", strconv.Itoa(dimension)},
		{"Unique found-optimal tours", strconv.Itoa(metrics.FoundOptimalTourCount)},
		{"Unique found-optimal edges", strconv.Itoa(metrics.UniqueFoundOptimalEdgeCount)},
		{"Raw MSA support found-optimal edge coverage [%]", formatPercent(metrics.OptimalEdgeMsaSupportCoverage)},
		{"Threshold MSA support threshold", formatFloat(metrics.HighMsaSupportThreshold)},
		{"Threshold MSA support edges", strconv.Itoa(metrics.HighMsaSupportEdgeCount)},
		{"Threshold MSA support precision [%]", formatPercent(metrics.HighMsaSupportPrecision)},
		{"Threshold MSA support recall [%]", formatPercent(metrics.HighMsaSupportRecall)},
		{"Threshold MSA support lift", formatFloat(metrics.HighMsaSupportLift)},
		{"Top MSA support edges", strconv.Itoa(metrics.TopMsaSupportEdgeCount)},
		{"Top MSA support precision [%]", formatPercent(metrics.TopMsaSupportPrecision)},
		{"Top MSA support recall [%]", formatPercent(metrics.TopMsaSupportRecall)},
		{"Top MSA support lift", formatFloat(metrics.TopMsaSupportLift)},
		{"Top MSA support avg tour coverage [%]", formatPercent(metrics.TopMsaSupportAverageTourCoverage)},
		{"Top MSA support best tour coverage [%]", formatPercent(metrics.TopMsaSupportMaxTourCoverage)},
	}

	for _, row := range rows {
		if err := writer.Write(row); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
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
		"High MSA support edges",
		"High MSA support density [%]",
		"Precision vs found-optimal edges [%]",
		"Recall vs found-optimal edges [%]",
		"Lift vs found-optimal edge density",
		"High-MSA support edges not seen in found-optimal tours",
	}
	if err := writer.Write(header); err != nil {
		return err
	}

	for _, metric := range metrics {
		row := []string{
			formatFloat(metric.Threshold),
			strconv.Itoa(metric.HighMsaSupportEdgeCount),
			formatPercent(metric.HighMsaSupportDensity),
			formatPercent(metric.Precision),
			formatPercent(metric.Recall),
			formatFloat(metric.Lift),
			strconv.Itoa(metric.FalsePositiveHighEdges),
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
		"Unique found-optimal tours",
		"Unique found-optimal edges",
		"Raw MSA support found-optimal edge coverage [%]",
		"Threshold MSA support threshold",
		"Threshold MSA support edges",
		"Threshold MSA support precision [%]",
		"Threshold MSA support recall [%]",
		"Threshold MSA support lift",
		"Top MSA support edges",
		"Top MSA support precision [%]",
		"Top MSA support recall [%]",
		"Top MSA support lift",
		"Top MSA support avg tour coverage [%]",
		"Top MSA support best tour coverage [%]",
	}
	if err := writer.Write(header); err != nil {
		return err
	}

	for _, analysis := range analyses {
		metrics := analysis.Metrics
		row := []string{
			analysis.Instance,
			strconv.Itoa(metrics.FoundOptimalTourCount),
			strconv.Itoa(metrics.UniqueFoundOptimalEdgeCount),
			formatPercent(metrics.OptimalEdgeMsaSupportCoverage),
			formatFloat(metrics.HighMsaSupportThreshold),
			strconv.Itoa(metrics.HighMsaSupportEdgeCount),
			formatPercent(metrics.HighMsaSupportPrecision),
			formatPercent(metrics.HighMsaSupportRecall),
			formatFloat(metrics.HighMsaSupportLift),
			strconv.Itoa(metrics.TopMsaSupportEdgeCount),
			formatPercent(metrics.TopMsaSupportPrecision),
			formatPercent(metrics.TopMsaSupportRecall),
			formatFloat(metrics.TopMsaSupportLift),
			formatPercent(metrics.TopMsaSupportAverageTourCoverage),
			formatPercent(metrics.TopMsaSupportMaxTourCoverage),
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

	instancesWithTours := make([]InstanceAnalysis, 0, len(analyses))
	for _, analysis := range analyses {
		if analysis.Metrics.FoundOptimalTourCount == 0 {
			continue
		}
		instancesWithTours = append(instancesWithTours, analysis)
	}

	count := float64(len(instancesWithTours))
	average := func(sum float64) float64 {
		if count == 0 {
			return 0
		}
		return sum / count
	}

	rawCoverageSum := 0.0
	thresholdEdgeSum := 0.0
	thresholdPrecisionSum := 0.0
	thresholdRecallSum := 0.0
	thresholdLiftSum := 0.0
	topEdgeSum := 0.0
	topPrecisionSum := 0.0
	topRecallSum := 0.0
	topLiftSum := 0.0
	topAverageTourCoverageSum := 0.0
	topMaxTourCoverageSum := 0.0
	for _, analysis := range instancesWithTours {
		metrics := analysis.Metrics
		rawCoverageSum += metrics.OptimalEdgeMsaSupportCoverage
		thresholdEdgeSum += float64(metrics.HighMsaSupportEdgeCount)
		thresholdPrecisionSum += metrics.HighMsaSupportPrecision
		thresholdRecallSum += metrics.HighMsaSupportRecall
		thresholdLiftSum += metrics.HighMsaSupportLift
		topEdgeSum += float64(metrics.TopMsaSupportEdgeCount)
		topPrecisionSum += metrics.TopMsaSupportPrecision
		topRecallSum += metrics.TopMsaSupportRecall
		topLiftSum += metrics.TopMsaSupportLift
		topAverageTourCoverageSum += metrics.TopMsaSupportAverageTourCoverage
		topMaxTourCoverageSum += metrics.TopMsaSupportMaxTourCoverage
	}

	var builder strings.Builder
	builder.WriteString("# MSA Support / Solution Analysis\n\n")
	builder.WriteString(fmt.Sprintf("- Instances analyzed: %d\n", len(analyses)))
	builder.WriteString(fmt.Sprintf("- Instances with found optimal tours: %d\n", len(instancesWithTours)))
	builder.WriteString(fmt.Sprintf("- Instances without found optimal tours: %d\n\n", len(analyses)-len(instancesWithTours)))

	builder.WriteString("## How To Read This\n\n")
	builder.WriteString("- `Precision` is the share of selected MSA support edges that were seen in found optimal tours.\n")
	builder.WriteString("- `Recall` is the share of the found-optimal edge union selected by the MSA support variant.\n")
	builder.WriteString("- `Lift` is precision divided by random edge-hit probability for that instance; higher means the selected edges are less random.\n")
	builder.WriteString("- `Tour coverage` is per-tour, not union-based: it measures how much of one complete found optimal tour is covered by the selected edges.\n\n")

	builder.WriteString("## Raw MSA Support Coverage\n\n")
	builder.WriteString("This checks whether MSA support contains found-optimal edges at all before selecting or filtering them.\n\n")
	builder.WriteString("| Metric | Average |\n")
	builder.WriteString("|---|---:|\n")
	builder.WriteString(fmt.Sprintf("| Found-optimal edges present in MSA support | %s%% |\n", formatPercent(average(rawCoverageSum))))
	builder.WriteString("\n")

	builder.WriteString("| Instance | Coverage % |\n")
	builder.WriteString("|---|---:|\n")
	for _, analysis := range instancesWithTours {
		metrics := analysis.Metrics
		builder.WriteString(fmt.Sprintf("| %s | %s |\n",
			analysis.Instance,
			formatPercent(metrics.OptimalEdgeMsaSupportCoverage),
		))
	}
	builder.WriteString("\n")

	builder.WriteString("## Threshold MSA support\n\n")
	builder.WriteString(fmt.Sprintf("This variant selects every edge with normalized MSA support at least %.1f.\n\n", highThreshold))
	builder.WriteString("| Metric | Average |\n")
	builder.WriteString("|---|---:|\n")
	builder.WriteString(fmt.Sprintf("| Selected edges | %.1f |\n", average(thresholdEdgeSum)))
	builder.WriteString(fmt.Sprintf("| Precision | %s%% |\n", formatPercent(average(thresholdPrecisionSum))))
	builder.WriteString(fmt.Sprintf("| Recall | %s%% |\n", formatPercent(average(thresholdRecallSum))))
	builder.WriteString(fmt.Sprintf("| Lift | %s |\n\n", formatFloat(average(thresholdLiftSum))))

	builder.WriteString("| Instance | Edges | Precision % | Recall % | Lift |\n")
	builder.WriteString("|---|---:|---:|---:|---:|\n")
	for _, analysis := range instancesWithTours {
		metrics := analysis.Metrics
		builder.WriteString(fmt.Sprintf("| %s | %d | %s | %s | %s |\n",
			analysis.Instance,
			metrics.HighMsaSupportEdgeCount,
			formatPercent(metrics.HighMsaSupportPrecision),
			formatPercent(metrics.HighMsaSupportRecall),
			formatFloat(metrics.HighMsaSupportLift),
		))
	}
	builder.WriteString("\n")

	builder.WriteString("## Top N-1 MSA support\n\n")
	builder.WriteString("This variant selects a fixed number of edges: the top `dimension - 1` positive MSA support edges, with deterministic tie-breaking by edge id.\n\n")
	builder.WriteString("| Metric | Average |\n")
	builder.WriteString("|---|---:|\n")
	builder.WriteString(fmt.Sprintf("| Selected edges | %.1f |\n", average(topEdgeSum)))
	builder.WriteString(fmt.Sprintf("| Precision | %s%% |\n", formatPercent(average(topPrecisionSum))))
	builder.WriteString(fmt.Sprintf("| Recall | %s%% |\n", formatPercent(average(topRecallSum))))
	builder.WriteString(fmt.Sprintf("| Lift | %s |\n", formatFloat(average(topLiftSum))))
	builder.WriteString(fmt.Sprintf("| Average tour coverage | %s%% |\n", formatPercent(average(topAverageTourCoverageSum))))
	builder.WriteString(fmt.Sprintf("| Best-tour coverage | %s%% |\n\n", formatPercent(average(topMaxTourCoverageSum))))

	builder.WriteString("| Instance | Edges | Precision % | Recall % | Avg tour coverage % | Best-tour coverage % |\n")
	builder.WriteString("|---|---:|---:|---:|---:|---:|\n")
	for _, analysis := range instancesWithTours {
		metrics := analysis.Metrics
		builder.WriteString(fmt.Sprintf("| %s | %d | %s | %s | %s | %s |\n",
			analysis.Instance,
			metrics.TopMsaSupportEdgeCount,
			formatPercent(metrics.TopMsaSupportPrecision),
			formatPercent(metrics.TopMsaSupportRecall),
			formatPercent(metrics.TopMsaSupportAverageTourCoverage),
			formatPercent(metrics.TopMsaSupportMaxTourCoverage),
		))
	}
	builder.WriteString("\n")

	if len(analyses) > len(instancesWithTours) {
		builder.WriteString("Instances without found optimal tours are omitted from the overlap tables because precision, recall, and tour coverage cannot be interpreted without a found-tour edge set.\n")
	}

	return os.WriteFile(path, []byte(builder.String()), 0644)
}

func buildTourEdgeSet(uniqueOptimalTours map[string][]int) map[Edge]struct{} {
	edges := make(map[Edge]struct{})
	for _, tourId := range sortedTourIds(uniqueOptimalTours) {
		for _, edge := range models.ConvertTourToEdges(uniqueOptimalTours[tourId]) {
			edges[edge] = struct{}{}
		}
	}
	return edges
}

func buildTopMsaSupportEdgeSet(msaSupport [][]float64, limit int) map[Edge]struct{} {
	if limit <= 0 {
		return map[Edge]struct{}{}
	}

	type msaSupportEdgeCandidate struct {
		edge   Edge
		signal float64
	}

	dimension := len(msaSupport)
	candidates := make([]msaSupportEdgeCandidate, 0, dimension*(dimension-1))
	for i := 0; i < dimension; i++ {
		for j := 0; j < len(msaSupport[i]); j++ {
			if i == j || msaSupport[i][j] <= 0 {
				continue
			}

			candidates = append(candidates, msaSupportEdgeCandidate{
				edge:   Edge{From: i, To: j},
				signal: normalizedMsaSupportValue(msaSupport[i][j], dimension),
			})
		}
	}

	sort.SliceStable(candidates, func(i, j int) bool {
		left := candidates[i]
		right := candidates[j]
		if left.signal != right.signal {
			return left.signal > right.signal
		}
		if left.edge.From != right.edge.From {
			return left.edge.From < right.edge.From
		}
		return left.edge.To < right.edge.To
	})

	if limit > len(candidates) {
		limit = len(candidates)
	}

	edges := make(map[Edge]struct{}, limit)
	for i := 0; i < limit; i++ {
		edges[candidates[i].edge] = struct{}{}
	}
	return edges
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

func calculateTourCoverages(uniqueOptimalTours map[string][]int, edgeSet map[Edge]struct{}) []float64 {
	coverages := make([]float64, 0, len(uniqueOptimalTours))
	for _, tourId := range sortedTourIds(uniqueOptimalTours) {
		tourEdges := models.ConvertTourToEdges(uniqueOptimalTours[tourId])
		if len(tourEdges) == 0 {
			coverages = append(coverages, 0)
			continue
		}

		matches := 0
		for _, edge := range tourEdges {
			if _, ok := edgeSet[edge]; ok {
				matches++
			}
		}
		coverages = append(coverages, float64(matches)/float64(len(tourEdges)))
	}
	return coverages
}

func validateSquareMatrix(matrix [][]float64) error {
	dimension := len(matrix)
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

func removeStaleTourPlots(config InstanceConfig) error {
	for _, path := range []string{
		config.ToursHeatmapPath,
		config.ToursHistogramPath,
		config.MsaSupportToursOverlapHeatmapPath,
	} {
		if err := removeIfExists(path); err != nil {
			return err
		}
	}

	return nil
}

func removeIfExists(path string) error {
	if path == "" {
		return nil
	}

	if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
		return err
	}

	return nil
}

func normalizeTour(tour []int) []int {
	if len(tour) == 0 {
		return tour
	}

	minIndex := 0
	for i, v := range tour {
		if v < tour[minIndex] {
			minIndex = i
		}
	}

	normalized := make([]int, len(tour))
	for i := range tour {
		normalized[i] = tour[(minIndex+i)%len(tour)]
	}

	return normalized
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

func average(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}

	sum := 0.0
	for _, value := range values {
		sum += value
	}
	return sum / float64(len(values))
}

func median(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}

	sortedValues := append([]float64(nil), values...)
	sort.Float64s(sortedValues)
	middle := len(sortedValues) / 2
	if len(sortedValues)%2 == 1 {
		return sortedValues[middle]
	}

	return (sortedValues[middle-1] + sortedValues[middle]) / 2
}

func minFloat(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}

	result := values[0]
	for _, value := range values[1:] {
		if value < result {
			result = value
		}
	}
	return result
}

func maxFloat(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}

	result := values[0]
	for _, value := range values[1:] {
		if value > result {
			result = value
		}
	}
	return result
}

func countUniqueValues(data []float64) int {
	unique := make(map[float64]struct{})
	for _, value := range data {
		unique[value] = struct{}{}
	}

	return len(unique)
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

func flattenMatrix(matrix [][]float64) []float64 {
	var result []float64
	for _, row := range matrix {
		result = append(result, row...)
	}
	return result
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
