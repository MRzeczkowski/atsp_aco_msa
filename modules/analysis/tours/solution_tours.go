package tours

import (
	"atsp_aco_msa/modules/algorithms/cyclecover"
	"atsp_aco_msa/modules/algorithms/msaHeuristic"
	"atsp_aco_msa/modules/models"
	"atsp_aco_msa/modules/utilities"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"sort"
)

type Edge = models.Edge

type InstanceConfig struct {
	Name                                string
	Dimension                           int
	MsaHeuristicDirectoryPath           string
	CycleCoverDirectoryPath             string
	OptimalToursCsvPath                 string
	ToursHeatmapPath                    string
	ToursHistogramPath                  string
	MsaHeuristicToursOverlapHeatmapPath string
	CycleCoverToursOverlapHeatmapPath   string
}

type Config struct {
	Instances []InstanceConfig
}

type TourStatistics struct {
	tourId                                                                                               string
	commonalityWithMsaHeuristic, minCommonalityWithMsa, averageCommonalityWithMsa, maxCommonalityWithMsa float64
}

func AnalyzeInstances(config Config) error {
	instances := append([]InstanceConfig(nil), config.Instances...)
	sort.SliceStable(instances, func(i, j int) bool {
		return instances[i].Name < instances[j].Name
	})

	for _, instance := range instances {
		if err := AnalyzeInstance(instance); err != nil {
			return err
		}
	}

	return nil
}

func AnalyzeInstance(config InstanceConfig) error {
	uniqueOptimalTours, err := ReadOptimal(config.OptimalToursCsvPath)
	if err != nil {
		return fmt.Errorf("%s: failed to read optimal tours: %w", config.Name, err)
	}

	msaHeuristicMatrix, err := msaHeuristic.Read(config.MsaHeuristicDirectoryPath)
	if err != nil {
		return fmt.Errorf("%s: failed to read MSA heuristic: %w", config.Name, err)
	}
	if err := validateSquareMatrix(msaHeuristicMatrix); err != nil {
		return fmt.Errorf("%s: invalid MSA heuristic: %w", config.Name, err)
	}
	if config.Dimension != 0 && config.Dimension != len(msaHeuristicMatrix) {
		return fmt.Errorf("%s: configured dimension %d does not match MSA heuristic dimension %d", config.Name, config.Dimension, len(msaHeuristicMatrix))
	}
	if err := validateTours(uniqueOptimalTours, len(msaHeuristicMatrix)); err != nil {
		return fmt.Errorf("%s: invalid found optimal tours: %w", config.Name, err)
	}

	toursCount := len(uniqueOptimalTours)
	if toursCount > 0 {
		cycleCoverMatrix, err := cyclecover.Read(config.CycleCoverDirectoryPath, len(msaHeuristicMatrix))
		if err != nil {
			return fmt.Errorf("%s: failed to read cycle cover: %w", config.Name, err)
		}

		toursMatrix := BuildMatrix(uniqueOptimalTours, config.Dimension)
		toursHeatmapPlotTitle := config.Name + " tours heatmap"
		if err := utilities.SaveHeatmapFromMatrix(toursMatrix, toursHeatmapPlotTitle, config.ToursHeatmapPath); err != nil {
			return err
		}

		overlapMatrix := BuildMsaOverlap(msaHeuristicMatrix, toursMatrix, toursCount)
		overlapHeatmapPlotTitle := config.Name + " MSA heuristic/tours overlap heatmap"
		if err := utilities.SaveHeatmapFromMatrix(overlapMatrix, overlapHeatmapPlotTitle, config.MsaHeuristicToursOverlapHeatmapPath); err != nil {
			return err
		}

		cycleCoverOverlapMatrix := BuildCycleCoverOverlap(cycleCoverMatrix, toursMatrix, toursCount)
		cycleCoverOverlapHeatmapPlotTitle := config.Name + " cycle cover/tours overlap heatmap"
		if err := utilities.SaveHeatmapFromMatrix(cycleCoverOverlapMatrix, cycleCoverOverlapHeatmapPlotTitle, config.CycleCoverToursOverlapHeatmapPath); err != nil {
			return err
		}

		if toursCount > 1 {
			dataForHistogram := filterZeroes(flattenMatrix(toursMatrix))
			if len(dataForHistogram) > 0 {
				toursHistogramPlotTitle := config.Name + " tours histogram"
				numberOfBins := countUniqueValues(dataForHistogram)

				if err := utilities.SaveHistogramFromData(dataForHistogram, numberOfBins, toursHistogramPlotTitle, config.ToursHistogramPath); err != nil {
					return err
				}
			}
		} else if err := removeIfExists(config.ToursHistogramPath); err != nil {
			return err
		}
	} else if err := removeStaleTourPlots(config); err != nil {
		return err
	}

	return nil
}

func ReadOptimal(optimalUniqueToursCsvPath string) (map[string][]int, error) {
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

func AddUnique(savedTours map[string][]int, tour []int) {
	normalizedTour := normalizeTour(tour)
	keyJson, _ := json.Marshal(normalizedTour)
	key := string(keyJson)

	if _, exists := savedTours[key]; exists {
		return
	}

	savedTours[key] = tour
}

func SaveStatistics(optimalUniqueToursCsvPath, msaHeuristicDir string, uniqueOptimalTours map[string][]int) error {
	toursStatistics, err := calculateToursStatistics(msaHeuristicDir, uniqueOptimalTours)
	if err != nil {
		return err
	}

	return saveOptimalToursStatistics(optimalUniqueToursCsvPath, toursStatistics)
}

func BuildMatrix(uniqueOptimalTours map[string][]int, dimension int) [][]float64 {
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

func BuildMsaOverlap(msaHeuristic, toursMatrix [][]float64, toursCount int) [][]float64 {
	dimension := len(msaHeuristic)
	overlapMatrix := make([][]float64, dimension)

	maxMsaHeuristicSelections := float64(dimension - 1)
	maxTourSelections := float64(toursCount)

	for i := 0; i < dimension; i++ {
		overlapMatrix[i] = make([]float64, dimension)
		for j := 0; j < dimension; j++ {
			msaHeuristicFrequency := 0.0
			if maxMsaHeuristicSelections > 0 {
				msaHeuristicFrequency = msaHeuristic[i][j] / maxMsaHeuristicSelections
			}

			toursFrequency := 0.0
			if maxTourSelections > 0 {
				toursFrequency = toursMatrix[i][j] / maxTourSelections
			}

			overlapMatrix[i][j] = msaHeuristicFrequency * toursFrequency
		}
	}

	return overlapMatrix
}

func BuildCycleCoverOverlap(cycleCoverMatrix, toursMatrix [][]float64, toursCount int) [][]float64 {
	dimension := len(cycleCoverMatrix)
	overlapMatrix := make([][]float64, dimension)
	maxTourSelections := float64(toursCount)

	for i := 0; i < dimension; i++ {
		overlapMatrix[i] = make([]float64, dimension)
		for j := 0; j < dimension; j++ {
			toursFrequency := 0.0
			if maxTourSelections > 0 {
				toursFrequency = toursMatrix[i][j] / maxTourSelections
			}

			overlapMatrix[i][j] = cycleCoverMatrix[i][j] * toursFrequency
		}
	}

	return overlapMatrix
}

func saveOptimalToursStatistics(optimalUniqueToursCsvPath string, toursStatistics []TourStatistics) error {
	header := []string{
		"Tour",
		"Commonality with MSA heuristic",
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
			fmt.Sprintf(floatFormat, tourStatistics.commonalityWithMsaHeuristic),
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

func calculateToursStatistics(msaHeuristicDir string, uniqueOptimalTours map[string][]int) ([]TourStatistics, error) {
	if len(uniqueOptimalTours) == 0 {
		return nil, nil
	}

	msaHeuristicMatrix, err := msaHeuristic.Read(msaHeuristicDir)
	if err != nil {
		return nil, err
	}

	msas, err := msaHeuristic.ReadMsas(msaHeuristicDir)
	if err != nil {
		return nil, err
	}

	toursStatistics := make([]TourStatistics, 0, len(uniqueOptimalTours))
	for _, tourId := range sortedTourIds(uniqueOptimalTours) {
		tour := uniqueOptimalTours[tourId]
		tourEdges := models.ConvertTourToEdges(tour)

		commonalityWithMsaHeuristic := calculateCommonalityWithMatrix(tourEdges, msaHeuristicMatrix)

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
			commonalityWithMsaHeuristic,
			minCommonalityWithMsa,
			averageCommonalityWithMsa,
			maxCommonalityWithMsa,
		})
	}

	sort.SliceStable(toursStatistics, func(i, j int) bool {
		if toursStatistics[i].commonalityWithMsaHeuristic != toursStatistics[j].commonalityWithMsaHeuristic {
			return toursStatistics[i].commonalityWithMsaHeuristic > toursStatistics[j].commonalityWithMsaHeuristic
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
		config.MsaHeuristicToursOverlapHeatmapPath,
		config.CycleCoverToursOverlapHeatmapPath,
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
