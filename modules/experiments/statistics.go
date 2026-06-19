package experiments

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
)

const (
	heuristicCycleCoverMsaPatching = "cycle-cover-msa-patching"
	heuristicRandomSparse          = "random-sparse"
	heuristicShuffledMsa           = "shuffled-msa"
	heuristicStrictShuffledMsa     = "strict-shuffled-msa"
	heuristicRootedShuffledMsa     = "rooted-shuffled-msa"
)

var StatisticsCsvHeader = []string{
	"Alpha",
	"Beta",
	"Rho",
	"Heuristic weight",
	"Iterations",
	"Min best at iteration",
	"Avg best at iteration",
	"Max best at iteration",
	"Min local search improvements",
	"Avg local search improvements",
	"Max local search improvements",
	"Min best deviation",
	"Avg best deviation",
	"Max best deviation",
	"Success rate [%]",
}

func StatisticsCsvHeaderForHeuristic(heuristic string) []string {
	if heuristic == heuristicCycleCoverMsaPatching {
		header := append([]string{}, StatisticsCsvHeader[:4]...)
		header = append(header, "MSA patch bias")
		header = append(header, StatisticsCsvHeader[4:]...)
		return header
	}

	if heuristicUsesRandomSeedColumn(heuristic) {
		header := append([]string{}, StatisticsCsvHeader[:4]...)
		header = append(header, "Random seed")
		header = append(header, StatisticsCsvHeader[4:]...)
		return header
	}

	return StatisticsCsvHeader
}

func heuristicUsesRandomSeedColumn(heuristic string) bool {
	return heuristic == heuristicRandomSparse ||
		heuristic == heuristicShuffledMsa ||
		heuristic == heuristicStrictShuffledMsa ||
		heuristic == heuristicRootedShuffledMsa
}

func ReadStatistics(csvFilePath string) ([]ExperimentsDataStatistics, error) {
	file, err := os.Open(csvFilePath)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)

	// Read the header and skip it
	header, err := reader.Read() // Read the first line (header)
	if err != nil {
		return nil, fmt.Errorf("failed to read header: %w", err)
	}
	includeMsaPatchBias := slices.Equal(header, StatisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching))
	includeRandomSeed := slices.Equal(header, StatisticsCsvHeaderForHeuristic(heuristicRandomSparse))
	if !slices.Equal(header, StatisticsCsvHeader) && !includeMsaPatchBias && !includeRandomSeed {
		return nil, fmt.Errorf("invalid statistics header")
	}
	defaultMsaPatchBias := defaultMsaPatchBiasForStatisticsPath(csvFilePath)

	var statistics []ExperimentsDataStatistics

	// Parse the remaining rows
	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("failed to read records: %w", err)
	}

	for _, record := range records {
		expectedRecordLength := len(StatisticsCsvHeader)
		if includeMsaPatchBias {
			expectedRecordLength = len(StatisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching))
		} else if includeRandomSeed {
			expectedRecordLength = len(StatisticsCsvHeaderForHeuristic(heuristicRandomSparse))
		}
		if len(record) != expectedRecordLength {
			return nil, fmt.Errorf("invalid record length: %d", len(record))
		}

		statistic, err := ParseStatisticsRecord(record, includeMsaPatchBias, includeRandomSeed, defaultMsaPatchBias)
		if err != nil {
			return nil, err
		}

		statistics = append(statistics, statistic)
	}

	return statistics, nil
}

func defaultMsaPatchBiasForStatisticsPath(csvFilePath string) float64 {
	if strings.Contains(filepath.Base(csvFilePath), "cycle_cover_msa_patching") {
		return 1.0
	}

	return 0.0
}

func ParseStatisticsRecord(record []string, includeMsaPatchBias, includeRandomSeed bool, defaultMsaPatchBias float64) (ExperimentsDataStatistics, error) {
	expectedRecordLength := len(StatisticsCsvHeader)
	if includeMsaPatchBias {
		expectedRecordLength = len(StatisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching))
	} else if includeRandomSeed {
		expectedRecordLength = len(StatisticsCsvHeaderForHeuristic(heuristicRandomSparse))
	}
	if len(record) != expectedRecordLength {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid record length: %d", len(record))
	}

	alpha, err := strconv.ParseFloat(record[0], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid Alpha: %w", err)
	}
	beta, err := strconv.ParseFloat(record[1], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid Beta: %w", err)
	}
	rho, err := strconv.ParseFloat(record[2], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid Rho: %w", err)
	}
	heuristicWeight, err := strconv.ParseFloat(record[3], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid heuristic weight: %w", err)
	}
	msaPatchBias := defaultMsaPatchBias
	indexOffset := 0
	if includeMsaPatchBias {
		var err error
		msaPatchBias, err = strconv.ParseFloat(record[4], 64)
		if err != nil {
			return ExperimentsDataStatistics{}, fmt.Errorf("invalid MSA patch bias: %w", err)
		}
		indexOffset = 1
	}
	randomSeed := int64(0)
	if includeRandomSeed {
		parsedRandomSeed, err := strconv.ParseInt(record[4+indexOffset], 10, 64)
		if err != nil {
			return ExperimentsDataStatistics{}, fmt.Errorf("invalid random seed: %w", err)
		}
		randomSeed = parsedRandomSeed
		indexOffset++
	}
	iterations, err := strconv.Atoi(record[4+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid Iterations: %w", err)
	}
	minBestAtIteration, err := strconv.Atoi(record[5+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min best at iteration: %w", err)
	}
	averageBestAtIteration, err := strconv.ParseFloat(record[6+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg best at iteration: %w", err)
	}
	maxBestAtIteration, err := strconv.Atoi(record[7+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max best at iteration: %w", err)
	}
	minThreeOptImprovementsCount, err := strconv.Atoi(record[8+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min local search improvements: %w", err)
	}
	averageThreeOptImprovementsCount, err := strconv.ParseFloat(record[9+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg local search improvements: %w", err)
	}
	maxThreeOptImprovementsCount, err := strconv.Atoi(record[10+indexOffset])
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max local search improvements: %w", err)
	}
	minBestDeviation, err := strconv.ParseFloat(record[11+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid min best deviation: %w", err)
	}
	averageBestDeviation, err := strconv.ParseFloat(record[12+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid avg best deviation: %w", err)
	}
	maxBestDeviation, err := strconv.ParseFloat(record[13+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid max best deviation: %w", err)
	}
	successRate, err := strconv.ParseFloat(record[14+indexOffset], 64)
	if err != nil {
		return ExperimentsDataStatistics{}, fmt.Errorf("invalid success rate: %w", err)
	}

	return ExperimentsDataStatistics{
		ExperimentParameters: ExperimentParameters{
			Alpha:           alpha,
			Beta:            beta,
			Rho:             rho,
			HeuristicWeight: heuristicWeight,
			MsaPatchBias:    msaPatchBias,
			RandomSeed:      randomSeed,
			Iterations:      iterations,
		},
		MinBestAtIteration:               minBestAtIteration,
		AverageBestAtIteration:           averageBestAtIteration,
		MaxBestAtIteration:               maxBestAtIteration,
		MinThreeOptImprovementsCount:     minThreeOptImprovementsCount,
		AverageThreeOptImprovementsCount: averageThreeOptImprovementsCount,
		MaxThreeOptImprovementsCount:     maxThreeOptImprovementsCount,
		MinBestDeviation:                 minBestDeviation,
		AverageBestDeviation:             averageBestDeviation,
		MaxBestDeviation:                 maxBestDeviation,
		SuccessRate:                      successRate,
	}, nil
}

func SaveStatistics(resultCsvPath, heuristic string, statistics []ExperimentsDataStatistics) {
	if err := os.MkdirAll(filepath.Dir(resultCsvPath), 0700); err != nil {
		fmt.Printf("Failed to create statistics directory: %v\n", err)
		return
	}

	file, err := os.Create(resultCsvPath)
	if err != nil {
		fmt.Printf("Failed to save Statistics: %v\n", err)
		return
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	includeMsaPatchBias := heuristic == heuristicCycleCoverMsaPatching
	includeRandomSeed := heuristicUsesRandomSeedColumn(heuristic)

	_ = writer.Write(StatisticsCsvHeaderForHeuristic(heuristic))

	for _, statistic := range statistics {
		writer.Write(StatisticsCsvRecord(statistic, includeMsaPatchBias, includeRandomSeed))
	}

	writer.Flush()
}

func SaveHeuristicStatistics(resultCsvPath string, statistics []HeuristicExperimentStatistics) error {
	if err := os.MkdirAll(filepath.Dir(resultCsvPath), 0700); err != nil {
		return err
	}

	file, err := os.Create(resultCsvPath)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	header := append([]string{"Heuristic"}, StatisticsCsvHeader...)
	if err := writer.Write(header); err != nil {
		return err
	}

	sortedStatistics := append([]HeuristicExperimentStatistics(nil), statistics...)
	sort.SliceStable(sortedStatistics, func(i, j int) bool {
		left, right := sortedStatistics[i].Statistics, sortedStatistics[j].Statistics
		if left.AverageBestDeviation != right.AverageBestDeviation {
			return left.AverageBestDeviation < right.AverageBestDeviation
		}
		if left.SuccessRate != right.SuccessRate {
			return left.SuccessRate > right.SuccessRate
		}
		return sortedStatistics[i].Heuristic < sortedStatistics[j].Heuristic
	})

	for _, statistic := range sortedStatistics {
		record := append([]string{statistic.Heuristic}, StatisticsCsvRecord(statistic.Statistics, false, false)...)
		if err := writer.Write(record); err != nil {
			return err
		}
	}

	writer.Flush()
	return writer.Error()
}

func ReadHeuristicStatistics(resultCsvPath string) ([]HeuristicExperimentStatistics, error) {
	file, err := os.Open(resultCsvPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		return nil, err
	}

	expectedHeader := append([]string{"Heuristic"}, StatisticsCsvHeader...)
	extendedHeader := append([]string{"Heuristic"}, StatisticsCsvHeaderForHeuristic(heuristicCycleCoverMsaPatching)...)
	extendedFormat := slices.Equal(header, extendedHeader)
	if !slices.Equal(header, expectedHeader) && !extendedFormat {
		return nil, fmt.Errorf("invalid heuristic statistics header")
	}

	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	statistics := make([]HeuristicExperimentStatistics, 0, len(records))
	for _, record := range records {
		expectedRecordLength := len(expectedHeader)
		if extendedFormat {
			expectedRecordLength = len(extendedHeader)
		}
		if len(record) != expectedRecordLength {
			return nil, fmt.Errorf("invalid heuristic statistics record length: got %d want %d", len(record), expectedRecordLength)
		}

		parsed, err := ParseStatisticsRecord(record[1:], extendedFormat, false, 0.0)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", record[0], err)
		}

		statistics = append(statistics, HeuristicExperimentStatistics{
			Heuristic:  record[0],
			Statistics: parsed,
		})
	}

	return statistics, nil
}
