package compositeMsa

import (
	"atsp_aco_msa/modules/algorithms/edmonds"
	"atsp_aco_msa/modules/models"
	"encoding/csv"
	"math"
	"os"
	"strconv"
)

type Edge = models.Edge

func ReadFromCsv(path string) ([][]float64, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)

	rows, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	data := make([][]float64, len(rows))
	for i, row := range rows {
		data[i] = make([]float64, len(row))
		for j, value := range row {
			data[i][j], err = strconv.ParseFloat(value, 64)
			if err != nil {
				return nil, err
			}
		}
	}

	return data, nil
}

func CreateFromData(matrix [][]float64) [][]float64 {
	vertices, edges, weights := models.ConvertToEdges(matrix)

	dimension := len(matrix)

	cmsa := make([][]float64, dimension)
	for i := range dimension {
		cmsa[i] = make([]float64, dimension)
	}

	msas := make([][]Edge, dimension)
	//msaStats := make([]GraphStats, dimension)
	//leafCounts := make([]float64, dimension)
	occurrences := make(map[Edge]float64)
	lengths := make([]float64, dimension)

	for i := 0; i < dimension; i++ {

		msa := edmonds.FindMSA(vertices, edges, i, weights)

		msas[i] = msa
		//msaStats[i] = calculateTreeStats(msa, weights)
		//leafCounts[i] = countLeaves(msa)

		for _, edge := range msa {
			occurrences[edge]++
			lengths[i] += weights[edge]
		}
	}

	// minLeafCount := math.MaxFloat64

	// for _, count := range leafCounts {
	// 	if count < minLeafCount {
	// 		minLeafCount = count
	// 	}
	// }

	// for i, msa := range msas {
	// 	for _, edge := range msa {
	// 		cmsa[edge.From][edge.To] += math.Pow(1/weights[edge], 5)
	// 		cmsa[edge.From][edge.To] += math.Pow(1/lengths[i], 5)
	// 	}
	// }

	for _, msa := range msas {

		// msaStats[i].stdDevWeight > 0.1*msaStats[i].avgWeight || msaStats[i].maxWeight > 1.5*msaStats[i].avgWeight

		// if math.Abs(msaStats[i].skewness) > 0.8 ||
		// 	msaStats[i].kurtosis > 5.0 {
		// 	continue
		// }

		for _, edge := range msa {
			cmsa[edge.From][edge.To] = math.Pow(occurrences[edge], 1)
			//cmsa[edge.From][edge.To] += math.Pow(1/msaStats[i].stdDevWeight, 1)
			//cmsa[edge.From][edge.To] += math.Pow(1/msaStats[i].skewness, 1)
			//cmsa[edge.From][edge.To] = math.Pow(1/msaStats[i].kurtosis, 1)
			//cmsa[edge.From][edge.To] += math.Pow(1/weights[edge], 1)
			//cmsa[edge.From][edge.To] += math.Pow(1/lengths[i], 1)
			//cmsa[edge.From][edge.To] += math.Pow(1/leafCounts[i], 1)

			//cmsa[edge.From][edge.To] /= occurrence[edge]
		}

		//fmt.Println(i, leafCounts[i])
	}

	// cmsaStats := calculateStats(cmsa)

	// fmt.Println("Min:", cmsaStats.minWeight)
	// fmt.Println("Max:", cmsaStats.maxWeight)
	// fmt.Println("Avg:", cmsaStats.avgWeight)
	// fmt.Println("StdDev:", cmsaStats.stdDevWeight)
	// fmt.Println("Skewness:", cmsaStats.skewness)
	// fmt.Println("Kurtosis:", cmsaStats.kurtosis)

	//return

	for i := 0; i < dimension; i++ {

		// sum := 0.0

		// for j := 0; j < dimension; j++ {
		// 	sum += cmsa[i][j]
		// }

		// for j := 0; j < dimension; j++ {
		// 	cmsa[i][j] /= sum
		// }

		//fmt.Println(cmsa[i])
	}

	return cmsa
}

func SaveToCsv(cmsa [][]float64, path string) error {

	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	for _, row := range cmsa {
		strRow := make([]string, len(row))
		for i, value := range row {
			strRow[i] = strconv.FormatFloat(value, 'f', -1, 64)
		}
		err := writer.Write(strRow)
		if err != nil {
			return err
		}
	}

	return nil
}
