package utilities

import (
	"fmt"
	"regexp"
	"strconv"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/palette/moreland"
	"gonum.org/v1/plot/plotter"
)

type HeatMapData struct {
	data       [][]float64
	cellWidth  float64
	cellHeight float64
	minX       float64
	minY       float64
}

func (h *HeatMapData) Dims() (c, r int) {
	return len(h.data[0]), len(h.data)
}

func (h *HeatMapData) X(c int) float64 {
	return h.minX + float64(c)*h.cellWidth
}

func (h *HeatMapData) Y(r int) float64 {
	// Reverse Y-axis to make the origin top-left
	return h.minY + float64(len(h.data)-1-r)*h.cellHeight
}

func (h *HeatMapData) Z(c, r int) float64 {
	return h.data[r][c]
}

var (
	plotWidth  font.Length = 800
	plotLength font.Length = 800
)

func SaveHeatmapFromMatrix(matrix [][]float64, plotTitle, plotPath string) error {
	heatmapData := &HeatMapData{
		data:       matrix,
		cellWidth:  1.0,
		cellHeight: 1.0,
		minX:       0.0,
		minY:       0.0,
	}

	heatmapPlot := plot.New()
	heatmapPlot.Title.Text = plotTitle
	heatmapPlot.HideAxes()

	palette := moreland.SmoothBlueRed().Palette(255)
	heatmap := plotter.NewHeatMap(heatmapData, palette)

	heatmapPlot.Add(heatmap)

	return heatmapPlot.Save(plotWidth, plotLength, plotPath)
}

func SaveHistogramFromData(data []float64, bins int, plotTitle, plotPath string) error {
	values := make(plotter.Values, len(data))
	copy(values, data)

	histogramPlot := plot.New()
	histogramPlot.Title.Text = plotTitle
	histogramPlot.X.Label.Text = "Value"
	histogramPlot.X.Min = 1.0
	histogramPlot.Y.Label.Text = "Frequency"

	hist, err := plotter.NewHist(values, bins)
	if err != nil {
		return err
	}

	histogramPlot.Add(hist)

	return histogramPlot.Save(plotWidth, plotLength, plotPath)
}

func ExtractNumber(input string) (int, error) {

	re := regexp.MustCompile(`\d+`)

	match := re.FindString(input)

	if match == "" {
		return 0, fmt.Errorf("no number found in the string")
	}

	number, err := strconv.Atoi(match)
	if err != nil {
		return 0, err
	}

	return number, nil
}

func FilterStrings(strings []string, condition func(string) bool) []string {
	result := []string{}

	for _, str := range strings {
		if condition(str) {
			result = append(result, str)
		}
	}

	return result
}

func GenerateRange(start, end, step float64) []float64 {
	var rangeSlice []float64

	for i := start; i <= end; i += step {
		rangeSlice = append(rangeSlice, i)
	}

	return rangeSlice
}

func TourLength(tour []int, distances [][]float64) float64 {
	n := len(tour)
	sum := 0.0

	for i := 0; i < n-1; i++ {
		start, end := tour[i], tour[i+1]
		sum += distances[start][end]
	}

	if n > 0 {
		last, first := tour[n-1], tour[0]
		sum += distances[last][first]
	}

	return sum
}
