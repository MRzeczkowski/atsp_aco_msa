package utilities

import (
	"fmt"
	"image/color"
	"math"
	"os"
	"path"
	"sort"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/palette/moreland"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
	"gonum.org/v1/plot/vg/vgimg"
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

type LinePlotData struct {
	Name   string
	Color  color.RGBA
	Values []float64
}

var (
	plotWidth  = vg.Length(800)
	plotHeight = vg.Length(800)
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

	return savePlotWithPadding(heatmapPlot, plotWidth, plotHeight, plotPath)
}

func SaveHistogramFromData(data []float64, bins int, plotTitle, plotPath string) error {
	values := make(plotter.Values, len(data))
	copy(values, data)

	histogramPlot := plot.New()
	histogramPlot.Title.Text = plotTitle
	histogramPlot.X.Label.Text = "Number of selections"
	histogramPlot.X.Min = 1.0
	histogramPlot.Y.Label.Text = "Count of edges"

	hist, err := plotter.NewHist(values, bins)
	if err != nil {
		return err
	}

	actualBins := len(hist.Bins)
	if actualBins == 0 {
		histogramPlot.Add(hist)
		return savePlotWithPadding(histogramPlot, plotWidth, plotHeight, plotPath)
	}

	firstBin := hist.Bins[0]
	middleBin := hist.Bins[actualBins/2]
	lastBin := hist.Bins[actualBins-1]

	histogramPlot.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for _, bin := range hist.Bins {
			center := (bin.Min + bin.Max) / 2

			label := ""
			if bin == firstBin || bin == middleBin || bin == lastBin {
				label = fmt.Sprintf("%d", int(math.Round(center)))
			}

			ticks = append(ticks, plot.Tick{
				Value: center,
				Label: label,
			})
		}

		return ticks
	})

	heightsSet := make(map[int]bool)
	for _, bin := range hist.Bins {
		heightsSet[int(bin.Weight)] = true
	}

	heights := make([]int, 0, len(heightsSet))
	for k := range heightsSet {
		heights = append(heights, k)
	}

	sort.Slice(heights, func(i, j int) bool {
		return heights[i] < heights[j]
	})

	heightsCount := len(heights)
	smallestHeight := heights[0]
	medianHeight := heights[heightsCount/2]
	biggestHeight := heights[heightsCount-1]

	histogramPlot.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick

		for i := 0; i <= biggestHeight; i++ {

			label := ""
			if i == smallestHeight || i == medianHeight || i == biggestHeight {
				label = fmt.Sprintf("%d", i)
			}

			ticks = append(ticks, plot.Tick{
				Value: float64(i),
				Label: label,
			})
		}

		return ticks
	})

	histogramPlot.Add(hist)

	return savePlotWithPadding(histogramPlot, plotWidth, plotHeight, plotPath)
}

func SaveHistogramFromDataWithLabels(data []float64, bins int, plotTitle, xLabel, yLabel, plotPath string) error {
	values := make(plotter.Values, len(data))
	copy(values, data)

	histogramPlot := plot.New()
	histogramPlot.Title.Text = plotTitle
	histogramPlot.X.Label.Text = xLabel
	histogramPlot.Y.Label.Text = yLabel
	if bins < 1 {
		bins = 1
	}

	if len(values) == 0 {
		return savePlotWithPadding(histogramPlot, plotWidth, plotHeight, plotPath)
	}

	minValue, maxValue := minMax(values)
	if minValue == maxValue {
		bars, err := plotter.NewBarChart(plotter.Values{float64(len(values))}, vg.Points(40))
		if err != nil {
			return err
		}
		histogramPlot.Add(bars)
		histogramPlot.NominalX(formatHistogramValue(minValue))
		return savePlotWithPadding(histogramPlot, plotWidth, plotHeight, plotPath)
	}

	hist, err := plotter.NewHist(values, bins)
	if err != nil {
		return err
	}

	histogramPlot.Add(hist)
	return savePlotWithPadding(histogramPlot, plotWidth, plotHeight, plotPath)
}

func minMax(values plotter.Values) (float64, float64) {
	minValue := values[0]
	maxValue := values[0]
	for _, value := range values[1:] {
		if value < minValue {
			minValue = value
		}
		if value > maxValue {
			maxValue = value
		}
	}
	return minValue, maxValue
}

func formatHistogramValue(value float64) string {
	rounded := math.Round(value)
	if math.Abs(value-rounded) < 1e-9 {
		return fmt.Sprintf("%.0f", rounded)
	}
	return fmt.Sprintf("%.2f", value)
}

func SaveLinePlotFromData(linePlotData []LinePlotData, plotTitle, plotPath string) error {
	plot := plot.New()
	plot.Title.Text = plotTitle
	plot.X.Label.Text = "Iteration"
	plot.Y.Label.Text = "Deviation"

	for _, data := range linePlotData {
		addLine(plot, data)
	}

	plot.Legend.Top = true

	return savePlotWithPadding(plot, plotWidth, plotHeight, plotPath)
}

func addLine(p *plot.Plot, data LinePlotData) {
	// Build the XYs data
	pts := make(plotter.XYs, len(data.Values))
	for i, v := range data.Values {
		pts[i].X = float64(i)
		pts[i].Y = v
	}

	// Create the line
	line, err := plotter.NewLine(pts)
	if err != nil {
		panic(err)
	}
	line.LineStyle.Width = vg.Points(2)
	line.LineStyle.Color = data.Color

	// Add the line and a legend entry
	p.Add(line)
	p.Legend.Add(data.Name, line)
}

func savePlotWithPadding(p *plot.Plot, width, height vg.Length, filePath string) error {
	// Increase the outer canvas size a bit more
	canvas := vgimg.New(width+4*vg.Centimeter, height+4*vg.Centimeter)
	dc := draw.New(canvas)

	// Increase the padding inside the canvas
	pad := vg.Centimeter / 5
	dc = draw.Crop(dc, 0, -pad, 0, 0)

	p.Draw(dc)

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		dir := path.Dir(filePath)
		os.MkdirAll(dir, 0700)
	}

	f, err := os.Create(filePath)
	if err != nil {
		return err
	}
	defer f.Close()

	png := vgimg.PngCanvas{Canvas: canvas}
	if _, err := png.WriteTo(f); err != nil {
		return err
	}

	return nil
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
