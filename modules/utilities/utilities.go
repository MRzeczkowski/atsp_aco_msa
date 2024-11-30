package utilities

import (
	"fmt"
	"regexp"
	"strconv"
)

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
	sum := 0.0
	p := len(tour)

	for i := 0; i < p-1; i++ {
		start, end := tour[i], tour[i+1]
		sum += distances[start][end]
	}

	if p > 0 {
		last, first := tour[p-1], tour[0]
		sum += distances[last][first]
	}

	return sum
}

func IndexOf(element int, data []int) int {
	for k, v := range data {
		if element == v {
			return k
		}
	}

	return -1 //not found.
}
