package utilities

import (
	"fmt"
	"math"
	"os"
	"regexp"
	"runtime/pprof"
	"strconv"
)

func FastPow(base, exp float64) float64 {
	switch exp {
	case 0.25:
		return math.Sqrt(math.Sqrt(base))
	case 0.5:
		return math.Sqrt(base)
	case 0.75:
		return math.Sqrt(base) * math.Sqrt(math.Sqrt(base))
	case 1:
		return base
	case 1.25:
		return base * math.Sqrt(math.Sqrt(base))
	case 1.5:
		return base * math.Sqrt(base)
	case 1.75:
		return base * math.Sqrt(base) * math.Sqrt(math.Sqrt(base))
	case 2:
		return base * base
	case 2.25:
		return base * base * math.Sqrt(math.Sqrt(base))
	case 2.5:
		return base * base * math.Sqrt(base)
	case 2.75:
		return base * base * math.Sqrt(base) * math.Sqrt(math.Sqrt(base))
	case 3:
		return base * base * base
	case 3.25:
		return base * base * base * math.Sqrt(math.Sqrt(base))
	case 3.5:
		return base * base * base * math.Sqrt(base)
	case 3.75:
		return base * base * base * math.Sqrt(base) * math.Sqrt(math.Sqrt(base))
	case 4:
		return base * base * base * base
	case 4.25:
		return base * base * base * base * math.Sqrt(math.Sqrt(base))
	case 4.5:
		return base * base * base * base * math.Sqrt(base)
	case 4.75:
		return base * base * base * base * math.Sqrt(base) * math.Sqrt(math.Sqrt(base))
	case 5:
		return base * base * base * base * base
	default:
		return math.Pow(base, exp)
	}
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

func StartProfiling() {
	f, err := os.Create("aco.prof")
	if err != nil {
		fmt.Println(err)
		return
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()
}
