package nearestNeighbors

import "sort"

// Build the nearest k neighbors list for each city
func BuildNearestNeighborsLists(distances [][]float64, k int) [][]int {
	n := len(distances)
	neighborsLists := make([][]int, n)

	for i := 0; i < n; i++ {

		type nodeDist struct {
			id       int
			distance float64
		}

		cityDistances := make([]nodeDist, n)
		for j := 0; j < n; j++ {
			cityDistances[j] = nodeDist{id: j, distance: distances[i][j]}
		}

		sort.Slice(cityDistances, func(x, y int) bool {
			return cityDistances[x].distance < cityDistances[y].distance
		})

		for j := 0; j < k; j++ {
			neighborsLists[i] = append(neighborsLists[i], cityDistances[j].id)
		}
	}

	return neighborsLists
}
