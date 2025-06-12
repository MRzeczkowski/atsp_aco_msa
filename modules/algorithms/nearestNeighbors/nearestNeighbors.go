package nearestNeighbors

import "sort"

// Build the nearest k neighbors list for each node
func BuildNearestNeighborsLists(distances [][]float64, k int) [][]int {
	n := len(distances)
	neighborsLists := make([][]int, n)

	for i := 0; i < n; i++ {

		type nodeDist struct {
			id       int
			distance float64
		}

		candidates := make([]nodeDist, 0, n-1)
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			candidate := nodeDist{id: j, distance: distances[i][j]}
			candidates = append(candidates, candidate)
		}

		sort.Slice(candidates, func(x, y int) bool {
			return candidates[x].distance < candidates[y].distance
		})

		for j := 0; j < k; j++ {
			neighborsLists[i] = append(neighborsLists[i], candidates[j].id)
		}
	}

	return neighborsLists
}
