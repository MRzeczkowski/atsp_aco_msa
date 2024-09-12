package edmonds

import (
	"atsp_aco_msa/modules/models"
	"math"

	"golang.org/x/exp/maps"
)

type Edge = models.Edge

func FindMSA(V []int, E []Edge, r int, w map[Edge]float64) []Edge {
	// Step 1: Removing all edges leading back to the root and adjusting edge set
	var _E []Edge
	_w := make(map[Edge]float64, len(w))
	for _, e := range E {
		if e.To != r {
			_E = append(_E, e)
			_w[e] = w[e]
		}
	}

	// Step 2: Finding minimum incoming edge for each vertex
	pi := make(map[int]int)
	for _, v := range V {
		if v == r {
			continue
		}
		minCost := math.MaxFloat64
		for _, e := range _E {
			if e.To == v && _w[e] < minCost {
				minCost = _w[e]
				pi[v] = e.From
			}
		}
	}

	// Step 3: Finding cycles
	cycleVertex := -1
	var visited map[int]bool
	for _, v := range V {
		if cycleVertex != -1 {
			break
		}

		visited = make(map[int]bool)
		next_v, ok := pi[v]

		for ok {
			if visited[next_v] {
				cycleVertex = next_v
				break
			}

			visited[next_v] = true
			next_v, ok = pi[next_v]
		}
	}

	var result []Edge

	// Step 4: No cycle
	if cycleVertex == -1 {
		for v, u := range pi {
			result = append(result, Edge{From: u, To: v})
		}

		return result
	}

	// Step 5: Handle cycle
	cycle := make(map[int]bool)
	cycle[cycleVertex] = true
	next_v := pi[cycleVertex]
	for next_v != cycleVertex {
		cycle[next_v] = true
		next_v = pi[next_v]
	}

	// Step 6: Contract the cycle into a new node v_c
	v_c := -(cycleVertex * cycleVertex) // Unique negative squared identifier
	V_prime := []int{}
	for _, v := range V {
		if !cycle[v] {
			V_prime = append(V_prime, v)
		}
	}

	V_prime = append(V_prime, v_c)
	E_prime := make(map[Edge]bool)
	w_prime := make(map[Edge]float64)
	correspondence := make(map[Edge]Edge)

	for _, uv := range _E {
		u := uv.From
		v := uv.To

		if !cycle[u] && cycle[v] {
			e := Edge{From: u, To: v_c}
			tmpEdge := Edge{From: pi[v], To: v}
			if E_prime[e] {
				if w_prime[e] < _w[uv]-_w[tmpEdge] {
					continue
				}
			}

			w_prime[e] = _w[uv] - _w[tmpEdge]
			correspondence[e] = uv
			E_prime[e] = true
		} else if cycle[u] && !cycle[v] {
			e := Edge{From: v_c, To: v}
			if E_prime[e] {
				old_u := correspondence[e].From

				tmpEdge := Edge{From: old_u, To: v}
				if _w[tmpEdge] < _w[uv] {
					continue
				}
			}

			E_prime[e] = true
			w_prime[e] = _w[uv]
			correspondence[e] = uv
		} else if !cycle[u] && !cycle[v] {
			e := uv
			E_prime[e] = true
			w_prime[e] = _w[uv]
			correspondence[e] = uv
		}
	}

	// Recursive call
	tree := FindMSA(V_prime, maps.Keys(E_prime), r, w_prime)

	// Step 8: Expanding back
	var cycle_edge Edge

	for _, e := range tree {
		u := e.From
		v := e.To

		if v == v_c {
			tmpEdge := Edge{From: u, To: v_c}
			old_v := correspondence[tmpEdge].To
			cycle_edge = Edge{From: pi[old_v], To: old_v}
			break
		}
	}

	resultSet := make(map[Edge]bool)

	for _, uv := range tree {
		resultSet[correspondence[uv]] = true
	}

	for v := range cycle {
		u := pi[v]
		tmpEdge := Edge{From: u, To: v}
		resultSet[tmpEdge] = true
	}

	delete(resultSet, cycle_edge)

	result = make([]Edge, 0)
	for e := range resultSet {
		result = append(result, e)
	}

	return result
}
