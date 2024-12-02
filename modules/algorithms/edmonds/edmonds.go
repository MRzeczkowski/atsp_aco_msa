package edmonds

import (
	"atsp_aco_msa/modules/models"
	"container/heap"
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
	pi := findMinimumEdges(_E, _w)

	// Step 3: Finding cycles using Union-Find
	cycleVertex := detectCycle(V, pi, r)

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
	V_prime, E_prime, w_prime, correspondence, v_c := contractCycle(V, _E, _w, pi, cycle, cycleVertex)

	// Recursive call
	tree := FindMSA(V_prime, E_prime, r, w_prime)

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

type EdgeHeap struct {
	edges []Edge
	w     map[Edge]float64 // Reference to the weights map
}

// Required heap.Interface methods
func (h EdgeHeap) Len() int           { return len(h.edges) }
func (h EdgeHeap) Less(i, j int) bool { return h.w[h.edges[i]] < h.w[h.edges[j]] }
func (h EdgeHeap) Swap(i, j int)      { h.edges[i], h.edges[j] = h.edges[j], h.edges[i] }

func (h *EdgeHeap) Push(x interface{}) {
	h.edges = append(h.edges, x.(Edge))
}

func (h *EdgeHeap) Pop() interface{} {
	old := h.edges
	n := len(old)
	x := old[n-1]
	h.edges = old[:n-1]
	return x
}

func findMinimumEdges(_E []Edge, w map[Edge]float64) map[int]int {
	// Initialize the priority queue with edges and weights
	edgeHeap := &EdgeHeap{
		edges: []Edge{},
		w:     w,
	}
	heap.Init(edgeHeap)

	// Populate the heap with edges
	for _, e := range _E {
		heap.Push(edgeHeap, e)
	}

	// Map to store parent pointers
	pi := make(map[int]int)
	selected := make(map[int]bool) // Track vertices that already have an incoming edge

	// Extract minimum incoming edges
	for edgeHeap.Len() > 0 {
		e := heap.Pop(edgeHeap).(Edge)
		if !selected[e.To] { // Process only if the vertex hasn't been assigned yet
			pi[e.To] = e.From
			selected[e.To] = true
		}
	}

	return pi
}

func detectCycle(V []int, pi map[int]int, root int) int {
	uf := newUnionFind(V) // Initialize Union-Find for |V| vertices

	cycleVertex := -1 // Default: no cycle detected
	for _, v := range V {
		if v == root {
			continue // Skip the root
		}

		parent, exists := pi[v]
		if !exists {
			continue // If there's no parent for this vertex, skip it
		}

		// Check if adding this edge forms a cycle
		if uf.find(v) == uf.find(parent) {
			cycleVertex = v // Cycle detected at this vertex
			break
		}

		// Otherwise, merge the two sets
		uf.union(v, parent)
	}

	return cycleVertex
}

func contractCycle(V []int, _E []Edge, _w map[Edge]float64, pi map[int]int, cycle map[int]bool, cycleVertex int) ([]int, []Edge, map[Edge]float64, map[Edge]Edge, int) {
	// Define the super-node for the cycle
	v_c := -(cycleVertex * cycleVertex) // Unique identifier for the super-node

	// Create new vertices list excluding the cycle vertices
	V_prime := []int{}
	for _, v := range V {
		if !cycle[v] {
			V_prime = append(V_prime, v)
		}
	}
	V_prime = append(V_prime, v_c) // Add the super-node

	// Process edges
	E_prime := []Edge{}
	w_prime := make(map[Edge]float64)
	correspondence := make(map[Edge]Edge)

	for _, uv := range _E {
		u := uv.From
		v := uv.To

		// Case 1: Incoming edge to the cycle
		if !cycle[u] && cycle[v] {
			e := Edge{From: u, To: v_c}
			if parent, exists := pi[v]; exists {
				tmpEdge := Edge{From: parent, To: v}
				if weight, exists := w_prime[e]; exists {
					// Keep the minimum weight edge
					if weight <= _w[uv]-_w[tmpEdge] {
						continue
					}
				}
				w_prime[e] = _w[uv] - _w[tmpEdge]
				correspondence[e] = uv
				E_prime = append(E_prime, e)
			}
		}

		// Case 2: Outgoing edge from the cycle
		if cycle[u] && !cycle[v] {
			e := Edge{From: v_c, To: v}
			if weight, exists := w_prime[e]; exists {
				// Keep the minimum weight edge
				if weight <= _w[uv] {
					continue
				}
			}
			w_prime[e] = _w[uv]
			correspondence[e] = uv
			E_prime = append(E_prime, e)
		}

		// Case 3: Edge outside the cycle
		if !cycle[u] && !cycle[v] {
			e := uv
			w_prime[e] = _w[uv]
			correspondence[e] = uv
			E_prime = append(E_prime, e)
		}
	}

	return V_prime, E_prime, w_prime, correspondence, v_c
}

type UnionFind struct {
	parent map[int]int
	size   map[int]int
	count  int
}

func newUnionFind(V []int) *UnionFind {
	numOfElements := len(V)
	parent := make(map[int]int, numOfElements)
	size := make(map[int]int, numOfElements)

	for _, v := range V {
		parent[v] = v
		size[v] = 1
	}

	return &UnionFind{
		parent: parent,
		size:   size,
		count:  numOfElements,
	}
}

// Time: O(logn) | Space: O(1)
func (uf *UnionFind) find(node int) int {
	for node != uf.parent[node] {
		// path compression
		uf.parent[node] = uf.parent[uf.parent[node]]
		node = uf.parent[node]
	}

	return node
}

// Time: O(1) | Space: O(1)
func (uf *UnionFind) union(node1, node2 int) {
	root1 := uf.find(node1)
	root2 := uf.find(node2)

	// already in the same set
	if root1 == root2 {
		return
	}

	if uf.size[root1] > uf.size[root2] {
		uf.parent[root2] = root1
		uf.size[root1] += 1
	} else {
		uf.parent[root1] = root2
		uf.size[root2] += 1
	}

	uf.count -= 1
}
