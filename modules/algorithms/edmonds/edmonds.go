package edmonds

import (
	"atsp_aco_msa/modules/models"
	"math"
)

type Edge = models.Edge

func FindMSA(root int, vertices []int, edges []Edge, weights map[Edge]float64) []Edge {

	// Step 1: Removing all edges leading back to the root and adjusting edge set
	filteredEdges := make([]Edge, 0, len(edges))
	filteredWeights := make(map[Edge]float64, len(edges))

	for _, edge := range edges {
		if edge.To != root {
			filteredEdges = append(filteredEdges, edge)
			filteredWeights[edge] = weights[edge]
		}
	}

	// Step 2: Finding minimum incoming edge for each vertex
	parentPointers := getParentPointers(root, vertices, filteredEdges, filteredWeights)

	// Step 3: Finding cycles
	cycleVertex, found := getCycleVertex(vertices, parentPointers)

	// Step 4: No cycle
	if !found {
		msaEdges := make([]Edge, 0, len(parentPointers))

		// Build the Minimum Spanning Arborescence from parentPointers
		for vertex, parent := range parentPointers {
			msaEdges = append(msaEdges, Edge{From: parent, To: vertex})
		}

		return msaEdges
	}

	// Step 5: Handle Cycle
	cycleVertices := make(map[int]bool) // Map to track vertices in the cycle
	cycleVertices[cycleVertex] = true   // Add the starting vertex of the cycle

	nextVertex := parentPointers[cycleVertex]
	for nextVertex != cycleVertex {
		cycleVertices[nextVertex] = true // Mark vertex as part of the cycle
		nextVertex = parentPointers[nextVertex]
	}

	// Step 6: Contract the cycle into a new node
	superNode := -(cycleVertex * cycleVertex) // Unique identifier for the super-node

	// Create new vertices list excluding the cycle vertices and including super-node
	contractedVerticesCount := len(vertices) - len(cycleVertices) + 1
	contractedVertices := make([]int, 0, contractedVerticesCount)
	for _, vertex := range vertices {
		if !cycleVertices[vertex] {
			contractedVertices = append(contractedVertices, vertex)
		}
	}
	contractedVertices = append(contractedVertices, superNode) // Add the super-node

	// Estimate the worst-case number of edges in the contracted graph
	maxEdges := contractedVerticesCount * (contractedVerticesCount - 1)

	// Initialize contracted edges, weights, and correspondence mapping
	contractedEdges := make([]Edge, 0, maxEdges) // Preallocate based on contracted vertices
	contractedWeights := make(map[Edge]float64, maxEdges)
	edgeMapping := make(map[Edge]Edge, maxEdges)

	for _, edge := range filteredEdges {
		sourceVertex := edge.From
		destinationVertex := edge.To

		// Case 1: Incoming edge to the cycle
		if !cycleVertices[sourceVertex] && cycleVertices[destinationVertex] {
			contractedEdge := Edge{From: sourceVertex, To: superNode}

			// Retrieve the minimum incoming edge to the cycle vertex
			if parentVertex, exists := parentPointers[destinationVertex]; exists {
				minCycleEdge := Edge{From: parentVertex, To: destinationVertex}

				// Calculate the relative weight adjustment
				relativeWeight := filteredWeights[edge] - filteredWeights[minCycleEdge]

				// Check if a contracted edge already exists and compare weights, keep the lighter one
				if weight, exists := contractedWeights[contractedEdge]; exists && weight <= relativeWeight {
					continue
				}

				contractedWeights[contractedEdge] = relativeWeight
				edgeMapping[contractedEdge] = edge
				contractedEdges = append(contractedEdges, contractedEdge)
			}
		}

		// Case 2: Outgoing edge from the cycle
		if cycleVertices[sourceVertex] && !cycleVertices[destinationVertex] {
			contractedEdge := Edge{From: superNode, To: destinationVertex}

			// Check if a contracted edge already exists and compare weights, keep the lighter one
			if weight, exists := contractedWeights[contractedEdge]; exists && weight <= filteredWeights[edge] {
				continue
			}

			contractedWeights[contractedEdge] = filteredWeights[edge]
			edgeMapping[contractedEdge] = edge
			contractedEdges = append(contractedEdges, contractedEdge)
		}

		// Case 3: Edge outside the cycle
		if !cycleVertices[sourceVertex] && !cycleVertices[destinationVertex] {

			// Keeping the edge and its weight unchanged, 1:1 mapping and adding as is to the contracted graph
			contractedWeights[edge] = filteredWeights[edge]
			edgeMapping[edge] = edge
			contractedEdges = append(contractedEdges, edge)
		}
	}

	// Step 7: Recursive call
	// Solving the same problem but without the cycle and using the super-node.
	contractedTree := FindMSA(root, contractedVertices, contractedEdges, contractedWeights)

	// Step 8: Expanding back
	var connectingEdgeToCycle Edge // Edge connecting the cycle to the rest of the graph

	// Identify the edge in the MSA solution that connects to the super-node
	for _, contractedEdge := range contractedTree {
		sourceVertex := contractedEdge.From
		targetVertex := contractedEdge.To

		if targetVertex == superNode {
			// Retrieve the original edge corresponding to this contracted edge
			originalTarget := edgeMapping[Edge{From: sourceVertex, To: superNode}].To

			// Determine the internal cycle edge to replace the contracted edge
			connectingEdgeToCycle = Edge{From: parentPointers[originalTarget], To: originalTarget}
			break
		}
	}

	// Translate the edges from the contracted graph back to the original graph
	expandedSetSize := len(contractedTree) + len(cycleVertices) - 1

	expandedEdgeSet := make(map[Edge]bool, expandedSetSize)
	for _, contractedEdge := range contractedTree {
		expandedEdgeSet[edgeMapping[contractedEdge]] = true
	}

	// Add the edges within the cycle to the result set
	for cycleVertex := range cycleVertices {
		parentVertex := parentPointers[cycleVertex]
		cycleEdge := Edge{From: parentVertex, To: cycleVertex}
		expandedEdgeSet[cycleEdge] = true
	}

	// Remove the edge that connected the cycle to the super-node
	delete(expandedEdgeSet, connectingEdgeToCycle)

	// Convert the edge set back into a slice for the final result
	finalMSA := make([]Edge, 0, expandedSetSize)
	for edge := range expandedEdgeSet {
		finalMSA = append(finalMSA, edge)
	}

	return finalMSA
}

func getParentPointers(root int, vertices []int, filteredEdges []Edge, filteredWeights map[Edge]float64) map[int]int {
	// This map tracks the parent of each vertex in the current minimum spanning tree under construction.
	parentPointers := make(map[int]int, len(vertices)-1)

	for _, vertex := range vertices {
		if vertex == root {
			continue
		}

		minCost := math.MaxFloat64

		// Find the minimum incoming edge for this vertex
		for _, edge := range filteredEdges {
			if edge.To == vertex && filteredWeights[edge] < minCost {
				minCost = filteredWeights[edge]
				parentPointers[vertex] = edge.From
			}
		}
	}

	return parentPointers
}

func getCycleVertex(vertices []int, parentPointers map[int]int) (int, bool) {
	var visited map[int]bool

	for _, vertex := range vertices {
		visited = make(map[int]bool, len(vertices)) // Track visited vertices for this traversal
		nextVertex, exists := parentPointers[vertex]

		for exists {
			if visited[nextVertex] {
				return nextVertex, true // Cycle detected
			}

			visited[nextVertex] = true
			nextVertex, exists = parentPointers[nextVertex]
		}
	}

	return -1, false // No cycle found
}
