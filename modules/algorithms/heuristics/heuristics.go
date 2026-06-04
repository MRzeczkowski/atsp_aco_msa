package heuristics

import (
	"math"
	"math/rand"
)

const msaHeuristicHighSignalThreshold = 1.0

type directedEdge struct {
	from, to int
}

type cyclePatch struct {
	fromA, toB   int
	fromB, toA   int
	costIncrease float64
	msaSupport   float64
	score        float64
	valid        bool
}

type patchCostCache struct {
	matrix, msaHeuristic [][]float64
	successors           []int
	cycleIndex           []int
	msaPatchBias         float64
	costs                [][]cyclePatch
	bestInRow            []cyclePatch
	dirtyRows            []bool
}

func BuildMsaHeuristicModifiers(msaHeuristic [][]float64, strength float64) [][]float64 {
	dimension := len(msaHeuristic)
	modifiers := BuildNeutralModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	maxMsaHeuristicSelections := float64(dimension - 1)
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i == j {
				continue
			}

			msaHeuristicSignal := msaHeuristic[i][j] / maxMsaHeuristicSelections
			if msaHeuristicSignal >= msaHeuristicHighSignalThreshold {
				modifiers[i][j] = 1.0 + msaHeuristicSignal*strength
			}
		}
	}

	return modifiers
}

func BuildRandomSparseModifiers(msaHeuristic [][]float64, strength float64, seed int64) [][]float64 {
	dimension := len(msaHeuristic)
	modifiers := BuildNeutralModifiers(dimension)
	if dimension <= 1 || strength == 0 {
		return modifiers
	}

	boostedEdgeCount := boostedModifierEdgeCount(BuildMsaHeuristicModifiers(msaHeuristic, 1.0))
	if boostedEdgeCount == 0 {
		return modifiers
	}

	edges := make([]directedEdge, 0, dimension*(dimension-1))
	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i != j {
				edges = append(edges, directedEdge{i, j})
			}
		}
	}

	random := rand.New(rand.NewSource(seed))
	for i := len(edges) - 1; i > 0; i-- {
		j := random.Intn(i + 1)
		edges[i], edges[j] = edges[j], edges[i]
	}

	for _, edge := range edges[:boostedEdgeCount] {
		modifiers[edge.from][edge.to] = 1.0 + strength
	}

	return modifiers
}

func BuildCycleCoverMsaPatchingModifiers(matrix, msaHeuristic, cycleCover [][]float64, patchingWeight, msaPatchBias float64) [][]float64 {
	if patchingWeight == 0 {
		return BuildNeutralModifiers(heuristicDimension(matrix, msaHeuristic, cycleCover))
	}

	patchingMatrix := BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(matrix, msaHeuristic, cycleCover, msaPatchBias)
	modifiers := BuildNeutralModifiers(len(patchingMatrix))

	for i := 0; i < len(patchingMatrix); i++ {
		for j := 0; j < len(patchingMatrix[i]); j++ {
			if i != j && patchingMatrix[i][j] > 0 {
				modifiers[i][j] = 1.0 + patchingMatrix[i][j]*patchingWeight
			}
		}
	}

	return modifiers
}

func BuildCycleCoverMsaPatchingMatrix(matrix, msaHeuristic, cycleCover [][]float64) [][]float64 {
	return BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(matrix, msaHeuristic, cycleCover, 1.0)
}

func BuildCycleCoverMsaPatchingMatrixWithMsaPatchBias(matrix, msaHeuristic, cycleCover [][]float64, msaPatchBias float64) [][]float64 {
	dimension := heuristicDimension(matrix, msaHeuristic, cycleCover)
	if dimension == 0 {
		return newZeroMatrix(dimension)
	}

	successors, ok := cycleCoverSuccessors(cycleCover, dimension)
	if !ok {
		return copyPositiveEdges(cycleCover, dimension)
	}

	cycles := buildCycles(successors)
	cycleIndex := vertexCycleIndexes(cycles, dimension)
	costCache := newPatchCostCache(matrix, msaHeuristic, successors, cycleIndex, msaPatchBias)

	activeCycles := len(cycles)
	for activeCycles > 1 {
		// GKS-style greedy patching process:
		// 1. cache patch costs for all vertex pairs,
		// 2. choose the best patch between two different current cycles,
		// 3. update only costs made stale by the merged cycles and swapped successors.
		patch := costCache.bestPatch()
		if !patch.valid {
			return copyPositiveEdges(cycleCover, dimension)
		}

		fromCycle := cycleIndex[patch.fromA]
		toCycle := cycleIndex[patch.fromB]
		fromCycleVertices := cycles[fromCycle]
		toCycleVertices := cycles[toCycle]

		applyCyclePatch(successors, patch)
		for _, vertex := range toCycleVertices {
			cycleIndex[vertex] = fromCycle
		}
		cycles[fromCycle] = append(cycles[fromCycle], toCycleVertices...)
		cycles[toCycle] = nil
		activeCycles--

		costCache.updateAfterPatch(patch, fromCycleVertices, toCycleVertices, fromCycle)
	}

	return buildSuccessorMatrix(successors)
}

func BuildCycleCoverModifiers(cycleCover [][]float64, strength float64) [][]float64 {
	dimension := len(cycleCover)
	modifiers := BuildNeutralModifiers(dimension)
	if dimension == 0 || strength == 0 {
		return modifiers
	}

	for i := 0; i < dimension; i++ {
		for j := 0; j < dimension; j++ {
			if i != j && cycleCover[i][j] != 0 {
				modifiers[i][j] = 1.0 + strength
			}
		}
	}

	return modifiers
}

func boostedModifierEdgeCount(modifiers [][]float64) int {
	count := 0
	for i := 0; i < len(modifiers); i++ {
		for j := 0; j < len(modifiers[i]); j++ {
			if i != j && modifiers[i][j] > 1.0 {
				count++
			}
		}
	}

	return count
}

func BuildNeutralModifiers(dimension int) [][]float64 {
	modifiers := make([][]float64, dimension)
	for i := range modifiers {
		modifiers[i] = make([]float64, dimension)
		for j := range modifiers[i] {
			modifiers[i][j] = 1.0
		}
	}
	return modifiers
}

func heuristicDimension(matrices ...[][]float64) int {
	for _, matrix := range matrices {
		if len(matrix) != 0 {
			return len(matrix)
		}
	}

	return 0
}

func newZeroMatrix(dimension int) [][]float64 {
	matrix := make([][]float64, dimension)
	for i := range matrix {
		matrix[i] = make([]float64, dimension)
	}

	return matrix
}

func copyPositiveEdges(matrix [][]float64, dimension int) [][]float64 {
	copyMatrix := newZeroMatrix(dimension)
	for i := 0; i < dimension && i < len(matrix); i++ {
		for j := 0; j < dimension && j < len(matrix[i]); j++ {
			if i != j && matrix[i][j] != 0 {
				copyMatrix[i][j] = 1.0
			}
		}
	}

	return copyMatrix
}

func cycleCoverSuccessors(cycleCover [][]float64, dimension int) ([]int, bool) {
	successors := make([]int, dimension)
	inDegree := make([]int, dimension)
	for vertex := range successors {
		successor := cycleCoverSuccessor(cycleCover, vertex)
		if successor < 0 || successor >= dimension {
			return nil, false
		}
		successors[vertex] = successor
		inDegree[successor]++
	}

	for _, degree := range inDegree {
		if degree != 1 {
			return nil, false
		}
	}

	return successors, true
}

func cycleCoverSuccessor(cycleCover [][]float64, vertex int) int {
	if vertex < 0 || vertex >= len(cycleCover) {
		return -1
	}

	for to, value := range cycleCover[vertex] {
		if vertex != to && value != 0 {
			return to
		}
	}

	return -1
}

func buildCycles(successors []int) [][]int {
	visited := make([]bool, len(successors))
	cycles := make([][]int, 0)

	for start := 0; start < len(successors); start++ {
		if visited[start] {
			continue
		}

		cycle := make([]int, 0)
		current := start
		for !visited[current] {
			visited[current] = true
			cycle = append(cycle, current)
			current = successors[current]
		}
		cycles = append(cycles, cycle)
	}

	return cycles
}

func bestGreedyCyclePatch(matrix, msaHeuristic [][]float64, successors []int, cycles [][]int, msaPatchBias float64) cyclePatch {
	best := cyclePatch{}
	cycleIndex := vertexCycleIndexes(cycles, len(successors))

	for fromA := range successors {
		for fromB := fromA + 1; fromB < len(successors); fromB++ {
			if cycleIndex[fromA] == cycleIndex[fromB] {
				continue
			}

			patch := newCyclePatch(matrix, msaHeuristic, successors, fromA, fromB, msaPatchBias)
			if betterCyclePatch(patch, best) {
				best = patch
			}
		}
	}

	return best
}

func newPatchCostCache(matrix, msaHeuristic [][]float64, successors, cycleIndex []int, msaPatchBias float64) *patchCostCache {
	dimension := len(successors)
	cache := &patchCostCache{
		matrix:       matrix,
		msaHeuristic: msaHeuristic,
		successors:   successors,
		cycleIndex:   cycleIndex,
		msaPatchBias: msaPatchBias,
		costs:        make([][]cyclePatch, dimension),
		bestInRow:    make([]cyclePatch, dimension),
		dirtyRows:    make([]bool, dimension),
	}

	for row := 0; row < dimension; row++ {
		cache.costs[row] = make([]cyclePatch, row)
		for col := 0; col < row; col++ {
			cache.costs[row][col] = cache.newCachedPatch(row, col)
			if betterCyclePatch(cache.costs[row][col], cache.bestInRow[row]) {
				cache.bestInRow[row] = cache.costs[row][col]
			}
		}
	}

	return cache
}

func (cache *patchCostCache) bestPatch() cyclePatch {
	cache.recalculateDirtyRows()

	best := cyclePatch{}
	for _, patch := range cache.bestInRow {
		if betterCyclePatch(patch, best) {
			best = patch
		}
	}
	return best
}

func (cache *patchCostCache) updateAfterPatch(patch cyclePatch, fromCycleVertices, toCycleVertices []int, mergedCycle int) {
	for _, from := range fromCycleVertices {
		for _, to := range toCycleVertices {
			cache.updatePair(from, to, cyclePatch{})
		}
	}

	for vertex := range cache.successors {
		if cache.cycleIndex[vertex] == mergedCycle {
			continue
		}

		cache.updatePair(patch.fromA, vertex, cache.newCachedPatchForPair(patch.fromA, vertex))
		cache.updatePair(patch.fromB, vertex, cache.newCachedPatchForPair(patch.fromB, vertex))
	}

	cache.recalculateDirtyRows()
}

func (cache *patchCostCache) updatePair(a, b int, patch cyclePatch) {
	if a == b {
		return
	}

	row, col := orderedPair(a, b)
	cache.costs[row][col] = patch
	if betterCyclePatch(patch, cache.bestInRow[row]) {
		cache.bestInRow[row] = patch
		return
	}

	if samePatchEndpoints(cache.bestInRow[row], a, b) {
		cache.dirtyRows[row] = true
	}
}

func (cache *patchCostCache) recalculateDirtyRows() {
	for row, dirty := range cache.dirtyRows {
		if !dirty {
			continue
		}

		cache.recalculateRow(row)
		cache.dirtyRows[row] = false
	}
}

func (cache *patchCostCache) recalculateRow(row int) {
	best := cyclePatch{}
	for col := range cache.costs[row] {
		if betterCyclePatch(cache.costs[row][col], best) {
			best = cache.costs[row][col]
		}
	}
	cache.bestInRow[row] = best
}

func (cache *patchCostCache) newCachedPatchForPair(a, b int) cyclePatch {
	row, col := orderedPair(a, b)
	return cache.newCachedPatch(row, col)
}

func (cache *patchCostCache) newCachedPatch(row, col int) cyclePatch {
	if cache.cycleIndex[row] == cache.cycleIndex[col] {
		return cyclePatch{}
	}

	return newCyclePatch(cache.matrix, cache.msaHeuristic, cache.successors, col, row, cache.msaPatchBias)
}

func orderedPair(a, b int) (int, int) {
	if a > b {
		return a, b
	}
	return b, a
}

func samePatchEndpoints(patch cyclePatch, a, b int) bool {
	if !patch.valid {
		return false
	}

	return (patch.fromA == a && patch.fromB == b) || (patch.fromA == b && patch.fromB == a)
}

func vertexCycleIndexes(cycles [][]int, dimension int) []int {
	indexes := make([]int, dimension)
	for i := range indexes {
		indexes[i] = -1
	}
	for cycleIndex, cycle := range cycles {
		for _, vertex := range cycle {
			if vertex >= 0 && vertex < dimension {
				indexes[vertex] = cycleIndex
			}
		}
	}

	return indexes
}

func newCyclePatch(matrix, msaHeuristic [][]float64, successors []int, fromA, fromB int, msaPatchBias float64) cyclePatch {
	toA := successors[fromA]
	toB := successors[fromB]
	costIncrease := matrixDistance(matrix, fromA, toB) +
		matrixDistance(matrix, fromB, toA) -
		matrixDistance(matrix, fromA, toA) -
		matrixDistance(matrix, fromB, toB)
	if math.IsInf(costIncrease, 0) || math.IsNaN(costIncrease) {
		return cyclePatch{}
	}

	msaSupport := msaHeuristicSignal(msaHeuristic, fromA, toB) +
		msaHeuristicSignal(msaHeuristic, fromB, toA)
	weightedMsaSupport := msaSupport * msaPatchBias
	score := costIncrease / (1.0 + weightedMsaSupport)
	if math.IsInf(score, 0) || math.IsNaN(score) {
		return cyclePatch{}
	}

	return cyclePatch{
		fromA:        fromA,
		toB:          toB,
		fromB:        fromB,
		toA:          toA,
		costIncrease: costIncrease,
		msaSupport:   weightedMsaSupport,
		score:        score,
		valid:        true,
	}
}

func applyCyclePatch(successors []int, patch cyclePatch) {
	successors[patch.fromA] = patch.toB
	successors[patch.fromB] = patch.toA
}

func msaHeuristicSignal(msaHeuristic [][]float64, from, to int) float64 {
	maxMsaHeuristicSelections := float64(len(msaHeuristic) - 1)
	if maxMsaHeuristicSelections <= 0 ||
		from < 0 || from >= len(msaHeuristic) ||
		to < 0 || to >= len(msaHeuristic[from]) {
		return 0
	}

	return msaHeuristic[from][to] / maxMsaHeuristicSelections
}

func betterCyclePatch(candidate, current cyclePatch) bool {
	if !candidate.valid {
		return false
	}
	if !current.valid {
		return true
	}
	if candidate.score != current.score {
		return candidate.score < current.score
	}
	if candidate.msaSupport != current.msaSupport {
		return candidate.msaSupport > current.msaSupport
	}
	if candidate.costIncrease != current.costIncrease {
		return candidate.costIncrease < current.costIncrease
	}
	if candidate.fromA != current.fromA {
		return candidate.fromA < current.fromA
	}
	if candidate.toB != current.toB {
		return candidate.toB < current.toB
	}
	if candidate.fromB != current.fromB {
		return candidate.fromB < current.fromB
	}
	return candidate.toA < current.toA
}

func matrixDistance(matrix [][]float64, from, to int) float64 {
	if from >= 0 && from < len(matrix) && to >= 0 && to < len(matrix[from]) {
		return matrix[from][to]
	}

	return math.Inf(1)
}

func buildSuccessorMatrix(successors []int) [][]float64 {
	matrix := newZeroMatrix(len(successors))
	for from, to := range successors {
		if from != to && to >= 0 && to < len(successors) {
			matrix[from][to] = 1.0
		}
	}

	return matrix
}
