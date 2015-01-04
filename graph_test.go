package graphs

import (
	//"math/rand"
	"reflect"
	"testing"
	//"fmt"
)

/*func createRandomGraph(edges int64, undirected bool) (ug *Graph) {
	ug = &Graph{
		VertexEdges: make(map[uint64]map[uint64]int64),
	}

	for i := int64(0); i < edges; i++ {
		from := uint64(rand.Int63() % edges)
		to := uint64(rand.Int63() % edges)
		if from != to {
			if _, ok := ug.VertexEdges[from]; ok {
				ug.VertexEdges[from][to] = int64(rand.Int63())
			} else {
				ug.VertexEdges[from] = map[uint64]int64{to: uint64(rand.Int63())}
			}

			if undirected {
				if _, ok := ug.VertexEdges[to]; ok {
					ug.VertexEdges[to][from] = int64(rand.Int63())
				} else {
					ug.VertexEdges[to] = map[uint64]int64{from: uint64(rand.Int63())}
				}
			}
		}
	}

	return
}*/

// This test checks if we can get by DFS the two paths that connects all the
// elements in two separate graphs without any connection between them
func TestUndDFS(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 1},
			[]uint64{0, 2},
			[]uint64{1, 2},
			[]uint64{2, 3},
			[]uint64{2, 4},

			[]uint64{5, 6},
			[]uint64{6, 7},
			[]uint64{6, 9},
			[]uint64{9, 5},
		},
		true,
	)

	expectedFromZero := map[uint64]bool{
		0: true,
		1: true,
		2: true,
		3: true,
		4: true,
	}
	expectedFromFive := map[uint64]bool{
		5: true,
		6: true,
		7: true,
		9: true,
	}
	if !reflect.DeepEqual(gr.Dfs(0), expectedFromZero) {
		t.Error("Expeceted path from Zero:", expectedFromZero, "but:", gr.Dfs(0), "obtained.")
	}
	if !reflect.DeepEqual(gr.Dfs(5), expectedFromFive) {
		t.Error("Expeceted path from Five:", expectedFromFive, "but:", gr.Dfs(5), "obtained.")
	}
}

func TestUndConnectedComponents(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 1},
			[]uint64{0, 2},
			[]uint64{1, 2},
			[]uint64{2, 3},
			[]uint64{2, 4},

			[]uint64{5, 6},
			[]uint64{6, 7},
			[]uint64{6, 9},
			[]uint64{9, 5},
		},
		true,
	)

	expected := []map[uint64]bool{
		map[uint64]bool{
			0: true,
			1: true,
			2: true,
			3: true,
			4: true,
		},
		map[uint64]bool{
			5: true,
			6: true,
			7: true,
			9: true,
		},
	}

	comps := gr.ConnectedComponents()
	if len(comps) != len(expected) {
		t.Error("We expected:", len(expected), "components, but:", len(comps), "found")
	}

compLoop:
	for _, c := range comps {
		for _, ec := range expected {
			if reflect.DeepEqual(c, ec) {
				continue compLoop
			}
		}

		t.Error("No component found:", c)
	}
}

func TestUndBFS(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 1},
			[]uint64{0, 2},
			[]uint64{0, 5},
			[]uint64{1, 2},
			[]uint64{2, 3},
			[]uint64{2, 4},
			[]uint64{4, 3},
			[]uint64{3, 5},
		},
		true,
	)

	expectedDistances := map[uint64]uint64{
		0: 0,
		1: 1,
		2: 1,
		3: 2,
		4: 2,
		5: 1,
	}
	expectedPaths := map[uint64]uint64{
		0: 0,
		1: 0,
		2: 0,
		3: 2,
		4: 2,
		5: 0,
	}
	path, dist := gr.Bfs(0)
	// We have moultiple paths with the same length, so we will check that
	// the path is not longer than one that we know that is one of the
	// bests
	if len(path) != len(expectedPaths) {
		t.Error("Expeceted paths from Zero:", expectedPaths, "but:", path, "obtained.")
	}

	if !reflect.DeepEqual(dist, expectedDistances) {
		t.Error("Expeceted distances from Zero:", expectedDistances, "but:", dist, "obtained.")
	}
}

func TestUndBipartite(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{1, 6},
			[]uint64{2, 8},
			[]uint64{3, 8},
			[]uint64{4, 6},
			[]uint64{4, 9},
			[]uint64{5, 8},
			[]uint64{5, 9},
			[]uint64{7, 2},
			[]uint64{7, 3},

			[]uint64{10, 11},
			[]uint64{11, 12},
			[]uint64{12, 13},
			[]uint64{10, 12},
		},
		true,
	)

	if !gr.IsBipartite(1) {
		t.Error("The graph:", gr.VertexEdges, "is bipartite from vertex 1, but was not detected")
	}

	if gr.IsBipartite(10) {
		t.Error("The graph:", gr.VertexEdges, "souldn't be bipartite from edge 10")
	}
}

func TestUndEulerianCycle(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 1},
			[]uint64{1, 4},
			[]uint64{4, 2},
			[]uint64{1, 5},
			[]uint64{5, 2},
			[]uint64{1, 2},
			[]uint64{2, 3},
			[]uint64{3, 0},
		},
		true,
	)

	tour, possible := gr.EulerianCycle(0)
	if !possible {
		t.Error("For the specified graph:", gr.VertexEdges, "the Eulerian cycle was not detected from vertex 0")
	}
	if tour[0] != tour[len(tour)-1] {
		t.Error("For the specified graph:", gr.VertexEdges, "the Eulerian cycle doesn't starts or ends in the same vertex:", tour)
	}
	if len(tour) != 9 { // The returned vertices has to be the number of edges + 1
		t.Error("For the specified graph:", gr.VertexEdges, "the Eulerian cycle doesn't walks through all the edges:", tour)
	}
}

func TestUndEulerianPath(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 1},
			[]uint64{1, 4},
			[]uint64{4, 2},
			[]uint64{1, 5},
			[]uint64{5, 2},
			[]uint64{1, 2},
			[]uint64{2, 3},
			[]uint64{3, 0},
			[]uint64{0, 2},
		},
		true,
	)

	tour, possible := gr.EulerianPath(0, 2)
	if !possible {
		t.Error("For the specified graph:", gr.VertexEdges, "the Eulerian Path was not detected from vertex 0 to vertex 2")
	}
	if tour[0] != 0 || tour[len(tour)-1] != 2 {
		t.Error("For the specified graph:", gr.VertexEdges, "the Eulerian Path doesn't starts on vertex 0, or ends on vertex 2:", tour)
	}
	if len(tour) != 10 { // The returned vertices has to be the number of edges + 1
		t.Error("For the specified graph:", gr.VertexEdges, "the Eulerian path doesn't walks through all the edges:", tour)
	}
}

func TestUndHamiltonPath(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{2, 3},
			[]uint64{3, 4},
			[]uint64{4, 5},
			[]uint64{5, 6},
			[]uint64{6, 7},
			[]uint64{7, 8},
			[]uint64{8, 9},
			[]uint64{9, 10},
			[]uint64{10, 11},
			[]uint64{11, 2},

			[]uint64{1, 2},
			[]uint64{4, 12},
			[]uint64{6, 13},
			[]uint64{8, 14},
			[]uint64{10, 15},

			[]uint64{1, 12},
			[]uint64{12, 13},
			[]uint64{13, 14},
			[]uint64{14, 15},
			[]uint64{15, 1},

			[]uint64{11, 16},
			[]uint64{9, 20},
			[]uint64{7, 19},
			[]uint64{5, 18},
			[]uint64{3, 17},

			[]uint64{16, 17},
			[]uint64{17, 18},
			[]uint64{18, 19},
			[]uint64{19, 20},
			[]uint64{20, 16},
		},
		true,
	)

	tour, possible := gr.HamiltonianPath(1, 2)

	if !possible {
		t.Error("Hamilton path not found for origin 1 and dest 2 on graph:", gr.VertexEdges)
	}
	if len(tour) != len(gr.VertexEdges) {
		t.Error("Hamilton path", tour, "doesn't covers all the vertices of the graph:", gr.VertexEdges)
	}

	tour, possible = gr.HamiltonianPath(1, 1)
	if !possible {
		t.Error("Hamilton tour not found for origin 1 and dest 1 on graph:", gr.VertexEdges)
	}
	if len(tour) != len(gr.VertexEdges)+1 {
		t.Error("Hamilton tour", tour, "doesn't covers all the vertices of the graph:", gr.VertexEdges)
	}
}

func TestTopologicalOrder(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 1},
			[]uint64{0, 5},
			[]uint64{0, 2},
			[]uint64{1, 4},
			[]uint64{5, 2},
			[]uint64{3, 2},
			[]uint64{3, 5},
			[]uint64{3, 4},
			[]uint64{3, 6},
			[]uint64{6, 0},
			[]uint64{6, 4},
		},
		false,
	)

	order, success := gr.TopologicalOrder()
	if !success {
		t.Error("Problem calculating topological order on graph:", gr.VertexEdges)
	}
	if len(order) != len(gr.Vertices) {
		t.Error("The number of vertices in the specified order:", order, "doesn't match with the total vertices on the graph:", gr.Vertices)
	}
}

func TestTopologicalCycle(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 1},
			[]uint64{0, 5},
			[]uint64{0, 2},
			[]uint64{1, 4},
			[]uint64{5, 2},
			[]uint64{3, 2},
			[]uint64{3, 5},
			[]uint64{3, 4},
			[]uint64{3, 6},
			[]uint64{6, 0},
			[]uint64{6, 4},
			[]uint64{2, 3},
		},
		false,
	)

	_, success := gr.TopologicalOrder()
	if success {
		t.Error("Problem calculating topological order on graph:", gr.VertexEdges, "a graph with a cycle can't have topological order")
	}
}

func TestStrongConnectedComponents(t *testing.T) {
	gr := GetUnWeightGraph(
		[][]uint64{
			[]uint64{0, 6},
			[]uint64{0, 2},
			[]uint64{1, 0},
			[]uint64{2, 3},
			[]uint64{2, 4},
			[]uint64{3, 2},
			[]uint64{3, 4},
			[]uint64{4, 5},
			[]uint64{4, 6},
			[]uint64{4, 11},
			[]uint64{5, 3},
			[]uint64{5, 0},
			[]uint64{6, 7},
			[]uint64{6, 8},
			[]uint64{8, 6},
			[]uint64{9, 7},
			[]uint64{9, 6},
			[]uint64{9, 12},
			[]uint64{10, 9},
			[]uint64{11, 9},
			[]uint64{12, 10},
			[]uint64{12, 11},
		},
		false,
	)

	comps, groups := gr.StronglyConnectedComponents()
	if len(groups) != 5 {
		t.Error("We have five strong components on the graph:", gr.VertexEdges, ", but:", len(groups), "was detected")
	}
	if comps[0] != comps[3] {
		t.Error("The components 0 and 3 should to be in the same group, but was not detected as it")
	}
	if comps[11] != comps[9] {
		t.Error("The components 0 and 3 should to be in the same group, but was not detected as it")
	}
}

func TestMst(t *testing.T) {
	gr := GetGraph(
		[]EdgeDefinition{
			EdgeDefinition{1, 7, 19},
			EdgeDefinition{0, 2, 26},
			EdgeDefinition{1, 3, 29},
			EdgeDefinition{2, 3, 17},
			EdgeDefinition{5, 7, 28},
			EdgeDefinition{2, 7, 34},
			EdgeDefinition{6, 4, 93},
			EdgeDefinition{4, 5, 35},
			EdgeDefinition{1, 5, 32},
			EdgeDefinition{1, 2, 36},
			EdgeDefinition{0, 4, 38},
			EdgeDefinition{4, 7, 37},
			EdgeDefinition{6, 2, 40},
			EdgeDefinition{3, 6, 52},
			EdgeDefinition{0, 7, 16},
			EdgeDefinition{6, 0, 58},
		},
		false,
	)

	expectedResult := []EdgeDefinition {
		EdgeDefinition{0, 7, 16},
		EdgeDefinition{2, 3, 17},
		EdgeDefinition{1, 7, 19},
		EdgeDefinition{0, 2, 26},
		EdgeDefinition{5, 7, 28},
		EdgeDefinition{4, 5, 35},
		EdgeDefinition{6, 2, 40},
	}
	result := gr.Mst()

	if !reflect.DeepEqual(result, expectedResult) {
		t.Error("Expeceted MST:", expectedResult, "but:", result, "obtained.")
	}
}
