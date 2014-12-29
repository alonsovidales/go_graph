package graphs

import (
	"sort"
)

// Graph Data strcture used to represent a graph, the VertexEdges var
// is a map where each key is a vertex, and the value a map where the keys are
// the vertices that can be reached from the main key vertex and the value the
// weight for this edge, the Vertices property is used as a set who contains
// all the available vertices in the graph
type Graph struct {
	Vertices    map[uint64]bool
	VertexEdges map[uint64]map[uint64]uint64
	Undirected  bool
}

// ByWeight Used to sort the graph edges by weight
type ByWeight [][]uint64

func (a ByWeight) Len() int           { return len(a) }
func (a ByWeight) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByWeight) Less(i, j int) bool { return a[i][2] < a[j][2] }

// GetGraph Returns an unweighted graph containing the specified edges,
// use the second boolean parameter in order to specify if the graph to be
// constructed is directed (true) or undirected (false)
func GetGraph(edges [][]uint64, undirected bool) (ug *Graph) {
	ug = &Graph{
		Vertices:    make(map[uint64]bool),
		VertexEdges: make(map[uint64]map[uint64]uint64),
		Undirected:  undirected,
	}

	weight := uint64(1)
	for _, edge := range edges {
		if len(edge) == 3 {
			weight = edge[2]
		}
		ug.Vertices[edge[0]] = true
		ug.Vertices[edge[1]] = true
		if _, ok := ug.VertexEdges[edge[0]]; ok {
			ug.VertexEdges[edge[0]][edge[1]] = weight
		} else {
			ug.VertexEdges[edge[0]] = map[uint64]uint64{edge[1]: weight}
		}
		if undirected {
			if _, ok := ug.VertexEdges[edge[1]]; ok {
				ug.VertexEdges[edge[1]][edge[0]] = weight
			} else {
				ug.VertexEdges[edge[1]] = map[uint64]uint64{edge[0]: weight}
			}
		}
	}

	return
}

// Mst Minimum Spanning Tree, Calculates the tree of edges who connects all the vertices with a minimun cost
// This method uses the Kruskal's algorithm:
//	- http://en.wikipedia.org/wiki/Kruskal%27s_algorithm
// An Union-find in order to detect cycles:
//	- http://en.wikipedia.org/wiki/Disjoint-set_data_structure
func (gr *Graph) Mst() (mst [][]uint64) {
	queue := [][]uint64{}
	mst = [][]uint64{}
	for v, e := range gr.VertexEdges {
		for tv, w := range e {
			queue = append(queue, []uint64{v, tv, w})
		}
	}

	var edgeToAdd, groupId uint64

	// Using union-find algorithm to detect cycles
	sort.Sort(ByWeight(queue))
	vertexByGroup := make(map[uint64][]uint64)
	vertexGroups := make(map[uint64]uint64)
	connect := make([]uint64, 2)
	lastUsedGroup := uint64(0)
queueLoop:
	for _, e := range queue {
		addToExistingGroup := false
		lG, lIn := vertexGroups[e[0]]
		rG, rIn := vertexGroups[e[1]]
		switch {
		case lIn && rIn:
			// We have a vertex :'(
			if lG == rG {
				continue queueLoop
			}

			// We will connect two groups
			if len(vertexByGroup[lG]) < len(vertexByGroup[rG]) {
				connect[0] = rG
				connect[1] = lG
			} else {
				connect[0] = lG
				connect[1] = rG
			}

			for _, v := range vertexByGroup[connect[1]] {
				vertexByGroup[connect[0]] = append(vertexByGroup[connect[0]], v)
				vertexGroups[v] = connect[0]
			}
			delete(vertexByGroup, connect[1])
		case lIn:
			groupId = lG
			edgeToAdd = e[1]
			addToExistingGroup = true
		case rIn:
			groupId = rG
			edgeToAdd = e[0]
			addToExistingGroup = true
		default:
			vertexByGroup[lastUsedGroup] = []uint64{e[0], e[1]}

			vertexGroups[e[0]] = lastUsedGroup
			vertexGroups[e[1]] = lastUsedGroup

			lastUsedGroup++
		}
		mst = append(mst, e)

		if addToExistingGroup {
			vertexGroups[edgeToAdd] = groupId
			if _, ok := vertexByGroup[groupId]; ok {
				vertexByGroup[groupId] = append(vertexByGroup[groupId], edgeToAdd)
			} else {
				vertexByGroup[groupId] = []uint64{edgeToAdd}
			}
		}
	}

	return
}

// NewReversedGraph Returns a copy allocated in a new memory space of the
// sorted graph, but with the edges in the opposite way
func (gr *Graph) NewReversedGraph() (rev *Graph) {
	rev = &Graph{
		Vertices:    gr.Vertices,
		VertexEdges: make(map[uint64]map[uint64]uint64),
		Undirected:  false,
	}

	for v, e := range gr.VertexEdges {
		for d, w := range e {
			if _, ok := rev.VertexEdges[d]; ok {
				rev.VertexEdges[d][v] = w
			} else {
				rev.VertexEdges[d] = map[uint64]uint64{v: w}
			}
		}
	}

	return
}

// StronglyConnectedComponents Detects all the strongly connected components in
// a directed graph and returns a map with the vertices as a key and the group
// where the vertex belongs as value, and another slice of maps where each map
// is used as a set who groups the vertices by groups, the keys of the maps
// are vertex number
// The algorithm used is the Kosaraju-Sharir's algorithm:
// 	- http://en.wikipedia.org/wiki/Kosaraju%27s_algorithm
func (gr *Graph) StronglyConnectedComponents() (components map[uint64]uint64, compGroups []map[uint64]bool) {
	currentGroup := uint64(0)
	components = make(map[uint64]uint64)
	compGroups = []map[uint64]bool{}

	topologOrder, _ := gr.TopologicalOrder()
	rev := gr.NewReversedGraph()
	for _, v := range topologOrder {
		if _, in := components[v]; !in {
			group := rev.Dfs(v)
			compGroups = append(compGroups, make(map[uint64]bool))
			for uv := range group {
				if _, in := components[uv]; !in {
					compGroups[len(compGroups)-1][uv] = true
					components[uv] = currentGroup
				}
			}
			currentGroup++
		}
	}

	return
}

// TopologicalOrder Calculates  the topological on directed graphs, where every
// directed edge uv from vertex u to vertex v, u comes before v in the ordering
// In case of doesn't be possible calculate the topological order of the graph
// because of a cycle, or any other cause, the second parameter will return
// false
func (gr *Graph) TopologicalOrder() (order []uint64, success bool) {
	if gr.Undirected {
		return nil, false
	}

	verticesToUse := make(map[uint64]bool)
	for v := range gr.Vertices {
		verticesToUse[v] = true
	}

	var orig uint64
	orderPos := uint64(0)
	order = make([]uint64, len(gr.Vertices))
	group := make(map[uint64]bool)
	for len(verticesToUse) > 0 {
		for orig = range verticesToUse {
			break
		}

		orderAux := make([]uint64, len(gr.Vertices))
		pos := uint64(0)
		gr.dfs(orig, group, orderAux, &pos)

		for i := uint64(0); i < pos; i++ {
			orderPos++
			delete(verticesToUse, orderAux[i])
			order[uint64(len(order))-orderPos] = orderAux[i]
		}
	}

	// Check if we have any cycle after sort the vertices, we have a cycle
	// if any of the vertices have a edge to a previously visited vertex
	mapPos := make(map[uint64]uint64)
	for i := uint64(0); i < uint64(len(order)); i++ {
		mapPos[order[i]] = i
		for edge := range gr.VertexEdges[order[i]] {
			if _, ok := mapPos[edge]; ok {
				return order, false
			}
		}
	}
	success = true

	return
}

// IsBipartite Checks if a graph is bipartite from the given vertex, in case of
// a graph composed by multiple components, checks if the component where this
// vertex is located is bipartite
func (gr *Graph) IsBipartite(origin uint64) bool {
	colours := map[uint64]bool{origin: false}
	queue := []uint64{origin}

	for len(queue) > 0 {
		current := queue[0]
		queue = queue[1:]
		for v := range gr.VertexEdges[current] {
			if _, visited := colours[v]; !visited {
				colours[v] = !colours[current]
				queue = append(queue, v)
			} else {
				if colours[v] == colours[current] {
					return false
				}
			}
		}
	}

	return true
}

// Copy Returns a copy allocated in a different memory space of the graph
func (gr *Graph) Copy() (cp *Graph) {
	cp = &Graph{
		VertexEdges: make(map[uint64]map[uint64]uint64),
	}
	for k, e := range gr.VertexEdges {
		for d, w := range e {
			if _, ok := cp.VertexEdges[k]; ok {
				cp.VertexEdges[k][d] = w
			} else {
				cp.VertexEdges[k] = map[uint64]uint64{d: w}
			}
		}
	}

	return
}

// EulerianPath Calculates a path that starting on the "orig" vertex and ending
// on the "end" vertex walks through all the edges on the graph.
// The second returned parameter specifies if existst or not a Eulerian path
// on the graph
func (gr *Graph) EulerianPath(orig uint64, end uint64) (path []uint64, success bool) {
	// For an Eulerian Path all the vertices but the origin and ending
	// vertices has to have a even degree, we will check in
	// EulerianCycle the even degree of all the vertices, so now we only
	// check the orig and end vertices
	if orig != end && (len(gr.VertexEdges[orig])%2 == 0 || len(gr.VertexEdges[end])%2 == 0) {
		return nil, false
	}

	// Remove the connection between the origin and end vertices
	newGr := gr.Copy()
	delete(newGr.VertexEdges[orig], end)
	delete(newGr.VertexEdges[end], orig)
	path, success = newGr.EulerianCycle(orig)
	if !success {
		return nil, false
	}
	path = append(path, end)

	return
}

// HamiltonianPath Calculates a path that visits each vertex exactly once. A
// same origin and destination can be specified in order to calculate a
// Hamilton tour
// This is a NP-complete problem
func (gr *Graph) HamiltonianPath(orig uint64, dest uint64) (path []uint64, success bool) {
	visited := make(map[uint64]bool)
	if orig != dest {
		visited[orig] = true
	}
	path = []uint64{orig}

	return gr.hamiltonianPath(orig, &dest, visited, path)
}

func (gr *Graph) hamiltonianPath(orig uint64, dest *uint64, visited map[uint64]bool, path []uint64) ([]uint64, bool) {
	if len(visited) == len(gr.VertexEdges) {
		if path[len(path)-1] == *dest {
			return path, true
		}

		return nil, false
	}

	for tv := range gr.VertexEdges[orig] {
		if _, ok := visited[tv]; !ok && (*dest != tv || len(visited) == len(gr.VertexEdges)-1) {
			visited[tv] = true
			path = append(path, tv)
			if path, found := gr.hamiltonianPath(tv, dest, visited, path); found {
				return path, true
			}
			path = path[:len(path)-1]
			delete(visited, tv)
		}
	}

	return nil, false
}

// EulerianCycle Calculates a cycle that starting and ending on the "orig"
// vertex walks through all the edges on the graph.
// The second returned parameter specifies if existst or not a Eulerian cycle
// on the graph
func (gr *Graph) EulerianCycle(orig uint64) (tour []uint64, success bool) {
	// For an Eulerian cirtuit all the vertices has to have a even degree
	for _, e := range gr.VertexEdges {
		if len(e)%2 != 0 {
			return nil, false
		}
	}

	// Hierholzer's algorithm
	var currentVertex, nextVertex uint64

	tour = []uint64{}
	stack := []uint64{orig}
	unvisitedEdg := gr.Copy().VertexEdges
	for len(stack) > 0 {
		currentVertex = stack[len(stack)-1]
		// Get an arbitrary edge from the current vertex
		if len(unvisitedEdg[currentVertex]) > 0 {
			for nextVertex = range unvisitedEdg[currentVertex] {
				break
			}
			delete(unvisitedEdg[currentVertex], nextVertex)
			delete(unvisitedEdg[nextVertex], currentVertex)
			stack = append(stack, nextVertex)
		} else {
			tour = append(tour, stack[len(stack)-1])
			stack = stack[:len(stack)-1]
		}
	}

	return tour, true
}

// ConnectedComponents Returns a slice of maps where the keys are the
// vertices, each element on the slice is a set of interconnected vertices but
// without connection with any other vertex in any other returned set of
// vertices
func (gr *Graph) ConnectedComponents() (groups []map[uint64]bool) {
	var groupToUse uint64
	usedVertex := make(map[uint64]uint64)
	currentGroup := uint64(0)
	for v := range gr.VertexEdges {
		if _, used := usedVertex[v]; !used {
			group := make(map[uint64]bool)
			gr.dfs(v, group, nil, nil)
			found := false
		groupSearch:
			for k := range group {
				if g, used := usedVertex[k]; used {
					groupToUse = g
					found = true

					break groupSearch
				}
			}
			if !found {
				groupToUse = currentGroup
				currentGroup++
				groups = append(groups, make(map[uint64]bool))
			}

			for k := range group {
				usedVertex[k] = groupToUse
				groups[groupToUse][k] = true
			}
		}
	}

	return
}

// Bfs Calculates the shortest path from the origin vertex to all the connected
// vertices and returns the list of edges and distances
func (gr *Graph) Bfs(origin uint64) (edgeTo map[uint64]uint64, distTo map[uint64]uint64) {
	queue := []uint64{origin}
	edgeTo = map[uint64]uint64{origin: origin}
	distTo = map[uint64]uint64{origin: 0}

	for len(queue) > 0 {
		current := queue[0]
		queue = queue[1:]
		deep := distTo[current] + 1
		for v := range gr.VertexEdges[current] {
			if _, visited := distTo[v]; !visited {
				distTo[v] = deep
				edgeTo[v] = current
				queue = append(queue, v)
			}
		}
	}

	return
}

// Dfs Finds all vertices connected to the "origin" vertex  and returns them as
// an slice of vertices.
// This method uses Depth-first search algorithm:
//	- http://en.wikipedia.org/wiki/Depth-first_search
// The Tremaux's algorithm is used to perform this search:
//	- http://en.wikipedia.org/wiki/Maze_solving_algorithm#Tr.C3.A9maux.27s_algorithm
func (gr *Graph) Dfs(root uint64) (usedVertex map[uint64]bool) {
	usedVertex = make(map[uint64]bool)
	gr.dfs(root, usedVertex, nil, nil)

	return
}

func (gr *Graph) dfs(origin uint64, usedVertex map[uint64]bool, order []uint64, pos *uint64) {
	usedVertex[origin] = true
	for v := range gr.VertexEdges[origin] {
		if _, visited := usedVertex[v]; !visited {
			gr.dfs(v, usedVertex, order, pos)
		}
	}
	if order != nil {
		order[*pos] = origin
		*pos++
	}
}
