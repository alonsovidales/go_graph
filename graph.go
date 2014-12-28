package graphs

//import "fmt"

// Graph Data strcture used to represent a graph, the VertexEdges var is a map
// where each key is a vertex, and the value a map where the keys are the
// vertices that can be reached from the main key vertex, the internal map is
// used as a set
type Graph struct {
	VertexEdges map[uint64]map[uint64]bool
}

// GetUndirected Returns an undirected graph containing the specified edges
func GetUndirected(edges [][2]uint64) (ug *Graph) {
	ug = &Graph{
		VertexEdges: make(map[uint64]map[uint64]bool),
	}

	for _, edge := range edges {
		if _, ok := ug.VertexEdges[edge[0]]; ok {
			ug.VertexEdges[edge[0]][edge[1]] = true
		} else {
			ug.VertexEdges[edge[0]] = map[uint64]bool{edge[1]: true}
		}
		if _, ok := ug.VertexEdges[edge[1]]; ok {
			ug.VertexEdges[edge[1]][edge[0]] = true
		} else {
			ug.VertexEdges[edge[1]] = map[uint64]bool{edge[0]: true}
		}
	}

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
		VertexEdges: make(map[uint64]map[uint64]bool),
	}
	for k, e := range gr.VertexEdges {
		for d := range e {
			if _, ok := cp.VertexEdges[k]; ok {
				cp.VertexEdges[k][d] = true
			} else {
				cp.VertexEdges[k] = map[uint64]bool{d: true}
			}
		}
	}

	return
}

// GetEulerianPath Calculates a path that starting on the "orig" vertex and ending
// on the "end" vertex walks through all the edges on the graph.
// The second returned parameter specifies if existst or not a Eulerian path
// on the graph
func (gr *Graph) GetEulerianPath(orig uint64, end uint64) (path []uint64, success bool) {
	// For an Eulerian Path all the vertices but the origin and ending
	// vertices has to have a even degree, we will check in
	// GetEulerianCycle the even degree of all the vertices, so now we only
	// check the orig and end vertices
	if orig != end && (len(gr.VertexEdges[orig])%2 == 0 || len(gr.VertexEdges[end])%2 == 0) {
		return nil, false
	}

	// Remove the connection between the origin and end vertices
	newGr := gr.Copy()
	delete(newGr.VertexEdges[orig], end)
	delete(newGr.VertexEdges[end], orig)
	path, success = newGr.GetEulerianCycle(orig)
	if !success {
		return nil, false
	}
	path = append(path, end)

	return
}

// GetEulerianCycle Calculates a cycle that starting and ending on the "orig"
// vertex walks through all the edges on the graph.
// The second returned parameter specifies if existst or not a Eulerian cycle
// on the graph
func (gr *Graph) GetEulerianCycle(orig uint64) (tour []uint64, success bool) {
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

// GetConnectedComponents Returns a slice of maps where the keys are the
// vertices, each element on the slice is a set of interconnected vertices but
// without connection with any other vertex in any other returned set of
// vertices
func (gr *Graph) GetConnectedComponents() (groups []map[uint64]bool) {
	usedVertex := make(map[uint64]bool)
	for v := range gr.VertexEdges {
		if _, used := usedVertex[v]; !used {
			group := make(map[uint64]bool)
			gr.dfs(v, group)
			groups = append(groups, group)
			for k := range group {
				usedVertex[k] = true
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
func (gr *Graph) Dfs(origin uint64) (usedVertex map[uint64]bool) {
	usedVertex = make(map[uint64]bool)
	gr.dfs(origin, usedVertex)

	return
}

func (gr *Graph) dfs(origin uint64, usedVertex map[uint64]bool) {
	usedVertex[origin] = true
	for v := range gr.VertexEdges[origin] {
		if _, visited := usedVertex[v]; !visited {
			gr.dfs(v, usedVertex)
		}
	}
}
