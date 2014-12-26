package graphs

//import "fmt"

type graph struct {
	vertexEdges map[uint64]map[uint64]bool
}

func GetUndirected(edges [][2]uint64) (ug *graph) {
	ug = &graph{
		vertexEdges: make(map[uint64]map[uint64]bool),
	}

	for _, edge := range edges {
		if _, ok := ug.vertexEdges[edge[0]]; ok {
			ug.vertexEdges[edge[0]][edge[1]] = true
		} else {
			ug.vertexEdges[edge[0]] = map[uint64]bool{edge[1]: true}
		}
		if _, ok := ug.vertexEdges[edge[1]]; ok {
			ug.vertexEdges[edge[1]][edge[0]] = true
		} else {
			ug.vertexEdges[edge[1]] = map[uint64]bool{edge[0]: true}
		}
	}

	return
}

func (gr *graph) IsBipartite() {
}

func (gr *graph) IsEulerian() {
}

func (gr *graph) GetConnectedComponents() (groups []map[uint64]bool) {
	usedVertex := make(map[uint64]bool)
	for v := range gr.vertexEdges {
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

// Calculates the shortest path from the origin vertex to all the connected
// vertices and returns the list of edges and distances
func (gr *graph) Bfs(origin uint64) (edgeTo map[uint64]uint64, distTo map[uint64]uint64) {
	queue := []uint64{origin}
	edgeTo = map[uint64]uint64{0: 0}
	distTo = map[uint64]uint64{0: 0}

	for len(queue) > 0 {
		current := queue[0]
		queue = queue[1:]
		deep := distTo[current] + 1
		for v := range gr.vertexEdges[current] {
			if _, visited := distTo[v]; !visited {
				distTo[v] = deep
				edgeTo[v] = current
				queue = append(queue, v)
			}
		}
	}

	return
}

// Finds all vertices connected to the "origin" vertex  and returns them as an
// slice of vertices.
// This method uses Depth-first search algorithm:
//	- http://en.wikipedia.org/wiki/Depth-first_search
// The Tremaux's algorithm is used to perform this search:
//	- http://en.wikipedia.org/wiki/Maze_solving_algorithm#Tr.C3.A9maux.27s_algorithm
func (gr *graph) Dfs(origin uint64) (usedVertex map[uint64]bool) {
	usedVertex = make(map[uint64]bool)
	gr.dfs(origin, usedVertex)

	return
}

func (gr *graph) dfs(origin uint64, usedVertex map[uint64]bool) {
	usedVertex[origin] = true
	for v := range gr.vertexEdges[origin] {
		if _, visited := usedVertex[v]; !visited {
			gr.dfs(v, usedVertex)
		}
	}
}
