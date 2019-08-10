from geopy.distance import great_circle
from copy import deepcopy


def tsp(data):
    # build a graph
    G = build_graph(data)
    # print("Graph: ", G)

    # build a minimum spanning tree
    MSTree = minimum_spanning_tree(G)
    MSTree_init = deepcopy(MSTree)
    # print("MSTree: ", MSTree)

    # find odd vertexes
    odd_vertexes = find_odd_vertexes(MSTree)
    odd_vertexes_init = deepcopy(odd_vertexes)
    # print("Odd vertexes in MSTree: ", odd_vertexes)

    # add minimum weight matching edges to MST
    new_added_matching = minimum_weight_matching(MSTree, G, odd_vertexes)
    united_MSTree_perfect_matching = deepcopy(MSTree)
    # print("Minimum weight matching: ", MSTree)

    # find an eulerian tour
    eulerian_tour = find_eulerian_tour(MSTree, G)

    # print("Eulerian tour: ", eulerian_tour)

    current = eulerian_tour[0]
    path = [current]
    visited = [False] * len(eulerian_tour)

    length = 0

    for v in eulerian_tour[1:]:
        if not visited[v]:
            path.append(v)
            visited[v] = True

            length += G[current][v]
            current = v

    # path.append(path[0])

    # print("Result path: ", path)
    # print("Result length of the path: ", length)

    return G, MSTree_init, odd_vertexes_init, new_added_matching, united_MSTree_perfect_matching, eulerian_tour, length, path



def get_length(x1, y1, x2, y2, name='great_circle'):
    '''
    x1: lat1
    y1: lng1
    x2: lat2
    y2: lng2
    '''
    if name == 'great_circle':
        return great_circle((x1,y1), (x2,y2)).km
    else:
        return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1 / 2)


def build_graph(data):
    graph = {}
    for this in range(len(data)):
        for another_point in range(len(data)):
            if this != another_point:
                if this not in graph:
                    graph[this] = {}

                graph[this][another_point] = get_length(data[this][0], 
                                                        data[this][1], 
                                                        data[another_point][0],
                                                        data[another_point][1],
                                                        name='great_circle')

    return graph


class UnionFind:
    def __init__(self):
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __iter__(self):
        return iter(self.parents)

    def union(self, *objects):
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r], r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest


def minimum_spanning_tree(G):
    tree = []
    subtrees = UnionFind()
    for W, u, v in sorted((G[u][v], u, v) for u in G for v in G[u]):
        if subtrees[u] != subtrees[v]:
            tree.append((u, v, W))
            subtrees.union(u, v)

    return tree


def find_odd_vertexes(MST):
    tmp_g = {}
    vertexes = []
    for edge in MST:
        if edge[0] not in tmp_g:
            tmp_g[edge[0]] = 0

        if edge[1] not in tmp_g:
            tmp_g[edge[1]] = 0

        tmp_g[edge[0]] += 1
        tmp_g[edge[1]] += 1

    for vertex in tmp_g:
        if tmp_g[vertex] % 2 == 1:
            vertexes.append(vertex)

    return vertexes


def minimum_weight_matching(MST, G, odd_vert):
    import random
    random.shuffle(odd_vert)

    new_added = []
    while odd_vert:
        v = odd_vert.pop()
        length = float("inf")
        u = 1
        closest = 0
        for u in odd_vert:
            if v != u and G[v][u] < length:
                length = G[v][u]
                closest = u

        MST.append((v, closest, length))
        new_added.append((v, closest, length))
        odd_vert.remove(closest)
    return new_added

def find_eulerian_tour(MatchedMSTree, G):
    # find neigbours
    neighbours = {}
    for edge in MatchedMSTree:
        if edge[0] not in neighbours:
            neighbours[edge[0]] = []

        if edge[1] not in neighbours:
            neighbours[edge[1]] = []

        neighbours[edge[0]].append(edge[1])
        neighbours[edge[1]].append(edge[0])

    # print("Neighbours: ", neighbours)

    # finds the hamiltonian circuit
    start_vertex = MatchedMSTree[0][0]
    EP = [neighbours[start_vertex][0]]

    while len(MatchedMSTree) > 0:
        for i, v in enumerate(EP):
            if len(neighbours[v]) > 0:
                break

        while len(neighbours[v]) > 0:
            w = neighbours[v][0]

            remove_edge_from_matchedMST(MatchedMSTree, v, w)

            del neighbours[v][(neighbours[v].index(w))]
            del neighbours[w][(neighbours[w].index(v))]

            i += 1
            EP.insert(i, w)

            v = w

    return EP


def remove_edge_from_matchedMST(MatchedMST, v1, v2):

    for i, item in enumerate(MatchedMST):
        if (item[0] == v2 and item[1] == v1) or (item[0] == v1 and item[1] == v2):
            del MatchedMST[i]

    return MatchedMST


def Euler_Tour(multigraph):
	""" Uses Fleury's algorithm to find the Euler Tour of the MultiGraph.
	"""
	tour = []
	temp_graph = nx.MultiGraph()
	graph_nodes = nx.nodes(multigraph)
	current_node = graph_nodes[0]
	tour.append(current_node)
	while nx.number_of_edges(multigraph) > 0: 	
		for edge in multigraph.edges(current_node):
			temp_graph = copy.deepcopy(multigraph)
			temp_graph.remove_edge(edge[0], edge[1], key=None)
			if nx.is_connected(temp_graph):
				tour.append(edge[1])
				current_node = edge[1]
				multigraph.remove_edge(edge[0], edge[1], key=None)
				break
		else:
			tour.append(edge[1])
			current_node = edge[1]
			multigraph.remove_edge(edge[0], edge[1], key=None)
			multigraph.remove_nodes_from(nx.isolates(multigraph))
	return tour


def shortcut_Euler_Tour(tour):
	"""Find's the shortcut of the Euler Tour to obtain the Approximation.
	"""
	Tour = []
	for vertex in tour:
		if vertex not in Tour:
			Tour.append(vertex)
	Tour.append(tour[0])
	return Tour
