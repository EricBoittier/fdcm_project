import networkx as nx
from anytree import Node, RenderTree
import numpy as np


def schedule_from_2drmsd(filename):
    # load data
    two_d_rmsd_data = load_adjacency_matrix(filename)
    # Sparsify matrix based on RMSD cut off
    two_d_rmsd_data_sparse = sparsify_matrix(two_d_rmsd_data)
    # Graph and visiting list
    G, visited = find_paths_from_adjacency(two_d_rmsd_data_sparse)
    # Using anytree package for Tree data struct.
    anytree = list_to_anytree(visited)
    # Using the tree and the entire graph, return the schedule of scan jobs
    scan_neighbours_schedule = scan_neighbours_schedule_from_anytree_graph(anytree, G)
    return scan_neighbours_schedule


def print_anytree(anytree):
    for pre, fill, node in RenderTree(anytree):
        print(f"{pre}{node.name}")


def list_to_anytree(lst):
    root_name = lst[0][0]
    root_node = Node(root_name)
    nodes = {root_name: root_node}  # keeping a dict of the nodes
    for branch in lst:
        assert branch[0] == root_name
        for parent_name, node_name in zip(branch, branch[1:]):
            node = nodes.setdefault(node_name, Node(node_name))
            parent_node = nodes[parent_name]
            if node.parent is not None:
                assert node.parent.name == parent_name
            else:
                node.parent = parent_node
    return root_node


def load_adjacency_matrix(two_d_rmsd):
    two_d_rmsd_load = np.loadtxt(two_d_rmsd)
    #  Remove the header line
    two_d_rmsd_data = two_d_rmsd_load[:, 1:]
    return two_d_rmsd_data


def sparsify_matrix(two_d_rmsd_data):
    #  Calculate cut off value
    #  i.e. the largest min. neighbour distance of all nodes
    CUTOFF = np.where(two_d_rmsd_data == 0, 100, two_d_rmsd_data).min(axis=1).max()
    #  Sparsify the graph by removing edges where RMSD > cut off
    two_d_rmsd_data_sparse = np.where(two_d_rmsd_data < CUTOFF, two_d_rmsd_data, 0)
    return two_d_rmsd_data_sparse


def find_paths_from_adjacency(two_d_rmsd_data_sparse):
    #  Create a weighted graph
    G = nx.from_numpy_matrix(two_d_rmsd_data_sparse)
    #  Find shortest paths
    shortest_paths = nx.single_source_dijkstra_path(G, 0)
    #  Find unique shortest paths (some paths may overlap)
    paths = list(dict(sorted(shortest_paths.items(), key=lambda item: len(item))).values())
    paths.reverse()
    path_visits = [paths[0]]
    for p in paths[1:]:
        if p not in [x[:-1] for x in path_visits]:
            path_visits.append(p)

    lengths = [len(x) for x in path_visits]
    path_lengths = {x: [] for x in range(np.min(lengths), np.max(lengths) + 1)}

    for i, path in enumerate(path_visits):
        path_lengths[len(path)].append(path)

    visited = []
    max_depth = max(lengths)
    min_depth = min(lengths)
    for i in range(max_depth, min_depth - 1, -1):
        for path in path_lengths[i]:
            if i == max_depth:
                visited.append(path)
            else:
                __sublist__ = False
                for v in visited:
                    if len(v) > i:
                        if v[i - 1] == path[-1]:
                            __sublist__ = True
                if not __sublist__:
                    visited.append(path)

    return G, visited


def scan_neighbours_schedule_from_anytree_graph(anytree, G):
    """
    Scan recursively through the tree, starting from the root
    """
    order = []

    def scan_children(node):
        if len(node.children) > 0:
            children = node.children
            for c in children:
                order.append([node.name, c.name])
                if len(c.children) > 0:
                    scan_children(c)
        else:
            pass

    scan_children(anytree)

    print_anytree(anytree)

    scan_neighbours_schedule = []
    visited = []
    for i in order:
        start, end = tuple(i)
        if len(visited) == 0:
            print("prev.:", order[0])
            previous = order[0][0]
        else:
            previous = [x for x in visited if x == start]
            print("prev.:", [x for x in visited if x == start])

        visited.append(previous)
        neighbours = [x for x in G.edges(start) if x[1] != end]
        schedule = [tuple(x) for x in order if x[0] == end]
        scan_neighbours_schedule.append(((previous, start, end), neighbours, schedule))

    return scan_neighbours_schedule

