import networkx as nx
from functools import lru_cache
from itertools import combinations
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


def simple_path_of_length_k(G: nx.Graph, k: int, start=None, end=None, only_in=None, only_out=None, leaves=None):
    def to_set(val):
        if val is None:
            return set()
        elif isinstance(val, (list, set)):
            return set(val)
        else:
            return {val}

    only_out = to_set(only_out)
    only_in = to_set(only_in)
    leaves = to_set(leaves)
    start = to_set(start)
    end = to_set(end)

    if k < 0 or k > G.number_of_edges():
        return None
    if k == 0:
        return [] if G.number_of_nodes() > 0 else None

    vpool = IDPool()
    cnf = CNF()

    allowed_positions = {}
    full = frozenset(range(k + 1))
    for v in G.nodes():
        if (len(start) == 1 and v in start) or v in only_out:
            allowed_positions[v] = {0}
        elif (len(end) == 1 and v in end) or v in only_in:
            allowed_positions[v] = {k}
        elif v in leaves:
            allowed_positions[v] = {0, k}
        else:
            allowed_positions[v] = full

    allowed_nodes = {i: [] for i in range(k + 1)}
    for v, positions in allowed_positions.items():
        for i in positions:
            allowed_nodes[i].append(v)

    # 1. Each position is occupied by exactly one node.
    for i in range(k + 1):
        lits = [vpool.id((v, i)) for v in allowed_nodes[i]]
        block = CardEnc.equals(lits=lits, bound=1, vpool=vpool, encoding=EncType.seqcounter)
        cnf.extend(block.clauses)

    # 2. Each node appears at most once.
    for v in G.nodes():
        lits = [vpool.id((v, i)) for i in allowed_positions[v]]
        block = CardEnc.atmost(lits=lits, bound=1, vpool=vpool, encoding=EncType.seqcounter)
        cnf.extend(block.clauses)

    # 3. Require that consecutive positions are connected by an edge.
    for u in G.nodes():
        for i in allowed_positions[u]:
            if i == k:
                continue
            clause = [-vpool.id((u, i))] + [vpool.id((v, i + 1)) for v in G.neighbors(u) if i + 1 in allowed_positions[v]]
            cnf.append(clause)
    '''# 3. Don't allow unconnected nodes between consecutive positions.
    if G.is_directed():
        def connected(u, v):
            return G.has_edge(u, v)
    else:
        def connected(u, v):
            return G.has_edge(u, v) or G.has_edge(v, u)

    for i in range(k):
        for u in G.nodes():
            for v in G.nodes():
                if u == v:
                    continue
                if not connected(u, v):
                    cnf.append([-vpool.id((u, i)), -vpool.id((v, i + 1))])'''

    # Set optional start/end nodes
    if start:
        cnf.append([vpool.id((s, 0)) for s in start])
    if end:
        cnf.append([vpool.id((e, k)) for e in end])

    with Solver(name="Cadical195", bootstrap_with=cnf.clauses) as solver:
        if solver.solve():
            model = set(solver.get_model())  # type: ignore
            assignment = [v for i in range(k + 1) for v in G.nodes() if vpool.id((v, i)) in model]
            return assignment
        return None


def longest_simple_path_linear_search(G: nx.Graph, start=None, end=None, only_in=None, only_out=None, leaves=None):
    longest_path = []

    for k in range(1, G.number_of_nodes()):
        path = simple_path_of_length_k(G, k, start, end, only_in, only_out, leaves)
        if path is not None:
            longest_path = path
        else:
            break

    return longest_path


def longest_simple_path_binary_search(G: nx.Graph, start=None, end=None, only_in=None, only_out=None, leaves=None):
    longest_path = []
    low = 0
    high = G.number_of_nodes() - 1

    while low <= high:
        mid = (low + high) // 2
        path = simple_path_of_length_k(G, mid, start, end, only_in, only_out, leaves)

        if path is not None:
            longest_path = path
            low = mid + 1
        else:
            high = mid - 1

    return longest_path


def longest_simple_path_components(C: nx.Graph):
    tree = C.graph["bridge_tree"]
    blocks = C.graph["bridge_components"]

    @lru_cache(None)
    def best_between(block_id: int, start=None, end=None):
        block = blocks[block_id]

        if block.number_of_nodes() == 1:
            return block.nodes()
        else:
            return longest_simple_path_binary_search(
                G=block,
                start=start,
                end=end,
                only_in=block.graph.get("only_in_nodes", []),
                only_out=block.graph.get("only_out_nodes", []),
                leaves=block.graph.get("leaves", [])
            )

    if tree.number_of_nodes() == 1:
        return best_between(block_id=0)

    longest_path = []

    # Only paths inside a block
    sorted_blocks = sorted(tree.nodes, key=lambda id: len(blocks[id]), reverse=True)
    for block_id in sorted_blocks:
        if len(blocks[block_id]) <= len(longest_path):
            break  # all subsequent blocks are even smaller
        path = best_between(block_id)
        if len(path) > len(longest_path):
            longest_path = path

    # Try out all paths in bridge_tree
    leaves = [n for n, d in tree.degree() if d <= 1]
    for s, t in combinations(leaves, 2):  # all possible start/end blocks
        path = nx.shortest_path(tree, s, t)
        path_nodes = []  # nodes of the current path
        path_splits = []  # paths that had to be split, because they can't forward a simple path (enter_node == exit_node)

        # First block
        first_edge = tree[path[0]][path[1]]
        first_node = first_edge["attach"][path[0]]
        path_nodes.extend(best_between(path[0], end=first_node))

        # Middle blocks
        for i in range(1, len(path) - 1):
            prev_block = path[i - 1]
            current_block = path[i]
            next_block = path[i + 1]

            enter_edge = tree[prev_block][current_block]
            exit_edge = tree[current_block][next_block]

            enter_node = enter_edge["attach"][current_block]
            exit_node = exit_edge["attach"][current_block]

            if blocks[current_block].number_of_nodes() > 1 and enter_node == exit_node:
                # Split the path if there is an invalid block inbetween
                path_nodes.extend(best_between(current_block, start=enter_node))
                path_splits.append(path_nodes)
                path_nodes = []
                path_nodes.extend(best_between(current_block, end=enter_node))
            else:
                path_nodes.extend(best_between(current_block, enter_node, exit_node))

        # Last block
        last_edge = tree[path[-2]][path[-1]]
        last_node = last_edge["attach"][path[-1]]
        path_nodes.extend(best_between(path[-1], start=last_node))

        if path_splits:
            path_splits.append(path_nodes)
            longest_path = max(path_splits, key=len)
        else:
            if len(path_nodes) > len(longest_path):
                longest_path = path_nodes

    return longest_path


def longest_simple_path(G: nx.Graph):
    longest_path = []

    for C in G.graph.get("connected_components", [G]):
        path = longest_simple_path_components(C)
        if len(path) > len(longest_path):
            longest_path = path

    return longest_path
