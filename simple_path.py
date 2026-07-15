import networkx as nx
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
            assignment = []
            for i in range(k + 1):
                for v in G.nodes():
                    if vpool.id((v, i)) in model:
                        assignment.append(v)
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
