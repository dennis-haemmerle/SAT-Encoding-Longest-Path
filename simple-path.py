import networkx as nx
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc
from pysat.solvers import Solver


def simple_path_of_length_k(G: nx.Graph, k: int) -> bool:
    if k < 0 or G.number_of_nodes() < k + 1:
        return False
    if k == 0:
        return G.number_of_nodes() > 0

    vpool = IDPool()
    cnf = CNF()

    # 1. Each position is occupied by exactly one node.
    for i in range(k + 1):
        lits = [vpool.id((v, i)) for v in G.nodes()]
        # block = CardEnc.equals(lits=lits, vpool=vpool, encoding=0)
        # cnf.extend(block.clauses)
        cnf.extend(CardEnc.atleast(lits=lits, vpool=vpool, encoding=0).clauses)
        cnf.extend(CardEnc.atmost(lits=lits, vpool=vpool, encoding=0).clauses)

    # 2. Each node appears at most once.
    for v in G.nodes():
        lits = [vpool.id((v, i)) for i in range(k + 1)]
        block = CardEnc.atmost(lits=lits, vpool=vpool, encoding=0)
        cnf.extend(block.clauses)

    '''
    # 3. Don't allow unconnected nodes between consecutive positions.
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
                    cnf.append([-vpool.id((u, i)), -vpool.id((v, i + 1))])
    '''
    # 3. Require that consecutive positions are connected by an edge.
    for i in range(k):
        for u in G.nodes():
            clause = [-vpool.id((u, i))] + [vpool.id((v, i + 1)) for v in G.neighbors(u)]
            cnf.append(clause)

    # solver = Solver(name="Cadical195")
    # for clause in cnf.clauses:
    #     solver.add_clause(clause)
    # return bool(solver.solve())
    with Solver(name="Cadical195", bootstrap_with=cnf.clauses) as solver:
        return bool(solver.solve())


def longest_simple_path_linear_search(G: nx.Graph):
    best = 0

    for k in range(1, G.number_of_nodes() + 1):
        if simple_path_of_length_k(G, k):
            best = k
        else:
            break

    return best


def longest_simple_path_binary_search(G: nx.Graph):
    low = 0
    high = G.number_of_nodes()
    best = 0

    while low <= high:
        mid = (low + high) // 2
        if simple_path_of_length_k(G, mid):
            best = mid
            low = mid + 1
        else:
            high = mid - 1

    return best


if __name__ == "__main__":
    G = nx.Graph()
    G.add_edges_from([
        (0, 1), (1, 4), (4, 5),
        (0, 2), (2, 6), (6, 7),
        (0, 3), (3, 8), (8, 9)
    ])
    assert simple_path_of_length_k(G, 6) is True
    assert simple_path_of_length_k(G, 7) is False
    assert longest_simple_path_linear_search(G) == 6
    assert longest_simple_path_binary_search(G) == 6

    G = nx.DiGraph()
    G.add_edges_from([
        (0, 1), (1, 4), (4, 5),
        (0, 2), (2, 6), (6, 7),
        (0, 3), (3, 8), (8, 9)
    ])
    assert simple_path_of_length_k(G, 3) is True
    assert simple_path_of_length_k(G, 4) is False
    assert longest_simple_path_linear_search(G) == 3
    assert longest_simple_path_binary_search(G) == 3
