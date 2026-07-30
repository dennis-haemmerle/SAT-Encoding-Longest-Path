import networkx as nx
from itertools import combinations
from collections import defaultdict
from pysat.formula import CNF, WCNF, IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
from pysat.examples.rc2 import RC2

from optimizations import optimize
from incremental_simple_path import IncrementalSimplePathEncoder


def read_from_file(filename: str):
    """
    Reads a graph from a file.

    Supported formats:
        p graph <num_nodes> <num_edges>
        p digraph <num_nodes> <num_edges>

        u v
        ...
    """

    G = None

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("c"):
                continue

            if line.startswith("p"):
                parts = line.split()

                if len(parts) != 4:
                    raise ValueError(f"Invalid header: {line}")

                _, graph_type, num_nodes, _ = parts
                num_nodes = int(num_nodes)

                if graph_type == "graph":
                    G = nx.Graph()
                elif graph_type == "digraph":
                    G = nx.DiGraph()
                else:
                    raise ValueError(f"Unknown graph type '{graph_type}'.")

                G.add_nodes_from(range(1, num_nodes + 1))
                continue

            if G is None:
                raise ValueError("Header must appear before any edge definitions.")

            u, v = map(int, line.split())
            G.add_edge(u, v)

    if G is None:
        raise ValueError("No header found.")

    return G


def to_set(val: list[int] | set[int] | int | None) -> set[int]:
    if val is None:
        return set()
    elif isinstance(val, (list, set)):
        return set(val)
    else:
        return {val}


def simple_path_of_length_k(G: nx.Graph, k: int, start=None, end=None, symmetry=None):
    if k < 0 or k > G.number_of_edges():
        return None
    if k == 0:
        return [] if G.number_of_nodes() > 0 else None

    start = to_set(start)
    end = to_set(end)

    vpool = IDPool()
    cnf = CNF()

    allowed_positions = {}
    for v in G.nodes():
        allowed = set(range(k + 1))

        if start:
            if (len(start) == 1 and v in start):
                allowed = {0}
            elif v not in start:
                allowed -= {0}

        if end:
            if (len(end) == 1 and v in end):
                allowed = {k}
            elif v not in end:
                allowed -= {k}

        allowed_positions[v] = allowed

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
        cnf.append([vpool.id((v, 0)) for v in start])
    if end:
        cnf.append([vpool.id((v, k)) for v in end])

    # Symmetry breaking
    if symmetry is not None:
        orbit_groups = symmetry.get("orbit_groups", {})
        for orbit in orbit_groups.values():
            if len(orbit) <= 1:
                continue

            # Only the representative of each orbit is allowed to be the start node.
            valid_starters = [v for v in orbit if 0 in allowed_positions[v]]
            if not valid_starters:
                continue

            starters = [v for v in valid_starters if v in start]
            representative = min(starters) if starters else min(valid_starters)

            for v in orbit:
                if v != representative:
                    cnf.append([-vpool.id((v, 0))])

    with Solver(name="Cadical195", bootstrap_with=cnf.clauses) as solver:
        if solver.solve():
            model = set(solver.get_model())  # type: ignore
            assignment = [v for i in range(k + 1) for v in G.nodes() if vpool.id((v, i)) in model]
            return assignment
        return None


def simple_path_of_length_atleast_k(G: nx.Graph, k: int, start=None, end=None, symmetry=None):
    if k < 0 or k > G.number_of_edges():
        return None
    if k == 0:
        return [] if G.number_of_nodes() > 0 else None

    start = to_set(start)
    end = to_set(end)

    vpool = IDPool()
    cnf = CNF()

    allowed_positions = {}
    for v in G.nodes():
        allowed = set(range(G.number_of_nodes()))

        if start:
            if (len(start) == 1 and v in start):
                allowed = {0}
            elif v not in start:
                allowed -= {0}

        allowed_positions[v] = allowed

    allowed_nodes = {i: [] for i in range(G.number_of_nodes())}
    for v, positions in allowed_positions.items():
        for i in positions:
            allowed_nodes[i].append(v)

    # 1. Each position is occupied by atmost one node.
    for i in range(G.number_of_nodes()):
        lits = [vpool.id((v, i)) for v in G.nodes()]
        cnf.extend(CardEnc.atmost(lits=lits, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)

    #    The first k positions are occupied by atleast one node.
    for i in range(k + 1):
        lits = [vpool.id((v, i)) for v in G.nodes()]
        cnf.extend(CardEnc.atleast(lits=lits, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 2. Each node appears at most once.
    for v in G.nodes():
        lits = [vpool.id((v, i)) for i in range(G.number_of_nodes())]
        cnf.extend(CardEnc.atmost(lits=lits, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 3. For each position i > k: if there is a node, then there has to be a node on position i - 1
    for i in range(k + 1, G.number_of_nodes()):
        prev_pos = [vpool.id((v, i - 1)) for v in G.nodes()]
        for u in G.nodes():
            clause = [-vpool.id((u, i))] + prev_pos
            cnf.append(clause)

    """ # 4. Require that consecutive positions are connected by an edge.
    if G.is_directed():
        predecessors = {v: list(G.predecessors(v)) for v in G.nodes()}  # type: ignore
    else:
        predecessors = {v: list(G.neighbors(v)) for v in G.nodes()}

    for v in G.nodes():
        for i in range(1, G.number_of_nodes()):
            clause = [-vpool.id((v, i))] + [vpool.id((u, i - 1)) for u in predecessors[v]]  # type: ignore
            cnf.append(clause) """
    # 4. Require that consecutive positions are connected by an edge.
    for i in range(G.number_of_nodes()):
        end_id = vpool.id(("end", i)) if i >= k else None

        if i < G.number_of_nodes() - 1:
            for u in G.nodes():
                clause = [-vpool.id((u, i))] + [vpool.id((v, i + 1)) for v in G.neighbors(u)] + ([end_id] if i >= k else [])
                cnf.append(clause)
        else:
            if end_id is not None:
                cnf.append([end_id])

    #   If the path ends at position i, there can be no more nodes at i+1
    for i in range(k, G.number_of_nodes() - 1):
        end_id = vpool.id(("end", i))
        for v in G.nodes():
            cnf.append([-end_id, -vpool.id((v, i + 1))])

    # Set optional start/end nodes
    if start:
        cnf.append([vpool.id((s, 0)) for s in start])
    if end:
        # Atleast one endpoint has to be at a valid position i >= k.
        cnf.append([vpool.id((e, i)) for e in end for i in range(k, G.number_of_nodes())])
        # If endpoint e is at position i, then there can not be any node at position i + 1.
        for i in range(k, G.number_of_nodes() - 1):
            for e in end:
                for v in G.nodes():
                    cnf.append([-vpool.id((e, i)), -vpool.id((v, i + 1))])

    # Symmetry breaking
    if symmetry is not None:
        orbit_groups = symmetry.get("orbit_groups", {})
        for orbit in orbit_groups.values():
            if len(orbit) <= 1:
                continue

            # Only the representative of each orbit is allowed to be the start node.
            valid_starters = [v for v in orbit if 0 in allowed_positions[v]]
            if not valid_starters:
                continue

            starters = [v for v in valid_starters if v in start]
            representative = min(starters) if starters else min(valid_starters)

            for v in orbit:
                if v != representative:
                    cnf.append([-vpool.id((v, 0))])

    with Solver(name="Cadical195", bootstrap_with=cnf.clauses) as solver:
        if solver.solve():
            model = set(solver.get_model())  # type: ignore
            assignment = [v for i in range(k + 1) for v in G.nodes() if vpool.id((v, i)) in model]
            return assignment
        return None


def simple_path_of_length_atleast_k_2(G: nx.Graph, k: int, start=None, end=None, symmetry=None):
    if k < 0 or k > G.number_of_edges():
        return None
    if k == 0:
        return [] if G.number_of_nodes() > 0 else None

    start = to_set(start)
    end = to_set(end)

    vpool = IDPool()
    cnf = CNF()

    allowed_positions = {}
    for v in G.nodes():
        allowed = set(range(G.number_of_nodes()))

        if start:
            if (len(start) == 1 and v in start):
                allowed = {0}
            elif v not in start:
                allowed -= {0}

        allowed_positions[v] = allowed

    allowed_nodes = {i: [] for i in range(G.number_of_nodes())}
    for v, positions in allowed_positions.items():
        for i in positions:
            allowed_nodes[i].append(v)

    # 0. Atleast k acitve variables
    for i in range(k + 1):
        cnf.append([vpool.id(i)])

    #    If position i+1 is active, then position i is also active
    for i in range(k, G.number_of_nodes() - 1):
        cnf.append([-vpool.id(i + 1), vpool.id(i)])

    # 1. Each position is occupied by exactly one node.
    for i in range(G.number_of_nodes()):
        lits = [vpool.id((v, i)) for v in G.nodes()]
        for clause in CardEnc.equals(lits=lits, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses:
            cnf.append([-vpool.id(i)] + clause)

        # If any node is located at position i, position i must be active.
        for v in G.nodes():
            cnf.append([-vpool.id((v, i)), vpool.id(i)])

    # 2. Each node appears at most once.
    for v in G.nodes():
        lits = [vpool.id((v, i)) for i in range(G.number_of_nodes())]
        cnf.extend(CardEnc.atmost(lits=lits, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 3. Require that consecutive positions are connected by an edge.
    for u in G.nodes():
        for i in range(G.number_of_nodes() - 1):
            clause = [-vpool.id((u, i)), -vpool.id(i + 1)] + [vpool.id((v, i + 1)) for v in G.neighbors(u)]
            cnf.append(clause)

    # Set optional start/end nodes
    if start:
        cnf.append([vpool.id((s, 0)) for s in start])
    if end:
        for i in range(k, G.number_of_nodes()):
            # if position i is active and position i+1 not, then atleast one endpoint has to be at position i
            clause = [-vpool.id(i)]
            if i < G.number_of_nodes() - 1:
                clause.append(vpool.id(i + 1))

            cnf.append(clause + [vpool.id((e, i)) for e in end])

    # Symmetry breaking
    if symmetry is not None:
        orbit_groups = symmetry.get("orbit_groups", {})
        for orbit in orbit_groups.values():
            if len(orbit) <= 1:
                continue

            # Only the representative of each orbit is allowed to be the start node.
            valid_starters = [v for v in orbit if 0 in allowed_positions[v]]
            if not valid_starters:
                continue

            starters = [v for v in valid_starters if v in start]
            representative = min(starters) if starters else min(valid_starters)

            for v in orbit:
                if v != representative:
                    cnf.append([-vpool.id((v, 0))])

    with Solver(name="Cadical195", bootstrap_with=cnf.clauses) as solver:
        if solver.solve():
            model = set(solver.get_model())  # type: ignore
            assignment = [v for i in range(k + 1) for v in G.nodes() if vpool.id((v, i)) in model]
            return assignment
        return None


def simple_path_of_length_k_edge_encoding(G: nx.Graph, k: int, start=None, end=None, symmetry=None, atleast_k=False):
    if k < 0 or k > G.number_of_edges():
        return None
    if k == 0:
        return [] if G.number_of_nodes() > 0 else None

    start = to_set(start)
    end = to_set(end)

    vpool = IDPool()
    cnf = CNF()

    def edge_var(e):
        u, v = e
        if G.is_directed():
            return vpool.id((u, v))
        return vpool.id((min(u, v), max(u, v)))

    # 1. Exactly k edges
    lits = [edge_var(e) for e in G.edges()]
    if atleast_k:
        cnf.extend(CardEnc.atleast(lits=lits, bound=k, vpool=vpool, encoding=EncType.seqcounter).clauses)
    else:
        cnf.extend(CardEnc.equals(lits=lits, bound=k, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 2. Exactly k+1 nodes
    lits = [vpool.id(v) for v in G.nodes()]
    if atleast_k:
        cnf.extend(CardEnc.atleast(lits=lits, bound=k + 1, vpool=vpool, encoding=EncType.seqcounter).clauses)
    else:
        cnf.extend(CardEnc.equals(lits=lits, bound=k + 1, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 3. Atleast one incident edge for each used node
    for v in G.nodes():
        if G.is_directed():
            incident = ([edge_var(e) for e in G.in_edges(v)] + [edge_var(e) for e in G.out_edges(v)])  # type: ignore
        else:
            incident = [edge_var(e) for e in G.edges(v)]

        for e in incident:
            cnf.append([-e, vpool.id(v)])  # e -> v
        cnf.append([-vpool.id(v)] + incident)  # v -> atleast one incident

    # 4. Node degree atmost 2
    for v in G.nodes():
        if G.is_directed():
            incoming = [edge_var(e) for e in G.in_edges(v)]  # type: ignore
            outgoing = [edge_var(e) for e in G.out_edges(v)]  # type: ignore

            cnf.extend(CardEnc.atmost(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)
            cnf.extend(CardEnc.atmost(lits=outgoing, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)
        else:
            incident = [edge_var(e) for e in G.edges(v)]
            cnf.extend(CardEnc.atmost(lits=incident, bound=2, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 5. Exactly two endpoints with degree 1
    if G.is_directed():
        start_vars = []
        end_vars = []

        for v in G.nodes():
            incoming = [edge_var(e) for e in G.in_edges(v)]  # type: ignore
            outgoing = [edge_var(e) for e in G.out_edges(v)]  # type: ignore

            # start_var = True <=> (indeg(v) = 0 ∧ outdeg(v) = 1)
            start_var = vpool.id(("start", v))
            start_vars.append(start_var)
            # end_var = True <=> (indeg(v) = 1 ∧ outdeg(v) = 0)
            end_var = vpool.id(("end", v))
            end_vars.append(end_var)

            # start_var -> outdeg(v) = 1
            if len(outgoing) == 0:
                cnf.append([-start_var])
            else:
                exactly_one_out = CardEnc.equals(lits=outgoing, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_out.clauses:
                    cnf.append([-start_var] + clause)

            # start_var -> indeg(v) = 0
            for e in incoming:
                cnf.append([-start_var, -e])

            # end_var -> indeg(v) = 1
            if len(incoming) == 0:
                cnf.append([-end_var])
            else:
                exactly_one_in = CardEnc.equals(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_in.clauses:
                    cnf.append([-end_var] + clause)

            # end_var -> outdeg(v) = 0
            for e in outgoing:
                cnf.append([-end_var, -e])

            out_deg_1 = vpool.id(("out_deg_1", v))
            if not outgoing:
                cnf.append([-out_deg_1])
            else:
                # out_deg_1 -> outdeg(v) = 1
                exactly_one_out = CardEnc.equals(lits=outgoing, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_out.clauses:
                    cnf.append([-out_deg_1] + clause)

                # outdeg(v) = 1 -> out_deg_1
                for i, e in enumerate(outgoing):
                    # exactly one outgoing -> out_deg_1
                    cnf.append([-e] + outgoing[:i] + outgoing[i + 1:] + [out_deg_1])

                for i in range(len(outgoing)):
                    for j in range(i + 1, len(outgoing)):
                        # atmost one outgoing -> out_deg_1
                        cnf.append([-outgoing[i], -outgoing[j], -out_deg_1])

            # (indeg(v) = 0 ∧ outdeg(v) = 1) -> start_var
            if incoming:
                cnf.append([-out_deg_1] + incoming + [start_var])
            else:
                cnf.append([-out_deg_1, start_var])

            in_deg_1 = vpool.id(("in_deg_1", v))
            if not incoming:
                cnf.append([-in_deg_1])
            else:
                # in_deg_1 -> indeg(v) = 1
                exactly_one_in = CardEnc.equals(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_in.clauses:
                    cnf.append([-in_deg_1] + clause)

                # indeg(v) = 1 -> in_deg_1
                for i, e in enumerate(incoming):
                    # exactly one incoming -> in_deg_1
                    cnf.append([-e] + incoming[:i] + incoming[i + 1:] + [in_deg_1])

                for i in range(len(incoming)):
                    for j in range(i + 1, len(incoming)):
                        # atmost one incoming -> in_deg_1
                        cnf.append([-incoming[i], -incoming[j], -in_deg_1])

            # (indeg(v) = 1 ∧ outdeg(v) = 0) -> end_var
            if outgoing:
                cnf.append([-in_deg_1] + outgoing + [end_var])
            else:
                cnf.append([-in_deg_1, end_var])

            # cnf.append([-start_var, -end_var])

            # (start_var v end_var) -> v
            cnf.append([-start_var, vpool.id(v)])
            cnf.append([-end_var, vpool.id(v)])

        if start:
            cnf.append([vpool.id(("start", v)) for v in start])
        if end:
            cnf.append([vpool.id(("end", v)) for v in end])

        # Symmetry breaking
        if symmetry is not None:
            orbit_groups = symmetry.get("orbit_groups", {})
            for orbit in orbit_groups.values():
                if len(orbit) <= 1:
                    continue

                # Only the representative of each orbit is allowed to be the start node.
                starters = [v for v in orbit if v in start] if start else list(orbit)
                if not starters:
                    continue

                representative = min(starters)

                for v in orbit:
                    if v != representative:
                        cnf.append([-vpool.id(("start", v))])
    else:
        endpoint_vars = []

        for v in G.nodes():
            incident = [edge_var(e) for e in G.edges(v)]

            # end_var = True <=> deg(v) = 1
            end_var = vpool.id(("end", v))
            endpoint_vars.append(end_var)

            if not incident:
                cnf.append([-end_var])
                continue

            # end_var -> deg(v) = 1
            exactly_one = CardEnc.equals(lits=incident, bound=1, vpool=vpool, encoding=EncType.seqcounter)
            for clause in exactly_one.clauses:
                cnf.append([-end_var] + clause)

            # deg(v) = 1 -> end_var
            for i, e in enumerate(incident):
                # atleast one incident -> end_var
                cnf.append([-e] + incident[:i] + incident[i + 1:] + [end_var])

            for i in range(len(incident)):
                for j in range(i + 1, len(incident)):
                    # atmost one incident -> end_var
                    cnf.append([-incident[i], -incident[j], -end_var])

            # end_var -> v
            cnf.append([-end_var, vpool.id(v)])

        if start:
            cnf.append([vpool.id(("end", v)) for v in start])
        if end:
            cnf.append([vpool.id(("end", v)) for v in end])

        # Symmetry breaking
        if symmetry is not None:
            orbit_groups = symmetry.get("orbit_groups", {})
            for orbit in orbit_groups.values():
                if len(orbit) <= 1:
                    continue

                # Only the representative of each orbit is allowed to be the start node.
                starters = [v for v in orbit if v in start] if start else list(orbit)
                if not starters:
                    continue

                representative = min(starters)

                for v in orbit:
                    if v != representative:
                        cnf.append([-vpool.id(("end", v))])

    if G.is_directed():
        cnf.extend(CardEnc.equals(lits=start_vars, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)
        cnf.extend(CardEnc.equals(lits=end_vars, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)
    else:
        cnf.extend(CardEnc.equals(lits=endpoint_vars, bound=2, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 6. Acyclicity via reachability
    reach = {}
    for u in G.nodes():
        for v in G.nodes():
            # Variables: r_{u,v} (u reaches v)
            reach[(u, v)] = vpool.id(("reach", (u, v)))

    if not G.is_directed():
        dir_vars = {}
        for u, v in G.edges():
            dir_vars[(u, v)] = vpool.id(("dir", (u, v)))
            dir_vars[(v, u)] = vpool.id(("dir", (v, u)))

    for e in G.edges():
        u, v = e
        if G.is_directed():
            # e_{u,v} -> r_{u,v}
            cnf.append([-edge_var(e), reach[(u, v)]])
        else:
            # e_{u,v} -> (dir_{u,v} v dir_{u,v})
            cnf.append([-edge_var(e), dir_vars[(u, v)], dir_vars[(v, u)]])
            # e_{u,v} -> not (dir_{u,v} ∧ dir_{u,v})
            cnf.append([-edge_var(e), -dir_vars[(u, v)], -dir_vars[(v, u)]])

            # (e ∧ dir(u,v)) -> reach(u,v)
            cnf.append([-edge_var(e), -dir_vars[(u, v)], reach[(u, v)]])
            # (e ∧ dir(v,u)) -> reach(v,u)
            cnf.append([-edge_var(e), -dir_vars[(v, u)], reach[(v, u)]])

        for x in G.nodes():
            if G.is_directed():
                # e_{u,v} ∧ r_{x,u} -> r_{x,v}
                cnf.append([-edge_var(e), -reach[(x, u)], reach[(x, v)]])
            else:
                # (e_{u,v} ∧ dir_{u,v} ∧ r_{x,u}) -> r_{x,v}
                cnf.append([-edge_var(e), -dir_vars[(u, v)], -reach[(x, u)], reach[(x, v)]])
                # (e_{u,v} ∧ dir_{v,u} ∧ r_{x,v}) -> r_{x,u}
                cnf.append([-edge_var(e), -dir_vars[(v, u)], -reach[(x, v)], reach[(x, u)]])

    """ # Ensure all active edges have a consistent orientation along the path
    if not G.is_directed():
        for e in G.edges():
            u, v = e
            for x in G.neighbors(v):
                if x == u:
                    continue
                # (e_{u,v} ∧ dir_{u,v} ∧ e_{v,x}) -> dir_{v,x}
                cnf.append([-edge_var(e), -dir_vars[(u, v)], -edge_var((v, x)), dir_vars[(v, x)]])
                pass
            for x in G.neighbors(u):
                if x == v:
                    continue
                # (e_{u,v} ∧ dir_{v,u} ∧ e_{u,x}) -> dir_{u,x}
                cnf.append([-edge_var(e), -dir_vars[(v, u)], -edge_var((u, x)), dir_vars[(u, x)]]) """
    # Ensure all active edges have a consistent orientation along the path
    if not G.is_directed():
        for v in G.nodes():
            incoming = [dir_vars[(u, v)] for u in G.neighbors(v)]
            cnf.extend(CardEnc.atmost(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)

    for v in G.nodes():
        cnf.append([-reach[(v, v)]])

    """ # 6. Subtour Elimination with DFJ
    for r in range(2, G.number_of_nodes()):
        for subset_nodes in combinations(G.nodes(), r):
            subset_nodes = set(subset_nodes)

            subset_edges = []
            for u, v in G.edges():
                if u in subset_nodes and v in subset_nodes:
                    subset_edges.append(edge_var((u, v)))

            diff = len(subset_edges) - len(subset_nodes)
            if diff >= 0:
                max_true = len(subset_edges) - (diff + 1)
                cnf.extend(CardEnc.atmost(lits=subset_edges, bound=max_true, vpool=vpool, encoding=EncType.seqcounter).clauses) """

    """ # 6. Subtour Elimination with DFJ (Lazy Cut Loop)
    def selected_subgraph(model):
        H = nx.Graph() if not G.is_directed() else nx.DiGraph()

        for v in G.nodes():
            if vpool.id(v) in model:
                H.add_node(v)

        for e in G.edges():
            if edge_var(e) in model:
                u, v = e
                H.add_edge(u, v)

        return H

    def add_dfj_cut(model, cycle_nodes):
        cycle_edges = [-edge_var((u, v)) for u, v in G.subgraph(cycle_nodes).edges() if edge_var((u, v)) in model]

        if cycle_edges:
            cnf.append(cycle_edges) """

    with Solver(name="Cadical195", bootstrap_with=cnf.clauses) as solver:
        """ # DFJ Lazy Cut Loop
        while solver.solve():
            model = solver.get_model()
            H = selected_subgraph(model)

            comps = list(nx.weakly_connected_components(H) if H.is_directed() else nx.connected_components(H))  # type: ignore
            cycles = []
            for comp in comps:
                used_edges = H.subgraph(comp).number_of_edges()
                used_nodes = H.subgraph(comp).number_of_nodes()
                if used_edges > 0 and used_nodes < k + 1 and used_edges >= used_nodes:
                    cycles.append(comp)

            if not cycles:
                model = set(solver.get_model())  # type: ignore
                assignment = [e for e in G.edges() if edge_var(e) in model]
                return assignment

            for cycle in cycles:
                add_dfj_cut(model, cycle)
                solver.add_clause(cnf.clauses[-1])

        return None """
        if solver.solve():
            model = set(solver.get_model())  # type: ignore
            assignment = [e for e in G.edges() if edge_var(e) in model]
            return extract_path(assignment, G.is_directed())
        return None


def longest_simple_path_edge_encoding_maxsat(G: nx.Graph):
    if G.number_of_nodes() <= 0:
        return []

    vpool = IDPool()
    wcnf = WCNF()

    def edge_var(e):
        u, v = e
        if G.is_directed():
            return vpool.id((u, v))
        return vpool.id((min(u, v), max(u, v)))

    # 1. Soft clauses: maximize number of selected edges
    for e in G.edges():
        wcnf.append([edge_var(e)], weight=1)

    # 4. Node degree atmost 2
    if G.is_directed():
        for v in G.nodes():
            incoming = [edge_var(e) for e in G.in_edges(v)]  # type: ignore
            outgoing = [edge_var(e) for e in G.out_edges(v)]  # type: ignore

            if incoming:
                wcnf.extend(CardEnc.atmost(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)
            if outgoing:
                wcnf.extend(CardEnc.atmost(lits=outgoing, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)
    else:
        for v in G.nodes():
            incident = [edge_var(e) for e in G.edges(v)]
            if incident:
                wcnf.extend(CardEnc.atmost(lits=incident, bound=2, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 5. Exactly two endpoints with degree 1.
    if G.is_directed():
        start_vars = []
        end_vars = []

        for v in G.nodes():
            incoming = [edge_var(e) for e in G.in_edges(v)]  # type: ignore
            outgoing = [edge_var(e) for e in G.out_edges(v)]  # type: ignore

            # start_var = True <=> (indeg(v) = 0 ∧ outdeg(v) = 1)
            start_var = vpool.id(("start", v))
            start_vars.append(start_var)
            # end_var = True <=> (indeg(v) = 1 ∧ outdeg(v) = 0)
            end_var = vpool.id(("end", v))
            end_vars.append(end_var)

            # start_var -> outdeg(v) = 1
            if len(outgoing) == 0:
                wcnf.append([-start_var])
            else:
                exactly_one_out = CardEnc.equals(lits=outgoing, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_out.clauses:
                    wcnf.append([-start_var] + clause)

            # start_var -> indeg(v) = 0
            for e in incoming:
                wcnf.append([-start_var, -e])

            # end_var -> indeg(v) = 1
            if len(incoming) == 0:
                wcnf.append([-end_var])
            else:
                exactly_one_in = CardEnc.equals(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_in.clauses:
                    wcnf.append([-end_var] + clause)

            # end_var -> outdeg(v) = 0
            for e in outgoing:
                wcnf.append([-end_var, -e])

            out_deg_1 = vpool.id(("out_deg_1", v))
            if not outgoing:
                wcnf.append([-out_deg_1])
            else:
                # out_deg_1 -> outdeg(v) = 1
                exactly_one_out = CardEnc.equals(lits=outgoing, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_out.clauses:
                    wcnf.append([-out_deg_1] + clause)

                # outdeg(v) = 1 -> out_deg_1
                for i, e in enumerate(outgoing):
                    # exactly one outgoing -> out_deg_1
                    wcnf.append([-e] + outgoing[:i] + outgoing[i + 1:] + [out_deg_1])

                for i in range(len(outgoing)):
                    for j in range(i + 1, len(outgoing)):
                        # atmost one outgoing -> out_deg_1
                        wcnf.append([-outgoing[i], -outgoing[j], -out_deg_1])

            # (indeg(v) = 0 ∧ outdeg(v) = 1) -> start_var
            if incoming:
                wcnf.append([-out_deg_1] + incoming + [start_var])
            else:
                wcnf.append([-out_deg_1, start_var])

            in_deg_1 = vpool.id(("in_deg_1", v))
            if not incoming:
                wcnf.append([-in_deg_1])
            else:
                # in_deg_1 -> indeg(v) = 1
                exactly_one_in = CardEnc.equals(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter)
                for clause in exactly_one_in.clauses:
                    wcnf.append([-in_deg_1] + clause)

                # indeg(v) = 1 -> in_deg_1
                for i, e in enumerate(incoming):
                    # exactly one incoming -> in_deg_1
                    wcnf.append([-e] + incoming[:i] + incoming[i + 1:] + [in_deg_1])

                for i in range(len(incoming)):
                    for j in range(i + 1, len(incoming)):
                        # atmost one incoming -> in_deg_1
                        wcnf.append([-incoming[i], -incoming[j], -in_deg_1])

            # (indeg(v) = 1 ∧ outdeg(v) = 0) -> end_var
            if outgoing:
                wcnf.append([-in_deg_1] + outgoing + [end_var])
            else:
                wcnf.append([-in_deg_1, end_var])

            # cnf.append([-start_var, -end_var])

            # (start_var v end_var) -> v
            # wcnf.append([-start_var, vpool.id(v)])
            # wcnf.append([-end_var, vpool.id(v)])
    else:
        endpoint_vars = []

        for v in G.nodes():
            incident = [edge_var(e) for e in G.edges(v)]

            # end_var = True <=> deg(v) = 1
            end_var = vpool.id(("end", v))
            # endpoint_vars.append(end_var)
            endpoint_vars.append(end_var)

            if not incident:
                wcnf.append([-end_var])
                continue

            # end_var -> deg(v) = 1
            exactly_one = CardEnc.equals(lits=incident, bound=1, vpool=vpool, encoding=EncType.seqcounter)
            for clause in exactly_one.clauses:
                wcnf.append([-end_var] + clause)

            # deg(v) = 1 -> end_var
            for i, e in enumerate(incident):
                # atleast one incident -> end_var
                wcnf.append([-e] + incident[:i] + incident[i + 1:] + [end_var])

            for i in range(len(incident)):
                for j in range(i + 1, len(incident)):
                    # atmost one incident -> end_var
                    wcnf.append([-incident[i], -incident[j], -end_var])

            # end_var -> v
            # wcnf.append([-end_var, vpool.id(v)])

    if G.is_directed():
        block = CardEnc.equals(lits=start_vars, bound=1, vpool=vpool, encoding=EncType.seqcounter)
        wcnf.extend(block.clauses)
        block = CardEnc.equals(lits=end_vars, bound=1, vpool=vpool, encoding=EncType.seqcounter)
        wcnf.extend(block.clauses)
    else:
        if len(endpoint_vars) >= 2:
            block = CardEnc.equals(lits=endpoint_vars, bound=2, vpool=vpool, encoding=EncType.seqcounter)
            wcnf.extend(block.clauses)

    """ # 6. Acyclity via reachability
    reach = {}
    for u in G.nodes():
        for v in G.nodes():
            # Variables: r_{u,v} (u reaches v)
            reach[(u, v)] = vpool.id(("reach", (u, v)))

    if not G.is_directed():
        dir_vars = {}
        for u, v in G.edges():
            dir_vars[(u, v)] = vpool.id(("dir", (u, v)))
            dir_vars[(v, u)] = vpool.id(("dir", (v, u)))

    for e in G.edges():
        u, v = e
        if G.is_directed():
            # e_{u,v} -> r_{u,v}
            wcnf.append([-edge_var(e), reach[(u, v)]])
        else:
            # e_{u,v} -> (dir_{u,v} v dir_{u,v})
            wcnf.append([-edge_var(e), dir_vars[(u, v)], dir_vars[(v, u)]])
            # e_{u,v} -> not (dir_{u,v} ∧ dir_{u,v})
            wcnf.append([-edge_var(e), -dir_vars[(u, v)], -dir_vars[(v, u)]])

            # (e ∧ dir(u,v)) -> reach(u,v)
            wcnf.append([-edge_var(e), -dir_vars[(u, v)], reach[(u, v)]])
            # (e ∧ dir(v,u)) -> reach(v,u)
            wcnf.append([-edge_var(e), -dir_vars[(v, u)], reach[(v, u)]])

        for x in G.nodes():
            if G.is_directed():
                # e_{u,v} ∧ r_{x,u} -> r_{x,v}
                wcnf.append([-edge_var(e), -reach[(x, u)], reach[(x, v)]])
            else:
                # (e_{u,v} ∧ dir_{u,v} ∧ r_{x,u}) -> r_{x,v}
                wcnf.append([-edge_var(e), -dir_vars[(u, v)], -reach[(x, u)], reach[(x, v)]])
                # (e_{u,v} ∧ dir_{v,u} ∧ r_{x,v}) -> r_{x,u}
                wcnf.append([-edge_var(e), -dir_vars[(v, u)], -reach[(x, v)], reach[(x, u)]])

        if not G.is_directed():
            for v in G.nodes():
                incoming = [dir_vars[(u, v)] for u in G.neighbors(v)]
                wcnf.extend(CardEnc.atmost(lits=incoming, bound=1, vpool=vpool, encoding=EncType.seqcounter).clauses)

    for v in G.nodes():
        wcnf.append([-reach[(v, v)]]) """

    # 6. Subtour Elimination with DFJ (Lazy Cut Loop)
    def selected_subgraph(model):
        H = nx.Graph() if not G.is_directed() else nx.DiGraph()

        for e in G.edges():
            if edge_var(e) in model:
                u, v = e
                H.add_edge(u, v)

        return H

    def add_dfj_cut(model, cycle_nodes):
        cycle_edges = [-edge_var((u, v)) for u, v in G.subgraph(cycle_nodes).edges() if edge_var((u, v)) in model]

        if cycle_edges:
            rc2.add_clause(cycle_edges)

    with RC2(wcnf) as rc2:
        # DFJ Lazy Cut Loop
        while True:
            model = rc2.compute()
            if model is None:
                return []

            H = selected_subgraph(model)

            comps = list(nx.weakly_connected_components(H) if H.is_directed() else nx.connected_components(H))  # type: ignore
            cycles = []
            for comp in comps:
                used_edges = H.subgraph(comp).number_of_edges()
                used_nodes = H.subgraph(comp).number_of_nodes()
                if used_edges > 0 and used_edges >= used_nodes:
                    cycles.append(comp)

            if not cycles:
                # return list(H.edges())
                return extract_path(list(H.edges()), G.is_directed())

            for cycle in cycles:
                add_dfj_cut(model, cycle)
        """ model = rc2.compute()
        if model:
            assignment = [e for e in G.edges() if edge_var(e) in model]
            return extract_path(assignment, G.is_directed())
        return [] """


def extract_path(edges: list[tuple], is_directed: bool) -> list:
    adj = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    total_degree = defaultdict(int)

    for u, v in edges:
        adj[u].append(v)
        out_degree[u] += 1
        in_degree[v] += 1
        total_degree[u] += 1
        total_degree[v] += 1

        if not is_directed:
            adj[v].append(u)

    if is_directed:
        for node in total_degree:
            if out_degree[node] - in_degree[node] == 1:
                start_node = node
                break
    else:
        for node, degree in total_degree.items():
            if degree == 1:
                start_node = node
                break

    path = [start_node]
    current = start_node
    visited_edges = set()

    for _ in range(len(edges)):
        next_node = None
        for next in adj[current]:
            edge_id = (current, next) if is_directed else tuple(sorted((current, next)))

            if edge_id not in visited_edges:
                visited_edges.add(edge_id)
                next_node = next
                break

        if next_node is None:
            break

        path.append(next_node)
        current = next_node

    return path


def longest_simple_path_linear_search(G: nx.Graph, start=None, end=None, symmetry=None, edge_encoding=False, atleast_k=False, atleast_k_variant=False):
    longest_path = []

    for k in range(1, G.number_of_nodes()):
        if edge_encoding:
            path = simple_path_of_length_k_edge_encoding(G, k, start, end, symmetry, atleast_k)
        else:
            if atleast_k:
                if atleast_k_variant:
                    path = simple_path_of_length_atleast_k_2(G, k, start, end, symmetry)
                else:
                    path = simple_path_of_length_atleast_k(G, k, start, end, symmetry)
            else:
                path = simple_path_of_length_k(G, k, start, end, symmetry)

        if path is not None:
            longest_path = path
        else:
            break

    return longest_path


def longest_simple_path_linear_search_top_down(G: nx.Graph, start=None, end=None, symmetry=None):
    for k in range(G.number_of_nodes() - 1, 0, -1):
        path = simple_path_of_length_k(G, k, start, end, symmetry)

        if path is not None:
            return path

    return []


def longest_simple_path_binary_search(G: nx.Graph, start=None, end=None, symmetry=None, edge_encoding=False, atleast_k=False, atleast_k_variant=False):
    longest_path = []
    low = 0
    high = G.number_of_nodes() - 1

    while low <= high:
        mid = (low + high) // 2
        if edge_encoding:
            path = simple_path_of_length_k_edge_encoding(G, mid, start, end, symmetry, atleast_k)
        else:
            if atleast_k:
                if atleast_k_variant:
                    path = simple_path_of_length_atleast_k_2(G, mid, start, end, symmetry)
                else:
                    path = simple_path_of_length_atleast_k(G, mid, start, end, symmetry)
            else:
                path = simple_path_of_length_k(G, mid, start, end, symmetry)

        if path is not None:
            longest_path = path
            low = mid + 1
        else:
            high = mid - 1

    return longest_path


def longest_simple_path_components(C: nx.Graph, dp: bool = True, encoding: int = 0):
    def longest_simple_path(subgraph: nx.Graph, reversed_subgraph: nx.Graph, enter_node=None, exit_node=None, symmetry=dict(), encoding: int = 0):
        """
        encoding [default: 0]:
            0 = position encoding (exactly k) - linear search (top-down)
            1 = position encoding (atleast k) variant 1 - linear search (bottom-up)
            2 = position encoding (atleast k) variant 2 - linear search (bottom-up)
            3 = position encoding (atleast k) variant 1 - binary search
            4 = position encoding (atleast k) variant 2 - binary search
            5 = edge encoding (atleast k) - linear search (bottom-up)
            6 = edge encoding (atleast k) - binary search
            7 = incremental position encoding
        """
        match encoding:
            case 1:
                if subgraph.is_directed() and not enter_node and exit_node:
                    return list(reversed(longest_simple_path_linear_search(
                        G=reversed_subgraph,
                        start=exit_node,
                        symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                        edge_encoding=False,
                        atleast_k=True,
                        atleast_k_variant=False
                    )))
                return longest_simple_path_linear_search(
                    G=subgraph,
                    start=enter_node,
                    end=exit_node,
                    # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                    edge_encoding=False,
                    atleast_k=True,
                    atleast_k_variant=False
                )
            case 2:
                if subgraph.is_directed() and not enter_node and exit_node:
                    return list(reversed(longest_simple_path_linear_search(
                        G=reversed_subgraph,
                        start=exit_node,
                        symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                        edge_encoding=False,
                        atleast_k=True,
                        atleast_k_variant=True
                    )))
                return longest_simple_path_linear_search(
                    G=subgraph,
                    start=enter_node,
                    end=exit_node,
                    # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                    edge_encoding=False,
                    atleast_k=True,
                    atleast_k_variant=True
                )
            case 3:
                if subgraph.is_directed() and not enter_node and exit_node:
                    return list(reversed(longest_simple_path_binary_search(
                        G=reversed_subgraph,
                        start=exit_node,
                        symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                        edge_encoding=False,
                        atleast_k=True,
                        atleast_k_variant=False
                    )))
                return longest_simple_path_binary_search(
                    G=subgraph,
                    start=enter_node,
                    end=exit_node,
                    # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                    edge_encoding=False,
                    atleast_k=True,
                    atleast_k_variant=False
                )
            case 4:
                if subgraph.is_directed() and not enter_node and exit_node:
                    return list(reversed(longest_simple_path_binary_search(
                        G=reversed_subgraph,
                        start=exit_node,
                        symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                        edge_encoding=False,
                        atleast_k=True,
                        atleast_k_variant=True
                    )))
                return longest_simple_path_binary_search(
                    G=subgraph,
                    start=enter_node,
                    end=exit_node,
                    # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                    edge_encoding=False,
                    atleast_k=True,
                    atleast_k_variant=True
                )
            case 5:
                return longest_simple_path_linear_search(
                    G=subgraph,
                    start=enter_node,
                    end=exit_node,
                    # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                    edge_encoding=True,
                    atleast_k=True
                )
            case 6:
                return longest_simple_path_binary_search(
                    G=subgraph,
                    start=enter_node,
                    end=exit_node,
                    # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                    edge_encoding=True,
                    atleast_k=True,
                )
            case 7:
                if subgraph.is_directed() and not enter_node and exit_node:
                    encoder = IncrementalSimplePathEncoder(reversed_subgraph)
                    path = list(reversed(encoder.longest_simple_path(
                        start=exit_node,
                        symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
                    )))
                else:
                    encoder = IncrementalSimplePathEncoder(subgraph)
                    path = encoder.longest_simple_path(
                        start=enter_node,
                        end=exit_node,
                        # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
                    )
                encoder.delete()
                return path
            case _:
                if subgraph.is_directed() and not enter_node and exit_node:
                    return list(reversed(longest_simple_path_linear_search_top_down(
                        G=reversed_subgraph,
                        start=exit_node,
                        symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
                    )))
                return longest_simple_path_linear_search_top_down(
                    G=subgraph,
                    start=enter_node,
                    end=exit_node,
                    # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
                )

    if C.is_directed():
        dag = C.graph["condensation_dag"]

        if dag.number_of_nodes() == 0:
            return []
        if dag.number_of_nodes() == 1:
            scc = next(iter(dag.nodes))
            subgraph = dag.nodes[scc]["subgraph"]
            return longest_simple_path_binary_search(subgraph)

        if not dp:
            # ==================================================
            # Brute Force Method
            # ==================================================

            longest_path = []

            # Find the longest path by checking paths between all pairs of blocks
            for start_scc, end_scc in combinations(dag.nodes, 2):
                # Check both directions for a path in topological order
                for src, dst in [(start_scc, end_scc), (end_scc, start_scc)]:
                    if not nx.has_path(dag, src, dst):
                        continue

                    # Try all simple paths in the condensations DAG
                    for dag_path in nx.all_simple_paths(dag, src, dst):
                        current_path = []
                        valid_dag_path = True
                        scc_connections = {}

                        # Collect possible connecting edges between consecutive SCCs
                        for i in range(len(dag_path) - 1):
                            u_scc = dag_path[i]
                            v_scc = dag_path[i + 1]
                            possible_edges = dag.edges[u_scc, v_scc]["original_edges"]

                            if not possible_edges:
                                valid_dag_path = False
                                break

                            scc_connections[i] = possible_edges

                        if not valid_dag_path:
                            continue

                        # Process SCCs along the DAG path
                        next_enter_node = None
                        for i, scc in enumerate(dag_path):
                            subgraph = dag.nodes[scc]["subgraph"]
                            symmetry = dag.nodes[scc]["symmetry"]
                            reversed_subgraph = subgraph.reverse()

                            enter_node = next_enter_node if i > 0 else None

                            # Determine exit nodes that connect current SCC to the next one
                            if i < len(dag_path) - 1:
                                valid_edges = [edge for edge in scc_connections[i] if edge[0] not in current_path and edge[1] not in current_path and (subgraph.number_of_nodes() == 1 or enter_node is None or edge[0] not in enter_node)]
                                if not valid_edges:
                                    valid_edges = [edge for edge in scc_connections[i] if edge[0] in enter_node and edge[1] not in current_path]
                                    if not valid_edges:
                                        valid_dag_path = False
                                        break
                                exit_node = [edge[0] for edge in valid_edges]
                            else:
                                exit_node = None

                            # Solve longest path inside current SCC
                            if subgraph.number_of_nodes() == 1:
                                scc_path = list(subgraph.nodes())
                            else:
                                if enter_node is not None and exit_node is not None:
                                    if enter_node[0] == exit_node[0]:
                                        scc_path = [enter_node[0]]
                                    else:
                                        # scc_path = longest_simple_path_linear_search_top_down(G=subgraph, start=enter_node, end=exit_node)
                                        scc_path = longest_simple_path(
                                            subgraph=subgraph,
                                            reversed_subgraph=reversed_subgraph,
                                            enter_node=enter_node,
                                            exit_node=exit_node,
                                            symmetry=symmetry,
                                            encoding=encoding
                                        )
                                else:
                                    # scc_path = longest_simple_path_binary_search(G=subgraph, start=enter_node, end=exit_node)
                                    scc_path = longest_simple_path(
                                        subgraph=subgraph,
                                        reversed_subgraph=reversed_subgraph,
                                        enter_node=enter_node,
                                        exit_node=exit_node,
                                        symmetry=symmetry,
                                        encoding=encoding
                                    )

                            if not scc_path:
                                current_path = []
                                break

                            current_path.extend(scc_path)

                            # Determine entry nodes for next SCC based on the exit node chosen by the sat solver
                            if i < len(dag_path) - 1:
                                chosen_exit = scc_path[-1]
                                corresponding_enter = [edge[1] for edge in valid_edges if edge[0] == chosen_exit]

                                if not corresponding_enter:
                                    current_path = []
                                    break

                                next_enter_node = corresponding_enter

                        if len(current_path) > len(longest_path):
                            longest_path = current_path

            return longest_path
        else:
            # ==================================================
            # Dynamic Programming with topological order
            # ==================================================

            DP = {}

            # Store all enter nodes for each scc
            enter_nodes = {scc: set() for scc in dag.nodes()}
            for u_scc, v_scc in dag.edges():
                for exit_node, next_enter_node in dag.edges[u_scc, v_scc]["original_edges"]:
                    enter_nodes[v_scc].add(next_enter_node)

            # Process SCCs in reverse topological order
            # Case: SCC is an intermediate or ending component of the global path
            for scc in reversed(list(nx.topological_sort(dag))):
                subgraph = dag.nodes[scc]["subgraph"]
                symmetry = dag.nodes[scc]["symmetry"]
                reversed_subgraph = subgraph.reverse()
                local_paths = {}

                # For each SCC compute the optimal path starting from every possible enter node
                for enter_node in enter_nodes[scc]:

                    if subgraph.number_of_nodes() == 1:
                        longest_path_from_enter_node = list(subgraph.nodes())
                    else:
                        # longest_path_from_enter_node = longest_simple_path_binary_search(G=subgraph, start=enter_node)
                        longest_path_from_enter_node = longest_simple_path(
                            subgraph=subgraph,
                            reversed_subgraph=reversed_subgraph,
                            enter_node=enter_node,
                            symmetry=symmetry,
                            encoding=encoding
                        )

                    # Explore transitions to successor SCCs in the DAG
                    for next_scc in dag.successors(scc):
                        edges_to_next = dag.edges[scc, next_scc]["original_edges"]

                        # Compute optimal path within the current SCC from an enter node to an exit node
                        for exit_node, next_enter_node in edges_to_next:
                            if enter_node == exit_node:
                                scc_path = [enter_node]
                            else:
                                if (enter_node, exit_node) not in local_paths:
                                    # local_paths[(enter_node, exit_node)] = longest_simple_path_linear_search_top_down(G=subgraph, start=enter_node, end=exit_node)
                                    local_paths[(enter_node, exit_node)] = longest_simple_path(
                                        subgraph=subgraph,
                                        reversed_subgraph=reversed_subgraph,
                                        enter_node=enter_node,
                                        exit_node=exit_node,
                                        symmetry=symmetry,
                                        encoding=encoding
                                    )
                                scc_path = local_paths[(enter_node, exit_node)]

                            # Extend path using previously computed DP results continuing from the next SCCs entry node
                            rest_path = DP.get((next_scc, next_enter_node), [])
                            total_path = scc_path + rest_path

                            if len(total_path) > len(longest_path_from_enter_node):
                                longest_path_from_enter_node = total_path

                    DP[(scc, enter_node)] = longest_path_from_enter_node

                # Case: SCC is the starting component of the global path
                longest_path_without_enter_node = []
                if len(enter_nodes[scc]) < subgraph.number_of_nodes():
                    if subgraph.number_of_nodes() == 1:
                        longest_path_without_enter_node = list(subgraph.nodes())
                    else:
                        # Find absolute the longest path within the scc
                        # longest_path_without_enter_node = longest_simple_path_binary_search(G=subgraph)
                        longest_path_without_enter_node = longest_simple_path(
                            subgraph=subgraph,
                            reversed_subgraph=reversed_subgraph,
                            symmetry=symmetry,
                            encoding=encoding
                        )
                        # paths that start with an enter node are already covered above
                        if longest_path_without_enter_node[0] in enter_nodes[scc]:
                            longest_path_without_enter_node = []

                    for next_scc in dag.successors(scc):
                        edges_to_next = dag.edges[scc, next_scc]["original_edges"]

                        # Find the longest path to all exit nodes
                        for exit_node, next_enter_node in edges_to_next:
                            if subgraph.number_of_nodes() == 1:
                                scc_path = list(subgraph.nodes())
                            else:
                                if exit_node not in local_paths:
                                    local_paths[exit_node] = list(reversed(longest_simple_path_binary_search(G=reversed_subgraph, start=exit_node)))
                                    local_paths[exit_node] = longest_simple_path(
                                        subgraph=subgraph,
                                        reversed_subgraph=reversed_subgraph,
                                        exit_node=exit_node,
                                        symmetry=symmetry,
                                        encoding=encoding
                                    )
                                scc_path = local_paths[exit_node]

                            if scc_path[0] not in enter_nodes[scc]:
                                rest_path = DP.get((next_scc, next_enter_node), [])
                                total_path = scc_path + rest_path

                                if len(total_path) > len(longest_path_without_enter_node):
                                    longest_path_without_enter_node = total_path

                    # Update DP table if a valid path starting from a non-enter node was found
                    if longest_path_without_enter_node:
                        if (scc, longest_path_without_enter_node[0]) not in DP:
                            DP[(scc, longest_path_without_enter_node[0])] = longest_path_without_enter_node
                        else:
                            DP[(scc, longest_path_without_enter_node[0])] = max(DP[(scc, longest_path_without_enter_node[0])], longest_path_without_enter_node, key=len)

            return max(DP.values(), key=len, default=[])
    else:
        block_cut_tree = C.graph["block_cut_tree"]
        blocks = [n for n, attr in block_cut_tree.nodes(data=True) if attr["type"] == "block"]

        if len(blocks) == 0:
            return []
        if len(blocks) == 1:
            subgraph = block_cut_tree.nodes[blocks[0]]["subgraph"]
            symmetry = block_cut_tree.nodes[blocks[0]]["symmetry"]
            """ return longest_simple_path_binary_search(
                G=subgraph,
                # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
            ) """
            return longest_simple_path(
                subgraph=subgraph,
                reversed_subgraph=subgraph,
                symmetry=symmetry,
                encoding=encoding
            )

        longest_path = []

        # Find the longest path by checking paths between all pairs of blocks
        for start_block, end_block in combinations(blocks, 2):
            # The path between two blocks in a tree is unique
            tree_path = nx.shortest_path(block_cut_tree, start_block, end_block)

            current_path = []
            for i, tree_node in enumerate(tree_path):
                if block_cut_tree.nodes[tree_node]["type"] == "block":
                    subgraph = block_cut_tree.nodes[tree_node]["subgraph"]
                    symmetry = block_cut_tree.nodes[tree_node]["symmetry"]

                    # The enter and exit nodes for a block are the adjacent cut nodes in the tree path
                    enter_node = tree_path[i - 1] if i > 0 else None
                    exit_node = tree_path[i + 1] if i < len(tree_path) - 1 else None

                    if enter_node is not None and exit_node is not None and subgraph.number_of_nodes() > 2:
                        """ block_path = longest_simple_path_linear_search_top_down(
                            G=subgraph,
                            start=enter_node,
                            end=exit_node,
                            # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
                        ) """
                        block_path = longest_simple_path(
                            subgraph=subgraph,
                            reversed_subgraph=subgraph,
                            enter_node=enter_node,
                            exit_node=exit_node,
                            symmetry=symmetry,
                            encoding=encoding
                        )
                        """ if block_path[0] == exit_node:
                            block_path = list(reversed(block_path)) """
                    else:
                        """ block_path = longest_simple_path_binary_search(
                            G=subgraph,
                            start=enter_node,
                            end=exit_node,
                            # symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None,
                        ) """
                        block_path = longest_simple_path(
                            subgraph=subgraph,
                            reversed_subgraph=subgraph,
                            enter_node=enter_node,
                            exit_node=exit_node,
                            symmetry=symmetry,
                            encoding=encoding
                        )
                        """ if enter_node is not None and block_path[-1] == enter_node:
                            block_path = list(reversed(block_path))
                        if exit_node is not None and block_path[0] == exit_node:
                            block_path = list(reversed(block_path)) """

                    if not current_path:
                        current_path.extend(block_path)
                    else:
                        # Skip the first node to prevent duplicate entries of the cut nodes
                        current_path.extend(block_path[1:])

            if len(current_path) > len(longest_path):
                longest_path = current_path

        return longest_path


def longest_simple_path(G: nx.Graph, dp: bool = True, encoding: int = 0):
    """
    encoding [default: 0]:
        0 = position encoding (exactly k) - linear search (top-down)
        1 = position encoding (atleast k) variant 1 - linear search (bottom-up)
        2 = position encoding (atleast k) variant 2 - linear search (bottom-up)
        3 = position encoding (atleast k) variant 1 - binary search
        4 = position encoding (atleast k) variant 2 - binary search
        5 = edge encoding (atleast k) - linear search (bottom-up)
        6 = edge encoding (atleast k) - binary search
        7 = incremental position encoding
    """
    H = optimize(G)
    longest_path = []

    for C in H.graph.get("connected_components", [H]):
        path = longest_simple_path_components(C, dp, encoding)
        if len(path) > len(longest_path):
            longest_path = path

    return longest_path
