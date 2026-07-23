import networkx as nx
from itertools import combinations
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from optimizations import optimize


def simple_path_of_length_k(G: nx.Graph, k: int, start=None, end=None, only_in=None, only_out=None, leaves=None, symmetry=None):
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


def simple_path_of_length_k_edge_encoding(G: nx.Graph, k: int):
    if k < 0 or k > G.number_of_edges():
        return None
    if k == 0:
        return [] if G.number_of_nodes() > 0 else None

    vpool = IDPool()
    cnf = CNF()

    def edge_var(e):
        u, v = e
        if G.is_directed():
            return vpool.id((u, v))
        return vpool.id((min(u, v), max(u, v)))

    # 1. Exactly k edges
    lits = [edge_var(e) for e in G.edges()]
    cnf.extend(CardEnc.equals(lits=lits, bound=k, vpool=vpool, encoding=EncType.seqcounter).clauses)

    # 2. Exactly k+1 nodes
    lits = [vpool.id(v) for v in G.nodes()]
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
    else:
        endpoint_vars = []

        for v in G.nodes():
            incident = [edge_var(e) for e in G.edges(v)]

            # end_var = True <=> deg(v) = 1
            end_var = vpool.id(("end", v))
            # endpoint_vars.append(end_var)
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
            assignment = [v for v in G.nodes() if vpool.id(v) in model]
            return assignment
        return None


def longest_simple_path_linear_search(G: nx.Graph, start=None, end=None, only_in=None, only_out=None, leaves=None, symmetry=None):
    longest_path = []

    for k in range(1, G.number_of_nodes()):
        path = simple_path_of_length_k(G, k, start, end, only_in, only_out, leaves, symmetry)

        if path is not None:
            longest_path = path
        else:
            break

    return longest_path


def longest_simple_path_linear_search_top_down(G: nx.Graph, start=None, end=None, only_in=None, only_out=None, leaves=None, symmetry=None):
    for k in range(G.number_of_nodes() - 1, 0, -1):
        path = simple_path_of_length_k(G, k, start, end, only_in, only_out, leaves, symmetry)

        if path is not None:
            return path

    return []


def longest_simple_path_binary_search(G: nx.Graph, start=None, end=None, only_in=None, only_out=None, leaves=None, symmetry=None):
    longest_path = []
    low = 0
    high = G.number_of_nodes() - 1

    while low <= high:
        mid = (low + high) // 2
        path = simple_path_of_length_k(G, mid, start, end, only_in, only_out, leaves, symmetry)

        if path is not None:
            longest_path = path
            low = mid + 1
        else:
            high = mid - 1

    return longest_path


def longest_simple_path_components(C: nx.Graph):
    if C.is_directed():
        dag = C.graph["condensation_dag"]

        if dag.number_of_nodes() == 0:
            return []
        if dag.number_of_nodes() == 1:
            scc = next(iter(dag.nodes))
            subgraph = dag.nodes[scc]["subgraph"]
            return longest_simple_path_binary_search(subgraph)

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
                                    scc_path = longest_simple_path_linear_search_top_down(G=subgraph, start=enter_node, end=exit_node)
                            else:
                                scc_path = longest_simple_path_binary_search(G=subgraph, start=enter_node, end=exit_node)

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
        block_cut_tree = C.graph["block_cut_tree"]
        blocks = [n for n, attr in block_cut_tree.nodes(data=True) if attr["type"] == "block"]

        if len(blocks) == 0:
            return []
        if len(blocks) == 1:
            subgraph = block_cut_tree.nodes[blocks[0]]["subgraph"]
            symmetry = block_cut_tree.nodes[blocks[0]]["symmetry"]
            return longest_simple_path_binary_search(
                G=subgraph,
                symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
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
                        block_path = longest_simple_path_linear_search_top_down(
                            G=subgraph,
                            start=enter_node,
                            end=exit_node,
                            symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
                        )
                    else:
                        block_path = longest_simple_path_binary_search(
                            G=subgraph,
                            start=enter_node,
                            end=exit_node,
                            symmetry=symmetry if symmetry.get("numorbits") < subgraph.number_of_nodes() else None
                        )

                    if not current_path:
                        current_path.extend(block_path)
                    else:
                        # Skip the first node to prevent duplicate entries of the cut nodes
                        current_path.extend(block_path[1:])

            if len(current_path) > len(longest_path):
                longest_path = current_path

        return longest_path


def longest_simple_path(G: nx.Graph):
    H = optimize(G)
    longest_path = []

    for C in H.graph.get("connected_components", [H]):
        path = longest_simple_path_components(C)
        if len(path) > len(longest_path):
            longest_path = path

    return longest_path
