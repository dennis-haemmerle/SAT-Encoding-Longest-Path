import networkx as nx
import pynauty as pn
from collections import defaultdict


# preprocessing
def optimize(G: nx.Graph):
    H = G.copy()

    H.remove_edges_from(list(nx.selfloop_edges(H)))
    H.remove_nodes_from(list(nx.isolates(H)))

    if H.is_directed():
        connected_components = [H.subgraph(c).copy() for c in nx.weakly_connected_components(H)]  # type: ignore

        for C in connected_components:
            strongly_connected_components = list(nx.strongly_connected_components(C))  # type: ignore
            # the condensation DAG stores each strongly connected component in a node
            condensation = nx.condensation(C, strongly_connected_components)  # type: ignore

            # map original edges to edges between SCC nodes
            for u, v in condensation.edges():
                condensation.edges[u, v]["original_edges"] = []

            mapping = condensation.graph["mapping"]
            for u, v in C.edges():
                u_scc = mapping[u]
                v_scc = mapping[v]
                if u_scc != v_scc:
                    condensation.edges[u_scc, v_scc]["original_edges"].append((u, v))

            for scc_id in condensation.nodes:
                # store subgraph for each SCC node
                scc_nodes = list(condensation.nodes[scc_id]["members"])
                subgraph = C.subgraph(scc_nodes).copy()
                condensation.nodes[scc_id]["type"] = "strongly_connected_component"
                condensation.nodes[scc_id]["subgraph"] = subgraph

                # compute symmetry for each SCC node
                condensation.nodes[scc_id]["symmetry"] = compute_symmetry_info(subgraph)

            C.graph["condensation_dag"] = condensation
    else:
        connected_components = [H.subgraph(c).copy() for c in nx.connected_components(H)]

        for C in connected_components:
            # compute_bridge_connectivity_info(C)
            cut_nodes = list(nx.articulation_points(C))
            biconnected_components = list(nx.biconnected_components(C))

            block_cut_tree = nx.Graph()

            # add cut nodes to the block cut tree
            for node in cut_nodes:
                block_cut_tree.add_node(node, type="cut_node")

            # add biconnected blocks as nodes and store their corresponding subgraphs
            for block_id, nodes in enumerate(biconnected_components):
                block_cut_tree.add_node(f"block{block_id}", type="block", subgraph=C.subgraph(nodes).copy())

                # connect the block to all its incident cut nodes
                for n in nodes:
                    if n in cut_nodes:
                        block_cut_tree.add_edge(f"block{block_id}", n)

            # compute symmetry for each block
            for block in block_cut_tree.nodes():
                if block_cut_tree.nodes[block]["type"] == "block":
                    subgraph = block_cut_tree.nodes[block]["subgraph"]
                    block_cut_tree.nodes[block]["symmetry"] = compute_symmetry_info(subgraph)

            C.graph["block_cut_tree"] = block_cut_tree

    H.graph["connected_components"] = connected_components

    return H


def compute_symmetry_info(G: nx.Graph, attachment_nodes=None):
    attachment_nodes = set(attachment_nodes or [])

    H, node_mapping, reverse_mapping = nx_to_pynauty(G)

    generators, grpsize1, grpsize2, orbits, numorbits = pn.autgrp(H)

    orbit_groups = defaultdict(list)
    for i, orbit in enumerate(orbits):
        v = reverse_mapping[i]
        orbit_groups[orbit].append(v)

    return {
        "orbits": orbits,  # Orbit representative for each node index
        "orbit_groups": orbit_groups,  # Groups of symmetric nodes
        "numorbits": numorbits,  # Number of orbit groups
        "automorphism_count": grpsize1 * (10 ** grpsize2),  # Number of graph automorphisms
    }


def nx_to_pynauty(G: nx.Graph) -> tuple[pn.Graph, dict, dict]:
    node_mapping = {v: i for i, v in enumerate(G.nodes())}
    reverse_mapping = {i: v for v, i in node_mapping.items()}

    adjacency_dict = {}
    for v in G.nodes():
        adjacency_dict[node_mapping[v]] = {node_mapping[u] for u in G.neighbors(v)}

    H = pn.Graph(number_of_vertices=G.number_of_nodes(), directed=G.is_directed(), adjacency_dict=adjacency_dict)

    return H, node_mapping, reverse_mapping
