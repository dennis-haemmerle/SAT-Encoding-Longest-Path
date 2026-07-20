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
    else:
        connected_components = [H.subgraph(c).copy() for c in nx.connected_components(H)]

    for C in connected_components:
        bridges = list(nx.bridges(C))

        # blocks after removing the bridges
        graph_without_bridges = C.copy()
        graph_without_bridges.remove_edges_from(bridges)
        bridge_components = [graph_without_bridges.subgraph(c).copy() for c in nx.connected_components(graph_without_bridges)]

        # node -> block_id
        node_to_block = {}
        for block_id, nodes in enumerate(nx.connected_components(graph_without_bridges)):
            for n in nodes:
                node_to_block[n] = block_id

        attachments = [set() for _ in bridge_components]

        # the bridge tree connects the blocks with the bridges
        bridge_tree = nx.Graph()
        bridge_tree.add_nodes_from(range(len(bridge_components)))

        for u, v in bridges:
            block_u = node_to_block[u]
            block_v = node_to_block[v]

            attachments[block_u].add(u)
            attachments[block_v].add(v)

            bridge_tree.add_edge(block_u, block_v, bridge=(u, v), attach={block_u: u, block_v: v})

        # length of the longest bridge chain
        bridge_chain_length = 0 if bridge_tree.number_of_nodes() <= 1 else nx.diameter(bridge_tree)

        C.graph["bridges"] = bridges
        C.graph["bridge_components"] = bridge_components  # 2_edge_connected_components
        C.graph["bridge_tree"] = bridge_tree
        C.graph["bridge_chain_length"] = bridge_chain_length

        for block_id, block in enumerate(bridge_components):
            block.graph["symmetry"] = compute_symmetry_info(block, attachments[block_id])

            if block.is_directed():
                block.graph["only_in_nodes"] = [v for v in block.nodes() if block.in_degree(v) == 0]  # type: ignore
                block.graph["only_out_nodes"] = [v for v in block.nodes() if block.out_degree(v) == 0]  # type: ignore
            else:
                block.graph["leaves"] = [v for v in block.nodes() if block.degree(v) == 1]

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
