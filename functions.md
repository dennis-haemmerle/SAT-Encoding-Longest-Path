# Functions

## Functions to determine a simple path of length atleast or exactly k

Position encoding (exactly k): simple_path_of_length_k(G: nx.Graph, k: int)
Position encoding (atleast k) variant 1: simple_path_of_length_atleast_k(G: nx.Graph, k: int)
Position encoding (atleast k) variant 2: simple_path_of_length_atleast_k_2(G: nx.Graph, k: int)
Edge encoding (exactly k or atleast k): simple_path_of_length_k_edge_encoding(G: nx.Graph, k: int, atleast_k=False)

### Notes:
- start=None, end=None, and symmetry=None are passed by higher-level functions and can be ignored
- atleast k is faster then exactly k (edge encoding)
- variant 1 and variant 2 have similar runtimes (position encoding atleast k)


## Functions to determine the longest simple path

Linear Search (bottom up 1...n): longest_simple_path_linear_search(G: nx.Graph, edge_encoding=False, atleast_k=False, atleast_k_variant=False)
- edge_encoding [standard: False]: True -> edge encoding / False -> position encoding
- atleast_k [standard: False]: True -> atleast k / False -> exactly k
- atleast_k_variant [standard: False] (only needed if not edge_encoding and atleast_k): True -> variant 2 / False -> variant 1

Linear Search (top down n...1): longest_simple_path_linear_search_top_down(G: nx.Graph)

Binary Search: longest_simple_path_binary_search(G: nx.Graph, edge_encoding=False, atleast_k=False, atleast_k_variant=False)
- edge_encoding [standard: False]: True -> edge encoding / False -> position encoding
- atleast_k [standard: False]: True -> atleast k / False -> exactly k
- atleast_k_variant [standard: False] (only relevant if not edge_encoding and atleast_k): True -> variant 2 / False -> variant 1

Edge encoding (MaxSAT): longest_simple_path_edge_encoding_maxsat(G: nx.Graph)

Optimized (preprocessed): longest_simple_path(G: nx.Graph, dp: bool = True, encoding: int = 0)
- dp (only relevant for DiGraph) [standard: True]: True -> dynamic programming with topological order / False -> Brute Force Method
- encoding = {0, ..., 7} [standard: 0]: all the function combinations from above

### Notes:
- start=None, end=None, and symmetry=None are passed by higher-level functions and can be ignored
- longest_simple_path_components(C: nx.Graph) is used by longest_simple_path(G: nx.Graph)
- the runtime of MaxSAT explodes really fast
- the preprocessed encoding differs between directed and undirected graphs


## Comparisons

- Position encoding (exactly k) vs Edge encoding (exactly k)
- Position encoding (atleast k) variant 1 vs Position encoding (atleast k) variant 2
- Edge encoding (exactly k) vs Edge encoding (atleast k)
- Linear search vs Binary search vs Incremental search (position encoding)
- Linear search vs Binary search vs MaxSAT (edge encoding)
- Optimized Encoding: Dynamic Programming vs Brute Force
- Optimized vs not Optimized ?