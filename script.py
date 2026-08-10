import argparse
import sys
import simple_path


def run_incremental_encoder(G):
    encoder = simple_path.IncrementalSimplePathEncoder(G)
    path = encoder.longest_simple_path()
    encoder.delete()
    return path


CONFIGURATIONS = {
    # Linear Search
    1: ("Linear Search (Position Encoding, Exactly k)",
        lambda G: simple_path.longest_simple_path_linear_search(G=G, edge_encoding=False, atleast_k=False)),
    2: ("Linear Search (Position Encoding, Atleast k, variant 1)",
        lambda G: simple_path.longest_simple_path_linear_search(G=G, edge_encoding=False, atleast_k=True, atleast_k_variant=False)),
    3: ("Linear Search (Position Encoding, Atleast k, variant 2)",
        lambda G: simple_path.longest_simple_path_linear_search(G=G, edge_encoding=False, atleast_k=True, atleast_k_variant=True)),
    4: ("Linear Search (Edge Encoding, Exactly k)",
        lambda G: simple_path.longest_simple_path_linear_search(G=G, edge_encoding=True, atleast_k=False)),
    5: ("Linear Search (Edge Encoding, Atleast k)",
        lambda G: simple_path.longest_simple_path_linear_search(G=G, edge_encoding=True, atleast_k=True)),

    # Binary Search
    6: ("Binary Search (Position Encoding, Exactly k)",
        lambda G: simple_path.longest_simple_path_binary_search(G=G, edge_encoding=False, atleast_k=False)),
    7: ("Binary Search (Position Encoding, Atleast k, variant 1)",
        lambda G: simple_path.longest_simple_path_binary_search(G=G, edge_encoding=False, atleast_k=True, atleast_k_variant=False)),
    8: ("Binary Search (Position Encoding, Atleast k, variant 2)",
        lambda G: simple_path.longest_simple_path_binary_search(G=G, edge_encoding=False, atleast_k=True, atleast_k_variant=True)),
    9: ("Binary Search (Edge Encoding, Exactly k)",
        lambda G: simple_path.longest_simple_path_binary_search(G=G, edge_encoding=True, atleast_k=False)),
    10: ("Binary Search (Edge Encoding, Atleast k)",
         lambda G: simple_path.longest_simple_path_binary_search(G=G, edge_encoding=True, atleast_k=True)),

    # Special Methods
    11: ("Incremental Position Encoding",
         run_incremental_encoder),
    12: ("MaxSAT Edge Encoding",
         lambda G: simple_path.longest_simple_path_edge_encoding_maxsat(G)),

    # Preprocessed Undirected Graphs
    13: ("Preprocessed (Position Encoding, Exactly k - Linear Search, Top-Down)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=0)),
    14: ("Preprocessed (Linear Search, Exactly k - Position Encoding, Atleast k, variant 1)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=1)),
    15: ("Preprocessed (Linear Search, Exactly k - Position Encoding, Atleast k, variant 2)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=2)),
    16: ("Preprocessed (Binary Search, Exactly k - Position Encoding, Atleast k, variant 1)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=3)),
    17: ("Preprocessed (Binary Search, Exactly k - Position Encoding, Atleast k, variant 2)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=4)),
    18: ("Preprocessed (Edge Encoding, Atleast k, Linear Search)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=5)),
    19: ("Preprocessed (Edge Encoding, Atleast k, Binary Search)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=6)),
    20: ("Preprocessed (Incremental Position Encoding)",
         lambda G: simple_path.longest_simple_path(G=G, encoding=7)),

    # Preprocessed Directed Graphs (Brute Force)
    21: ("Preprocessed (Position Encoding, Exactly k - Linear Search, Top-Down)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=0)),
    22: ("Preprocessed (Linear Search, Exactly k - Position Encoding, Atleast k, variant 1)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=1)),
    23: ("Preprocessed (Linear Search, Exactly k - Position Encoding, Atleast k, variant 2)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=2)),
    24: ("Preprocessed (Binary Search, Exactly k - Position Encoding, Atleast k, variant 1)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=3)),
    25: ("Preprocessed (Binary Search, Exactly k - Position Encoding, Atleast k, variant 2)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=4)),
    26: ("Preprocessed (Edge Encoding, Atleast k, Linear Search)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=5)),
    27: ("Preprocessed (Edge Encoding, Atleast k, Binary Search)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=6)),
    28: ("Preprocessed (Incremental Position Encoding)",
         lambda G: simple_path.longest_simple_path(G=G, dp=False, encoding=7)),

    # Preprocessed Directed Graphs (Dynamic Programming)
    29: ("Preprocessed (Position Encoding, Exactly k - Linear Search, Top-Down)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=0)),
    30: ("Preprocessed (Linear Search, Exactly k - Position Encoding, Atleast k, variant 1)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=1)),
    31: ("Preprocessed (Linear Search, Exactly k - Position Encoding, Atleast k, variant 2)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=2)),
    32: ("Preprocessed (Binary Search, Exactly k - Position Encoding, Atleast k, variant 1)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=3)),
    33: ("Preprocessed (Binary Search, Exactly k - Position Encoding, Atleast k, variant 2)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=4)),
    34: ("Preprocessed (Edge Encoding, Atleast k, Linear Search)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=5)),
    35: ("Preprocessed (Edge Encoding, Atleast k, Binary Search)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=6)),
    36: ("Preprocessed (Incremental Position Encoding)",
         lambda G: simple_path.longest_simple_path(G=G, dp=True, encoding=7))
}


BASIC_IDS = [1, 5, 6, 10, 11, 12]
PREPROCESSED_IDS = list(range(13, 37))


def print_result(description, path):
    print(description)
    if path is not None:
        print(f"len {len(path) - 1 if len(path) > 0 else 0}")
        for node in path:
            print(node)
    else:
        print("no path found")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph_file", type=str)
    parser.add_argument("func", type=int, choices=list(CONFIGURATIONS.keys()), nargs="?", default=None)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--basic", action="store_true")
    parser.add_argument("--preprocessed", action="store_true")

    args = parser.parse_args()

    try:
        G = simple_path.read_from_file(args.graph_file)
    except Exception as e:
        print(f"Error loading '{args.graph_file}': {e}", file=sys.stderr)
        sys.exit(1)

    if args.all:
        for func_id, (desc, func) in CONFIGURATIONS.items():
            path = func(G)
            print_result(f"[{func_id}] {desc}", path)
    elif args.basic:
        for func_id in BASIC_IDS:
            desc, func = CONFIGURATIONS[func_id]
            path = func(G)
            print_result(f"[{func_id}] {desc}", path)
    elif args.preprocessed:
        for func_id in PREPROCESSED_IDS:
            desc, func = CONFIGURATIONS[func_id]
            path = func(G)
            print_result(f"[{func_id}] {desc}", path)
    elif args.func:
        desc, func = CONFIGURATIONS[args.func]
        path = func(G)
        print_result(f"[{args.func}] {desc}", path)


if __name__ == "__main__":
    main()
